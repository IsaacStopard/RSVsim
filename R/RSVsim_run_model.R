#' Runs the RSV model using odin and dust
#'
#' Function to run the transmission model with cohort aging. Prevalence, incidence between the given time-steps and the incidence rate per day are also calculated.
#'
#' @param parameters List of parameters from \code{RSVsim_parameters} function.
#' @param times Simulation times. Default: 0 to 365.25 days with intervals of 0.25 days.
#' @param cohort_step_size Time steps to run the model over before adjusting the ages of all cohorts. Default: 1 month. If \code{is.numeric(cohort_step_size) == FALSE} then cohort ageing is not applied. Can have a maximum of 3 decimal places.
#' @param warm_up Length of time-points to exclude before calculating the likelihood. Default: \code{NULL}.
#' @return Simulation output (dataframe). In the dataframe, age refers to the lowest age in the age group.
#' @export
RSVsim_run_model <- function(parameters,
                             times = seq(0, 365.25*1, 0.25),
                             cohort_step_size = 1/12 * 365.25,
                             warm_up = NULL
                             ){

  ##########################################################

  size_cohorts <- parameters$size_cohorts

  rel_sizes <- parameters$rel_sizes

  # running for the time of the smallest cohort
  if(is.numeric(cohort_step_size) & round(max(diff(times)), digits = 5) >= round(cohort_step_size, digits = 5)){
    stop("The maximum time difference is greater than or equal to the cohort step size")
  }

  if(is.numeric(cohort_step_size) & round(min(size_cohorts), digits = 5) < round(cohort_step_size, digits = 5)){
    stop("The smallest cohort age size is smaller than the cohort_step_size: increase the differences in age limits or decrease the cohort_step_size")
  }

  max_t <- max(times)

  # runs the model with cohort ageing
  RSV_dust <- dust2::dust_system_create(generator = RSV_ODE,
                                        pars = parameters,
                                        n_particles = 1,
                                        n_groups = 1,
                                        time = 0,
                                        deterministic = TRUE
                                        ) |> invisible()


  dust2::dust_system_set_state_initial(RSV_dust)

  # do not change this
  # the order is determined by the order in the odin RSV_ODE.R file
  states <- dust2::dust_unpack_index(RSV_dust)

  states_order <- names(states)

  # creating the infection state, vaccination state and age corresponding to each position in dust
  states <- rep(states_order, lengths(states))

  # full_groups <- c("Sp", "Ss", "Ep", "Es", "Ip", "Is", "R", "Incidence", "DetIncidence", "Incidence_rate", "DetIncidence_rate", "prev", "prev_p", "prev_s","doses")
  # age_groups <- c()
  # vac_groups <- c()

  # states_order_group <- dplyr::case_when(states_order %in% full_groups ~ "full", states_order %in% age_groups ~ "age", states_order %in% vac_groups ~ "vac", .default = "other")

  # age_lookup <- list("full" = rep(parameters$age.limits, parameters$nVaccStates), "age" = parameters$age.limits, "vac" = rep(NA, parameters$nVaccStates), "other" = NA)
  # vac_lookup <- list("full"  = rep(1:parameters$nVaccStates, each = parameters$nAges), "age" = rep(NA, parameters$nAges), "vac" = 1:parameters$nVaccStates, "other" = NA)

  # ages <- unlist(age_lookup[states_order_group], use.names = FALSE)
  # vacc_states <- unlist(vac_lookup[states_order_group], use.names = FALSE)

  n_states <- length(states_order)

  # important for maintaining the same order as RSV_ODE
  ages <- rep(rep(parameters$age.limits, parameters$nVaccStates), n_states)
  vacc_states <- rep(rep(1:parameters$nVaccStates, each = parameters$nAges), n_states)

  # no cohort aging in these variables
  output_variable_names <- states_order[!(states_order %in% c("Sp", "Ep", "Ip", "Ss", "Es", "Is", "R"))]

  dust_df <- data.frame("state" = factor(states, levels = states_order),
                        "age" = ages,
                        "vacc_state" = vacc_states
                        )

  if(is.numeric(cohort_step_size)){

    transition_rate <- cohort_step_size/size_cohorts

    n_steps <- ceiling(max_t / cohort_step_size)

    # running the model with cohort aging (run for a single cohort, move cohort, change initial states, repeat)
    out_list <- vector(mode = "list", length = n_steps)

    times_all <- sort(unique(round(c(times, 1:n_steps * cohort_step_size), digits = 10)))

    times_in <- lapply(1:n_steps, FUN = function(i){
      c(times_all[times_all >= ((i - 1) * cohort_step_size) & times_all < (i * cohort_step_size)])
      })

    for(i in 1:n_steps){

      out <- dust2::dust_system_simulate(RSV_dust, times = times_in[[i]])

      out_list[[i]] <- out |> as.data.frame()

      colnames(out_list[[i]]) <- times_in[[i]]

      out_list[[i]] <- out_list[[i]] |> cbind(dust_df) |>
        tidyr::pivot_longer(cols = 1:length(times_in[[i]]),
                            values_to = "value",
                            names_to = "time") |>
        dplyr::mutate(time = as.numeric(time))

      # calculating the initial states with cohort ageing in the S, E, I, R compartments only
      # assumes that ages are sequential

      next_state <- out_list[[i]] |>
        dplyr::filter(time == max(time)) |>
        dplyr::arrange(state, vacc_state, age)

      # filling in births to keep the population size constant
      next_state_age <- next_state |>
        dplyr::filter(!(state %in% output_variable_names)) |>
        dplyr::group_by(state, vacc_state) |>
        dplyr::mutate(transition_ct = value * transition_rate,
                      lag_transition = dplyr::lag(transition_ct, 1, order_by = state),
                      lag_transition = # everyone is born unvaccinated
                        ifelse(is.na(lag_transition) & age == 0 & state == "Sp" & vacc_state == 1,
                               rel_sizes[1] * transition_rate[1] * parameters$total_population,
                               ifelse(is.na(lag_transition) & age == 0, 0, lag_transition)),
                      next_value = value - transition_ct + lag_transition
                      ) |>
        dplyr::ungroup()

      next_state_output_variables <- next_state |> dplyr::filter(state %in% output_variable_names) |> dplyr::rename(next_value = value)

      next_state_in <- dust_df |>
        dplyr::left_join(rbind(next_state_age[,colnames(next_state_output_variables)], next_state_output_variables),
                         by = c("state", "vacc_state", "age")) |>
        dplyr::select(next_value) |>
        unlist() |> unname() |> as.vector()

      dust2::dust_system_set_state(sys = RSV_dust, state = next_state_in)

      # if(all(dust2::dust_system_state(RSV_dust) != next_state)){
      #   stop("RSVsim_run_model: states have not been updated")
      # }
      }

    out <- invisible(
      dplyr::bind_rows(out_list) |> tidyr::pivot_wider(names_from = state, values_from = value) |>
        dplyr::left_join(data.frame(age = parameters$age.limits, age_chr = parameters$age_chr), by = dplyr::join_by(age))
      )

  } else{

    out <- dust2::dust_system_simulate(RSV_dust, times = times) |> as.data.frame()
    colnames(out) <- times
    out <- out |> cbind(dust_df) |>
      tidyr::pivot_longer(cols = 1:length(times),
                          values_to = "value",
                          names_to = "time") |>
      dplyr::mutate(time = as.numeric(time)) |>
      dplyr::arrange(factor(state,
                            levels = states_order), vacc_state, age) |>
      tidyr::pivot_wider(names_from = state, values_from = value)
  }

  # incidence calculation
  out_checkout <- out |> dplyr::group_by(age, vacc_state) |> dplyr::mutate(cumulativeIncidence = Incidence,
                                                                           cumulativeDetIncidence = DetIncidence,
                                                                           cumulativeDoses = doses,
                                                                           Incidence = tidyr::replace_na(Incidence - dplyr::lag(Incidence, 1), 0),
                                                                           DetIncidence = tidyr::replace_na(DetIncidence - dplyr::lag(DetIncidence, 1), 0),
                                                                           Doses = tidyr::replace_na(doses - dplyr::lag(doses, 1), 0)
                                                                           ) |>
    dplyr::ungroup()

  if(!is.null(warm_up)){
      out_checkout <- out_checkout |> dplyr::filter(time >= warm_up) |> dplyr::mutate(time = time - warm_up)
  }

  # checking the total population is correct


  if(any(abs(out_checkout |> dplyr::group_by(time) |> dplyr::summarise(total = sum(Sp) + sum(Ep) + sum(Ip) + sum(Ss) + sum(Es) + sum(Is) + sum(R)) |> dplyr::select(total) - parameters$total_population) > 1E-1)){
    stop("RSVsim_run_model: population does not sum to the correct number")
  }

  return(out_checkout |> as.data.frame())
}

