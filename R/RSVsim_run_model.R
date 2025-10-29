#' Runs the RSV model using odin and dust
#'
#' Function to run the transmission model with cohort aging. Prevalence, incidence between the given time-steps and the incidence rate per day are also calculated.
#'
#' @param parameters List of parameters from \code{RSVsim_parameters} function.
#' @param times Simulation times. Default: 0 to 365 days with intervals of 0.25 days.
#' @param cohort_step_size Time steps to run the model over before adjusting the ages of all cohorts. Default: 60 days. If \code{is.numeric(cohort_step_size) == FALSE} then cohort ageing is not applied. Can have a maximum of 3 decimal places.
#' If \code{NULL}: 0.1% RSV prevalence is assumed for people during the primary infection, which is seeded at the beginning of the simulation.
#' All other people are assumed to be susceptible to their primary infection.
#' @param warm_up Length of time-points to exclude before calculating the likelihood. Default: \code{NULL}.
#' @return Simulation output (dataframe). In the dataframe, age refers to the lowest age in the age group.
#' @export
RSVsim_run_model <- function(parameters,
                             times = seq(0, 365*1, 0.25),
                             cohort_step_size = 60,
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

  states <- lapply(1:length(states), function(i){rep(names(states[i]), length(states[[i]]))}) |> unlist()

  n_states <- length(states_order)

  ages <- rep(parameters$age.limits, n_states)

  # no cohort aging in these variables
  output_variable_names <- c("Incidence", "DetIncidence", "prev", "prev_p", "prev_s", "Incidence_rate", "DetIncidence_rate")

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

      out_list[[i]] <- out_list[[i]] |> dplyr::mutate(state = states, age = ages) |>
        tidyr::pivot_longer(cols = 1:length(times_in[[i]]),
                            values_to = "value",
                            names_to = "time") |>
        dplyr::mutate(time = as.numeric(time)) |>
        dplyr::arrange(factor(state,
                              levels = states_order), age)

      # calculating the initial states with cohort ageing in the S, E, I, R compartments only
      next_state <- out_list[[i]] |>
        dplyr::filter(time == max(time) & !(state %in% output_variable_names)) |>
        dplyr::group_by(state) |>
        dplyr::mutate(transition_ct = value * transition_rate,
                      lag_transition = dplyr::lag(transition_ct, 1, order_by = state))

    # check NAs are in the correct place
    #if(!all(
    #  next_state[which(is.na(next_state$lag_transition)), "age"] == 0) |
    #   sum(is.na(next_state$lag_transition)) != n_states){
    #  stop("RSVsim_run_model: error when lagging the cohort aging")
    #}

      # filling in births to keep the population size constant
      next_state <- next_state |> dplyr::ungroup() |>
        dplyr::mutate(lag_transition =
                        ifelse(is.na(lag_transition) & age == 0 & state == "Sp",
                               rel_sizes[1] * transition_rate[1] * parameters$total_population,
                               ifelse(is.na(lag_transition) & age == 0, 0, lag_transition)
                               ),
                      next_value = value - transition_ct + lag_transition
                      )

      next_state_output_variables = out_list[[i]] |> dplyr::filter(time == max(time) & state %in% output_variable_names) |> dplyr::rename(next_value = value)

      next_state <- rbind(next_state[,colnames(next_state_output_variables)], next_state_output_variables) |>
        dplyr::arrange(factor(state,
                              levels = states_order), age) |>
        dplyr::select(next_value) |>
        unlist() |> unname() |> as.vector()

      dust2::dust_system_set_state(sys = RSV_dust, state = next_state)

      # if(all(dust2::dust_system_state(RSV_dust) != next_state)){
      #   stop("RSVsim_run_model: states have not been updated")
      # }

      out_list[[i]] <- out_list[[i]] |> tidyr::pivot_wider(names_from = state, values_from = value)
      }

    out <- invisible(
      dplyr::bind_rows(out_list) |>
        dplyr::left_join(data.frame(age = parameters$age.limits, age_chr = parameters$age_chr), by = dplyr::join_by(age))
      )

  } else{

    out <- dust2::dust_system_simulate(RSV_dust, times = times) |> as.data.frame()
    colnames(out) <- times
    out <- out |> dplyr::mutate(state = states, age = ages) |>
      tidyr::pivot_longer(cols = 1:length(times),
                          values_to = "value",
                          names_to = "time") |>
      dplyr::mutate(time = as.numeric(time)) |>
      dplyr::arrange(factor(state,
                            levels = states_order), age) |>
      tidyr::pivot_wider(names_from = state, values_from = value)
  }

  # incidence calculation
  out_checkout <- out |> dplyr::group_by(age) |> dplyr::mutate(Incidence = tidyr::replace_na(Incidence - dplyr::lag(Incidence, 1), 0),
                                                               DetIncidence = tidyr::replace_na(DetIncidence - dplyr::lag(DetIncidence, 1), 0)
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

