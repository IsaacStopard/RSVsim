#' Runs the RSV model using odin and dust
#'
#' Function to run the transmission model with cohort aging.
#'
#' @param parameters List of parameters from \code{get_params} function.
#' @param times Simulation times. Default: 0 - 3650 days with intervals of 0.25 days.
#' @param cohort_step_size Time steps to run the model over before adjusting the ages of all cohorts. Default: 10 days.
#' @param init_conds Initial conditions to run the model. List. Default: \code{NULL}. If \code{NULL}: 1% RSV prevalence is assumed for people during the primary infection.
#' All other people are assumed to be susceptible to their primary infection.
#' @param warm_up Length of time-points to exclude before calculating the likelihood. Default: \code{NULL}.

#' @return Simulation output (dataframe). In the dataframe, age refers to the lowest age in the age group.
#' @export
RSVsim_run_model_dust <- function(parameters,
                                  times = seq(0, 3650, 0.25),
                                  cohort_step_size = 10,
                                  init_conds = NULL,
                                  warm_up = NULL
                             ){

  ##########################################################
  # labels for the ages

  model <- stan_package_model(name = "RSV_ODE_stan", package = "RSVsim")

  age_chr <- c()

  for(i in 1:(parameters$nAges - 1)){
    age_chr <- c(age_chr, c(paste0("[",round(parameters$age.limits[i], digits = 2),",", round(parameters$age.limits[i+1], digits = 2),")")))
  }

  age_chr <- c(age_chr, paste0("[",round(parameters$age.limits[parameters$nAges], digits = 2),",", round(parameters$max_age, digits = 2),")"))

  # age differences in days
  size_cohorts <- with(parameters,
                       c(diff(age.limits * 365.25), max_age*365.25 - age.limits[length(age.limits)]*365.25)
                       )

  # running for the time of the smallest cohort
  if(max(diff(times)) >= cohort_step_size){
    stop("The maximum time difference is greater than or equal to the cohort step size")
  }

  max_t <- max(times)

  transition_rate <- cohort_step_size/size_cohorts

  rel_sizes <- size_cohorts/sum(size_cohorts)

  n_steps <- ceiling(max_t / cohort_step_size)

  # initial conditions: choose whether to reset here
  if(is.null(init_conds)) {

    rep_z <- rep(0, parameters$nAges)

    parameters <- purrr::list_modify(
      parameters,
        Sp0 = rel_sizes * parameters$total_population * 0.99,
        Ep0 = rep_z,
        Ip0 = rel_sizes * parameters$total_population * 0.01,
        Ss0 = rep_z,
        Es0 = rep_z,
        Is0 = rep_z,
        R0 = rep_z
      )
  } else{
    parameters <- with(init_conds,{
      for(name in c("Sp0", "Ep0", "Ip0", "Ss0", "Es0", "Is0", "R0")){
        if(length(get(name)) != parameters$nAges){
          stop(paste("RSVsim_run_model: initial conditions for", name,"are not all the same length as the number of age categories", sep = " "))
        }
      }
      purrr::list_modify(
        parameters,
        Sp0 = Sp0,
        Ep0 = Ep0,
        Ip0 = Ip0,
        Ss0 = Ss0,
        Es0 = Es0,
        Is0 = Is0,
        R0 = R0
      )
      }
    )
  }

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

  states <- lapply(1:length(states), function(i){rep(names(states[i]), length(states[[i]]))}) |> unlist()

  n_states <- length(unique(states))

  ages <- rep(parameters$age.limits, n_states)

  incidence_i <- which(states %in% c("Incidence", "DetIncidence", "prev", "prev_p", "prev_s"))

  # running the model with cohort aging (run for a single cohort, move cohort, change initial states, repeat)
  out_list <- vector(mode = "list", length = n_steps)

  for(i in 1:n_steps){

    times_in <- if(i == 1){
      sort(unique(c(0,
                    times[times >= 0 & times <= i * cohort_step_size],
                    i * cohort_step_size))
           )
    } else{
      sort(unique(c(times[times > ((i - 1) * cohort_step_size) & times <= (i * cohort_step_size)],
                    i * cohort_step_size)
                  )
           )
    }

    out <- dust2::dust_system_simulate(RSV_dust,
                                       times = times_in)

    # checking the total population is correct
    if(!all(abs(out[-incidence_i,] |> colSums() - parameters$total_population) < 1E-5)){
      stop("RSVsim_run_model: population does not sum to the correct number")
    }

    out_list[[i]] <- out |> as.data.frame()

    colnames(out_list[[i]]) <- times_in

    out_list[[i]] <- out_list[[i]] |> dplyr::mutate(state = states, age = ages) |>
      tidyr::pivot_longer(cols = 1:length(times_in),
                          values_to = "value",
                          names_to = "time") |>
      dplyr::mutate(time = as.numeric(time)) |>
      dplyr::arrange(factor(state,
                            levels = c("Sp", "Ep", "Ip", "Ss", "Es", "Is", "R", "Incidence", "DetIncidence", "prev", "prev_p", "prev_s")),
                            age)

    # calculating the initial states
    next_state <- out_list[[i]] |>
      dplyr::filter(time == max(time)) |>
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
                             ifelse(is.na(lag_transition) & age == 0, 0,
                                    lag_transition)
                             )
                    )

    next_state <- next_state |>
      dplyr::mutate(next_value = value - transition_ct + lag_transition) |>
      dplyr::select(next_value) |>
      unlist() |> unname() |> as.vector()

    dust2::dust_system_set_state(sys = RSV_dust,
                                 state = next_state)

    if(all(dust2::dust_system_state(RSV_dust) != next_state)){
      stop("RSVsim_run_model: states have not been updated")
    }

    out_list[[i]] <- out_list[[i]] |> tidyr::pivot_wider(names_from = state, values_from = value)

  }

  out <- invisible(
    dplyr::bind_rows(out_list) |>
      dplyr::left_join(data.frame(age = parameters$age.limits,
                                  age_chr = age_chr), by = dplyr::join_by(age)) |>
      as.data.frame()
      )

  if(!is.null(warm_up)){
    if(!is.numeric(warm_up)){
      stop("warm_up is not numeric")
    } else{
      out <- out |> dplyr::filter(time >= warm_up) |> dplyr::mutate(time = time - warm_up)
    }
  }

  return(out)
}

