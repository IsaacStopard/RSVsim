#' Runs the RSV model
#'
#' Function to run the transmission model with cohort aging.
#' Cohorts are aged at time intervals of the smallest age group given in \code{parameters[["age.limits"]]}.
#'
#' @param parameters List of parameters from \code{get_params} function.
#' @param max_t Simulation maximum time. Default: 2000 days.
#' @param dt Time steps to run the model over. Default: 0.25 days.
#' @param init_conds Initial conditions to run the model. List. Default: \code{NULL}. If \code{NULL}: 1% RSV prevalence is assumed for people during the primary infection.
#' All other people are assumed to be susceptible to their primary infection.
#' @return Simulation output (dataframe). In the dataframe, age refers to the lowest age in the age group.
#' @export
run_model <- function(parameters,
                      max_t = 2000,
                      dt = 0.25,
                      init_conds = NULL
                      ){

  ##########################################################
  # labels for the ages
  age_chr <- c(paste0("[",round(parameters$age.limits[1], digits = 2),",", round(parameters$age.limits[2], digits = 2),"]"))

  for(i in 2:(parameters$nAges - 1)){
    age_chr <- c(age_chr, c(paste0("[",round(parameters$age.limits[i], digits = 2),",", round(parameters$age.limits[i+1], digits = 2),")")))
  }

  age_chr <- c(age_chr, paste0("[",round(parameters$age.limits[parameters$nAges], digits = 2),",", round(parameters$max_age, digits = 2),")"))

  # age differences in days
  size_cohorts <- with(parameters,
                       c(diff(age.limits * 365.25), max_age*365.25 - age.limits[length(age.limits)]*365.25)
                       )

  # running for the time of the smallest cohort
  cohort_step_size <- floor(min(size_cohorts))

  transition_rate <- cohort_step_size/size_cohorts

  rel_sizes <- size_cohorts/sum(size_cohorts)

  n_steps <- floor(max_t / cohort_step_size)

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
        R0 = rep_z,
        DetIncidence0 = rep_z,
        Incidence0 = rep_z
      )
  } else{
    parameters <- with(init_conds,{
      for(name in c("Sp0", "Ep0", "Ip0", "Ss0", "Es0", "Is0",
                    "R0", "DetIncidence0", "Incidence0")){
        if(length(get(name)) != parameters$nAges){
          stop(paste("Initial conditions for", name,"are not all the same length as the number of age categories", sep = " "))
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
        R0 = R0,
        DetIncidence0 = DetIncidence0,
        Incidence0 = Incidence0
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

  incidence_i <- which(states %in% c("Incidence", "DetIncidence"))

  # running the model with cohort aging (run for a single cohort, move cohort, change initial states, repeat)
  out_list <- vector(mode = "list", length = n_steps)

  for(i in 1:n_steps){

    times_in <- if(i == 1){
      seq(0, i * cohort_step_size, dt)
    } else{
      seq((i - 1) * cohort_step_size + dt, i * cohort_step_size, dt)
    }

    out <- dust2::dust_system_simulate(RSV_dust,
                                       times = times_in)

    # checking the total population is correct
    if(!all(abs(out[-incidence_i,] |> colSums() - parameters$total_population) < 1E-5)){
      stop("Population does not sum to the correct number")
    }

    out_list[[i]] <- out |> as.data.frame()

    colnames(out_list[[i]]) <- times_in

    out_list[[i]] <- out_list[[i]] |> dplyr::mutate(state = states, age = ages) |>
      tidyr::pivot_longer(cols = 1:length(times_in),
                          values_to = "value",
                          names_to = "time") |>
      dplyr::mutate(time = as.numeric(time)) |>
      dplyr::arrange(factor(state,
                            levels = c("Sp", "Ep", "Ip", "Ss", "Es", "Is", "R", "Incidence", "DetIncidence")),
                            age)

    # calculating the initial states
    next_state <- out_list[[i]] |>
      dplyr::filter(time == max(time)) |>
      dplyr::group_by(state) |>
      dplyr::mutate(transition_ct = value * transition_rate,
                    lag_transition = dplyr::lag(transition_ct, 1, order_by = state))

    # check NAs are in the correct place
    if(!all(
      next_state[which(is.na(next_state$lag_transition)), "age"] == 0) |
       sum(is.na(next_state$lag_transition)) != n_states){
      stop("error when lagging the cohort aging")
    }

    # filling in births to keep the population size constant
    next_state <- next_state |> dplyr::ungroup() |>
      dplyr::mutate(lag_transition =
                      ifelse(is.na(lag_transition) & age == 0 & state == "Sp",
                             rel_sizes[1] * transition_rate[1] * parameters$total_population,
                             ifelse(is.na(lag_transition) & age == 0, 0,
                                    lag_transition)
                             )
                    )

    if(
      abs(
        next_state |> dplyr::filter(age == max(age) &
                                      state %in% c("Sp", "Ep", "Ip", "Ss", "Es", "Is", "R")
                                      ) |> dplyr::ungroup() |>
          dplyr::select(transition_ct) |> as.vector() |> unlist() |> sum() -
        rel_sizes[1] * transition_rate[1] * parameters$total_population
        ) > 1E-5){
      stop("births are not correct to keep the population size constant")
    }

    next_state <- next_state |>
      dplyr::mutate(next_value = value - transition_ct + lag_transition) |>
      dplyr::select(next_value) |>
      unlist() |> unname() |> as.vector()

    if(abs(sum(next_state[-incidence_i]) - parameters$total_population) > 1E-7){
      stop("total population for the next cohort states is not correct")
    }

    dust2::dust_system_set_state(sys = RSV_dust,
                                 state = next_state)

    if(all(dust2::dust_system_state(RSV_dust) != next_state)){
      stop("states have not been updated")
    }

  }

  return(
    invisible(
      dplyr::bind_rows(out_list) |>
      dplyr::left_join(data.frame(age = parameters$age.limits,
                                  age_chr = age_chr), by = dplyr::join_by(age)) |>
      dplyr::select(-age) |>
      dplyr::rename(age = age_chr) |> as.data.frame()
      )
    )
}

