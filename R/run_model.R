#' Runs the RSV model
#'
#' Function to run the transmission model with cohort aging.
#'
#' @param parameters List of parameters from \code{get_params} function.
#' @param max_t Simulation maximum time. Default: 2000 days.
#' @param dt Time steps to run the model over. Default: 0.25 days.
#' @param init_conds Initial conditions to run the model. List. Default: \code{NULL}. If \code{NULL} 1% RSV prevalence is assumed for people during the primary infection.
#' All other people are assumed to be susceptible to their primary infection.
#' @param save_final_states Boolean. Choose whether to save final model state as initial conditions for next simulation. Default: \code{FALSE}.
#' @return Simulation output
#' @export
run_model <- function(parameters,
                      max_t = 2000,
                      dt = 0.25,
                      init_conds = NULL,
                      save_final_states = FALSE #
                      ){

  ##########################################################

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

  # do not change this
  # the order is determined by the order in the odin RSV_ODE.R file
  states <- dust2::dust_unpack_index(RSV_dust)
  states <- lapply(1:length(states), function(i){rep(names(states[i]), length(states[[i]]))}) |> unlist()

  ages <- rep(parameters$age.limits, length(unique(states)))

  dust2::dust_system_set_state_initial(RSV_dust)

  out_df <- vector(mode = "list", length = n_steps)

  for(i in 1:n_steps){

    times_in <- if(i == 1){
      seq(0, i * cohort_step_size, dt)
    } else{
      seq((i - 1) * cohort_step_size + dt, i * cohort_step_size, dt)
    }

    out <- dust2::dust_system_simulate(RSV_dust,
                                times = times_in)

    out_df[[i]] <- out |> as.data.frame()

    colnames(out_df[[i]]) <- times_in

    out_df[[i]] <- out_df[[i]] |> dplyr::mutate(state = states, age = ages) |>
      tidyr::pivot_longer(cols = 1:length(times_in),
                          values_to = "value",
                          names_to = "time") |>
      dplyr::mutate(time = as.numeric(time)) |>
      dplyr::arrange(factor(state,
                            levels = c("Sp", "Ep", "Ip", "Ss", "Es", "Is", "Incidence", "DetIncidence")),
                            age)

    next_state <- out_df[[i]] |>
      dplyr::filter(time == max(time)) |>
      dplyr::group_by(state) |>
      dplyr::mutate(transition_ct = value * transition_rate,
                    lag_transition = dplyr::lag(transition_ct, 1, order_by = state))

    # check NAs are in the correct place
    if(!all(
      next_state[which(is.na(next_state$lag_transition)), "age"] == 0) |
       sum(is.na(next_state$lag_transition)) != 9){
      stop("error calculating the cohort aging")
    }

    next_state <- next_state |> dplyr::ungroup() |>
      dplyr::mutate(lag_transition =
                      ifelse(is.na(lag_transition) & age == 0 & state == "Sp",
                             rel_sizes[1] * transition_rate[1] * parameters$total_population,
                             ifelse(is.na(lag_transition) & age == 0,
                                    0, lag_transition)))

    if(
      abs(next_state |> dplyr::filter(age == max(age) & state %in% c("Sp", "Ep", "Ip", "Ss", "Es", "Is", "R")) |>
          dplyr::ungroup() |> dplyr::select(transition_ct) |> as.vector() |> unlist() |> sum() -
      rel_sizes[1] * transition_rate[1] * parameters$total_population) > 1E-5){
      stop("births are not correct to keep the population size constant")
    }

    next_state <- next_state |>
      dplyr::mutate(next_value = value - transition_ct + lag_transition) |>
      dplyr::ungroup() |> dplyr::select(next_value) |>
      unlist() |> unname() |> as.vector()

    dust2::dust_system_set_state(sys = RSV_dust,
                                 state = next_state)
  }

  while (T0 <= max_t){

    # solve the odes first
    t <- seq(from = T0, to = T0 + 1, by = dt)
    m <- mod$run(t)
    pop <- mod$transform_variables(m)

    if (T0 == 0){
      pop_out <- pop
    } else {
      pop_out$time <- c(pop_out$time, pop$time[5])
      pop_out$S <- rbind(pop_out$S, pop$S[5,])
      pop_out$E <- rbind(pop_out$E, pop$E[5,])
      pop_out$I <- rbind(pop_out$I, pop$I[5,])
      pop_out$R <- rbind(pop_out$R, pop$R[5,])
      pop_out$Incidence <- rbind(pop_out$Incidence, pop$Incidence[5,])
      pop_out$DetIncidence <- rbind(pop_out$DetIncidence, pop$DetIncidence[5,])
    }

    # cohort aging

    # extract the final state from pop
    S <- as.vector(t(data.table::last(pop$S)))
    E <- as.vector(t(data.table::last(pop$E)))
    I <- as.vector(t(data.table::last(pop$I)))
    R <- as.vector(t(data.table::last(pop$R)))
    Incidence <- as.vector(t(data.table::last(pop$Incidence)))
    DetIncidence <- as.vector(t(data.table::last(pop$DetIncidence)))

    # initialise the new initial condition vectors
    I0 <- rel_sizes*0
    S0 <- rel_sizes*0
    E0 <- rel_sizes*0
    R0 <- rel_sizes*0
    Incidence <- rel_sizes*0
    DetIncidence <- rel_sizes*0

    # then fill them in
    S0[1] <- perth_pop * rel_sizes[1]

    for(i in c(2:nAges)){
      S0[i] = S[(i - 1)] * trans_rate[(i-1)] + S[i] - S[i] * trans_rate[i]
      E0[i] = E[(i - 1)] * trans_rate[(i-1)] + E[i] - E[i] * trans_rate[i]
      I0[i] = I[(i - 1)] * trans_rate[(i-1)] + I[i] - I[i] * trans_rate[i]
      R0[i] = R[(i - 1)] * trans_rate[(i-1)] + R[i] - R[i] * trans_rate[i]
      Incidence0[i] = Incidence[(i - 1)] * trans_rate[(i-1)] + Incidence[i] - Incidence[i] * trans_rate[i]
      DetIncidence0[i] = DetIncidence[(i - 1)] * trans_rate[(i-1)] + DetIncidence[i] - DetIncidence[i] * trans_rate[i]
    }

    pars <- list(
      b0 = b0,
      b1 = b1,
      phi = phi,
      delta = delta,
      gamma = gamma,
      nu = nu,
      prop_detected_vect = prop_detected_vect,
      sigma_vect = sigma_vect,
      omega_vect = omega_vect,
      mixing = mixing,
      S0 = S0,
      E0 = E0,
      I0 = I0,
      R0 = R0,
      Incidence0 = Incidence0,
      DetIncidence0 = DetIncidence0
    )

    mod <- x(user = pars)
    T0 <- T0 + 1

  }
  pop_out <- pop_out[2:8]
  return(pop_out)
}

