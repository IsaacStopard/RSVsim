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

  # initial conditions: choose whether to reset here, or read in from csv file
  if(is.null(init_conds)) {
    parameters <- purrr::list_modify(
      parameters,
        Sp0 = rel_sizes * parameters$population * 0.99,
        Ep0 = 0,
        Ip0 = rel_sizes * parameters$population * 0.01,
        Ss0 = 0,
        Es0 = 0,
        Is0 = 0,
        R0 = 0,
        DetIncidence0 = 0,
        Incidence0 = 0
      )
  } else{
    parameters <- purrr::list_modify(
      parameters,
      Sp0 = init_conds$Sp0,
      Ep0 = init_conds$Ep0,
      Ip0 = init_conds$Ip0,
      Ss0 = init_conds$Ss0,
      Es0 = init_conds$Es0,
      Is0 = init_conds$Is0,
      R0 = init_conds$R0,
      DetIncidence0 = init_conds$DetIncidence0,
      Incidence0 = init_conds$Incidence0
    )
  }

  # times
  T0 <- 0

  # runs the model with cohort ageing
  mod <- RSV_ODE$new(pars = parameters,
                     time = T0,
                     seed = 123,
                     deterministic = TRUE)


  # aging related parameters
  size_cohorts <- c(diff(parameters$age.limits), parameters$max_age - parameters$age.limits[length(parameters$age.limits)])

  transition_rate <- 1/size_cohorts
  rel_sizes <- size_cohorts/sum(size_cohorts)

  pop_out <-

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

    # cohort ageing

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

