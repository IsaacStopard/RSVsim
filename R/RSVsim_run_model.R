#' Runs the RSV model using Stan
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
RSVsim_run_model <- function(parameters,
                             times = seq(0, 365 * 4, 0.25),
                             cohort_step_size = 10,
                             init_conds = NULL,
                             warm_up = NULL
){

  ##########################################################

  # https://discourse.mc-stan.org/t/exposing-stan-user-defined-functions-using-cmdstanr-and-rcpp/27104/9
  #
  # m <- instantiate::stan_package_model(name = "RSV_ODE_stan", package = "RSVsim", compile_standalone = TRUE)
  #
  cmdstanr::cmdstan_model(stan_file = "src/stan/RSV_ODE_stan.stan", compile_standalone = TRUE)

  a <- cmdstanr::cmdstan_model(stan_file = "src/stan/RSV_ODE_stan.stan", compile_standalone = TRUE, compile = TRUE)

  a$functions$existing_exe <- FALSE

  # labels for the ages
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
  if(is.null(init_conds)){
    rep_z <- rep(0, parameters$nAges)
    Sp0 <- rel_sizes * parameters$total_population * 0.99
    Ep0 <- rep_z
    Ip0 <- rel_sizes * parameters$total_population * 0.01
    Ss0 <- rep_z
    Es0 <- rep_z
    Is0 <- rep_z
    R0 <- rep_z
    init_conds <- c(Sp0, Ep0, Ip0, Ss0, Es0, Is0, R0)
  } else{
    with(init_conds,{
      for(name in c("Sp0", "Ep0", "Ip0", "Ss0", "Es0", "Is0", "R0")){
        if(length(get(name)) != parameters$nAges){
          stop(paste("RSVsim_run_model: initial conditions for", name,"are not all the same length as the number of age categories", sep = " "))
        }
      }
    }
  )
    init_conds <- as.vector(c(init_conds[["Sp0"]], init_conds[["Ep0"]], init_conds[["Ip0"]], init_conds[["Ss0"]], init_conds[["Es0"]], init_conds[["Is0"]], init_conds[["R0"]]))

  }

  # running the model with cohort ageing (run for a single cohort, move cohort, change initial states, repeat)
  times_vec <- times
  n_times_vec <- length(times_vec)
  t0 <- times[1]
  times <- times[-1]

  n_times <- length(times)

  times_array <- lapply(1:n_steps, function(i, times, cohort_step_size){
    return(
      sort(unique(c(times[times > ((i - 1) * cohort_step_size) & times <= (i * cohort_step_size)])))
      )
    }, times = times, cohort_step_size = cohort_step_size
    )

  n_times_array <- as.vector(unlist(lapply(times_array, length)))
  cumn_times_array <- as.integer(cumsum(n_times_array))
  cumn_times_array <- c(0, cumn_times_array[-length(cumn_times_array)])

  Sp_index <- 0
  Ep_index <- 1
  Ip_index <- 2
  Ss_index <- 3
  Es_index <- 4
  Is_index <- 5
  R_index <- 6

  out <- cohort_ageing_stan(n_times = length(times),
                            n_steps = n_steps,
                            t0 = t0,
                            times_array = times_array,
                            n_times_array = n_times_array,
                            cumn_times_array = cumn_times_array,
                            Sp_index = Sp_index,
                            Ep_index = Ep_index,
                            Ip_index = Ip_index,
                            Ss_index = Ss_index,
                            Es_index = Es_index,
                            Is_index = Is_index,
                            R_index = R_index,
                            init_conds = init_conds,
                            nAges = parameters$nAges,
                            total_population = parameters$total_population,
                            b0 = parameters$b0,
                            b1 = parameters$b1,
                            phi = parameters$phi,
                            delta = parameters$delta,
                            gammas = parameters$gamma_s,
                            gammap = parameters$gamma_p,
                            nu = parameters$nu,
                            omega_vect = parameters$omega_vect,
                            sigma_vect = parameters$sigma_vect,
                            alpha_vect = parameters$alpha_vect,
                            matrix_mean = parameters$matrix_mean,
                            transition_rate = transition_rate,
                            rel_sizes = rel_sizes)

  inc <- calc_incidence_stan(out = out,
                             times_vec = times_vec,
                             n_times_vec = n_times_vec,
                             nAges = parameters$nAges,
                             Sp_index = Sp_index,
                             Ep_index = Ep_index,
                             Ip_index = Ip_index,
                             Ss_index = Ss_index,
                             Es_index = Es_index,
                             Is_index = Is_index,
                             R_index = R_index,
                             b0 = parameters$b0,
                             b1 = parameters$b1,
                             phi = parameters$phi,
                             omega_vect = parameters$omega_vect,
                             sigma_vect = parameters$sigma_vect,
                             alpha_vect = parameters$alpha_vect,
                             matrix_mean = parameters$matrix_mean)

  n_states <- 7
  states <- rep(NA, n_states)
  states[Sp_index + 1] <- "Sp"
  states[Ep_index + 1] <- "Ep"
  states[Ip_index + 1] <- "Ip"
  states[Ss_index + 1] <- "Ss"
  states[Es_index + 1] <- "Es"
  states[Is_index + 1] <- "Is"
  states[R_index + 1] <- "R"
  states <- rep(states, each = length(parameters$age.limits))

  age_vec <- rep(parameters$age.limits, n_states)
  age_chr_vec <- rep(age_chr, n_states)

  out <- lapply(1:length(out), function(i, out, times_vec, states, age_vec, age_chr_vec){
    return(data.frame(time = times_vec[i],
                      state = states,
                      value = out[[i]],
                      age = age_vec,
                      age_chr = age_chr_vec))
  }, out = out, times_vec = times_vec, states = states, age_vec = age_vec, age_chr_vec = age_chr_vec) |>
    dplyr::bind_rows() |> tidyr::pivot_wider(names_from = state, values_from = value) |> as.data.frame() |>
    dplyr::left_join(
      lapply(1:length(inc), function(i, inc, times_vec = times_vec, age, age_chr){
        return(data.frame(
          time = times_vec[i],
          incidence = inc[[i]],
          age = age,
          age_chr = age_chr)
          )
        }, inc = inc, age = parameters$age.limits, age_chr = age_chr, times_vec = times_vec) |> dplyr::bind_rows(),
      by = c("time", "age", "age_chr")
      )

  out[, "prev"] = (out[, "Ip"] + out[, "Is"])/(out[, "Sp"] + out[, "Ep"] + out[, "Ip"] + out[, "Ss"] + out[, "Es"] + out[, "Is"] + out[, "R"])

  if(!is.null(warm_up)){
    if(!is.numeric(warm_up)){
      stop("warm_up is not numeric")
    } else{
      out <- out |> dplyr::filter(time >= warm_up) |> dplyr::mutate(time = time - warm_up)
    }
  }

  return(out)
}

