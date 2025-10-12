#' Summary function to calculate the age-specific timing of peak incidence.
#'
#' Function calculates the time when incidence peaks stratified by age over the whole simulation (excluding warm up).
#'
#' @param out \code{RSVsim_run_model} function output.
#' @return Summary statistics (peak incidence time ordered by age-group) from the simulation output (dataframe).
#' @export
RSVsim_peak <- function(out){

  peak_times <- out |> dplyr::group_by(age) |> dplyr::filter(Incidence == max(Incidence)) |>
    dplyr::slice_head(n = 1) |> dplyr::select(age, time) |> dplyr::arrange(age) |> dplyr::pull(time)

  return(peak_times)

}

#' Summary function to calculate the age-specific total incidence over the simulation.
#'
#' Function calculates the total incidence stratified by age over the whole simulation (excluding warm up).
#' @param out \code{RSVsim_run_model} function output.
#' @return Summary statistics (total incidence ordered by age-group) from the simulation output (dataframe).
#' @export
RSVsim_total_incidence <- function(out){

  total_incidence <- out |> dplyr::group_by(age) |> dplyr::summarise(total_incidence = sum(Incidence)) |>
    dplyr::arrange(age) |> dplyr::pull(total_incidence)

  return(total_incidence)

}

#' Summary function to calculate the age-specific amplitude of changes in incidence over the simulation.
#'
#' Function calculates the difference between maximum and minimum incidence stratified by age over the whole simulation (excluding warm up).
#'
#' @param out \code{RSVsim_run_model} function output.
#' @return Summary statistics (amplitude by age-group) from the simulation output (dataframe).
#' @export
RSVsim_amplitude <- function(out){

  amplitude <- out |> dplyr::group_by(age) |> dplyr::summarise(amplitude = (max(Incidence) - min(Incidence))) |>
    dplyr::arrange(age) |> dplyr::pull(amplitude)

  return(amplitude)

}

#' Distance function to calculate the absolute difference between two vectors.
#'
#' Calculates the absolute distance. Vectorised.
#'
#' @param target First vector of metrics.
#' @param target_star Second vector of metrics.
#' @return Absolute difference.
#' @export
RSVsim_abs_dist_fun <- function(target, target_star){
  return(abs(target - target_star))
}

#' Distance function to calculate the shortest periodic difference between two vectors.
#'
#' Calculates the absolute distance between two times when the times are circular. Vectorised.
#'
#' @param target First vector of metrics.
#' @param target_star Second vector of metrics.
#' @param period Cosine wave period.
#' @return Absolute difference.
#' @export
RSVsim_shortest_periodic_dist_fun <- function(target, target_star, period){
  return(abs(pmin(target - target_star, period - (target - target_star))))
}

#' Function to run a Approximate Bayesian Computation (ABC) rejection algorithm
#'
#' Runs the ABC-rejection algorithm.
#'
#' @param target Values to fit to.
#' @param epsilon Acceptable error for each target.
#' @param summary_fun Function to calculate the summary statistics equivalent to the target values.
#' @param dist_fun Function to the calculate the error between the target and \code{summary_fun} outputs.
#' @param prior_fun Function to sample from the priors for all parameters. Must return a vector.
#' @param n_prior_attempts Number of random samples from the prior to attempt for each accepted particle.
#' @param nparticles Integer. Number of samples from the approximate posterior.
#' @param used_seeds_all Vector. Seeds used when generating the prior samples for each accepted particle.
#' @param ncores Number of cores. If greater than one then it is run in parallel.
#' @param fitted_parameter_names Vector of names of the parameters that are being estimated.
#' @param fixed_parameter_list List of parameter values to run the model excluding the fitted parameters.
#' @inheritParams RSVsim_run_model
#' @return List of fixed parameters, max_t and warm_up.
#' @export
RSVsim_ABC_rejection <- function(target,
                                 epsilon,
                                 summary_fun,
                                 dist_fun,
                                 prior_fun,
                                 n_prior_attempts,
                                 nparticles,
                                 used_seeds_all,
                                 ncores=1,
                                 fitted_parameter_names,
                                 fixed_parameter_list,
                                 times = seq(0, 365*5, 0.25), # maximum time to run the model for
                                 cohort_step_size = 0.2*365, # time at which to age people
                                 warm_up = 365 * 4){

  # checking the inputs
  if(!all(is.numeric(target))){
    stop("target must be numeric")
  }

  if(!all(is.numeric(epsilon))){
    stop("epsilon must be numeric")
  }

  if(!(length(target) == length(epsilon) | length(epsilon) == 1)){
    stop("if a specific tolerance is used for each summary output length(epsilon) must be equal to length(target) else length(epsilon) must equal 1")
  }

  if(any(epsilon <= 0)){
    stop("increase epsilon above zero")
  }

  if(length(used_seeds_all)!=nparticles){
    stop('the length used_seeds_all must be nparticles')
  }

  nparams <- length(fitted_parameter_names)

  nAges <- fixed_parameter_list$nAges

  while_fun <- function(particle){

    i <- 1 # initialise the number of accepted particles
    j <- 1 # initialise the number of proposed particles

    used_seed <- used_seeds_all[particle]
    set.seed(used_seed)
    fitted_parameters_all <- prior_fun(n_prior_attempts)

    while(i <= 1){

      if(j > n_prior_attempts){
        stop("No particle accepted: increase n_prior_attempts or the tolerance (epsilon)")
      }

      fitted_parameters <- fitted_parameters_all[j,]

      parameters_in <- c(setNames(unlist(fitted_parameters, recursive = FALSE), fitted_parameter_names),
                         fixed_parameter_list)

      out <- RSVsim_run_model(parameters = parameters_in,
                              times = times,
                              cohort_step_size = cohort_step_size,
                              warm_up = warm_up)

      summary_stats <- summary_fun(out)

      distance <- dist_fun(target, summary_stats)

      if(all(distance <= epsilon)){

        acc_params <- c(unlist(fitted_parameters, recursive = FALSE), j, particle, used_seed)
        names(acc_params) <- c(fitted_parameter_names, "attempts", "particle_number", "prior_function_seed")

        i <- i + 1
      }

      j <- j + 1
    }

    return(acc_params)
  }


  if(ncores > 1){
    cl <- parallel::makePSOCKcluster(ncores) #not to overload your computer
    parallel::clusterExport(cl, varlist = c("RSV_ODE",
                                            "while_fun",
                                            "RSVsim_run_model",
                                            "RSVsim_total_incidence",
                                            "RSVsim_amplitude",
                                            "RSVsim_peak",
                                            "RSVsim_abs_dist_fun",
                                            "RSVsim_shortest_periodic_dist_fun",
                                            "fixed_parameter_list",
                                            "fitted_parameter_names",
                                            "nAges",
                                            "times",
                                            "cohort_step_size",
                                            "warm_up",
                                            "target",
                                            "epsilon",
                                            "prior_fun",
                                            "n_prior_attempts",
                                            "summary_fun",
                                            "dist_fun",
                                            "nparams"),
                            envir = environment()
                            )

    parallel::clusterEvalQ(cl, library(RSVsim))

  } else{
    cl <- NULL
  }

  pbapply::pboptions(type = "timer", char = "-")
  res <- pbapply::pblapply(cl = cl, X = 1:nparticles, FUN = while_fun)

  res <- do.call(rbind, res)

  if(ncores > 1){

    #stop cluster
    parallel::stopCluster(cl)

  }

  return(res)
}

#' #' Function to run a Approximate Bayesian Computation Sequential Monte Carlo (ABC-SMC) algorithm
#' #'
#' #' ----- NOT FINISHED ------ implementation of an ABC-SMC algorithm.
#' #'
#' #' @param target Values to fit to.
#' #' @param epsilon_matrix Matrix of acceptable error for each target (columns) for each generation (rows). The row number of is used to determine the number of generations.
#' #' @param summary_fun Function to calculate the summary statistics equivalent to the target values. Should take one argument: the model outputs.
#' #' @param dist_fun Function to the calculate the error between the target and \code{summary_fun} outputs. Should take one argument: the model outputs.
#' #' @param prior_fun Function to sample from the priors for all parameters. Must return a vector.
#' #' @param prior_dens_fun Function that calculates the joint probability density of all parameters given the prior distributions.
#' #' @param nparticles Number of samples from the approximate posterior for each generation.
#' #' @param ncores Number of cores. If greater than one then it is run in parallel.
#' #' @param fitted_parameter_names Vector of names of the parameters that are being estimated.
#' #' @param fixed_parameter_list List of parameter values to run the model excluding the fitted parameters.
#' #' @inheritParams RSVsim_run_model
#' #' @return List of fixed parameters, max_t and warm_up.
#' #' @export
#' RSVsim_ABC_SMC <- function(target,
#'                            epsilon_matrix,
#'                            summary_fun,
#'                            dist_fun,
#'                            prior_fun,
#'                            prior_dens_fun,
#'                            particle_low,
#'                            particle_up,
#'                            nparticles,
#'                            ncores=1,
#'                            fitted_parameter_names,
#'                            fixed_parameter_list,
#'                            times = seq(0, 365*5, 0.25), # maximum time to run the model for
#'                            cohort_step_size = 0.2*365, # time at which to age people
#'                            warm_up = 365 * 4){
#'
#'   nparams <- length(fitted_parameter_names)
#'
#'   ntargets <- length(target)
#'
#'   nAges <- fixed_parameter_list$nAges
#'
#'   # checking the inputs
#'   if(!all(is.numeric(target))){
#'     stop("target must be numeric")
#'   }
#'
#'   if(!all(is.numeric(epsilon_matrix))){
#'     stop("epsilon_matrix must be numeric")
#'   }
#'
#'   if(any(epsilon_matrix <= 0)){
#'     stop("increase all epsilon_matrix values above zero")
#'   }
#'
#'   if(ncol(epsilon_matrix)!=nparams){
#'     stop("incorrect number of columns in epsilon_matrix")
#'   }
#'
#'   if(length(particle_up) != length(particle_low) | length(particle_low) != nparams){
#'     stop("incorrect number of lower or upper particle boundaries")
#'   }
#'
#'   # number of generations
#'   G <- nrow(epsilon_matrix)
#'
#'   # Number of simulations for each parameter set
#'   n <- 1
#'
#'   # Empty matrices to store results (5 model parameters)
#'   res_old <- matrix(ncol = nparams, nrow = nparticles)
#'   res_new <- matrix(ncol = nparams, nrow = nparticles)
#'
#'   # Empty vectors to store weights
#'   w_old <- matrix(ncol = 1, nrow = nparticles)
#'   w_new <- matrix(ncol = 1, nrow = nparticles)
#'
#'   for(g in 1:G){
#'
#'     while(i <= nparticles){
#'
#'       ### selecting parameters (particles)
#'       if(g == 1){
#'
#'         fitted_parameters <- prior_fun(1)
#'
#'       } else{
#'
#'         # sample particle set from previously fitted parameters
#'         p <- sample(seq(1, nparticles), 1, prob = w_old)
#'
#'         # perturb the particle to obtain theta**
#'         fitted_parameters <- tmvtnorm::rtmvnorm(1, mean = res_old[p,], sigma = sigma, lower = particle_low, upper = particle_up)
#'
#'       }
#'
#'       parameters <- c(as.list(setNames(fitted_parameters, fitted_parameter_names)),
#'                       fixed_parameter_list)
#'
#'       p_non_zero <- as.numeric(prior_dens_fun(fitted_parameters) > 0)
#'
#'       if(p_non_zero){
#'         m <- 0
#'         distance <- matrix(nrow = n, ncol = ntargets)
#'
#'         for(j in 1:n){
#'
#'           out <- RSVsim_run_model(parameters = parameters,
#'                               times = times,
#'                               cohort_step_size = cohort_step_size,
#'                               warm_up = warm_up)
#'
#'           target_star <- summary_fun(out)
#'
#'           distance[j,] <- dist_fun(target, target_star)
#'
#'
#'         }
#'
#'
#'
#'       }
#'
#'
#'     }
#'
#'
#'   }
#'
#'
#'   if(ncores > 1){
#'
#'     chunk_size <- ceiling(nparticles / ncores)
#'
#'     start_rows <- seq(1, nparticles, by = chunk_size)
#'     end_rows <- c(start_rows[-1] - 1, nparticles)
#'
#'     res_chunks <-
#'       lapply(1:ncores, function(i){
#'         as.matrix(res[start_rows[i]:end_rows[i], ])
#'       }
#'       )
#'
#'     cl <- parallel::makeCluster(ncores) #not to overload your computer
#'     parallel::clusterExport(cl, varlist = c("while_fun",
#'                                             "RSVsim_run_model",
#'                                             "RSVsim_total_incidence",
#'                                             "RSVsim_amplitude",
#'                                             "RSVsim_peak",
#'                                             "RSVsim_abs_dist_fun",
#'                                             "RSVsim_shortest_periodic_dist_fun",
#'                                             "fixed_parameter_list",
#'                                             "fitted_parameter_names",
#'                                             "nAges",
#'                                             "times", "cohort_step_size",
#'                                             "warm_up",
#'                                             "target",
#'                                             "epsilon",
#'                                             "prior_fun",
#'                                             "summary_fun",
#'                                             "dist_fun",
#'                                             "res_chunks",
#'                                             "nparams")
#'     )
#'
#'     res_out <- parallel::parLapply(cl = cl, X = res_chunks, fun = while_fun)
#'
#'     #stop cluster
#'     parallel::stopCluster(cl)
#'
#'     res_out <- do.call(rbind, res_out)
#'
#'   } else{
#'     res_out <- while_fun(res_in = res)
#'   }
#'
#'   return()
#' }
#'

