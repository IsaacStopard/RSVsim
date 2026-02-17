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

#' Helper function to update the parameter list by names
#'
#' Changes the values in the list of parameters to run the model.
#' To fit parameters stored within vectors or matrices, \code{fitted_parameter_names} can also include indexing such as "parameter_name\[i\]".
#'
#' @param fixed_parameter_list List of parameters, should be the output of the \code{RSVsim_parameters} function.
#' @param fitted_parameter_names Vector of parameter names to be updated.
#' @param fitted_parameter_values Vector of updated parameter values.
#' @return Parameter list.
#' @export
RSVsim_update_parameters <- function(fixed_parameter_list, fitted_parameter_names, fitted_parameter_values){

  updated_parameters <- fixed_parameter_list

  update_string <- paste(
    paste0("updated_parameters$", fitted_parameter_names, " <- ", fitted_parameter_values),
    collapse = "\n")

  # 3. Evaluate once
  eval(parse(text = update_string))

  return(updated_parameters)
}

#' Function to run an Approximate Bayesian Computation (ABC) rejection algorithm
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
#' @param fitted_parameter_names Vector of names of the parameters that are being estimated. To fit vectors or matrices the individual element must be identified in the fitted_parameter_names.
#' @param fixed_parameter_list List of the original parameter values to run the model.
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
                                 times = seq(0, 365.25*5, 0.25), # maximum time to run the model for
                                 cohort_step_size = 1/12*365.25, # time at which to age people
                                 warm_up = 365.25 * 4){

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

      parameters_in <- RSVsim_update_parameters(fixed_parameter_list, fitted_parameter_names, fitted_parameters)

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
                                            "RSVsim_update_parameters",
                                            "RSVsim_run_model",
                                            "RSVsim_total_incidence",
                                            "RSVsim_amplitude",
                                            "RSVsim_peak",
                                            "RSVsim_abs_dist_fun",
                                            "RSVsim_shortest_periodic_dist_fun",
                                            "RSVsim_update_parameters",
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
                                            "nparams",
                                            "used_seeds_all"),
                            envir = environment()
                            )

    parallel::clusterEvalQ(cl, {
      library(RSVsim)
      library(tidyr)
      library(dplyr)
      }
      )

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

#' Function to run an Approximate Bayesian Computation Sequential Monte Carlo (ABC-SMC) algorithm
#'
#' Runs the ABC-SMC algorithm.
#'
#' @param epsilon_matrix Matrix of tolerance values. Different columns correspond to the values for different data points and different rows correspond to the values for the different generations.
#' @param n_param_attempts_per_accept Number of samples to try for each accepted particle.
#' @param used_seed_matrix Matrix of seeds: number of rows must be equal to the number of generations and number of columns must be equal to nparticles.
#' @param prior_dens_fun Function that calculates the probability density of all parameters given the prior distributions. The joint probability is the product of the values returned by this function.
#' @param particle_low Lower bounds on the parameters.
#' @param particle_up Upper bounds on the parameters.
#' @inheritParams RSVsim_run_model
#' @inheritParams RSVsim_ABC_rejection
#' @return List of fixed parameters, max_t and warm_up.
#' @export
RSVsim_ABC_SMC <- function(target,
                           epsilon_matrix,
                           summary_fun,
                           dist_fun,
                           prior_fun,
                           n_param_attempts_per_accept,
                           used_seed_matrix,
                           prior_dens_fun,
                           particle_low,
                           particle_up,
                           nparticles,
                           ncores=1,
                           fitted_parameter_names,
                           fixed_parameter_list,
                           times = seq(0, 365.25*5, 0.25), # maximum time to run the model for
                           cohort_step_size = 1/12*365.25, # time at which to age people
                           warm_up = 365.25 * 4){

  nparams <- length(fitted_parameter_names)

  ntargets <- length(target)

  nAges <- fixed_parameter_list$nAges

  # checking the inputs
  if(!all(is.numeric(target))){
    stop("target must be numeric")
  }

  if(!all(is.numeric(epsilon_matrix))){
    stop("epsilon_matrix must be numeric")
  }

  if(any(epsilon_matrix <= 0)){
    stop("increase all epsilon_matrix values above zero")
  }

  if(ncol(epsilon_matrix)!=ntargets){
    stop("incorrect number of columns in epsilon_matrix")
  }

  if(length(particle_up) != length(particle_low) | length(particle_low) != nparams){
    stop("incorrect number of lower or upper particle boundaries")
  }

  # number of generations
  G <- nrow(epsilon_matrix)

  # Number of simulations for each parameter set
  n <- 1

  while_fun_SMC <- function(particle, g, w_old, res_old, nparticles, sigma, n_param_attempts_per_accept, n, target, epsilon_matrix, fitted_parameter_names, fixed_parameter_list,
                            times, cohort_step_size, warm_up, prior_fun, prior_dens_fun, particle_low, particle_up, used_seed_matrix, ntargets){

    used_seed <- used_seed_matrix[g, particle]
    set.seed(used_seed)

    ### selecting parameters (particles)
      if(g == 1){

        fitted_parameters <- prior_fun(n_param_attempts_per_accept)

      } else{

        p <- sample(seq(1, nparticles), n_param_attempts_per_accept, prob = w_old, replace = TRUE)

        # sample particle set from previously fitted parameters
        fitted_parameters <- as.data.frame(
          t(
            sapply(1:n_param_attempts_per_accept, function(a){
              return(
                tmvtnorm::rtmvnorm(1,
                                   mean = unlist(as.vector(unname(res_old[p[a],, drop = FALSE]))),
                                   sigma = sigma, lower = particle_low, upper = particle_up))
              },
              simplify = TRUE)
            )
          )
      }

    i <- 1  # initialise the number of accepted particles
    k <- 1 # initialise the number of attempted particles

    while(i <= 1){

      parameters <- fitted_parameters[k, ]

      parameters_ODE <- RSVsim_update_parameters(fixed_parameter_list, fitted_parameter_names, parameters)

      p_non_zero <- as.numeric(prod(prior_dens_fun(parameters)) > 0)

      if(p_non_zero){
        m <- 0
        distance <- matrix(nrow = n, ncol = ntargets)

        for(j in 1:n){

          out <- RSVsim_run_model(parameters = parameters_ODE,
                                  times = times,
                                  cohort_step_size = cohort_step_size,
                                  warm_up = warm_up)

          target_star <- summary_fun(out)

          distance[j,] <- dist_fun(target, target_star)

          if(all(distance[j,] <= epsilon_matrix[g,])){
            m <- m + 1
          }
        }

        if(m > 0){
          res_new_out <- parameters

          w1 <- prod(prior_dens_fun(parameters))

          if(g == 1){
            w2 <- 1
          } else{
            w2 <- sum(sapply(1:nparticles, function(a){
              w_old[a]*tmvtnorm::dtmvnorm(as.vector(unname(unlist(res_new_out))), mean = as.vector(unname(unlist(res_old[a,]))), sigma = sigma, lower = particle_low, upper = particle_up)
            }))
          }
          w_new_out <- (m/n)*w1/w2
          i <- i + 1
        }
      }

      k <- k + 1

    }

    return(list("res_new" = res_new_out,
                "w_new" = w_new_out,
                "seed" = used_seed))
  }


  w_list <- vector(mode = "list", length = G)
  res_list <- vector(mode = "list", length = G)
  sigma_list <- vector(mode = "list", length = G)

  # Empty matrices to store results (5 model parameters)
  res_old <- matrix(nrow = nparticles, ncol = nparams)
  res_new <- matrix(nrow = nparticles, ncol = nparams)

  # Empty vectors to store weights
  w_old <- matrix(nrow = nparticles, ncol = 1)
  w_new <- matrix(nrow = nparticles, ncol = 1)

  sigma <- matrix(nrow = nparams, ncol = nparams)

  for(g in 1:G){

    print(paste("Generation", g, "of", G))

    if(ncores > 1){
      cl <- parallel::makePSOCKcluster(ncores) #not to overload your computer
      parallel::clusterExport(cl, varlist = c("RSV_ODE",
                                              "while_fun_SMC",
                                              "RSVsim_run_model",
                                              "RSVsim_update_parameters",
                                              "fixed_parameter_list",
                                              "fitted_parameter_names",
                                              "nAges",
                                              "times",
                                              "cohort_step_size",
                                              "warm_up",
                                              "target",
                                              "epsilon_matrix",
                                              "prior_fun",
                                              "prior_dens_fun",
                                              "n_param_attempts_per_accept",
                                              "summary_fun",
                                              "dist_fun",
                                              "nparams",
                                              "used_seed_matrix",
                                              "particle_low",
                                              "particle_up",
                                              "ntargets"),
                              envir = environment()
      )

      parallel::clusterEvalQ(cl, {
        library(RSVsim)
        library(tidyr)
        library(dplyr)
        library(tmvtnorm)
      }
      )

    } else{
      cl <- NULL
    }

    pbapply::pboptions(type = "timer", char = "-")
    res <- pbapply::pblapply(cl = cl, X = 1:nparticles, FUN = while_fun_SMC,
                             g = g,
                             w_old = w_old,
                             res_old = res_old,
                             nparticles = nparticles,
                             sigma = sigma,
                             n_param_attempts_per_accept = n_param_attempts_per_accept,
                             n = n,
                             target = target,
                             epsilon_matrix = epsilon_matrix,
                             fitted_parameter_names = fitted_parameter_names,
                             fixed_parameter_list = fixed_parameter_list,
                             times = times,
                             cohort_step_size = cohort_step_size,
                             warm_up = warm_up,
                             prior_fun = prior_fun,
                             prior_dens_fun = prior_dens_fun,
                             particle_low = particle_low,
                             particle_up = particle_up,
                             used_seed_matrix = used_seed_matrix,
                             ntargets = ntargets)

    if(ncores > 1){

      #stop cluster
      parallel::stopCluster(cl)

    }

    w_new <- sapply(res, '[[', "w_new")
    res_new <- do.call(rbind, lapply(res, '[[', 'res_new'))

    sigma <- cov(res_new)
    res_old <- res_new
    w_old <- w_new / sum(w_new)

    w_list[[g]] <- w_old
    res_list[[g]] <- res_old
    sigma_list[[g]] <- sigma
  }

  return(list("fitted_parameters" = res_list,
              "weights" = w_list,
              "sigma" = sigma_list))
}


