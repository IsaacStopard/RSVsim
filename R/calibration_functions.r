#' Get fixed parameters for the RSVsim_calibrate function from the data
#'
#' Helper function to assess the ages in the data, and return a list of suitable fixed parameters and
#' the model control option values given as arguments to the function.
#'
#' @param data Dataframe of data to fit to: "age_chr" columns must be included.
#' @inheritParams RSVsim_parameters
#' @inheritParams RSVsim_contact_matrix
#' @inheritParams RSVsim_run_model
#' @return List of fixed parameters, max_t and warm_up.
#' @export
RSVsim_fixed_parameters <- function(fitted_parameter_names = c("b0"),
                                    data,
                                    overrides = list(),
                                    country = "United Kingdom"){


  if(!all(is.character(data[,"age_chr"]))){
    stop("Not all values in age_chr_data are characters")
  }

  age_chr <- unique(data$age_chr)

  age.limits <- lapply(age_chr, function(i){
    split <- strsplit(i, ",") |> unlist()
    number <- readr::parse_number(split)}) |>
    unlist() |> unique() |> sort()

  contact_population_list <- RSVsim_contact_matrix(country = country, age.limits = age.limits)

  fixed_parameters <- RSVsim_parameters(overrides = overrides,
                                        contact_population_list = contact_population_list,
                                        fitted_parameter_names = fitted_parameter_names
                                        )

  return(fixed_parameters)
}

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
  return(pmin(target - target_star, period - (target - target_star)))
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
#' @param nparticles Number of samples from the approximate posterior.
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
                                 nparticles,
                                 ncores=1,
                                 fitted_parameter_names,
                                 fixed_parameter_list,
                                 times = seq(0, 365*5, 0.25), # maximum time to run the model for
                                 cohort_step_size = 0.2*365, # time at which to age people\
                                 init_conds = NULL,
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

  nparams <- length(fitted_parameter_names)

  # res <- matrix(ncol = nparams + 1, nrow = nparticles) # Empty matrix to store results
  # colnames(res) <- c(fitted_parameter_names)

  nAges <- fixed_parameter_list$nAges

  #pb <- utils::txtProgressBar(min = i, max = n_particles_to_fit, style = 3)

  while_fun <- function(particle){

    i <- 1 # initialise the number of accepted particles
    j <- 1 # initialise the number of proposed particles

    while(i <= 1){

      used_seed <- MQMF::getseed()

      fitted_parameters <- prior_fun()

      parameters <- c(as.list(setNames(fitted_parameters, fitted_parameter_names)),
                      fixed_parameter_list)

      out <- RSVsim_run_model(parameters = parameters,
                              times = times,
                              cohort_step_size = cohort_step_size,
                              init_conds = init_conds,
                              warm_up = warm_up)

      summary_stats <- summary_fun(out)

      distance <- dist_fun(target, summary_stats)

      if(all(distance <= epsilon)){

        acc_params <- c(fitted_parameters, i/j, particle, used_seed)
        names(acc_params) <- c(fitted_parameter_names, "acc_rate", "particle_number", "seed")

        i <- i + 1
      }

      j <- j + 1
    }

    return(acc_params)
  }


  if(ncores > 1){

    # chunk_size <- ceiling(nparticles / ncores)
    #
    # start_rows <- seq(1, nparticles, by = chunk_size)
    # end_rows <- c(start_rows[-1] - 1, nparticles)
    #
    # res_chunks <-
    #   lapply(1:ncores, function(i){
    #     as.matrix(res[start_rows[i]:end_rows[i], ])
    #   }
    #   )

    cl <- parallel::makeCluster(ncores) #not to overload your computer
    parallel::clusterExport(cl, varlist = c("while_fun",
                                            "RSVsim_run_model",
                                            "RSVsim_total_incidence",
                                            "RSVsim_amplitude",
                                            "RSVsim_peak",
                                            "RSVsim_abs_dist_fun",
                                            "RSVsim_shortest_periodic_dist_fun",
                                            "fixed_parameter_list",
                                            "fitted_parameter_names",
                                            "nAges",
                                            "times", "cohort_step_size",
                                            "init_conds", "warm_up",
                                            "target",
                                            "epsilon",
                                            "prior_fun",
                                            "summary_fun",
                                            "dist_fun",
                                            "nparams"),
                            envir = environment()
                            )

  } else{
    cl <- NULL
  }

  pbapply::pboptions(type = "timer", char = "-")
  #pboptions(type = "tk")
  res <- pbapply::pblapply(cl = cl, X = 1:nparticles, FUN = while_fun)
  res <- do.call(rbind, res)

  if(ncores > 1){

    #stop cluster
    parallel::stopCluster(cl)

  }

  return(res)
}

#' Function to run a Approximate Bayesian Computation Sequential Monte Carlo (ABC-SMC) algorithm
#'
#' ----- NOT FINISHED ------ implementation of an ABC-SMC algorithm.
#'
#' @param target Values to fit to.
#' @param epsilon Acceptable error for each target - must be list of equal length to the generations.
#' @param summary_fun Function to calculate the summary statistics equivalent to the target values.
#' @param dist_fun Function to the calculate the error between the target and \code{summary_fun} outputs.
#' @param prior_fun Function to sample from the priors for all parameters. Must return a vector.
#' @param nparticles Number of samples from the approximate posterior for each generation.
#' @param nparticles Number of generations.
#' @param ncores Number of cores. If greater than one then it is run in parallel.
#' @param fitted_parameter_names Vector of names of the parameters that are being estimated.
#' @param fixed_parameter_list List of parameter values to run the model excluding the fitted parameters.
#' @inheritParams RSVsim_run_model
#' @return List of fixed parameters, max_t and warm_up.
#' @export
RSVsim_ABC_SMC <- function(target,
                           epsilon,
                           summary_fun,
                           dist_fun,
                           prior_fun,
                           nparticles,
                           ngenerations,
                           ncores=1,
                           fitted_parameter_names,
                           fixed_parameter_list,
                           times = seq(0, 365*5, 0.25), # maximum time to run the model for
                           cohort_step_size = 0.2*365, # time at which to age people\
                           init_conds = NULL,
                           warm_up = 365 * 4){

  nparams <- length(fitted_parameter_names)

  res <- matrix(ncol = nparams, nrow = nparticles) # Empty matrix to store results
  colnames(res) <- c(fitted_parameter_names)

  nAges <- fixed_parameter_list$nAges

  # checking the inputs
  if(!all(is.numeric(target))){
    stop("target must be numeric")
  }

  if(!all(is.numeric(epsilon))){
    stop("epsilon must be numeric")
  }


  while_fun <- function(res_in, cl){

    n_particles_to_fit <- nrow(res_in)

    i <- 1 # initialise the number of accepted particles
    j <- 1 # initialise the number of proposed particles

    pb <- utils::txtProgressBar(min = i, max = n_particles_to_fit, style = 3)

    while(i <= n_particles_to_fit){

      # look into the latin hypercube

      fitted_parameters <- prior_fun()

      parameters <- c(as.list(setNames(fitted_parameters, fitted_parameter_names)),
                      fixed_parameter_list)

      out <- RSVsim_run_model(parameters = parameters,
                              times = times,
                              cohort_step_size = cohort_step_size,
                              init_conds = init_conds,
                              warm_up = warm_up)

      summary_stats <- summary_fun(out)

      distance <- dist_fun(target, summary_stats)

      if(all(distance <= epsilon)){

        res_in[i, 1:nparams] <- c(fitted_parameters)


        i <- i + 1
      }

      j <- j + 1

      utils::setTxtProgressBar(pb, i)
    }

    close(pb)

    print(paste("acc_rate:", i/j*100, "%"))

    return(res_in)
  }


  if(ncores > 1){

    chunk_size <- ceiling(nparticles / ncores)

    start_rows <- seq(1, nparticles, by = chunk_size)
    end_rows <- c(start_rows[-1] - 1, nparticles)

    res_chunks <-
      lapply(1:ncores, function(i){
        as.matrix(res[start_rows[i]:end_rows[i], ])
      }
      )

    cl <- parallel::makeCluster(ncores) #not to overload your computer
    parallel::clusterExport(cl, varlist = c("while_fun",
                                            "RSVsim_run_model",
                                            "RSVsim_total_incidence",
                                            "RSVsim_amplitude",
                                            "RSVsim_peak",
                                            "RSVsim_abs_dist_fun",
                                            "RSVsim_shortest_periodic_dist_fun",
                                            "fixed_parameter_list",
                                            "fitted_parameter_names",
                                            "nAges",
                                            "times", "cohort_step_size",
                                            "init_conds", "warm_up",
                                            "target",
                                            "epsilon",
                                            "prior_fun",
                                            "summary_fun",
                                            "dist_fun",
                                            "res_chunks",
                                            "nparams")
    )

    res_out <- parallel::parLapply(cl = cl, X = res_chunks, fun = while_fun)

    #stop cluster
    parallel::stopCluster(cl)

    res_out <- do.call(rbind, res_out)

  } else{
    res_out <- while_fun(res_in = res)
  }

  return(res_out)
}


