#' Get parameters for the RSVsim_log_likelihood function from the data
#'
#' Helper function to assess the times and ages in the data, and return a list of suitable fixed parameters and
#' the model control option values given as arguments to the function.
#' @param data Data frame of incidence data: "incidence", "time" and "age_chr" columns must be included.
#' @inheritParams RSVsim_parameters
#' @inheritParams RSVsim_contact_matrix
#' @inheritParams RSVsim_run_model_dust
#' @return List of fixed parameters, max_t and warm_up.
#' @export
RSVsim_calibration_parameters <- function(fitted_parameter_names = c("b0", "b1", "phi"),
                                          data,
                                          overrides = list(),
                                          country = "United Kingdom",
                                          cohort_step_size = 10,
                                          warm_up = NULL){

  if(!all(is.integer(data[,"incidence"]))){
    stop("Not all values in incidence_data are integers")
    }

  if(!all(is.numeric(data[,"time"]))){
    stop("Not all values in time_data are numeric")
    }

  if(!all(is.character(data[,"age_chr"]))){
    stop("Not all values in age_chr_data are characters")
    }

   max_t <- max(data$time)

   times <- sort(unique(data$time))

   if(!is.null(warm_up)){
     max_t <- max_t + warm_up
     times <- times + warm_up
   }

   times <- sort(unique(c(times, seq(0, max_t, cohort_step_size / 20))))

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

   out <- list("fixed_parameters" = fixed_parameters,
               "fitted_parameter_names" = fitted_parameter_names,
               "times" = times,
               "cohort_step_size" = cohort_step_size,
               "warm_up" = warm_up)

   return(out)
}

#' Calculate the likelihood of observing incidence data given the model outputs
#'
#' Function to run the model and calculate the likelihood. A poisson likelihood is used for incidence.
#' The parameters must be constrained, such that b0 is greater than 0, b1 is between 0 and 1 and phi is between 0 and 365.25.
#' The incidence data must be provided as a count.
#' @param fitted_parameters index 1: b0 - Transmission rate coefficient parameter, index 2: b1 - amplitude of seasonal forcing parameter and
#'  index 3: phi - phase shift of seasonality parameter.
#' @inheritParams RSVsim_parameters
#' @inheritParams RSVsim_run_model_dust
#' @return Log-likelihood assuming the incidence data is Poisson distributed.
#' @export
RSVsim_log_likelihood_dust <- function(fitted_parameters,
                                       fitted_parameter_names,
                                       fixed_parameters,
                                       data,
                                       times,
                                       cohort_step_size = 10,
                                       init_conds = NULL,
                                       warm_up = NULL
                                  ){

  for(i in 1:length(fitted_parameters)){
    fixed_parameters[[fitted_parameter_names[i]]] <- fitted_parameters[i]
  }

  parameters <- fixed_parameters

  out <- RSVsim_run_model_dust(parameters = parameters,
                               times = times,
                               warm_up = warm_up,
                               cohort_step_size = cohort_step_size,
                               init_conds = init_conds
                               )

  data <- dplyr::left_join(data, out[,c("time", "Incidence", "age_chr")], by = c("time", "age_chr")) |>
    dplyr::mutate(log_likelihood = dpois(x = incidence, lambda = Incidence, log = TRUE))

  log_likelihood <- sum(data$log_likelihood, na.rm = FALSE)

  if (!is.finite(log_likelihood)){
    warning("The log-likelihood is not finite")
    return(-Inf)
    } else{
      return(log_likelihood)
      }
}

#' Constrained maximum-likelihood estimation of the b0, b1 and phi parameters from incidence data.
#'
#' Wrapper function for the \code{stats::nlminb function} (non-linear box-constrained quasi-Newton optimisation).
#' Used to maximise the \code{RSVsim_log_likelihood_dust} function with data and parameters formatted by the \code{RSVsim_calibration_parameters} function.
#' @inheritParams RSVsim_log_likelihood_dust
#' @param scale_parameters Boolean operator indicating whether to scale the parameters between 0 and 1 to help fitting multiple parameters on difference scales. Default: \code{FALSE}.
#' @param lower_ll Vector of lower and boundaries on the fitted parameters.
#' @param upper_ll Vector of upper boundaries on the fitted parameters.
#' @return List of fitted parameter values.
#' @export
RSVsim_max_likelihood_dust <- function(data,
                                       fitted_parameter_names,
                                       fixed_parameters,
                                       times,
                                       cohort_step_size = 10,
                                       init_conds = NULL,
                                       warm_up = NULL,
                                       scale_parameters = FALSE,
                                       lower_ll,
                                       upper_ll){

  np <- length(fitted_parameter_names)

  if(length(lower_ll) != length(upper_ll) |
     length(upper_ll) != np){
    stop("The number of fitted parameters does not equal the lower or upper limit lengths")
  }

  # function to provide the parameters as a single vector and return the negative log-likelihood
  RSVsim_negative_log_likelihood <- function(fitted_parameters,
                                             fitted_parameter_names,
                                             fixed_parameters,
                                             data,
                                             times,
                                             cohort_step_size,
                                             init_conds,
                                             warm_up,
                                             scale_parameters,
                                             lower_ll,
                                             upper_ll){

    if(scale_parameters){
      fitted_parameters <- fitted_parameters * (upper_ll - lower_ll) + lower_ll
      }

    return(
      -RSVsim_log_likelihood_dust(fitted_parameters = fitted_parameters, # must be negative because nlminb minimises the output
                                  fitted_parameter_names = fitted_parameter_names,
                                  fixed_parameters = fixed_parameters,
                                  data = data,
                                  times = times,
                                  cohort_step_size = cohort_step_size,
                                  init_conds = init_conds,
                                  warm_up = warm_up
                                  )
      )
    }

  if(scale_parameters){
    lower_optim <- rep(0, np)
    upper_optim <- rep(1, np)
  } else{
    lower_optim <- lower_ll
    upper_optim <- upper_ll
  }

  start_optim <- lower_optim

  out <- stats::nlminb(start = start_optim,
                       objective = RSVsim_negative_log_likelihood,
                       lower = lower_optim,
                       upper = upper_optim,
                       control = list("trace" = 0, "rel.tol" = 1e-10),
                       fitted_parameter_names = fitted_parameter_names,
                       fixed_parameters = fixed_parameters,
                       data = data,
                       times = times,
                       cohort_step_size = cohort_step_size,
                       init_conds = init_conds,
                       warm_up = warm_up,
                       scale_parameters = scale_parameters,
                       lower_ll = lower_ll,
                       upper_ll = upper_ll)

  if(scale_parameters){
    out$par <- out$par * (upper_ll - lower_ll) + lower_ll
  }

  return(out)
}
