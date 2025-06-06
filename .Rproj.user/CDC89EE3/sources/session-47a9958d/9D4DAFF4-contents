#' Get parameters for the calibration_likelihood function from the data
#'
#' Helper function to assess the times and ages in the data, and return suitable parameters to run the model.
#' @param data List of data frames of incidence data. For all data frames "incidence", "time" and "age_chr" columns must be included.
#' For incidence data frames a "value" column is needed to contain the incidence values.
#' For the prevalence data frame a "sample" column is required with the number of people tested and a "positive" column with the number of people that test positive.
#' @param data_populations Vector of total population sizes for each dataframe in \code{data}.
#' @inheritParams create_contact_matrix
#' @inheritParams run_model
#' @return List of fixed parameters, max_t and warm_up.
#' @export
get_calibration_parameters <- function(data,
                                       data_populations,
                                       country = "United Kingdom",
                                       warm_up = NULL){

  n_d <- length(names(data))

  out <- vector(mode = "list", length = n_d)

  for(name in names(data)){
    for(column_name in c("time", "age_chr", "incidence")){

      if(!column_name %in% colnames(data[[name]])){
        stop(paste("Cannot find the", column_name, "column in the", name, "data frame.", sep = " "))
        }

      if(column_name == "incidence"){
        if(!all(is.integer(data[[name]][[column_name]]))){
          stop(paste("Not all values in the", column_name , "column in the", name, "data frame are integers", sep = " "))
        }
        }

      if(column_name == "time"){
        if(!all(is.numeric(data[[name]][[column_name]]))){
          stop(paste("Not all values in the", column_name , "column in the", name, "data frame are numeric", sep = " "))
        }
        }

      if(column_name == "age_chr"){
        if(!all(is.character(data[[name]][[column_name]]))){
          stop(paste("Not all values in the", column_name, "column in the", name, "data frame are characters or factors", sep = " "))
        }
        }

    }
  }

  for(i in 1:n_d){

   df <- data[[i]]

   max_t <- max(df$time)
   if(!is.null(warm_up)){
     max_t <- max_t + warm_up
   }

   age_chr <- unique(df$age_chr)

   age.limits <- lapply(age_chr, function(i){
     split <- strsplit(i, ",") |> unlist()
     number <- readr::parse_number(split)}
     ) |> unlist() |> unique() |> sort()

   contact_population_list <- create_contact_matrix(country = country, age.limits = age.limits)

   fixed_parameters <- get_parameters(overrides = list("total_population" = data_populations[i]),
                                      contact_population_list = contact_population_list,
                                      fitted = c("b0", "b1", "phi"))

   out_t <- list("fixed_parameters" = fixed_parameters,
                 "max_t" = max_t,
                 "warm_up" = warm_up)

   out[[i]] <- out_t

  }

  return(out)
}

#' Calculate the likelihood of observing incidence data given the model outputs
#'
#' Function to run the model and calculate the likelihood. A poisson likelihood is used for incidence.
#' The parameters must be constrained, such that b0 is greater than 0, b1 is between 0 and 1 and phi is between 0 and 365.25.
#' The incidence data must be provided as a count.
#' @param fitted_parameters Vector of parameters that we want to fit. Names must include "b0", "b1" and "phi".
#' @param fixed_parameter_list List of parameters that are fixed.
#' @inheritParams get_calibration_parameters
#' @param minimise Boolean operator indicating whether the log-likelihood should be multiplied by -1.
#' This is useful is trying to estimate the maximum likelihood with an optimisation algorithm that estimates the minimum.
#' @param scale_parameters List of lower and upper boundaries on the fitted parameters.
#' Default: \code{NULL}. If \code{NULL} the parameters are not transformed.
#' If a list is provided the parameters log-likelihood is calculated assuming the \code{fitted_parameters}
#' have already been scaled between 0 and 1 with the given lower and upper limits.
#' @inheritDotParams run_model max_t cohort_step_size dt init_conds warm_up
#' @return Log-likelihood.
#' @export
calibration_likelihood <- function(fitted_parameters,
                                   fixed_parameter_list,
                                   data,
                                   minimise = FALSE,
                                   scale_parameters = NULL,
                                   ...
                                   ){

  if(!is.null(scale_parameters)){
    fitted_parameters <- fitted_parameters * (scale_parameters$upper - scale_parameters$lower) + scale_parameters$lower
  }

  b0 <- fitted_parameters[1]
  b1 <- fitted_parameters[2]
  phi <- fitted_parameters[3]

  if(length(data) != length(fixed_parameter_list)){
    stop("data list length does not equal fixed parameter list length")
  }

  n_d <- length(data)

  param_names <- c("b0", "b1", "phi", "delta", "gamma_s", "gamma_p", "nu",
                   "omega_vect", "prop_detected_vect", "sigma_vect", "alpha_vect", "nAges", "age.limits", "matrix_mean", "max_age", "total_population")

  log_likelihood <- c()

  for(i in 1:n_d){
    parameters <- purrr::list_modify(fixed_parameter_list[[i]]$fixed_parameters,
                                     b0 = b0,
                                     b1 = b1,
                                     phi = phi)

    if(!all(param_names %in% names(parameters))){
      stop(paste("calibration_likelihood: not all required parameters are included for dataset at position", i,"in fixed_parameter_list", sep = " "))
    }

    out <- run_model(parameters = parameters,
                     max_t = fixed_parameter_list[[i]]$max_t,
                     warm_up = fixed_parameter_list[[i]]$warm_up,
                     ...)

    data[[i]] <- dplyr::left_join(data[[i]], out[,c("time", "Incidence", "age_chr")], by = c("time", "age_chr")) |>
      dplyr::mutate(log_likelihood = dpois(x = incidence, lambda = Incidence, log = TRUE))

    if(any(is.na(data[[i]]$log_likelihood))){
      stop(paste("There were", sum(is.na(data[[i]]$log_likelihood)), "NAs in the log-likelihood for dataset", i, sep = " "))
    }

    log_likelihood <- c(log_likelihood, sum(data[[i]]$log_likelihood))

  }

  if(minimise){
    return(-sum(log_likelihood))
  } else{
    return(sum(log_likelihood))
  }

}

#' Constrained maximum-likelihood estimation of the b0, b1 and phi parameters from incidence data.
#'
#' Helper function to assess the times and ages in the data, and return suitable parameters to run the model. Initial parameter values are set to 0.
#' @inheritParams calibration_likelihood
#' @return List of fitted parameter values for b0, b1 and phi.
#' @export
constrained_max_likelihood <- function(fixed_parameter_list,
                                       data,
                                       scale_parameters,
                                       cohort_step_size,
                                       dt){

  if(!is.null(scale_parameters)){
    lower_optim <- rep(0, 3)
    upper_optim <- rep(1, 3)
    start_optim <- c(0.25, 0.25, 0.25)
  } else{
    lower_optim <- c(0.01, 0, 0)
    upper_optim <- c(10, 1, 365.25)
    start_optim <- c(0.1, 0.5, 100)
  }

  out <- nlminb(start = start_optim,
                objective = calibration_likelihood,
                lower = lower_optim,
                upper = upper_optim,
                control = list("trace" = 1, "rel.tol" = 1e-10),
                fixed_parameter_list = fixed_parameter_list,
                data = data,
                minimise = TRUE,
                scale_parameters = scale_parameters,
                cohort_step_size = cohort_step_size,
                dt = dt)

  if(!is.null(scale_parameters)){
    out$par <- out$par * (scale_parameters$upper - scale_parameters$lower) + scale_parameters$lower
  }

  return(out)
}
