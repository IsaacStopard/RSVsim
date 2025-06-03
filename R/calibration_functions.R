#' Get parameters for the calibration_likelihood function from the data
#'
#' Function to assess the times and ages in the data, and return suitable parameters to run the model.
#' @param data
#' @param warm_up Length of time-points to exclude before calculating the likelihood. Default: \code{NULL}.
#' @return List of parameters, max_t and warm_up.
#' @export

get_calibration_variables <- function(data,
                                      warm_up = NULL){

  max_t <- max(sapply(data, function(df){max(df$time)}))

  if(!is.null(warm_up)){
    max_t <- max_t + warm_up
  }

  age_chr <- unique(sapply(data, function(df){unique(df$age_chr)}, simplify = FALSE)) |> unlist()


  lapply(age_chr, function(i){
    split <- strsplit(i, ",") |> unlist()
    number <- readr::parse_number(split)
    bracket <- grepl("(", split)
    }
    )
}


#' Calculate the likelihood of observing incidence or prevalence data given the model outputs
#'
#' Function to run the model and calculate the likelihood. A poisson likelihood is used for incidence.
#' The incidence data must be provided as a count in the
#' @param b0 Vector of parameter names that we want to fit.
#' @param fixed_parameters List of parameters that are fixed.
#' @inheritParams data
#' @inheritDotParams run_model max_t cohort_step_size dt init_conds.
#' @return Log-likelihood.
#' @export
calibration_likelihood <- function(b0,
                                   fixed_parameters,
                                   data,
                                   ...
                                   ){

  param_names <- c("b0", "b1", "phi", "delta", "gamma_s", "gamma_p", "nu", "omega_vect", "prop_detected_vect", "sigma_vect", "alpha_vect", "nAges", "age.limits", "matrix_mean", "max_age", "total_population")

  parameters <- purrr::list_modify(fixed_parameters,
                                   b0 = b0)

  if(!all(param_names %in% names(parameters))){
    stop(paste("calibration_likelihood: not all required parameters are included", sep = " "))
  }

  out <- run_model(parameters = parameters, max_t = max_t, ...)

  if("prevalence" %in% names(data)){
    dplyr::left_join(data[["prevalence"]], out[,c("time", "prev", "age_chr")], by = c("time", "age_chr")) <-

  }

  if("incidence" %in% names(data)){

  }


  #if()

}

incidence <- data.frame(time = seq(0, 365)) |>
  dplyr::mutate(value = rpois(1, 1 + 0.5 * cos(2 * pi * time / 365.25)),
                age_chr = ifelse(time < 182, "[0, 0.2]", "(20, 40]"))

prevalence <- data.frame(time = seq(0, 365)) |>
  dplyr::mutate(sample = rpois(1, 100),
                positive = rbinom(1, size = sample, prob = 0.1 + 0.1 * cos(2 * pi * time / 365.25)),
                age_chr = ifelse(time < 182, "[0, 0.2]", "(20, 40]"))

data <- list("incidence" = incidence,
             "prevalence" = prevalence)
