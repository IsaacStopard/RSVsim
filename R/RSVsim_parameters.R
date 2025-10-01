#' Get parameters for RSVsim
#'
#' Provides the default parameters to run the RSV model.
#'
#' @param overrides List of default parameters to change.
#' @param contact_population_list List of outputs from the \code{create_contact_matrix} function.
#' @param fitted_parameter_names Vector of parameter names to exclude from the list because they are fitted. Default: \code{NULL}.
#' @return Parameter list. Default values: \cr
#' \code{b0}: 0.08. Transmission rate coefficient. \cr
#' \code{b1}: 0. Amplitude of seasonal forcing in the transmission rate. \cr
#' \code{phi}: 0. Phase shift of the seasonal forcing. \cr
#' \code{delta}: 1/4. Inverse of the latent period. \cr
#' \code{gamma_s}: 1/8. Inverse of the infectious period of subsequent (secondary) infections. \cr
#' \code{gamma_p}: 1/9. Inverse of the infectious period of first (primary) infections. \cr
#' \code{nu}: 1/200. Inverse of the duration of natural immunity. \cr
#' \code{omega_vect}: 0.35 if less than 5 years old else 1. Reduced infectiousness due to immunity in early life. \cr
#' \code{prop_detected_vect}: age-specific proportion detected. \cr
#' \code{sigma_vect}: 0.7 if 1 months old, 0.8 if 2 months old, 0.9 if 3 months old, else 1. Reduced susceptibility due to immunity in early life. \cr
#' \code{alpha_vect} = 0.4 if less than 10 years old else 0.3. Reduced susceptibility to secondary infections in each age group. \cr
#' \code{nAges}: must be obtained from \code{contact_population_list}. Number of age groups - same length as age.limits. \cr
#' \code{age.limits}: must be obtained from \code{contact_population_list}. Lower age limits in years. \cr
#' \code{matrix_mean}: must be obtained from \code{contact_population_list}. Contact matrix of mean contacts per person. \cr
#' \code{max_age}: must be obtained from \code{contact_population_list}. Maximum age in years. \cr
#' \code{age_chr}: Age limits as a character vector. \cr
#' \code{size_cohorts}: the differences in ages in days - used to calculate the size of the cohorts. \cr
#' \code{total_population}: 1861923. Total population size. \cr
#' @export
RSVsim_parameters <- function(overrides = list(),
                              contact_population_list,
                              fitted_parameter_names = NULL
                       ){

  # override parameters with any client specified ones
  if (!is.list(overrides)) {
    stop('RSVsim_parameters: overrides must be a list')
  }

  if(!is.list(contact_population_list)){
    stop("RSVsim_parameters: contact_population_list must be a list")
  }

  if(is.list(contact_population_list)){
    for(name in c("matrix_mean", "age.limits", "age_distribution")){
      if(!name %in% names(contact_population_list)){
        stop(paste("RSVsim_parameters: contact_population_list must contain", name, sep = " "))
        }
    }
  }

  nAges <- length(contact_population_list$age.limits)



  parameters <- with(contact_population_list,
       {
         if(!is.matrix(matrix_mean)){
           stop("RSVsim_parameters: contact_population_list$matrix_mean must be a square matrix with the same age groups specified in age.limits")
         }

         if(is.matrix(matrix_mean) & ncol(matrix_mean) != nAges |
            is.matrix(matrix_mean) & nrow(matrix_mean) != nAges){
           stop("RSVsim_parameters: contact_population_list$matrix_mean must be a square matrix with the same age groups specified in age.limits")
         }

         alpha_vect <- sigma_vect <- prop_detected_vect <- omega_vect <- rep(NA, nAges)

         omega_vect[age.limits >= 5] <- 0.35
         omega_vect[age.limits < 5] <- 1

         sigma_vect[age.limits <= 1/12] <- 0.7
         sigma_vect[age.limits > 1/12 & age.limits <= 2/12] <- 0.8
         sigma_vect[age.limits > 2/12 & age.limits < 3/12] <- 0.9
         sigma_vect[age.limits >= 3/12] <- 1

         prop_detected_vect[age.limits <= 3 * 1/12] <- 0.424
         prop_detected_vect[age.limits > 3 * 1/12 & age.limits <= 6 * 1/12] <- 0.088
         prop_detected_vect[age.limits > 6 * 1/12 & age.limits <= 1] <- 0.047
         prop_detected_vect[age.limits > 1 & age.limits < 2] <- 0.02
         prop_detected_vect[age.limits >= 2] <- 0.01 # very little detection in > 2 year olds

         alpha_vect[nAges >= 10] <- 0.3
         alpha_vect[nAges < 10] <- 0.4

         parameters <- list(
           "b0" = 0.08, # transmission rate coefficient 0.087
           "b1" = 0,#-0.193, # amplitude of seasonal forcing
           "phi" = 0, # phase shift of seasonal forcing
           "delta" = 1/4, # inverse of the latent period
           "gamma_s" = 1/8, # inverse of the infectious period of subsequent (secondary +) infections
           "gamma_p" = 1/9, # inverse of the infectious period of first (primary) infections
           "nu" = 1/200, # inverse of the duration of natural immunity
           "omega_vect" = omega_vect,
           "prop_detected_vect" = prop_detected_vect,
           "sigma_vect" = sigma_vect,
           "alpha_vect" = alpha_vect, # reduced susceptibility for infection in age group i
           "nAges" = nAges,
           "age.limits" = age.limits,
           "matrix_mean" = matrix_mean, # contact matrix
           "max_age" = max_age,
           "age_chr" = age_chr,
           "size_cohorts" = size_cohorts,
           "total_population" = 1861923
         )
         }
  )

  for (name in names(overrides)) {
    if (!(name %in% names(parameters))) {
      stop(paste('RSVsim_parameters: unknown parameter:', name, sep=' '))
    }

    if(name %in% c("alpha_vect", "prop_detected_vect", "sigma_vect", "omega_vect") &
       length(overrides[[name]]) != nAges){
      stop(paste("RSVsim_parameters:", name, 'is not correct length', sep = ' '))
    }

    parameters[[name]] <- overrides[[name]]
  }

  if(!is.null(fitted_parameter_names)){
    for(name in fitted_parameter_names){
      names_all <- names(parameters)
      if (!name %in% names_all) {
        stop(paste('RSVsim_parameters: unknown fitted parameter:', name, sep=' '))
      } else{
        parameters <- parameters[-which(names_all == name)]
      }
    }

  }

  return(parameters)
}
