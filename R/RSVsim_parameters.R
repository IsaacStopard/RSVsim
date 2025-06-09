#' Get parameters for RSVsim
#'
#' Provides the default parameters to run the RSV model.
#'
#' @param overrides List of default parameters to change.
#' @param contact_population_list List of outputs from the \code{create_contact_matrix} function.
#' @param fitted Vector of parameter names to exclude from the list because they are fitted. Default: \code{NULL}.
#' @return Parameter list.
#' @export

RSVsim_parameters <- function(overrides = list(),
                           contact_population_list,
                           fitted = NULL
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
         prop_detected_vect[age.limits >= 2] <- 0

         alpha_vect[nAges >= 10] <- 0.3
         alpha_vect[nAges < 10] <- 0.4

         parameters <- list(
           "b0" = 0.15, # transmission rate coefficient 0.087
           "b1" = -0.193, # amplitude of seasonal forcing
           "phi" = 1.536 * 30.436875, # phase shift of seasonal forcing
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

  if(!is.null(fitted)){
    for(name in fitted){
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
