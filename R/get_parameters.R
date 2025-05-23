#' get parameters for RSVsim
#'
#' Provides the default parameters to run the RSV model
#'
#' @param overrides List of default parameters to change.
#' @param country Country for use in the \code{contact_matrix} function in the socialmixr package. Can be given as country name or 2 digit ISO code.
#' @param age.limits Lower limits of the age groups to run the simulation for (must be in years).
#' @param mixing Contact matrix.
#' @return Parameter list.
#' @export

get_parameters <- function(overrides = list(),
                           age.limits = c(seq(0,5,1/12), seq(10,70,5)),
                           contact_matrix
                       ){

  # override parameters with any client specified ones
  if (!is.list(overrides)) {
    stop('overrides must be a list')
  }

  if(!is.numeric(age.limits)){
    stop("age.limits must be numeric")
  }

  nAges <- length(age.limits)

  mixing <- if("adjusted_matrix_mean" %in% names(contact_matrix)){
    contact_matrix$adjusted_matrix_mean
  } else{
    contact_matrix$default_matrix_mean
  }

  if(!is.matrix(mixing) | ncol(mixing) != nAges | nrow(mixing) != nAges){
    stop("mixing must be a square matrix with the same age groups specified in age.limits")
  }

  alpha_vect <- sigma_vect <- prop_detected_vect <- omega_vect <- rep(NA, nAges)

  omega_vect[age.limits >= 5] <- 0.35
  omega_vect[age.limits < 5] <- 1

  sigma_vect[age.limits <= 1/12] <- 0.7
  sigma_vect[age.limits > 1/12 & age.limits <= 2/12] <- 0.8
  sigma_vect[age.limits > 2/12 & age.limits <= 3/12] <- 0.9
  sigma_vect[age.limits > 3/12] <- 1

  prop_detected_vect[age.limits <= 3 * 1/12] <- 0.424
  prop_detected_vect[age.limits > 3 * 1/12 & age.limits <= 6 * 1/12] <- 0.088
  prop_detected_vect[age.limits > 6 * 1/12 & age.limits <= 1] <- 0.047
  prop_detected_vect[age.limits > 1 & age.limits <= 2] <- 0.02
  prop_detected_vect[age.limits > 2] <- 0

  alpha_vect[nAges >= 10] <- 0.3
  alpha_vect[nAges < 10] <- 0.4

  parameters <- list(
    "b0" = 0.087, # transmission rate coefficient
    "b1" = -0.193, # amplitude of seasonal forcing
    "phi" = 1.536, # phase shift of seasonal forcing
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
    "mixing" = mixing, # contact matrix
    "max_age" = contact_matrix$max_age,
    "population" = 1861923
  )

  for (name in names(overrides)) {
    if (!(name %in% names(parameters))) {
      stop(paste('unknown parameter', name, sep=' '))
    }

    if(name %in% c("alpha_vect", "prop_detected_vect", "sigma_vect", "omega_vect") &
       length(overrides[[name]]) != nAges){
      stop(paste(name, 'is not correct length', sep = ' '))
    }

    parameters[[name]] <- overrides[[name]]
  }

  return(parameters)
}
