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
#' \code{rel_sizes}: the relative size of each cohort (age category duration in days divided by the total number of days). \cr
#' \code{total_population}: 1861923. Total population size. \cr
#' \code{Sp0, Ss0, Ep0, Es0, Ip0, Is0, R0, Incidence0}: initial conditions to run the model for each compartment - these are given as prevalence and the initial conditions calculated in this function. List. Default: \code{NULL}.
#' @export
RSVsim_parameters <- function(overrides = list(),
                              contact_population_list,
                              fitted_parameter_names = NULL
                       ){

  # override parameters with any client specified ones
  if(!is.list(overrides)) {
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

         sigma_vect[age.limits < 1/12] <- 0.7
         sigma_vect[age.limits >= 1/12 & age.limits < 2/12] <- 0.8
         sigma_vect[age.limits >= 2/12 & age.limits < 3/12] <- 0.9
         sigma_vect[age.limits > 3/12] <- 1

         prop_detected_vect[age.limits < 3 * 1/12] <- 0.424
         prop_detected_vect[age.limits >= 3 * 1/12 & age.limits < 6 * 1/12] <- 0.088
         prop_detected_vect[age.limits >= 6 * 1/12 & age.limits < 1] <- 0.047
         prop_detected_vect[age.limits >= 1 & age.limits < 2] <- 0.02
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
           "rel_sizes" = rel_sizes,
           "total_population" = 1861923
         )
         }
  )




  for (name in names(overrides)) {

    if (!(name %in% c(names(parameters), "Sp0", "Ep0", "Ip0", "Ss0", "Es0", "Is0", "R0", "Incidence0"))){
      stop(paste('RSVsim_parameters: unknown parameter:', name, sep=' '))
    }

    if(name %in% c("alpha_vect", "prop_detected_vect", "sigma_vect", "omega_vect", "Sp0", "Ep0", "Ip0", "Ss0", "Es0", "Is0", "R0", "Incidence0") &
       length(overrides[[name]]) != nAges){
      stop(paste("RSVsim_parameters:", name, 'is not correct length', sep = ' '))
    }

    parameters[[name]] <- overrides[[name]]
  }

  # default initial conditions
  rep_z <- rep(0, parameters$nAges)


  if(is.null(parameters$Sp0)){
    parameters$Sp0 = contact_population_list$rel_sizes * parameters$total_population * 0.999
  }

  if(is.null(parameters$Ep0)){
    parameters$Ep0 <- rep_z
    }

  if(is.null(parameters$Ip0)){
    parameters$Ip0 <- contact_population_list$rel_sizes * parameters$total_population * 0.001
    }

    if(is.null(parameters$Ss0)){
      parameters$Ss0 <- rep_z
      }

    if(is.null(parameters$Es0)){
      parameters$Es0 <- rep_z
      }

    if(is.null(parameters$Is0)){
      parameters$Is0 <- rep_z
    }

    if(is.null(parameters$R0)){
      parameters$R0 <- rep_z
    }

    if(is.null(parameters$Incidence0)){
      parameters$Incidence0 <- rep_z
    }

  with(parameters,
       {
         if((round(sum(Sp0) + sum(Ss0) + sum(Es0) + sum(Ep0) + sum(Ip0) + sum(Is0) + sum(R0), digits = 2)) != round(total_population, digits = 2)){
           stop("The total number of people in the initial conditions does not sum to the total population")
         }
       }
       )

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
