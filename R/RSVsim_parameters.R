#' Get parameters for RSVsim
#'
#' Provides the default parameters to run the RSV model.
#'
#' @param overrides List of default parameters to change.
#' @param contact_population_list List of outputs from the \code{create_contact_matrix} function.
#' @return Parameter list. Default values: \cr
#' \code{b0}: 0.08. Transmission rate coefficient. \cr
#' \code{b1}: 0. Amplitude of seasonal forcing in the transmission rate. \cr
#' \code{phi}: 0. Phase shift of the seasonal forcing. \cr
#' \code{delta}: 1/4. Inverse of the latent period. \cr
#' \code{gamma_s}: 1/8. Inverse of the infectious period of subsequent (secondary) infections. \cr
#' \code{gamma_p}: 1/9. Inverse of the infectious period of first (primary) infections. \cr
#' \code{nu}: 1/200. Inverse of the duration of natural immunity. \cr
#' \code{omega_vect}: 0.35 if greater than 5 years old else 1. Reduced infectiousness due to immunity in early life. \cr
#' \code{prop_detected_vect}: age-specific proportion detected. \cr
#' \code{sigma_vect}: 0.7 if 1 months old, 0.8 if 2 months old, 0.9 if 3 months old, else 1. Reduced susceptibility due to immunity in early life. \cr
#' \code{alpha_vect} = 0.4 if less than 10 years old else 0.3. Reduced susceptibility to secondary infections in each age group. \cr
#' \code{nAges}: must be obtained from \code{contact_population_list}. Number of age groups - same length as age_limits. \cr
#' \code{age_limits}: must be obtained from \code{contact_population_list}. Lower age limits in years. \cr
#' \code{matrix_per_person}: must be obtained from \code{contact_population_list}. Contact matrix of mean contacts per person. \cr
#' \code{max_age}: must be obtained from \code{contact_population_list}. Maximum age in years. \cr
#' \code{age_chr}: Age limits as a character vector. \cr
#' \code{size_cohorts}: the differences in ages in days - used to calculate the size of the cohorts. \cr
#' \code{rel_sizes}: the relative size of each cohort (age category duration in days divided by the total number of days). \cr
#' \code{total_population}: 1861923. Total population size. \cr
#' \code{nVaccStates}: 2. Number of vaccination states. Default: 1: unvaccinated, 2: vaccinated. Minimum is 2. \cr
#' \code{gamma_vaccine}: 1 / (365*2). Rate of waning between vaccinated states. Size of vector must be nVaccStates minus one. \cr
#' \code{nVaccTimes}: 5. Length of start times of vaccination distributions. \cr
#' \code{vaccine_times}: c(0, 365.25, 365.25 + 30, 730.5, 730.5 + 30). Times of vaccination distributions. First time should be 0. \cr
#' \code{vaccine_period}: rep(30, 5). Duration of vaccination distributions.
#' \code{vaccine_cov}: matrix with the rows corresponding to the vaccine_times and the columns corresponding to the ages. \cr
#' Age-specific proportion of unvaccinated people to have been vaccinated by the end of the vaccination distribution for each \code{vaccine_time},
#' assuming no waning of vaccination. The vaccination rate is calculated as \code{-log(1 - vaccine_cov) / vaccine_period}. Default is 0 coverage for all ages.
#' Changes in effective coverage with are given as a model output.  \cr
#' \code{VE}: 0.85 for all ages and vaccinated states. Vaccine efficacy for each age group and vaccinated state. Number of rows must be equal to the number of age groups,
#' and the number of columns should be equal to nVaccStates - 1. \cr
#' \code{Sp0, Ss0, Ep0, Es0, Ip0, Is0, R0, Incidence0, doses0}: initial conditions to run the model for each compartment -
#' these are given as prevalence and the initial conditions calculated in this function. List. Default: \code{NULL}.
#' If \code{NULL}: 0.1% RSV prevalence is assumed for people during the primary infection, which is seeded at the beginning of the simulation.
#' All other people are assumed to be susceptible to their primary infection.
#' @export
RSVsim_parameters <- function(overrides = list(),
                              contact_population_list
                       ){

  # override parameters with any user specified ones
  if(!is.list(overrides)) {
    stop('RSVsim_parameters: overrides must be a list')
  }

  if(!is.list(contact_population_list)){
    stop("RSVsim_parameters: contact_population_list must be a list")
  } else {
    for(name in c("matrix_per_person", "age_limits", "age_distribution")){
      if(!name %in% names(contact_population_list)){
        stop(paste("RSVsim_parameters: contact_population_list must contain", name, sep = " "))
      }
      }
  }

  if(any(diff(contact_population_list$age_limits) <= 0)){
    stop("RSVsim_parameters: age_limits must be in ascending order")
  }

  nAges <- length(contact_population_list$age_limits)

  parameters <- with(contact_population_list,
       {
         if(!is.matrix(matrix_per_person)){
           stop("RSVsim_parameters: contact_population_list$matrix_per_person must be a square matrix with the same age groups specified in age_limits")
         }

         if(is.matrix(matrix_per_person) & ncol(matrix_per_person) != nAges |
            is.matrix(matrix_per_person) & nrow(matrix_per_person) != nAges){
           stop("RSVsim_parameters: contact_population_list$matrix_per_person must be a square matrix with the same age groups specified in age_limits")
         }

         alpha_vect <- sigma_vect <- prop_detected_vect <- omega_vect <- rep(NA, nAges)

         omega_vect[age_limits >= 5] <- 0.35
         omega_vect[age_limits < 5] <- 1

         sigma_vect[age_limits < 1/12] <- 0.7
         sigma_vect[age_limits >= 1/12 & age_limits < 2/12] <- 0.8
         sigma_vect[age_limits >= 2/12 & age_limits < 3/12] <- 0.9
         sigma_vect[age_limits >= 3/12] <- 1

         prop_detected_vect[age_limits < 3 * 1/12] <- 0.424
         prop_detected_vect[age_limits >= 3 * 1/12 & age_limits < 6 * 1/12] <- 0.088
         prop_detected_vect[age_limits >= 6 * 1/12 & age_limits < 1] <- 0.047
         prop_detected_vect[age_limits >= 1 & age_limits < 2] <- 0.02
         prop_detected_vect[age_limits >= 2] <- 0.01 # very little detection in > 2 year olds

         alpha_vect[age_limits >= 10] <- 0.3
         alpha_vect[age_limits < 10] <- 0.4

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
           "age_limits" = age_limits,
           "matrix_per_person" = matrix_per_person, # contact matrix
           "max_age" = max_age,
           "age_chr" = age_chr,
           "size_cohorts" = size_cohorts,
           "rel_sizes" = rel_sizes,
           "total_population" = 1861923,
           "nVaccStates" = 2,
           "gamma_vaccine" = 1/(365*2),
           "nVaccTimes" = 5,
           "vaccine_times" = c(0, 365.25, 365.25 + 30, 730.5, 730.5 + 30),
           "vaccine_period" = rep(30, 5),
           "vaccine_cov" = matrix(rep(0, nAges * 5), nrow = nAges, ncol = 5),
           "VE" = matrix(rep(0.85, nAges), ncol = 1, nrow = nAges)
         )
         }
  )

  nVaccStates <- if(!is.null(overrides[["nVaccStates"]])){overrides[["nVaccStates"]]} else{parameters[["nVaccStates"]]}

  if(nVaccStates < 2){
    stop("nVaccStates must be at least 2 (unvaccinated and vaccinated)")
  }

  nVaccTimes <- if(!is.null(overrides[["nVaccTimes"]])){overrides[["nVaccTimes"]]} else{parameters[["nVaccTimes"]]}

  for (name in names(overrides)) {

    if (!(name %in% c(names(parameters), "Sp0", "Ep0", "Ip0", "Ss0", "Es0", "Is0", "R0", "Incidence0", "doses0"))){
      stop(paste("RSVsim_parameters: unknown parameter:", name, sep=' '))
    }

    if(name %in% c("alpha_vect", "prop_detected_vect", "sigma_vect", "omega_vect") &
       length(overrides[[name]]) != nAges){
      stop(paste("RSVsim_parameters:", name, 'is not correct length', sep = ' '))
    }

    if(name %in% c("Sp0", "Ep0", "Ip0", "Ss0", "Es0", "Is0", "R0", "Incidence0", "doses0") &&
       length(overrides[[name]]) != nAges * nVaccStates){
      stop(paste("RSVsim_parameters:", name, "is not the correct length", sep = ' '))
    }

    if(name == "VE" & length(overrides[[name]]) != nAges * (nVaccStates - 1)){
      stop("RSVsim_parameters: VE is not the correct dimensions")
    }

    if(name == "max_cov" &&
       length(overrides[[name]]) != nAges){
      stop(paste("RSVsim_parameters:", name, "is not the correct length", sep = ' '))
    }

    if(name == "vaccine_times" && length(overrides[[name]]) != nVaccTimes){
      stop("RSVsim_parameters: vaccine_times is not the correct length")
      }

    if(name %in% c("gamma_vaccine") && length(overrides[[name]]) != (nVaccStates - 1)){
      stop(paste("RSVsim_parameters:", name, "is not the correct length", sep = ' '))
    }

    if(name %in% c("vaccine_cov") &&
       length(overrides[[name]]) != nAges * nVaccTimes){
      stop(paste("RSVsim_parameters:", name, "is not the correct size", sep = ' '))
    }

    parameters[[name]] <- overrides[[name]]
  }

  # default initial conditions
  rep_z <- rep(0, parameters$nAges)
  mat_z <- matrix(rep(rep_z, nVaccStates),
                  nrow = nAges, ncol = nVaccStates)

  if(is.null(parameters$Sp0)){
    parameters$Sp0 = matrix(c(contact_population_list$rel_sizes * parameters$total_population * 0.999,
                              rep(rep_z, nVaccStates - 1)),
                            nrow = nAges, ncol = nVaccStates)

    }

  if(is.null(parameters$Ep0)){
    parameters$Ep0 <- mat_z
    }

  if(is.null(parameters$Ip0)){
    parameters$Ip0 <- matrix(c(contact_population_list$rel_sizes * parameters$total_population * 0.001,
                               rep(rep_z, nVaccStates - 1)),
                             nrow = nAges, ncol = nVaccStates)
    }

  if(is.null(parameters$Ss0)){
      parameters$Ss0 <- mat_z
    }

    if(is.null(parameters$Es0)){
      parameters$Es0 <- mat_z
      }

    if(is.null(parameters$Is0)){
      parameters$Is0 <- mat_z
    }

    if(is.null(parameters$R0)){
      parameters$R0 <- mat_z
    }

    if(is.null(parameters$Incidence0)){
      parameters$Incidence0 <- mat_z
    }

  if(is.null(parameters$doses0)){
    parameters$doses0 <- mat_z
  }

  if(max(parameters$vaccine_cov) > 1 | min(parameters$vaccine_cov)){
    stop("vaccine_cov values must be between 0 and 1")
  }

  if(max(parameters$VE) > 1 | min(parameters$VE) < 0){
    stop("VE must be between 0 and 1")
  }

  with(parameters,
       {
         if((round(sum(Sp0) + sum(Ss0) + sum(Es0) + sum(Ep0) + sum(Ip0) + sum(Is0) + sum(R0), digits = 2)) != round(total_population, digits = 2)){
           stop("The total number of people in the initial conditions does not sum to the total population")
         }
       }
       )

  return(parameters)
}
