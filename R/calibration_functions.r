#' Runs the RSV model using odin and dust
#'
#' Function to run the transmission model with cohort aging. Prevalence, incidence between the given time-steps and the incidence rate per day are also calculated.
#'
#' @param parameters List of parameters from \code{get_params} function.
#' @param times Simulation times. Default: 0 - 3650 days with intervals of 0.25 days.
#' @param cohort_step_size Time steps to run the model over before adjusting the ages of all cohorts. Default: 10 days.
#' @param init_conds Initial conditions to run the model. List. Default: \code{NULL}.
#' If \code{NULL}: 1% RSV prevalence is assumed for people during the primary infection, which is seeded at the beginning of the simulation.
#' All other people are assumed to be susceptible to their primary infection.
#' @param warm_up Length of time-points to exclude before calculating the likelihood. Default: \code{NULL}.

#' @return Simulation output (dataframe). In the dataframe, age refers to the lowest age in the age group.
#' @export
#'
RSVsim_ <- function(){

}
