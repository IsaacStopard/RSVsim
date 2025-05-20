#' Runs the RSV model
#'
#' @param params List of parameters from \code{get_params} function
#'
#' @return Simulation output
#' @export
run_model <- function(params){

  # times
  T0 <- 0
  max_t <- 2000
  dt <- 0.25

  # runs the model with cohort ageing

}

# move different sizes of the population proportional - so same number of population size would move each age group - so always balanced.
# have a look at the burden studies - referenced in the online document
# need to calibrate better in terms of adults vs children
# proportional roles of each - check literature
# other models have not validated the outputs in adults
# time series data in children - incidence of detected disease (no prevalence information)
# how to work out the annual burden
# Trish Campbell individual based model of RSV - explicitly linking immune status of mothers and infants
# need to assume about hospitalisation and how that varies with age

# some models also stratify by primary and secondary infections
# calculate Rt as an output - look at covid modelling - nimue model could feed in R0 - convert between beta

