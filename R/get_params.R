#' Provides the default parameters to run the RSV model
#'
#' @param overrides list of default parameters to change: b0 (),
#' @param country country for use in the \code{contact_matrix} function in the socialmixr package. Can be given as country name or 2 digit ISO code
#' @param age.limits lower limits of the age groups to run the simulation for (must be in years)
#' @return Simulation output
#' @export

get_params <- function(overrides = list(),

                       ){

  nAges <- length(age.limits)

  params <- list(
    "b0" = 0.015421, #
    "b1" = 0.39722,
    "omega" = 0.6, # reduced infectiousness from older age groups
    "delta" = 7.632, # latency rate
    "gamma_p" = 3.342, # infectious rate
    "gamma_s" = ,
    "nu" = 0.132, # rate of waning immunity
    "phi" = 0.98456,

    "nAges" = nAges,
    "mixing" = mixing, # contact matrix

  )
  init_conds_from_file <- 0 # choose whether to read in some existing ICs
  save_init_conds <- 1 # choose whether to save final model state as ICs for next time
  max_t <- 2000

  prop_detected_vect <- rep(0, parameters[["nAges"]])

  for(i in 1)

    as.vector(c(rep(0.424, 3),
                                     rep(0.088, 3),
                                     rep(0.047, 6),
                                     rep(0.020, 12),
                                     rep(0, (nAges-24))))

  # override parameters with any client specified ones
  if (!is.list(overrides)) {
    stop('overrides must be a list')
  }

  for (name in names(overrides)) {
    if (!(name %in% names(parameters))) {
      stop(paste('unknown parameter', name, sep=' '))
    }
    parameters[[name]] <- overrides[[name]]
  }

}
