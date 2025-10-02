test_that("calibration functions", {

  skip_if(Sys.getenv("RUN_TESTS") != "true", message = "Skipping tests (set environment to RUN_TESTS to true to run)")

  contact_population_list <- RSVsim_contact_matrix()

  parameters <- RSVsim_parameters(contact_population_list = contact_population_list)

  sim <- RSVsim_run_model(parameters = parameters,
                          times = seq(0, 365*4, 0.25),
                          cohort_step_size = min(parameters$size_cohorts),
                          init_conds = NULL,
                          warm_up = 365 * 3)

  total_incidence <- RSVsim_total_incidence(sim)
  peak <- RSVsim_peak(sim)
  amplitude <- RSVsim_amplitude(sim)

  expect_true(length(total_incidence) == parameters$nAges)
  expect_true(length(peak) == parameters$nAges)
  expect_true(length(amplitude) == parameters$nAges)

  expect_true(sum(is.na(total_incidence)) == 0)
  expect_true(sum(is.na(peak)) == 0)
  expect_true(sum(is.na(amplitude)) == 0)

  expect_true(sum(RSVsim_shortest_periodic_dist_fun(peak, peak, 365.25)) == 0)
  expect_true(sum(RSVsim_abs_dist_fun(total_incidence, total_incidence)) == 0)

  prior_fun <- function(n_prior_attempts){
    return(as.matrix(runif(n_prior_attempts), 0.01, 0.1))
  }

  fitted_parameter_names <- c("b0")

  fixed_parameter_list <- RSVsim_parameters(contact_population_list = contact_population_list,
                                            fitted_parameter_names = fitted_parameter_names)

  check <- RSVsim_ABC_rejection(target = total_incidence,
                                epsilon = total_incidence * 0.5,
                                summary_fun = RSVsim_total_incidence,
                                dist_fun = RSVsim_abs_dist_fun,
                                prior_fun = prior_fun,
                                n_prior_attempts = 10000,
                                nparticles = 1,
                                used_seeds_all = c(1),
                                ncores = 1,
                                fitted_parameter_names = fitted_parameter_names,
                                fixed_parameter_list = fixed_parameter_list,
                                times = seq(0, 365*4, 0.25), # maximum time to run the model for
                                cohort_step_size = min(parameters$size_cohorts), # time at which to age people\
                                init_conds = NULL,
                                warm_up = 365 * 3)

  expect_true(sum(is.na(check)) == 0)
  expect_true(all(dim(check) == c(1, 4)))
  expect_true(check[1, "particle_number"] == 1)

})
