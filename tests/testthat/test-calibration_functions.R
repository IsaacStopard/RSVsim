test_that("calibration functions", {

  #skip_if(Sys.getenv("RUN_TESTS") != "true", message = "Skipping tests (set environment to RUN_TESTS to true to run)")

  testthat::expect_warning(contact_population_list <- RSVsim_contact_matrix())

  parameters <- RSVsim_parameters(contact_population_list = contact_population_list)

  sim <- RSVsim_run_model(parameters = parameters,
                          times = seq(0, 365*4, 0.5),
                          cohort_step_size = min(parameters$size_cohorts),
                          warm_up = 365 * 3)

  total_incidence <- RSVsim_total_incidence(sim)
  peak <- RSVsim_peak(sim)
  amplitude <- RSVsim_amplitude(sim)

  testthat::expect_true(length(total_incidence) == parameters$nAges)
  testthat::expect_true(length(peak) == parameters$nAges)
  expect_true(length(amplitude) == parameters$nAges)

  testthat::expect_true(sum(is.na(total_incidence)) == 0)
  testthat::expect_true(sum(is.na(peak)) == 0)
  testthat::expect_true(sum(is.na(amplitude)) == 0)

  testthat::expect_true(sum(RSVsim_shortest_periodic_dist_fun(peak, peak, 365.25)) == 0)
  testthat::expect_true(sum(RSVsim_abs_dist_fun(total_incidence, total_incidence)) == 0)

  prior_fun <- function(n_prior_attempts){
    return(matrix(c(runif(n_prior_attempts, 0.05, 0.1)), ncol = 1))
  }

  fitted_parameter_names <- c("b0")

  fixed_parameter_list <- RSVsim_parameters(contact_population_list = contact_population_list)

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
                                warm_up = 365 * 3)

  check_p <- RSVsim_ABC_rejection(target = total_incidence,
                                epsilon = total_incidence * 0.5,
                                summary_fun = RSVsim_total_incidence,
                                dist_fun = RSVsim_abs_dist_fun,
                                prior_fun = prior_fun,
                                n_prior_attempts = 10000,
                                nparticles = 2,
                                used_seeds_all = c(1,2),
                                ncores = 2,
                                fitted_parameter_names = fitted_parameter_names,
                                fixed_parameter_list = fixed_parameter_list,
                                times = seq(0, 365*4, 0.25), # maximum time to run the model for
                                cohort_step_size = min(parameters$size_cohorts), # time at which to age people\
                                warm_up = 365 * 3)

  testthat::expect_true(sum(is.na(check)) == 0)
  testthat::expect_true(all(dim(check) == c(1, 4)))
  testthat::expect_true(check[1, "particle_number"] == 1)

  prior_dens_fun <- function(x){
    return(c(dunif(x[1], 0.05, 0.1, log = FALSE)))
  }

  check_smc_p <- RSVsim_ABC_SMC(
    target = total_incidence,
    epsilon_matrix = matrix(c(total_incidence * 0.75, total_incidence * 0.5), nrow = 2),
    summary_fun = RSVsim_total_incidence,
    dist_fun = RSVsim_abs_dist_fun,
    prior_fun = prior_fun,
    n_prior_attempts = 10000,
    nparticles = 10,
    used_seed_matrix = matrix(seq(1,25*2), nrow = 2),
    prior_dens_fun = prior_dens_fun,
    particle_low = c(0.05),
    particle_up = c(0.1),
    ncores = 2,
    fitted_parameter_names = fitted_parameter_names,
    fixed_parameter_list = fixed_parameter_list,
    times = seq(0, 365*4, 0.5), # maximum time to run the model for
    cohort_step_size = min(parameters$size_cohorts), # time at which to age people\
    warm_up = 365 * 3
  )

  check_smc <- RSVsim_ABC_SMC(
    target = total_incidence,
    epsilon_matrix = matrix(c(total_incidence * 0.75, total_incidence * 0.5), nrow = 2),
    summary_fun = RSVsim_total_incidence,
    dist_fun = RSVsim_abs_dist_fun,
    prior_fun = prior_fun,
    n_prior_attempts = 10000,
    nparticles = 5,
    used_seed_matrix = matrix(seq(1,25*2), nrow = 2),
    prior_dens_fun = prior_dens_fun,
    particle_low = c(0.05),
    particle_up = c(0.1),
    ncores = 1,
    fitted_parameter_names = fitted_parameter_names,
    fixed_parameter_list = fixed_parameter_list,
    times = seq(0, 365*4, 0.5), # maximum time to run the model for
    cohort_step_size = min(parameters$size_cohorts), # time at which to age people\
    warm_up = 365 * 3
  )

})
