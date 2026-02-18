test_that("RSVsim_run_model function works", {

  #skip_if(Sys.getenv("RUN_TESTS") != "true", message = "Skipping tests (set environment to RUN_TESTS to true to run)")

  testthat::expect_warning(contact_population_list <- RSVsim_contact_matrix())
  parameters <- RSVsim_parameters(overrides = list("b0" = 0.11, "b1" = 0.3, "vaccine_times" = c(0.00, c(365.25, 395.25, 730.50, 760.50) + 365.25)),
                                  contact_population_list = contact_population_list)

  vacc_cov_old <- rbind(matrix(rep(0, (length(contact_population_list$age_limits)-1) * 5),
                               nrow = (length(contact_population_list$age_limits)-1)),
                        c(0, 0.99, 0, 0.99, 0))

  #c_fun <- cinterpolate::interpolation_function(x = parameters_vac$vaccine_times, y = vaccine_rate, type = "constant")
  #c_fun(425)

  parameters_vac <- RSVsim_parameters(overrides = list("b0" = 0.11, "b1" = 0.3, "vaccine_cov" = vacc_cov_old,
                                                       "vaccine_times" = c(0.00, c(365.25, 395.25, 730.50, 760.50) + 365.25)),
                                      contact_population_list = contact_population_list)

  s_time <- Sys.time()
  sim <- RSVsim_run_model(parameters = parameters,
                          times = seq(0, 365.25*5, 0.25),
                          cohort_step_size = 1/12 * 365.25,
                          warm_up = 365.25 * 1)
  e_time <- Sys.time()
  e_time - s_time

  testthat::expect_true(sum(is.na(subset(sim, vacc_state == 1))) == 0)
  testthat::expect_true(max(sim$time) <= 3650)

  sim_vac <- RSVsim_run_model(parameters = parameters_vac,
                              times = seq(0, 365.25*5, 0.25),
                              cohort_step_size = 1/12 * 365.25,
                              warm_up = 365.25 * 1)

  testthat::expect_true(sum(sim_vac$doses) > 0)

  testthat::expect_error(RSVsim_run_model(parameters = parameters,
                                times = seq(0, 3650, 0.25),
                                cohort_step_size = min(parameters$size_cohorts) + 10,
                                warm_up = 365 * 3))

  testthat::expect_true(max(na.omit(sim$prev)) <= 1 & min(na.omit(sim$prev)) >= 0)

  # checking the total population is correct
  testthat::expect_true(all(abs(sim |> dplyr::group_by(time) |> dplyr::summarise(Sp = sum(Sp), Ep = sum(Ep), Ip = sum(Ip), Ss = sum(Ss), Es = sum(Es), Is = sum(Is), R = sum(R)) |>
              dplyr::select(Sp, Ep, Ip, Ss, Es, Is, R) |> rowSums() - parameters$total_population) < 1E-5)
              )

})
