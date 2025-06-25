test_that("calibration functions", {
  age.limits <- seq(0, 80, 10)
  contact_population_list <- RSVsim_contact_matrix(age.limits = age.limits)
  parameters <- RSVsim_parameters(contact_population_list = contact_population_list,
                                  overrides = list("b1" = 0.25))
  model_simulation <- RSVsim_run_model_dust(parameters,
                                            times = seq(0, 365*4, 0.25),
                                            cohort_step_size = 10,
                                            init_conds = NULL,
                                            warm_up = 365 * 3)

  df <- data.frame(time = model_simulation$time,
                   incidence = as.integer(rnbinom(n = nrow(model_simulation), mu = model_simulation$Incidence, size = 100)),
                   age_chr = model_simulation$age_chr) |>
    subset(time %% 5 == 0)

  fixed_parameter_list <- RSVsim_calibration_parameters(data = df,
                                                        fitted_parameter_names = c("b0", "b1"),
                                                        overrides = list(),
                                                        country = "United Kingdom",
                                                        warm_up = 365 * 3,
                                                        cohort_step_size = 10)

  expect_true(sum(is.na(fixed_parameter_list)) == 0)

  ll <- RSVsim_log_likelihood_dust(fitted_parameters = c("b0" = 0.15, "b1" = 0.25),
                                   data = df,
                                   fitted_parameter_names = fixed_parameter_list$fitted_parameter_names,
                                   fixed_parameters = fixed_parameter_list$fixed_parameters,
                                   times = fixed_parameter_list$times,
                                   cohort_step_size = fixed_parameter_list$cohort_step_size,
                                   init_conds = NULL,
                                   warm_up = fixed_parameter_list$warm_up)

  expect_true(is.numeric(ll))

  fitted_values <- RSVsim_max_likelihood_dust(data = df,
                                              fitted_parameter_names = fixed_parameter_list$fitted_parameter_names,
                                              fixed_parameters = fixed_parameter_list$fixed_parameters,
                                              times = fixed_parameter_list$times,
                                              cohort_step_size = fixed_parameter_list$cohort_step_size,
                                              init_conds = NULL,
                                              warm_up = fixed_parameter_list$warm_up,
                                              scale_parameters = FALSE,
                                              lower_ll = c(0.01, 0),
                                              upper_ll = c(10, 1))

  expect_true(is.numeric(fitted_values$par))

  })
