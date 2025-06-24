test_that("multiplication works", {
  age.limits <- seq(0, 80, 10)
  contact_population_list <- RSVsim_contact_matrix(age.limits = age.limits)
  parameters <- RSVsim_parameters(contact_population_list = contact_population_list)
  model_simulation <- RSVsim_run_model_dust(parameters,
                                            warm_up = 365 * 9,
                                            max_t = 365 * 10)

  df <- data.frame(time = model_simulation$time,
                   incidence = as.integer(rnbinom(n = nrow(model_simulation), mu = model_simulation$Incidence, size = 100)),
                   age_chr = model_simulation$age_chr) |>
    subset(time %% 5 == 0)

  data <- list("sim" = df)

  fixed_parameter_list <- RSVsim_calibration_parameters(data = data,
                                                        data_populations = rep(parameters$total_population, length(data)),
                                                        warm_up = 365 * 9,
                                                        cohort_step_size = 10,
                                                        init_conds = NULL)

  expect_equal(length(data), length(fixed_parameter_list))
  expect_true(sum(is.na(fixed_parameter_list)) == 0)

  ll <- RSVsim_log_likelihood_dust(b0 = 0.15,
                                   b1 = 0.25,
                                   phi = 10,
                                   fixed_parameter_list = fixed_parameter_list,
                                   data = data)

  expect_true(is.numeric(ll))

  fitted_values <- RSVsim_max_likelihood_dust(fixed_parameter_list = fixed_parameter_list,
                                              data = data,
                                              scale_parameters = list(lower = c(0.01, 0, 0), upper = c(10, 1, 365.25)))

  expect_true(is.numeric(fitted_values$par))

  })
