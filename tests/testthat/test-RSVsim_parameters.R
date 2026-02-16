test_that("RSVsim_parameters_function", {

  age.limits <- c(seq(0,5,1/12), seq(10,70,5))

  contact_population_list <- RSVsim_contact_matrix(country = "United Kingdom",
                                                   age.limits = age.limits)

  expect_error(RSVsim_parameters(overrides = list("vaccine_cov" = 1), contact_population_list = contact_population_list))

  expect_error(RSVsim_parameters(overrides = list("t" = 10), contact_population_list = contact_population_list))

  expect_error(RSVsim_parameters(overrides = list("omega_vect" = 1), contact_population_list = contact_population_list))

  expect_true(sum(is.na(RSVsim_parameters(contact_population_list = contact_population_list))) == 0)

  expect_true(RSVsim_parameters(overrides = list("b0" = 1), contact_population_list = contact_population_list)$b0 == 1)
})
