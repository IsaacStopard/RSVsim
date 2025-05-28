test_that("get_parameters_function", {

  age.limits <- c(seq(0,5,1/12), seq(10,70,5))

  contact_population_list <- create_contact_matrix(country = "United Kingdom",
                                                   age.limits = age.limits)

  expect_error(get_parameters(overrides = list("t" = 10), contact_population_list = contact_population_list))

  expect_error(get_parameters(overrides = list("omega_vect" = 1), contact_population_list = contact_population_list))

  expect_true(sum(is.na(get_parameters(contact_population_list = contact_population_list))) == 0)

  expect_true(get_parameters(overrides = list("b0" = 1), contact_population_list = contact_population_list)$b0 == 1)
})
