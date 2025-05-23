test_that("get_parameters_function", {

  age.limits <- c(seq(0,5,1/12), seq(10,70,5))

  mixing <- create_contact_matrix(country = "United Kingdom",
                                  age.limits = age.limits)$adjusted_matrix_mean

  expect_error(get_parameters(overrides = list("t" = 10), age.limits = age.limits, mixing = mixing))

  expect_error(get_parameters(overrides = list("omega_vect" = 1), age.limits = age.limits, mixing = mixing))

  expect_true(sum(is.na(params)) == 0)

  expect_true(get_parameters(overrides = list("b0" = 1), age.limits = age.limits, mixing = mixing)$b0 == 1)
})
