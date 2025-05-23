test_that("get_parameters_function", {

  age.limits <- c(seq(0,5,1/12), seq(10,70,5))

  contact_matrix <- create_contact_matrix(country = "United Kingdom",
                                  age.limits = age.limits)

  expect_error(get_parameters(overrides = list("t" = 10), age.limits = age.limits, contact_matrix = contact_matrix))

  expect_error(get_parameters(overrides = list("omega_vect" = 1), age.limits = age.limits, contact_matrix = contact_matrix))

  expect_true(sum(is.na(get_parameters(contact_matrix = contact_matrix))) == 0)

  expect_true(get_parameters(overrides = list("b0" = 1), age.limits = age.limits, contact_matrix = contact_matrix)$b0 == 1)
})
