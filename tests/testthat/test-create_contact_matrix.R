test_that("contact matrix", {
  t_mat <- create_contact_matrix(country = "United Kingdom",
                                  age.limits = c(seq(0, 5, 1/12), seq(10, 75, 5)))

  # checking the matrix is numeric
  expect_true(
    is.numeric(t_mat$adjusted_matrix_per_capita)
  )

  # testing the matrix is symmetrical
  expect_true(
    round(t_mat$adjusted_matrix_per_capita, digits = 10) == t(round(t_mat$adjusted_matrix_per_capita, digits = 10)) |> all()
    )


  expect_true()
})
