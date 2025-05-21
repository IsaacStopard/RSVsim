test_that("contact matrix", {

  list_age_limits <- list(c(seq(0, 5, 1), seq(10, 75, 5)),
                          c(seq(0, 5, 1/12), seq(10, 90, 5)),
                          c(seq(0, 75, 5)))

  countries <- c("United Kingdom", "Belgium", "Finland")

  expect_error(create_contact_matrix(country = "France",
                                     age.limits = seq(0, 70, 5)))

  for(i in 1:length(list_age_limits)){
    t_mat <- create_contact_matrix(country = countries[i],
                                   age.limits = list_age_limits[[i]])

    # checking the matrix is numeric
    testthat::expect_true(
      is.numeric(t_mat$adjusted_matrix)
    )

    testthat::expect_true(
      is.numeric(t_mat$default_matrix)
    )

    testthat::expect_equal(nrow(t_mat$adjusted_matrix), length(list_age_limits[[i]]))

    testthat::expect_true(sum(is.na(t_mat$adjusted_matrix)) == 0)

    # checking matrix is symmetric
    testthat::expect_true(
      all(t(round(t_mat$default_matrix, digits = 10)) == round(t_mat$default_matrix, digits = 10))
    )

    # checking the rows sums are consistent
    age.limits.default <- colnames(t_mat$default_matrix)

    testthat::expect_true(
      all(
        round(rowSums(t_mat$adjusted_matrix[which(rownames(t_mat$adjusted_matrix) %in% age.limits.default),]), digits = 10) ==
          round(rowSums(t_mat$default_matrix), digits = 10))
    )
  }

})
