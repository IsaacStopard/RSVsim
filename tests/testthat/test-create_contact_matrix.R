test_that("contact matrix", {

  expect_error(create_contact_matrix(country = "France",
                                     age.limits = seq(0, 70, 5)))

  expect_no_warning(create_contact_matrix(country = "United Kingdom",
                                          age.limits = seq(0, 70, 5)))

  list_age_limits <- list(c(seq(0, 5, 1), seq(10, 75, 5)),
                          c(seq(0, 5, 1/12), seq(10, 80, 5)),
                          c(seq(0, 75, 5)))

  countries <- c("United Kingdom", "Belgium", "Finland")

  for(i in 1:length(list_age_limits)){

    t_mat <- create_contact_matrix(country = countries[i],
                                   age.limits = list_age_limits[[i]])

    # checking the matrix is numeric
    testthat::expect_true(
      all(c(is.numeric(t_mat$default_matrix_mean),
            is.numeric(t_mat$default_matrix_contacts),
            is.numeric(t_mat$default_population)))
    )

    testthat::expect_true(sum(c(sum(is.na(t_mat$default_matrix_mean)),
                                sum(is.na(t_mat$default_matrix_contacts)),
                                sum(is.na(t_mat$default_population)))) == 0
                          )

    # checking matrix is symmetric
    testthat::expect_true(
      isSymmetric(t_mat$default_matrix_contacts, check.attributes = FALSE)
    )

    if(any(age.limits %% 1 != 0)){
      testthat::expect_true(
        all(c(is.numeric(t_mat$adjusted_matrix_mean),
              is.numeric(t_mat$adjusted_matrix_contacts),
              is.numeric(t_mat$adjusted_population)))
      )

       testthat::expect_true(sum(c(sum(is.na(t_mat$adjusted_matrix_mean)),
                                sum(is.na(t_mat$adjusted_matrix_contacts)),
                                sum(is.na(t_mat$adjusted_population)))) == 0
                             )

      testthat::expect_equal(nrow(t_mat$adjusted_matrix_mean), length(list_age_limits[[i]]))

      testthat::expect_true(
      isSymmetric(t_mat$adjusted_matrix_contacts, check.attributes = FALSE)
      )
    }



  }

})
