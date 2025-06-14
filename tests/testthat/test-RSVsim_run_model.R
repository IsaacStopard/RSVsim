test_that("RSVsim_run_model function works", {

  contact_population_list <- RSVsim_contact_matrix()
  parameters <- RSVsim_parameters(contact_population_list = contact_population_list)

  sim <- RSVsim_run_model(parameters = parameters,
                          max_t = 3650,
                          cohort_step_size = 10,
                          dt = 0.25,
                          init_conds = NULL)

  expect_true(sum(is.na(sim)) == 0)
  expect_true(max(sim$time) <= 3650)

})
