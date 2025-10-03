test_that("RSVsim_run_model function works", {

  #skip_if(Sys.getenv("RUN_TESTS") != "true", message = "Skipping tests (set environment to RUN_TESTS to true to run)")

  contact_population_list <- RSVsim_contact_matrix()
  parameters <- RSVsim_parameters(contact_population_list = contact_population_list)

  sim <- RSVsim_run_model(parameters = parameters,
                          times = seq(0, 3650, 0.25),
                          cohort_step_size = 10,
                          init_conds = NULL,
                          warm_up = 365 * 3)

  expect_true(sum(is.na(sim)) == 0)
  expect_true(max(sim$time) <= 3650)

  expect_error(RSVsim_run_model(parameters = parameters,
                                times = seq(0, 3650, 0.25),
                                cohort_step_size = min(parameters$size_cohorts) + 10,
                                init_conds = NULL,
                                warm_up = 365 * 3))

  # checking the total population is correct
  if(!all(abs(sim |> dplyr::group_by(time) |> dplyr::summarise(Sp = sum(Sp), Ep = sum(Ep), Ip = sum(Ip), Ss = sum(Ss), Es = sum(Es), Is = sum(Is), R = sum(R)) |>
              dplyr::select(Sp, Ep, Ip, Ss, Es, Is, R) |> rowSums() - parameters$total_population) < 1E-5)){
    stop("RSVsim_run_model: population does not sum to the correct number")
  }

})
