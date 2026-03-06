test_that("RSVsim_hospitalisation works", {

  hosp_rate <- c(5338.91375, 1751.145228, 715.9402014, 208.4621667, 10.86702733, 2.58129891, 2.321354938, 1.868628613,
                 1.320440418, 1.834551244, 1.883359955, 2.861263391, 4.744902927, 6.645799647, 9.348531561, 14.3538506,
                 22.62406397, 32.94296492, 87.63366855) / 100000

  age_limits_hosp_rate <- c(0, 0.5, 1, 2, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75)

  step <- 1/12
  cohort_step_size <- step * 365.25
  age_times <- seq(0, 365.25 - cohort_step_size, cohort_step_size)
  age_limits <- c(seq(0, 5, step), seq(10, 75, 5))

  testthat::expect_warning(contact_population_list <- RSVsim_contact_matrix(country = "United Kingdom", age_limits = age_limits))

  time_step <- 0.25
  times_in <- seq(0, 365.25*11)
  warm_up_in <- 365.25 * 10

  parameters <- RSVsim_parameters(contact_population_list = contact_population_list)

  sim_status_quo <- RSVsim_run_model(parameters = parameters,
                                     times = times_in,
                                     cohort_step_size = cohort_step_size,
                                     warm_up = warm_up_in)

  max_time <- max(times_in)

  vaccine_times <- c(0, c(90, 120) + warm_up_in)
  vaccine_period <- diff(c(vaccine_times, max_time))
  nVaccTimes <- length(vaccine_times)
  nAges <- length(contact_population_list$age_limits)
  vaccine_cov <- rbind(matrix(rep(0, (nAges-1) * nVaccTimes),
                              nrow = (nAges-1), ncol = nVaccTimes),
                       c(0, 0.80, 0))

  # not necessary but included for completeness
  nVaccStates <- 2
  VE_inf <- matrix(rep(0.85, nAges * (nVaccStates-1)), nrow = nAges) # vaccine efficacy
  gamma_vaccine <- 1 / (365.25 * 2)

  sim_vaccination <- RSVsim_run_model(parameters = RSVsim_parameters(overrides = list("vaccine_times" = vaccine_times,
                                                                                      "vaccine_period" = vaccine_period,
                                                                                      "nVaccTimes" = nVaccTimes,
                                                                                      "vaccine_cov" = vaccine_cov,
                                                                                      "nVaccStates" = nVaccStates,
                                                                                      "VE_inf" = VE_inf,
                                                                                      "gamma_vaccine" = gamma_vaccine
                                                                                      ),
                                                                     contact_population_list = contact_population_list),
                                      times = times_in,
                                      cohort_step_size = cohort_step_size, # time at which to age people
                                      warm_up = warm_up_in)

  hosp_rate_min_time <- 0
  hosp_rate_max_time <- 365.25

  VE_hosp_given_inf <- matrix(rep(0.1, length(age_limits)), ncol = 1)

  testthat::expect_error(RSVsim_prob_hosp_given_inf(hosp_rate = hosp_rate,
                                                    age_limits_hosp_rate = c(1),
                                                    hosp_rate_min_time = hosp_rate_min_time,
                                                    hosp_rate_max_time = hosp_rate_max_time,
                                                    sim_status_quo = sim_status_quo,
                                                    VE_hosp_given_inf = VE_hosp_given_inf,
                                                    VE_hosp_vacc_states = c(2)
                                                    ))

  p_h_g_i <- RSVsim_prob_hosp_given_inf(hosp_rate = hosp_rate,
                                        age_limits_hosp_rate = age_limits_hosp_rate,
                                        hosp_rate_min_time = hosp_rate_min_time,
                                        hosp_rate_max_time = hosp_rate_max_time,
                                        sim_status_quo = sim_status_quo,
                                        VE_hosp_given_inf = VE_hosp_given_inf,
                                        VE_hosp_vacc_states = c(2))

  age_limits <- round(age_limits, digits = 5)
  age_limits_hosp_rate <- round(age_limits_hosp_rate, digits = 5)

  age_limits_index <- data.frame(age_index = base::findInterval(age_limits, age_limits_hosp_rate),
                                 age = age_limits)
  age_limits_index[, "age_hosp"] <- age_limits_hosp_rate[age_limits_index[, "age_index"]]

  pop <- p_h_g_i$sim_status_quo |> dplyr::filter(time == 0) |>
    dplyr::left_join(as.data.frame(age_limits_index), by = c("age")) |>
    dplyr::group_by(age_hosp) |>
    dplyr::summarise(tot_pop = sum(Sp + Ep + Ip + Ss + Es + Is + R))

  hosp_rate_check <- p_h_g_i$sim_status_quo |>
    dplyr::left_join(as.data.frame(age_limits_index), by = c("age")) |>
    dplyr::group_by(age_hosp) |>
    dplyr::summarise(tot_hosp = sum(Hospitalisations)) |>
    dplyr::left_join(as.data.frame(pop), by = c("age_hosp")) |>
    dplyr::mutate(hosp_rate = tot_hosp / tot_pop) |>
    dplyr::pull(hosp_rate) |>
    round(digits = 5)

  testthat::expect_true(all(hosp_rate_check == round(hosp_rate, digits = 5)))

  testthat::expect_true(all(c("prob_hosp_given_inf_uv", "sim_status_quo") %in% names(p_h_g_i)))

  testthat::expect_true(is.numeric(p_h_g_i[["sim_status_quo"]][, "Hospitalisations"]))

  testthat::expect_true(all(p_h_g_i[["sim_status_quo"]][, "Hospitalisations"] <= p_h_g_i[["sim_status_quo"]][, "Incidence"]))

  sim_vacc_hosp <- RSVsim_hospitalisations(prob_hosp_given_inf = p_h_g_i[["prob_hosp_given_inf_uv"]],
                                           age_limits_hosp_rate = age_limits_hosp_rate,
                                           sim_scenario = sim_vaccination,
                                           VE_hosp_given_inf = VE_hosp_given_inf,
                                           VE_hosp_vacc_states = c(2))

  testthat::expect_true(all(is.numeric(sim_vacc_hosp[, "Hospitalisations"])))

  sim_status_quo_hosp <- RSVsim_hospitalisations(prob_hosp_given_inf = p_h_g_i[["prob_hosp_given_inf_uv"]],
                                                 age_limits_hosp_rate = age_limits_hosp_rate,
                                                 sim_scenario = sim_status_quo,
                                                 VE_hosp_given_inf = VE_hosp_given_inf,
                                                 VE_hosp_vacc_states = c(2))

  testthat::expect_true(all(p_h_g_i[["sim_status_quo"]][, "Hospitalisations"] == sim_status_quo_hosp[,"Hospitalisations"]))

  })
