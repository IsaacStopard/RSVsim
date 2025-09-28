# age.limits <- c(seq(0, 2, 0.5), seq(10, 60, 5))
# contact_population_list <- RSVsim_contact_matrix(country = "United Kingdom", age.limits = age.limits)
#
# parameters <- RSVsim_parameters(overrides = list("b0" = 0.08),
#                                 contact_population_list = contact_population_list,
#                                 fitted = NULL)
# nAges <- length(age.limits)
#
# out <- RSVsim_run_model(parameters = parameters,
#                         times = seq(0, 365*5, 0.25), # maximum time to run the model for
#                         cohort_step_size = 0.2*365, # time at which to age people\
#                         init_conds = list("Sp0" = rep(0, nAges),
#                                           "Ep0" = rep(0, nAges),
#                                           "Ip0" = rep(0.01, nAges),
#                                           "Ss0" = rep(1 - 0.01 - 0.1, nAges),
#                                           "Es0" = rep(0, nAges),
#                                           "Is0" = rep(0, nAges),
#                                           "R0" = rep(0.1, nAges),
#                                           "Incidence0" = rep(0, nAges)),
#                         warm_up = 365*4)
#
# ggplot(data = out, aes(x = time, y = Incidence, col = age, group = age_chr)) +
#   geom_line() +
#   facet_wrap(~age_chr, scales = "free_y")
#
# fitted_parameter_names <- c("b0")
#
# fixed_parameter_list <- RSVsim_fixed_parameters(fitted_parameter_names = fitted_parameter_names, data = out)
#
# calibration_targets <- RSVsim_peak_plus_total_incidence(fitted_parameters = c(0.15),
#                                                         fitted_parameter_names = fitted_parameter_names,
#                                                         fixed_parameter_list = fixed_parameter_list,
#                                                         init_conds = list("Sp0" = rep(0, nAges),
#                                                                           "Ep0" = rep(0, nAges),
#                                                                           "Ip0" = rep(0.01, nAges),
#                                                                           "Ss0" = rep(1 - 0.01 - 0.1, nAges),
#                                                                           "Es0" = rep(0, nAges),
#                                                                           "Is0" = rep(0, nAges),
#                                                                           "R0" = rep(0.1, nAges),
#                                                                           "Incidence0" = rep(0, nAges)))
#
# RSVsim_peak_plus_total_incidence_wrapper <- function(x){
#   return(
#     RSVsim_peak_plus_total_incidence(fitted_parameters = x,
#                                      fitted_parameter_names = fitted_parameter_names,
#                                      fixed_parameter_list = fixed_parameter_list,
#                                      times = seq(0, 365*1, 0.25), # maximum time to run the model for
#                                      cohort_step_size = 10, # time at which to age people\
#                                      init_conds = list("Sp0" = rep(0, nAges),
#                                                        "Ep0" = rep(0, nAges),
#                                                        "Ip0" = rep(0.01, nAges),
#                                                        "Ss0" = rep(1 - 0.01 - 0.1, nAges),
#                                                        "Es0" = rep(0, nAges),
#                                                        "Is0" = rep(0, nAges),
#                                                        "R0" = rep(0.1, nAges),
#                                                        "Incidence0" = rep(0, nAges)),
#                                      warm_up = NULL)
#   )
# }
#
# n_t <- length(calibration_targets)/2
#
# # Run the ABC simulation
# # nb_simul: The number of simulations to run
# # summary_stat_target: Your observed summary statistics
# # tol: The proportion of simulations to keep (e.g., 10%)
# abc_results_t0.3 <- EasyABC::ABC_rejection(model = RSVsim_peak_plus_total_incidence_wrapper,
#                                        prior = list(c("unif", 0, 5)),
#                                        nb_simul = 5000,
#                                        summary_stat_target = calibration_targets,
#                                        tol = 0.3,
#                                        progress_bar = TRUE,
#                                        n_cluster=1
#                                        )
#
# hist(abc_results_t0.3$param)
# quantile(abc_results_t0.3$param, probs = c(0.025, 0.5, 0.975))
#
# abc_results_t0.1 <- EasyABC::ABC_rejection(model = RSVsim_peak_plus_total_incidence_wrapper,
#                                            prior = list(c("unif", 0, 5)),
#                                            nb_simul = 5000,
#                                            summary_stat_target = calibration_targets,
#                                            tol = 0.1,
#                                            progress_bar = TRUE,
#                                            n_cluster=1
# )
#
# hist(abc_results_t0.1$param)
# quantile(abc_results_t0.1$param, probs = c(0.025, 0.5, 0.975))
#
# tolerance=0.5
# c_drov=0.7
#
# abc_smc_results <- EasyABC::ABC_sequential(method = "Drovandi",
#                                           model = RSVsim_peak_plus_total_incidence_wrapper,
#                                           prior = list(c("unif", 0, 5)),
#                                           nb_simul = 500,
#                                           summary_stat_target = calibration_targets,
#                                           tolerance_tab = 0.5,
#                                           progress_bar = TRUE,
#                                           n_cluster = 1,
#                                           c = 0.7
#                                           )
#
# # simulating the data
# age.limits <- c(seq(0, 2, 0.2), seq(10, 60, 20))
#
