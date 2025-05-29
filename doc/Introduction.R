## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(RSVsim)
library(ggplot2)

# specifying the some age limits to use with the model
age.limits <- c(seq(0, 1, 0.2), seq(10, 60, 20))

## -----------------------------------------------------------------------------
contact_population_list <- RSVsim::create_contact_matrix(country = "United Kingdom", age.limits = age.limits)

matrix_mean_plot <- contact_population_list$matrix_mean |> 
  as.data.frame() |> 
  tibble::rownames_to_column("age_y") |> 
  tidyr::pivot_longer(-age_y, names_to = "age_x", values_to = "mean_contacts")

ggplot(data = matrix_mean_plot, 
aes(x = age_x, y = age_y, fill = mean_contacts)) + 
geom_tile() +
  xlab("Age") + ylab("Age") + theme_bw()

## -----------------------------------------------------------------------------
parameters <- RSVsim::get_parameters(overrides = list("b0" = 1), 
                                     contact_population_list = contact_population_list)

## -----------------------------------------------------------------------------
out <- RSVsim::run_model(parameters = parameters,
                         max_t = 3650, # maximum time to run the model for
                         cohort_step_size = 365.25/12,
                         dt = 0.5, # time step to get model outputs for
                         init_conds = NULL)

ggplot(data = out |> dplyr::filter(state %in% c("Sp", "Ip", "Ss", "Is") & time > (3650 - 365 * 3)), 
       aes(x = time, y = value, col = state)) +
  geom_line() +
  facet_wrap(~age_chr, scales = "free_y") + theme_bw()

## -----------------------------------------------------------------------------
ggplot(data = out |> dplyr::filter(state %in% c("Incidence", "DetIncidence") & time > (3650 - 365 * 3)), 
       aes(x = time, y = value, col = state)) +
  geom_line() +
  facet_wrap(~age_chr, scales = "free_y") + theme_bw()

