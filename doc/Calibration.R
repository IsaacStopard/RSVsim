## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(RSVsim)
library(ggplot2)


## -----------------------------------------------------------------------------
# getting the model parameters
age.limits <- c(seq(0, 5, 1), seq(10, 70, 20))

contact_population_list <- create_contact_matrix(country = "United Kingdom", age.limits = age.limits)

parameters <- get_parameters(overrides = list("b0" = 0.15, 
                                              "b1" = 0.25, 
                                              "phi" = 10),
                             contact_population_list = contact_population_list)

# running the model
model_simulation <- run_model(parameters = parameters, 
                              warm_up = 365 * 9, 
                              max_t = 365 * 10)

# adding noise to the incidence values and selecting samples every 5 days
set.seed(123)
df <- data.frame(time = model_simulation$time,
                 incidence = as.integer(rnbinom(n = nrow(model_simulation), mu = model_simulation$Incidence, size = 50)),
                 age_chr = model_simulation$age_chr) |> 
  subset(time %% 5 == 0)



## -----------------------------------------------------------------------------

data <- list("sim" = df)

fixed_parameter_list <- get_calibration_parameters(data = data,
                                                   data_populations = rep(parameters$total_population, length(data)),
                                                   warm_up = 365 * 9)

calibration_likelihood(fitted_parameters = c("b0" = 0.15, "b1" = 0.25, "phi" = 10),
                       fixed_parameter_list = fixed_parameter_list,
                       minimise = FALSE,
                       data = data,
                       cohort_step_size = 10,
                       dt = 0.25) |> print()

fitted_values <- constrained_max_likelihood(fixed_parameter_list = fixed_parameter_list,
                                            data = data,
                                            scale_parameters = list(lower = c(0.01, 0, 0), upper = c(10, 1, 365.25)),
                                            cohort_step_size = 10,
                                            dt = 0.25)


## -----------------------------------------------------------------------------

parameters_fit <- get_parameters(overrides = list("b0" = fitted_values$par[1], 
                                                  "b1" = fitted_values$par[2], 
                                                  "phi" = fitted_values$par[3]),
                                 contact_population_list = contact_population_list)

model_simulation_fit <- run_model(parameters_fit, 
                                  warm_up = 365 * 9, 
                                  max_t = 365 * 10)

# visualising the simulated data and fitted model
ggplot(data = df |> 
         dplyr::mutate(age_chr = 
                         factor(df$age_chr, 
                                levels = c("[0,1)", "[1,2)", "[2,3)", "[3,4)", "[4,5)", "[5,10)", "[10,30)", "[30,50)", "[50,70)", "[70,90)"))), 
       aes(x = time, y = incidence, col = age_chr)) + 
  geom_point() + 
  theme_classic() +
  geom_line(data = model_simulation_fit, inherit.aes = FALSE, 
            aes(x = time, y = Incidence, group = age_chr, col = age_chr)
            )



