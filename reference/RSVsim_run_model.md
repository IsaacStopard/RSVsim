# Runs the RSV model using odin and dust

Function to run the transmission model with cohort aging. Prevalence,
incidence between the given time-steps and the incidence rate per day
are also calculated.

## Usage

``` r
RSVsim_run_model(
  parameters,
  times = seq(0, 365.25 * 1, 0.25),
  cohort_step_size = 1/12 * 365.25,
  warm_up = NULL
)
```

## Arguments

- parameters:

  List of parameters from `RSVsim_parameters` function.

- times:

  Simulation times. Default: 0 to 365.25 days with intervals of 0.25
  days. The minimum interval size is 1E-5.

- cohort_step_size:

  Time steps to run the model over before adjusting the ages of all
  cohorts. Default: 1 month. If `is.numeric(cohort_step_size) == FALSE`
  then cohort ageing is not applied. Can have a maximum of 7 decimal
  places.

- warm_up:

  Length of time-points to exclude before calculating the likelihood.
  Default: `NULL`.

## Value

Simulation output (dataframe). In the dataframe, age refers to the
lowest age in the age group.
