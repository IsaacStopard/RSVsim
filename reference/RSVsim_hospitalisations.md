# Calculate number of hospitalisations for `RSVsim_run_model` output

Function to estimate the number of hospitalisations from the modelled
incidence and estimated probability of hospitalisation given
unvaccinated infection.

## Usage

``` r
RSVsim_hospitalisations(
  prob_hosp_given_inf,
  age_limits_hosp_rate,
  sim_scenario,
  VE_HOSP,
  VE_HOSP_vacc_states
)
```

## Arguments

- prob_hosp_given_inf:

  Age-specific probability of hospitalisation given unvaccinated
  infection. Ages must correspond to those in `age_limits_hosp_rate`.

- age_limits_hosp_rate:

  Vector of ages corresponding to the first age of each age category for
  the hospitalisation rates and probability of hospitalisation given
  infection. These ages must be present in the age_limits used to run
  the simulations.

- sim_scenario:

  `RSVsim_run_model` output. Hospitalisations are calculated for this
  scenario using the estimated probability of hospitalisation given
  infection.

- VE_HOSP:

  Matrix of age-specific vaccine efficacy against hospitalisation. Rows
  must correspond to the ages in `age_limits_hosp_rate`, and columns
  must correspond to the vaccinated states used in `RSVsim_run_model`.

- VE_HOSP_vacc_states:

  vector of integers indicating the vaccinated states used in
  `RSVsim_run_model`.

## Value

sim_scenario dataframe with hospitalisations column
