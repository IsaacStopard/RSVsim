# Probability of hospitalisation given unvaccinated infection

Calculates the age-specific probability of hospitalisation given
infection, using observed hospitalisation rates per person among the
population as a whole and RSVsim to estimate number of infections per
person.

## Usage

``` r
RSVsim_prob_hosp_given_inf(
  hosp_rate,
  age_limits_hosp_rate,
  hosp_rate_min_time,
  hosp_rate_max_time,
  sim_status_quo,
  VE_HOSP,
  VE_HOSP_vacc_states
)
```

## Arguments

- hosp_rate:

  Vector of age-specific RSV hospitalisation rates (per person) in the
  status quo scenario.

- age_limits_hosp_rate:

  Vector of ages corresponding to the first age of each age category for
  the hospitalisation rates and probability of hospitalisation given
  infection. These ages must be present in the age_limits used to run
  the simulations.

- hosp_rate_min_time:

  Minimum simulation time to select when summing the hospitalisations to
  calculate the hospitalisation rates.

- hosp_rate_max_time:

  Maximum simulation time to select when summing the hospitalisations to
  calculate the hospitalisation rates.

- sim_status_quo:

  `RSVsim_run_model` output used to estimate the number of infections
  when calculating the probability of hospitalisation given infection.

- VE_HOSP:

  Matrix of age-specific vaccine efficacy against hospitalisation. Rows
  must correspond to the ages in `age_limits_hosp_rate`, and columns
  must correspond to the vaccinated states used in `RSVsim_run_model`.

- VE_HOSP_vacc_states:

  vector of integers indicating the vaccinated states used in
  `RSVsim_run_model`.

## Value

List of probability of hospitalisation given infection and
`sim_status_quo` model outputs with hospitalisations.
