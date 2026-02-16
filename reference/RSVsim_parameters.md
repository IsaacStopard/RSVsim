# Get parameters for RSVsim

Provides the default parameters to run the RSV model.

## Usage

``` r
RSVsim_parameters(overrides = list(), contact_population_list)
```

## Arguments

- overrides:

  List of default parameters to change.

- contact_population_list:

  List of outputs from the `create_contact_matrix` function.

## Value

Parameter list. Default values:  
`b0`: 0.08. Transmission rate coefficient.  
`b1`: 0. Amplitude of seasonal forcing in the transmission rate.  
`phi`: 0. Phase shift of the seasonal forcing.  
`delta`: 1/4. Inverse of the latent period.  
`gamma_s`: 1/8. Inverse of the infectious period of subsequent
(secondary) infections.  
`gamma_p`: 1/9. Inverse of the infectious period of first (primary)
infections.  
`nu`: 1/200. Inverse of the duration of natural immunity.  
`omega_vect`: 0.35 if greater than 5 years old else 1. Reduced
infectiousness due to immunity in early life.  
`prop_detected_vect`: age-specific proportion detected.  
`sigma_vect`: 0.7 if 1 months old, 0.8 if 2 months old, 0.9 if 3 months
old, else 1. Reduced susceptibility due to immunity in early life.  
`alpha_vect` = 0.4 if less than 10 years old else 0.3. Reduced
susceptibility to secondary infections in each age group.  
`nAges`: must be obtained from `contact_population_list`. Number of age
groups - same length as age.limits.  
`age.limits`: must be obtained from `contact_population_list`. Lower age
limits in years.  
`matrix_per_person`: must be obtained from `contact_population_list`.
Contact matrix of mean contacts per person.  
`max_age`: must be obtained from `contact_population_list`. Maximum age
in years.  
`age_chr`: Age limits as a character vector.  
`size_cohorts`: the differences in ages in days - used to calculate the
size of the cohorts.  
`rel_sizes`: the relative size of each cohort (age category duration in
days divided by the total number of days).  
`total_population`: 1861923. Total population size.  
`nVaccStates`: 2. Number of vaccination states. Default: 1:
unvaccinated, 2: vaccinated. Minimum is 2.  
`gamma_vaccine`: 1 / (365\*2). Rate of waning between vaccinated states.
Size of vector must be nVaccStates minus one.  
`nVaccTimes`: 5. Length of start times of vaccination distributions.  
`vaccine_times`: c(0, 365.25, 365.25 + 30, 730.5, 730.5 + 30). Times of
vaccination distributions. First time should be 0.  
`vaccine_period`: rep(30, 5). Duration of vaccination distributions.
`vaccine_cov`: matrix with the rows corresponding to the vaccine_times
and the columns corresponding to the ages. Age-specific proportion of
unvaccinated people to have been vaccinated by the end of the
vaccination distribution for each `vaccine_time`, assuming no waning of
vaccination. The vaccination rate is calculated as
`-log(1 - vaccine_cov) / vaccine_period`. Default is 0 coverage for all
ages. Changes in effective coverage with are given as a model output.  
`VE`: 0.85. Vaccine efficacy for each vaccinated state. Length should be
equal to nVaccStates - 1.
`Sp0, Ss0, Ep0, Es0, Ip0, Is0, R0, Incidence0`: initial conditions to
run the model for each compartment - these are given as prevalence and
the initial conditions calculated in this function. List. Default:
`NULL`. If `NULL`: 0.1% RSV prevalence is assumed for people during the
primary infection, which is seeded at the beginning of the simulation.
All other people are assumed to be susceptible to their primary
infection.
