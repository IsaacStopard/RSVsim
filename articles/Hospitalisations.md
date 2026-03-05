# Hospitalisations

``` r
library(RSVsim)
```

We do not explicitly include hospitalisation in the ODE model, but can
estimate the hospitalisations from the simulated incidence. Note the
time until hospitalisation is therefore ignored.

## Probability of hospitalisation given infection

The probability of hospitalisation given unvaccinated infection for age
group $i$ is calculated as
$p\left( hosp|inf,i \right) = \frac{H_{i}}{I_{i}}$, where $H_{i}$ is the
number of hospitalisations in age group $i$ and $I_{i}$ is the
corresponding incidence. The age-specific hospitalisation rates per
person is $h_{i} = \frac{H_{i}}{N_{i}}$, where $N_{i}$ is the population
size, meaning $p\left( hosp|inf \right) = \frac{h_{i}}{I_{i}/N_{i}}$.
Similarly, we can calculate $p\left( hosp|inf \right)$ in unvaccinated
people assuming a given vaccine efficacy against hospitalisation for
each vaccinated state ($VE_{v}$) and hospitalisation rates, thus,

$$p\left( hosp|inf,i \right) = \frac{h_{i}}{\frac{I_{1i}}{N_{1i}} + \sum\limits_{2}^{V}\frac{I_{vi}}{N_{vi}}\left( 1 - VE_{vi} \right)}$$

In reality, hospitalisation rates are typically sampled, but the total
incidence is unknown. We can therefore use to estimate the incidence in
the scenario corresponding to setting in which the age specific
hospitalisation rates were estimated. In this example we assume no
existing vaccination. The default vaccination states are 1
(unvaccinated) and 2 (vaccinated). First, we use to approximate the
incidence.

``` r

age_limits <- c(seq(0, 5, 1/12), seq(10, 75, 5))

contact_population_list <- RSVsim_contact_matrix(age_limits = age_limits)

parameters <- RSVsim_parameters(contact_population_list = contact_population_list)

sim_status_quo <- RSVsim_run_model(parameters = parameters,
                                   times = seq(0, 365.25*11, 0.25),
                                   cohort_step_size = 1/12 * 365.25,
                                   warm_up = 365.25 * 10)
```

We can then use the function to estimate $p\left( hosp|inf \right)$ for
given values of $h_{i}$ (note these must be per person). We must also
provide the times over which the hospitalisations are estimated and the
vaccine efficacy against hospitalisation for each vaccinated state.

``` r

# age specific number of hospitalisations per person

hosp_rate <- c(5338.91375, 1751.145228, 715.9402014, 208.4621667, 10.86702733, 2.58129891, 2.321354938, 1.868628613,
                 1.320440418, 1.834551244, 1.883359955, 2.861263391, 4.744902927, 6.645799647, 9.348531561, 14.3538506,
                 22.62406397, 32.94296492, 87.63366855) / 100000

age_limits_hosp_rate <- c(0, 0.5, 1, 2, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75) # lower age limits of the hospitalisation rate data

VE_HOSP <- matrix(rep(0.1, length(age_limits_hosp_rate)), ncol = 1)

prob_hosp_given_inf_list <- RSVsim_prob_hosp_given_inf(hosp_rate = hosp_rate,
                                                       age_limits_hosp_rate = age_limits_hosp_rate,
                                                       hosp_rate_min_time = 0,
                                                       hosp_rate_max_time = 365.25,
                                                       sim_status_quo = sim_status_quo,
                                                       VE_HOSP = VE_HOSP,
                                                       VE_HOSP_vacc_states = c(2) # vector identifying the vaccination states that are vaccinated (position corresponds to the column position in VE_HOSP)
                                                       )
```

## Hospitalisations

The function returns a list of $p\left( hosp|inf,i \right)$ and the
model outputs with hospitalisations. We can then use these probabilities
to calculate the hospitalisations from the simulated incidence in
different scenarios.

``` r

# single 30 day vaccination campaign in the oldest age group
vaccine_times <- c(0, c(90, 120) + 365.25 * 10)
vaccine_period <- diff(c(vaccine_times, max(vaccine_times)))
nVaccTimes <- length(vaccine_times)

nAges <- length(contact_population_list$age_limits)

# vaccination campaign target coverage: 80%
vaccine_cov <- rbind(matrix(rep(0, (nAges-1) * nVaccTimes),
                            nrow = (nAges-1), ncol = nVaccTimes),
                     c(0, 0.70, 0))

# not necessary due to default settings but included for completeness
nVaccStates <- 2
VE <- matrix(rep(0.85, nAges * (nVaccStates-1)), nrow = nAges) # vaccine efficacy against infection
gamma_vaccine <- 1 / (365.25 * 2)

sim_vaccination <- RSVsim_run_model(parameters = RSVsim_parameters(overrides = list("vaccine_times" = vaccine_times,
                                                                                    "vaccine_period" = vaccine_period,
                                                                                    "nVaccTimes" = nVaccTimes,
                                                                                    "vaccine_cov" = vaccine_cov,
                                                                                    "nVaccStates" = nVaccStates,
                                                                                    "VE" = VE,
                                                                                    "gamma_vaccine" = gamma_vaccine),
                                                                   contact_population_list = contact_population_list),
                                    times = seq(0, 365.25*11, 0.25),
                                    cohort_step_size = 1/12 * 365.25, # time at which to age people
                                    warm_up = 365.25 * 10)
```
