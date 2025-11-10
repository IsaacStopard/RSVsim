# Function to run an Approximate Bayesian Computation Sequential Monte Carlo (ABC-SMC) algorithm

—– NOT FINISHED —— implementation of an ABC-SMC algorithm. This function
will not work when fitting the initial conditions.

## Usage

``` r
RSVsim_ABC_SMC(
  target,
  epsilon_matrix,
  summary_fun,
  dist_fun,
  prior_fun,
  n_param_attempts_per_accept,
  used_seed_matrix,
  prior_dens_fun,
  particle_low,
  particle_up,
  nparticles,
  ncores = 1,
  fitted_parameter_names,
  fixed_parameter_list,
  times = seq(0, 365.25 * 5, 0.25),
  cohort_step_size = 1/12 * 365.25,
  warm_up = 365.25 * 4
)
```

## Arguments

- target:

  Values to fit to.

- epsilon_matrix:

  Matrix of tolerance values. Different columns correspond to the values
  for different data points and different rows correspond to the values
  for the different generations.

- summary_fun:

  Function to calculate the summary statistics equivalent to the target
  values.

- dist_fun:

  Function to the calculate the error between the target and
  `summary_fun` outputs.

- prior_fun:

  Function to sample from the priors for all parameters. Must return a
  vector.

- n_param_attempts_per_accept:

  Number of samples to try for each accepted particle.

- used_seed_matrix:

  Matrix of seeds: number of rows must be equal to the number of
  generations and number of columns must be equal to nparticles.

- prior_dens_fun:

  Function that calculates the probability density of all parameters
  given the prior distributions. The joint probability is the product of
  the values returned by this function.

- particle_low:

  Lower bounds on the parameters.

- particle_up:

  Upper bounds on the parameters.

- nparticles:

  Integer. Number of samples from the approximate posterior.

- ncores:

  Number of cores. If greater than one then it is run in parallel.

- fitted_parameter_names:

  Vector of names of the parameters that are being estimated.

- fixed_parameter_list:

  List of parameter values to run the model excluding the fitted
  parameters.

- times:

  Simulation times. Default: 0 to 365.25 days with intervals of 0.25
  days.

- cohort_step_size:

  Time steps to run the model over before adjusting the ages of all
  cohorts. Default: 1 month. If `is.numeric(cohort_step_size) == FALSE`
  then cohort ageing is not applied. Can have a maximum of 3 decimal
  places.

- warm_up:

  Length of time-points to exclude before calculating the likelihood.
  Default: `NULL`.

## Value

List of fixed parameters, max_t and warm_up.
