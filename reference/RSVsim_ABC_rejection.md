# Function to run an Approximate Bayesian Computation (ABC) rejection algorithm

Runs the ABC-rejection algorithm. This will also work with fitting the
initial conditions.

## Usage

``` r
RSVsim_ABC_rejection(
  target,
  epsilon,
  summary_fun,
  dist_fun,
  prior_fun,
  n_prior_attempts,
  nparticles,
  used_seeds_all,
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

- epsilon:

  Acceptable error for each target.

- summary_fun:

  Function to calculate the summary statistics equivalent to the target
  values.

- dist_fun:

  Function to the calculate the error between the target and
  `summary_fun` outputs.

- prior_fun:

  Function to sample from the priors for all parameters. Must return a
  vector.

- n_prior_attempts:

  Number of random samples from the prior to attempt for each accepted
  particle.

- nparticles:

  Integer. Number of samples from the approximate posterior.

- used_seeds_all:

  Vector. Seeds used when generating the prior samples for each accepted
  particle.

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
