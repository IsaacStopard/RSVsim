# Calibration

``` r
library(RSVsim)
```

In this vignette, we give an overview of how to use Approximate Bayesian
Computation (ABC) rejection and Sequential Monte Carlo ABC algorithm
(ABC-SMC) algorithms to calibrate the mathematical model of RSV
transmission to simulated data. The code to run these samplers was
adapted from code provided in Minter, Amanda, and Renata Retkute.
“Approximate Bayesian Computation for infectious disease modelling.”
Epidemics 29 (2019): 100368
<https://doi.org/10.1016/j.epidem.2019.100368>.

First, we run the model to generate some data to fit to.

``` r

# specify age categories to run the model with
age_limits <- c(seq(0, 2, 0.2), seq(10, 60, 5))

# get contact matrix
contact_population_list <- RSVsim_contact_matrix(country = "United Kingdom", age_limits = age_limits)

# get the default model parameters with a b0 value of 0.08, b1 value of 0, phi value of 0 and change the initial conditions 
overrides <- list("b0" = 0.11, "b1" = 0.3, "phi" = 10)

parameters <- fixed_parameter_list <- RSVsim_parameters(overrides = overrides, 
                                                        contact_population_list = contact_population_list)

out <- RSVsim_run_model(parameters = parameters,
                        times = seq(0, 365 * 2, 0.25), # maximum time to run the model for
                        cohort_step_size = 0.2 * 365, # time at which to age people\
                        warm_up = 365)

fitted_parameter_names <- c("b0", "b1", "phi")

nAges <- parameters[["nAges"]]
```

Next, we define a function to calculate some summary statistics from the
model output. Specifically, we calculate the age-specific peak time of
incidence infections, the amplitude of incidence (difference between
maximum and minimum incidence) and the total yearly incidence. RSVsim
has built in functions to calculate these metrics from the model output.
We define a function `summary_fun` to combine these metrics. We will use
these as the target for the ABC-rejection algorithm. We must also
specify a function that randomly samples from the prior distributions of
all fitted parameters and returns a vector of the parameters. We use
latin hypercube sampling to efficiently sample the prior space. In this
example we fit the `b0`, `b1` and `phi` parameters. A vector of the
parameter names that are being fitted in the same order as the prior
samples must also be provided to set up the parameters to run the model.

``` r

#########################
##### ABC functions #####
#########################

# calculate the summary statistics
summary_fun <- function(out){
  return(c(RSVsim_total_incidence(out), RSVsim_amplitude(out), RSVsim_peak(out)))
}

target <- summary_fun(out)

# a function that samples from the priors is required
# we use latin hypercube sampling to efficiently sample from the prior distributions
prior_fun <- function(n_prior_attempts){
  
  x <- lhs::randomLHS(n_prior_attempts, 3)
  
  # adjusting the prior distributions
  x[,1] <- qunif(x[,1], min = 0.05, max = 0.25)
  x[,2] <- qunif(x[,2], min = 0, max = 1) # not necessary but added for completeness
  x[,3] <- qunif(x[,3], min = -365.25/2, max = 365.25/2)
  
  return(as.matrix(x, nrow = n_prior_attempts))
}

# a function that calculates the distance between the summary statistics is required
dist_fun <- function(target, target_star, n = nAges){
  return(
    c(
      RSVsim_abs_dist_fun(target[1:n], target_star[1:n]),
      RSVsim_abs_dist_fun(target[(n+1):(n*2)], target_star[(n+1):(n*2)]),
      RSVsim_shortest_periodic_dist_fun(target[(2*n+1):(n*3)], target_star[(2*n+1):(n*3)], period = 365.25)
    )
  )

}
```

To implement the ABC-rejection algorithm we must specify the tolerance.
We specify a specific tolerance for each metric. To calculate the
tolerance values that means approximately 0.1% of all simulations are
accepted we run the model 1000 times with different parameter
combinations drawn from the priors. We then calculate the error
(distances) between the summary statistics from each simulation and the
target summary statistics. We set the a range of tolerances using
different percentiles of the distances and check the number of these
simulations that are accepted for each tolerance, and use the smallest
tolerance with at least 1 simulation accepted. **This method of
calculating suitable tolerance values is optional and custom values can
be specified by the user.**

``` r
#################################
##### setting the tolerance #####
#################################

# this code is optional but is used to approximate epsilon values that will give us an acceptance rate of 0.1%

# calculating a tolerances that means at least 1 particle combination is accepted every 100 simulations
# getting 1000 samples from the priors
set.seed(123)
n_check <- 1000
prior_params <- prior_fun(n_check)

# simulating the summary statistics for each particle combination
prior_distances <- sapply(1:n_check, function(i){
  parameters_in <- RSVsim_update_parameters(fixed_parameter_list, fitted_parameter_names, prior_params[i, ])
  
  out <- RSVsim_run_model(parameters = parameters_in,
                   times = seq(0, 365*2, 0.25), # maximum time to run the model for
                   cohort_step_size = 0.2*365, # time at which to age people
                   warm_up = 365)
  
  return(dist_fun(target, 
                  summary_fun(out)
                  )
         )
})

# calculating the number of particles for which all summary statistics are within the tolerance given different percentiles of the summary statistics for all the particles 
# selecting the tolerance with at least 1 model simulation using the prior samples that is within the tolerance 
nsuccess <- rep(NA, n_check)
q <- seq(0.01, 1, 0.01)

for(i in 1:100){
  epsilon_check <- round(apply(prior_distances, 1, quantile, probs = c(q[i])), digits = 2)
  nsuccess[i] <- sum(sapply(1:100, function(j){all(prior_distances[,j] <= epsilon_check)}))
}

# selecting a tolerance where at least 1 particle is accepted for the 1000 simulations
epsilon <- round(apply(prior_distances, 1, quantile, probs = c(q[min(which(nsuccess >= 1))])), digits = 2)
```

To run the model we specify the number of particles to fit, and the
seeds for each particle. The number of cores can also be specified - if
greater than one then parallel processing is implemented using a
parallel socket cluster with the parallel R package.

``` r
###############################################
##### running the ABC-rejection algorithm #####
###############################################

nparticles = 500

used_seeds_all <- seq(1, nparticles)

ncores <- 8

fit <- RSVsim_ABC_rejection(target = target,
                            epsilon = epsilon,
                            summary_fun = summary_fun,
                            dist_fun = dist_fun,
                            prior_fun = prior_fun,
                            n_prior_attempts = 10000,
                            nparticles = nparticles,
                            used_seeds_all = used_seeds_all,
                            ncores = ncores,
                            fitted_parameter_names = fitted_parameter_names,
                            fixed_parameter_list = fixed_parameter_list,
                            times = seq(0, 365, 0.25), # maximum time to run the model for
                            cohort_step_size = 0.2*365, # time at which to age people\
                            warm_up = NULL)

# to visualise the posterior distributions
fit_b0 <- unlist(fit[,"b0"])
fit_b1 <- unlist(fit[,"b1"])
fit_phi <- unlist(fit[,"phi"])

hist(fit_b0)
hist(fit_b1)
hist(fit_phi)
```

We also provide a function to fit the model using ABC-SMC. The ABC-SMC
algorithm requires a function to return the prior density and a
sequential reduction in the tolerances.

``` r

# function that returns the prior likelihood of each parameter (prior density function)
# returns a vector of the same length as the number of parameters
prior_dens_fun <- function(x){
  
  # adjusting the prior distributions
  return(c(dunif(x[1], min = 0.05, max = 0.25),
           dunif(x[2], min = 0, max = 1),
           dunif(x[3], min = -365.25/2, max = 365.25/2)
           )
         )
}
```

Similar to the ABC-rejection algorithm, we use the previously simulated
summary statistics to calculate the tolerances corresponding to a given
number of particle combinations that are accepted. We do this for a
decreasing acceptance rate; the tolerances are stored in a matrix with
each column corresponding to each summary statistic and the rows
corresponding to each tolerance level (the lowest tolerance in the last
row) (`G`; generation). **Again, this method of calculating suitable
tolerance values is optional and custom values can be specified by the
user.**

``` r
nsuccess <- rep(NA, n_check)
q <- seq(1/n_check, 1, 1/n_check)
for(i in 1:n_check){
  epsilon_check <- round(apply(prior_distances, 1, quantile, probs = c(q[i])), digits = 5)
  nsuccess[i] <- sum(sapply(1:n_check, function(j){all(prior_distances[,j] <= epsilon_check)}))
}

acceptance_rate <- sort(c(1, 10, 50, 75), decreasing = TRUE)

q_percentile <- vector(mode = "list", length = length(acceptance_rate))

q_percentile <- lapply(acceptance_rate, function(ar){q[min(which(nsuccess/n_check * 100 > ar))]}) # 1

q_percentile_ar1 <- q_percentile[[which(acceptance_rate == 1)]]

# selecting a tolerance where at least 1 particle is accepted for the 1000 simulations
epsilon_matrix <- as.matrix(do.call(rbind, lapply(q_percentile, function(q){round(apply(prior_distances, 1, quantile, probs = q), digits = 3)})))

G <- nrow(epsilon_matrix)
```

We can then run the ABC-SMC using the code below. The function also
requires: (1) a matrix of random seeds for each generation and accepted
particle (`used_seed_matrix`), (2) the minimum and maximum values for
each fitted parameter used in the prior distributions (`particle_low`
and `particle_up` respectively), (3) the number of accepted particles
for each generation (`nparticles`), (4) the number of prior samples to
try for each accepted particle combination (`n_prior_attempts`), (5) the
names of the fitted parameters (`fitted_parameter_names`, these must be
characters and can include indexing into parameters that are stored as
vectors or matrices), (6) a list of all parameters equivalent to those
used in `RSVsim_run_model` `parameters` (`fixed_parameter_list`) and (7)
the times and cohort step sizes used in `RSVsim_run_model`.

``` r
used_seed_matrix <- matrix(seq(1, nparticles * G), nrow = G)

fit_smc <- RSVsim_ABC_SMC(target = target,
                          epsilon_matrix = epsilon_matrix,
                          summary_fun = summary_fun,
                          dist_fun = dist_fun,
                          prior_fun = prior_fun,
                          n_prior_attempts = 10000,
                          used_seed_matrix = used_seed_matrix,
                          prior_dens_fun = prior_dens_fun,
                          particle_low = c(0.05, 0, -365.25/2),
                          particle_up = c(0.25, 1, 365.25/2),
                          nparticles = nparticles,
                          ncores = ncores,
                          fitted_parameter_names = fitted_parameter_names,
                          fixed_parameter_list = fixed_parameter_list,
                          times = seq(0, 365, 0.25), # maximum time to run the model for
                          cohort_step_size = 0.2*365, # time at which to age people
                          warm_up = NULL)
```

**Running the model with the fitted parameters.** The ABC algorithms
return matrices of the fitted parameter values. For this is the main
output and for a list of the fitted parameter values is returned for
each generation. The fitted parameter columns correspond to the
parameter and the rows correspond to each accepted parameter
combination. To run the model with one particle (combination of fitted
parameters) we can therefore select the correct matrix of fitted
parameters and use the to adjust the parameters. For example, if I want
to run the *first* accepted particle:

``` r

# ABC rejection
parameters_ABC_rejection <- RSVsim_update_parameters(parameters, fitted_parameter_names, c(fit[1, c("b0", "b1", "phi")]))

# ABC SMC
# for the ABC SMC algorithm I must select the fitted parameters from the last generation

fit_params_SMC <- fit_smc$fitted_parameters[[G]]
colnames(fit_params_SMC) <- fitted_parameter_names

parameters_ABC_SMC <- RSVsim_update_parameters(parameters, fitted_parameter_names, c(fit_params_SMC[1, c("b0", "b1", "phi")]))
```
