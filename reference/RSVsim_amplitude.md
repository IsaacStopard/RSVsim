# Summary function to calculate the age-specific amplitude of changes in incidence over the simulation.

Function calculates the difference between maximum and minimum incidence
stratified by age over the whole simulation (excluding warm up).

## Usage

``` r
RSVsim_amplitude(out)
```

## Arguments

- out:

  `RSVsim_run_model` function output.

## Value

Summary statistics (amplitude by age-group) from the simulation output
(dataframe).
