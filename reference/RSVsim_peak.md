# Summary function to calculate the age-specific timing of peak incidence.

Function calculates the time when incidence peaks stratified by age over
the whole simulation (excluding warm up).

## Usage

``` r
RSVsim_peak(out)
```

## Arguments

- out:

  `RSVsim_run_model` function output.

## Value

Summary statistics (peak incidence time ordered by age-group) from the
simulation output (dataframe).
