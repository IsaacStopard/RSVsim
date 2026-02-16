# Helper function to update the parameter list by names

Calculates the absolute distance between two times when the times are
circular. Vectorised.

## Usage

``` r
RSVsim_update_parameters(
  fixed_parameter_list,
  fitted_parameter_names,
  fitted_parameter_values
)
```

## Arguments

- fixed_parameter_list:

  List of parameters, should be the output of the `RSVsim_parameters`
  function.

- fitted_parameter_names:

  Vector of parameter names to be updated. Can include indexing such as
  "parameter_name\[i\]".

- fitted_parameter_values:

  Vector of updated parameter values.

## Value

Parameter list.
