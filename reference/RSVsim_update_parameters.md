# Helper function to update the parameter list by names

Changes the values in the list of parameters to run the model. To fit
parameters stored within vectors or matrices, `fitted_parameter_names`
can also include indexing such as "parameter_name\[i\]".

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

  Vector of parameter names to be updated.

- fitted_parameter_values:

  Vector of updated parameter values.

## Value

Parameter list.
