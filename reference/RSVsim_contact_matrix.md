# Create contact matrix for RSVsim

Uses the `socialmixr` package to access the POLYMOD data and calculate
per-person mean contact matrix. If age groupings smaller than 1-year are
required, the function calculates the total contacts in the population
with the nearest integer age groupings (required for `socialmixr`
package) and splits the total contacts into smaller age groupings by
dividing them uniformly.

## Usage

``` r
RSVsim_contact_matrix(
  country = "United Kingdom",
  age_limits = c(seq(0, 5, 2/12), seq(10, 70, 5))
)
```

## Arguments

- country:

  Country for use in the `contact_matrix` function in the `socialmixr`
  package. Can be given as country name or 2 digit ISO code. United
  Kingdom default.

- age_limits:

  Lower limits of the age groups to run the simulation. Ages must be in
  years. Note the smallest age limits allowed are 1E-8.

## Value

List including the contact matrix (`matrix_mean` is the mean per person,
`matrix_contacts` is the total contacts), the age limits, the age
distribution used to calculate the total contacts and the maximum age in
the contact data.
