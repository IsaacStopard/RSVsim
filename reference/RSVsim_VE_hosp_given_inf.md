# Vaccine efficacy against hospitalisation conditional on infection

Calculates the vaccine efficacy against hospitalisation given infection
using the overall vaccine efficacy against hospitalisation and vaccine
efficacy against infection. If the number ages in VE_inf are greater
than the number of ages in VE_hosp we just select the

## Usage

``` r
RSVsim_VE_hosp_given_inf(VE_hosp, VE_inf)
```

## Arguments

- VE_hosp:

  Vaccine efficacy against hospitalisation.

- VE_inf:

  Vaccine efficacy against infection.

## Value

Vaccine efficacy against hospitalisation given infection
(`VE_hosp_given_inf`).
