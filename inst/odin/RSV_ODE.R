# Time Steps
nAges <- parameter()

## Model Parameters
b0 <- parameter()
b1 <- parameter()
phi <- parameter()
delta <- parameter()
gamma_s <- parameter()
gamma_p <- parameter()
nu <- parameter()
omega_vect <- parameter()
prop_detected_vect <- parameter()
sigma_vect <- parameter()
alpha_vect <- parameter() # reduced susceptibility for infection in age group i
matrix_per_person <- parameter()

dim(matrix_per_person) <- c(nAges, nAges)
dim(sigma_vect) <- nAges
dim(omega_vect) <- nAges
dim(prop_detected_vect) <- nAges
dim(alpha_vect) <- nAges

## Derivatives for Flows Between Compartments
##------------------------------------------------------------------------------

# prevalence of infection in age group j multiplied by reduced infectiousness in age group j (due to immunity early in life)
temp[] <- omega_vect[i] * (Is[i] + Ip[i]) / N[i]
# contacts multiplied by prevalence
s_ij[,] <- matrix_per_person[i,j] * temp[j]
lambda[] <- b0 * (1 - b1 * cos(2 * 3.14159265358979323846 / 365.25 * (time - phi))) * sum(s_ij[i,])

incidence_rate_p[] <- lambda[i] * sigma_vect[i] * Sp[i]
incidence_rate_s[] <- lambda[i] * sigma_vect[i] * alpha_vect[i] * Ss[i]

# primary infection
deriv(Sp[1:nAges]) <- -incidence_rate_p[i]
deriv(Ep[1:nAges]) <- incidence_rate_p[i] - delta * Ep[i]
deriv(Ip[1:nAges]) <- delta * Ep[i] - gamma_p * Ip[i]
# secondary infection
deriv(Ss[1:nAges]) <- -incidence_rate_s[i] + nu * R[i]
deriv(Es[1:nAges]) <- incidence_rate_s[i] - delta * Es[i]
deriv(Is[1:nAges]) <- delta * Es[i] - gamma_s * Is[i]

deriv(R[1:nAges]) <- gamma_p * Ip[i] + gamma_s * Is[i] - nu * R[i]

# used to tracks the incidence in a given age-group (the incidence of the age when infected)
deriv(Incidence[1:nAges]) <- incidence_rate_s[i] + incidence_rate_p[i]

N[1:nAges] <- Ss[i] + Es[i] + Is[i] + Sp[i] + Ep[i] + Ip[i] + R[i]

output(DetIncidence[1:nAges]) <- prop_detected_vect[i] * Incidence[i]
output(Incidence_rate[1:nAges]) <- incidence_rate_s[i] + incidence_rate_p[i]
output(DetIncidence_rate[1:nAges]) <- prop_detected_vect[i] * Incidence_rate[i]
output(prev[1:nAges]) <- (Ip[i] + Is[i])/N[i]
output(prev_p[1:nAges]) <- Ip[i]/N[i]
output(prev_s[1:nAges]) <- Is[i]/N[i]

## Initial states:
initial(Sp[1:nAges]) <- Sp0[i]
initial(Ep[1:nAges]) <- Ep0[i]
initial(Ip[1:nAges]) <- Ip0[i]
initial(Ss[1:nAges]) <- Ss0[i]
initial(Es[1:nAges]) <- Es0[i]
initial(Is[1:nAges]) <- Is0[i]
initial(R[1:nAges]) <- R0[i]
initial(Incidence[1:nAges]) <- Incidence0[i]

##Initial vectors
Sp0 <- parameter()
Ep0 <- parameter()
Ip0 <- parameter()
Ss0 <- parameter()
Es0 <- parameter()
Is0 <- parameter()
R0 <- parameter()
Incidence0 <- parameter()

##Dimensions of the different "vectors" used
# For the State Variables
dim(Sp) <- nAges
dim(Ep) <- nAges
dim(Ip) <- nAges
dim(R) <- nAges
dim(Ss) <- nAges
dim(Es) <- nAges
dim(Is) <- nAges
dim(Incidence) <- nAges

dim(N) <- nAges
dim(DetIncidence) <- nAges
dim(Incidence_rate) <- nAges
dim(DetIncidence_rate) <- nAges
dim(prev) <- nAges
dim(prev_s) <- nAges
dim(prev_p) <- nAges
dim(incidence_rate_s) <- nAges
dim(incidence_rate_p) <- nAges
dim(Sp0) <- nAges
dim(Ep0) <- nAges
dim(Ip0) <- nAges
dim(Ss0) <- nAges
dim(Es0) <- nAges
dim(Is0) <- nAges
dim(R0) <- nAges
dim(Incidence0) <- nAges

dim(lambda) <- nAges
dim(s_ij) <- c(nAges,nAges)
dim(temp) <- nAges

## Calculate Vaccination Rate
##------------------------------------------------------------------------------

# Interpolation of vaccination rate over time
# mv is vaccines per day
# tt_vaccine is the timepoints at which vacc rollouts commence
# max_vaccine is the doses delivered at those timepoints
# max_doses is the maximum number of doses distributed on any given day
# vacc_coverage is the vector of coverages of dimension nAges (input by the user)
# target_pop is the simulation size (input by the user)
# vaccine_period is the period over which vaccines are distributed in a single campaign (e.g. single year) (input by the user)

vacc_coverage[] <- user() # dimension nAges
target_pop[] <- user() # integer
vaccine_period[] <- user() # integer
vacc_t[] <- user() # vector, can be any length specified by user

max_doses <- max(vacc_coverage) * target_pop / vaccine_period
vacc_t <- c(1, 360, 720) # example so we can see how it works - can remove

tt_vaccine <- 0
max_vaccine <- 0


for (i in 1:length(vacc_t)){
  tt_vaccine <- c(tt_vaccine, vacc_t[i], vacc_t[i]+vaccine_period)
  max_vaccine <- c(max_vaccine, max_doses, 0)
}

# creates the interpolation function (f you want to test outside odin, try cinterpolate::interpolation_function())
mv <- interpolate(tt_vaccine, max_vaccine, "constant")

# Number of people available for vaccination at time point
vr_temp[] <- Sp[i,1] * vacc_coverage[i] + Ep[i,1] * vacc_coverage[i] + Ss[i,1] * vacc_coverage[i] + Es[i,1] * vacc_coverage[i] + R[i,1] * vacc_coverage[i]
dim(vr_temp) <- nAges

# Catch so vaccination rate does not exceed 1 if the number of people available for vaccination < rate
vr_den <- if(sum(vr_temp) < mv) mv else sum(vr_temp)
vr <- if(mv==0) 0 else mv / vr_den # Vaccination rate to achieve capacity given number in vaccine-eligible population

