# number of age groups
nAges <- parameter()

# number of vaccination groups including unvaccinated
nVaccStates <- parameter()

## Model Parameters
b0 <- parameter()
b1 <- parameter()
phi <- parameter()
delta <- parameter()
gamma_s <- parameter()
gamma_p <- parameter()
gamma_vaccine <- parameter()
nu <- parameter()
omega_vect <- parameter() # reduced infectiousness due to immunity in early life
prop_detected_vect <- parameter()
sigma_vect <- parameter() # reduced susceptibility due to immunity in early life
alpha_vect <- parameter() # reduced susceptibility to secondary infections in age group i
matrix_per_person <- parameter()

dim(matrix_per_person) <- c(nAges, nAges)

# allowing these parameters to vary by both age and vaccination status
dim(sigma_vect) <- c(nAges)
dim(omega_vect) <- c(nAges)
dim(alpha_vect) <- c(nAges)

dim(prop_detected_vect) <- c(nAges)

dim(gamma_vaccine) <- nVaccStates - 1

## Derivatives for Flows Between Compartments
##------------------------------------------------------------------------------

# FOI
# number of infections in age group i and vaccination status j multiplied by reduced infectiousness (due to immunity early in life and vaccination)

N[1:nAges, 1:nVaccStates] <- Ss[i,j] + Es[i,j] + Is[i,j] + Sp[i,j] + Ep[i,j] + Ip[i,j] + R[i,j]

temp_rel[,] <- omega_vect[i] * (Is[i,j] + Ip[i,j])

dim(temp_rel) <- c(nAges, nVaccStates)

temp[] <- sum(temp_rel[i,]) / N_age[i] # sum over all the vaccination states, divide by population size to get the prevalence
dim(temp) <- c(nAges)

# contacts multiplied by prevalence of infectious people
s_ij[,] <- matrix_per_person[i,j] * temp[j]
dim(s_ij) <- c(nAges, nAges)

# age-specific force of infection
lambda[] <- b0 * (1 - b1 * cos(2 * 3.14159265358979323846 / 365.25 * (time - phi))) * sum(s_ij[i,])
dim(lambda) <- c(nAges)

incidence_rate_p[1:nAges,1] <- lambda[i] * sigma_vect[i] * Sp[i,1]
incidence_rate_p[1:nAges,2:nVaccStates] <- lambda[i] * sigma_vect[i] * (1 - VE[i,(j-1)]) * Sp[i,j]

incidence_rate_s[1:nAges,1] <- lambda[i] * sigma_vect[i] * alpha_vect[i] * Ss[i,1]
incidence_rate_s[1:nAges,2:nVaccStates] <- lambda[i] * sigma_vect[i] * (1 - VE[i,(j-1)]) * alpha_vect[i] * Ss[i,j]

VE <- parameter()
dim(VE) <- c(nAges, nVaccStates - 1)

# vr should be by age

# primary infection
deriv(Sp[1:nAges, 1]) <- -incidence_rate_p[i,1] - (VaccState_rate[i,1] * Sp[i,1]) + (VaccState_rate[i, nVaccStates] * Sp[i,nVaccStates])
deriv(Sp[1:nAges, 2:nVaccStates]) <- -incidence_rate_p[i,j] - (VaccState_rate[i,j] * Sp[i,j]) + (VaccState_rate[i,(j-1)] * Sp[i,(j-1)])

deriv(Ep[1:nAges, 1]) <- incidence_rate_p[i,1] - (delta + VaccState_rate[i, 1]) * Ep[i,1] + VaccState_rate[i, nVaccStates] * Ep[i,nVaccStates]
deriv(Ep[1:nAges, 2:nVaccStates]) <- incidence_rate_p[i,j] - (delta + VaccState_rate[i,j]) * Ep[i,j] + (VaccState_rate[i,(j-1)] * Ep[i,(j-1)])

# no vaccination of infectious people
deriv(Ip[1:nAges, 1]) <- delta * Ep[i,1] - gamma_p * Ip[i,1] + gamma_vaccine_I[nVaccStates] * Ip[i,nVaccStates]
deriv(Ip[1:nAges, 2:nVaccStates]) <- delta * Ep[i,j] - (gamma_p + gamma_vaccine_I[j]) * Ip[i,j] + gamma_vaccine_I[(j-1)] * Ip[i, (j-1)]

# secondary infection
deriv(Ss[1:nAges, 1]) <- -incidence_rate_s[i,1] + (nu * R[i,1]) - (VaccState_rate[i,1] * Ss[i,1]) + (VaccState_rate[i, nVaccStates] * Ss[i, nVaccStates])
deriv(Ss[1:nAges, 2:nVaccStates]) <- -incidence_rate_s[i,j] + (nu * R[i,j]) + (VaccState_rate[i, (j-1)] * Ss[i,(j-1)]) - (VaccState_rate[i, j] * Ss[i,j])

deriv(Es[1:nAges, 1]) <- incidence_rate_s[i,1] - (delta + VaccState_rate[i,1]) * Es[i,1] + VaccState_rate[i,nVaccStates] * Es[i,nVaccStates]
deriv(Es[1:nAges, 2:nVaccStates]) <- incidence_rate_s[i,j] - (delta + VaccState_rate[i,j]) * Es[i,j] + (VaccState_rate[i,(j-1)] * Es[i, (j-1)])

deriv(Is[1:nAges, 1]) <- delta * Es[i,1] - gamma_s * Is[i,1] + gamma_vaccine_I[nVaccStates] * Is[i,nVaccStates]
deriv(Is[1:nAges, 2:nVaccStates]) <- delta * Es[i,j] - (gamma_s + gamma_vaccine_I[j]) * Is[i,j] + gamma_vaccine_I[(j-1)] * Is[i,(j-1)]

deriv(R[1:nAges, 1]) <- gamma_p * Ip[i,1] + gamma_s * Is[i,1] - (nu + VaccState_rate[i,1]) * R[i,1] + VaccState_rate[i, nVaccStates] * R[i, nVaccStates]
deriv(R[1:nAges, 2:nVaccStates]) <- gamma_p * Ip[i,j] + gamma_s * Is[i,j] - (nu + VaccState_rate[i,j]) * R[i,j] + (VaccState_rate[i,(j-1)] * R[i,(j-1)])

# used to tracks the incidence in a given age-group and vaccination status (the incidence of the age when infected)
deriv(Incidence[1:nAges, 1:nVaccStates]) <- incidence_rate_s[i,j] + incidence_rate_p[i,j]
output(DetIncidence[1:nAges, 1:nVaccStates]) <- prop_detected_vect[i] * Incidence[i,j]

output(Incidence_rate[1:nAges, 1:nVaccStates]) <- incidence_rate_s[i,j] + incidence_rate_p[i,j]
output(DetIncidence_rate[1:nAges, 1:nVaccStates]) <- prop_detected_vect[i] * Incidence_rate[i,j]

output(prev[1:nAges, 1:nVaccStates]) <- (Ip[i,j] + Is[i,j])/N[i,j]
output(prev_p[1:nAges, 1:nVaccStates]) <- Ip[i,j]/N[i,j]
output(prev_s[1:nAges, 1:nVaccStates]) <- Is[i,j]/N[i,j]

N_age[1:nAges] <- sum(N[i,]) # number of people in each age group

# used to calculate the number of doses given
deriv(doses[1:nAges, 1:nVaccStates]) <- VaccState_doses[i,j] * (Sp[i,j] + Ss[i,j] + Ep[i,j] + Es[i,j] + R[i,j])
dim(doses) <- c(nAges, nVaccStates)
# vac_cov[1:nAges] <- 1 - (Sp[i,1] + Ep[i,1] + Ip[i,1] + Ss[i,1] + Es[i,1] + Is[i,1] + R[i,1])/N_age[i]
# vac_cov_all <-  1 - (sum(Sp[,1]) + sum(Ep[,1]) + sum(Ip[,1]) + sum(Ss[,1]) + sum(Es[,1]) + sum(Is[,1]) + sum(R[,1]))/sum(N[,])

## Initial states:
initial(Sp[1:nAges, 1:nVaccStates]) <- Sp0[i,j]
initial(Ep[1:nAges, 1:nVaccStates]) <- Ep0[i,j]
initial(Ip[1:nAges, 1:nVaccStates]) <- Ip0[i,j]
initial(Ss[1:nAges, 1:nVaccStates]) <- Ss0[i,j]
initial(Es[1:nAges, 1:nVaccStates]) <- Es0[i,j]
initial(Is[1:nAges, 1:nVaccStates]) <- Is0[i,j]
initial(R[1:nAges, 1:nVaccStates]) <- R0[i,j]
initial(Incidence[1:nAges, 1:nVaccStates]) <- Incidence0[i,j]
initial(doses[1:nAges, 1:nVaccStates]) <- doses0[i,j]

##Initial vectors
Sp0 <- parameter()
Ep0 <- parameter()
Ip0 <- parameter()
Ss0 <- parameter()
Es0 <- parameter()
Is0 <- parameter()
R0 <- parameter()
Incidence0 <- parameter()
doses0 <- parameter()

##Dimensions of the different "vectors" used
# For the State Variables
dim(Sp) <- c(nAges, nVaccStates)
dim(Ep) <- c(nAges, nVaccStates)
dim(Ip) <- c(nAges, nVaccStates)
dim(R) <- c(nAges, nVaccStates)
dim(Ss) <- c(nAges, nVaccStates)
dim(Es) <- c(nAges, nVaccStates)
dim(Is) <- c(nAges, nVaccStates)
dim(Incidence) <- c(nAges, nVaccStates)

dim(N) <- c(nAges, nVaccStates)
dim(N_age) <- nAges
dim(DetIncidence) <- c(nAges, nVaccStates)
dim(Incidence_rate) <- c(nAges, nVaccStates)
dim(DetIncidence_rate) <- c(nAges, nVaccStates)
dim(prev) <- c(nAges, nVaccStates)
dim(prev_s) <- c(nAges, nVaccStates)
dim(prev_p) <- c(nAges, nVaccStates)
dim(incidence_rate_s) <- c(nAges, nVaccStates)
dim(incidence_rate_p) <- c(nAges, nVaccStates)
dim(Sp0) <- c(nAges, nVaccStates)
dim(Ep0) <- c(nAges, nVaccStates)
dim(Ip0) <- c(nAges, nVaccStates)
dim(Ss0) <- c(nAges, nVaccStates)
dim(Es0) <- c(nAges, nVaccStates)
dim(Is0) <- c(nAges, nVaccStates)
dim(R0) <- c(nAges, nVaccStates)
dim(Incidence0) <- c(nAges, nVaccStates)
dim(doses0) <- c(nAges, nVaccStates)
# dim(vac_cov) <- c(nAges)

## Calculate Vaccination Rate
##------------------------------------------------------------------------------

nVaccTimes <- parameter()

vaccine_times <- parameter()
dim(vaccine_times) <- c(nVaccTimes)

vaccine_cov <- parameter()
dim(vaccine_cov) <- c(nAges, nVaccTimes)

vaccine_period <- parameter()
dim(vaccine_period) <- nVaccTimes

vaccine_rate[1:nAges, 1:nVaccTimes] <- -log(1 - vaccine_cov[i,j]) / vaccine_period[j]
dim(vaccine_rate) <- c(nAges, nVaccTimes)

vr <- interpolate(vaccine_times, vaccine_rate, "constant")
dim(vr) <- nAges

output(VaccState_rate[1:nAges, 1:nVaccStates]) <- if(j == 1) vr[i] else gamma_vaccine[(j-1)]
dim(VaccState_rate) <- c(nAges, nVaccStates)

gamma_vaccine_I[1:nVaccStates] <- if(i==1) 0 else gamma_vaccine[(i-1)]
dim(gamma_vaccine_I) <- nVaccStates

VaccState_doses[1:nAges, 1:nVaccStates] <- if(j == 1) vr[i] else 0
dim(VaccState_doses) <- c(nAges, nVaccStates)
