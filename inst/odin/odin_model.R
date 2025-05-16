# Time Steps
output(time) <- TRUE
nAges <- user()

## Model Parameters
b0 <- user()
b1 <- user()
phi <- user()
delta <- user()
gamma <- user()
nu <- user()
omega_vect[] <- user()
prop_detected_vect[] <- user()
sigma_vect[] <- user()
mixing[,] <- user()

dim(mixing) <- c(nAges, nAges)
dim(sigma_vect) <- nAges
dim(omega_vect) <- nAges
dim(prop_detected_vect) <- nAges

## Derivatives for Flows Between Compartments
##------------------------------------------------------------------------------
temp[] <- omega_vect[i] * I[i] / N[i]
s_ij[,] <- mixing[i,j] * temp[j]
lambda[] <- b0 * (1 + b1 * cos(2 * 3.14159265358979323846 * time / 12 + phi)) * sum(s_ij[i,])
infect[] <- lambda[i] * sigma_vect[i] * S[i]

deriv(S[1:nAges]) <- -infect[i] + nu * R[i]
deriv(E[1:nAges]) <- infect[i] - delta * E[i]
deriv(I[1:nAges]) <- delta*E[i] - gamma * I[i]
deriv(R[1:nAges]) <- gamma * I[i] - nu * R[i]
N[1:nAges] <- S[i] + E[i] + I[i] + R[i]
deriv(Incidence[1:nAges]) <- infect[i]
deriv(DetIncidence[1:nAges]) <- prop_detected_vect[i] * infect[i]

## Initial states:
initial(S[1:nAges]) <- S0[i]
initial(E[1:nAges]) <- E0[i]
initial(I[1:nAges]) <- I0[i]
initial(R[1:nAges]) <- R0[i]
initial(Incidence[1:nAges]) <- Incidence0[i]
initial(DetIncidence[1:nAges]) <- DetIncidence0[i]

##Initial vectors
S0[] <- user()
E0[] <- user()
I0[] <- user()
R0[] <- user()
Incidence0[] <- user()
DetIncidence0[] <- user()

##Dimensions of the different "vectors" used
# For the State Variables
dim(S) <- nAges
dim(E) <- nAges
dim(I) <- nAges
dim(R) <- nAges
dim(N) <- nAges
dim(Incidence) <- nAges
dim(DetIncidence) <- nAges
dim(S0) <- nAges
dim(E0) <- nAges
dim(I0) <- nAges
dim(R0) <- nAges
dim(Incidence0) <- nAges
dim(DetIncidence0) <- nAges
dim(lambda) <- nAges
dim(s_ij) <- c(nAges,nAges)
dim(temp) <- nAges
dim(infect) <- nAges
