# Time Steps
output(time) <- TRUE
nAges <- parameter()

## Model Parameters
alpha <- parameter() # reduced susceptibility for infection in age group i
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
mixing <- parameter()

dim(mixing) <- c(nAges, nAges)
dim(sigma_vect) <- nAges
dim(omega_vect) <- nAges
dim(prop_detected_vect) <- nAges

## Derivatives for Flows Between Compartments
##------------------------------------------------------------------------------

# prevalence of infection in age group j multiplied by reduced infectiousness in age group j (due to immunity early in life)
temp[] <- omega_vect[j] * (Is[j] + Ip[j]) / N[j]
# contacts multiplied by prevalence
s_ij[,] <- mixing[i,j] * temp[j]
lambda[] <- b0 * (1 + b1 * cos(2 * 3.14159265358979323846 * time / 12 + phi)) * sum(s_ij[i,])

infect_p[] <- lambda[i] * sigma_vect[i] * Sp[i]
infect_s[] <- lambda[i] * sigma_vect[i] * alpha[i] * Ss[i]

# primary infection
deriv(Sp[1:nAges]) <- -infect_p[i]
deriv(Ep[1:nAges]) <- infect_p[i] - delta * Ep[i]
deriv(Ip[1:nAges]) <- delta * Ep[i] - gamma_p * Ip[i]
# secondary infection
deriv(Ss[1:nAges]) <- -infect_s[i] + nu * R[i]
deriv(Es[1:nAges]) <- infect_s[i] - delta * Es[i]
deriv(Is[1:nAges]) <- delta * Es[i] - gamma_s * Is[i]

deriv(R[1:nAges]) <- gamma_p * Ip[i] + gamma_s * Is[i] - nu * R[i]

N[1:nAges] <- Ss[i] + Es[i] + Is[i] + Sp[i] + Ep[i] + Ip[i] + R[i]

deriv(Incidence[1:nAges]) <- infect_s[i] + infect_p[i]
deriv(DetIncidence[1:nAges]) <- prop_detected_vect[i] * infect_s[i] + prop_detected_vect[i] * infect_p[i]

## Initial states:
initial(Sp[1:nAges]) <- Sp0[i]
initial(Ep[1:nAges]) <- Ep0[i]
initial(Ip[1:nAges]) <- Ip0[i]
initial(Sp[1:nAges]) <- Ss0[i]
initial(Es[1:nAges]) <- Es0[i]
initial(Is[1:nAges]) <- Is0[i]
initial(R[1:nAges]) <- R0[i]
initial(Incidence[1:nAges]) <- Incidence0[i]
initial(DetIncidence[1:nAges]) <- DetIncidence0[i]

##Initial vectors
Sp0 <- parameter()
Ep0 <- parameter()
Ip0 <- parameter()
Ss0 <- parameter()
Es0 <- parameter()
Is0 <- parameter()
R0 <- parameter()
Incidence0 <- parameter()
DetIncidence0 <- parameter()

##Dimensions of the different "vectors" used
# For the State Variables
dim(Sp) <- nAges
dim(Ep) <- nAges
dim(Ip) <- nAges
dim(R) <- nAges
dim(Ss) <- nAges
dim(Es) <- nAges
dim(Is) <- nAges

dim(N) <- nAges
dim(Incidence) <- nAges
dim(DetIncidence) <- nAges
dim(Sp0) <- nAges
dim(Ep0) <- nAges
dim(Ip0) <- nAges
dim(Ss0) <- nAges
dim(Es0) <- nAges
dim(Is0) <- nAges
dim(R0) <- nAges

dim(Incidence0) <- nAges
dim(DetIncidence0) <- nAges
dim(lambda) <- nAges
dim(s_ij) <- c(nAges,nAges)
dim(temp) <- nAges
dim(infect) <- nAges
