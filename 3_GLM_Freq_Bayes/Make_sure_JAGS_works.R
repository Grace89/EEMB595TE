
# Install packages
install.packages("rjags", dep =  TRUE)
install.packages("R2jags", dep =  TRUE)
install.packages("coda", dep =  TRUE)
install.packages("jagsUI", dep =  TRUE)

# Open the library
library("rjags")
library("R2jags")
library("coda")
library("jagsUI")

# Simulate the data
nSites <- 100

z <- rbinom(nSites, 1, p = 0.6)

# Make the model

sink("model.txt")
cat("
model{

#--- Priors

p ~ dunif(0, 1)

# Ecological model

for(i in 1:nSites){
  z[i] ~ dbern(p)
}

}
", fill = TRUE)
sink()

# Bundle the data

jags.data <- list(nSites = nSites,
                  z = z)

# Initial values
inits <- function(){list(p = runif(1, 0, 1))} 

# Parameters to monitor
params <- c("p")

# Set MCMC settings
ni <- 1000
nb <- 10
nt <- 1
nc <- 3

# Run the model
out <- jags(data = jags.data, inits = inits, parameters.to.save = params, model.file = "model.txt", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)

# Look at model output
out
