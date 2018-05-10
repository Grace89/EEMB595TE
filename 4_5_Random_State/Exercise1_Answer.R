# Scenario
# You are interested in the predictor variables of amphibian presense/absence.
# You do surveys for amphibians at 100 sites with varying forest cover.

# First, we simulate the data with the forest cover.

# Your job: Add the line(s) in the JAGS model that will analyze the data this way.


# Hint1: The way we simulate the data is almost identical to how we analyze the data.
# Hint2: Use the "logit" function in the JAGS model.
# Hint3: You will need to change the priors
# Hint4: We need to add something into the "jags.data" bundle.
# Hint5: We are no longer monitoring "p"- there are two other parameters we are interested in.


# Open the library
library("jagsUI")

# Simulate the data
nSites <- 100

# Parameters
alpha <- -3
beta <- 0.1

# Covariate- Percent forest cover
Forest_cover <- runif(nSites, min = 0, max = 100)

Forest_cover <- Forest_cover[order(Forest_cover)]

p <- plogis(alpha + beta * Forest_cover)

z <- rbinom(nSites, 1, p = p)

# Visulaize the data
plot(Forest_cover, jitter(z, 0.1), ylab = "Amphibian presence/absence", xlab = "Forest Cover percent", las= 1, pch = 21, col = "black", bg = "dodgerblue")

lines(Forest_cover, p, col = "red", lwd = 3)

# Make the model
sink("model.txt")
cat("
model{

#--- Priors

alpha ~ dnorm(0, 0.001)
beta ~ dnorm(0, 0.001)

# Ecological model

for(i in 1:nSites){
  z[i] ~ dbern(p[i])
    logit(p[i]) = alpha + beta * Forest_cover[i]
}

}
", fill = TRUE)
sink()

# Bundle the data

jags.data <- list(nSites = nSites,
                  z = z,
                  Forest_cover = Forest_cover)

# Initial values
inits <- function(){list(alpha = runif(1, -4, -2),
                         beta = runif(1, 0, 1))} 

# Parameters to monitor
params <- c("alpha", "beta")

# Set MCMC settings
ni <- 10000
nb <- 1000
nt <- 10
nc <- 3

# Run the model
out <- jags(data = jags.data, inits = inits, parameters.to.save = params, model.file = "model.txt", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)

# Look at model output
out

# Remeber
# alpha = -3
# beta= 0.1