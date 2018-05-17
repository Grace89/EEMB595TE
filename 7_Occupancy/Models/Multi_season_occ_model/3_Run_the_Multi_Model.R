# Bundle data
win.data <- list(y = data$y, nsite = dim(data$y)[1], nrep = dim(data$y)[2], nyear = dim(data$y)[3])

# Initial values
zst <- apply(win.data$y, c(1, 3), max)	# Observed occurrence as inits for z

inits <- function(){ list(z = zst)}

# Parameters monitored
params <- c("psi", "phi", "gamma", "p") 

# MCMC settings
ni <- 2500
nt <- 4
nb <- 500
nc <- 3

# Call JAGS from R (BRT 3 min)
out <- jags(win.data, inits, params, "Dynocc.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)

# Summarize posteriors
print(out, dig = 2)
