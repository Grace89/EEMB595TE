# Simulate data

# Define parameters
lambda <- 4 # Average number of animals per site
p <- 0.3    # Individual detection probability

# Define survey conditions
nSites <- 1000
nSurveys <- 3

# Simulate the truth
N <- rpois(nSites, lambda)

# Simulate the observed data
x <- array(NA, dim = c(nSites, nSurveys))

for(i in 1:nSites){
  for(j in 1:nSurveys){
    x[i,j] <- rbinom(1, N[i], p)
  }
}

# y = observed data
# z = truth

# To save simulated data:
write.csv(x, file = "/Users/Cici/EEMB595TE/7_Occupancy/ant_data_abund.csv")
