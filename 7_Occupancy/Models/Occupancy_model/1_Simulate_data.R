# Simulate data

# Define parameters
psi <- 0.6 # Occupancy probability
p <- 0.3   # Detection probability

# Define survey conditions
nSites <- 1000
nSurveys <- 3

# Simulate the truth
z <- rbinom(nSites, 1, psi)

# Simulate the observed data
y <- array(NA, dim = c(nSites, nSurveys))

for(i in 1:nSites){
  for(j in 1:nSurveys){
    y[i,j] <- rbinom(1, 1, z[i]*p)
  }
}

# y = observed data
# z = truth

# To save simulated data:
#write.csv(y, file = "/Users/Cici/EEMB595TE/7_Occupancy/ant_data.csv")
