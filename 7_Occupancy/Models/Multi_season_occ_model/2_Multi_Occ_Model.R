sink("Dynocc.jags")
cat("
model {

# Priors
psi1 ~ dunif(0, 1)

for (k in 1:(nyear-1)){
  phi[k] ~ dunif(0, 1)
  gamma[k] ~ dunif(0, 1)
  p[k] ~ dunif(0, 1) 
}

p[nyear] ~ dunif(0, 1)

# Ecological submodel: Define state conditional on parameters
for (i in 1:nsite){
  z[i,1] ~ dbern(psi1)
    
  for (k in 2:nyear){
    
    muZ[i,k]<- z[i,k-1]*phi[k-1] + (1-z[i,k-1])*gamma[k-1]
    
    z[i,k] ~ dbern(muZ[i,k])
  } #k

} #i
    
# Observation model
for (i in 1:nsite){

  for (j in 1:nrep){

    for (k in 1:nyear){

      muy[i,j,k] <- z[i,k]*p[k]

      y[i,j,k] ~ dbern(muy[i,j,k])

    } #k

  } #j

} #i
    
}
",fill = TRUE)
sink()