sink("model.txt")
cat("
model 
{
# Priors

psi ~ dunif(0, 1)
p ~ dunif(0, 1)

# (1) Ecological model

for(i in 1:nSites){

  z[i] ~ dbern(psi)

  for(j in 1:nSurveys){

    y[i,j] ~ dbern(p.eff[i])

    p.eff[i] <-  z[i] * p

  } # j

} # i

}
",fill = TRUE)
sink()