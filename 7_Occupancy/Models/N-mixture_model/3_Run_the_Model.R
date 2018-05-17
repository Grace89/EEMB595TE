dat <- read.csv("~/EEMB595TE/7_Occupancy/ant_data_abund.csv")[,-1]

jags.data <- list( x = dat,
                   nSites = nrow(dat),
                   nSurveys = ncol(dat))

Ninit <- apply(dat, 1, max)

inits <- function(){list(lambda = runif(1, 0, 10),
                         p = runif(1, 0, 1),
                         N = Ninit
)}

#------- Parameters monitored

params <- c("lambda", 
            "p")

#------- MCMC settings
ni <- 1000
nb <- 100
nt <- 10
nc <- 3
na <- 100

#------ call Library
library("jagsUI")

#------- Call JAGS from R

out <- jags(data = jags.data, inits = inits, parameters.to.save = params, model.file = "model.txt", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, n.adapt = na, parallel = TRUE)

print(out, dig = 3)

plot(out)
