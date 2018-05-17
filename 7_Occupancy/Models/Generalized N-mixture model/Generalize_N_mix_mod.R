#Generalized N-mixture model assumptions (Dail & Madsen #2011):
#  (1) Sites and individual animals detections are #independent
#(2) The abundance at site i at time t only depends on the #abundance at site i at time t-1 (i.e., Markov property)
#(3) All animals present at site i at time t are assumed #to have the same detection probability


# Number of sites
M <- 500

# Number of repeated samples
T <- 5

# Average number of hosts per site
lambda1 <- 3

# Average recruitment rate
gamma1 <- 3

# Average survival rate
omega1 <- 0.9

# Average detection probability
p1 <- 0.4

# Days
days <- matrix(1, nrow = M, ncol = T)
  
# Data matrix + true number of hosts
y1 <- N1 <- matrix(NA, M, T)

# Survival + recruitment matrix
S1 <- G1 <- matrix(NA, M, T-1)

# First year= poisson process
N1[,1] <- rpois(M, lambda = lambda1)

# Second + years; survival + recruitment = abundance
for(t in 1:(T-1)) {
  S1[,t] <- rbinom(M, N1[,t], omega1^days[,t])

    # Discrete probability distribution describing the number of sucesses in a sequence of independent events, each with a probability p of sucess
  G1[,t] <- rpois(M, lambda = gamma1 * days[,t])

    # Discrete probability distribution describing the probability of a given number of events occurring in a fixed time interval if the events have an average rate + are idendentent of the time since last event
  N1[,t+1] <- S1[,t] + G1[,t]

    # Total abundance[next year] = survivors [last sampling] + gains[since last sampling]
}

y1[] <- rbinom(M*T, N1, p1)

sink("model.txt")
cat("
model{

#------------ Priors    

# Initial abundance
    lambda ~ dunif(-10, 10)

# Immigration rates
    gamma ~ dunif(-10, 10)

# Annual survivorship
    omega1 ~ dunif(0, 1)

# Detection probability
    p1 ~ dunif(0, 1)

#------------ Likelihood model   

for(j in 1:nSites){
    
# Year 1
    
    N[j, 1] ~ dpois(lambda1[j])
          log(lambda1[j]) <- lambda

    for(t in 2:nReps){
    
# Estimate Survivorship
    
    S[j, t-1] ~ dbin(omega1^days[j,t], N[j, t-1])

# Estimate recruitment
    G[j, t-1] ~ dpois(gamma1[j, t-1] * days[j,t])
          log(gamma1[j, t-1]) <- gamma

# Sum to get total N at location j in year t
    N[j, t] <- S[j, t-1] + G[j, t-1]  

    } #ts
    
 } # js

for(j in 1:nSites){

  for(t in 1:nReps){

    n[j, t] ~ dbin(p1, N[j, t])

  }
}
    
} 
    
    ", fill = TRUE)
sink()

#--------------------------- MCMC settings
ni <- 30000
nb <- 3000
nt <- 30
nc <- 3		


#--------------------------- Bundle the data

Ni <- array(NA, dim = c(M, T, 1))

y <- y1

win.data <- list(nSites = M, nReps = T, 
                 n = y, days = days)

#--------------------------- Initial values

Si <- Gi <- matrix(NA, nrow = M, ncol = T-1)

Ni <- y1+1

Si[] <- 1
Gi <- Ni[,-1] - Si

Ni[ ,-1] <- NA


inits <- function(){
  list(
    omega1 = 0.9, 
    p1 = 0.4, 
    gamma = 1,
    lambda = 1,
    N = Ni,
    S = Si,
    G = Gi
  )
}

#--------------------------- Parameters to monitor

params <- c("lambda", 
            "omega1",
            "p1", 
            "gamma")


#--------------------------- Run the model
library("jagsUI")

output2 <- jags(win.data, inits, params, "model.txt", 
                n.chains = nc, n.thin = nt, n.iter = ni, 
                n.burnin = nb, parallel = TRUE)

print(output2, dig = 3)

plot(output2)

#------------ Compare true to simulated values

sim <- c(log(lambda1), omega1, p1, log(gamma1))

names <- c("lambda", 
           "omega1",
           "p1",
           "gamma")

dd <- data.frame(names = names, TRUTH = sim)

MOD <- cbind(names, unlist(output2$mean[names]), 
             unlist(output2$q2.5[names]), 
             unlist(output2$q97.5[names]))

colnames(MOD)[1] <- "names" 

MER <- merge(dd, MOD, by = "names")

fores <- as.data.frame(MER)

for (i in 2:ncol(fores)){
  fores[,i] <- as.numeric(as.character(fores[,i]))
}

colnames(fores) <- c("x", "TRUTH", "y", "ylo", "yhi")

credplot.gg <- function(d){
  require(ggplot2)
  ggplot(d, aes(x=x, y=y, ymin=ylo, ymax=yhi))+geom_pointrange()+ xlab('Parameter Name') +   
    coord_flip() + ylab('Probability') +
   theme_bw() + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
}


credplot.gg(fores) + geom_point(aes(y = TRUTH, x = x), size = 2, col = "red")

