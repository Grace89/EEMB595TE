# Analysis of the JS model as a multistate model

# Augment data approach allows us to estimate the total population size
# Idea = add a large number of 0 capture histories, and then the model will determine if those individuals were missed or do not exit using information from the detection probability and recruitment rates
nz <- 500

# Format the data with an added dummy column
# Dummy column = 1st column is all 1s
# Allows the model to calculate recrutiment into the population for the first season

CH.ms <- array(0, dim = c(dim(CH)[1]+nz, dim(CH)[2]+1))

for(i in 1:dim(CH)[1]){
  for(j in 1:dim(CH)[2]){
    CH.ms[i,j+1] <- CH[i, j]
  }
}

CH.du <- CH.ms[1:dim(CH)[1],]

# Recode CH matrix: a 0 is not allowed in JAGS!
CH.ms[CH.ms==0] <- 3    # Not seen = 3, seen = 1 or 2

# Bundle data
jags.data <- list(y = CH.ms, 
                  M = dim(CH.ms)[1],
                  K = dim(CH.ms)[2],
                  n.surv = 1,
                  state = 3
)

# Initial values
ch <- CH.du

js.multistate.init <- function(ch, nz){
  # 4 state system
  # 1 = Not entered
  # 2 = uninfected
  # 3 = infected
  # 4 = Dead 
  
  # 3 observation states
  # seen uninfected
  # seen infected
  # not seen
  
  # Put an NA when an individual was not seen ( = 0 or 3)
  ch[ch==0] <- NA
  ch[ch==3] <- NA
  
  state <- ch
  colnames(state) <- 1:dim(state)[2]
  
  # When the individual is known to be alive between first and last capture fill it in with 2s
  for(i in 1:dim(ch)[1]){    # For each individual
    
    n1 <- min(which(ch[i,] < 3))
    n2 <- max(which(ch[i,] < 3))
    
    fill <- which(is.na(state[i,n1:n2]) == TRUE)
    state[i, names(fill)] <- 2
  }
  state <- state + 1
  
  f <- array(NA, dim = c(dim(ch)[1]))
  
  for(i in 1:dim(ch)[1]){
    f[i] <- min(which(!is.na(ch[i, ])))
  }
  
  l <- array(NA, dim = c(dim(ch)[1]))
  
  for(i in 1:dim(ch)[1]){
    l[i] <- max(which(!is.na(ch[i, ])))
  }
  
  for (i in 1:dim(ch)[1]){
    # Before initial observation- not observed
    state[i,1:(f[i]-1)] <- 1
    
    # If the last time the animal was seen != the last survey date
    # Then add 1 to the survey date, and fill the rest to the last occasion with 4
    if(l[i]!= dim(ch)[2]){state[i, (l[i]+1):dim(ch)[2]] <- 4}
    
    if(ch[i, f[i]] == 1){state[i, f[i]] <- 2}
    if(ch[i, f[i]] == 2){state[i, f[i]] <- 3} 
  }
  
  state2 <- array(NA, dim = c(dim(ch)[1]+nz, dim(ch)[2]))
  state3 <- array(NA, dim = c(dim(ch)[1]+nz, dim(ch)[2]))
  
  # Copy over the matrix you generated to the one with the right dimensions    
  for(i in 1:dim(ch)[1]){
    for(j in 1:(dim(ch)[2])){
      state2[i,j] <- state[i, j]
    }
  }
  
  # For all the data augmented individuals- their state == 1
  for(i in (dim(state)[1]+1):dim(state2)[1]){ 
    for(j in 1:(dim(state)[2])){ # For all occasions
      state2[i,j] <- 1
    }
  }  
  
  return(state2)
}

n.occasions <- dim(CH.ms)[2]

zinit <- js.multistate.init(CH.du, nz)

zinit <- apply(zinit, c(1, 2), max)
zinit[,1] <- NA


#--- Bundle the inits 

inits <- function(){list(phi_U = phi_U, 
                         phi_I = phi_I, 
                         
                         p_U = p_U, 
                         p_I = p_I, 
                         
                         beta_UI = beta_UI,
                         beta_IU = beta_IU,
                         
                         gamma_U = rep(gamma_U, times = n.occasions-1), 
                         gamma_I = rep(gamma_I, times = n.occasions-1),
                         
                         z = zinit
)}    

#------- Parameters monitored

params <- c("phi_U", 
            "phi_I", 
            "p_U", 
            "p_I",
            "beta_UI",
            "beta_IU",
            "gamma_U",
            "gamma_I")

#------- MCMC settings

ni <- 10
nb <- 1
nt <- 1
nc <- 3

#------ call Library
library("jagsUI")

#------- Call JAGS from R

js.ms <- jags(data = jags.data, inits = inits, parameters.to.save = params, model.file = "model_JS.txt", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)

print(js.ms, dig = 3)