#------ Model host system
#--- Host = Multi-state model with capture-mark-recapture analysis
        # Capture-mark-recapture data

#--- Host states

# 4 true states
  # Not entered
  # uninfected
  # infected
  # Dead

# 3 observed states
  # seen uninfected
  # seen infected
  # not seen

#--- Survey conditions
n.occasions <- 5   # Number of primary sampling occasions
n.surveys <- 1    # Secondary surveys

#---- Population parameters
# Super population size
N_U <- 100  # Number uninfected
N_I <- 100  # Number infected

N <- N_U + N_I  # Total population size

# Number of true states
n.states <- 4

# Number of observed states
n.obs <- 3

#--- Define parameter values
# Uninfected survival probability
phi_U <- 0.9

# Infected survival probability
phi_I <- 0.8

# Entry probability     
  # Must sum to one across all sampling ocassions
  # With n.occasions, n.occasions entry probability must be defined
gamma_U <- (1/(n.occasions))/4
gamma_I <- (1/(n.occasions))/4
 
# Transition probability 
beta_UI <- 0.5  # Going from Uninfected to Infected

beta_IU <- 0.2 # Going from Infected to Uninfected

# Detection probability
p_U <- 0.8  # Uninfected host
p_I <- 0.6  # Infected host

#---- Define function to simulate capture-recapture data 
simul.js <- function(
  
  n.occasions,  # Number of primary seasons
  n.surveys,    # Number of secondary surveys within primary seasons
  gamma_U,      # Uninfected entry probability 
  gamma_I,      # Infected entry probability 
  N,            # Total abundance
  N_U,          # Uninfected abundance
  N_I,          # Infected abundance
  p_U,          # Uninfected detection probability
  p_I           # Infected detection probability
  
){

########################################################   
#############   Table of Contents ######################       
########################################################  
  
  #---- This function has 5 parts
  # 1. Set up model parameters
  # 2. Assign true infection state
  # 3. Create observed capture history
  # 4.  Cleaning up the data

########################################################   
######   1.  Set up parameters for the model ###########    
########################################################  
#---- 1. Define matrix for state process
    # 4-dimensional matrix
        # 1. Leaving state
        # 2. Arriving state
        # 3. Host abundance
        # 4. Number of transitions
PSI.state <- array(NA, dim = c(n.states, n.states, N, n.occasions-1))

for(i in 1:N){
  for(j in 1:n.occasions-1){
    PSI.state[,,i,j] <- matrix(c(1 - gamma_U - gamma_I, gamma_U,  gamma_I,  0,
                                 0,     phi_U * (1- beta_UI), phi_U * beta_UI, 1-phi_U,
                                 0,     phi_I * beta_IU, phi_I * (1- beta_IU), 1-phi_I,
                                 0,     0,        0,       1), 
                               nrow = 4, ncol = 4, byrow = T)
  }
}
  
# Set the number of sampling occasions  
  n.occasions <- dim(PSI.state)[4] + 1
  
# Generate no. of entering hosts between occasions
  B_U <- rmultinom(1, N_U, rep(gamma_U, times = n.occasions))
  B_I <- rmultinom(1, N_I, rep(gamma_I, times = n.occasions))

# Calculate total number of entering hosts each season
  B <- B_U + B_I

# Create empty matricies to store true and survey data 
  CH.sur <- matrix(0, ncol = n.occasions, nrow = N)
  
# Define a vector with the occasion when individuals enter the population
  ent.occ_U <- ent.occ_I <- numeric()

# Populate vector with the occasion each individual enters the population  
for (t in 1:n.occasions){
    ent.occ_U <- c(ent.occ_U, rep(t, B_U[t]))
    ent.occ_I <- c(ent.occ_I, rep(t, B_I[t]))
}

# This is 1 vector for all individuals - regardless of infection status
ent.occ <- c(ent.occ_U, ent.occ_I)  

########################################################   
######   2.  Assign true infection state  ##############    
########################################################  

#-- Each uninfected hosts enters the population as uninfected (= 2)
for (i in 1:length(ent.occ_U)){     # For each individual in uninfected
  CH.sur[i, ent.occ_U[i]] <- 2   # 2 = uninfected
}

#-- Each infected hosts enters the population as infected (= 3)
for(i in (length(ent.occ_U)+1):length(ent.occ)){    
  CH.sur[i, ent.occ[i]] <- 3   # 3 = infected
}

#-- Now determine what happens to that individual after it enters the population
  # Does it survive?
  # Does it transistion disease states?
for (i in 1:N){    # For each individual in the super population       
  # If the entry occasion = the last occasion, go next
  if (ent.occ[i] == n.occasions) next   
  # For each occasion from when an individuals enters 2 the last sampling occasion    
    for(t in (ent.occ[i]+1):n.occasions){  
      # Determine individual state history
      sur <- which(rmultinom(1, 1, PSI.state[CH.sur[i, t-1], , i, t-1]) == 1)
      CH.sur[i, t] <- sur
    } #t
} #i

# Calculate the total number of individuals in the matricies created
Nt <- numeric(n.occasions)
for(i in 1:n.occasions){
  Nt[i] <- length(which(CH.sur[,i] > 0))
}
 
########################################################   
######   3.  Create observed capture history  ##########   
########################################################  

# Calculate the probability of misdiagnosing individuals
# This occurs when the swab or qPCR sample results in 0 zoospores
# You would use ySp and Spore to create the misclassification probability

#---- Define matrix for observation process
  # 4-dimensional matrix
    # 1. true state
    # 2. Observation state
    # 3. Host abundance
    # 4. Number of occasions
PSI.obs <- array(NA, dim = c(n.states, n.obs, N, n.occasions))

for(i in 1:N){
  for(j in 1:(n.occasions)){
    PSI.obs[,,i,j] <- matrix(c(0,       0,   1,
                              p_U,      0,   1-p_U,
                              0,      p_I,   1-p_I,
                              0,        0,   1), 
                             nrow = 4, ncol = 3, byrow = T)
    }
  }


#-- Create empty matrix to hold the data
# Dimensions of data
  # row = individual
  # column = primary period
  # sheets = secondary surveys
CH.p <- array(0, dim = c(N, n.occasions))


#-- 
for(i in 1:N){    # For each individual in the super population       
  # For each occasion from when an individuals enters 2 last occasion    
    for(j in ent.occ[i]:n.occasions){  
      # Determine if the individual is captured
        event <- which(rmultinom(1, 1, PSI.obs[CH.sur[i, j], , i, j]) == 1)
        CH.p[i, j] <- event
    } #js
} #is


########################################################   
######   4.  Cleaning up the data  ##################### 
########################################################  

# Dimensions of data
  # row = individual
  # column = primary period
  # sheets = secondary surveys

##########
# Remove individuals never captured
##########

never_cap <- numeric()

for(i in 1:dim(CH.p)[1]){

  # This only works if the entire capture history is == to 3 or ==0  
  n_seen <- grep("[03]", CH.p[i,])
  # Need to add capture history with 3s and 0s only
 # which(ch[i,,j] != 1)
  
  if(length(n_seen) == dim(CH.p)[2]){
    never_cap <- c(never_cap, i)
  }
  }

  
if(length(never_cap) > 0){
  CH.p <- CH.p[-never_cap,]
  CH.sur2 <- CH.sur[-never_cap,]
}
# Individuals never captured = 1
CH.sur[CH.sur == 0] <- 1
CH.sur2[CH.sur2 == 0] <- 1


return(list(CH.p = CH.p, 
            CH.sur = CH.sur, 
            CH.sur2 = CH.sur2,
            B = B,
            Nt = Nt
            ))
}

# Execute simulation function

sim <- simul.js(n.occasions = n.occasions, n.surveys = n.surveys,
                gamma_U = gamma_U, gamma_I = gamma_I, N = N,
                p_U = p_U, p_I = p_I,
                N_U = N_U, N_I = N_I 
)

#------------

dim(sim$CH.p)
dim(sim$CH.sur)
dim(sim$CH.sur2)

sim1 <- sim$CH.p
sim1[sim1 == 3] <- 0

write.csv(sim1, file = "/Users/Cici/EEMB595TE/6_Multi-state_Mark_Recapture_Models/Data/data_amphibian.csv")
