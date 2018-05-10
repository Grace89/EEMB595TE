
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
beta_IU <- 0.2  # Going from Infected to Uninfected

# Detection probability
p_U <- 0.8  # Uninfected host
p_I <- 0.6  # Infected host

#--- Define parameter values
true <- c(phi_U, 
          phi_I, 
          p_U,
          p_I,
          beta_UI,
          beta_IU,
          
          rep(gamma_U, times = n.occasions-1), 
          rep(gamma_I, times = n.occasions-1)
)

names <- c("Survival Uninfected", 
           "Survival Infected", 
           "Detection Uninfected",
           "Detection Infected",
           "Transmission: U to I",
           "Recovery: I to U",
           paste("Recruit U T", 1:(n.occasions-1), sep = ""), 
           paste("Recruit I T", 1:(n.occasions-1), sep = "")
)

names <- factor(names, levels = names)

mod.mean <- c(
  js.ms$mean$phi_U,
  js.ms$mean$phi_I,
  js.ms$mean$p_U,
  js.ms$mean$p_I,
  js.ms$mean$beta_UI,
  js.ms$mean$beta_IU,
  js.ms$mean$gamma_U,
  js.ms$mean$gamma_I
)

mod.q2.5 <- c(
  js.ms$q2.5$phi_U,
  js.ms$q2.5$phi_I,
  js.ms$q2.5$p_U,
  js.ms$q2.5$p_I,
  js.ms$q2.5$beta_UI,
  js.ms$q2.5$beta_IU,
  js.ms$q2.5$gamma_U,
  js.ms$q2.5$gamma_I
)

mod.q97.5 <- c(
  js.ms$q97.5$phi_U,
  js.ms$q97.5$phi_I,
  js.ms$q97.5$p_U,
  js.ms$q97.5$p_I,
  js.ms$q97.5$beta_UI,
  js.ms$q97.5$beta_IU,
  js.ms$q97.5$gamma_U,
  js.ms$q97.5$gamma_I
)

dat <- data.frame(names = names, true = true, mod.mean = mod.mean, mod.q2.5 = mod.q2.5, mod.q97.5 = mod.q97.5)


library(ggplot2)

cols <- c("Truth" = "red", "Estimated" = "black")

ggplot(dat, aes(x= names, y=mod.mean, ymin=mod.q2.5, ymax=mod.q97.5))+ 
  geom_linerange(size = 1) +
  geom_point(size = 3, aes(x = names, y = mod.mean, col = "Estimated")) +
  geom_point(size = 3, aes(x = names, y = true, col = "Truth")) +
  scale_colour_manual("Values", values=cols)+
  geom_hline(yintercept = 0, lty=2) +
  coord_flip() + ylab('Parameter estimates') +
  xlab("Parameter names") +
  theme_bw()+ 
  theme(axis.text.x = element_text(size = 17, color = "black"), 
        axis.text.y = element_text(size = 17, color = "black"), 
        axis.title.y = element_text(size = 17, color = "black"), 
        axis.title.x =element_text(size = 17, color = "black"),
        legend.title =element_text(size = 17, color = "black"),
        legend.text =element_text(size = 17, color = "black"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) 