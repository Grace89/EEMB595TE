
# Define parameters
psi <- 0.6 # Occupancy probability
p <- 0.3   # Detection probability

#--- Define parameter values
true <- c(psi,
          p
)

names <- c("Occupancy",
           "Detection"
)

names <- factor(names, levels = names)

mod.mean <- c(
  out$mean$psi,
  out$mean$p
)

mod.q2.5 <- c(
  out$q2.5$psi,
  out$q2.5$p
)

mod.q97.5 <- c(
  out$q97.5$psi,
  out$q97.5$p
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