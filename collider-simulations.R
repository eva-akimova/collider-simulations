# "Heritable environments: bias due to conditioning on a collider in models with polygenic scores"
# Code for simulation analyses and figures
# 23-Sep-2020

# Packages needed
install.packages("dplyr")
install.packages("broom")
install.packages("purrr")
install.packages("mvtnorm")
install.packages("ggplot2")
install.packages("cowplot")
install.packages("tidyr")

library(dplyr)
library(broom)
library(purrr)
library(mvtnorm)
library(ggplot2)
library(tidyr)

# Additive model
n <- 1000000
sim_data <- function(r, n = 1000000) {
  epsilon <- rmvnorm(n, c(0, 0), matrix(c(1, 0.2, 0.2, 1), nrow = 2))
  G <- rnorm(n, 0, 1) # Genes
  E <- r*G + epsilon[, 1] #E = exposure
  Y <- 0.6 * G + 0.6 * E + epsilon[, 2]
  tibble(PGScore = G,
         E = E,
         Outcome = Y)
}

run_model <- function(df) {
  lm(Outcome ~ PGScore + E, df)
}

get_rsquared <- function(ml) {
  glance(ml)$r.squared[[1]]
}

get_pgscore_beta <- function(ml) {
  tidy(ml)[[2, 2]]
}

get_e_beta <- function(ml) {
  tidy(ml)[[3, 2]]
}

r_range <- seq(0, 0.5, 0.05)
dat <- map(r_range, sim_data) %>%
  map(run_model)

rsqs <- map(dat, get_rsquared) %>% unlist()
pgs_betas <- map(dat, get_pgscore_beta) %>% unlist()
e_betas <- map(dat, get_e_beta) %>% unlist()

# R-square | Fifure 1C
plot_rsq_inflation <- 
  tibble(r = r_range,
         'R-Squared' = rsqs) %>%
  gather('key', 'value', 'R-Squared')
plot_rsq_inflation
# Following the discussion and furmulas in Supplemental Data, the true share of the variance in Y explained by G and E 
# in the context of this simulation is 0.461178,
plot_rsq_inflation <- transform( plot_rsq_inflation,
                                  base = 0.461178 
)
plot_rsq_inflation <- transform( plot_rsq_inflation,
                                  inflation = value*100/base-100
)

plot_rsq_inflation

plt <- ggplot(plot_rsq_inflation, aes(x=r, y=inflation, group=key))+ geom_line(size = 1, color="violetred4")+
  geom_point(size = 4, color="violetred4", shape=15) +
  geom_hline(yintercept = 0, linetype = "longdash", size = 1, color="violetred4") +
  scale_y_continuous(limits=c(0, 40)) +
  cowplot::theme_cowplot()
plt  + theme(legend.position = "bottom") +
  labs(title="R-square change depending on rGE values",
       x ="Gene-environment correlation (rGE)", y = "% change") +
  theme(
    panel.background = element_rect(fill = "white",
                                    colour = "white",
                                    size = 0.45, linetype = "solid"),
    panel.grid.major = element_line(size = 0.45, linetype = 'solid',
                                    colour = "gray92"), 
    panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                    colour = "gray92")
  ) 

# Coefficients | Figure 1B
plot_ebeta_inflation <- 
  tibble(r = r_range,
         'E Beta' = e_betas) %>%
  gather('key', 'value', 'E Beta')
plot_ebeta_inflation
plot_ebeta_inflation <- transform( plot_ebeta_inflation,
                                    base = 0.6 
)
plot_ebeta_inflation <- transform( plot_ebeta_inflation,
                                    inflation = value*100/base-100
)

plot_ebeta_inflation

plot_ebetg_inflation <- 
  tibble(r = r_range,
         'PGS Beta' = pgs_betas) %>%
  gather('key', 'value', 'PGS Beta')
plot_ebetg_inflation

plot_ebetg_inflation <- transform( plot_ebetg_inflation,
                                    base = 0.6 
)
plot_ebetg_inflation <- transform( plot_ebetg_inflation,
                                    inflation = value*100/base-100
)

plot_ebetg_inflation

combinebetas = rbind(plot_ebeta_inflation, plot_ebetg_inflation)

plt <- ggplot(combinebetas, aes(x=r, y=inflation, group=key))+ geom_line(size = 1, aes(color=key))+
  geom_point(size = 4, aes(color=key, shape=key)) +
  geom_hline(yintercept = 0, linetype = "longdash", size = 1, color="gray54") +
  scale_y_continuous(limits=c(-20, 40)) + 
  cowplot::theme_cowplot()  + 
  theme(legend.position = "bottom") +
  labs(title="Coefficients inflation depending on rGE values",
       x ="Gene-environment correlation (rGE)", y = "% inflation") +
  scale_colour_discrete(name  ="",
                        breaks=c("E Beta", "PGS Beta"),
                        labels=c("Environment", "PGScore")) +
  scale_shape_discrete(name  ="",
                       breaks=c("E Beta", "PGS Beta"),
                       labels=c("Environment", "PGScore")) +
  theme(
    panel.background = element_rect(fill = "white",
                                    colour = "white",
                                    size = 0.45, linetype = "solid"),
    panel.grid.major = element_line(size = 0.45, linetype = 'solid',
                                    colour = "gray92"), 
    panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                    colour = "gray92")
  ) 
plt
##############################################################################################
# Gene-environment model
n <- 1000000
sim_data <- function(r, n = 1000000) {
  epsilon <- rmvnorm(n, c(0, 0, 0), matrix(c(1, 0.2, 0.2, 0.2, 1, 0.2, 0.2, 0.2, 1), nrow = 3))
  G <- rnorm(n, 0, 1) # Genes
  E <- r*G + epsilon[, 1] #E = exposure
  GxE <- E*G # Interaction term
  Y <- 0.6 * G + 0.6 * E + 0.1 * GxE + epsilon[, 2]
  tibble(PGScore = G,
         E = E,
         GxE = GxE,
         Outcome = Y)
}

run_model <- function(df) {
  lm(Outcome ~ PGScore + E + GxE, df)
}

get_rsquared <- function(ml) {
  glance(ml)$r.squared[[1]]
}

get_pgscore_beta <- function(ml) {
  tidy(ml)[[2, 2]]
}

get_e_beta <- function(ml) {
  tidy(ml)[[3, 2]]
}

get_gxe_beta <- function(ml) {
  tidy(ml)[[4, 2]]
}

r_range <- seq(0, 0.5, 0.05)
dat <- map(r_range, sim_data) %>%
  map(run_model)

rsqs <- map(dat, get_rsquared) %>% unlist()
pgs_betas <- map(dat, get_pgscore_beta) %>% unlist()
e_betas <- map(dat, get_e_beta) %>% unlist()
gxe_betas <- map(dat, get_gxe_beta) %>% unlist()


# R-square | Fifure 2C
plot_rsq_inflation <- 
  tibble(r = r_range,
         'R-Squared' = rsqs) %>%
  gather('key', 'value', 'R-Squared')
plot_rsq_inflation
# Following the discussion and furmulas in Supplimental Data, the true share of the variance in Y explained by G and E 
# in the context of this simulation is 0.4627.
plot_rsq_inflation <- transform( plot_rsq_inflation,
                                 base = 0.4627 
)
plot_rsq_inflation <- transform( plot_rsq_inflation,
                                 inflation = value*100/base-100
)

plot_rsq_inflation

plt <- ggplot(plot_rsq_inflation, aes(x=r, y=inflation, group=key))+ geom_line(size = 1, color="violetred4")+
  geom_point(size = 6, color="violetred4", shape=18) +
  geom_hline(yintercept = 0, linetype = "longdash", size = 1, color="violetred4") +
  scale_y_continuous(limits=c(0, 35)) +
  cowplot::theme_cowplot()
plt  + theme(legend.position = "bottom") +
  labs(title="R-square inflation depending on rGE values",
       x ="Gene-environment correlation (rGE)", y = "% inflation") +
  theme(
    panel.background = element_rect(fill = "white",
                                    colour = "white",
                                    size = 0.45, linetype = "solid"),
    panel.grid.major = element_line(size = 0.45, linetype = 'solid',
                                    colour = "gray92"), 
    panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                    colour = "gray92")
  ) 

# Coefficients | Figure 2B
plot_ebeta_inflation <- 
  tibble(r = r_range,
         'E Beta' = e_betas) %>%
  gather('key', 'value', 'E Beta')
plot_ebeta_inflation
plot_ebeta_inflation <- transform( plot_ebeta_inflation,
                                   base = 0.6
)
plot_ebeta_inflation <- transform( plot_ebeta_inflation,
                                   inflation = value*100/base-100
)
plot_ebeta_inflation <- transform( plot_ebeta_inflation,
                                   case = 'beta E = beta G'
)
plot_ebeta_inflation

plot_ebetg_inflation <- 
  tibble(r = r_range,
         'PGS Beta' = pgs_betas) %>%
  gather('key', 'value', 'PGS Beta')
plot_ebetg_inflation
plot_ebetg_inflation <- transform( plot_ebetg_inflation,
                                   base = 0.6
)
plot_ebetg_inflation <- transform( plot_ebetg_inflation,
                                   inflation = value*100/base-100
)
plot_ebetg_inflation <- transform( plot_ebetg_inflation,
                                   case = 'beta E = beta G'
)
plot_ebetg_inflation

plot_gxe_inflation <- 
  tibble(r = r_range,
         'GxE Beta' = gxe_betas) %>%
  gather('key', 'value', 'GxE Beta')
plot_gxe_inflation
plot_gxe_inflation <- transform( plot_gxe_inflation,
                                 base = 0.1 
)
plot_gxe_inflation <- transform( plot_gxe_inflation,
                                 inflation = value*100/base-100
)
plot_gxe_inflation <- transform( plot_gxe_inflation,
                                 case = 'beta E = beta G'
)
plot_gxe_inflation

combinebetas = rbind(plot_ebeta_inflation, plot_ebetg_inflation, plot_gxe_inflation)

plt <- ggplot(combinebetas, aes(x=r, y=inflation, group=key))+ geom_line(size = 1, aes(color=key))+
  geom_point(size = 4, aes(color=key, shape=key)) +
  geom_hline(yintercept = 0, linetype = "longdash", size = 1, color="gray54") +
  scale_y_continuous(limits=c(-20, 40)) + 
  cowplot::theme_cowplot()  + 
  theme(legend.position = "bottom") +
  labs(title="Coefficients inflation depending on rGE values",
       x ="Gene-environment correlation (rGE)", y = "% inflation") +
  scale_colour_discrete(name  ="",
                        breaks=c("E Beta", "PGS Beta", "GxE Beta"),
                        labels=c("Environment", "PGScore", "GxE interaction")) +
  scale_shape_discrete(name  ="",
                       breaks=c("E Beta", "PGS Beta", "GxE Beta"),
                       labels=c("Environment", "PGScore", "GxE interaction")) +
  theme(
    panel.background = element_rect(fill = "white",
                                    colour = "white",
                                    size = 0.45, linetype = "solid"),
    panel.grid.major = element_line(size = 0.45, linetype = 'solid',
                                    colour = "gray92"), 
    panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                    colour = "gray92")
  ) 
plt
