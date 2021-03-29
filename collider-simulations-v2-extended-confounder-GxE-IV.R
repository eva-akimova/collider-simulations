# "Heritable environments: bias due to conditioning on a collider in models with polygenic scores"
# Code for simulation analyses and figures
# 23-Sep-2020
# Additional analyses from 8-Mar-2021

# Packages needed
install.packages("dplyr")
install.packages("broom")
install.packages("purrr")
install.packages("mvtnorm")
install.packages("ggplot2")
install.packages("cowplot")
install.packages("tidyr")
install.packages("AER")
install.packages("forcats")

library(dplyr)
library(broom)
library(purrr)
library(mvtnorm)
library(ggplot2)
library(tidyr)
library(forcats)

######################### Figure 1B: bias in additive models ################################
# Scenario 1: Modest confounder, U
n <- 1000000
sim_data1 <- function(r, n = 1000000) {
  epsilon <- rmvnorm(n, c(0, 0), matrix(c(1, 0.12, 0.12, 1), nrow = 2))
  G <- rnorm(n, 0, 1) # Genes
  E <- r*G + epsilon[, 1] #E = exposure
  Y1 <- 0.6 * G + 0.6 * E + epsilon[, 2]
  tibble(PGScore = G,
         E = E,
         Outcome1 = Y1)
}

run_model1 <- function(df) {
  lm(Outcome1 ~ PGScore + E, df)
}

get_rsquared1 <- function(ml) {
  glance(ml)$r.squared[[1]]
}

get_pgscore_beta1 <- function(ml) {
  tidy(ml)[[2, 2]]
}

get_e_beta1 <- function(ml) {
  tidy(ml)[[3, 2]]
}

r_range <- seq(0, 0.5, 0.05)
dat <- map(r_range, sim_data1) %>%
  map(run_model1)

rsqs1 <- map(dat, get_rsquared1) %>% unlist()
pgs_betas1 <- map(dat, get_pgscore_beta1) %>% unlist()
e_betas1 <- map(dat, get_e_beta1) %>% unlist()

# Scenario 2: Moderate confounder, U
n <- 1000000
sim_data2 <- function(r, n = 1000000) {
  epsilon <- rmvnorm(n, c(0, 0), matrix(c(1, 0.25, 0.25, 1), nrow = 2))
  G <- rnorm(n, 0, 1) # Genes
  E <- r*G + epsilon[, 1] #E = exposure
  Y2 <- 0.6 * G + 0.6 * E + epsilon[, 2]
  tibble(PGScore = G,
         E = E,
         Outcome2 = Y2)
}

run_model2 <- function(df) {
  lm(Outcome2 ~ PGScore + E, df)
}

get_rsquared2 <- function(ml) {
  glance(ml)$r.squared[[1]]
}

get_pgscore_beta2 <- function(ml) {
  tidy(ml)[[2, 2]]
}

get_e_beta2 <- function(ml) {
  tidy(ml)[[3, 2]]
}

r_range <- seq(0, 0.5, 0.05)
dat <- map(r_range, sim_data2) %>%
  map(run_model2)

rsqs2 <- map(dat, get_rsquared2) %>% unlist()
pgs_betas2 <- map(dat, get_pgscore_beta2) %>% unlist()
e_betas2 <- map(dat, get_e_beta2) %>% unlist()

# Scenario 3: Strong confounder, U
n <- 1000000
sim_data3 <- function(r, n = 1000000) {
  epsilon <- rmvnorm(n, c(0, 0), matrix(c(1, 0.38, 0.38, 1), nrow = 2))
  G <- rnorm(n, 0, 1) # Genes
  E <- r*G + epsilon[, 1] #E = exposure
  Y3 <- 0.6 * G + 0.6 * E + epsilon[, 2]
  tibble(PGScore = G,
         E = E,
         Outcome3 = Y3)
}

run_model3 <- function(df) {
  lm(Outcome3 ~ PGScore + E, df)
}

get_rsquared3 <- function(ml) {
  glance(ml)$r.squared[[1]]
}

get_pgscore_beta3 <- function(ml) {
  tidy(ml)[[2, 2]]
}

get_e_beta3 <- function(ml) {
  tidy(ml)[[3, 2]]
}

r_range <- seq(0, 0.5, 0.05)
dat <- map(r_range, sim_data3) %>%
  map(run_model3)

rsqs3 <- map(dat, get_rsquared3) %>% unlist()
pgs_betas3 <- map(dat, get_pgscore_beta3) %>% unlist()
e_betas3 <- map(dat, get_e_beta3) %>% unlist()

## Betas inflation figure
## Environment
# Scenario 1
plot_ebeta_inflation1 <- 
  tibble(r = r_range,
         'E Beta' = e_betas1) %>%
  gather('key', 'value', 'E Beta')
plot_ebeta_inflation1
plot_ebeta_inflation1 <- transform( plot_ebeta_inflation1,
                                    base = 0.6 
)
plot_ebeta_inflation1 <- transform( plot_ebeta_inflation1,
                                    inflation = value*100/base-100
)
plot_ebeta_inflation1 <- transform( plot_ebeta_inflation1,
                                    case = '(1) Modest confounder, U'
)
plot_ebeta_inflation1 <- transform( plot_ebeta_inflation1,
                                    case2 = 'Coefficients inflation'
)
plot_ebeta_inflation1

# Scenario 2
plot_ebeta_inflation2 <- 
  tibble(r = r_range,
         'E Beta' = e_betas2) %>%
  gather('key', 'value', 'E Beta')
plot_ebeta_inflation2
plot_ebeta_inflation2 <- transform( plot_ebeta_inflation2,
                                    base = 0.6 
)
plot_ebeta_inflation2 <- transform( plot_ebeta_inflation2,
                                    inflation = value*100/base-100
)
plot_ebeta_inflation2 <- transform( plot_ebeta_inflation2,
                                    case = '(2) Moderate confounder, U'
)
plot_ebeta_inflation2 <- transform( plot_ebeta_inflation2,
                                    case2 = 'Coefficients inflation'
)
plot_ebeta_inflation2

# Scenario 3
plot_ebeta_inflation3 <- 
  tibble(r = r_range,
         'E Beta' = e_betas3) %>%
  gather('key', 'value', 'E Beta')
plot_ebeta_inflation3
plot_ebeta_inflation3 <- transform( plot_ebeta_inflation3,
                                    base = 0.6 
)
plot_ebeta_inflation3 <- transform( plot_ebeta_inflation3,
                                    inflation = value*100/base-100
)
plot_ebeta_inflation3 <- transform( plot_ebeta_inflation3,
                                    case = '(3) Strong confounder, U'
)
plot_ebeta_inflation3 <- transform( plot_ebeta_inflation3,
                                    case2 = 'Coefficients inflation'
)
plot_ebeta_inflation3

## Genes
# Scenario 1
plot_ebetg_inflation1 <- 
  tibble(r = r_range,
         'PGS Beta' = pgs_betas1) %>%
  gather('key', 'value', 'PGS Beta')
plot_ebetg_inflation1
plot_ebetg_inflation1 <- transform( plot_ebetg_inflation1,
                                    base = 0.6 
)
plot_ebetg_inflation1 <- transform( plot_ebetg_inflation1,
                                    inflation = value*100/base-100
)
plot_ebetg_inflation1 <- transform( plot_ebetg_inflation1,
                                    case = '(1) Modest confounder, U'
)
plot_ebetg_inflation1 <- transform( plot_ebetg_inflation1,
                                    case2 = 'Coefficients inflation'
)
plot_ebetg_inflation1

# Scenario 2
plot_ebetg_inflation2 <- 
  tibble(r = r_range,
         'PGS Beta' = pgs_betas2) %>%
  gather('key', 'value', 'PGS Beta')
plot_ebetg_inflation2
plot_ebetg_inflation2 <- transform( plot_ebetg_inflation2,
                                    base = 0.6 
)
plot_ebetg_inflation2 <- transform( plot_ebetg_inflation2,
                                    inflation = value*100/base-100
)
plot_ebetg_inflation2 <- transform( plot_ebetg_inflation2,
                                    case = '(2) Moderate confounder, U'
)
plot_ebetg_inflation2 <- transform( plot_ebetg_inflation2,
                                    case2 = 'Coefficients inflation'
)
plot_ebetg_inflation2

# Scenario 3
plot_ebetg_inflation3 <- 
  tibble(r = r_range,
         'PGS Beta' = pgs_betas3) %>%
  gather('key', 'value', 'PGS Beta')
plot_ebetg_inflation3
plot_ebetg_inflation3 <- transform( plot_ebetg_inflation3,
                                    base = 0.6 
)
plot_ebetg_inflation3 <- transform( plot_ebetg_inflation3,
                                    inflation = value*100/base-100
)
plot_ebetg_inflation3 <- transform( plot_ebetg_inflation3,
                                    case = '(3) Strong confounder, U'
)
plot_ebetg_inflation3 <- transform( plot_ebetg_inflation3,
                                    case2 = 'Coefficients inflation'
)
plot_ebetg_inflation3

#R-square
# Scenario 1: Modest confounder, U
plot_rsq_inflation1 <- 
  tibble(r = r_range,
         'R-Squared' = rsqs1) %>%
  gather('key', 'value', 'R-Squared')
plot_rsq_inflation1
# Following the discussion and formulas in Supplemental Data, the true 
# share of the variance in Y explained by G and E 
# in the context of this simulation is 0.461178,
plot_rsq_inflation1 <- transform( plot_rsq_inflation1,
                                  base = 0.461178 
)
plot_rsq_inflation1 <- transform( plot_rsq_inflation1,
                                  inflation = value*100/base-100
)
plot_rsq_inflation1 <- transform( plot_rsq_inflation1,
                                  case = '(1) Modest confounder, U'
)
plot_rsq_inflation1 <- transform( plot_rsq_inflation1,
                                  case2 = 'R-square change'
)
plot_rsq_inflation1

# Scenario 2: Moderate confounder, U
plot_rsq_inflation2 <- 
  tibble(r = r_range,
         'R-Squared' = rsqs2) %>%
  gather('key', 'value', 'R-Squared')
plot_rsq_inflation2
# Following the discussion and formulas in Supplemental Data, the true 
# share of the variance in Y explained by G and E 
# in the context of this simulation is 0.461178,
plot_rsq_inflation2 <- transform( plot_rsq_inflation2,
                                  base = 0.461178 
)
plot_rsq_inflation2 <- transform( plot_rsq_inflation2,
                                  inflation = value*100/base-100
)
plot_rsq_inflation2 <- transform( plot_rsq_inflation2,
                                  case = '(2) Moderate confounder, U'
)
plot_rsq_inflation2 <- transform( plot_rsq_inflation2,
                                  case2 = 'R-square change'
)
plot_rsq_inflation2

# Scenario 3: Moderate confounder, U
plot_rsq_inflation3 <- 
  tibble(r = r_range,
         'R-Squared' = rsqs3) %>%
  gather('key', 'value', 'R-Squared')
plot_rsq_inflation3
# Following the discussion and formulas in Supplemental Data, the true 
# share of the variance in Y explained by G and E 
# in the context of this simulation is 0.461178,
plot_rsq_inflation3 <- transform( plot_rsq_inflation3,
                                  base = 0.461178 
)
plot_rsq_inflation3 <- transform( plot_rsq_inflation3,
                                  inflation = value*100/base-100
)
plot_rsq_inflation3 <- transform( plot_rsq_inflation3,
                                  case = '(3) Strong confounder, U'
)
plot_rsq_inflation3 <- transform( plot_rsq_inflation3,
                                  case2 = 'R-square change'
)
plot_rsq_inflation3

combinebetas1 = rbind(plot_ebeta_inflation1, plot_ebetg_inflation1)
combinebetas2 = rbind(plot_ebeta_inflation2, plot_ebetg_inflation2)
combinebetas3 = rbind(plot_ebeta_inflation3, plot_ebetg_inflation3)

combinebetas = rbind(combinebetas1, combinebetas2, combinebetas3, 
                     plot_rsq_inflation1, plot_rsq_inflation2, 
                     plot_rsq_inflation3)
combinebetas <- mutate(combinebetas, key = as_factor(key))

plt <- ggplot(combinebetas, aes(x=r, y=inflation, group=key, color = key, shape = key))+ 
  geom_line(size = 1)+
  geom_point(size = 3.5) +
  scale_shape_manual(values=c(16,15)) +
  scale_y_continuous(limits=c(-40, 70)) +
  geom_hline(yintercept = 0, color = "gray54", linetype = "longdash", size = 0.9) +
  cowplot::theme_cowplot()
plt + theme(legend.position = "bottom") + 
  facet_grid(cols = vars(case), rows = vars(case2), scales = "free")  +
  theme(
    panel.background = element_rect(fill = "white",
                                    colour = "white",
                                    size = 0.45, linetype = "solid"),
    panel.grid.major = element_line(size = 0.45, linetype = 'solid',
                                    colour = "gray92"), 
    panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                    colour = "gray92"),
    panel.border=element_rect(colour="gray70",size=0.65)
  ) + labs( #title="Coefficients inflation depending on rGE values and unobserved confounder, U",
    x ="Gene-environment correlation (rGE)", y = "% inflation") +
  scale_colour_discrete(name  ="Coefficients",
                        breaks=c("E Beta", "PGS Beta", "R-square"),
                        labels=c("Environment", "PGScore", "R-square"), 
                        type = c("cornflowerblue", "springgreen3", "orchid4")) +
  scale_shape_discrete(name  ="Coefficients",
                       breaks=c("E Beta", "PGS Beta", "R-square"),
                       labels=c("Environment", "PGScore", "R-square"))

#####################################################################################################
######################### Figure 2B: bias in gene-environment models ################################
# Scenario 1: Modest confounder, U
n <- 1000000
sim_data1 <- function(r, n = 1000000) {
  epsilon <- rmvnorm(n, c(0, 0, 0), matrix(c(1, 0.12, 0.12, 0.12, 1, 0.12, 0.12, 0.12, 1), nrow = 3))
  G <- rnorm(n, 0, 1) # Genes
  E <- r*G + epsilon[, 1] #E = exposure
  GxE <- E*G # Interaction term
  Y1 <- 0.6 * G + 0.6 * E + 0.1 * GxE + epsilon[, 2]
  tibble(PGScore = G,
         E = E,
         GxE = GxE,
         Outcome1 = Y1)
}

run_model1 <- function(df) {
  lm(Outcome1 ~ PGScore + E + GxE, df)
}

get_rsquared1 <- function(ml) {
  glance(ml)$r.squared[[1]]
}

get_pgscore_beta1 <- function(ml) {
  tidy(ml)[[2, 2]]
}

get_e_beta1 <- function(ml) {
  tidy(ml)[[3, 2]]
}

get_gxe_beta1 <- function(ml) {
  tidy(ml)[[4, 2]]
}

r_range <- seq(0, 0.5, 0.05)
dat <- map(r_range, sim_data1) %>%
  map(run_model1)

rsqs1 <- map(dat, get_rsquared1) %>% unlist()
pgs_betas1 <- map(dat, get_pgscore_beta1) %>% unlist()
e_betas1 <- map(dat, get_e_beta1) %>% unlist()
gxe_betas1 <- map(dat, get_gxe_beta1) %>% unlist()

# Scenario 2: Moderate confounder, U
n <- 1000000
sim_data2 <- function(r, n = 1000000) {
  epsilon <- rmvnorm(n, c(0, 0, 0), matrix(c(1, 0.25, 0.25, 0.25, 1, 0.25, 0.25, 0.25, 1), nrow = 3))
  G <- rnorm(n, 0, 1) # Genes
  E <- r*G + epsilon[, 1] #E = exposure
  GxE <- E*G # Interaction term
  Y2 <- 0.6 * G + 0.6 * E + 0.1 * GxE + epsilon[, 2]
  tibble(PGScore = G,
         E = E,
         GxE = GxE,
         Outcome2 = Y2)
}

run_model2 <- function(df) {
  lm(Outcome2 ~ PGScore + E + GxE, df)
}

get_rsquared2 <- function(ml) {
  glance(ml)$r.squared[[1]]
}

get_pgscore_beta2 <- function(ml) {
  tidy(ml)[[2, 2]]
}

get_e_beta2 <- function(ml) {
  tidy(ml)[[3, 2]]
}

get_gxe_beta2 <- function(ml) {
  tidy(ml)[[4, 2]]
}

r_range <- seq(0, 0.5, 0.05)
dat <- map(r_range, sim_data2) %>%
  map(run_model2)

rsqs2 <- map(dat, get_rsquared2) %>% unlist()
pgs_betas2 <- map(dat, get_pgscore_beta2) %>% unlist()
e_betas2 <- map(dat, get_e_beta2) %>% unlist()
gxe_betas2 <- map(dat, get_gxe_beta2) %>% unlist()

# Scenario 3: Strong confounder, U
n <- 1000000
sim_data3 <- function(r, n = 1000000) {
  epsilon <- rmvnorm(n, c(0, 0, 0), matrix(c(1, 0.38, 0.38, 0.38, 1, 0.38, 0.38, 0.38, 1), nrow = 3))
  G <- rnorm(n, 0, 1) # Genes
  E <- r*G + epsilon[, 1] #E = exposure
  GxE <- E*G # Interaction term
  Y3 <- 0.6 * G + 0.6 * E + 0.1 * GxE + epsilon[, 2]
  tibble(PGScore = G,
         E = E,
         GxE = GxE,
         Outcome3 = Y3)
}

run_model3 <- function(df) {
  lm(Outcome3 ~ PGScore + E + GxE, df)
}

get_rsquared3 <- function(ml) {
  glance(ml)$r.squared[[1]]
}

get_pgscore_beta3 <- function(ml) {
  tidy(ml)[[2, 2]]
}

get_e_beta3 <- function(ml) {
  tidy(ml)[[3, 2]]
}

get_gxe_beta3 <- function(ml) {
  tidy(ml)[[4, 2]]
}

r_range <- seq(0, 0.5, 0.05)
dat <- map(r_range, sim_data3) %>%
  map(run_model3)

rsqs3 <- map(dat, get_rsquared3) %>% unlist()
pgs_betas3 <- map(dat, get_pgscore_beta3) %>% unlist()
e_betas3 <- map(dat, get_e_beta3) %>% unlist()
gxe_betas3 <- map(dat, get_gxe_beta3) %>% unlist()

## Betas inflation figure
## Environment
# Scenario 1
plot_ebeta_inflation1 <- 
  tibble(r = r_range,
         'E Beta' = e_betas1) %>%
  gather('key', 'value', 'E Beta')
plot_ebeta_inflation1
plot_ebeta_inflation1 <- transform( plot_ebeta_inflation1,
                                    base = 0.6 
)
plot_ebeta_inflation1 <- transform( plot_ebeta_inflation1,
                                    inflation = value*100/base-100
)
plot_ebeta_inflation1 <- transform( plot_ebeta_inflation1,
                                    case = '(1) Modest confounder, U'
)
plot_ebeta_inflation1 <- transform( plot_ebeta_inflation1,
                                    case2 = 'Coefficients inflation'
)
plot_ebeta_inflation1

# Scenario 2
plot_ebeta_inflation2 <- 
  tibble(r = r_range,
         'E Beta' = e_betas2) %>%
  gather('key', 'value', 'E Beta')
plot_ebeta_inflation2
plot_ebeta_inflation2 <- transform( plot_ebeta_inflation2,
                                    base = 0.6 
)
plot_ebeta_inflation2 <- transform( plot_ebeta_inflation2,
                                    inflation = value*100/base-100
)
plot_ebeta_inflation2 <- transform( plot_ebeta_inflation2,
                                    case = '(2) Moderate confounder, U'
)
plot_ebeta_inflation2 <- transform( plot_ebeta_inflation2,
                                    case2 = 'Coefficients inflation'
)
plot_ebeta_inflation2

# Scenario 3
plot_ebeta_inflation3 <- 
  tibble(r = r_range,
         'E Beta' = e_betas3) %>%
  gather('key', 'value', 'E Beta')
plot_ebeta_inflation3
plot_ebeta_inflation3 <- transform( plot_ebeta_inflation3,
                                    base = 0.6 
)
plot_ebeta_inflation3 <- transform( plot_ebeta_inflation3,
                                    inflation = value*100/base-100
)
plot_ebeta_inflation3 <- transform( plot_ebeta_inflation3,
                                    case = '(3) Strong confounder, U'
)
plot_ebeta_inflation3 <- transform( plot_ebeta_inflation3,
                                    case2 = 'Coefficients inflation'
)
plot_ebeta_inflation3

## Genes
# Scenario 1
plot_ebetg_inflation1 <- 
  tibble(r = r_range,
         'PGS Beta' = pgs_betas1) %>%
  gather('key', 'value', 'PGS Beta')
plot_ebetg_inflation1
plot_ebetg_inflation1 <- transform( plot_ebetg_inflation1,
                                    base = 0.6 
)
plot_ebetg_inflation1 <- transform( plot_ebetg_inflation1,
                                    inflation = value*100/base-100
)
plot_ebetg_inflation1 <- transform( plot_ebetg_inflation1,
                                    case = '(1) Modest confounder, U'
)
plot_ebetg_inflation1 <- transform( plot_ebetg_inflation1,
                                    case2 = 'Coefficients inflation'
)
plot_ebetg_inflation1

# Scenario 2
plot_ebetg_inflation2 <- 
  tibble(r = r_range,
         'PGS Beta' = pgs_betas2) %>%
  gather('key', 'value', 'PGS Beta')
plot_ebetg_inflation2
plot_ebetg_inflation2 <- transform( plot_ebetg_inflation2,
                                    base = 0.6 
)
plot_ebetg_inflation2 <- transform( plot_ebetg_inflation2,
                                    inflation = value*100/base-100
)
plot_ebetg_inflation2 <- transform( plot_ebetg_inflation2,
                                    case = '(2) Moderate confounder, U'
)
plot_ebetg_inflation2 <- transform( plot_ebetg_inflation2,
                                    case2 = 'Coefficients inflation'
)
plot_ebetg_inflation2

# Scenario 3
plot_ebetg_inflation3 <- 
  tibble(r = r_range,
         'PGS Beta' = pgs_betas3) %>%
  gather('key', 'value', 'PGS Beta')
plot_ebetg_inflation3
plot_ebetg_inflation3 <- transform( plot_ebetg_inflation3,
                                    base = 0.6 
)
plot_ebetg_inflation3 <- transform( plot_ebetg_inflation3,
                                    inflation = value*100/base-100
)
plot_ebetg_inflation3 <- transform( plot_ebetg_inflation3,
                                    case = '(3) Strong confounder, U'
)
plot_ebetg_inflation3 <- transform( plot_ebetg_inflation3,
                                    case2 = 'Coefficients inflation'
)
plot_ebetg_inflation3

## GxE
# Scenario 1
plot_gxe_inflation1 <- 
  tibble(r = r_range,
         'GxE Beta' = gxe_betas1) %>%
  gather('key', 'value', 'GxE Beta')
plot_gxe_inflation1
plot_gxe_inflation1 <- transform( plot_gxe_inflation1,
                                  base = 0.1 
)
plot_gxe_inflation1 <- transform( plot_gxe_inflation1,
                                  inflation = value*100/base-100
)
plot_gxe_inflation1 <- transform( plot_gxe_inflation1,
                                  case = '(1) Modest confounder, U'
)
plot_gxe_inflation1 <- transform( plot_gxe_inflation1,
                                  case2 = 'Coefficients inflation'
)
plot_gxe_inflation1

# Scenario 2
plot_gxe_inflation2 <- 
  tibble(r = r_range,
         'GxE Beta' = gxe_betas2) %>%
  gather('key', 'value', 'GxE Beta')
plot_gxe_inflation2
plot_gxe_inflation2 <- transform( plot_gxe_inflation2,
                                  base = 0.1 
)
plot_gxe_inflation2 <- transform( plot_gxe_inflation2,
                                  inflation = value*100/base-100
)
plot_gxe_inflation2 <- transform( plot_gxe_inflation2,
                                  case = '(2) Moderate confounder, U'
)
plot_gxe_inflation2 <- transform( plot_gxe_inflation2,
                                  case2 = 'Coefficients inflation'
)
plot_gxe_inflation2

# Scenario 3
plot_gxe_inflation3 <- 
  tibble(r = r_range,
         'GxE Beta' = gxe_betas3) %>%
  gather('key', 'value', 'GxE Beta')
plot_gxe_inflation3
plot_gxe_inflation3 <- transform( plot_gxe_inflation3,
                                  base = 0.1 
)
plot_gxe_inflation3 <- transform( plot_gxe_inflation3,
                                  inflation = value*100/base-100
)
plot_gxe_inflation3 <- transform( plot_gxe_inflation3,
                                  case = '(3) Strong confounder, U'
)
plot_gxe_inflation3 <- transform( plot_gxe_inflation3,
                                  case2 = 'Coefficients inflation'
)
plot_gxe_inflation3

#R-square
# Scenario 1: Modest confounder, U
plot_rsq_inflation1 <- 
  tibble(r = r_range,
         'R-Squared' = rsqs1) %>%
  gather('key', 'value', 'R-Squared')
plot_rsq_inflation1
# Following the discussion and formulas in Supplemental Data, the true share 
# of the variance in Y explained by G and E 
# in the context of this simulation is 0.4627,
plot_rsq_inflation1 <- transform( plot_rsq_inflation1,
                                  base = 0.4627 
)
plot_rsq_inflation1 <- transform( plot_rsq_inflation1,
                                  inflation = value*100/base-100
)
plot_rsq_inflation1 <- transform( plot_rsq_inflation1,
                                  case = '(1) Modest confounder, U'
)
plot_rsq_inflation1 <- transform( plot_rsq_inflation1,
                                  case2 = 'R-square change'
)
plot_rsq_inflation1

# Scenario 2: Moderate confounder, U
plot_rsq_inflation2 <- 
  tibble(r = r_range,
         'R-Squared' = rsqs2) %>%
  gather('key', 'value', 'R-Squared')
plot_rsq_inflation2
# Following the discussion and formulas in Supplemental Data, the true share 
# of the variance in Y explained by G and E 
# in the context of this simulation is 0.4627,
plot_rsq_inflation2 <- transform( plot_rsq_inflation2,
                                  base = 0.4627 
)
plot_rsq_inflation2 <- transform( plot_rsq_inflation2,
                                  inflation = value*100/base-100
)
plot_rsq_inflation2 <- transform( plot_rsq_inflation2,
                                  case = '(2) Moderate confounder, U'
)
plot_rsq_inflation2 <- transform( plot_rsq_inflation2,
                                  case2 = 'R-square change'
)
plot_rsq_inflation2

# Scenario 3: Moderate confounder, U
plot_rsq_inflation3 <- 
  tibble(r = r_range,
         'R-Squared' = rsqs3) %>%
  gather('key', 'value', 'R-Squared')
plot_rsq_inflation3
# Following the discussion and formulas in Supplemental Data, the true share 
# of the variance in Y explained by G and E 
# in the context of this simulation is 0.4627,
plot_rsq_inflation3 <- transform( plot_rsq_inflation3,
                                  base = 0.4627 
)
plot_rsq_inflation3 <- transform( plot_rsq_inflation3,
                                  inflation = value*100/base-100
)
plot_rsq_inflation3 <- transform( plot_rsq_inflation3,
                                  case = '(3) Strong confounder, U'
)
plot_rsq_inflation3 <- transform( plot_rsq_inflation3,
                                  case2 = 'R-square change'
)
plot_rsq_inflation3

combinebetas1 = rbind(plot_ebeta_inflation1, plot_ebetg_inflation1)
combinebetas2 = rbind(plot_ebeta_inflation2, plot_ebetg_inflation2)
combinebetas3 = rbind(plot_ebeta_inflation3, plot_ebetg_inflation3)

combinebetas = rbind(combinebetas1, combinebetas2, combinebetas3, 
                     plot_rsq_inflation1, plot_rsq_inflation2, plot_rsq_inflation3,
                     plot_gxe_inflation1, plot_gxe_inflation2, plot_gxe_inflation3)
combinebetas <- mutate(combinebetas, key = as_factor(key))

plt <- ggplot(combinebetas, aes(x=r, y=inflation, group=key))+ geom_line(aes(color=key), size = 1)+
  geom_point(aes(color=key, shape=key), size = 3.5) +
  scale_shape_manual(values=c(16,15)) +
  scale_y_continuous(limits=c(-40, 70)) +
  geom_hline(yintercept = 0, color = "gray54", linetype = "longdash", size = 0.9) +
  cowplot::theme_cowplot()
plt + theme(legend.position = "bottom") + 
  facet_grid(cols = vars(case), rows = vars(case2), scales = "free")  +
  theme(
    panel.background = element_rect(fill = "white",
                                    colour = "white",
                                    size = 0.45, linetype = "solid"),
    panel.grid.major = element_line(size = 0.45, linetype = 'solid',
                                    colour = "gray92"), 
    panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                    colour = "gray92"),
    panel.border=element_rect(colour="gray70",size=0.65),
    plot.title = element_text(size=14)
  ) + labs( #title="Coefficients inflation depending on rGE values and unobserved confounder, U",
    x ="Gene-environment correlation (rGE)", y = "% inflation") +
  scale_colour_discrete(name  ="Coefficients",
                        breaks=c("E Beta", "PGS Beta", "GxE Beta", "R-square"),
                        labels=c("Environment", "PGScore", "GxE interaction", "R-square"), 
                        type = c("cornflowerblue", "springgreen3", "orchid4", "hotpink1")) +
  scale_shape_discrete(name  ="Coefficients",
                       breaks=c("E Beta", "PGS Beta", "GxE Beta", "R-square"),
                       labels=c("Environment", "PGScore", "GxE interaction", "R-square"))

#####################################################################################################
######################### Figure S1: bias in gene-environment models ################################
# Scenario 1: Modest confounder, U
n <- 4000000
sim_data1 <- function(r, n = 4000000) {
  epsilon <- rmvnorm(n, c(0, 0, 0), matrix(c(1, 0.11, 0.11, 0.11, 1, 0.11, 0.11, 0.11, 1), nrow = 3))
  G <- rnorm(n, 0, 1) # Genes
  Z <- rnorm(n) #confounder that interacts with E
  E <- r*G + 0.1*Z + epsilon[, 1] #E = exposure
  GxE <- E*G # Interaction term
  Y1 <- 0.6 * G + 0.6 * E + 0.1 * GxE + 0.2*E*Z + 0.1*Z + epsilon[, 2]
  tibble(PGScore = G,
         E = E,
         GxE = GxE,
         Outcome1 = Y1)
}

run_model1 <- function(df) {
  lm(Outcome1 ~ PGScore + E + GxE, df)
}

get_rsquared1 <- function(ml) {
  glance(ml)$r.squared[[1]]
}

get_pgscore_beta1 <- function(ml) {
  tidy(ml)[[2, 2]]
}

get_e_beta1 <- function(ml) {
  tidy(ml)[[3, 2]]
}

get_gxe_beta1 <- function(ml) {
  tidy(ml)[[4, 2]]
}

r_range <- seq(0, 0.5, 0.05)
dat <- map(r_range, sim_data1) %>%
  map(run_model1)

rsqs1 <- map(dat, get_rsquared1) %>% unlist()
pgs_betas1 <- map(dat, get_pgscore_beta1) %>% unlist()
e_betas1 <- map(dat, get_e_beta1) %>% unlist()
gxe_betas1 <- map(dat, get_gxe_beta1) %>% unlist()

# Scenario 2: Moderate confounder, U
n <- 4000000
sim_data2 <- function(r, n = 4000000) {
  epsilon <- rmvnorm(n, c(0, 0, 0), matrix(c(1, 0.19, 0.19, 0.19, 1, 0.19, 0.19, 0.19, 1), nrow = 3))
  G <- rnorm(n, 0, 1) # Genes
  Z <- rnorm(n) #confounder that interacts with E
  E <- r*G + 0.3*Z + epsilon[, 1] #E = exposure
  GxE <- E*G # Interaction term
  Y2 <- 0.6 * G + 0.6 * E + 0.1 * GxE + 0.2*E*Z + 0.3*Z + epsilon[, 2]
  tibble(PGScore = G,
         E = E,
         GxE = GxE,
         Outcome2 = Y2)
}

run_model2 <- function(df) {
  lm(Outcome2 ~ PGScore + E + GxE, df)
}

get_rsquared2 <- function(ml) {
  glance(ml)$r.squared[[1]]
}

get_pgscore_beta2 <- function(ml) {
  tidy(ml)[[2, 2]]
}

get_e_beta2 <- function(ml) {
  tidy(ml)[[3, 2]]
}

get_gxe_beta2 <- function(ml) {
  tidy(ml)[[4, 2]]
}

r_range <- seq(0, 0.5, 0.05)
dat <- map(r_range, sim_data2) %>%
  map(run_model2)

rsqs2 <- map(dat, get_rsquared2) %>% unlist()
pgs_betas2 <- map(dat, get_pgscore_beta2) %>% unlist()
e_betas2 <- map(dat, get_e_beta2) %>% unlist()
gxe_betas2 <- map(dat, get_gxe_beta2) %>% unlist()

# Scenario 3: Strong confounder, U
n <- 4000000
sim_data3 <- function(r, n = 4000000) {
  epsilon <- rmvnorm(n, c(0, 0, 0), matrix(c(1, 0.25, 0.25, 0.25, 1, 0.25, 0.25, 0.25, 1), nrow = 3))
  G <- rnorm(n, 0, 1) # Genes
  Z <- rnorm(n) #confounder that interacts with E
  E <- r*G + 0.5*Z + epsilon[, 1] #E = exposure
  GxE <- E*G # Interaction term
  Y3 <- 0.6 * G + 0.6 * E + 0.1 * GxE + 0.2*E*Z + 0.5*Z + epsilon[, 2]
  tibble(PGScore = G,
         E = E,
         GxE = GxE,
         Outcome3 = Y3)
}

run_model3 <- function(df) {
  lm(Outcome3 ~ PGScore + E + GxE, df)
}

get_rsquared3 <- function(ml) {
  glance(ml)$r.squared[[1]]
}

get_pgscore_beta3 <- function(ml) {
  tidy(ml)[[2, 2]]
}

get_e_beta3 <- function(ml) {
  tidy(ml)[[3, 2]]
}

get_gxe_beta3 <- function(ml) {
  tidy(ml)[[4, 2]]
}

r_range <- seq(0, 0.5, 0.05)
dat <- map(r_range, sim_data3) %>%
  map(run_model3)

rsqs3 <- map(dat, get_rsquared3) %>% unlist()
pgs_betas3 <- map(dat, get_pgscore_beta3) %>% unlist()
e_betas3 <- map(dat, get_e_beta3) %>% unlist()
gxe_betas3 <- map(dat, get_gxe_beta3) %>% unlist()

## Betas inflation figure
## Environment
# Scenario 1
plot_ebeta_inflation1 <- 
  tibble(r = r_range,
         'E Beta' = e_betas1) %>%
  gather('key', 'value', 'E Beta')
plot_ebeta_inflation1
plot_ebeta_inflation1 <- transform( plot_ebeta_inflation1,
                                    base = 0.6 
)
plot_ebeta_inflation1 <- transform( plot_ebeta_inflation1,
                                    inflation = value*100/base-100
)
plot_ebeta_inflation1 <- transform( plot_ebeta_inflation1,
                                    case = '(1) Modest confounder, U'
)
plot_ebeta_inflation1 <- transform( plot_ebeta_inflation1,
                                    case2 = 'Coefficients inflation'
)
plot_ebeta_inflation1

# Scenario 2
plot_ebeta_inflation2 <- 
  tibble(r = r_range,
         'E Beta' = e_betas2) %>%
  gather('key', 'value', 'E Beta')
plot_ebeta_inflation2
plot_ebeta_inflation2 <- transform( plot_ebeta_inflation2,
                                    base = 0.6 
)
plot_ebeta_inflation2 <- transform( plot_ebeta_inflation2,
                                    inflation = value*100/base-100
)
plot_ebeta_inflation2 <- transform( plot_ebeta_inflation2,
                                    case = '(2) Moderate confounder, U'
)
plot_ebeta_inflation2 <- transform( plot_ebeta_inflation2,
                                    case2 = 'Coefficients inflation'
)
plot_ebeta_inflation2

# Scenario 3
plot_ebeta_inflation3 <- 
  tibble(r = r_range,
         'E Beta' = e_betas3) %>%
  gather('key', 'value', 'E Beta')
plot_ebeta_inflation3
plot_ebeta_inflation3 <- transform( plot_ebeta_inflation3,
                                    base = 0.6 
)
plot_ebeta_inflation3 <- transform( plot_ebeta_inflation3,
                                    inflation = value*100/base-100
)
plot_ebeta_inflation3 <- transform( plot_ebeta_inflation3,
                                    case = '(3) Strong confounder, U'
)
plot_ebeta_inflation3 <- transform( plot_ebeta_inflation3,
                                    case2 = 'Coefficients inflation'
)
plot_ebeta_inflation3

## Genes
# Scenario 1
plot_ebetg_inflation1 <- 
  tibble(r = r_range,
         'PGS Beta' = pgs_betas1) %>%
  gather('key', 'value', 'PGS Beta')
plot_ebetg_inflation1
plot_ebetg_inflation1 <- transform( plot_ebetg_inflation1,
                                    base = 0.6 
)
plot_ebetg_inflation1 <- transform( plot_ebetg_inflation1,
                                    inflation = value*100/base-100
)
plot_ebetg_inflation1 <- transform( plot_ebetg_inflation1,
                                    case = '(1) Modest confounder, U'
)
plot_ebetg_inflation1 <- transform( plot_ebetg_inflation1,
                                    case2 = 'Coefficients inflation'
)
plot_ebetg_inflation1

# Scenario 2
plot_ebetg_inflation2 <- 
  tibble(r = r_range,
         'PGS Beta' = pgs_betas2) %>%
  gather('key', 'value', 'PGS Beta')
plot_ebetg_inflation2
plot_ebetg_inflation2 <- transform( plot_ebetg_inflation2,
                                    base = 0.6 
)
plot_ebetg_inflation2 <- transform( plot_ebetg_inflation2,
                                    inflation = value*100/base-100
)
plot_ebetg_inflation2 <- transform( plot_ebetg_inflation2,
                                    case = '(2) Moderate confounder, U'
)
plot_ebetg_inflation2 <- transform( plot_ebetg_inflation2,
                                    case2 = 'Coefficients inflation'
)
plot_ebetg_inflation2

# Scenario 3
plot_ebetg_inflation3 <- 
  tibble(r = r_range,
         'PGS Beta' = pgs_betas3) %>%
  gather('key', 'value', 'PGS Beta')
plot_ebetg_inflation3
plot_ebetg_inflation3 <- transform( plot_ebetg_inflation3,
                                    base = 0.6 
)
plot_ebetg_inflation3 <- transform( plot_ebetg_inflation3,
                                    inflation = value*100/base-100
)
plot_ebetg_inflation3 <- transform( plot_ebetg_inflation3,
                                    case = '(3) Strong confounder, U'
)
plot_ebetg_inflation3 <- transform( plot_ebetg_inflation3,
                                    case2 = 'Coefficients inflation'
)
plot_ebetg_inflation3

## GxE
# Scenario 1
plot_gxe_inflation1 <- 
  tibble(r = r_range,
         'GxE Beta' = gxe_betas1) %>%
  gather('key', 'value', 'GxE Beta')
plot_gxe_inflation1
plot_gxe_inflation1 <- transform( plot_gxe_inflation1,
                                  base = 0.1 
)
plot_gxe_inflation1 <- transform( plot_gxe_inflation1,
                                  inflation = value*100/base-100
)
plot_gxe_inflation1 <- transform( plot_gxe_inflation1,
                                  case = '(1) Modest confounder, U'
)
plot_gxe_inflation1 <- transform( plot_gxe_inflation1,
                                  case2 = 'Coefficients inflation'
)
plot_gxe_inflation1

# Scenario 2
plot_gxe_inflation2 <- 
  tibble(r = r_range,
         'GxE Beta' = gxe_betas2) %>%
  gather('key', 'value', 'GxE Beta')
plot_gxe_inflation2
plot_gxe_inflation2 <- transform( plot_gxe_inflation2,
                                  base = 0.1 
)
plot_gxe_inflation2 <- transform( plot_gxe_inflation2,
                                  inflation = value*100/base-100
)
plot_gxe_inflation2 <- transform( plot_gxe_inflation2,
                                  case = '(2) Moderate confounder, U'
)
plot_gxe_inflation2 <- transform( plot_gxe_inflation2,
                                  case2 = 'Coefficients inflation'
)
plot_gxe_inflation2

# Scenario 3
plot_gxe_inflation3 <- 
  tibble(r = r_range,
         'GxE Beta' = gxe_betas3) %>%
  gather('key', 'value', 'GxE Beta')
plot_gxe_inflation3
plot_gxe_inflation3 <- transform( plot_gxe_inflation3,
                                  base = 0.1 
)
plot_gxe_inflation3 <- transform( plot_gxe_inflation3,
                                  inflation = value*100/base-100
)
plot_gxe_inflation3 <- transform( plot_gxe_inflation3,
                                  case = '(3) Strong confounder, U'
)
plot_gxe_inflation3 <- transform( plot_gxe_inflation3,
                                  case2 = 'Coefficients inflation'
)
plot_gxe_inflation3

#R-square
# Scenario 1: Modest confounder, U
plot_rsq_inflation1 <- 
  tibble(r = r_range,
         'R-Squared' = rsqs1) %>%
  gather('key', 'value', 'R-Squared')
plot_rsq_inflation1
# Following the discussion and formulas in Supplemental Data, the true share 
# of the variance in Y explained by G and E 
# in the context of this simulation is 0.4627,
plot_rsq_inflation1 <- transform( plot_rsq_inflation1,
                                  base = 0.4627 
)
plot_rsq_inflation1 <- transform( plot_rsq_inflation1,
                                  inflation = value*100/base-100
)
plot_rsq_inflation1 <- transform( plot_rsq_inflation1,
                                  case = '(1) Modest confounder, U'
)
plot_rsq_inflation1 <- transform( plot_rsq_inflation1,
                                  case2 = 'R-square change'
)
plot_rsq_inflation1

# Scenario 2: Moderate confounder, U
plot_rsq_inflation2 <- 
  tibble(r = r_range,
         'R-Squared' = rsqs2) %>%
  gather('key', 'value', 'R-Squared')
plot_rsq_inflation2
# Following the discussion and formulas in Supplemental Data, the true share 
# of the variance in Y explained by G and E 
# in the context of this simulation is 0.4627,
plot_rsq_inflation2 <- transform( plot_rsq_inflation2,
                                  base = 0.4627 
)
plot_rsq_inflation2 <- transform( plot_rsq_inflation2,
                                  inflation = value*100/base-100
)
plot_rsq_inflation2 <- transform( plot_rsq_inflation2,
                                  case = '(2) Moderate confounder, U'
)
plot_rsq_inflation2 <- transform( plot_rsq_inflation2,
                                  case2 = 'R-square change'
)
plot_rsq_inflation2

# Scenario 3: Moderate confounder, U
plot_rsq_inflation3 <- 
  tibble(r = r_range,
         'R-Squared' = rsqs3) %>%
  gather('key', 'value', 'R-Squared')
plot_rsq_inflation3
# Following the discussion and formulas in Supplemental Data, the true share 
# of the variance in Y explained by G and E 
# in the context of this simulation is 0.4627,
plot_rsq_inflation3 <- transform( plot_rsq_inflation3,
                                  base = 0.4627 
)
plot_rsq_inflation3 <- transform( plot_rsq_inflation3,
                                  inflation = value*100/base-100
)
plot_rsq_inflation3 <- transform( plot_rsq_inflation3,
                                  case = '(3) Strong confounder, U'
)
plot_rsq_inflation3 <- transform( plot_rsq_inflation3,
                                  case2 = 'R-square change'
)
plot_rsq_inflation3

combinebetas1 = rbind(plot_ebeta_inflation1, plot_ebetg_inflation1)
combinebetas2 = rbind(plot_ebeta_inflation2, plot_ebetg_inflation2)
combinebetas3 = rbind(plot_ebeta_inflation3, plot_ebetg_inflation3)

combinebetas = rbind(combinebetas1, combinebetas2, combinebetas3, 
                     plot_rsq_inflation1, plot_rsq_inflation2, plot_rsq_inflation3,
                     plot_gxe_inflation1, plot_gxe_inflation2, plot_gxe_inflation3)
combinebetas <- mutate(combinebetas, key = as_factor(key))

plt <- ggplot(combinebetas, aes(x=r, y=inflation, group=key))+ geom_line(aes(color=key), size = 1)+
  geom_point(aes(color=key, shape=key), size = 3.5) +
  scale_shape_manual(values=c(16,15)) +
  scale_y_continuous(limits=c(-40, 70)) +
  geom_hline(yintercept = 0, color = "gray54", linetype = "longdash", size = 0.9) +
  cowplot::theme_cowplot()
plt + theme(legend.position = "bottom") + 
  facet_grid(cols = vars(case), rows = vars(case2), scales = "free")  +
  theme(
    panel.background = element_rect(fill = "white",
                                    colour = "white",
                                    size = 0.45, linetype = "solid"),
    panel.grid.major = element_line(size = 0.45, linetype = 'solid',
                                    colour = "gray92"), 
    panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                    colour = "gray92"),
    panel.border=element_rect(colour="gray70",size=0.65),
    plot.title = element_text(size=14)
  ) + labs( #title="Coefficients inflation depending on rGE values and unobserved confounder, U",
    x ="Gene-environment correlation (rGE)", y = "% inflation") +
  scale_colour_discrete(name  ="Coefficients",
                        breaks=c("E Beta", "PGS Beta", "GxE Beta", "R-square"),
                        labels=c("Environment", "PGScore", "GxE interaction", "R-square"), 
                        type = c("cornflowerblue", "springgreen3", "orchid4", "hotpink1")) +
  scale_shape_discrete(name  ="Coefficients",
                       breaks=c("E Beta", "PGS Beta", "GxE Beta", "R-square"),
                       labels=c("Environment", "PGScore", "GxE interaction", "R-square"))

#####################################################################################################
# Instrumental variable solution example
library(AER)
n <- 1000000
epsilon <- rmvnorm(n, c(0, 0), matrix(c(1, 0.2, 0.2, 1), nrow = 2))
G <- rnorm(n, 0, 1) # Genes
# independent random error terms for D and Y
uD <- rnorm(n)
uY <- lm(rnorm(n) ~ uD)$residual
Z <- rnorm(n)  # instrumental variable predictor of E, environment
E <- 0.5*G + 0.5*Z + epsilon[, 1] + uD #E = exposure
Y <- 0.6 * G + 0.6 * E + epsilon[, 2] + uY

xHat <- lm(E ~ Z + G)$fitted.values
x <- lm(Y ~ xHat + G)
summary(x)





