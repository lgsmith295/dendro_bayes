#####################################################################
# Hierarchical reconstruction of temperature from Tornetrask Scots pine tree rings
#
# Daniel J. Hocking
# Based on Schofield et al. 2016, 2017 and Steinscheider et al. 2017
#####################################################################

###### Load libraries #####
library(rjags) # for running JAGS
library(MASS) # for boxcox
library(lattice) # for xyplot
library(dplyr)
# install.packages("devtools")
# devtools::install_github("mjskay/tidybayes")
library(tidybayes)
library(tidyverse)
library(ggplot2)
library(ggfan)
library(parallel)
require(gridExtra)

source("Code/functions.R")

load(file = "Results/JAGS/model_prep.RData") # climate_df, years, M, Tea, log_y, x_full, a_use, 







#### RCS + climate spline 25 year knots #####





###### 25-yr spline #####
load(file = "Results/JAGS/rcs_spline_25.RData")

# validation
x_id_25 = which(substr(varnames(m_rcs_spline_25),1,2)=="x[") # finds the indices of the x variables
post_climate_25 = colMeans(as.matrix(m_rcs_spline_25[,x_id_25])) # finds the posterior mean of the x variables
x_valid_post_m <- post_climate_25[hold_out]*x_sd + x_mean
x_valid <- x_full[hold_out]
plot(x_valid, x_valid_post_m, type = "p")
abline(0, 1, col = "red")
cor(x_valid, x_valid_post_m)

# R2

# RMSE
res2 <- (x_valid_post_m - x_valid)^2
sqrt(sum(res2) / length(x_valid_post_m))
sqrt(mean(res2))

# Reconstruction
recon <- plot_recon(m_rcs_spline_25, obs = data.frame(year = years, value = x_full))
recon2 <- recon + theme_bw_poster() + ylim(0, 20) # + ggtitle("NegExp RCS with 25-yr Spline Climate")
ggsave(plot = recon2, filename = "Results/Figures/JAGS/rcs_spline_25_poster.pdf", dpi = 300)
recon2 <- recon + theme_bw_journal()
ggsave(plot = recon2, filename = "Results/Figures/JAGS/rcs_spline_25_paper.pdf", dpi = 1000)

recon_valid <- plot_recon(m_rcs_spline_25, valid_yrs = year[hold_out], obs = climate_df)
ggsave(filename = "Results/Figures/JAGS/torn_recon_rcs_spline_25_valid_back.tiff", plot = recon_valid, width = 8, height = 4, units = "in") 
recon_valid2 <- recon_valid + theme_bw_poster() + ylim(0, 20) # + ggtitle("NegExp RCS with 25-yr Spline Climate")
ggsave(plot = recon_valid2, filename = "Results/Figures/JAGS/rcs_spline_25_poster_valid.pdf", dpi = 300)

rm(m_rcs_spline_25)
