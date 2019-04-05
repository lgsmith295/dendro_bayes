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
library(splines)

source("Code/functions.R")

load(file = "Results/NM/JAGS/model_prep_nm.RData")

##### negexp detrend linear climate M2-non-centered ####
load(file = "Results/NM/JAGS/negexp_norm.RData")

x_id_50 = which(substr(varnames(m_negexp_norm),1,2)=="x[")
post_climate_50 = colMeans(as.matrix(m_negexp_norm[,x_id_50]))
x_valid_post_m <- post_climate_50[hold_out]*x_sd + x_mean
x_valid <- x_full[hold_out]

performance_df <- perf_stats(est_valid = x_valid_post_m, observed = x_full, valid_id = hold_out, cal_id = cal_ids, mod_id = "negexp_linear")

# reconstruction plot
recon <- plot_recon(m_negexp_norm, obs = data.frame(year = years, value = x_full) , y_lab = "Mean Jan-Jul Precip. (cm)")
recon2 <- recon + theme_bw_poster()
ggsave(plot = recon2, filename = "Results/Figures/NM/negexp_norm_poster.pdf", dpi = 300)
recon2 <- recon + theme_bw_journal()
ggsave(plot = recon2, filename = "Results/Figures/NM/negexp_norm_paper.pdf", dpi = 1000)
ggsave(plot = recon2, filename = "Results/Figures/NM/negexp_norm_paper.png")

rm(m_negexp_norm)


##### negexp detrend 1 changepoint climate ####
x_min <- min(x_use, na.rm = TRUE)
x_max <- max(x_use, na.rm = TRUE)
load(file = "Results/NM/JAGS/negexp_1change.RData")

x_id_50 = which(substr(varnames(m_negexp_1change),1,2)=="x[")
post_climate_50 = colMeans(as.matrix(m_negexp_1change[,x_id_50]))
x_valid_post_m <- post_climate_50[hold_out]*x_sd + x_mean
x_valid <- x_full[hold_out]

performance_df[2, ] <- perf_stats(est_valid = x_valid_post_m, observed = x_full, valid_id = hold_out, cal_id = cal_ids, mod_id = "negexp_1change")

recon <- plot_recon(m_negexp_1change, obs = data.frame(year = years, value = x_full) , y_lab = "Mean Jan-Jul Precip. (cm)")

# consider doing observed values as points to better see where they fall within the credible interval
recon_valid <- plot_recon(m_negexp_1change, obs = data.frame(year = years, value = x_full) , y_lab = "Mean Jan-Jul Precip. (cm)", valid_yrs = years[hold_out])

# reconstruction plot
recon <- plot_recon(m_negexp_1change, obs = data.frame(year = years, value = x_full) , y_lab = "Mean Jan-Jul Precip. (cm)")
recon2 <- recon + theme_bw_poster()
ggsave(plot = recon2, filename = "Results/Figures/NM/negexp_1change_poster.pdf", dpi = 300)
recon2 <- recon + theme_bw_journal()
ggsave(plot = recon2, filename = "Results/Figures/NM/negexp_1change_paper.pdf", dpi = 1000)
ggsave(plot = recon2, filename = "Results/Figures/NM/negexp_1change_paper.png")

rm(m_negexp_1change)




