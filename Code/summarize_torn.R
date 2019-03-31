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


##### negexp detrend linear climate M2-non-centered ####
load(file = "Results/JAGS/M2_nc.RData")
# Reconstruction
recon <- plot_recon(outM2_nc, obs = climate_df)
ggsave(filename = "Results/Figures/JAGS/torn_recon_post_negexp_linear.tiff", plot = recon, width = 8, height = 4, units = "in") # , dpi = 1000) # poster vs paper formatting - see past work and make package or github source
recon <- plot_recon(outM2_nc, obs = data.frame(year = years, value = x_full))
recon2 <- recon + theme_bw_poster() + ylim(0, 20) # + ggtitle("NegExp with Linear Climate")
ggsave(plot = recon2, filename = "Results/Figures/JAGS/negexp_norm_poster.pdf", dpi = 300)
recon2 <- recon + theme_bw_journal()
ggsave(plot = recon2, filename = "Results/Figures/JAGS/negexp_norm_paper.pdf", dpi = 1000)

recon_valid <- plot_recon(outM2_nc, valid_yrs = years[hold_out], obs = climate_df)
ggsave(filename = "Results/Figures/JAGS/torn_recon_negexp_linear_valid_back.tiff", plot = recon_valid, width = 8, height = 4, units = "in") 
recon_valid2 <- recon_valid + theme_bw_poster() + ylim(0, 20)
ggsave(plot = recon_valid2, filename = "Results/Figures/JAGS/negexp_norm_poster_valid.pdf", dpi = 300) # + ggtitle("NegExp with Linear Climate")

if(FALSE) {
# Validation plot
temp_valid <- temp_df_long %>%
  dplyr::filter(year %in% years[hold_out]) %>%
  dplyr::mutate(Value = "estimated")

climate_valid <- climate_df %>%
  dplyr::mutate(Value = "observed") %>%
  dplyr::rename(temp = x_full) %>%
  dplyr::filter(year %in% years[hold_out])

# temp_valid <- bind_rows(temp_valid, climate_df)

g <- ggplot(temp_valid, aes(x = year, y = temp))
g + geom_fan() + geom_line(data = climate_valid, aes(year, temp), colour = "red") + ylab("Temperature (C)") + xlab("Year") + scale_fill_distiller() + theme_bw()
ggsave(filename = "Results/Figures/torn_recon_negexp_linear_valid_back.tiff", width = 8, height = 4, units = "in") 
}

##### spline on age standardization #####
load(file = "Results/JAGS/detrend_spl.RData")

if(!dir.exists("Results/Figures/Detrend/Splines/")) dir.create("Results/Figures/Detrend/Splines/", recursive = TRUE)

if(FALSE) {
for(i in 1:M) {
  #tiff(paste0("Results/Figures/Detrend/Splines/detrend_tree_spl_", i, ".tif"), height = 12, width = 17, units = 'cm', compression = "lzw", res = 600)
  pdf(paste0("Results/Figures/Detrend/Splines/detrend_tree_spl_", i, ".pdf"), width = 8, height = 6, family = "Helvetica")
  # postscript(paste0("Results/Figures/Detrend/Splines/detrend_tree_spl_", i, ".eps"), width = 8, height = 6, horizontal = FALSE, onefile = FALSE, paper = "special", colormodel = "cmyk", family = "Courier")
  xidx = which(substr(varnames(out_detrend_spl),1,(8 + nchar(i)))==paste0("alpha0[", i, ","))
  postxm = colMeans(as.matrix(out_detrend_spl[,xidx])) 
  plot(year[f[i]]:year[l[i]], y_orig[i, f[i]:l[i]]/1000, xlab = "year", ylab = "Tree ring width (mm)", xlim = c(min(year, na.rm = TRUE), max(year, na.rm = TRUE)), ylim = c(0, 3))
  lines(year[f[i]]:year[l[i]], exp(postxm)/1000, col = "red", lwd = 2)
  dev.off()
}
dev.off()
}

if(FALSE) {
# reconstruction plot
temp_df <- as_tibble(t(as.matrix(out_detrend_spl[,x_id])))
temp_df <- temp_df %>% 
  mutate(year = year)

temp_df_long <- temp_df %>% 
  gather(key = sim, value = temp, -year) %>%
  dplyr::mutate(temp = temp*x_sd + x_mean)

# solid credible interval
x_post_quantiles <- data.frame(year = year, t(apply(x_post, 2, quantile, probs = c(0.025, 0.5, 0.975))))
names(x_post_quantiles) <- c("year", "lcri", "median", "ucri")
x_post_quantiles <- x_post_quantiles %>%
  dplyr::mutate(lcri = lcri*x_sd + x_mean,
                median = median*x_sd + x_mean,
                ucri = ucri*x_sd + x_mean) %>%
  as.data.frame()
ggplot(x_post_quantiles, aes(x = year, y = median)) + geom_ribbon(data = x_post_quantiles, aes(ymin = lcri, ymax = ucri), alpha=0.3) + geom_line(aes(year, median)) + theme_bw()

# scaled posterior interval
g <- ggplot(temp_df_long, aes(x = year, y = temp))
g + geom_fan() + geom_line(data = data.frame(year = year, median = x_post_median*x_sd + x_mean), aes(x = year, y = median), colour="black", size = 0.2) + ylab("Temperature (C)") + xlab("Year") + theme_bw() + scale_fill_distiller()
ggsave(filename = "Results/Figures/torn_recon_post_detrend_spl.tiff", width = 8, height = 4, units = "in") # , dpi = 1000) # poster vs paper formatting - see past work and make package or github source

# compare Schofield and my reconstructions
xidx_m2 = which(substr(varnames(outM2),1,2)=="x[")
postxm_m2 = colMeans(as.matrix(outM2[,xidx_m2])) 
xidx = which(substr(varnames(out_detrend_spl),1,2)=="x[")
postxm = colMeans(as.matrix(out_detrend_spl[,xidx])) 
xidx_con = which(substr(varnames(outM2_con),1,2)=="x[")
postxm_con = colMeans(as.matrix(outM2_con[,xidx_con])) 

par(mfrow = c(2, 2))
plot(year,postxm_m2*x_sd + x_mean,type="l", ylim = c(8, 17), main = "M2 non-centered")
plot(year,postxm_con*x_sd + x_mean,type="l", ylim = c(8, 17), main = "M2_CON")
plot(year,postxm*x_sd + x_mean,type="l", ylim = c(8, 17), main = "Detrend Spline")
plot(year, x_post_median*x_sd + x_mean, type="l", ylim = c(8, 17), main = "Spline detrend median")
par(mfrow = c(1, 1))

plot(year, post_climate_25*x_sd + x_mean,type="l")
plot(year, post_climate_50*x_sd + x_mean,type="l")
plot(year, post_climate_100*x_sd + x_mean,type="l")

plot(postxm_m2, postxm)
abline(0, 1, col = "red")

plot(postxm_con, postxm)
abline(0, 1, col = "red")

plot(x_full[hold_out], postxm[hold_out]*x_sd + x_mean)
abline(0, 1, col = "red")
}

# Reconstruction
recon <- plot_recon(out_detrend_spl, obs = climate_df)
ggsave(filename = "Results/Figures/JAGS/torn_recon_post_spline_linear.tiff", plot = recon, width = 8, height = 4, units = "in") # , dpi = 1000) # poster vs paper formatting - see past work and make package or github source
recon <- plot_recon(out_detrend_spl, obs = data.frame(year = years, value = x_full))
recon2 <- recon + theme_bw_poster() + ylim(0, 20) # # + ggtitle("Cubic Spline with Linear Climate")
ggsave(plot = recon2, filename = "Results/Figures/JAGS/spline_norm_poster.pdf", dpi = 300)
recon2 <- recon + theme_bw_journal()
ggsave(plot = recon2, filename = "Results/Figures/JAGS/spline_norm_paper.pdf", dpi = 1000)

recon_valid <- plot_recon(out_detrend_spl, valid_yrs = years[hold_out], obs = climate_df)
ggsave(filename = "Results/Figures/JAGS/torn_recon_spline_linear_valid_back.tiff", plot = recon_valid, width = 8, height = 4, units = "in") 
recon_valid2 <- recon_valid + theme_bw_poster() + ylim(0, 20)
ggsave(plot = recon_valid2, filename = "Results/Figures/JAGS/spline_norm_poster_valid.pdf", dpi = 300)

###### linear climate AR1 #####

load(file = "Results/JAGS/our_ar.RData")

# Reconstruction
recon <- plot_recon(our_ar, obs = climate_df)
ggsave(filename = "Results/Figures/JAGS/torn_recon_post_negexp_linear_ar.tiff", plot = recon, width = 8, height = 4, units = "in") # , dpi = 1000) # poster vs paper formatting - see past work and make package or github source
recon <- plot_recon(our_ar, obs = data.frame(year = years, value = x_full))
recon2 <- recon + theme_bw_poster() + ylim(0, 20)
ggsave(plot = recon2, filename = "Results/Figures/JAGS/negexp_linear_ar_poster.pdf", dpi = 300)
recon2 <- recon + theme_bw_journal()
ggsave(plot = recon2, filename = "Results/Figures/JAGS/negexp_linear_ar_paper.pdf", dpi = 1000)

recon_valid <- plot_recon(our_ar, valid_yrs = years[hold_out], obs = climate_df)
ggsave(filename = "Results/Figures/JAGS/torn_recon_spline_linear_valid_back.tiff", plot = recon_valid, width = 8, height = 4, units = "in") 
recon_valid2 <- recon_valid + theme_bw_poster() + ylim(0, 20)
ggsave(plot = recon_valid2, filename = "Results/Figures/JAGS/negexp_linear_ar_poster_valid.pdf", dpi = 300)





##### Climate Spline #####
library(splines)
load(file = "Results/JAGS/climate_spline_25.RData")

# Reconstruction
recon <- plot_recon(out_m_climate_spline_25, obs = climate_df)
ggsave(filename = "Results/Figures/JAGS/torn_recon_post_negexp_spl25.tiff", plot = recon, width = 8, height = 4, units = "in") # , dpi = 1000) # poster vs paper formatting - see past work and make package or github source
recon <- plot_recon(out_m_climate_spline_25, obs = data.frame(year = years, value = x_full))
recon2 <- recon + theme_bw_poster() + ylim(0, 20)
ggsave(plot = recon2, filename = "Results/Figures/JAGS/negexp_spl25_poster.pdf", dpi = 300)
recon2 <- recon + theme_bw_journal()
ggsave(plot = recon2, filename = "Results/Figures/JAGS/negexp_spl25_paper.pdf", dpi = 1000)

recon_valid <- plot_recon(out_m_climate_spline_25, valid_yrs = years[hold_out], obs = climate_df)
ggsave(filename = "Results/Figures/JAGS/torn_recon_negexp_spl25_valid_back.tiff", plot = recon_valid, width = 8, height = 4, units = "in") 
recon_valid2 <- recon_valid + theme_bw_poster() + ylim(0, 20)
ggsave(plot = recon_valid2, filename = "Results/Figures/JAGS/negexp_spl25_poster_valid.pdf", dpi = 300)

##### negexp detrend 1 changepoint climate - seems good ####
load(file = "Results/JAGS/negexp_1change.RData")

# Reconstruction
recon <- plot_recon(m_negexp_1change, obs = climate_df)
ggsave(filename = "Results/Figures/JAGS/torn_recon_post_negexp_1change.tiff", plot = recon, width = 8, height = 4, units = "in") # , dpi = 1000) # poster vs paper formatting - see past work and make package or github source
recon <- plot_recon(m_negexp_1change, obs = data.frame(year = years, value = x_full))
recon2 <- recon + theme_bw_poster() + ylim(0, 20)
ggsave(plot = recon2, filename = "Results/Figures/JAGS/negexp_1change_poster.pdf", dpi = 300)
recon2 <- recon + theme_bw_journal()
ggsave(plot = recon2, filename = "Results/Figures/JAGS/negexp_1change_paper.pdf", dpi = 1000)

recon_valid <- plot_recon(m_negexp_1change, valid_yrs = years[hold_out], obs = climate_df)
ggsave(filename = "Results/Figures/JAGS/torn_recon_negexp_1change_valid_back.tiff", plot = recon_valid, width = 8, height = 4, units = "in") 
recon_valid2 <- recon_valid + theme_bw_poster() + ylim(0, 20)
ggsave(plot = recon_valid2, filename = "Results/Figures/JAGS/negexp_1change_poster_valid.pdf", dpi = 300)

## any better than random points around the mean?



#### RCS (invariant) + climate spline 25 year knots #####
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

recon_valid <- plot_recon(m_rcs_spline_25, valid_yrs = years[hold_out], obs = climate_df)
ggsave(filename = "Results/Figures/JAGS/torn_recon_rcs_spline_25_valid_back.tiff", plot = recon_valid, width = 8, height = 4, units = "in") 
recon_valid2 <- recon_valid + theme_bw_poster() + ylim(0, 20) # + ggtitle("NegExp RCS with 25-yr Spline Climate")
ggsave(plot = recon_valid2, filename = "Results/Figures/JAGS/rcs_spline_25_poster_valid.pdf", dpi = 300)

rm(m_rcs_spline_25)
