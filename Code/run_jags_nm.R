#####################################################################
# Hierarchical reconstruction of precip from southest mountains NM multispecies tree rings
# NOAA Climate Division CONUS 2904
# Annual precip to start
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
library(splines)

source("Code/functions.R")

##### Set conditions #####
testing <- FALSE


if(testing) {
  nb = 1000
  ni = 1000
  nt = 1
  nc = 3
} else {
  nb = 10000
  ni = 10000
  nt = 10
  nc = 4
}

if(!dir.exists("Results/NM/JAGS/")) dir.create("Results/NM/JAGS/", recursive = TRUE)
if(!dir.exists("Results/NM/Figures/")) dir.create("Results/NM/Figures/", recursive = TRUE)

#####

#### Load and Prep Data #####
# load(file = "Data/itrdb_pilo_mount_washington/pilo_rwl_climate.RData")
raw <- readRDS(file = "Data/az_nm/raw_correlated_2904.RData")

# Based on sample depth only do reconstruction of last 1000 years (and before 1500 is based completely on 1 species: PSME)
if(TRUE) {
cores <- raw %>%
  group_by(core) %>%
  dplyr::filter(!is.na(rwl)) %>%
  dplyr::select(core, year, rwl) %>%
  # dplyr::summarise(n()) %>%
  # dplyr::filter()
  summarise(f_yr = max(year),
            rwl = min(rwl)) %>%
  dplyr::filter(f_yr >= 1000) # 1510

raw <- raw %>%
  ungroup() %>%
  dplyr::filter(core %in% cores$core) # ,
                # year >= 1710) # not sure why this still isn't working with creating first and last year vectors

raw <- raw %>%
  filter(year >= 1000)
}

  # make yij core x year table
y_ij <- raw %>%
  group_by(study, noaa_id, sp_code, core) %>%
  dplyr::select(-rwl_norm, -core1) %>%
  tidyr::spread(key = year, value = rwl)

# make core-species table
core_df <- raw %>%
  group_by(study, noaa_id, sp_code, core) %>%
  dplyr::select(study, noaa_id, sp_code, core) %>%
  distinct() %>%
  mutate(tree = str_sub(core, end = -2), # separate cores from the same core -isn't going to work for trees with only 1 core
         tree_core = str_sub(core, end = 1)) %>%
  mutate(tree = if_else(purrr::is_character(tree_core), tree, core)) # assume any cores from the same tree have identical IDs except ending with a letter

unique(core_df$tree_core)
unique(core_df$tree)

# Climate
load(file = "Data/az_nm/jan_jul_ppt_cm.RData") # jan-jul mean precip
climate2 <- climate2 %>%
  dplyr::filter(year <= max(raw$year, na.rm = TRUE))

years <- sort(unique(raw$year))

# standardize climate data
x_obs <- as.numeric(climate2$precip)
x_full <- c(rep(NA_real_, times = length(years)-length(x_obs)), x_obs) # augment data with NA for estimation

x_mean <- mean(x_full, na.rm=TRUE)
x_sd <- sd(x_full, na.rm=TRUE);
x_s <- (x_full - x_mean) / x_sd # standardize temperatures to z-scores for numerical stability
x_use <- x_s

# Hold out years for validation (first half)
climate_rows <- which(!is.na(x_full))  # the indices that we hold out
hold_out <- climate_rows[1:floor(length(climate_rows) / 2)]
x_use[hold_out] <- NA
cal_ids <- (max(hold_out)+1):(max(climate_rows))

years[hold_out]

# Get tree ring data into tree x year matrix
y_orig <- y_ij %>%
  ungroup() %>%
  dplyr::select(-study, -noaa_id, -sp_code, -core) %>%
  as.data.frame() %>%
  as.matrix()

# For now for the sake of time and memory just do last 500 years of climate
# y_orig <- y_orig[ , (ncol(y_orig)-1000):ncol(y_orig)]
  
y_std <- (y_orig - mean(y_orig, na.rm = T) / sd(y_orig, na.rm = T))
log_y <- log(y_orig + 0.001) # add tiny increase for years with zero growth (missing rings?)

M <- nrow(y_orig)
Tea <- ncol(y_orig)

l <- f <- rep(NA, M) # set up blank first and last year vectors for each tree
for(i in 1:M) {
  tmp <- which(!is.na(y_orig[i, ])) # column (year) with first non-NA
  f[i] <- min(tmp)
  l[i] <- max(tmp)
}

# Age to use
# Can't know age because of growth form of these trees so can't do RCS but can just do individual standardizations with year index as age
age <- matrix(NA, nrow = nrow(y_orig), ncol = ncol(y_orig))
for(i in 1:M) {
  columns <- f[i]:l[i]
  ages <- 1:length(columns)
  age[i, columns] <- ages
}
a_use <- (age - mean(age, na.rm = TRUE)) / sd(age, na.rm = TRUE)

K <- length(unique(y_ij$sp_code))

# testing using tornetrask from other script
if(FALSE) {
  K <- 1
  species <- rep(1, times = M)
}

species <- as.integer(as.factor(y_ij$sp_code))

save(climate2, years, M, Tea, log_y, x_full, a_use, K, cores, x_mean, x_sd, x_use, species, nb, ni, nt, hold_out, cal_ids, file = "Results/NM/JAGS/model_prep_nm.RData")

########## Run Models ###########

##### negexp detrend linear climate M2-non-centered ####
initialize_m2_nc = function(){
  alpha0 = rnorm(M, rowMeans(log_y, na.rm=TRUE), 0.25)
  alpha1 = -rlnorm(M, -1, 0.25)
  sd_y = rlnorm(M, 0, 1) 
  mu_a0 = mean(alpha0)
  mu_a1 = mean(alpha1)
  sd_a0 = sd(alpha0)
  sd_a1 = sd(alpha1)
  eta = matrix(rnorm(Tea*K, 0, 0.25), Tea, K)
  # for(k in 1:K) {
  #   sd_eta[k] = sd(eta[ , k])
  # }
  sd_eta = runif(K, 0.2, 0.3)
  beta0 = rnorm(K, 0.5, 0.1)
  sd_x = rlnorm(1, 0, 1) 
  return(list(sd_y = sd_y,
              alpha0 = alpha0,
              alpha1 = alpha1,
              # mu_a0 = mu_a0,
              # mu_a1 = mu_a1,
              # sd_a0 = sd_a0,
              # sd_a1 = sd_a1,
              sd_eta = sd_eta,
              beta0 = beta0,
              sd_x = sd_x))
}

m2_nc_data <- list(y = log_y, 
                   f = f, 
                   l = l, 
                   M = M, 
                   Tea = Tea, 
                   a = a_use,
                   K = K,
                  species = species,
                   # v = 410, # or 320 
                   x = x_use)

params <- c("x", "beta0", "sd_eta", "sd_x", "sd_y")

cl <- makeCluster(nc)                       # Request # cores
clusterExport(cl, c("m2_nc_data", "initialize_m2_nc", "params", "M", "f", "l", "K", "Tea", "a_use", "log_y", "x_use", "nb", "ni", "nt", "K", "species"))
clusterSetRNGStream(cl = cl, 8675301)
system.time({ 
  out <- clusterEvalQ(cl, {
    library(rjags)
    jm <- jags.model("Code/JAGS/m2_nc_multisp.txt", m2_nc_data, initialize_m2_nc, n.adapt=nb, n.chains=1) # Compile model and run burnin
    out <- coda.samples(jm, params, n.iter=ni, thin=nt) # Sample from posterior distribution
    return(as.mcmc(out))
  })
}) # 

stopCluster(cl)

# Results
m_negexp_norm <- mcmc.list(out)
save(m_negexp_norm, file = "Results/NM/JAGS/negexp_norm.RData")

plot(m_negexp_norm[ ,c("sd_eta[1]", "sd_eta[2]", "sd_eta[3]", "sd_x", "beta0[1]", "beta0[2]", "beta0[3]")])
par(mfrow = c(1,1))

gelman.diag(m_negexp_norm[ ,c("sd_eta[1]", "sd_eta[2]", "sd_eta[3]", "sd_x", "beta0[1]", "beta0[2]", "beta0[3]")])

summary(m_negexp_norm[ ,c("sd_eta[1]", "sd_eta[2]", "sd_eta[3]", "sd_x", "beta0[1]", "beta0[2]", "beta0[3]")])

x_id_50 = which(substr(varnames(m_negexp_norm),1,2)=="x[") # finds the indices of the x variables
post_climate_50 = colMeans(as.matrix(m_negexp_norm[,x_id_50])) # finds the posterior mean of the x variables
# plot(post_climate_50,type="l") # plots the posterior mean of the x variables
plot(years, post_climate_50*x_sd + x_mean,type="l") # plots the posterior mean of the x variables on the original scale (degrees C) and year on x-axis

# validation
x_valid_post_m <- post_climate_50[hold_out]*x_sd + x_mean
x_valid <- x_full[hold_out]
plot(x_valid, x_valid_post_m, type = "p")
abline(0, 1, col = "red")
cor(x_valid, x_valid_post_m)

# Validation Performance Statistics

# performance statistics just using median of posterior prob for valid data
performance_df <- perf_stats(est_valid = x_valid_post_m, observed = x_full, valid_id = hold_out, cal_id = cal_ids, mod_id = "negexp_linear")

# reconstruction plot
recon <- plot_recon(m_negexp_norm, obs = data.frame(year = years, value = x_full) , y_lab = "Mean Jan-Jul Precip. (cm)")

recon2 <- recon + theme_bw_poster()
ggsave(plot = recon2, filename = "Results/Figures/NM/negexp_norm_poster.pdf", dpi = 300)

recon2 <- recon + theme_bw_journal()
ggsave(plot = recon2, filename = "Results/Figures/NM/negexp_norm_paper.pdf", dpi = 1000)

# consider doing observed values as points to better see where they fall within the credible interval
recon_valid <- plot_recon(m_negexp_norm, obs = data.frame(year = years, value = x_full) , y_lab = "Mean Jan-Jul Precip. (cm)", valid_yrs = years[hold_out])
recon_valid2 <- recon_valid + theme_bw_poster()
ggsave(plot = recon_valid2, filename = "Results/Figures/NM/negexp_norm_poster_valid.pdf", dpi = 300)

rm(out)
rm(m_negexp_norm)

## any better than random points around the mean?


# recon + geom_smooth(se = FALSE) # do not do this with the MCMC chains. Maybe smooth through the mean or something else because this will kill the computer.

#####


##### negexp detrend linear climate without holdout ####
initialize_m2_nc = function(){
  alpha0 = rnorm(M, rowMeans(log_y, na.rm=TRUE), 0.25)
  alpha1 = -rlnorm(M, -1, 0.25)
  sd_y = rlnorm(M, 0, 1) 
  mu_a0 = mean(alpha0)
  mu_a1 = mean(alpha1)
  sd_a0 = sd(alpha0)
  sd_a1 = sd(alpha1)
  eta = matrix(rnorm(Tea*K, 0, 0.25), Tea, K)
  # for(k in 1:K) {
  #   sd_eta[k] = sd(eta[ , k])
  # }
  sd_eta = runif(K, 0.2, 0.3)
  beta0 = rnorm(K, 0.5, 0.1)
  sd_x = rlnorm(1, 0, 1) 
  return(list(sd_y = sd_y,
              alpha0 = alpha0,
              alpha1 = alpha1,
              # mu_a0 = mu_a0,
              # mu_a1 = mu_a1,
              # sd_a0 = sd_a0,
              # sd_a1 = sd_a1,
              sd_eta = sd_eta,
              beta0 = beta0,
              sd_x = sd_x))
}

m2_nc_data <- list(y = log_y, 
                   f = f, 
                   l = l, 
                   M = M, 
                   Tea = Tea, 
                   a = a_use,
                   K = K,
                   species = species,
                   # v = 410, # or 320 
                   x = x_s)

params <- c("x", "beta0", "sd_eta", "sd_x", "sd_y")

cl <- makeCluster(nc)                       # Request # cores
clusterExport(cl, c("m2_nc_data", "initialize_m2_nc", "params", "M", "f", "l", "K", "Tea", "a_use", "log_y", "x_use", "nb", "ni", "nt", "K", "species"))
clusterSetRNGStream(cl = cl, 8675301)
system.time({ 
  out <- clusterEvalQ(cl, {
    library(rjags)
    jm <- jags.model("Code/JAGS/m2_nc_multisp.txt", m2_nc_data, initialize_m2_nc, n.adapt=nb, n.chains=1) # Compile model and run burnin
    out <- coda.samples(jm, params, n.iter=ni, thin=nt) # Sample from posterior distribution
    return(as.mcmc(out))
  })
}) # 

stopCluster(cl)

# Results
m_negexp_norm_full <- mcmc.list(out)
save(m_negexp_norm_full, file = "Results/NM/JAGS/negexp_norm_full.RData")

plot(m_negexp_norm_full[ ,c("sd_eta[1]", "sd_eta[2]", "sd_eta[3]", "sd_x", "beta0[1]", "beta0[2]", "beta0[3]")])
par(mfrow = c(1,1))

gelman.diag(m_negexp_norm_full[ ,c("sd_eta[1]", "sd_eta[2]", "sd_eta[3]", "sd_x", "beta0[1]", "beta0[2]", "beta0[3]")])

summary(m_negexp_norm_full[ ,c("sd_eta[1]", "sd_eta[2]", "sd_eta[3]", "sd_x", "beta0[1]", "beta0[2]", "beta0[3]")])

x_id_50 = which(substr(varnames(m_negexp_norm_full),1,2)=="x[") # finds the indices of the x variables
post_climate_50 = colMeans(as.matrix(m_negexp_norm_full[,x_id_50])) # finds the posterior mean of the x variables
# plot(post_climate_50,type="l") # plots the posterior mean of the x variables
plot(years, post_climate_50*x_sd + x_mean,type="l") # plots the posterior mean of the x variables on the original scale (degrees C) and year on x-axis

# validation
x_valid_post_m <- post_climate_50[hold_out]*x_sd + x_mean
x_valid <- x_full[hold_out]
plot(x_valid, x_valid_post_m, type = "p")
abline(0, 1, col = "red")
cor(x_valid, x_valid_post_m)

# Validation Performance Statistics

# performance statistics just using median of posterior prob for valid data
# performance_df <- perf_stats(est_valid = x_valid_post_m, observed = x_full, valid_id = hold_out, cal_id = cal_ids, mod_id = "negexp_linear")

# reconstruction plot
recon <- plot_recon(m_negexp_norm_full, obs = data.frame(year = years, value = x_full) , y_lab = "Mean Jan-Jul Precip. (cm)")

recon2 <- recon + theme_bw_poster()
ggsave(plot = recon2, filename = "Results/Figures/NM/negexp_norm_full_poster.pdf", dpi = 300)

recon2 <- recon + theme_bw_journal()
ggsave(plot = recon2, filename = "Results/Figures/NM/negexp_norm_full_paper.pdf", dpi = 1000)

# consider doing observed values as points to better see where they fall within the credible interval
recon_valid <- plot_recon(m_negexp_norm_full, obs = data.frame(year = years, value = x_full) , y_lab = "Mean Jan-Jul Precip. (cm)", valid_yrs = years[hold_out])
recon_valid2 <- recon_valid + theme_bw_poster()
ggsave(plot = recon_valid2, filename = "Results/Figures/NM/negexp_norm_full_poster_valid.pdf", dpi = 300)

rm(out)
rm(m_negexp_norm_full)

## any better than random points around the mean?


# recon + geom_smooth(se = FALSE) # do not do this with the MCMC chains. Maybe smooth through the mean or something else because this will kill the computer.

#####


##### negexp detrend 1 changepoint climate ####

# x_min <- (0 - x_mean) / x_sd # value of x on standardized scale when climate = 0
x_min <- min(x_use, na.rm = TRUE)
x_max <- max(x_use, na.rm = TRUE)

initialize_m2_nc = function(){
  alpha0 = rnorm(M, rowMeans(log_y, na.rm=TRUE), 0.25)
  alpha1 = -rlnorm(M, -1, 0.25)
  sd_y = rlnorm(M, 0, 1) 
  mu_a0 = mean(alpha0)
  mu_a1 = mean(alpha1)
  sd_a0 = sd(alpha0)
  sd_a1 = sd(alpha1)
  x_1 = runif(K, -1.1, -0.5)
  eta = matrix(rnorm(Tea*K, 0, 0.25), Tea, K)
  beta0 <- matrix(NA, nrow = K, ncol = 2)
  for(k in 1:K) {
    beta0[k, 1] = rnorm(1, 1, 0.1)
    beta0[k, 2] = rnorm(1, 0.5, 0.1)
  #   sd_eta[k] = sd(eta[ , k])
  }
  sd_eta = runif(K, 0.2, 0.3)
  sd_x = rlnorm(1, 0, 1) 
  return(list(sd_y = sd_y,
              alpha0 = alpha0,
              alpha1 = alpha1,
              # mu_a0 = mu_a0,
              # mu_a1 = mu_a1,
              # sd_a0 = sd_a0,
              # sd_a1 = sd_a1,
              sd_eta = sd_eta,
              x_1 = x_1,
              beta0 = beta0,
              sd_x = sd_x))
}

m2_nc_data <- list(y = log_y, 
                   f = f, 
                   l = l, 
                   M = M, 
                   Tea = Tea, 
                   a = a_use,
                   K = K,
                   species = species,
                   x_min = x_min,
                   x_max = x_max,
                   # v = 410, # or 320 
                   x = x_use)

params <- c("x", "beta0", "sd_eta", "sd_x", "sd_y", "x_1")

cl <- makeCluster(nc)                       # Request # cores
clusterExport(cl, c("m2_nc_data", "initialize_m2_nc", "params", "M", "f", "l", "K", "Tea", "a_use", "log_y", "x_use", "x_min", "nb", "ni", "nt", "K", "species"))
clusterSetRNGStream(cl = cl, 8675301)
system.time({ 
  out <- clusterEvalQ(cl, {
    library(rjags)
    jm <- jags.model("Code/JAGS/negexp_1changept_multispp.txt", m2_nc_data, initialize_m2_nc, n.adapt=nb, n.chains=1) # Compile model and run burnin
    out <- coda.samples(jm, params, n.iter=ni, thin=nt) # Sample from posterior distribution
    return(as.mcmc(out))
  })
}) # 

stopCluster(cl)

# Results
m_negexp_1change <- mcmc.list(out)
save(m_negexp_1change, file = "Results/NM/JAGS/negexp_1change.RData")

m_negexp_1change <- mcmc.list(m_negexp_1change[[2]], m_negexp_1change[[3]], m_negexp_1change[[4]])

plot(m_negexp_1change[ ,c("sd_eta[1]", "sd_eta[2]", "sd_eta[3]", "sd_x", "beta0[1,1]", "beta0[2,1]", "beta0[3,1]", "beta0[1,2]", "beta0[2,2]", "beta0[3,2]", "x_1[1]", "x_1[2]", "x_1[3]")])
par(mfrow = c(1,1))

gelman.diag(m_negexp_1change[ ,c("sd_eta[1]", "sd_eta[2]", "sd_eta[3]", "sd_x", "beta0[1,1]", "beta0[2,1]", "beta0[3,1]", "beta0[1,2]", "beta0[2,2]", "beta0[3,2]", "x_1[1]", "x_1[2]", "x_1[3]")])

effectiveSize(m_negexp_1change[ ,c("sd_eta[1]", "sd_eta[2]", "sd_eta[3]", "sd_x", "beta0[1,1]", "beta0[2,1]", "beta0[3,1]", "beta0[1,2]", "beta0[2,2]", "beta0[3,2]", "x_1[1]", "x_1[2]", "x_1[3]")])

x_id_50 = which(substr(varnames(m_negexp_1change),1,2)=="x[") # finds the indices of the x variables
x_id_est <- x_id_50[!(x_id_50 %in% years[cal_ids])] # not working

N_eff_x <- effectiveSize(as.matrix(m_negexp_1change[ , x_id_est])) #
range(N_eff_x[which(N_eff_x > 0)])
hist(N_eff_x[which(N_eff_x > 0)])

post_climate_50 = colMeans(as.matrix(m_negexp_1change[ , x_id_50])) # finds the posterior mean of the x variables
# plot(post_climate_50,type="l") # plots the posterior mean of the x variables
plot(years, post_climate_50*x_sd + x_mean,type="l") # plots the posterior mean of the x variables on the original scale (degrees C) and year on x-axis

# validation
x_valid_post_m <- post_climate_50[hold_out]*x_sd + x_mean
x_valid <- x_full[hold_out]
plot(x_valid, x_valid_post_m, type = "p")
abline(0, 1, col = "red")
cor(x_valid, x_valid_post_m)

# Validation Performance Statistics

# performance statistics just using median of posterior prob for valid data
performance_df[2, ] <- perf_stats(est_valid = x_valid_post_m, observed = x_full, valid_id = hold_out, cal_id = cal_ids, mod_id = "negexp_1change")

# reconstruction plot
recon <- plot_recon(m_negexp_1change, obs = data.frame(year = years, value = x_full) , y_lab = "Mean Jan-Jul Precip. (cm)")

# consider doing observed values as points to better see where they fall within the credible interval
recon_valid <- plot_recon(m_negexp_1change, obs = data.frame(year = years, value = x_full) , y_lab = "Mean Jan-Jul Precip. (cm)", valid_yrs = years[hold_out])

# reconstruction plot
recon <- plot_recon(m_negexp_1change, obs = data.frame(year = years, value = x_full) , y_lab = "Mean Jan-Jul Precip. (cm)")

recon2 <- recon + theme_bw_poster()
ggsave(plot = recon2, filename = "Results/Figures/NM/negexp_1change_poster.pdf", dpi = 300)

recon2 <- recon + theme_bw_journal()
ggsave(plot = recon2, filename = "Results/Figures/NM/negexp_1change_paper.pdf", dpi = 1000)

rm(out)
rm(m_negexp_1change)

#####

###### No Detrending Mean Climate #####

x_min <- (0 - x_mean) / x_sd # value of x on standardized scale when climate = 0

initialize_m2_nc = function(){
  alpha0 = rnorm(M, rowMeans(log_y, na.rm=TRUE), 0.25)
 # alpha1 = -rlnorm(M, -1, 1)
  sd_y = rlnorm(M, 0, 0.25) 
  sd_a0 = sd(alpha0)
  # sd_a1 = sd(alpha1)
  eta = matrix(rnorm(Tea*K, 0, 1), Tea, K)
  beta0 <- matrix(NA, nrow = K, ncol = 2)
  for(k in 1:K) {
    beta0[k, 1] = 0
    beta0[k, 2] = rnorm(1, 0.5, 0.1)
    #   sd_eta[k] = sd(eta[ , k])
  }
  sd_eta = runif(K, 0.2, 0.3)
  sd_x = rlnorm(1, 0, 0.25) 
  return(list(sd_y = sd_y,
              alpha0 = alpha0,
              # alpha1 = alpha1,
              # mu_a0 = mu_a0,
              # mu_a1 = mu_a1,
              # sd_a0 = sd_a0,
              # sd_a1 = sd_a1,
              sd_eta = sd_eta,
              beta0 = beta0,
              sd_x = sd_x))
}

m2_nc_data <- list(y = log_y, 
                   f = f, 
                   l = l, 
                   M = M, 
                   Tea = Tea, 
                   # a = a_use,
                   K = K,
                   species = species,
                   x_min = x_min,
                   # v = 410, # or 320 
                   x = x_use)

params <- c("x", "beta0", "sd_eta", "sd_x", "sd_y", "x_1")

cl <- makeCluster(nc)                       # Request # cores
clusterExport(cl, c("m2_nc_data", "initialize_m2_nc", "params", "M", "f", "l", "K", "Tea", "log_y", "x_use", "x_min", "nb", "ni", "nt", "K", "species"))
clusterSetRNGStream(cl = cl, 9675303)
system.time({ 
  out <- clusterEvalQ(cl, {
    library(rjags)
    jm <- jags.model("Code/JAGS/no_detrend_1changept_multispp.txt", m2_nc_data, initialize_m2_nc, n.adapt=nb, n.chains=1) # Compile model and run burnin
    out <- coda.samples(jm, params, n.iter=ni, thin=nt) # Sample from posterior distribution
    return(as.mcmc(out))
  })
}) # 

stopCluster(cl)

# Results
m_none_1change <- mcmc.list(out)
save(m_none_1change, file = "Results/NM/JAGS/none_1change.RData")

plot(m_none_1change[ ,c("sd_eta[1]", "sd_eta[2]", "sd_eta[3]", "sd_x", "beta0[1,1]", "beta0[2,1]", "beta0[3,1]", "beta0[1,2]", "beta0[2,2]", "beta0[3,2]", "x_1")])
par(mfrow = c(1,1))

summary(m_none_1change[ ,c("sd_eta[1]", "sd_eta[2]", "sd_eta[3]", "sd_x", "beta0[1,1]", "beta0[2,1]", "beta0[3,1]", "beta0[1,2]", "beta0[2,2]", "beta0[3,2]", "x_1")])

gelman.diag(m_none_1change[ ,c("sd_eta[1]", "sd_eta[2]", "sd_eta[3]", "sd_x", "beta0[1,1]", "beta0[2,1]", "beta0[3,1]", "beta0[1,2]", "beta0[2,2]", "beta0[3,2]", "x_1")])

x_id_50 = which(substr(varnames(m_none_1change),1,2)=="x[") 
post_climate_50 = colMeans(as.matrix(m_none_1change[,x_id_50]))
x_valid_post_m <- post_climate_50[hold_out]*x_sd + x_mean
x_valid <- x_full[hold_out]
plot(x_valid, x_valid_post_m, type = "p")
abline(0, 1, col = "red")
cor(x_valid, x_valid_post_m)

# Validation Performance Statistics

# performance statistics just using median of posterior prob for valid data
performance_df[3, ] <- perf_stats(est_valid = x_valid_post_m, observed = x_full, valid_id = hold_out, cal_id = cal_ids, mod_id = "mean_1change")

rm(out)
rm(m_none_1change)

#####

#### climate spline 50 year knots - not getting reasonable results w/25 or 50 year knots #####

knots <- seq(0, Tea, by = 25)
B <- bs(1:Tea, knots = knots) # neet to sort x_use? then doesn't line up with rest.
H <- length(knots) + 3

# run model
initialize_spl = function(){
  alpha0 = rnorm(M, rowMeans(log_y, na.rm=TRUE), 0.25)
  alpha1 = -rlnorm(M, -1, 0.25)
  sd_y = rlnorm(M, 0, 1) 
  mu_a0 = mean(alpha0)
  mu_a1 = mean(alpha1)
  sd_a0 = sd(alpha0)
  sd_a1 = sd(alpha1)
  mu_gamma = rnorm(1, 0, 0.2)
  sd_g = runif(1, 0.1, 0.3)
  eta = matrix(rnorm(Tea*K, 0, 0.25), Tea, K)
  sd_eta = runif(K, 0.2, 0.3) 
  beta0 = rnorm(K, 0.5, 0.1)
  sd_x = rlnorm(1, 0, 1) 
  return(list(sd_y = sd_y,
              alpha0 = alpha0,
              alpha1 = alpha1,
              # mu_a0 = mu_a0,
              # mu_a1 = mu_a1,
              # sd_a0 = sd_a0,
              # sd_a1 = sd_a1,
              mu_gamma = mu_gamma,
              sd_g = sd_g,
              sd_eta = sd_eta,
              beta0 = beta0,
              sd_x = sd_x))
}

spl_data <- list(y = log_y, 
                 f = f, 
                 l = l, 
                 M = M, 
                 Tea = Tea, 
                 a = a_use, 
                 x = x_use,
                 B = B,
                 H = H,
                 K = K,
                 species = species)

params <- c("x", "sd_eta", "sd_x", "sd_y", "mu_gamma", "sd_g")

cl <- makeCluster(nc)                       # Request # cores
clusterExport(cl, c("spl_data", "initialize_spl", "params", "M", "f", "l", "Tea", "a_use", "log_y", "x_use", "nb", "ni", "nt", "B", "K", "species"))
clusterSetRNGStream(cl = cl, 8675302)
system.time({ 
  out <- clusterEvalQ(cl, {
    library(rjags)
    jm <- jags.model("Code/JAGS/climate_spline_multisp.txt", spl_data, initialize_spl, n.adapt=nb, n.chains=1) # Compile model and run burnin
    out <- coda.samples(jm, params, n.iter=ni, thin=nt) # Sample from posterior distribution
    return(as.mcmc(out))
  })
}) # 
stopCluster(cl)

# Results
m_spline_25 <- mcmc.list(out)
save(m_spline_25, file = "Results/NM/JAGS/negexp_spline25.RData")

plot(m_spline_25[ ,c("sd_eta[1]", "sd_eta[2]", "sd_eta[3]", "sd_x", "mu_gamma", "sd_g")])
par(mfrow = c(1,1))

x_id_25 = which(substr(varnames(m_spline_25),1,2)=="x[") 
post_climate_25 = colMeans(as.matrix(m_spline_25[,x_id_25]))
x_valid_post_m <- post_climate_25[hold_out]*x_sd + x_mean
x_valid <- x_full[hold_out]
plot(x_valid, x_valid_post_m, type = "p")
abline(0, 1, col = "red")
cor(x_valid, x_valid_post_m)

# Validation Performance Statistics

# performance statistics just using median of posterior prob for valid data
performance_df[4, ] <- perf_stats(est_valid = x_valid_post_m, observed = x_full, valid_id = hold_out, cal_id = cal_ids, mod_id = "negexp_spline25")

# reconstruction plots
recon <- plot_recon(m_spline_25)
ggsave(filename = "Results/Figures/torn_recon_post_spl25.tiff", plot = recon, width = 8, height = 4, units = "in") # , dpi = 1000) # poster vs paper formatting - see past work and make package or github source

recon_valid <- plot_recon(m_spline_25, valid_yrs = year[hold_out])
ggsave(filename = "Results/Figures/torn_recon_spl25_valid_back.tiff", plot = recon_valid, width = 8, height = 4, units = "in") 


# reconstruction plot
recon <- plot_recon(m_spline_25, obs = data.frame(year = years, value = x_full) , y_lab = "Mean Jan-Jul Precip. (cm)")

recon2 <- recon + theme_bw_poster()
ggsave(plot = recon2, filename = "Results/Figures/NM/negexp_spline25_poster.pdf", dpi = 300)

recon2 <- recon + theme_bw_journal()
ggsave(plot = recon2, filename = "Results/Figures/NM/negexp_spline25_paper.pdf", dpi = 1000)

rm(out)
rm(m_spline_25)

write.csv(performance_df, file = "Results/NM/nm_performance_stats.csv", row.names = FALSE)

#####




####### pre-whitening and re-reddening - AR model? #####

##### Extra ######

if(FALSE) {
  
  ##### linear detrend 1 changepoint climate - linear error - not working and probably not a great model ####
  
  x_min <- (0 - x_mean) / x_sd # value of x on standardized scale when climate = 0
  
  initialize_m2_nc = function(){
    alpha0 = rnorm(M, rowMeans(y_std, na.rm=TRUE), 0.5)
    alpha1 = -rlnorm(M, -1, 1)
    sd_y = rlnorm(M, 0, 0.5) 
    mu_a0 = mean(alpha0)
    mu_a1 = mean(alpha1)
    sd_a0 = sd(alpha0)
    sd_a1 = sd(alpha1)
    eta = matrix(rnorm(Tea*K, 0, 1), Tea, K)
    beta0 <- matrix(NA, nrow = K, ncol = 2)
    for(k in 1:K) {
      beta0[k, 1] = 0
      beta0[k, 2] = rnorm(1, 0.5, 0.1)
      #   sd_eta[k] = sd(eta[ , k])
    }
    sd_eta = runif(K, 0.2, 0.3)
    sd_x = rlnorm(1, 0, 1) 
    return(list(sd_y = sd_y,
                alpha0 = alpha0,
                alpha1 = alpha1,
                # mu_a0 = mu_a0,
                # mu_a1 = mu_a1,
                # sd_a0 = sd_a0,
                # sd_a1 = sd_a1,
                sd_eta = sd_eta,
                beta0 = beta0,
                sd_x = sd_x))
  }
  
  m2_nc_data <- list(y = y_std, 
                     f = f, 
                     l = l, 
                     M = M, 
                     Tea = Tea, 
                     a = a_use,
                     K = K,
                     species = species,
                     x_min = x_min,
                     # v = 410, # or 320 
                     x = x_use)
  
  params <- c("x", "beta0", "sd_eta", "sd_x", "sd_y", "x_1")
  
  cl <- makeCluster(nc)                       # Request # cores
  clusterExport(cl, c("m2_nc_data", "initialize_m2_nc", "params", "M", "f", "l", "K", "Tea", "a_use", "y_std", "x_use", "x_min", "nb", "ni", "nt", "K", "species"))
  clusterSetRNGStream(cl = cl, 8675301)
  system.time({ 
    out <- clusterEvalQ(cl, {
      library(rjags)
      jm <- jags.model("Code/JAGS/negexp_1changept_multispp.txt", m2_nc_data, initialize_m2_nc, n.adapt=nb, n.chains=1) # Compile model and run burnin
      out <- coda.samples(jm, params, n.iter=ni, thin=nt) # Sample from posterior distribution
      return(as.mcmc(out))
    })
  }) # 
  
  stopCluster(cl)
  
  # Results
  m_linear_1change <- mcmc.list(out)
  
  plot(m_linear_1change[ ,c("sd_eta[1]", "sd_eta[2]", "sd_eta[3]", "sd_x", "beta0[1,1]", "beta0[2,1]", "beta0[3,1]", "beta0[1,2]", "beta0[2,2]", "beta0[3,2]", "x_1")])
  par(mfrow = c(1,1))
  x_id_50 = which(substr(varnames(m_negexp_1change),1,2)=="x[") # finds the indices of the x variables
  post_climate_50 = colMeans(as.matrix(m_negexp_1change[,x_id_50])) # finds the posterior mean of the x variables
  # plot(post_climate_50,type="l") # plots the posterior mean of the x variables
  plot(years, post_climate_50*x_sd + x_mean,type="l") # plots the posterior mean of the x variables on the original scale (degrees C) and year on x-axis
  
  # validation
  x_valid_post_m <- post_climate_50[hold_out]*x_sd + x_mean
  x_valid <- x_full[hold_out]
  plot(x_valid, x_valid_post_m, type = "p")
  abline(0, 1, col = "red")
  cor(x_valid, x_valid_post_m)
  
  res2 <- (x_valid_post_m - x_valid)^2
  sqrt(sum(res2) / length(x_valid_post_m))
  sqrt(mean(res2))
  
  # reconstruction plot
  plot_recon <- function(mcmc, obs = climate, mean = x_mean, sd = x_sd, valid_yrs = NULL) {
    xidx_m2 = which(substr(varnames(mcmc),1,2)=="x[")
    temp_df <- as_tibble(t(as.matrix(mcmc[ , xidx_m2])))
    temp_df <- temp_df %>% 
      mutate(year = obs$year)
    
    temp_df_long <- temp_df %>% 
      gather(key = sim, value = temp, -year) %>%
      dplyr::mutate(temp = temp*sd + mean)
    
    # Validation plot
    if(!is.null(valid_yrs)) {
      temp_valid <- temp_df_long %>%
        dplyr::filter(year %in% valid_yrs) %>%
        dplyr::mutate(Value = "estimated")
      
      climate_valid <- obs %>%
        dplyr::mutate(Value = "observed") %>%
        dplyr::rename(temp = value) %>% # x_full) %>%
        dplyr::filter(year %in% valid_yrs)
      
      g <- ggplot(temp_valid, aes(x = year, y = temp))
      g <- g + geom_fan() + geom_line(data = climate_valid, aes(year, temp), colour = "red") + ylab("Mean Jan-Jul Precip. (cm)") + xlab("Year") + scale_fill_distiller() + theme_bw()
    } else {
      # scaled posterior interval
      g <- ggplot(temp_df_long, aes(x = year, y = temp))
      g <- g + geom_fan() + geom_line(data = obs, aes(x = year, y = value), colour="black", size = 0.2) + ylab("Mean Jan-Jul Precip. (cm)") + xlab("Year") + theme_bw() + scale_fill_distiller() #x_full),
    }
    return(g)
  }
  
  recon <- plot_recon(m_negexp_norm, obs = data.frame(year = years, value = x_full) , y_lab = "Mean Jan-Jul Precip. (cm)")
  
  # consider doing observed values as points to better see where they fall within the credible interval
  recon_valid <- plot_recon(m_negexp_norm, obs = data.frame(year = years, value = x_full) , y_lab = "Mean Jan-Jul Precip. (cm)", valid_yrs = years[hold_out])
  
  ## any better than random points around the mean?
  
  
  # recon + geom_smooth(se = FALSE) # do not do this with the MCMC chains. Maybe smooth through the mean or something else because this will kill the computer.
  
  #####
  
  
  ##### negexp detrend with power law climate ####
  
  # negExp (log-linear)
  m2_nc_data <- list(y = log_y, # log_y, # try with linear scale first
                     f = f, 
                     l = l, 
                     M = M, 
                     Tea = Tea, 
                     a = a_use,
                     K = K,
                     species = species,
                     # v = 410, # or 320 
                     x = log(x_full)) # power law with decimals needs positive values of X
  
  initialize_m2_nc = function(){
    alpha0 = rnorm(M, rowMeans(log_y, na.rm=TRUE), 0.25) # log_y
    alpha1 = -rlnorm(M, -1, 0.25)
    sd_y = rlnorm(M, 0, 1) 
    mu_a0 = mean(alpha0)
    mu_a1 = mean(alpha1)
    sd_a0 = sd(alpha0)
    sd_a1 = sd(alpha1)
    eta = matrix(rnorm(Tea*K, 0, 0.25), Tea, K)
    # for(k in 1:K) {
    #   sd_eta[k] = sd(eta[ , k])
    # }
    sd_eta = runif(K, 0.2, 0.3)
    beta0 = rnorm(K, 0.5, 0.1)
    theta = runif(K, 0.3, 0.7)
    sd_x = rlnorm(1, 0, 1) 
    return(list(sd_y = sd_y,
                alpha0 = alpha0,
                alpha1 = alpha1,
                # mu_a0 = mu_a0,
                # mu_a1 = mu_a1,
                # sd_a0 = sd_a0,
                # sd_a1 = sd_a1,
                sd_eta = sd_eta,
                beta0 = beta0,
                theta = theta,
                sd_x = sd_x))
  }
  
  # Linear
  
  library(nlme)
  lx <- log(x_full)
  df1 <- data.frame(y1 = as.numeric(y_orig[1, ]), y20 = as.numeric(y_orig[2, ]), log_y = log(y_orig[1, ]+0.001), lx = lx) %>%
    drop_na()
  
  lm1 <- lm(y1 ~ 1 + lx, data = df1)
  summary(lm1)
  plot(fitted(lm1), df1$y1)
  abline(0, 1)
  
  pl1 <- nls(y1 ~ k + a * lx ^ b, df1) #, start=c(a=1, b=1, k=5))
  summary(pl1)
  plot(pl1)
  plot(fitted(pl1), df1$y)
  abline(0, 1)
  
  plot(exp(df1$lx), fitted(pl1))
  
  pl20 <- nls(y20 ~ k + a * lx ^ b, df1) #, start=c(a=1, b=1, k=5))
  summary(pl20)
  
  pl1_log <- nls(log_y ~ k + a * lx ^ b, df1) #, start=c(a=1, b=1, k=5))
  summary(pl1)
  plot(pl1)
  plot(fitted(pl1), df1$y)
  abline(0, 1)
  
  m2_nc_data <- list(y = y_orig, # log_y, # try with linear scale first
                     f = f, 
                     l = l, 
                     M = M, 
                     Tea = Tea, 
                     a = a_use,
                     K = K,
                     species = species,
                     # v = 410, # or 320 
                     x = log(x_full)) # power law with decimals needs positive values of X
  
  initialize_m2_nc = function(){
    alpha0 = rnorm(M, rowMeans(log_y, na.rm=TRUE), 0.25) # log_y
    alpha1 = runif(M, 0, 0.1)
    sd_y = rlnorm(M, 0, 1) 
    mu_a0 = mean(alpha0)
    mu_a1 = mean(alpha1)
    sd_a0 = sd(alpha0)
    sd_a1 = sd(alpha1)
    eta = matrix(rnorm(Tea*K, 0, 0.25), Tea, K)
    # for(k in 1:K) {
    #   sd_eta[k] = sd(eta[ , k])
    # }
    sd_eta = runif(K, 0.2, 0.3)
    beta0 = rnorm(K, 0.2, 0.05)
    # theta = runif(K, 1, 1.5)
    theta = runif(K, 1, 1.5)
    sd_x = rlnorm(1, log(x_obs), 0.1) 
    # x = c(rep(mean(log(x_full), na.rm = T), times = 1895), log(x_obs))
    return(list(
      alpha0 = alpha0,
      alpha1 = alpha1,
      # mu_a0 = mu_a0,
      # mu_a1 = mu_a1,
      # sd_a0 = sd_a0,
      # sd_a1 = sd_a1,
      sd_eta = sd_eta,
      beta0 = beta0,
      theta = theta,
      # x = x,
      # sd_x = sd_x,
      sd_y = sd_y))
  }
  
  params <- c("x", "beta0", "theta", "sd_eta", "sd_x", "sd_y")
  
  cl <- makeCluster(nc)                       # Request # cores
  clusterExport(cl, c("m2_nc_data", "initialize_m2_nc", "params", "M", "f", "l", "K", "Tea", "a_use", "log_y", "x_use", "nb", "ni", "nt", "K", "species", "x_obs", "x_full", "y_orig"))
  clusterSetRNGStream(cl = cl, 8675301)
  system.time({ 
    out <- clusterEvalQ(cl, {
      library(rjags)
      jm <- jags.model("Code/JAGS/negexp_power_multispp.txt", m2_nc_data, initialize_m2_nc, n.adapt=nb, n.chains=1) # Compile model and run burnin
      out <- coda.samples(jm, params, n.iter=ni, thin=nt) # Sample from posterior distribution
      return(as.mcmc(out))
    })
  }) # 
  
  stopCluster(cl)
  
  # Results
  m_negexp_power <- mcmc.list(out)
  
  plot(m_negexp_power[ ,c("sd_eta", "sd_x", "beta0")])
  plot(m_negexp_power[ ,c("sd_eta[1]", "sd_x", "beta0[1]")])
  plot(m_negexp_power[ ,c("sd_eta[1]", "sd_eta[2]", "sd_eta[3]", "sd_x", "beta0[1]", "beta0[2]", "beta0[3]")])
  
  plot(m_negexp_power[ ,c("sd_eta", "sd_x", "beta0[1]", "beta0[2]", "beta0[3]")])
  par(mfrow = c(1,1))
  x_id_50 = which(substr(varnames(m_negexp_power),1,2)=="x[") # finds the indices of the x variables
  post_climate_50 = colMeans(as.matrix(m_negexp_power[,x_id_50])) # finds the posterior mean of the x variables
  plot(post_climate_50,type="l") # plots the posterior mean of the x variables
  plot(years, exp(post_climate_50),type="l") # plots the posterior mean of the x variables on the original scale (degrees C) and year on x-axis
  
  # validation
  x_valid_post_m <- post_climate_50[hold_out]*x_sd + x_mean
  x_valid <- x_full[hold_out]
  plot(x_valid, x_valid_post_m, type = "p")
  abline(0, 1, col = "red")
  cor(x_valid, x_valid_post_m)
  
  res2 <- (x_valid_post_m - x_valid)^2
  sqrt(sum(res2) / length(x_valid_post_m))
  sqrt(mean(res2))
  
  # reconstruction plot
  plot_recon <- function(mcmc, obs = climate, mean = x_mean, sd = x_sd, valid_yrs = NULL) {
    xidx_m2 = which(substr(varnames(mcmc),1,2)=="x[")
    temp_df <- as_tibble(t(as.matrix(mcmc[ , xidx_m2])))
    temp_df <- temp_df %>% 
      mutate(year = obs$year)
    
    temp_df_long <- temp_df %>% 
      gather(key = sim, value = temp, -year) %>%
      dplyr::mutate(temp = temp*sd + mean)
    
    # Validation plot
    if(!is.null(valid_yrs)) {
      temp_valid <- temp_df_long %>%
        dplyr::filter(year %in% valid_yrs) %>%
        dplyr::mutate(Value = "estimated")
      
      climate_valid <- obs %>%
        dplyr::mutate(Value = "observed") %>%
        dplyr::rename(temp = value) %>% # x_full) %>%
        dplyr::filter(year %in% valid_yrs)
      
      g <- ggplot(temp_valid, aes(x = year, y = temp))
      g <- g + geom_fan() + geom_line(data = climate_valid, aes(year, temp), colour = "red") + ylab("Mean Jan-Jul Precip. (cm)") + xlab("Year") + scale_fill_distiller() + theme_bw()
    } else {
      # scaled posterior interval
      g <- ggplot(temp_df_long, aes(x = year, y = temp))
      g <- g + geom_fan() + geom_line(data = obs, aes(x = year, y = value), colour="black", size = 0.2) + ylab("Mean Jan-Jul Precip. (cm)") + xlab("Year") + theme_bw() + scale_fill_distiller() #x_full),
    }
    return(g)
  }
  
  recon <- plot_recon(m_negexp_norm, obs = data.frame(year = years, value = x_full) , y_lab = "Mean Jan-Jul Precip. (cm)")
  
  # consider doing observed values as points to better see where they fall within the credible interval
  recon_valid <- plot_recon(m_negexp_norm, obs = data.frame(year = years, value = x_full) , y_lab = "Mean Jan-Jul Precip. (cm)", valid_yrs = years[hold_out])
  
  ## any better than random points around the mean?
  
  
  # recon + geom_smooth(se = FALSE) # do not do this with the MCMC chains. Maybe smooth through the mean or something else because this will kill the computer.
  
  #####
  
  ###### tree response #####
  # tree response to precip nonlinear - logistic? - currently assuming exponential (log-linear)?
  ###########
  
}


