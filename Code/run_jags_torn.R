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

source("Code/functions.R")

testing <- FALSE

#### Load and Prep Data #####
# load data
a_use <- as.matrix(read.table('Data/Tornetrask/age.txt', header = FALSE))
year <- c(read.table("Data/Tornetrask/year.txt"), recursive = TRUE) 
x_full <- c(read.table("Data/Tornetrask/xtorn.txt"), recursive=TRUE) # climate data

# Combine year and climate data
climate_df <- data.frame(year, x_full, stringsAsFactors = FALSE) %>%
  dplyr::rename(value = x_full)
# instramental temperature record from 1816-1995 (179 years)

# standardize climate data
x_mean <- mean(x_full, na.rm=TRUE)
x_sd <- sd(x_full, na.rm=TRUE);
x_s <- (x_full - x_mean) / x_sd # standardize temperatures to z-scores for numerical stability
x_use <- x_s

# Hold out years for validation (first half)
hold_out <- 321:410  # the indices that we hold out
cal_ids <- 411:500
x_use[hold_out] <- NA

# Get tree ring data
y_orig <- as.matrix(read.table('Data/Tornetrask/ytorn.txt', header = FALSE))
log_y <- log(y_orig)

M <- nrow(y_orig)
Tea <- ncol(y_orig)

l <- f <- rep(NA, M) # set up blank first and last year vectors for each tree
for(i in 1:M) {
  tmp <- which(!is.na(y_orig[i, ])) # column (year) with first non-NA
  f[i] <- min(tmp)
  l[i] <- max(tmp)
}

years <- unique(year)

save(climate_df, years, M, Tea, log_y, x_full, a_use, x_sd, x_mean, hold_out, x_use, cal_ids, file = "Results/JAGS/model_prep.RData")

#### Does creating a standard chronology FORCE the climate reconstruction to be (normally) distributed with fluctuations around a mean????? Over flatten??? #####

##### dplR chronology #####

library(dplR)
y_detrend <- detrend(as.data.frame(t(y_orig)), method = "ModNegExp", return.info = TRUE)
y_detrend$model.info$V1$ModNegExp$coefs

y_detrend <- detrend(as.data.frame(t(y_orig)), method = "Spline")

y_crn <- chron(y_detrend, prefix = "HUR", prewhiten = FALSE)
y_crn$year <- year

y_crn_std <- chron(y_detrend, prefix = "HUR", prewhiten = TRUE)
y_crn_std$year <- year

y_rwl <- as.data.frame(t(y_orig))
po <- data.frame(series = colnames(y_rwl), pith.offset = as.integer(1+ colSums(!is.na(y_rwl)))) # 1)
  
y_rcs <- rcs(y_rwl, po = po)

########## Run Models ###########

if(testing) {
  nb = 1000
  ni = 1000
  nt = 1
  nc = 3
} else {
  nb = 10000
  ni = 20000
  nt = 10
  nc = 4
}

if(!dir.exists("Results/JAGS")) dir.create("Results/JAGS", recursive = TRUE)
if(!dir.exists("Results/Figures/JAGS")) dir.create("Results/Figures/JAGS", recursive = TRUE)

##### negexp detrend linear climate M2-non-centered ####

# my code should have the same output as M2 from Schofield I think

initialize_m2_nc = function(){
  alpha0 = rnorm(M, rowMeans(log_y, na.rm=TRUE), 0.25)
  alpha1 = -rlnorm(M, -1, 0.25)
  sd_y = rlnorm(M, 0, 1) 
  mu_a0 = mean(alpha0)
  mu_a1 = mean(alpha1)
  sd_a0 = sd(alpha0)
  sd_a1 = sd(alpha1)
  eta = rnorm(Tea, 0, 0.25)
  sd_eta = sd(eta) 
  beta0 = rnorm(1, 0, 0.1)
  sd_x = rlnorm(1, 0, 1) 
  return(list(sd_y = sd_y,
              mu_a0 = mu_a0,
              mu_a1 = mu_a1,
              sd_a0 = sd_a0,
              sd_a1 = sd_a1,
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
                   # v = 410, # or 320 
                   x = x_use)

params <- c("x", "mu_a0", "mu_a1", "sd_a0", "sd_a1", "beta0", "sd_eta", "sd_x", "sd_y")

cl <- makeCluster(nc)                       # Request # cores
clusterExport(cl, c("m2_nc_data", "initialize_m2_nc", "params", "M", "f", "l", "Tea", "a_use", "log_y", "x_use", "nb", "ni", "nt"))
clusterSetRNGStream(cl = cl, 8675301)
system.time({ 
  out <- clusterEvalQ(cl, {
    library(rjags)
    jm <- jags.model("Code/JAGS/m2_nc.txt", m2_nc_data, initialize_m2_nc, n.adapt=nb, n.chains=1) # Compile model and run burnin
    out <- coda.samples(jm, params, n.iter=ni, thin=nt) # Sample from posterior distribution
    return(as.mcmc(out))
  })
}) # 

stopCluster(cl)

# Results
outM2_nc <- mcmc.list(out)
plot(outM2_nc[ , c("mu_a0", "mu_a1", "sd_a0", "sd_a1", "beta0", "sd_eta", "sd_x")])
par(mfrow = c(1,1))
# summary(outM2_nc[,c("mu_a0", "mu_a1", "sd_a0", "sd_a1", "beta0", "sd_eta", "sd_x")])

save(outM2_nc, file = "Results/JAGS/M2_nc.RData")

xidx_m2 = which(substr(varnames(outM2_nc),1,2)=="x[")
postxm_m2 = colMeans(as.matrix(outM2_nc[,xidx_m2])) 
x_post = as.matrix(outM2_nc[ , xidx_m2])
x_post_median <- apply(x_post, 2, median)

# validation
x_valid_post_m <- postxm_m2[hold_out]*x_sd + x_mean
x_valid <- x_full[hold_out]
plot(x_valid, x_valid_post_m, type = "p")
abline(0, 1, col = "red")
cor(x_valid, x_valid_post_m)

res2 <- (x_valid_post_m - x_valid)^2
sqrt(sum(res2) / length(x_valid_post_m))
sqrt(mean(res2))

#+ geom_line(data = data.frame(year = year, median = x_post_median*x_sd + x_mean), aes(x = year, y = median), colour="black", size = 0.2) + ylab("Temperature (C)") + xlab("Year") + theme_bw() + scale_fill_distiller()

# m2_climate <- outM2_nc[, grep("x[", colnames(m2_nc_data), fixed=T)]
# 
# foo <- outM2_nc %>% 
#   gather_draws(x) %>%
#   head(15)
# 
# m2_climate <- as.matrix(outM2_nc[, c("x")])


rm(out)
rm(outM2_nc)

#####

##### spline on age standardization #####

library(splines)

# penalty matrix function
makeQ = function(degree, K, epsilon=1e-3){
  x <- diag(K)
  E <- diff(x, differences=degree)
  return( t(E) %*% E + x*epsilon)
} 

# Set up basis function for each tree ring series
K = 6 # for Q it needs to be the 4 betas for cubic plus start and end node I think

B <- array(NA, dim = c(500, 6, M))
for(i in 1:M) {
  y_i <- y_orig[i, f[i]:l[i]]
  x <- f[i]:l[i]
  tmp <- bs(x, knots = c(min(x, na.rm = T), floor(length(x) * 0.67), max(x, na.rm = T)))
  B[f[i]:l[i], 1:K, i] <- tmp
}

Bt <- array(NA, c(K, 500, M))
for(i in 1:M) {
  for(j in 1:500) {
    for(k in 1:K) {
      Bt[k, j, i] <- B[j, k, i] # t(B[ , , i])
    }
  }
} 

Q <- makeQ(2, K=K) # set up penalty matrix

# test to see if calcs work
# alpha0 <- matrix(NA, M, 500)
# foo <- array(NA, c(500, K, M))
# beta <- matrix(rnorm(K*M, 0, 1), K, M)
# for(i in 1:M){
#   for(j in f[i]:l[i]) {
#     for(k in 1:K) {
#       foo[j, k, i] <- beta[k, i] * Bt[k, j, i]
#     }
#     alpha0[i, j] <- sum(foo[j, , i])
#     
#   }
# }

initialize_detrend_spl = function(){
  # alpha0 = rnorm(M, rowMeans(log_y, na.rm=TRUE), 0.25)
  # alpha1 = -rlnorm(M, -1, 0.25)
  sd_y = rlnorm(M, 0, 1) 
  eta = rnorm(Tea, 0, 0.25)
  sd_eta = sd(eta) 
  beta0 = rnorm(1, 0, 0.1)
  sd_x = rlnorm(1, 0, 1) 
  return(list(sd_y = sd_y,
              eta = eta,
              sd_eta = sd_eta,
              beta0 = beta0,
              sd_x = sd_x))
}

data_detrend_spl <- list(y = log_y, 
                f = f, 
                l = l, 
                M = M, 
                Tea = Tea, 
                # a = a_use, 
                x = x_use,
                B = Bt,
                K = K,
                B_pred = Bt,
                Q = Q)

m_detrend_spl = jags.model('Code/JAGS/detrend_spline.txt', 
                  data = data_detrend_spl, 
                  inits = initialize_detrend_spl, 
                  n.chains = nc, 
                  n.adapt = nb)
out_detrend_spl = coda.samples(m_detrend_spl, c("x", "beta0", "sd_eta", "sd_x", "sd_y"), ni) # "alpha0", can monitor alpha0 if want to see the spline for each series but it takes a massive amount of memory (6 GB RAM held in R for 3 chains at 1000 iterations). can't monitor mu, y_rep, or mu_rep without massive memory

save(out_detrend_spl, file = "Results/JAGS/detrend_spl.RData")

plot(out_detrend_spl[ , c("x[1]", "sd_eta", "sd_x", "beta0")])
par(mfrow = c(1,1))

# Quick plot reconstruction
par(mfrow = c(1,1))
x_id = which(substr(varnames(out_detrend_spl),1,2)=="x[")
x_post = as.matrix(out_detrend_spl[ ,x_id])
x_post_mean <- colMeans(x_post)
plot(year, x_post_mean*x_sd + x_mean, type="l")

x_post_median <- apply(x_post, 2, median)
plot(year, x_post_median*x_sd + x_mean, type="l")

plot(x, y_i)
lines(x, postxm)

# validation
x_valid_post_m <- postxm_m2[hold_out]*x_sd + x_mean
x_valid <- x_full[hold_out]
plot(x_valid, x_valid_post_m, type = "p")
abline(0, 1, col = "red")
cor(x_valid, x_valid_post_m)

res2 <- (x_valid_post_m - x_valid)^2
sqrt(sum(res2) / length(x_valid_post_m))
sqrt(mean(res2))


# plot(outM_spline)
par(mfrow = c(1,1))

gamidx = which(substr(varnames(out_detrend_spl),1,6)=="beta.0") # finds the indices of the gamma variables
xyplot(out_detrend_spl[,gamidx]) # traceplot of the gamma variables
summary(out_detrend_spl[,gamidx]) # summary of the gamma variables

rm(out)
rm(out_detrend_spl)

#####

###### linear climate AR1 #####
initialize_ar = function(){
  alpha0 = rnorm(M, rowMeans(log_y, na.rm=TRUE), 0.25)
  alpha1 = -rlnorm(M, -1, 0.25)
  sd_y = rlnorm(M, 0, 1) 
  mu_a0 = mean(alpha0)
  mu_a1 = mean(alpha1)
  sd_a0 = sd(alpha0)
  sd_a1 = sd(alpha1)
  eta = rnorm(Tea, 0, 0.25)
  sd_eta = sd(eta) 
  beta0 = rnorm(1, 0, 0.1)
  sd_x = rlnorm(1, 0, 1) 
  delta = runif(M, -0.7, 0.7)
  return(list(sd_y = sd_y,
              mu_a0 = mu_a0,
              mu_a1 = mu_a1,
              sd_a0 = sd_a0,
              sd_a1 = sd_a1,
              sd_eta = sd_eta,
              beta0 = beta0,
              sd_x = sd_x,
              delta = delta))
}

ar_data <- list(y = log_y, 
                f = f, 
                l = l, 
                M = M, 
                Tea = Tea, 
                a = a_use, 
                x = x_use)

m_ar = jags.model('Code/JAGS/linear_ar.txt', 
                data = ar_data, 
                inits = initialize_ar, 
                n.chains = 3, 
                n.adapt = 5000)
our_ar = coda.samples(m_ar, c("x", "mu_a0", "alpha1", "sd_a0", "beta0", "sd_eta", "sd_x", "sd_y", "mean_delta"), 3000)

save(our_ar, file = "Results/JAGS/our_ar.RData")

plot(our_ar[ , c("mu_a0", "sd_a0", "sd_eta", "sd_x", "beta0", "mean_delta")])
par(mfrow = c(1,1))

summary(our_ar[ , c("mu_a0", "sd_a0", "sd_eta", "sd_x", "beta0", "mean_delta")])

xidx = which(substr(varnames(our_ar),1,2)=="x[") # finds the indices of the x variables
postxm = colMeans(as.matrix(our_ar[,xidx])) # finds the posterior mean of the x variables
plot(postxm,type="l") # plots the posterior mean of the x variables
plot(year,postxm*x_sd + x_mean,type="l") # plots the posterior mean of the x variables on the original scale (degrees C) and year on x-axis


rm(out)
rm(our_ar)

#####

##### Climate Spline #####
library(splines)

knots <- seq(0, 500, by = 25)
B <- bs(1:Tea, knots = knots)
H <- length(knots) + 3

matplot(1:Tea, B, type="l",xlab="Time",ylab="Basis function, Bj(X)",cex.lab=1.5,cex.axis=1.5,lwd=2)
# d <- max(4, floor(Tea/35))
# K <- floor(Tea/d - 1)

# visualize priors (D draws from prior distribution)
D <- 10
g <- matrix(0, Tea, D)

for(j in 1:D){
  beta  <- rnorm(length(knots) + 3, 0, 10)
  g[ , j] <- B%*%beta
}

matplot(1:500, g, lwd=2, type="l", cex.lab=1.5, cex.axis=1.5)

# run model

initialize_spl = function(){
  alpha0 = rnorm(M, rowMeans(log_y, na.rm=TRUE), 0.25)
  alpha1 = -rlnorm(M, -1, 0.25)
  sd_y = rlnorm(M, 0, 1) 
  mu_a0 = mean(alpha0)
  mu_a1 = mean(alpha1)
  sd_a0 = sd(alpha0)
  sd_a1 = sd(alpha1)
  eta = rnorm(Tea, 0, 0.25)
  sd_eta = sd(eta) 
  sd_x = rlnorm(1, 0, 1) 
  return(list(sd_y = sd_y,
              mu_a0 = mu_a0,
              mu_a1 = mu_a1,
              sd_a0 = sd_a0,
              sd_a1 = sd_a1,
              sd_eta = sd_eta,
              sd_x = sd_x))
}

m_spline_25 = jags.model('Code/JAGS/climate_spline.txt', 
                         list(y = log_y, 
                              f = f, 
                              l = l, 
                              M = M, 
                              Tea = Tea, 
                              a = a_use, 
                              x = x_use,
                              B = B,
                              H = H), 
                         inits = initialize_spl, 
                         n.chains = 3, 
                         n.adapt = 3000)

out_m_climate_spline_25 = coda.samples(m_spline_25, c("x", "mu_a0", "mu_a1", "sd_a0", "sd_a1", "sd_eta", "sd_x", "sd_y", "mu_gamma", "sd_g"), 2000)

save(out_m_climate_spline_25, file = "Results/JAGS/climate_spline_25.RData")

plot(out_m_climate_spline_25[ ,c("mu_a0", "mu_a1", "sd_a0", "sd_a1", "sd_eta", "sd_x")])
par(mfrow = c(1,1))
x_idx_25 = which(substr(varnames(out_m_climate_spline_25),1,2)=="x[") # finds the indices of the x variables
post_climate_25 = colMeans(as.matrix(out_m_climate_spline_25[ ,x_idx_25])) # finds the posterior mean of the x variables
plot(postxm,type="l") # plots the posterior mean of the x variables
plot(year, post_climate_25*x_sd + x_mean,type="l") # plots the posterior mean of the x variables on the original scale (degrees C) and year on x-axis

# validation
x_valid_post_m <- post_climate_25[hold_out]*x_sd + x_mean
x_valid <- x_full[hold_out]
plot(x_valid, x_valid_post_m, type = "p")
abline(0, 1, col = "red")
cor(x_valid, x_valid_post_m)

res2 <- (x_valid_post_m - x_valid)^2
sqrt(sum(res2) / length(x_valid_post_m))
sqrt(mean(res2))


rm(out)
rm(out_m_climate_spline_25)

##### negexp detrend 1 changepoint climate - seems good ####

x_min <- min(x_use, na.rm = T) # value of x on standardized scale when climate = 0
x_max <- max(x_use, na.rm = T)

initialize_m2_nc = function(){
  alpha0 = rnorm(M, rowMeans(log_y, na.rm=TRUE), 0.25)
  alpha1 = -rlnorm(M, -1, 0.25)
  sd_y = rlnorm(M, 0, 1) 
  mu_a0 = mean(alpha0)
  mu_a1 = mean(alpha1)
  sd_a0 = sd(alpha0)
  sd_a1 = sd(alpha1)
  eta = rnorm(Tea, 0, 0.25)
  beta0 <- rep(NA_real_, times = 2)
  beta0[1] = rnorm(2, -10, 0.5)
  beta0[2] = rnorm(2, 0.5, 0.1)
  #   sd_eta[k] = sd(eta[ , k])
  sd_eta = runif(1, 0.2, 0.3)
  x_1 = runif(1, x_min, 0)
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
              x_1 = x_1,
              sd_x = sd_x))
}

m2_nc_data <- list(y = log_y, 
                   f = f, 
                   l = l, 
                   M = M, 
                   Tea = Tea, 
                   a = a_use,
                   # K = K,
                   # species = species,
                   x_min = x_min,
                   x_max = x_max,
                   # v = 410, # or 320 
                   x = x_use)

params <- c("x", "beta0", "sd_eta", "sd_x", "sd_y", "x_1")

cl <- makeCluster(nc)                       # Request # cores
clusterExport(cl, c("m2_nc_data", "initialize_m2_nc", "params", "M", "f", "l", "Tea", "a_use", "log_y", "x_use", "x_min", "x_max", "nb", "ni", "nt"))
clusterSetRNGStream(cl = cl, 98708761)
system.time({ 
  out <- clusterEvalQ(cl, {
    library(rjags)
    jm <- jags.model("Code/JAGS/negexp_1changept.txt", m2_nc_data, initialize_m2_nc, n.adapt=nb, n.chains=1) # Compile model and run burnin
    out <- coda.samples(jm, params, n.iter=ni, thin=nt) # Sample from posterior distribution
    return(as.mcmc(out))
  })
}) # 

stopCluster(cl)

# Results
m_negexp_1change <- mcmc.list(out)

save(m_negexp_1change, file = "Results/JAGS/negexp_1change.RData")

plot(m_negexp_1change[ ,c("sd_eta", "sd_x", "beta0[1]", "beta0[2]", "x_1")])
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


rm(out)
rm(m_negexp_1change)

# recon + geom_smooth(se = FALSE) # do not do this with the MCMC chains. Maybe smooth through the mean or something else because this will kill the computer.

#####


#### RCS + climate spline 25 year knots #####
knots <- seq(0, 250, by = 25)
B <- bs(1:Tea, knots = knots) # neet to sort x_use? then doesn't line up with rest.
H <- length(knots) + 3

# run model
initialize_spl = function(){
  alpha0 = rnorm(1, rowMeans(log_y, na.rm=TRUE), 0.25)
  alpha1 = -rlnorm(1, -1, 0.25)
  sd_y = rlnorm(M, 0, 1) 
  eta = rnorm(Tea, 0, 0.25)
  sd_eta = sd(eta) 
  sd_x = rlnorm(1, 0, 1) 
  return(list(sd_y = sd_y,
              alpha0 = alpha0,
              alpha1 = alpha1,
              sd_eta = sd_eta,
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
                 H = H)

params <- c("x", "alpha0", "alpha1", "sd_eta", "sd_x", "sd_y", "mu_gamma", "sd_g")

cl <- makeCluster(nc)                       # Request # cores
clusterExport(cl, c("spl_data", "initialize_spl", "params", "M", "f", "l", "Tea", "a_use", "log_y", "x_use", "nb", "ni", "nt"))
clusterSetRNGStream(cl = cl, 8675302)
system.time({ 
  out <- clusterEvalQ(cl, {
    library(rjags)
    jm <- jags.model("Code/JAGS/rcs_climate_spl.txt", spl_data, initialize_spl, n.adapt=nb, n.chains=1) # Compile model and run burnin
    out <- coda.samples(jm, params, n.iter=ni, thin=nt) # Sample from posterior distribution
    return(as.mcmc(out))
  })
}) # 
stopCluster(cl)

# Results
m_rcs_spline_25 <- mcmc.list(out)

save(m_rcs_spline_25, file = "Results/JAGS/rcs_spline_25.RData")

plot(m_rcs_spline_25[ ,c("alpha0", "alpha1", "sd_eta", "sd_x", "mu_gamma", "sd_g")])
par(mfrow = c(1,1))

# reconstruction plots
recon <- plot_recon(m_rcs_spline_25)
ggsave(filename = "Results/Figures/torn_recon_post_rcs_spl25.tiff", plot = recon, width = 8, height = 4, units = "in") # , dpi = 1000) # poster vs paper formatting - see past work and make package or github source

y_crn$scaled <- (y_crn$HURstd * x_sd) + x_mean
recon + geom_line(data = y_crn, aes(year, scaled), color = "pink") # compare to chronology

p1 <- recon + geom_line(data = y_crn, aes(year, HURstd + 1), color = "pink") + geom_line(data = y_crn_std, aes(year, HURres), color = "gray") 

p2 <- ggplot(data = y_crn, aes(year, samp.depth)) + geom_step() + theme_bw()

# Expressed Population Size (EPS)
# library(dplR)
# rwi.stats()

require(gridExtra)
grid.arrange(p1 + theme(legend.position="top"), p2, ncol = 1, heights = c(2, 1))

recon_valid <- plot_recon(m_rcs_spline_25, valid_yrs = year[hold_out])
ggsave(filename = "Results/Figures/torn_recon_rcs_spl25_valid_back.tiff", plot = recon_valid, width = 8, height = 4, units = "in") 

rm(out)
rm(m_rcs_spline_25)

######



if(FALSE) {
  
  ##### climate p-spline - not mixing #####
  library(splines)
  
  # penalty matrix function
  makeQ = function(degree, K, epsilon=1e-3){
    x <- diag(K)
    E <- diff(x, differences=degree)
    return( t(E) %*% E + x*epsilon)
  } 
  
  # Set up basis function for each tree ring series
  K = 6 # for Q it needs to be the 4 betas for cubic plus start and end node I think
  
  knots <- seq(0, 500, by = 25)
  B <- bs(1:Tea, knots = knots)
  H <- length(knots) + 3
  
  Q <- makeQ(2, K = H) # not sure what the degree is
  # run model
  
  
  # sam clifford blog
  d <- max(4, floor(Tea/35))
  K <- floor(Tea/d - 1)
  
  bspline <- function(x, K, bdeg=3, cyclic=FALSE, xl=min(x), xr=max(x)){
    x <- as.matrix(x,ncol=1)
    
    ndx <- K - bdeg
    
    # as outlined in Eilers and Marx (1996)
    dx <- (xr - xl) / ndx
    t <- xl + dx * (-bdeg:(ndx+bdeg))
    T <- (0 * x + 1) %*% t
    X <- x %*% (0 * t + 1)
    P <- (X - T) / dx
    B <- (T <= X) & (X < (T + dx))
    r <- c(2:length(t), 1)
    
    for (k in 1:bdeg){
      B <- (P * B + (k + 1 - P) * B[ ,r]) / k; 
    }
    
    B <- B[,1:(ndx+bdeg)]
    
    if (cyclic == 1){
      for (i in 1:bdeg){
        B[ ,i] <- B[ ,i] + B[ ,K-bdeg+i]    
      }
      B <- B[ , 1:(K-bdeg)]
    }
    
    return(B)
  }
  
  B <- bspline(x = 1:500, K = H) #, xl=0, xr=1)
  
  makeQ = function(degree, K, epsilon=1e-3){
    x <- diag(K)
    E <- diff(x, differences=degree)
    return( t(E) %*% E + x*epsilon)
  } 
  
  Q <- makeQ(2, K)
  
  round(eigen(makeQ(2, K))$values, 4)
  
  
  initialize_spl = function(){
    alpha0 = rnorm(M, rowMeans(log_y, na.rm=TRUE), 0.25)
    alpha1 = -rlnorm(M, -1, 0.25)
    sd_y = rlnorm(M, 0, 1) 
    mu_a0 = mean(alpha0)
    mu_a1 = mean(alpha1)
    sd_a0 = sd(alpha0)
    sd_a1 = sd(alpha1)
    eta = rnorm(Tea, 0, 0.25)
    sd_eta = sd(eta) 
    sd_x = rlnorm(1, 0, 1) 
    return(list(sd_y = sd_y,
                mu_a0 = mu_a0,
                mu_a1 = mu_a1,
                sd_a0 = sd_a0,
                sd_a1 = sd_a1,
                sd_eta = sd_eta,
                sd_x = sd_x))
  }
  
  m_p_spline_25 = jags.model('Code/JAGS/climate_p_spline.txt', 
                             list(y = log_y, 
                                  f = f, 
                                  l = l, 
                                  M = M, 
                                  Tea = Tea, 
                                  a = a_use, 
                                  x = x_use,
                                  Q = Q,
                                  B = B,
                                  K = K), 
                             inits = initialize_spl, 
                             n.chains = 3, 
                             n.adapt = 1000)
  
  out_p_spline_25 = coda.samples(m_p_spline_25, c("x", "gamma_0", "gamma_00", "lambda_x", "sd_eta", "sd_x", "sd_y"), 1000) # "alpha0", can monitor alpha0 if want to see the spline for each series but it takes a massive amount of memory (6 GB RAM held in R for 3 chains at 1000 iterations). can't monitor mu, y_rep, or mu_rep without massive memory
  
  plot(out_p_spline_25[ , c("x[100]", "sd_eta", "sd_x", "gamma_00", "lambda_x")])
  par(mfrow = c(1,1))
  
  # Quick plot reconstruction
  par(mfrow = c(1,1))
  x_id = which(substr(varnames(out_p_spline_25),1,2)=="x[")
  x_post = as.matrix(out_p_spline_25[ ,x_id])
  x_post_mean <- colMeans(x_post)
  plot(year, x_post_mean*x_sd + x_mean, type="l")
  
  x_post_median <- apply(x_post, 2, median)
  plot(year, x_post_median*x_sd + x_mean, type="l")
  
  plot(x, y_i)
  lines(x, postxm)
  
  # validation
  x_valid_post_m <- postxm_m2[hold_out]*x_sd + x_mean
  x_valid <- x_full[hold_out]
  plot(x_valid, x_valid_post_m, type = "p")
  abline(0, 1, col = "red")
  cor(x_valid, x_valid_post_m)
  
  res2 <- (x_valid_post_m - x_valid)^2
  sqrt(sum(res2) / length(x_valid_post_m))
  sqrt(mean(res2))
  
  
  
  
  #####
  
#### climate spline 100 year knots #####

knots <- seq(0, 500, by = 100)
B <- bs(1:Tea, knots = knots) # neet to sort x_use? then doesn't line up with rest.
H <- length(knots) + 3

matplot(1:Tea, B, type="l",xlab="Time",ylab="Basis function, Bj(X)",cex.lab=1.5,cex.axis=1.5,lwd=2)
# d <- max(4, floor(Tea/35))
# K <- floor(Tea/d - 1)

# visualize priors (D draws from prior distribution)
D <- 10
g <- matrix(0, Tea, D)

for(j in 1:D){
  beta  <- rnorm(length(knots) + 3, 0, 10)
  g[ , j] <- B%*%beta
}

matplot(1:500, g, lwd=2, type="l", cex.lab=1.5, cex.axis=1.5)

# run model

initialize_spl = function(){
  alpha0 = rnorm(M, rowMeans(log_y, na.rm=TRUE), 0.25)
  alpha1 = -rlnorm(M, -1, 0.25)
  sd_y = rlnorm(M, 0, 1) 
  mu_a0 = mean(alpha0)
  mu_a1 = mean(alpha1)
  sd_a0 = sd(alpha0)
  sd_a1 = sd(alpha1)
  eta = rnorm(Tea, 0, 0.25)
  sd_eta = sd(eta) 
  sd_x = rlnorm(1, 0, 1) 
  return(list(sd_y = sd_y,
              mu_a0 = mu_a0,
              mu_a1 = mu_a1,
              sd_a0 = sd_a0,
              sd_a1 = sd_a1,
              sd_eta = sd_eta,
              sd_x = sd_x))
}

m_spline_100 = jags.model('Code/JAGS/climate_spline.txt', 
                          list(y = log_y, 
                               f = f, 
                               l = l, 
                               M = M, 
                               Tea = Tea, 
                               a = a_use, 
                               x = x_use,
                               B = B,
                               H = H), 
                          inits = initialize_spl, 
                          n.chains = 3, 
                          n.adapt = 3000)

out_m_climate_spline_100 = coda.samples(m_spline_100, c("x", "mu_a0", "mu_a1", "sd_a0", "sd_a1", "sd_eta", "sd_x", "sd_y", "mu_gamma", "sd_g"), 2000)

plot(out_m_climate_spline_100[ ,c("mu_a0", "mu_a1", "sd_a0", "sd_a1", "sd_eta", "sd_x")])
par(mfrow = c(1,1))
x_id_100 = which(substr(varnames(out_m_climate_spline_100),1,2)=="x[") # finds the indices of the x variables
post_climate_100 = colMeans(as.matrix(out_m_climate_spline_100[,x_id_100])) # finds the posterior mean of the x variables
plot(post_climate_100,type="l") # plots the posterior mean of the x variables
plot(year, post_climate_100*x_sd + x_mean,type="l") # plots the posterior mean of the x variables on the original scale (degrees C) and year on x-axis

# validation
x_valid_post_m <- post_climate_100[hold_out]*x_sd + x_mean
x_valid <- x_full[hold_out]
plot(x_valid, x_valid_post_m, type = "p")
abline(0, 1, col = "red")
cor(x_valid, x_valid_post_m)

res2 <- (x_valid_post_m - x_valid)^2
sqrt(sum(res2) / length(x_valid_post_m))
sqrt(mean(res2))


#####

##### negexp detrend 2 changepoints climate  - no convergence ####
# probably not enough observed values when trees stop growing, would need strong prior or to fix it as Schofield did based on other studies

x_min <- min(x_use, na.rm = T)
x_max <- max(x_use, na.rm = T)

initialize_m2_nc = function(){
  alpha0 = rnorm(M, rowMeans(log_y, na.rm=TRUE), 0.25)
  alpha1 = -rlnorm(M, -1, 0.25)
  sd_y = rlnorm(M, 0, 1) 
  mu_a0 = mean(alpha0)
  mu_a1 = mean(alpha1)
  sd_a0 = sd(alpha0)
  sd_a1 = sd(alpha1)
  eta = rnorm(Tea, 0, 0.25)
  beta0 <- rep(NA_real_, times = 2)
  # beta0[1] = rnorm(2, -10, 0.5)
  beta0[2] = rnorm(2, 0.5, 0.1)
  beta0[3] = rnorm(2, 0.5, 0.1)
  #   sd_eta[k] = sd(eta[ , k])
  sd_eta = runif(1, 0.2, 0.3)
  x_1 = runif(2, x_min, x_max)
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
              x_1 = x_1,
              sd_x = sd_x))
}

m2_nc_data <- list(y = log_y, 
                   f = f, 
                   l = l, 
                   M = M, 
                   Tea = Tea, 
                   a = a_use,
                   # K = K,
                   # species = species,
                   x_min = x_min,
                   x_max = x_max,
                   # v = 410, # or 320 
                   x = x_use)

params <- c("x", "beta0", "sd_eta", "sd_x", "sd_y", "x_1", "x_change")

cl <- makeCluster(nc)                       # Request # cores
clusterExport(cl, c("m2_nc_data", "initialize_m2_nc", "params", "M", "f", "l", "Tea", "a_use", "log_y", "x_use", "x_min", "x_max", "nb", "ni", "nt"))
clusterSetRNGStream(cl = cl, 98708761)
system.time({ 
  out <- clusterEvalQ(cl, {
    library(rjags)
    jm <- jags.model("Code/JAGS/negexp_2changept.txt", m2_nc_data, initialize_m2_nc, n.adapt=nb, n.chains=1) # Compile model and run burnin
    out <- coda.samples(jm, params, n.iter=ni, thin=nt) # Sample from posterior distribution
    return(as.mcmc(out))
  })
}) # 

stopCluster(cl)

# Results
m_negexp_2change <- mcmc.list(out)

plot(m_negexp_2change[ ,c("sd_eta", "sd_x", "beta0[1]", "beta0[2]", "x_1[1]")])
par(mfrow = c(1,1))
x_id_50 = which(substr(varnames(m_negexp_2change),1,2)=="x[") # finds the indices of the x variables
post_climate_50 = colMeans(as.matrix(m_negexp_2change[,x_id_50])) # finds the posterior mean of the x variables
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
    g <- g + geom_fan() + geom_line(data = climate_valid, aes(year, temp), colour = "red") + ylab("Annual Precipitation (cm)") + xlab("Year") + scale_fill_distiller() + theme_bw()
  } else {
    # scaled posterior interval
    g <- ggplot(temp_df_long, aes(x = year, y = temp))
    g <- g + geom_fan() + geom_line(data = obs, aes(x = year, y = value), colour="black", size = 0.2) + ylab("Annual Precipitation (cm)") + xlab("Year") + theme_bw() + scale_fill_distiller() #x_full),
  }
  return(g)
}

recon <- plot_recon(m_negexp_norm, obs = data.frame(year = years, value = x_full))

# consider doing observed values as points to better see where they fall within the credible interval
recon_valid <- plot_recon(m_negexp_norm, obs = data.frame(year = years, value = x_full), valid_yrs = years[hold_out])

## any better than random points around the mean?


# recon + geom_smooth(se = FALSE) # do not do this with the MCMC chains. Maybe smooth through the mean or something else because this will kill the computer.

#####

#### Run splines or neg exp on series and if not negative then remove from data before passing to model?? ####


#################
# models of TS using year instead of age with partial pooling and not
# - Tornetrask data don't have actual age but are assumed to get near pith and that first ring if first year of growth
################

#### climate spline 50 year knots #####

knots <- seq(0, 500, by = 50)
B <- bs(1:Tea, knots = knots) # neet to sort x_use? then doesn't line up with rest.
H <- length(knots) + 3

matplot(1:Tea, B, type="l",xlab="Time",ylab="Basis function, Bj(X)",cex.lab=1.5,cex.axis=1.5,lwd=2)
# d <- max(4, floor(Tea/35))
# K <- floor(Tea/d - 1)

# visualize priors (D draws from prior distribution)
D <- 10
g <- matrix(0, Tea, D)

for(j in 1:D){
  beta  <- rnorm(length(knots) + 3, 0, 10)
  g[ , j] <- B%*%beta
}

matplot(1:500, g, lwd=2, type="l", cex.lab=1.5, cex.axis=1.5)

# run model
initialize_spl = function(){
  alpha0 = rnorm(M, rowMeans(log_y, na.rm=TRUE), 0.25)
  alpha1 = -rlnorm(M, -1, 0.25)
  sd_y = rlnorm(M, 0, 1) 
  mu_a0 = mean(alpha0)
  mu_a1 = mean(alpha1)
  sd_a0 = sd(alpha0)
  sd_a1 = sd(alpha1)
  eta = rnorm(Tea, 0, 0.25)
  sd_eta = sd(eta) 
  sd_x = rlnorm(1, 0, 1) 
  return(list(sd_y = sd_y,
              mu_a0 = mu_a0,
              mu_a1 = mu_a1,
              sd_a0 = sd_a0,
              sd_a1 = sd_a1,
              sd_eta = sd_eta,
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
                 H = H)

params <- c("x", "mu_a0", "mu_a1", "sd_a0", "sd_a1", "sd_eta", "sd_x", "sd_y", "mu_gamma", "sd_g")

cl <- makeCluster(nc)                       # Request # cores
clusterExport(cl, c("spl_data", "initialize_spl", "params", "M", "f", "l", "Tea", "a_use", "log_y", "x_use", "nb", "ni", "nt"))
clusterSetRNGStream(cl = cl, 8675302)
system.time({ 
  out <- clusterEvalQ(cl, {
    library(rjags)
    jm <- jags.model("Code/JAGS/climate_spline.txt", spl_data, initialize_spl, n.adapt=nb, n.chains=1) # Compile model and run burnin
    out <- coda.samples(jm, params, n.iter=ni, thin=nt) # Sample from posterior distribution
    return(as.mcmc(out))
  })
}) # 
stopCluster(cl)

# Results
m_spline_50 <- mcmc.list(out)
# 
# m_spline_50 = jags.model('Code/JAGS/climate_spline.txt', 
#                          list(y = log_y, 
#                               f = f, 
#                               l = l, 
#                               M = M, 
#                               Tea = Tea, 
#                               a = a_use, 
#                               x = x_use,
#                               B = B,
#                               H = H), 
#                          inits = initialize_spl, 
#                          n.chains = 3, 
#                          n.adapt = 3000)
# 
# out_m_climate_spline_50 = coda.samples(m_spline_50, c("x", "mu_a0", "mu_a1", "sd_a0", "sd_a1", "sd_eta", "sd_x", "sd_y", "mu_gamma", "sd_g"), 2000)

plot(m_spline_50[ ,c("mu_a0", "mu_a1", "sd_a0", "sd_a1", "sd_eta", "sd_x", "mu_gamma", "sd_g")])
par(mfrow = c(1,1))
x_id_50 = which(substr(varnames(out_m_climate_spline_50),1,2)=="x[") # finds the indices of the x variables
post_climate_50 = colMeans(as.matrix(out_m_climate_spline_50[,x_id_50])) # finds the posterior mean of the x variables
plot(post_climate_50,type="l") # plots the posterior mean of the x variables
plot(year, post_climate_50*x_sd + x_mean,type="l") # plots the posterior mean of the x variables on the original scale (degrees C) and year on x-axis

# validation
x_valid_post_m <- post_climate_50[hold_out]*x_sd + x_mean
x_valid <- x_full[hold_out]
plot(x_valid, x_valid_post_m, type = "p")
abline(0, 1, col = "red")
cor(x_valid, x_valid_post_m)

res2 <- (x_valid_post_m - x_valid)^2
sqrt(sum(res2) / length(x_valid_post_m))
sqrt(mean(res2))

# reconstruction plots
recon <- plot_recon(m_spline_50)
ggsave(filename = "Results/Figures/torn_recon_post_spl50.tiff", plot = recon, width = 8, height = 4, units = "in") # , dpi = 1000) # poster vs paper formatting - see past work and make package or github source

recon_valid <- plot_recon(m_spline_50, valid_yrs = year[hold_out])
ggsave(filename = "Results/Figures/torn_recon_spl50_valid_back.tiff", plot = recon_valid, width = 8, height = 4, units = "in") 

#### climate spline 100 year knots #####

knots <- seq(0, 500, by = 100)
B <- bs(1:Tea, knots = knots) # neet to sort x_use? then doesn't line up with rest.
H <- length(knots) + 3

matplot(1:Tea, B, type="l",xlab="Time",ylab="Basis function, Bj(X)",cex.lab=1.5,cex.axis=1.5,lwd=2)
# d <- max(4, floor(Tea/35))
# K <- floor(Tea/d - 1)

# visualize priors (D draws from prior distribution)
D <- 10
g <- matrix(0, Tea, D)

for(j in 1:D){
  beta  <- rnorm(length(knots) + 3, 0, 10)
  g[ , j] <- B%*%beta
}

matplot(1:500, g, lwd=2, type="l", cex.lab=1.5, cex.axis=1.5)

# run model

initialize_spl = function(){
  alpha0 = rnorm(M, rowMeans(log_y, na.rm=TRUE), 0.25)
  alpha1 = -rlnorm(M, -1, 0.25)
  sd_y = rlnorm(M, 0, 1) 
  mu_a0 = mean(alpha0)
  mu_a1 = mean(alpha1)
  sd_a0 = sd(alpha0)
  sd_a1 = sd(alpha1)
  eta = rnorm(Tea, 0, 0.25)
  sd_eta = sd(eta) 
  sd_x = rlnorm(1, 0, 1) 
  return(list(sd_y = sd_y,
              mu_a0 = mu_a0,
              mu_a1 = mu_a1,
              sd_a0 = sd_a0,
              sd_a1 = sd_a1,
              sd_eta = sd_eta,
              sd_x = sd_x))
}

m_spline_100 = jags.model('Code/JAGS/climate_spline.txt', 
                          list(y = log_y, 
                               f = f, 
                               l = l, 
                               M = M, 
                               Tea = Tea, 
                               a = a_use, 
                               x = x_use,
                               B = B,
                               H = H), 
                          inits = initialize_spl, 
                          n.chains = 3, 
                          n.adapt = 3000)

out_m_climate_spline_100 = coda.samples(m_spline_100, c("x", "mu_a0", "mu_a1", "sd_a0", "sd_a1", "sd_eta", "sd_x", "sd_y", "mu_gamma", "sd_g"), 2000)

plot(out_m_climate_spline_100[ ,c("mu_a0", "mu_a1", "sd_a0", "sd_a1", "sd_eta", "sd_x")])
par(mfrow = c(1,1))
x_id_100 = which(substr(varnames(out_m_climate_spline_100),1,2)=="x[") # finds the indices of the x variables
post_climate_100 = colMeans(as.matrix(out_m_climate_spline_100[,x_id_100])) # finds the posterior mean of the x variables
plot(post_climate_100,type="l") # plots the posterior mean of the x variables
plot(year, post_climate_100*x_sd + x_mean,type="l") # plots the posterior mean of the x variables on the original scale (degrees C) and year on x-axis

# validation
x_valid_post_m <- post_climate_100[hold_out]*x_sd + x_mean
x_valid <- x_full[hold_out]
plot(x_valid, x_valid_post_m, type = "p")
abline(0, 1, col = "red")
cor(x_valid, x_valid_post_m)

res2 <- (x_valid_post_m - x_valid)^2
sqrt(sum(res2) / length(x_valid_post_m))
sqrt(mean(res2))


#####

##### Tree-specific climate response - not identifiable?  ####
# variation in climate response by tree - not identifiable - could try informative priors

initialize_tree_clim = function(){
  alpha0 = rnorm(M, rowMeans(log_y, na.rm=TRUE), 0.25)
  alpha1 = -rlnorm(M, -1, 0.25)
  sd_y = rlnorm(M, 0, 1) 
  mu_a0 = mean(alpha0)
  mu_a1 = mean(alpha1)
  sd_a0 = sd(alpha0)
  sd_a1 = sd(alpha1)
  eta = rnorm(Tea, 0, 0.25)
  sd_eta = sd(eta) 
  beta0 = rnorm(M, 0, 0.1)
  mu_b0 = mean(beta0)
  sd_b0 = sd(beta0)
  sd_x = rlnorm(1, 0, 1) 
  return(list(sd_y = sd_y,
              mu_a0 = mu_a0,
              mu_a1 = mu_a1,
              sd_a0 = sd_a0,
              sd_a1 = sd_a1,
              sd_eta = sd_eta,
              mu_b0 = mu_b0,
              sd_b0 = sd_b0,
              sd_x = sd_x))
}

m_tree_clim = jags.model('Code/JAGS/tree_clim_nc.txt', 
                   list(y = log_y, 
                        f = f, 
                        l = l, 
                        M = M, 
                        Tea = Tea, 
                        a = a_use, 
                        x = x_use), 
                   inits = initialize_tree_clim, 
                   n.chains = 3, 
                   n.adapt = 1000)

out_m_tree_clim = coda.samples(m_tree_clim, c("x", "mu_a0", "mu_a1", "sd_a0", "sd_a1", "mu_b0", "sd_b0", "sd_eta", "sd_x", "sd_y"), 1000)

plot(out_m_tree_clim[ ,c("mu_a0", "mu_a1", "sd_a0", "sd_a1", "mu_b0", "sd_b0", "sd_eta", "sd_x")])
par(mfrow = c(1,1))
xidx = which(substr(varnames(out_m_tree_clim),1,2)=="x[") # finds the indices of the x variables
postxm = colMeans(as.matrix(out_m_tree_clim[,xidx])) # finds the posterior mean of the x variables
plot(postxm,type="l") # plots the posterior mean of the x variables
plot(year,postxm*x_sd + x_mean,type="l") # plots the posterior mean of the x variables on the original scale (degrees C) and year on x-axis





########## Climate and Detrending Splines - not working right - need tensor splines for decorrelating GAM? #########

# Set up detrending splines for each series independently
# penalty matrix function
makeQ = function(degree, K, epsilon=1e-3){
  x <- diag(K)
  E <- diff(x, differences=degree)
  return( t(E) %*% E + x*epsilon)
} 

# Set up basis function for each tree ring series
K = 6 # for Q it needs to be the 4 betas for cubic plus start and end node I think

B <- array(NA, dim = c(500, 6, M))
for(i in 1:M) {
  y_i <- y_orig[i, f[i]:l[i]]
  x <- f[i]:l[i]
  tmp <- bs(x, knots = c(min(x, na.rm = T), floor(length(x) * 0.67), max(x, na.rm = T)))
  B[f[i]:l[i], 1:K, i] <- tmp
}

Bt <- array(NA, c(K, 500, M))
for(i in 1:M) {
  for(j in 1:500) {
    for(k in 1:K) {
      Bt[k, j, i] <- B[j, k, i] # t(B[ , , i])
    }
  }
} 

Q <- makeQ(2, K=K) # set up penalty matrix


# Climate Basis Function Setup
knots <- seq(0, 500, by = 50)
B_c <- bs(1:Tea, knots = knots) # neet to sort x_use? then doesn't line up with rest.
H <- length(knots) + 3

# run model
initialize_spl_spl = function(){
  sd_y = rlnorm(M, 0, 1) 
  eta = rnorm(Tea, 0, 0.25)
  sd_eta = sd(eta) 
  sd_x = rlnorm(1, 0, 1) 
  return(list(sd_y = sd_y,
              sd_eta = sd_eta,
              sd_x = sd_x))
}

m_spline_spline_50 = jags.model('Code/JAGS/climate_detrend_splines.txt', 
                         list(y = log_y, 
                              f = f, 
                              l = l, 
                              M = M, 
                              Tea = Tea, 
                              # a = a_use, 
                              x = x_use,
                              Q = Q,
                              B = Bt,
                              K = K,
                              B_c = B_c,
                              H = H), 
                         inits = initialize_spl_spl, 
                         n.chains = 3, 
                         n.adapt = 3000)

out_m_spline_spline_50 = coda.samples(m_spline_spline_50, c("x", "sd_eta", "sd_x", "sd_y", "mu_gamma", "sd_g"), 2000)

plot(out_m_spline_spline_50[ ,c("sd_y[1]", "x[1]", "sd_eta", "sd_x")])
par(mfrow = c(1,1))
x_id_spl_50 = which(substr(varnames(out_m_spline_spline_50),1,2)=="x[") # finds the indices of the x variables
post_spline_spline_50 = colMeans(as.matrix(out_m_spline_spline_50[ ,x_id_spl_50])) # finds the posterior mean of the x variables
plot(post_spline_spline_50,type="l") # plots the posterior mean of the x variables
plot(year, post_spline_spline_50*x_sd + x_mean,type="l") # plots the posterior mean of the x variables on the original scale (degrees C) and year on x-axis












params <- c("x", "mu_a0", "mu_a1", "sd_a0", "sd_a1", "beta0", "sd_eta", "sd_x", "sd_y")

cl <- makeCluster(3)                       # Request # cores
clusterExport(cl, c("m2_data", "initizalize_spl", "params")) # , "Nst", "ni", "nb", "nt")) # Make these available
clusterSetRNGStream(cl = cl, 13459784)

system.time({ # no status bar (% complete) when run in parallel
  out <- clusterEvalQ(cl, {
    library(rjags)
    jm <- jags.model("Code/Jags_Models/final_elev_od.txt", pjor.od.data, inits, n.adapt=nb, n.chains=1) # Compile model and run burnin
    out <- coda.samples(jm, params, n.iter = 1000, thin = 1) # Sample from posterior distribution
    return(as.mcmc(out))
  })
}) # 

stopCluster(cl)

# Results
pjor_od <- mcmc.list(out)


#####

##### M2 hierarchical centering - currently not identifiable #####
initialize_m2c = function(){
  alpha0 = rnorm(M, rowMeans(log_y, na.rm=TRUE), 0.25)
  alpha1 = -rlnorm(M, -1, 0.25)
  sd_y = rlnorm(M, 0, 1) 
  mu_a0 = mean(alpha0)
  mu_a1 = mean(alpha1)
  sd_a0 = sd(alpha0)
  sd_a1 = sd(alpha1)
  eta = rnorm(Tea, 0, 0.25)
  sd_eta = sd(eta) 
  beta0 = rnorm(1, 0, 0.1)
  sd_x = rlnorm(1, 0, 1) 
  return(list(sd_y = sd_y,
              mu_a0 = mu_a0,
              mu_a1 = mu_a1,
              sd_a0 = sd_a0,
              sd_a1 = sd_a1,
              sd_eta = sd_eta,
              beta0 = beta0,
              sd_x = sd_x))
}

m2_data <- list(y = log_y, 
                f = f, 
                l = l, 
                M = M, 
                Tea = Tea, 
                a = a_use, 
                x = x_use)

M2 = jags.model('Code/JAGS/modM2_rm.txt', 
                list(y = log_y, 
                     f = f, 
                     l = l, 
                     M = M, 
                     Tea = Tea, 
                     a = a_use, 
                     x = x_use), 
                inits = initialize_m2c, 
                n.chains = 3, 
                n.adapt = 1000)
outM2 = coda.samples(M2, c("x", "mu_a0", "mu_a1", "sd_a0", "sd_a1", "beta0", "sd_eta", "sd_x", "sd_y"), 1000)

outM2_list <- mcmc.list(outM2)
plot(outM2[ ,c("mu_a0", "mu_a1", "sd_a0", "sd_a1", "beta0", "sd_eta", "sd_x")])



#####

##### Models from other papers #####

##### Model M2: stable climate #####
initialize_m2 = function(){ # function that generates initial values for the model
  inita0 = rnorm(M,rowMeans(log_y,na.rm=TRUE),0.25) # initial values for alpha_0 -- centered around the mean yuse value
  inita1 = -rlnorm(M,-1,0.25) # initial values for alpha_1 -- note that these must be negative
  sdy = rlnorm(M) # initial values for sigma
  mua = c(mean(inita0),mean(inita1)) # use the initial values for alpha_0 and alpha_1 to find initial value for mu_a0 and mu_a1
  sda = c(sd(inita0),sd(inita1)) # as above to find initial value for sig_a0 and sig_a1
  eta = rnorm(Tea,0,0.25) # find initial values for eta
  sdeta = sd(eta) # find initial value for sig_eta
  gamma = rnorm(2,c(0,0.4),0.1) # initial values for gamma.  The initial value for gamma 2 is set to be positive
  sdx = rlnorm(1) # initial value for sig_x
  return(list(alpha = cbind(inita0,inita1),eta = eta, mua = mua, sda = sda, sdeta = sdeta ,sdy = sdy,gamma=gamma,sdx = sdx))
}

M2 = jags.model('Code/JAGS/modM2.txt', list(y = log_y, f = f, l = l, M = M, Tee = Tea, a = a_use, x=x_use), inits = initialize_m2, n.chains = 3, n.adapt = 500)
outM2 = coda.samples(M2,c('x','eta','alpha','sdy','sdeta','mua','sda','gamma','sdx'),1000)

gamidx = which(substr(varnames(outM2),1,5)=="gamma") # finds the indices of the gamma variables
xyplot(outM2[,gamidx]) # traceplot of the gamma variables
summary(outM2[,gamidx]) # summary of the gamma variables

xidx = which(substr(varnames(outM2),1,2)=="x[") # finds the indices of the x variables
postxm = colMeans(as.matrix(outM2[,xidx])) # finds the posterior mean of the x variables
plot(postxm,type="l") # plots the posterior mean of the x variables
plot(year,postxm*x_sd + x_mean,type="l") # plots the posterior mean of the x variables on the original scale (degrees C) and year on x-axis

##### RUNNING MODEL M6 #####
initfunM6 = function(){ # function that generates initial values for the model
  inita0 = rnorm(M,rowMeans(yuse,na.rm=TRUE),0.25) # initial values for alpha_0 -- centered around the mean yuse value
  inita1 = -rlnorm(M,-1,0.25) # initial values for alpha_1 -- note that these must be negative
  sdy = rlnorm(M) # initial values for sigma
  mua = c(mean(inita0),mean(inita1)) # use the initial values for alpha_0 and alpha_1 to find initial value for mu_a0 and mu_a1
  sda = c(sd(inita0),sd(inita1)) # as above to find initial value for sig_a0 and sig_a1
  eta = rnorm(Y,0,0.25) # find initial values for eta
  sdeta = sd(eta) # find initial value for sig_eta
  gamma = rnorm(2,c(0,0.4),0.1) # initial values for gamma.  The initial value for gamma 2 is set to be positive
  sdx = rlnorm(1) # initial value for sig_x
  rho = 2*rbeta(1,2,1)-1 # initial value for AR1 parameter -- "triangle"
  return(list(alpha = cbind(inita0,inita1),eta = eta, mua = mua, sda = sda, sdeta = sdeta ,sdy = sdy,gamma=gamma,sdx = sdx,rho = rho))
}

M6 = jags.model('modM6.txt', list(y = yuse, f = f, l = l, M = M, Tee = Tee, a = ause, x=xuse), inits = initfunM6, n.chains = 3, n.adapt = 1000)
outM6 = coda.samples(M6,c('x','eta','alpha','sdy','sdeta','mua','sda','gamma','sdx','rho'),1000)

rhoidx = which(substr(varnames(outM6),1,3)=="rho") # finds the indices of the gamma variables
xyplot(outM6[,rhoidx]) # traceplot of the gamma variables

##### M7 hugershoff #####
M7 = jags.model('hugershoff_autocorr.txt', list(y = yuse, f = f, l = l, M = M, Tee = Tee, a = ause, x=xuse), inits = initfunM6, n.chains = 3, n.adapt = 1000)
outM7 = coda.samples(M6,c('x','eta','alpha','sdy','sdeta','mua','sda','gamma','sdx','rho'),1000)

#####

#### Extra #####

M2_con = jags.model('Examples/Schofield_2016/MB_TS_CON.txt', 
                    list(y = log_y, 
                         f = f, 
                         l = l, 
                         m = 410, 
                         k = M,
                         n = Tea, 
                         a = a_use, 
                         x = x_use), 
                    # inits = initialize_m2_nc, 
                    n.chains = 3, 
                    n.adapt = 500)

outM2_con = coda.samples(M2_con, c("x", "gamma", "sdeta", "sdx", "sdy"), 1000)

plot(outM2_con[ ,c("mu_a0", "mu_a1", "sd_a0", "sd_a1", "beta0", "sd_eta", "sd_x")])

# compare Schofield and my reconstructions
xidx_m2 = which(substr(varnames(outM2),1,2)=="x[")
postxm_m2 = colMeans(as.matrix(outM2[,xidx_m2])) 
xidx = which(substr(varnames(outM2_nc),1,2)=="x[")
postxm = colMeans(as.matrix(outM2_nc[,xidx])) 
xidx_con = which(substr(varnames(outM2_con),1,2)=="x[")
postxm_con = colMeans(as.matrix(outM2_con[,xidx_con])) 

par(mfrow = c(3, 1))
plot(year,postxm_m2*x_sd + x_mean,type="l")
plot(year,postxm_con*x_sd + x_mean,type="l")
plot(year,postxm*x_sd + x_mean,type="l")
par(mfrow = c(1, 1))

plot(postxm_m2, postxm)
abline(0, 1, col = "red")

plot(postxm_con, postxm)
abline(0, 1, col = "red")

plot(x_full[hold_out], postxm[hold_out]*x_sd + x_mean)
abline(0, 1, col = "red")

lm1 <- lm(as.numeric(x_full[hold_out]) ~ as.numeric(postxm[hold_out]*x_sd + x_mean))
summary(lm1)

}