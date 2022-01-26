#####################################################################
# Hierarchical reconstruction of precip from Mount Washinton bristlecone pine tree rings
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

# Move to source or package
# reconstruction plot
plot_recon <- function(mcmc, obs = climate, mean = x_mean, sd = x_sd, valid_yrs = NULL) {
  xidx_m2 = which(substr(varnames(mcmc),1,2)=="x[")
  temp_df <- as_tibble(t(as.matrix(mcmc[ , xidx_m2])))
  temp_df <- temp_df %>% 
    mutate(year = year)
  
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
    g <- g + geom_fan() + geom_line(data = climate_valid, aes(year, temp), colour = "red") + ylab("Temperature (C)") + xlab("Year") + scale_fill_distiller() + theme_bw()
  } else {
    # scaled posterior interval
    g <- ggplot(temp_df_long, aes(x = year, y = temp))
    g <- g + geom_fan() + geom_line(data = obs, aes(x = year, y = value), colour="black", size = 0.2) + ylab("Temperature (C)") + xlab("Year") + theme_bw() + scale_fill_distiller() #x_full),
  }
  return(g)
}

#### Load and Prep Data #####
load(file = "Data/itrdb_pilo_mount_washington/pilo_rwl_climate.RData")

# Filter to just summer precip
tree_climate <- tree_climate %>%
  dplyr::filter(type == "ppt" | is.na(type)) # %>%
 # dplyr::filter(month %in% c(NA, 5, 6, 7, 8))

# make year vector
year <- unique(tree_climate$year)

# standardize climate data
climate <- tree_climate %>%
  group_by(year) %>%
  summarise(value = mean(value)) %>%
  ungroup() # %>%
  # dplyr::select(ppt)

x_full <- as.numeric(climate$value)

x_mean <- mean(x_full, na.rm=TRUE)
x_sd <- sd(x_full, na.rm=TRUE);
x_s <- (x_full - x_mean) / x_sd # standardize temperatures to z-scores for numerical stability
x_use <- x_s

# Hold out years for validation (first half)
climate_rows <- which(!is.na(x_full))  # the indices that we hold out
hold_out <- climate_rows[1:floor(length(climate_rows) / 2)]
x_use[hold_out] <- NA

# Get tree ring data into tree x year matrix
y1 <- tree_climate %>%
  dplyr::select(tree, year, rwi) %>%
  group_by(tree, year) %>%
  summarise(rwi = mean(rwi)) # mean across cores from same tree because of lack of time to deal with in a more nuanced way
  
y_orig <- y1 %>%
  ungroup() %>%
  tidyr::spread(key = year, value = rwi) %>%
  dplyr::select(-tree) %>%
  as.data.frame() %>%
  as.matrix()

y_orig <- y_orig[1:(nrow(y_orig)-1), ] # last row is a problem

y <- y_orig + 0.00001 # add tiny increase for years with zero growth (missing rings?)

log_y <- log(y)

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

########## Run Models ###########
testing <- TRUE
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

if(!dir.exists("Results/JAGS")) dir.create("Results/JAGS", recursive = TRUE)

##### negexp detrend linear climate M2-non-centered - not working ####

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

m2_nc_data <- list(y = y_orig, 
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
m_negexp_norm <- mcmc.list(out)


plot(m_negexp_norm[ ,c("mu_a0", "mu_a1", "sd_a0", "sd_a1", "sd_eta", "sd_x", "beta0")])
par(mfrow = c(1,1))
x_id_50 = which(substr(varnames(m_negexp_norm),1,2)=="x[") # finds the indices of the x variables
post_climate_50 = colMeans(as.matrix(m_negexp_norm[,x_id_50])) # finds the posterior mean of the x variables
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
recon <- plot_recon(m_negexp_norm, obs = tree_climate)

ggsave(filename = "Results/Figures/wash_recon_post_negexp_const.tiff", plot = recon, width = 8, height = 4, units = "in") # , dpi = 1000) # poster vs paper formatting - see past work and make package or github source

recon_valid <- plot_recon(m_negexp_norm, valid_yrs = year[hold_out])
ggsave(filename = "Results/Figures/wash_recon_negexp_const_valid_back.tiff", plot = recon_valid, width = 8, height = 4, units = "in") 




#####

##### spline on age standardization - not working  #####

library(splines)

# penalty matrix function
makeQ = function(degree, K, epsilon=1e-3){
  x <- diag(K)
  E <- diff(x, differences=degree)
  return( t(E) %*% E + x*epsilon)
} 

# Set up basis function for each tree ring series
K = 6 # for Q it needs to be the 4 betas for cubic plus start and end node I think

B <- array(NA, dim = c(Tea, 6, M))
for(i in 1:M) {
  y_i <- y_orig[i, f[i]:l[i]]
  x <- f[i]:l[i]
  tmp <- bs(x, knots = c(min(x, na.rm = T), floor(length(x) * 0.67), max(x, na.rm = T)))
  B[f[i]:l[i], 1:K, i] <- tmp
}

Bt <- array(NA, c(K, Tea, M))
for(i in 1:M) {
  for(j in 1:Tea) {
    for(k in 1:K) {
      Bt[k, j, i] <- B[j, k, i] # t(B[ , , i])
    }
  }
} 

Q <- makeQ(2, K=K) # set up penalty matrix

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

params <- c("x", "beta0", "sd_eta", "sd_x", "sd_y")

library(parallel)
cl <- makeCluster(nc)                       # Request # cores
clusterExport(cl, c("data_detrend_spl", "initialize_detrend_spl", "params", "M", "f", "l", "Tea", "a_use", "log_y", "x_use", "nb", "ni", "nt", "Bt", "Q", "K"))
clusterSetRNGStream(cl = cl, 8675303)
system.time({ 
  out <- clusterEvalQ(cl, {
    library(rjags)
    jm <- jags.model("Code/JAGS/detrend_spline.txt", data_detrend_spl, initialize_detrend_spl, n.adapt=nb, n.chains=1) # Compile model and run burnin
    out <- coda.samples(jm, params, n.iter=ni, thin=nt) # Sample from posterior distribution
    return(as.mcmc(out))
  })
}) # 

stopCluster(cl)

# Results
m_spline_norm <- mcmc.list(out)

plot(m_spline_norm[ , c("x[1]", "sd_eta", "sd_x", "beta0")])
par(mfrow = c(1,1))

# Quick plot reconstruction
par(mfrow = c(1,1))
x_id = which(substr(varnames(m_spline_norm),1,2)=="x[")
x_post = as.matrix(m_spline_norm[ ,x_id])
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


if(!dir.exists("Results/Figures/Detrend/Splines/")) dir.create("Results/Figures/Detrend/Splines/", recursive = TRUE)

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

#####

#### climate spline 50 year knots - not working  #####

knots <- seq(0, Tea, by = 50)
B <- bs(1:Tea, knots = knots) # neet to sort x_use? then doesn't line up with rest.
H <- length(knots) + 3

matplot(1:Tea, B, type="l",xlab="Time",ylab="Basis function, Bj(X)",cex.lab=1.5,cex.axis=1.5,lwd=2)
# d <- max(4, floor(Tea/35))
# K <- floor(Tea/d - 1)

# visualize priors (D draws from prior distribution)
D <- 10
g <- matrix(0, Tea, D)

for(j in 1:D){
  beta  <- rnorm(length(knots) + 3, 0, 1)
  g[ , j] <- B%*%beta
}

matplot(1:Tea, g, lwd=2, type="l", cex.lab=1.5, cex.axis=1.5)

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
clusterExport(cl, c("spl_data", "initialize_spl", "params", "M", "f", "l", "Tea", "a_use", "log_y", "x_use", "nb", "ni", "nt", "B", "H"))
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
x_id_50 = which(substr(varnames(m_spline_50),1,2)=="x[") # finds the indices of the x variables
post_climate_50 = colMeans(as.matrix(m_spline_50[,x_id_50])) # finds the posterior mean of the x variables
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




#### linear detrend climate spline 50 year knots #####

knots <- seq(0, Tea, by = 50)
B <- bs(1:Tea, knots = knots) # neet to sort x_use? then doesn't line up with rest.
H <- length(knots) + 3

matplot(1:Tea, B, type="l",xlab="Time",ylab="Basis function, Bj(X)",cex.lab=1.5,cex.axis=1.5,lwd=2)
# d <- max(4, floor(Tea/35))
# K <- floor(Tea/d - 1)

# visualize priors (D draws from prior distribution)
D <- 10
g <- matrix(0, Tea, D)

for(j in 1:D){
  beta  <- rnorm(length(knots) + 3, 0, 1)
  g[ , j] <- B%*%beta
}

matplot(1:Tea, g, lwd=2, type="l", cex.lab=1.5, cex.axis=1.5)

# run model
initialize_spl = function(){
  alpha0 = rlnorm(M, rowMeans(log_y, na.rm=TRUE), 0.15)
  alpha1 = runif(M, -0.1, 0)
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

spl_data <- list(y = y_orig, 
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
clusterExport(cl, c("spl_data", "initialize_spl", "params", "M", "f", "l", "Tea", "a_use", "log_y", "x_use", "nb", "ni", "nt", "B", "H"))
clusterSetRNGStream(cl = cl, 86753091)
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
m_linear_spline_50 <- mcmc.list(out)
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
x_id_50 = which(substr(varnames(m_spline_50),1,2)=="x[") # finds the indices of the x variables
post_climate_50 = colMeans(as.matrix(m_spline_50[,x_id_50])) # finds the posterior mean of the x variables
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
recon <- plot_recon(m_spline_50, obs = climate)
ggsave(filename = "Results/Figures/wash_recon_post_spl50.tiff", plot = recon, width = 8, height = 4, units = "in") # , dpi = 1000) # poster vs paper formatting - see past work and make package or github source

recon_valid <- plot_recon(m_spline_50, obs = climate, valid_yrs = year[hold_out])
ggsave(filename = "Results/Figures/wash_recon_spl50_valid_back.tiff", plot = recon_valid, width = 8, height = 4, units = "in") 




##### linear detrend linear climate ####

initialize_m2_nc = function(){
  alpha0 = rnorm(M, rowMeans(y_orig, na.rm=TRUE), 0.25)
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

m2_nc_data <- list(y = y_orig, 
                   f = f, 
                   l = l, 
                   M = M, 
                   Tea = Tea, 
                   a = a_use, 
                   # v = 410, # or 320 
                   x = x_use)

params <- c("x", "mu_a0", "mu_a1", "sd_a0", "sd_a1", "beta0", "sd_eta", "sd_x", "sd_y")

cl <- makeCluster(nc)                       # Request # cores
clusterExport(cl, c("m2_nc_data", "initialize_m2_nc", "params", "M", "f", "l", "Tea", "a_use", "y_orig", "x_use", "nb", "ni", "nt"))
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
m_linear_linear <- mcmc.list(out)


plot(m_linear_linear[ ,c("mu_a0", "mu_a1", "sd_a0", "sd_a1", "sd_eta", "sd_x", "beta0")])
par(mfrow = c(1,1))
x_id_50 = which(substr(varnames(m_linear_linear),1,2)=="x[") # finds the indices of the x variables
post_climate_50 = colMeans(as.matrix(m_linear_linear[,x_id_50])) # finds the posterior mean of the x variables
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
recon <- plot_recon(m_linear_linear, obs = tree_climate)

ggsave(filename = "Results/Figures/wash_recon_post_linear_const.tiff", plot = recon, width = 8, height = 4, units = "in") # , dpi = 1000) # poster vs paper formatting - see past work and make package or github source

recon_valid <- plot_recon(m_linear_linear, valid_yrs = year[hold_out])
ggsave(filename = "Results/Figures/wash_recon_linear_const_valid_back.tiff", plot = recon_valid, width = 8, height = 4, units = "in") 




#####