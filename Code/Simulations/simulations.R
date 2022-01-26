
###### Load libraries #####
library(rjags) # for running JAGS
library(MASS) # for boxcox
library(lattice) # for xyplot
library(dplyr)
# install.packages("devtools")
# devtools::install_github("mjskay/tidybayes")
# library(tidybayes)
library(tidyverse)
library(ggplot2)
library(ggfan)
library(parallel)

library(devtools)
install_github(repo = "lgsmith295/simDendro")
library(simDendro)

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

##### NegExp Growth - Mean Climate #####
years <- 1:500

# Growth
growth <- negexp_growth(length(years), 1, 0.03, 0.001)
plot(growth, type = "b")

# climate random, uncorrelated around mean
mu_climate <- 12
climate <- rnorm(length(years), mu_climate, 3)
plot(climate, type = "l")

# Climate effect on growth
beta_mu <- 0.8 # real scale slope - 1 degree C increase in T results in this change in mm tree ring growth
# log_beta <- log(beta_mu)
beta <- rnorm(length(years), beta_mu, 0.05)
eta <- beta * climate
plot(eta, type = "l")
plot(climate, eta)

# IID error
eps <- rnorm(length(years), -2.5, 0.5)
error <- exp(eps)
summary(error)
plot(error, type = "l")

# Put growth model together
rwl <- growth * eta * error
plot(rwl, type = "b")


# loop over multiple trees
M <- 100 # number of trees
Tea <- 500
years <- 1500 + 1:Tea

# Growth
tgrowth <- matrix(NA, nrow = Tea, ncol = M)
# assume each tree lives 200 years for simplicity and germinated in a random year between 1350 and 1950 - assumes can't get to pith in older trees
fy <- floor(runif(M, 1351, 1951))
k <- rlnorm(M, log(0.2), 0.2) # tree biological growth curves vary
b <- rlnorm(M, log(0.04), 0.1) # tree biological growth curves vary
a <- rlnorm(M, log(1), 0.2) # tree biological growth curves vary - maybe should be multivariate normal between the parameters
hist(a)
hist(b)
hist(k)
sd(b)

for(i in 1:M) {
  bio <- negexp_growth(200, a[i], b[i], k[i])
  yrs <- fy[i]:(fy[i] + 199)
  df_tmp <- data.frame(tree = i, year = yrs, bio)
  if(i == 1) {
    df <- df_tmp
  } else {
    df <- bind_rows(df, df_tmp)
  }
}

str(df)

# climate random, uncorrelated around mean
mu_climate <- 12
climate <- rnorm(Tea, mu_climate, 3)
plot(climate, type = "l")

# Climate effect on growth
mu_beta <- 0.8 # real scale slope - 1 degree C increase in T results in this change in mm tree ring growth
# log_beta <- log(beta_mu)
beta <- rnorm(Tea, mu_beta, 0.05) # this should probably vary by tree but I guess that's going into other variance components
eta <- beta * climate
plot(eta, type = "l")
plot(climate, eta)
eta <- data.frame(year = years, eta, climate)

# IID error
eps <- rnorm(Tea, -2.5, 0.5)
error <- exp(eps)
summary(error)
plot(error, type = "l")

# Random IID tree variance among years
sigma_eps_i <- rlnorm(M, log(0.5), 0.05)
hist(sigma_eps_i)
eps_ij <- data.frame(expand.grid(tree = 1:M, year = years, eps = NA_real_)) %>%
  dplyr::arrange(tree, year)
for(i in 1:M) {
  eps_ij[which(eps_ij$tree == i), ]$eps <- rlnorm(Tea, log(0.08), sigma_eps_i[i])
}
hist(eps_ij$eps)
# tree_sd <- data.frame(tree = 1:M, eps_i = eps_i)

# Put model together
rwl_long <- df %>%
    group_by(tree) %>%
    left_join(eta) %>%
    left_join(eps_ij) %>%
    mutate(rwl = bio * eta * eps,
           rwi = rwl - bio,
           climate_std = (climate - mean(climate, na.rm = TRUE)) / sd(climate, na.rm = TRUE)) %>%
    dplyr::filter(year %in% years)

ggplot(rwl_long, aes(year, rwl)) + geom_line(alpha = 0.3, color = "blue") + theme(legend.position = "none") + theme_bw()

ggplot(rwl_long, aes(year, rwi)) + geom_point(alpha = 0.1, color = "turquoise") + theme(legend.position = "none") + theme_bw()

ggplot(rwl_long, aes(year, rwi)) + 
  geom_point(alpha = 0.1, color = "turquoise", aes(group = tree)) + 
  geom_smooth() + 
  geom_line(aes(year, climate_std), alpha = 0.1) + 
  geom_smooth(aes(year, climate_std), color = "red", alpha = 0.5) + 
  theme(legend.position = "none") + 
  theme_bw()

# just show 1st 10 trees
rwl_small <- rwl_long %>%
  dplyr::filter(tree %in% 11:20)

ggplot(rwl_small, aes(year, rwl, group = tree, color = as.factor(tree))) + geom_line(alpha = 0.3) + theme_bw() + theme(legend.position = "none")

# Fit model
y_ij <- rwl_long %>%
  group_by(tree) %>%
  dplyr::select(tree, year, rwl) %>%
  spread(year, rwl) %>%
  ungroup() %>%
  dplyr::select(-tree)

# standardize climate data
x_mean <- mean(climate, na.rm=TRUE)
x_sd <- sd(climate, na.rm=TRUE);
x_s <- (climate - x_mean) / x_sd # standardize temperatures to z-scores for numerical stability
x_use <- x_s

# Hold out years for validation (first half)
hold_out <- 1:400  # the indices that we hold out
x_use[hold_out] <- NA

# Get tree ring data
log_y <- log(as.matrix(y_ij))

M <- nrow(log_y)
Tea <- ncol(log_y)

l <- f <- rep(NA, M) # set up blank first and last year vectors for each tree
for(i in 1:M) {
  tmp <- which(!is.na(log_y[i, ])) # column (year) with first non-NA
  f[i] <- min(tmp)
  l[i] <- max(tmp)
}

# Age to use
# Can't know age because of growth form of these trees so can't do RCS but can just do individual standardizations with year index as age
age <- matrix(NA, nrow = nrow(log_y), ncol = ncol(log_y))
for(i in 1:M) {
  columns <- f[i]:l[i]
  ages <- 1:length(columns)
  age[i, columns] <- ages
}
a_use <- (age - mean(age, na.rm = TRUE)) / sd(age, na.rm = TRUE)

##### negexp detrend linear climate M2-non-centered ####

# my code should have the same output as M2 from Schofield I think

initialize = function(){
  alpha0 = rnorm(M, rowMeans(log_y, na.rm=TRUE), 0.25)
  alpha1 = -rlnorm(M, -1, 0.25)
  sd_y = rlnorm(M, 0, 1) 
  eta = rnorm(Tea, 0, 0.25)
  sd_eta = sd(eta) 
  beta0 = rnorm(1, 0, 0.1)
  sd_x = rlnorm(1, 0, 1) 
  return(list(sd_y = sd_y,
              alpha0 = alpha0,
              alpha1 = alpha1,
              sd_eta = sd_eta,
              beta0 = beta0,
              sd_x = sd_x))
}

m_data <- list(y = log_y, 
                   f = f, 
                   l = l, 
                   M = M, 
                   Tea = Tea, 
                   a = a_use, 
                   x = x_use)

params <- c("x", "beta0", "sd_eta", "sd_x", "sd_y")

cl <- makeCluster(nc)                       # Request # cores
clusterExport(cl, c("m_data", "initialize", "params", "M", "f", "l", "Tea", "a_use", "log_y", "x_use", "nb", "ni", "nt"))
clusterSetRNGStream(cl = cl, 867541)
system.time({ 
  out <- clusterEvalQ(cl, {
    library(rjags)
    jm <- jags.model("Code/JAGS/ind_negexp_mean_climate.txt", m_data, initialize, n.adapt=nb, n.chains=1) # Compile model and run burnin
    out <- coda.samples(jm, params, n.iter=ni, thin=nt) # Sample from posterior distribution
    return(as.mcmc(out))
  })
}) # 

stopCluster(cl)

# Results
outM2_nc <- mcmc.list(out)
plot(outM2_nc[ , c("beta0", "sd_eta", "sd_x")])
par(mfrow = c(1,1))
# summary(outM2_nc[,c("mu_a0", "mu_a1", "sd_a0", "sd_a1", "beta0", "sd_eta", "sd_x")])

save(outM2_nc, file = "Results/JAGS/sim_mean_climate.RData")

xidx_m2 = which(substr(varnames(outM2_nc),1,2)=="x[")
postxm_m2 = colMeans(as.matrix(outM2_nc[,xidx_m2])) 
x_post = as.matrix(outM2_nc[ , xidx_m2])
x_post_median <- apply(x_post, 2, median)

plot(x_post_median)

plot(years, postxm_m2*x_sd + x_mean,type="l") # plots the posterior mean of the x variables on the original scale (degrees C) and year on x-axis

plot(postxm_m2*x_sd + x_mean, climate)




##### climate with linear trend and random noise #####

# Growth
growth <- negexp_growth(length(years), 1, 0.03, 0.1)
plot(growth, type = "b")

# climate
mu_climate <- 10 + 0.1 * years
climate <- rnorm(length(years), mu_climate, 3)
plot(climate, type = "l")

# Climate effect on growth
beta_mu <- 0.8 # real scale slope - 1 degree C increase in T results in this change in mm tree ring growth
# log_beta <- log(beta_mu)
beta <- rnorm(length(years), beta_mu, 0.05)
eta <- beta * climate
plot(eta, type = "l")
plot(climate, eta)

# IID error
eps <- rnorm(length(years), -2.5, 0.5)
error <- exp(eps)
summary(error)
plot(error, type = "l")

# Put growth model together
rwl <- growth * eta * error
plot(rwl, type = "b")


##### Allow climate to vary over time #####

# Growth
growth <- negexp_growth(length(years), 1, 0.03, 0.1)
plot(growth, type = "b")

# Climate: autoregressive model around set mean
mu_climate <- 12
yt <- arima.sim(n=500, list(order=c(1,0,0),ar=0.9), sd = 0.5) + mu_climate  # adding intercept and trend to the autoregressive 
plot(yt)
acf(yt)

### Climate: autoregressive model with linear trend
mu_climate <- 12
yt <- arima.sim(n=500, list(order=c(1,0,0),ar=0.9), sd = 0.5) + mu_climate + years*0.01  # adding intercept and trend to the autoregressive 
plot(yt)
acf(yt)


### Climate: basis splines
mu_climate <- 12

library(mgcv)
x <- seq(0,500,length=501)
knots1 <- seq(0, 500, by = 100)
sm <- smoothCon(s(x,bs="bs"), data.frame(x), knots = data.frame(knots1))[[1]] # addition of knots function not working and not throwing warning or error
matplot(1:501, sm$X, type="l",xlab="Time",ylab="Basis function, Bj(X)",cex.lab=1.5,cex.axis=1.5,lwd=2)

set.seed(416923)
beta  <- rnorm(ncol(sm$X), 0, 3)
g <- (sm$X %*% beta) + mu_climate
plot(g, type = "l")

# add random noise around trend
climate <- g + rnorm(length(g), 0, 2)
plot(climate, type = "l")
acf(climate)

# loop over multiple trees
M <- 100 # number of trees
Tea <- 500
years <- 1500 + 1:Tea

# Growth

# assume each tree lives 200 years for simplicity and germinated in a random year between 1350 and 1950 - assumes can't get to pith in older trees
fy <- floor(runif(M, 1351, 1951))
k <- rlnorm(M, log(0.2), 0.2) # tree biological growth curves vary
b <- rlnorm(M, log(0.04), 0.1) # tree biological growth curves vary
a <- rlnorm(M, log(1), 0.2) # tree biological growth curves vary - maybe should be multivariate normal between the parameters
hist(a)
hist(b)
hist(k)
sd(b)

for(i in 1:M) {
  bio <- negexp_growth(200, a[i], b[i], k[i])
  yrs <- fy[i]:(fy[i] + 199)
  df_tmp <- data.frame(tree = i, year = yrs, bio)
  if(i == 1) {
    df <- df_tmp
  } else {
    df <- bind_rows(df, df_tmp)
  }
}

str(df)

# climate spline
mu_climate <- 12

library(mgcv)
x <- years
knots1 <- seq(0, max(years), by = 100)
sm <- smoothCon(s(x,bs="bs"), data.frame(x), knots = data.frame(knots1))[[1]] # addition of knots function not working and not throwing warning or error

set.seed(416923)
beta  <- rnorm(ncol(sm$X), 0, 3)
g <- (sm$X %*% beta) + mu_climate
# plot(g, type = "l")

# add random noise around trend
climate <- g + rnorm(length(g), 0, 0.5)
# plot(climate, type = "l")
# acf(climate)

# Climate effect on growth
mu_beta <- 0.8 # real scale slope - 1 degree C increase in T results in this change in mm tree ring growth
# log_beta <- log(beta_mu)
beta <- rnorm(Tea, mu_beta, 0.05) # this should probably vary by tree but I guess that's going into other variance components
eta <- beta * climate
# plot(eta, type = "l")
# plot(climate, eta)
eta <- data.frame(year = years, eta, climate)

# IID error
eps <- rnorm(Tea, -2.5, 0.5)
error <- exp(eps)
# summary(error)
# plot(error, type = "l")

# Random IID tree variance among years
sigma_eps_i <- rlnorm(M, log(0.5), 0.05)
hist(sigma_eps_i)
eps_ij <- data.frame(expand.grid(tree = 1:M, year = years, eps = NA_real_)) %>%
  dplyr::arrange(tree, year)
for(i in 1:M) {
  eps_ij[which(eps_ij$tree == i), ]$eps <- rlnorm(Tea, log(0.08), sigma_eps_i[i])
}
# hist(eps_ij$eps)
# tree_sd <- data.frame(tree = 1:M, eps_i = eps_i)

# Put model together
rwl_long <- df %>%
  group_by(tree) %>%
  left_join(eta) %>%
  left_join(eps_ij) %>%
  mutate(rwl = bio * eta * eps,
         rwi = rwl - bio,
         rwi_div = rwl / bio,
         climate_std = (climate - mean(climate, na.rm = TRUE)) / sd(climate, na.rm = TRUE)) %>%
  dplyr::filter(year %in% years)

ggplot(rwl_long, aes(year, rwl)) + geom_line(alpha = 0.3, color = "blue") + theme(legend.position = "none") + theme_bw()

ggplot(rwl_long, aes(year, rwi)) + geom_point(alpha = 0.1, color = "turquoise") + theme(legend.position = "none") + theme_bw()

ggplot(rwl_long, aes(year, rwi)) + 
  geom_point(alpha = 0.1, color = "turquoise", aes(group = tree)) + 
  geom_smooth() + 
  geom_line(aes(year, climate_std), alpha = 0.1) + 
  geom_smooth(aes(year, climate_std), color = "red", alpha = 0.5) + 
  theme(legend.position = "none") + 
  theme_bw()

ggplot(rwl_long, aes(year, rwi_div)) + 
  geom_point(alpha = 0.1, color = "turquoise", aes(group = tree)) + 
  geom_smooth() + 
  geom_line(aes(year, climate_std), alpha = 0.1) + 
  geom_smooth(aes(year, climate_std), color = "red", alpha = 0.5) + 
  theme(legend.position = "none") + 
  theme_bw()

# just show 1st 10 trees
rwl_small <- rwl_long %>%
  dplyr::filter(tree %in% 11:20)

ggplot(rwl_small, aes(year, rwl, group = tree, color = as.factor(tree))) + geom_line(alpha = 0.3) + theme_bw() + theme(legend.position = "none")

# Fit model
y_ij <- rwl_long %>%
  group_by(tree) %>%
  dplyr::select(tree, year, rwl) %>%
  spread(year, rwl) %>%
  ungroup() %>%
  dplyr::select(-tree)

# standardize climate data
x_mean <- mean(climate, na.rm=TRUE)
x_sd <- sd(climate, na.rm=TRUE);
x_s <- (climate[ , 1] - x_mean) / x_sd # standardize temperatures to z-scores for numerical stability
x_use <- x_s

# Hold out years for validation (first half)
hold_out <- 1:400  # the indices that we hold out
x_use[hold_out] <- NA

# Get tree ring data
log_y <- log(as.matrix(y_ij))

M <- nrow(log_y)
Tea <- ncol(log_y)

l <- f <- rep(NA, M) # set up blank first and last year vectors for each tree
for(i in 1:M) {
  tmp <- which(!is.na(log_y[i, ])) # column (year) with first non-NA
  f[i] <- min(tmp)
  l[i] <- max(tmp)
}

# Age to use
# Can't know age because of growth form of these trees so can't do RCS but can just do individual standardizations with year index as age
age <- matrix(NA, nrow = nrow(log_y), ncol = ncol(log_y))
for(i in 1:M) {
  columns <- f[i]:l[i]
  ages <- 1:length(columns)
  age[i, columns] <- ages
}
a_use <- (age - mean(age, na.rm = TRUE)) / sd(age, na.rm = TRUE)

##### negexp detrend linear climate M2-non-centered ####

# my code should have the same output as M2 from Schofield I think

initialize = function(){
  alpha0 = rnorm(M, rowMeans(log_y, na.rm=TRUE), 0.25)
  alpha1 = -rlnorm(M, -1, 0.25)
  sd_y = rlnorm(M, 0, 1) 
  eta = rnorm(Tea, 0, 0.25)
  sd_eta = sd(eta) 
  beta0 = rnorm(1, 0, 0.1)
  sd_x = rlnorm(1, 0, 1) 
  return(list(sd_y = sd_y,
              alpha0 = alpha0,
              alpha1 = alpha1,
              sd_eta = sd_eta,
              beta0 = beta0,
              sd_x = sd_x))
}

m_data <- list(y = log_y, 
               f = f, 
               l = l, 
               M = M, 
               Tea = Tea, 
               a = a_use, 
               x = x_use)

params <- c("x", "beta0", "sd_eta", "sd_x", "sd_y")

cl <- makeCluster(nc)                       # Request # cores
clusterExport(cl, c("m_data", "initialize", "params", "M", "f", "l", "Tea", "a_use", "log_y", "x_use", "nb", "ni", "nt"))
clusterSetRNGStream(cl = cl, 867541)
system.time({ 
  out <- clusterEvalQ(cl, {
    library(rjags)
    jm <- jags.model("Code/JAGS/ind_negexp_mean_climate.txt", m_data, initialize, n.adapt=nb, n.chains=1) # Compile model and run burnin
    out <- coda.samples(jm, params, n.iter=ni, thin=nt) # Sample from posterior distribution
    return(as.mcmc(out))
  })
}) # 

stopCluster(cl)

# Results
outM2_nc <- mcmc.list(out)
plot(outM2_nc[ , c("beta0", "sd_eta", "sd_x")])
par(mfrow = c(1,1))
# summary(outM2_nc[,c("mu_a0", "mu_a1", "sd_a0", "sd_a1", "beta0", "sd_eta", "sd_x")])

save(outM2_nc, file = "Results/JAGS/sim_mean_climate.RData")

xidx_m2 = which(substr(varnames(outM2_nc),1,2)=="x[")
postxm_m2 = colMeans(as.matrix(outM2_nc[,xidx_m2])) 
x_post = as.matrix(outM2_nc[ , xidx_m2])
x_post_median <- apply(x_post, 2, median)

plot(x_post_median)

plot(years, postxm_m2*x_sd + x_mean,type="l") # plots the posterior mean of the x variables on the original scale (degrees C) and year on x-axis

plot(postxm_m2*x_sd + x_mean, climate)


# run with climate spline
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

m_spline_25 = jags.model('Code/JAGS/ind_negexp_climate_spline.txt', 
                         list(y = log_y, 
                              f = f, 
                              l = l, 
                              M = M, 
                              Tea = Tea, 
                              a = a_use, 
                              x = x_use,
                              B = B,
                              H = H), 
                         inits = initialize, 
                         n.chains = 3, 
                         n.adapt = 1000)

out_m_climate_spline_25 = coda.samples(m_spline_25, c("x", "sd_eta", "sd_x", "sd_y", "mu_gamma", "sd_g"), 1000)

plot(out_m_climate_spline_25[ ,c("sd_eta", "sd_x")])
par(mfrow = c(1,1))
x_idx_25 = which(substr(varnames(out_m_climate_spline_25),1,2)=="x[") # finds the indices of the x variables
post_climate_25 = colMeans(as.matrix(out_m_climate_spline_25[ ,x_idx_25])) # finds the posterior mean of the x variables
plot(post_climate_25,type="l") # plots the posterior mean of the x variables
plot(years, post_climate_25*x_sd + x_mean, type="l") # plots the posterior mean of the x variables on the original scale (degrees C) and year on x-axis
plot(years, climate, type = "l")

plot(post_climate_25*x_sd + x_mean, climate)
abline(0, 1, col = "red")

# validation
x_valid_post_m <- post_climate_25[hold_out]*x_sd + x_mean
x_valid <- x_full[hold_out]
plot(x_valid, x_valid_post_m, type = "p")
abline(0, 1, col = "red")
cor(x_valid, x_valid_post_m)

res2 <- (x_valid_post_m - x_valid)^2
sqrt(sum(res2) / length(x_valid_post_m))
sqrt(mean(res2))

##### negexp detrend 1 changepoint climate - seems good ####

x_min <- min(x_use, na.rm = TRUE) # value of x on standardized scale when climate = 0
x_max <- max(x_use, na.rm = TRUE)

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
  x_1 = runif(1, x_min, max(x_use, na.rm = TRUE))
  sd_x = rlnorm(1, 0, 1) 
  return(list(sd_y = sd_y,
              alpha0 = alpha0,
              alpha1 = alpha1,
              # mu_a0 = mu_a0,
              # mu_a1 = mu_a1,
              # sd_a0 = sd_a0,
              # sd_a1 = sd_a1,
              sd_eta = sd_eta,
              # beta0 = beta0,
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

plot(m_negexp_1change[ ,c("sd_eta", "sd_x", "beta0[1]", "beta0[2]", "x_1")])
par(mfrow = c(1,1))
x_id_50 = which(substr(varnames(m_negexp_1change),1,2)=="x[") # finds the indices of the x variables
post_climate_50 = colMeans(as.matrix(m_negexp_1change[,x_id_50])) # finds the posterior mean of the x variables
# plot(post_climate_50,type="l") # plots the posterior mean of the x variables
plot(years, post_climate_50*x_sd + x_mean,type="l") # plots the posterior mean of the x variables on the original scale (degrees C) and year on x-axis


## Create series for multiple trees with the same climate


## Partial Samples (not near pith) - maybe move this to own repository using the package and separate from the vignette

# Simulate all tree growth using a biological growth model but then cut off the samples to more recent years due to things like heart rot, alternative growth forms, and missing the pith. See how different models assuming different growth (e.g. negex, spline, linear, mean) recover the data and reconstruct climate.

