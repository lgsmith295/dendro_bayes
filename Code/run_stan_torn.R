#####################################################################
# Hierarchical reconstruction of temperature from Tornetrask Scots pine tree rings
#
# Daniel J. Hocking
# Based on Schofield et al. 2016, 2017 and Steinscheider et al. 2017
#####################################################################

###### Load libraries #####
library(dplyr)
# install.packages("devtools")
# devtools::install_github("mjskay/tidybayes")
# devtools::install_github("tidyverse/tidyr")
library(tidybayes)
library(tidyverse)
library(ggplot2)
library(ggfan)
library(rstan)

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

#### Load and Prep Data #####
# load data
a_use <- as.matrix(read.table('Data/Tornetrask/age.txt', header = FALSE))
year <- c(read.table("Data/Tornetrask/year.txt"), recursive = TRUE) 
x_full <- c(read.table("Data/Tornetrask/xtorn.txt"), recursive=TRUE) # climate data

# Combine year and climate data
climate_df <- data.frame(year, x_full, stringsAsFactors = FALSE) %>%
  dplyr::rename(temp = x_full)
# instramental temperature record from 1816-1995 (179 years)

# Get tree ring data
y_orig <- as.matrix(read.table('Data/Tornetrask/ytorn.txt', header = FALSE))

y_long <- t(y_orig) %>%
  as.data.frame() %>%
  gather(key = tree, value = rwl) %>%
  mutate(year = rep(year, times = nrow(y_orig)),
         year_index = rep(1:500, times = nrow(y_orig))) %>%
  dplyr::filter(!is.na(rwl))

a_long <- t(a_use) %>%
  as.data.frame() %>%
  gather(key = tree, value = age) %>%
  mutate(year = rep(year, times = nrow(a_use)))

df_full <- y_long %>%
  left_join(a_long) %>%
  left_join(climate_df) %>%
  tidyr::separate(col = tree, into = c(NA, "tree_id"), sep = "V", remove = FALSE) %>%
  mutate(x_s = as.numeric(scale(temp)),
         x_obs = ifelse(is.na(.$temp), 0, 1),
         tree_id = as.integer(tree_id),
         log_y = log(rwl)) %>%
  dplyr::arrange(x_obs, tree_id, year)
  
# standardize climate data
x_mean <- mean(df_full$x_s, na.rm=TRUE)
x_sd <- sd(df_full$x_s, na.rm=TRUE);
# x_s <- (x_full - x_mean) / x_sd # standardize temperatures to z-scores for numerical stability
x_use <- df_full[which(df_full$x_obs == 1), ]$x_s

# Hold out years for validation (first half)
hold_out <- 321:410  # the indices that we hold out
# x_use[hold_out] <- NA

log_y <- log(df_full$rwl)

M <- nrow(y_orig)
Tea <- ncol(y_orig)

l <- f <- rep(NA, M) # set up blank first and last year vectors for each tree
for(i in 1:M) {
  tmp <- which(!is.na(y_orig[i, ])) # column (year) with first non-NA
  f[i] <- min(tmp)
  l[i] <- max(tmp)
}

#### Does creating a standard chronology FORCE the climate reconstruction to be (normally) distributed with fluctuations around a mean????? Over flatten??? #####

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

if(!dir.exists("Results/Stan")) dir.create("Results/Stan", recursive = TRUE)

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

df_obs <- df_full %>%
  filter(x_obs == 1)

df_miss <- df_full %>%
  filter(x_obs != 1)

dat_stan <- list(N_obs = nrow(df_obs),
                 x_obs = df_obs$x_s,
                 N_miss = nrow(df_miss), # problem in estimating the same values of climate many times
                 log_y = df_full$log_y,
                 M = M,
                 N = nrow(df_full),
                 tree_id = df_full$tree_id,
                 year = df_full$year_index,
                 age = df_full$age)


params <- c("x_miss", "mu_a0", "mu_a1", "sd_a0", "sd_a1", "beta0", "sd_eta", "sd_x", "sd_y")

m_negexp_linear <- stan(file = "Code/Stan/negexp_linear.stan", 
                        data = dat_stan, 
                        pars = params,
                        iter = 1000,
                        warmup = 500,
                        thin = 1)
