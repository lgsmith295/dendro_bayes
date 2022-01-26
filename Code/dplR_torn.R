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