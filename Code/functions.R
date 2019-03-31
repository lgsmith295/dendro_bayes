theme_bw_journal <- function (base_family = "") {
  theme_grey(base_family = base_family) %+replace%
    theme(
      axis.text = element_text(size = rel(0.9)),
      axis.title = element_text(size = rel(0.9)),
      axis.ticks = element_line(colour = "black"),
      legend.key = element_rect(colour = "grey80"),
      panel.background = element_rect(fill = "white",
                                      colour = NA),
      panel.border = element_rect(fill = NA,
                                  colour = "grey50"),
      panel.grid.major = element_line(colour = "grey90",
                                      size = 0.2),
      panel.grid.minor = element_line(colour = "grey98",
                                      size = 0.5),
      strip.background = element_rect(
        fill = "grey80",
        colour = "grey50",
        size = 0.2
      )
    )
}

theme_set(theme_bw_journal())

theme_bw_poster <- function (base_family = "") {
  theme_grey(base_family = base_family) %+replace%
    theme(
      # legend.position = "none",
      title = element_text(size = 20),
      axis.text = element_text(size = 20),
      axis.title = element_text(size = 25),
      axis.ticks = element_line(colour = "black"),
      legend.key = element_rect(colour = "grey80"),
      legend.title=element_text(size=20), 
      legend.text=element_text(size=15),
      panel.background = element_rect(fill = "white",
                                      colour = NA),
      panel.border = element_rect(fill = NA,
                                  colour = "grey50"),
      panel.grid.major = element_line(colour = "grey90",
                                      size = 0.2),
      panel.grid.minor = element_line(colour = "grey98",
                                      size = 0.5),
      strip.background = element_rect(
        fill = "grey80",
        colour = "grey50",
        size = 0.2
      )
    )
}

theme_set(theme_bw_poster())

plot_recon <- function(mcmc, obs = climate, mean = x_mean, sd = x_sd, valid_yrs = NULL, y_lab = "Mean Temperature (C)") {
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
    g <- g + geom_fan() + geom_line(data = climate_valid, aes(year, temp), colour = "red") + ylab(y_lab) + xlab("Year") + scale_fill_distiller() + theme_bw()
  } else {
    # scaled posterior interval
    g <- ggplot(temp_df_long, aes(x = year, y = temp))
    g <- g + geom_fan() + geom_line(data = obs, aes(x = year, y = value), colour="black", size = 0.2) + ylab(y_lab) + xlab("Year") + theme_bw() + scale_fill_distiller() #x_full),
  }
  
  return(g)
}


# reconstruction plot
plot_recon_temp <- function(mcmc, obs = climate_df, x_mean = x_mean, x_sd = x_sd, valid_yrs = NULL) {
  xidx_m2 = which(substr(varnames(mcmc),1,2)=="x[")
  temp_df <- as_tibble(t(as.matrix(mcmc[ , xidx_m2])))
  temp_df <- temp_df %>% 
    mutate(year = year)
  
  temp_df_long <- temp_df %>% 
    gather(key = sim, value = temp, -year) %>%
    dplyr::mutate(temp = temp*x_sd + x_mean)
  
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

# Performance statistics
perf_stats <- function(est_valid, observed, valid_id = hold_out, cal_id = cal_ids, mod_id = "m") {
  x_valid <- observed[valid_id]
  R2 <- sum((est_valid - mean(x_valid))^2) / sum((x_valid - mean(x_valid))^2)
  Pearson_R2 <- (sum((x_valid - mean(x_valid)) * (est_valid - mean(est_valid))))^2 / (sum((x_valid - mean(x_valid))^2) * sum((est_valid - mean(est_valid))^2))
  RE <- 1 - sum((x_valid - est_valid) ^ 2) / sum((x_valid - mean(x_full[cal_id], na.rm = TRUE))^2)
  CE <- 1 - sum((x_valid - est_valid) ^ 2) / sum((x_valid - mean(x_full[valid_id])) ^ 2)
  
  df <- data.frame(model = mod_id, R2, Pearson_R2, RE, CE, stringsAsFactors = FALSE)
  return(df)
}


# create observed vs expected plot with CRI?????
if(FALSE) {
xidx_m2 = which(substr(varnames(mcmc),1,2)=="x[")
temp_df <- as_tibble(t(as.matrix(mcmc[ , xidx_m2])))
temp_df <- temp_df %>% 
  mutate(year = obs$year)

temp_df_long <- temp_df %>% 
  gather(key = sim, value = temp, -year) %>%
  dplyr::mutate(temp = temp*sd + mean)

# Validation plot
  temp_valid <- temp_df_long %>%
    dplyr::filter(year %in% valid_yrs) %>%
    dplyr::mutate(Value = "estimated")
  
  climate_valid <- obs %>%
    dplyr::mutate(Value = "observed") %>%
    dplyr::rename(temp = value) %>% # x_full) %>%
    dplyr::filter(year %in% valid_yrs)
  
  df <- left_join(climate_valid, temp_valid, by = c("year"))
  
  ggplot(df, aes(temp.x, temp.y)) + geom_point(alpha = 0.01) # + geom_fan()
  
  g <- ggplot(df, aes(x = as.factor(temp.x), y = temp.y)) + geom_boxplot()
  
}