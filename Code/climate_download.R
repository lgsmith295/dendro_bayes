############################################
# Import Monthly Climate Data from PRISM
############################################

# Adapted from http://eremrah.com/articles/How-to-extract-data-from-PRISM-raster/

library(tidyr)
library(stringr)
library(prism)
library(raster)
library(dplyr)
library(ggplot2)

if(!dir.exists("Output")) {
  dir.create("Output")
}

##### Set Conditions #####
# Start with growing season streamflow and precip
mo_start <- 5 # start month - integer
mo_end <- 10 # end moneth - integer
res <- "4km" # resolution 800m or 4 km  - only 4km monthly available for free

##### Use Metadata from ITRDB to get coordinates #####
library(rjson)

# Check metadata
meta_files <- list.files("Data/itrdb_pilo_mount_washington/metadata/")
meta_list <- list()

for(i in 1:length(meta_files)) {
  meta_list[[i]] <- fromJSON(file = paste0("Data/itrdb_pilo_mount_washington/metadata/", meta_files[i]))
  names(meta_list)[[i]] <- meta_list[[i]]$studyCode
}

str(meta_list)

coords <- meta_list$NV520$site[[1]]$geo$geometry$coordinates

##### Import Monthly Climate #####
ptm <- proc.time() # time how long this takes - start timer
options(prism.path = "Data/Prism") # set location for downloaded data
get_prism_monthlys(type = "ppt", years = 1895:2017, mon = mo_start:mo_end, keepZip = FALSE)
proc.time() - ptm # 28 minutes

temp <- proc.time() # time how long this takes - start timer
options(prism.path = "Data/Prism") # set location for downloaded data
get_prism_monthlys(type = "tmean", years = 1895:2017, mon = mo_start:mo_end, keepZip = FALSE)
proc.time() - temp # 28 minutes

# check list
# ls_prism_data(name = TRUE)
str(ls_prism_data(name = TRUE))

# Create a Raster Stack
precip_stack <- ls_prism_data() %>%  
  prism_stack(.)  # Stack files

# Get proj from raster stack
prism_crs <- precip_stack@crs@projargs

# Study Locations
locs <- data.frame(id = c(1),
                       lat = c(as.numeric(coords[1])),
                       lon = c(as.numeric(coords[2])))

# Convert points to spatial points data frame
coordinates(locs) <- c('lon', 'lat')
proj4string(locs) <- CRS(prism_crs)

# Extract data from raster (very slow for huge stack)
df <- data.frame(coordinates(locs), locs$id, extract(precip_stack, locs))

# Reshape data
climate <- df %>%  
  gather(date, value, 4:ncol(df))

# Split date into type, year, month
climate <- separate(climate, "date", c("t1", "type", "stable", "scale", "YearMonth", "end"), sep = c("_"))
climate <- separate(climate, "YearMonth", c("year", "month"), sep = c(4)) %>%
  dplyr::select(-t1, -stable, -scale, -end) %>%
  dplyr::mutate(year = as.integer(year),
         month = as.integer(month))

str(climate)
unique(climate$type)

# Reshape data again
climate2 <- climate %>%
  spread(type, value)

panel.cor <- function(x, y, digits=2, prefix="", cex.cor, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- cor(x, y, use="complete.obs")
  txt <- format(c(r, 0.123456789), digits=digits)[1]
  txt <- paste(prefix, txt, sep="")
  if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt, cex = cex.cor*0.5, col=c("gray60", "black")[(abs(r)>0.65)+1])
}

panel.hist <- function(x, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(usr[1:2],0,1.5) )
  h <- hist(x, plot = FALSE)
  breaks <- h$breaks; nB <- length(breaks)
  y <- h$counts; y <- y/max(y)
  rect(breaks[-nB], 0, breaks[-1], y, col="cyan", ...)
}

pairs(climate2[ , c("ppt", "tmean")],
      upper.panel = panel.cor, diag.panel = panel.hist)

# 
# # Order data
# df <- df[order(df$locs.id), ]


# Save
save(climate2, file = paste0("Data/prism_", names(meta_list)[[1]], ".RData"))






