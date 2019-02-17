################################################################################
# Hierarchical Modeling for Dendroclimatology
# Daniel J. Hocking and Laura G. Smith
# 
# 12 December 2018
# 
# Building off of Schofield et al. 2016
################################################################################

#------------ Load libraries ------------------

library(dplyr)
library(lubridate)
library(ggplot2)
library(ggmap)
library(maps)
library(mapdata)
library(devtools)
# devtools::install_github("ropensci/FedData")
library(FedData)

#------------ Pull Data from ITRDB ------------

# Create data directory if doesn't exist
if(!file.exists("Data")) {
  dir.create("Data", showWarnings = FALSE)
}

# Create bounding box of area from which to pull the data
vepPolygon <- polygon_from_extent(raster::extent(-80.925, -73.235, 37.195, 42.634), proj4string='+proj=longlat +datum=WGS84 +ellps=WGS84 +zone=17')

norway_polygon <- polygon_from_extent(raster::extent(4.99207807783, 58.0788841824, 31.29341841, 80.6571442736), proj4string='+proj=longlat +datum=WGS84 +ellps=WGS84 +zone=17')
  

# Get the ITRDB records
ITRDB <- get_itrdb(template = norway_polygon, label=NULL, species = c("PISY"), makeSpatial = TRUE, force.redo = TRUE, measurement.type = "Ring Width", chronology.type = "Measurements Only") # , species = "JUVI", - none found but maybe JUSP?, , measurement.type = "Ring Width", template = vepPolygon, , extraction.dir = "Data/ITRDB/"

ITRDB <- get_itrdb(template = NULL, label='NORW015', makeSpatial = TRUE, force.redo = TRUE)


ITRDB <- get_itrdb(template = vepPolygon, label='VA017', makeSpatial = TRUE, force.redo = FALSE) # , species = "JUVI", - none found but maybe JUSP?, , measurement.type = "Ring Width", template = vepPolygon, , extraction.dir = "Data/ITRDB/"

str(ITRDB)

# check what species extracted
unique(ITRDB$metadata@data$SPECIES)

# check who contributed (Cook?)
unique(ITRDB$metadata@data$CONTRIBUTOR)

# check series
unique(ITRDB$metadata@data$SERIES)

# check names
unique(ITRDB$metadata@data$NAME)


# Create map
usa <- map_data("usa")
states <- map_data("state")
tva_states <- subset(states, region %in% c("tennessee", "kentucky"))

ggplot(data = states) + 
  geom_polygon(aes(x = long, y = lat, group = group), fill = "gray90", color = "black") +
  geom_point(data = ITRDB$metadata@data, aes(x = LON, y = LAT, group = SERIES), color = "red") + 
  coord_fixed(1.3) +
  theme_bw() +
  guides(fill=FALSE)  # do this to leave off the color legend