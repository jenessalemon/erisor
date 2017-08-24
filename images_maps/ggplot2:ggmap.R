setwd('/Users/jimblotter/Desktop/Grad_School/Data_Analysis/erisor/QC/') #set path
library("ggplot2")
library("ggmap")
library("Imap")
library("maps")
library("mapdata")
library("plyr")
library("magick")
library("ggrepel")
gps <- read.csv("gps_erisor.csv", header = TRUE)            #read in csv with population coordinates

lat_up <- max(gps$Latitude) +1                                 #map boundaries
lat_down <- min(gps$Latitude) -1
long_left <- min(gps$Longitude) -1
long_right <- max(gps$Longitude) +1

#Get blank map
map <- get_map(location = c(long_left, lat_down, long_right, lat_up),
          color = "color",
          source = "google",
          maptype = "terrain", #roadmap? hybrid? terrain
          zoom = 6)

not_familiar <- image_read(path = "/Users/jimblotter/Desktop/Grad_School/Data_Analysis/erisor/images_maps/not_familiar.png")
not_fam <- as.raster(not_familiar)

#Plot points
ggmap(map) +
          geom_point(data = gps, aes(x = gps$Longitude, y = gps$Latitude,
          colour = gps$Species,
          fill = gps$Species,
          size = 0.5, shape = 21)) + scale_shape_identity() + 
          coord_cartesian(ylim = c(lat_down, lat_up), xlim=c(long_left, long_right)) +
          annotation_raster(not_fam,ymin = lat_down,ymax= 37.2,xmin = long_left,xmax = -113)

######################################################################################################
#Map for those unfamiliar with USA
usa <- map_data("usa")
ggblank <- ggplot() +
geom_polygon(data = usa, aes(x=long, y = lat, group = group), fill = NA, color = "black") +
coord_fixed(1.3)

bounds_notfamiliar <- data.frame(
long = c(-119.0975,-106.9835,-106.9835,-119.0975),
lat = c(43.89989,43.89989,34.18568,34.18568),
stringsAsFactors = FALSE
)

ggblank +
geom_polygon(data = bounds_notfamiliar, aes(x = long, y = lat), color = "blue", size = 1, fill = NA) +
borders("state", colour = "black")

#######################################################################################################
#San Francisco Mountain Range Enlarged
#Get map for insert: (but I've decided not to insert it)
gg1 <- ggplot() +
  geom_polygon(data = usa, aes(x=long, y = lat, group = group), fill = NA, color = "black") +
  coord_fixed(1.3) + borders("state", colour = "black")

bounds_SFMR <- data.frame(
  long = c(-113.9,-112.9,-112.9,-113.9),
  lat = c(39,39,38.2,38.2),
  stringsAsFactors = FALSE
)

gg1 +
  coord_fixed(xlim = c(-118.5, -106.5), ylim = c(35, 43), ratio = 1.3) +
  geom_point(data = gps, aes(x = gps2$Longitude, y = gps2$Latitude,
                             colour = gps2$Species,
                             fill = gps2$Species,
                             size = 1, shape = 20)) + 
                             scale_shape_identity() +
  geom_polygon(data = bounds_SFMR, aes(x = long, y = lat), color = "blue", size = 1, fill = NA) +
  borders("state", colour = "black")

where_SFMR <- image_read(path = "/Users/jimblotter/Desktop/Grad_School/Data_Analysis/erisor/images_maps/SFMR_where_in_distribution.png")
where_SFMR_is <- as.raster(where_SFMR)

#Get blank map
SFMRmap <- get_map(location = c(-113.77, 38.35, -113, 38.87),
               color = "color",
               source = "google",
               maptype = "terrain", #roadmap? hybrid? terrain
               zoom = 10)

#Plot points
gpsSFMR <- read.csv("gps_erisor_SFMR.csv", header = TRUE, stringsAsFactors = FALSE)

ggmap(SFMRmap) +
  geom_point(data = gpsSFMR, aes(x = gpsSFMR$Longitude, y = gpsSFMR$Latitude,
  colour = gpsSFMR$Species,
  fill = gpsSFMR$Species,
  size = 0.5, shape = 20)) + scale_shape_identity() +
  geom_label_repel(aes(x = gpsSFMR$Longitude, y = gpsSFMR$Latitude,
                       label = gpsSFMR$Population))

########################################################################################
#Another way
usa <- map_data("usa")
gps2 <- read.csv("gps_erisor.csv", header = TRUE)
#get rid of the SFMR labels
gps2[1,"Population"] <- NA
gps2[2,"Population"] <- NA
gps2[3,"Population"] <- NA
gps2[5,"Population"] <- NA
gps2[6,"Population"] <- NA
gps2[7,"Population"] <- NA
gps2[8,"Population"] <- NA
gps2[9,"Population"] <- NA
gps2[29,"Population"] <- NA
gps2[31,"Population"] <- NA
gps2[33,"Population"] <- NA
gps2[34,"Population"] <- NA
#define bounds for blue box to show the SFMR
bounds_SFMR <- data.frame(
  long = c(-113.9,-112.9,-112.9,-113.9),
  lat = c(39,39,38.2,38.2),
  stringsAsFactors = FALSE
)
#get blank map
gg1 <- ggplot() +
  geom_polygon(data = usa, aes(x=long, y = lat, group = group), fill = NA, color = "black") +
  coord_fixed(1.3) + borders("state", colour = "black") +
  geom_polygon(data = bounds_SFMR, aes(x = long, y = lat), color = "blue", size = 1, fill = NA) +
    borders("state", colour = "black")

gg1 +
  coord_fixed(xlim = c(-119, -107.0), ylim = c(34, 44), ratio = 1.3) +
  geom_point(data = gps, aes(x = gps2$Longitude, y = gps2$Latitude,
    colour = gps2$Species,
    fill = gps2$Species,
    size = 1, shape = 20)) + 
  scale_shape_identity() +
  geom_label_repel(aes(x = gps2$Longitude, y = gps2$Latitude,
    label=gps2$Population)) +
  annotation_raster(not_fam,ymin = 33.5,ymax= 36.2,xmin = -119,xmax = -113)

