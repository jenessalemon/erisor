setwd('/Users/jimblotter/Desktop/Grad_School/Data_Analysis/erisor/QC/') #set path
library("ggplot2")
library("ggmap")
library("Imap")
library("maps")
library("mapdata")
library("plyr")
library("magick")
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
gg1 <- ggplot() +
geom_polygon(data = usa, aes(x=long, y = lat, group = group), fill = NA, color = "black") +
coord_fixed(1.3)

bounds <- data.frame(
long = c(-119.0975,-106.9835,-106.9835,-119.0975),
lat = c(43.89989,43.89989,34.18568,34.18568),
stringsAsFactors = FALSE
)

gg1 +
geom_polygon(data = bounds, aes(x = long, y = lat), color = "blue", size = 1, fill = NA) +
borders("state", colour = "black")

#map of states collected in
c <- c("UTAH", "NEVADA", "IDAHO", "ARIZONA", "COLORADO", "WYOMING", "NEW MEXICO", "CALIFORNIA")
map(database = "state")
map(database = "state",regions = c,col = "blue",fill=T,add=TRUE)

#######################################################################################################
#San Francisco Mountain Range Enlarged
#Get blank map
map <- get_map(location = c(-113.7, 38.3, -113.1, 38.81),
color = "color",
source = "google",
maptype = "terrain", #roadmap? hybrid? terrain
zoom = 10)
#Plot points
ggmap(map) +
geom_point(data = gps, aes(x = gps$Longitude, y = gps$Latitude,
colour = gps$Species,
fill = gps$Species,
size = 0.5, shape = 20)) + scale_shape_identity()

########################################################################################
#Another way
usa <- map_data("usa")
gg1 <- ggplot() +
geom_polygon(data = usa, aes(x=long, y = lat, group = group), fill = NA, color = "black") +
coord_fixed(1.3) + borders("state", colour = "black")

gg1 +
coord_fixed(xlim = c(-119, -107.0), ylim = c(34, 44), ratio = 1.3) +
geom_point(data = gps, aes(x = gps$Longitude, y = gps$Latitude,
colour = gps$Species,
fill = gps$Species,
size = 1, shape = 20)) + scale_shape_identity() +
geom_text(aes(x = gps$Longitude, y = gps$Latitude,
label=gps$Population), #overlapping
hjust=0, vjust=1)

