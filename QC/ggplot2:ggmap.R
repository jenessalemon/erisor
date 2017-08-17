setwd('/Users/jimblotter/Desktop/Grad_School/Data_Analysis/erisor/QC/') #set path
library("ggplot2")
library("ggmap")
library("Imap")
library("maps")
gps <- read.csv("gps_erisor.csv", header = TRUE)            #read in csv with population coordinates

lat_up <- max(gps$Latitude) +1                                 #map boundaries
lat_down <- min(gps$Latitude) -1
long_left <- min(gps$Longitude) -1
long_right <- max(gps$Longitude) +1

#Get blank map
map <- get_map(location = c(long_left, lat_down, long_right, lat_up),
               color = "color",
               source = "google",
               maptype = "hybrid", #roadmap? hybrid? terrain
               zoom = 6)

ggmap(map) +                                         
  geom_point(data = gps, aes(x = gps$Longitude, y = gps$Latitude, 
                             colour = gps$Species,
                             fill = gps$Species,
                             size = 0.5, shape = 21)) + scale_shape_identity()

ggmap(map) +                                         
  geom_point(data = gps, aes(x = gps$Longitude, y = gps$Latitude, 
                             colour = gps$Region, 
                             fill = gps$Region, 
                             alpha = 0.8, 
                             size = 2))

#map of regions
c <- c("UTAH", "NEVADA", "IDAHO", "ARIZONA", "COLORADO", "WYOMING", "NEW MEXICO", "CALIFORNIA")
map(database = "state")
map(database = "state",regions = c,col = "blue",fill=T,add=TRUE)

#map of region
gg1 <- ggplot() + 
  geom_polygon(data = usa, aes(x=long, y = lat, group = group), fill = NA, color = "black") + 
  coord_fixed(1.3)

labs <- data.frame(
  long = c(-119.0975,-106.9835,-106.9835,-119.0975),
  lat = c(43.89989,43.89989,34.18568,34.18568),
  stringsAsFactors = FALSE
)  

gg1 + 
  #geom_path(data = usa, aes(x=long, y = lat, group = group), fill = NA, color = "blue") + 
  geom_polygon(data = labs, aes(x = long, y = lat), color = "blue", size = 1, fill = NA)



