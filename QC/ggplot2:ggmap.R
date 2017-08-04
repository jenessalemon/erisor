setwd('/Users/jimblotter/Desktop/Grad_School/Data_Analysis/erisor/QC/') #remember to put the .str file there!
library("ggplot2")
library("ggmap")

gps <- read.csv("gps_erisor.csv", header = TRUE)            #read in csv with population coordinates
lat_up <- max(gps$Latitude)                                 #map boundaries
lat_down <- min(gps$Latitude)
long_left <- min(gps$Longitude)
long_right <- max(gps$Longitude)

#Get blank map
map <- get_map(location = c(long_left, lat_down, long_right, lat_up),
               color = "color",
               source = "google",
               maptype = "hybrid", #roadmap? hybrid? terrain
               zoom = 6)
#Add points! Need to investigate coloring the dots. I've seen Aaron Duffy do it in a code somewhere.
ggmap(map) +                                         
  geom_point(data = gps, aes(x = gps$Longitude, y = gps$Latitude, fill = "red", alpha = 0.8), size = 2, shape = 21) +
  guides(col = gps$Population) #I think Aaron has written code for this, check it out!
