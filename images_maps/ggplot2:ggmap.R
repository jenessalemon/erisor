setwd('/Users/jimblotter/Desktop/Grad_School/Data_Analysis/erisor/images_maps/') #set path
library("ggplot2")
library("ggmap")
library("Imap")
library("maps")
library("mapdata")
library("plyr")
library("magick")
library("ggrepel")
library("ggsn")
gps <- read.csv("gps_erisor.csv", header = TRUE)            #read in csv with population coordinates

###################### GGMAP attempt (main distribution map) #####################################
#determine map boundaries
lat_up <- max(gps$Latitude) +1                                 
lat_down <- min(gps$Latitude) -1
long_left <- min(gps$Longitude) -1
long_right <- max(gps$Longitude) +1

#set boudnaries for SFMR box
bounds_SFMR <- data.frame(
  long = c(-113.9,-112.95,-112.95,-113.9),
  lat = c(39,39,38.25,38.25),
  stringsAsFactors = FALSE)

#Get blank map
map <- get_map(location = c(long_left, lat_down, long_right, lat_up),
          color = "color",
          source = "google",
          maptype = "terrain", #roadmap? hybrid? terrain
          zoom = 6) + coord_fixed(1.3)

#read in insert image
not_familiar <- image_read(path = "/Users/jimblotter/Desktop/Grad_School/Data_Analysis/erisor/images_maps/not_familiar.png")
not_fam <- as.raster(not_familiar) #convert to raster

#Plot points over the blank map previously created (map)
ggmap(map) +
          geom_point(data = gps, x = gps$Longitude, y = gps$Latitude,
            colour = "black",
            size = 2, shape = 20) + scale_shape_identity() + 
          geom_label_repel(aes(x = gps2$Longitude, y = gps2$Latitude, label = gps2$Population),
                   data = gps2, nudge_x = -0.1) +
          geom_polygon(data = bounds_SFMR, aes(x = long, y = lat), color = "black", size = 0.75, fill = NA) + 
          coord_cartesian(ylim = c(lat_down, lat_up), xlim=c(long_left, long_right)) +
          annotation_raster(not_fam,ymin = 34,ymax= 36.8,xmin = long_left,xmax = -113.5)


############### Map for those unfamiliar with USA (ggpolot2) ###############################
#get blank map
usa <- map_data("usa")
  ggblank <- ggplot() +
  geom_polygon(data = usa, aes(x=long, y = lat, group = group), fill = "white", color = "black") +
  coord_fixed(1.3)
  
#set bounds for the box enclosing the range
bounds_notfamiliar <- data.frame(
long = c(-119.0975,-106.9835,-106.9835,-119.0975),
lat = c(43.89989,43.89989,34.18568,34.18568),
stringsAsFactors = FALSE)

#plot bounds over map
ggblank +
geom_polygon(data = bounds_notfamiliar, aes(x = long, y = lat), color = "blue", size = 1, fill = NA) +
borders("state", colour = "black") + theme_nothing()

################### San Francisco Mountain Range Enlarged (ggmap) #######################

#Get blank map
SFMRmap <- get_map(location = c(-113.6, 38.46, -113.3, 38.75),
               color = "color",
               source = "google",
               maptype = "terrain", #roadmap? hybrid? terrain
               zoom = 10)
ggmap(SFMRmap)


#Needed to plot points
gpsSFMR <- read.csv("gps_erisor_SFMR.csv", header = TRUE, stringsAsFactors = FALSE)
gpslake <- read.csv("gps_erisor_SFMR.csv", header = TRUE, stringsAsFactors = FALSE)
gpslake[13, "Population"] <- "Sevier Lake"
gpslake[13, "Latitude"] <- 38.9
gpslake[13, "Longitude"] <- -113.1
gpslake[13, "Species"] <- "Sevier Lake"   #I can remove this label, I was just messing around.
gpslake

#Needed for scalebar:
bb <- attr(SFMRmap, "bb")
bb
bb2 <- data.frame(long = unlist(bb[c(2, 4)]), lat = unlist(bb[c(1,3)]))

mymap <- ggmap(SFMRmap) + 
  geom_point(data = gpsSFMR, aes(x = gpsSFMR$Longitude, y = gpsSFMR$Latitude,
            size = 0.5, 
            fill= gpsSFMR$Species,
            shape = gpsSFMR$Species),
            show.legend = TRUE) +
  guides(color = FALSE, size = FALSE) +
  theme(legend.title=element_blank(), + legend.text=element_text(size=5))  +
  scale_shape_manual(values = c(18,111)) +
  scale_shape_identity() +
  geom_label_repel(aes(x = gpslake$Longitude, y = gpslake$Latitude, label = gpslake$Population), data = gpslake) +
  scalebar(data = bb2, dist = 10, dd2km = TRUE, model  = "WGS84", 
           location = "topleft", 
           anchor = c(
             x = bb$ll.lon + 0.08 * (bb$ur.lon - bb$ll.lon), 
             y = bb$ll.lat + 0.95 * (bb$ur.lat - bb$ll.lat)
           ) 
  )
print(mymap + scale_shape_manual(values = c(111,18)))
#fill = guide_legend(title = "Species", override.aes = list(alpha = 1))
  #coord_fixed(ratio = 1:1) +       #in my opinion the N is redundant because the y axis is increasing in latitude.
  #annotation_raster(north,ymin = 38.85,ymax= 38.94,xmin = -113.81,xmax = -113.7)

########################## GGPLOT2 attempt (main distribution map) ###########################################
#read in data, assign usa
usa <- map_data("usa")
gps2 <- read.csv("gps_erisor.csv", header = TRUE)

#get insert image
north_arrow <- image_read(path = "/Users/jimblotter/Desktop/Grad_School/Data_Analysis/erisor/images_maps/north-arrow-2.png")
north <- as.raster(north_arrow) #convert to raster

#get north symbol
not_familiar <- image_read(path = "/Users/jimblotter/Desktop/Grad_School/Data_Analysis/erisor/images_maps/not_familiar.png")
not_fam <- as.raster(not_familiar) #convert to raster

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
  geom_polygon(data = usa, aes(x=long, y = lat, group = group), fill = "white", color = "black") +
  coord_fixed(1.3) + borders("state", colour = "black") +
  geom_polygon(data = bounds_SFMR, aes(x = long, y = lat), color = "grey", size = 1, fill = NA) +
    borders("state", colour = "black")
  
#For scalebar

#plot points, add repelled labels, add image for those unfamiliar in the corner
gg1 +
  coord_fixed(xlim = c(-119.39, -107.0), ylim = c(34, 43), ratio = 1.3) +
  geom_point(data = gps, aes(x = gps2$Longitude, y = gps2$Latitude,
    colour = gps2$Species,
    fill = gps2$Species), size = 3.5, shape = 20) + 
  scale_shape_identity() +
  geom_label_repel(aes(x = gps2$Longitude, y = gps2$Latitude,
    label=gps2$Population), nudge_x = -0.2) +
  annotation_raster(not_fam,ymin = 33.4,ymax= 36.1,xmin = -120,xmax = -113.69) +
  annotation_raster(north,ymin = 42.15,ymax= 43.4,xmin = -120,xmax = -118.5) +
  scalebar(data = usa, dist = 10, dd2km = TRUE, model  = "WGS84", 
           location = "bottomright")

########################## Left to fix ################################
#a distance scale on ggplot distribution map
#edit the legend? (both)
#change color of the points to stand out? (both)
#a North arrow (on ggmap if necessary)
#a border (maybe with lat long markers?) I have this?
#some major cities?
http://egallic.fr/scale-bar-and-north-arrow-on-a-ggplot2-map/

library("ggmap")
library("ggsn")
map <- get_googlemap(center =c(23.89, 52.74), zoom = 12, maptype = "hybrid")
bb <- attr(map, "bb")
bb2 <- data.frame(long = unlist(bb[c(2, 4)]), lat = unlist(bb[c(1,3)]))
ggmap(map) + 
  scalebar(data = bb2, dist = 50, dd2km = TRUE, model  = "WGS84", 
           location = "topleft", 
           anchor = c(
             x = bb$ll.lon + 0.1 * (bb$ur.lon - bb$ll.lon), 
             y = bb$ll.lat + 0.9 * (bb$ur.lat - bb$ll.lat)
           )
  )
