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

############ Map for those unfamiliar with USA (ggpolot2) ###############################
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
geom_polygon(data = bounds_notfamiliar, aes(x = long, y = lat), color = "grey", size = 1, fill = NA) +
borders("state", colour = "black") + theme_nothing()

################### San Francisco Mountain Range Enlarged (ggmap) #######################

#Get blank map
SFMRmap <- get_map(location = c(-113.6, 38.55, -113.4, 38.7),
               color = "color",
               source = "google",
               maptype = "terrain", #roadmap? hybrid? terrain
               zoom = 10)

#Needed to plot points
gpsSFMR <- read.csv("gps_erisor_SFMR.csv", header = TRUE, stringsAsFactors = FALSE)
gpsSFMR[11,] <- NA
gpsSFMR[12,] <- NA
gpslake <- read.csv("gps_erisor_SFMR.csv", header = TRUE, stringsAsFactors = FALSE)
gpslake[13, "Population"] <- "Sevier Lake"
gpslake[13, "Latitude"] <- 38.9
gpslake[13, "Longitude"] <- -113.1
gpslake[13, "Species"] <- "Sevier Lake"   #I can remove this label, I was just messing around.
gpslake[11,] <- NA
gpslake[12,] <- NA
gpslake

#Needed for scalebar:
bb <- attr(SFMRmap, "bb")
bb
bb2 <- data.frame(long = unlist(bb[c(2, 4)]), lat = unlist(bb[c(1,3)]))

mymap <- ggmap(SFMRmap) + 
  geom_point(data = gpsSFMR, aes(x = gpsSFMR$Longitude, y = gpsSFMR$Latitude,
            size = 0.3, 
            shape = 20,
            fill= gpsSFMR$Species,
            color = gpsSFMR$Species),
            show.legend = FALSE) +
  #scale_shape_manual(values = c(18,20)) + REMOVED SO THATI AM CONSISTANT
  scale_color_manual(values = c("#66CDAA", "#FF7F50")) +
  scale_shape_identity() +
  geom_label_repel(data = gpslake, aes(x = gpslake$Longitude, y = gpslake$Latitude, 
          label = gpslake$Population)) +
  scalebar(data = bb2, dist = 10, dd2km = TRUE, model  = "WGS84", 
           location = "topleft", 
           anchor = c(
             x = bb$ll.lon + 0.08 * (bb$ur.lon - bb$ll.lon), 
             y = bb$ll.lat + 0.95 * (bb$ur.lat - bb$ll.lat)))

print(mymap + scale_fill_manual(values = c("#66CDAA", "#FF7F50")))

#fill = guide_legend(title = "Species", override.aes = list(alpha = 1))
  #coord_fixed(ratio = 1:1) +       #in my opinion the N is redundant because the y axis is increasing in latitude.
  #annotation_raster(north,ymin = 38.85,ymax= 38.94,xmin = -113.81,xmax = -113.7)

########################## GGPLOT2 (main distribution map) ###########################################
#read in data, assign usa
usa <- map_data("usa")
gps2 <- read.csv("gps_erisor.csv", header = TRUE)

#get north insert image
north_arrow <- image_read(path = "/Users/jimblotter/Desktop/Grad_School/Data_Analysis/erisor/images_maps/north-arrow-2.png")
north <- as.raster(north_arrow) #convert to raster

#get not familiar
not_familiar <- image_read(path = "/Users/jimblotter/Desktop/Grad_School/Data_Analysis/erisor/images_maps/box_black.png")
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

#plot points, add repelled labels, add image for those unfamiliar in the corner
gg1 +
  coord_fixed(xlim = c(-119.39, -107.0), ylim = c(34, 43), ratio = 1.3) +
  geom_point(data = gps, aes(x = gps2$Longitude, y = gps2$Latitude,
    colour = gps2$Species,
    fill = gps2$Species), size = 3.5, shape = 20) + 
  #scale_color_manual(values = wes_palette("FantasticFox", 3, type = "discrete")) +
  scale_color_manual(values = c("#66CDAA", "#66CDAA", "#FF7F50")) +
  scale_shape_identity() +
  geom_label_repel(aes(x = gps2$Longitude, y = gps2$Latitude,
    label=gps2$Population), nudge_x = -0.2) +
  annotation_raster(north,ymin = 42.15,ymax= 43.4,xmin = -119.95,xmax = -118.5) +
  annotation_raster(not_fam,ymin = 33.9,ymax= 36.1,xmin = -120,xmax = -114.19) +
  scalebar(data = usa, dist = 10, dd2km = TRUE, model  = "WGS84", 
           location = "bottomright") +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1.5)) +
  theme(legend.position="none") +
  scale_bar(lon = -112.5, lat = 42.5, 
              distance_lon = 200, distance_lat = 20, distance_legend = 40, 
              dist_unit = "km", orientation = FALSE)  

#shockleyi color: "#009E73"
#soredium color: "#E69F00"
########################## Left to fix ################################
#a distance scale on ggplot distribution map
#edit the legend? (both)
#change color of the points to stand out? (both)
#a North arrow (on ggmap if necessary)
#a border (maybe with lat long markers?) I have this?
#some major cities?
http://egallic.fr/scale-bar-and-north-arrow-on-a-ggplot2-map/

library(maps)
library(maptools)
library(ggplot2)
library(grid)

create_scale_bar <- function(lon,lat,distance_lon,distance_lat,distance_legend, dist_units = "km"){
  # First rectangle
  bottom_right <- gcDestination(lon = lon, lat = lat, bearing = 90, dist = distance_lon, dist.units = dist_units, model = "WGS84")
  
  topLeft <- gcDestination(lon = lon, lat = lat, bearing = 0, dist = distance_lat, dist.units = dist_units, model = "WGS84")
  rectangle <- cbind(lon=c(lon, lon, bottom_right[1,"long"], bottom_right[1,"long"], lon),
                     lat = c(lat, topLeft[1,"lat"], topLeft[1,"lat"],lat, lat))
  rectangle <- data.frame(rectangle, stringsAsFactors = FALSE)
  
  # Second rectangle t right of the first rectangle
  bottom_right2 <- gcDestination(lon = lon, lat = lat, bearing = 90, dist = distance_lon*2, dist.units = dist_units, model = "WGS84")
  rectangle2 <- cbind(lon = c(bottom_right[1,"long"], bottom_right[1,"long"], bottom_right2[1,"long"], bottom_right2[1,"long"], bottom_right[1,"long"]),
                      lat=c(lat, topLeft[1,"lat"], topLeft[1,"lat"], lat, lat))
  rectangle2 <- data.frame(rectangle2, stringsAsFactors = FALSE)
  
  # Now let's deal with the text
  on_top <- gcDestination(lon = lon, lat = lat, bearing = 0, dist = distance_legend, dist.units = dist_units, model = "WGS84")
  on_top2 <- on_top3 <- on_top
  on_top2[1,"long"] <- bottom_right[1,"long"]
  on_top3[1,"long"] <- bottom_right2[1,"long"]
  
  legend <- rbind(on_top, on_top2, on_top3)
  legend <- data.frame(cbind(legend, text = c(0, distance_lon, distance_lon*2)), stringsAsFactors = FALSE, row.names = NULL)
  return(list(rectangle = rectangle, rectangle2 = rectangle2, legend = legend))
}

create_orientation_arrow <- function(scale_bar, length, distance = 1, dist_units = "km"){
  lon <- scale_bar$rectangle2[1,1]
  lat <- scale_bar$rectangle2[1,2]
  
  # Bottom point of the arrow
  beg_point <- gcDestination(lon = lon, lat = lat, bearing = 0, dist = distance, dist.units = dist_units, model = "WGS84")
  lon <- beg_point[1,"long"]
  lat <- beg_point[1,"lat"]
  
  # Let us create the endpoint
  on_top <- gcDestination(lon = lon, lat = lat, bearing = 0, dist = length, dist.units = dist_units, model = "WGS84")
  
  left_arrow <- gcDestination(lon = on_top[1,"long"], lat = on_top[1,"lat"], bearing = 225, dist = length/5, dist.units = dist_units, model = "WGS84")
  
  right_arrow <- gcDestination(lon = on_top[1,"long"], lat = on_top[1,"lat"], bearing = 135, dist = length/5, dist.units = dist_units, model = "WGS84")
  
  res <- rbind(
    cbind(x = lon, y = lat, xend = on_top[1,"long"], yend = on_top[1,"lat"]),
    cbind(x = left_arrow[1,"long"], y = left_arrow[1,"lat"], xend = on_top[1,"long"], yend = on_top[1,"lat"]),
    cbind(x = right_arrow[1,"long"], y = right_arrow[1,"lat"], xend = on_top[1,"long"], yend = on_top[1,"lat"]))
  
  res <- as.data.frame(res, stringsAsFactors = FALSE)
  
  # Coordinates from which "N" will be plotted
  coords_n <- cbind(x = lon, y = (lat + on_top[1,"lat"])/2)
  
  return(list(res = res, coords_n = coords_n))
}

scale_bar <- function(lon, lat, distance_lon, distance_lat, distance_legend, dist_unit = "km", rec_fill = "white", rec_colour = "black", rec2_fill = "black", rec2_colour = "black", legend_colour = "black", legend_size = 3, orientation = TRUE, arrow_length = 500, arrow_distance = 300, arrow_north_size = 6){
  the_scale_bar <- create_scale_bar(lon = lon, lat = lat, distance_lon = distance_lon, distance_lat = distance_lat, distance_legend = distance_legend, dist_unit = dist_unit)
  # First rectangle
  rectangle1 <- geom_polygon(data = the_scale_bar$rectangle, aes(x = lon, y = lat), fill = rec_fill, colour = rec_colour)
  
  # Second rectangle
  rectangle2 <- geom_polygon(data = the_scale_bar$rectangle2, aes(x = lon, y = lat), fill = rec2_fill, colour = rec2_colour)
  
  # Legend
  scale_bar_legend <- annotate("text", label = paste(the_scale_bar$legend[,"text"], dist_unit, sep=""), x = the_scale_bar$legend[,"long"], y = the_scale_bar$legend[,"lat"], size = legend_size, colour = legend_colour)
  
  res <- list(rectangle1, rectangle2, scale_bar_legend)
  
  if(orientation){# Add an arrow pointing North
    coords_arrow <- create_orientation_arrow(scale_bar = the_scale_bar, length = arrow_length, distance = arrow_distance, dist_unit = dist_unit)
    arrow <- list(geom_segment(data = coords_arrow$res, aes(x = x, y = y, xend = xend, yend = yend)), annotate("text", label = "N", x = coords_arrow$coords_n[1,"x"], y = coords_arrow$coords_n[1,"y"], size = arrow_north_size, colour = "black"))
    res <- c(res, arrow)
  }
  return(res)
}

usa_map <- map_data("state")
P <- ggplot() + geom_polygon(data = usa_map, aes(x = long, y = lat, group = group)) + coord_map()

P + scale_bar(lon = -130, lat = 26, 
              distance_lon = 500, distance_lat = 100, distance_legend = 200, 
              dist_unit = "km", orientation = FALSE)
