setwd('/Users/jimblotter/Desktop/Grad_School/Data_Analysis/erisor/QC/') #set path
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

########### San Francisco Mountain Range Enlarged (ggmap) #######################

#Get blank map
SFMRmap <- get_map(location = c(-113.77, 38.35, -113, 38.87),
               color = "color",
               source = "google",
               maptype = "terrain", #roadmap? hybrid? terrain
               zoom = 10)

#Plot points
gpsSFMR <- read.csv("gps_erisor_SFMR.csv", header = TRUE, stringsAsFactors = FALSE)
gpslake <- read.csv("gps_erisor_SFMR.csv", header = TRUE, stringsAsFactors = FALSE)
gpslake[13, "Population"] <- "Sevier Lake"
gpslake[13, "Latitude"] <- 38.9
gpslake[13, "Longitude"] <- -113.1
gpslake[13, "Species"] <- "Sevier Lake"   #I can remove this label, I was just messing around.
gpslake

ggmap(SFMRmap) +
  geom_point(data = gpsSFMR, aes(x = gpsSFMR$Longitude, y = gpsSFMR$Latitude,
  colour = gpsSFMR$Species,
  fill = gpsSFMR$Species,
  size = 0.5, shape = 20)) + scale_shape_identity() +
  geom_label_repel(aes(x = gpslake$Longitude, y = gpslake$Latitude, label = gpslake$Population),
            data = gpslake) 

########################## GGPLOT2 attempt (main distribution map) ###########################################
#read in data, assign usa
usa <- map_data("usa")
gps2 <- read.csv("gps_erisor.csv", header = TRUE)

#get insert image
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
  

#plot points, add repelled labels, add image for those unfamiliar in the corner
erisor_dist <- gg1 +
  coord_fixed(xlim = c(-119.39, -107.0), ylim = c(34, 43), ratio = 1.3) +
  geom_point(data = gps, aes(x = gps2$Longitude, y = gps2$Latitude,
    colour = gps2$Species,
    fill = gps2$Species), size = 3.5, shape = 20) + 
  scale_shape_identity() +
  geom_label_repel(aes(x = gps2$Longitude, y = gps2$Latitude,
    label=gps2$Population), nudge_x = -0.2) +
  annotation_raster(not_fam,ymin = 33.4,ymax= 36.1,xmin = -120,xmax = -113.69) 
  scalebar(data = gps, location = "topright", dist = 50, st.dist = 0.2, 
           st.bottom = TRUE, st.size = 5, ddkm = TRUE, model = "WGS84", 
           x.min = , x.max, y.min, y.max, anchor = NULL, 
           facet.var = NULL, facet.lev = NULL)

########################## Left to fix ################################
#a border (maybe with lat long markers?)
#a distance scale
#a North arrow
#some major cities?
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

http://egallic.fr/scale-bar-and-north-arrow-on-a-ggplot2-map/
