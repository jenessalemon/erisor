setwd('/Users/jimblotter/Desktop/Grad_School/Data_Analysis/erisor/QC/') #set path
library("ggplot2")
library("ggmap")
library("Imap")

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

#Geographic distance matrix
places_names <- c("Museum of Modern Art New York, NY",
                  "Smithsonian Museum of American Art Washington, DC",
                  "Brooklyn Museum Brooklyn, NY",
                  "Walker Art Center Minneapolis, MN",
                  "Fralin Museum of Art Charlottesville, VA")
dist_list <- lapply(places_names, 
                    function(z) 
                      sapply(z, 
                             function(x) 
                               mapdist(from=places_names, to=x)$miles)
)
dist_list
#BUT this is limited to 250 queries a day. Use Imap:


places_names <- c("Museum of Modern Art New York, NY",
                  "Smithsonian Museum of American Art Washington, DC",
                  "Brooklyn Museum Brooklyn, NY",
                  "Walker Art Center Minneapolis, MN",
                  "Fralin Museum of Art Charlottesville, VA")

# geocode place names
places_lat <- geocode(places_names, source="google")$lat
places_lon <- geocode(places_names, source="google")$lon

# create a data frame to store all variables
places_df <- data.frame(names = places_names,
                        lat = places_lat,
                        lon = places_lon)

# calculate geodesic distance with gdist() from Imap package

# create an empty list
dist_list <- list()

# iterate through data frame placing calculated distance next to place place names
for (i in 1:nrow(places_df)) {
  
  dist_list[[i]] <- gdist(lon.1 = places_df$lon[i], 
                          lat.1 = places_df$lat[i], 
                          lon.2 = places_df$lon, 
                          lat.2 = places_df$lat, 
                          units="miles")
  
}

# view results as list
dist_list