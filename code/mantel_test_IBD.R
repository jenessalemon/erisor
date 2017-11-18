# Input: .str file, and a .csv file with two columns, "lon" and "lat" for each individual in the .str file.
# In this code, we accomplish the following:
# 1 Convert .str file -> genind object -> dist object (Euclidean genetic distances),
# 2 text file of geographic coordinates -> geographic distance matrix
# 3 Mantel test of isolation by distance.
#   -histogram of simulations with a p-value
#   -heat map

## 1. ##
setwd('') #remember to put the .str file there!
library("ape")
library("genetics")    
library("seqinr")
library("ggplot2")
library("adegenet")

# Read in data.
obj_155 <- read.structure("erisor_shock.str", n.ind = 118, onerowperind= FALSE, n.loc = 331, col.lab = 1, col.pop = 0, col.others = NULL, row.marknames = 0)
indNames(obj_155) # I only get 247 because (24) low coverage individuals get filtered out.
ploidy(obj_155) # should return 2 since we gave it 2 alleles for each marker.

# Get Euclidean genetic distance matrix.
D_155 <- dist(tab(obj_155))

## 2 ##
# Geographic matrix function, obtained from https://eurekastatistics.com/calculating-a-distance-matrix-for-geographic-points-using-r/
ReplaceLowerOrUpperTriangle <- function(m, triangle.to.replace){
  # If triangle.to.replace="lower", replaces the lower triangle of a square matrix with its upper triangle.
  # If triangle.to.replace="upper", replaces the upper triangle of a square matrix with its lower triangle.
  
  if (nrow(m) != ncol(m)) stop("Supplied matrix must be square.")
  if      (tolower(triangle.to.replace) == "lower") tri <- lower.tri(m)
  else if (tolower(triangle.to.replace) == "upper") tri <- upper.tri(m)
  else stop("triangle.to.replace must be set to 'lower' or 'upper'.")
  m[tri] <- t(m)[tri]
  return(m)
}

GeoDistanceInMetresMatrix <- function(df.geopoints){
  # Returns a matrix (M) of distances between geographic points.
  # M[i,j] = M[j,i] = Distance between (df.geopoints$lat[i], df.geopoints$lon[i]) and
  # (df.geopoints$lat[j], df.geopoints$lon[j]).
  # The row and column names are given by df.geopoints$name.
  
  GeoDistanceInMetres <- function(g1, g2){
    # Returns a vector of distances. (But if g1$index > g2$index, returns zero.)
    # The 1st value in the returned vector is the distance between g1[[1]] and g2[[1]].
    # The 2nd value in the returned vector is the distance between g1[[2]] and g2[[2]]. Etc.
    # Each g1[[x]] or g2[[x]] must be a list with named elements "index", "lat" and "lon".
    # E.g. g1 <- list(list("index"=1, "lat"=12.1, "lon"=10.1), list("index"=3, "lat"=12.1, "lon"=13.2))
    DistM <- function(g1, g2){
      require("Imap")
      return(ifelse(g1$index > g2$index, 0, gdist(lat.1=g1$lat, lon.1=g1$lon, lat.2=g2$lat, lon.2=g2$lon, units="m")))
    }
    return(mapply(DistM, g1, g2))
  }
  
  n.geopoints <- nrow(df.geopoints)
  
  # The index column is used to ensure we only do calculations for the upper triangle of points
  df.geopoints$index <- 1:n.geopoints
  
  # Create a list of lists
  list.geopoints <- by(df.geopoints[,c("index", "lat", "lon")], 1:n.geopoints, function(x){return(list(x))})
  
  # Get a matrix of distances (in metres)
  mat.distances <- ReplaceLowerOrUpperTriangle(outer(list.geopoints, list.geopoints, GeoDistanceInMetres), "lower")
  
  # Set the row and column names
  rownames(mat.distances) <- df.geopoints$name
  colnames(mat.distances) <- df.geopoints$name
  
  return(mat.distances)
}
# Get geographic distance matrix.
setwd('')
geo <- read.csv("erisor_shock_coords.csv", header = TRUE, sep=",") 
geo
geo.df <- data.frame(geo)
geo.df
geo_notdist <- round(GeoDistanceInMetresMatrix(geo.df) / 1000)
Dgeo_155 <- as.dist(geo_notdist)

## 3 ##

# Confirm that the dimensions of matrices match.
str(Dgeo_155)
str(D_155)

# Mantel test for isolation by distance.
isobd <- mantel.randtest(D_155,Dgeo_155)
isobd

# Results will vary due to the simulation! But here's the results I got from this example:
"Monte-Carlo test
Call: mantel.randtest(m1 = D_155, m2 = Dgeo_155)

Observation: 0.2179472 

Based on 999 replicates
Simulated p-value: 0.001 
Alternative hypothesis: greater 

Std.Obs   Expectation      Variance 
5.2805989947 -0.0002957777  0.0017081039"

# Plot histrogram.
plot(isobd)
plot(Dgeo_155,D_155, xlab="Geographic Distance (m)", ylab="Genetic Distance")
abline(lm(D_155~Dgeo_155), col="red",lty=2)

# Plot heat map.
library(MASS)
dens <- kde2d(Dgeo_155,D_155, n=300)
myPal <- colorRampPalette(c("white","blue","gold", "orange", "red"))
plot(Dgeo_155, D_155, pch=20,cex=.5, xlab="Geographic Distance (m)", ylab="Genetic Distance")
image(dens, col=transp(myPal(300),.7), add=TRUE)
abline(lm(D_155~Dgeo_155))
title("Isolation by Distance")

