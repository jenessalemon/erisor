# Regression of geographic and genetic distance. I found code to make the
# geographic distance matrix on the internet (below). The gentic distance matrix
# was generated in "str_analysis".
######################## Geographic Distance Matrix ####################################
# Got this code from: https://eurekastatistics.com/calculating-a-distance-matrix-for-geographic-points-using-r/
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

#Now to call on my own data:
setwd('/Users/jimblotter/Desktop/Grad_School/Data_Analysis/erisor/QC/') #remember to put the .str file there!
geo <- read.csv("geo.csv", header = TRUE, sep=",")
geo.df <- data.frame(geo)
geo_dist_matrix <- round(GeoDistanceInMetresMatrix(geo.df) / 1000)

##################### Genetic Distance Matrix ###################################
#Read in data
obj1 <- read.structure("p1.str", n.ind = 89, n.loc = 372, col.lab = 1, col.pop = 0, col.others = NULL, row.marknames = 0) #place cursor in console
# It will prompt for info:
#   genotypes = 244  (number of samples) This number can be found in the ipyrad _stats file, I had 266 but I threw out two samples with no data.
#   markers = 1886 (number of loci) Also find in ipyrad _stats file.
#   column with labels for genotypes = 1 (the names of the samples)
#   column with population factor = 0 (if we were specifying populations a priori but we aren't)
#   other optional columns - just hit enter, there are none in our data
#   row with the marker names = 0 (We don't have a marker names row)
#   genotypes coded by a single row = n (We have 2 lines per sample in the dataset.)
indNames(obj1) # I only get 247 because (24) low coverage individuals get filtered out.
ploidy(obj1) # should return 2 since we gave it 2 alleles for each marker.


#Neighbor joining euclidian distance tree
D <- dist(tab(obj1))               #super hard to read, create a distance matrix!
D                                  #0 means they are identical...and should be expected across the diagonal as each sample is identical to itself. 0s other than the diagonal better be replicates... Nas just mean that there are no loci in common between the two samples (we have no info ont their relatedness.)

#Convert to matrix for indexing.
M <- as(D, "matrix")               #Makes a normal matrix of the distance values.
M
str(M)

############################# Regression #####################################


