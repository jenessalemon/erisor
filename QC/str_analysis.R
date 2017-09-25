## In this script I read in the relevant .str file, and convert to a genind object.
## I create a distance matrix and find the 5 closest individuals to any given
## individual (used for replicate analysis). I run Nei's genetic distance to get
## a different perspective. I generate a geographic distance matrix to perform a 
## regressiong against the genetic distance matrix and look for an association.

setwd('/Users/jimblotter/Desktop/Grad_School/Data_Analysis/erisor/QC/') #remember to put the .str file there!
#install.packages("genetics")
library("ape")
library("genetics")        #there is no package called ‘genetics’ = install.packages("genetics")
library("pegas")
library("seqinr")
library("ggplot2")
library("adegenet")
#library("poppr")
#library("ggmap")
#SUBTRACT ONE FROM NUMBER OF LOCI
#Read in data
obj1 <- read.structure("shock_51_noNA.str", n.ind = 34, onerowperind= FALSE, n.loc = 6935, col.lab = 1, col.pop = 0, col.others = NULL, row.marknames = 0) #place cursor in console
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
str(D)

#Convert to matrix for indexing.
M <- as(D, "matrix")               #Makes a normal matrix of the distance values.
M
str(M)

#Convert to dataframe:
df <- as.data.frame(M)
df

#Indexing Distance Matrix (M)
#All of the distance values for each individual are in each row, and each column, due to nature of the matrix.
row1 <- M["p_001s_02",]
row1                                #row 1 has the distance values between p_001s_01 and every other sample in this library.
str(row1)
row1["p_001s_03"]                   #index by sample name
row1[3]                             #by index number


#Now we identify the 5 closest relatives to a replicate sample:

#Find 5 closest relatives for just one sample (one row of the matrix):
find_relatives <-function(row){
    relatives <- list()                                    #initialize an empty list
    decreasing_index <- order(row, decreasing = FALSE)     #sort, smallest distance values first
    for (i in 1:5){
        relatives[i] <- decreasing_index[i]                  #first item in the list is the closest relative
    }
    return(relatives)
}
#call the function with one row of the matrix (one individual).
row1 <- M["p_026s_12",]                  #when looking for my replicates all I need to do is enter the sample name here and run to line 74.
output <- find_relatives(row1)
output

#Now I convert the output of find_relatives (index -> sample name).
names_list <- indNames(obj1)
names_list
index_to_samples <- function(find_relatives_output){
    samples <- list()
    for(i in find_relatives_output){
        #print(i)
        print(names_list[i])
        #samples <- names_list[i]
    }
    #return(samples)
}
index_to_samples(output)    #call the function with the output from "find_relatives"

#Now let's see if we can get the relatives from the entire list of replciate samples.
reps <- read.csv("list_of_replicates_first_library.csv", header = FALSE, sep = ",")
reps.list <- as.list(as.data.frame(t(reps)))
reps.list[3]
reps16 <- c("p_001s_02", "p_001s_04", "p_026s_12", "p_027s_06", "p_026s_14", "p_002s_06", "p_030s_11", "p_009s_13", "p_1027s_12", "p_013s_01", "p_012s_05", "p_015s_15", "p_016s_14", "p_019s_13", "p_013s_05", "p_024s_04")
reps16

for(i in reps.list){
    print(i)
}

for(row in M){
    if(row == index){
        R[i] <- M[i]
    }
    reps.list
    
    #So that I don't have to do this individually with 38 replicates:
    #Create a matrix that only retains line of replicate individuals, and apply to that entire matrix.
    
    replicate_matrix <- function(reps_list, matrix){       #create a matrix with only rows representing replicate samples
        R <- matrix(nrow = length(reps_list), ncol = length(indNames(obj1))    #Initialize an empty matrix
        print(str(R))
        for(row in row.names(matrix)){                       #M is the entire distance matrix
            if(row %in% reps_list == TRUE){                    #If the sample name is in the list of replciates
                R[row] <- M[row]                                 #add that row to our new matrix
            }
        }
        return(R)
    }
    rep_mat <- replicate_matrix(reps16, M)
    rep_mat
    str(rep_mat)
    ##
    length(indNames(obj1)
    if("p_019s_13" %in% reps16 == TRUE){
        print("hi")
    }
    #works
    for(row in row.names(M)){                       #M is the entire distance matrix
        if(row %in% reps16 == TRUE){                    #If the sample name is in the list of replciates
            print("hello")                                 #add that row to our new matrix
        }
    }
    #works
    
    
    find_my_replicates <- function(rep_matrix){              #take the output from replciate_matrix (line 97)
        for(row in rep_matrix){
            print(row)
            index <- find_relatives()
            index_to_samples(indices)
        }
    }
    
    indices <- find_relatives(reps.list)
    index_to_samples(indices)
    
#Nei's Distance
nei <- nei.dist(obj1)
aboot(nei)               #dendrogram using Nei's distance. Passing in indices rather than sample names
#need to call sample names not indexes.

############## Genetic Distance Matrix - Population Level (INCOMPLETE) ############################
M            #this is the genetic distance matrix at the individual level.
             #I would like to get a genetic distance matrix at the population level.
# First, I converted to a dataframe.
df <- as.data.frame(M)
df
# I will need to make pairwise comparisons across all indiviudals in the matrix.
# I would like to write a function that will do parwise comparisons for pop1 and pop2, pop1 and pop 3, etc.

# Let's start small with a dummy dataframe:
n = c(2, 3, 5) 
s = c(3, 4, 1) 
b = c(1, 2, 6) 
df2 = data.frame(n, s, b)

############# Isolation by Distance #############################
#1-genetic distance matrix (JUST USE D)
#pop <- c(1,1,2.5,2.5,2.5,2,2,2,2,2,2,2,2,2,2,4.5,5,5,6,6,6,9,9,11,11,11,12,12,13,13,15,15,16,16,16,16,16,16,17,17,17,26,26,26,27,27,27,29,31,1027,1027)
#popyo <- as.data.frame(pop)    #converting to a dataframe because apparently that is all strata() will take
#strata(obj1, combine = TRUE, value = popyo, name = "pop") #can't get the name to stick
#toto <- genind2genpop(obj1)          #but I ended up using D...
#Dgen <- dist.genpop(toto,method=2)

#2-geographic distance matrix
###Need to read in the coordinates...from a file in the same order and with the same samples.
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

setwd('/Users/jimblotter/Desktop/Grad_School/Data_Analysis/erisor/QC/') #remember to put the .str file there!
geo <- read.csv("geo_for_shock51.csv", header = TRUE, sep=",")
geo.df <- data.frame(geo)
geo_notdist <- round(GeoDistanceInMetresMatrix(geo.df) / 1000)
Dgeo <- as.dist(geo_notdist) #convert to dist object

# Ok, now I can perform the mantel test
isobd <- mantel.randtest(D,Dgeo)
isobd
# Result for shock_51:
"Monte-Carlo test
Call: mantel.randtest(m1 = D, m2 = Dgeo)

Observation: 0.5836126 

Based on 999 replicates
Simulated p-value: 0.001 
Alternative hypothesis: greater 

Std.Obs  Expectation     Variance 
12.606791936 -0.001904555  0.002157099"
# plotting
plot(isobd)
plot(Dgeo,D)
abline(lm(D~Dgeo), col="red",lty=2)
# heat map
library(MASS)
dens <- kde2d(Dgeo,D, n=300)
myPal <- colorRampPalette(c("white","blue","gold", "orange", "red"))
plot(Dgeo, D, pch=20,cex=.5)
image(dens, col=transp(myPal(300),.7), add=TRUE)
abline(lm(D~Dgeo))
title("Isolation by distance plot")

