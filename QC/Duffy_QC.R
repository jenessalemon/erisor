setwd('/Users/jimblotter/Desktop/Grad_School/Data_Analysis/QC/') #remember to put the .str file there!
#install.packages("genetics")
library("ape")
library("genetics")        #there is no package called ‘genetics’ = install.packages("genetics")
library("pegas")
library("seqinr")
library("ggplot2")
library("adegenet")

#Read in data
obj1 <- read.structure("erisor_reps_9sout.str", n.ind = 262, n.loc = 1886, col.lab = 1, col.pop = 0, col.others = NULL, row.marknames = 0) #place cursor in console
# It will prompt for info:
#   genotypes = 244  (number of samples) This number can be found in the ipyrad _stats file, I had 266 but I threw out two samples with no data. 
#   markers = 1886 (number of loci) Also find in ipyrad _stats file.
#   column with labels for genotypes = 1 (the names of the samples)
#   column with population factor = 0 (if we were specifying populations a priori but we aren't)
#   other optional columns - just hit enter, there are none in our data
#   row with the marker names = 0 (We don't have a marker names row)
#   genotypes coded by a single row = n (We have 2 lines per sample in the dataset.)
indNames(obj1) # I only get 238 because (24) low coverage individuals get filtered out.
ploidy(obj1) # should return 2 since we gave it 2 alleles for each marker.


#Neighbor joining euclidian distance tree
D <- dist(tab(obj1))               #super hard to read, create a distance matrix! 

#What does D look like?
D                                  #0 means they are identical...and should be expected across the diagonal as each sample is identical to itself. 0s other than the diagonal better be replicates... Nas just mean that there are no loci in common between the two samples (we have no info ont their relatedness.)
class(D)                           #dist
attributes(D)
str(D)
D[100]                             #can index by row number
D["p_019s_11"]                     #can index by sample name

#Convert to matrix for indexing.
M <- as(D, "matrix")               #Makes a normal matrix of the distance values.
M                                  
str(M)

#Indexing Distance Matrix (M)
#All of the distance values for each individual are in each row, and each column, due to nature of the matrix.
row1 <- M["p_001s_01",]
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
row1 <- M["p_023s_13",]                  #when looking for my replicates all I need to do is enter the sample name here and run to line 74.
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
  R <- matrix(nrow = length(reps_list), ncol = 238)    #Initialize an empty matrix
  for(row in row.names(matrix)){                       #M is the entire distance matrix
    if(row %in% reps_list == TRUE){                    #If the sample name is in the list of replciates
      R[row] <- M[row]                                 #add that row to our new matrix
    }
  }
  return(R)
}
R <- replicate_matrix(reps.list, M)
str(R)
  
find_my_replicates <- function(rep_matrix){              #take the output from replciate_matrix (line 97)
  for(row in rep_matrix){
    print(row)
    index <- find_relatives()
    index_to_samples(indices)
  }
}

indices <- find_relatives(reps.list)
index_to_samples(indices)

#Making a Tree
tre <- njs(D)
jpeg(height=960, width=960)
par(xpd=TRUE, mar=c(0,0,0,0))
plot(tre, type="phylogram", edge.w=1)
dev.off()