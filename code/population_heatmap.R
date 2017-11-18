## This code generates a heat map showing distances between populations.
## The input is a structure file (.str), and a .csv file of the populations
## each sample belongs to, in the order which they appear in the the .str file.

# Load required software.
setwd("")
library("ape")
library("genetics")    
library("seqinr")
library("ggplot2")
library("adegenet")

# Read in the data.
obj_155_pop <- read.structure("erisor.str", n.ind = 155, onerowperind= FALSE, n.loc = 211, col.lab = 1, col.pop = 0, col.others = NULL, row.marknames = 0)
indNames(obj_155_pop) 
ploidy(obj_155_pop)

# Add Population data.
popFac <- read.csv("popFac.csv", header = FALSE)
popFac <- factor(unlist(popFac))
pop(obj_155_pop) <- popFac
pop(obj_155_pop)

# Convert to genpop.
genind_155 <- genind2genpop(obj_155_pop, pop = popFac, quiet = FALSE, process.other = FALSE,
              other.action = mean)

# Generate the distance matrix.
pop_dist <- dist.genpop(genind_155, method = 1, diag = TRUE, upper = TRUE)

# Draw the heatmap.
heatmap(as.matrix(pop_dist), symm=TRUE, keep.dendro = FALSE)
