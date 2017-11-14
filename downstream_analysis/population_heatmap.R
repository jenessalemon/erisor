setwd('/Users/jimblotter/Desktop/Grad_School/Data_Analysis/erisor/downstream_analysis') #remember to put the .str file there!
library("ape")
library("genetics")    
library("seqinr")
library("ggplot2")
library("adegenet")
#SUBTRACT ONE FROM NUMBER OF LOCI
#read in data
obj_155_pop <- read.structure('/Users/jimblotter/Desktop/Grad_School/Data_Analysis/erisor/QC/mantel/str_files/final163.str', n.ind = 155, onerowperind= FALSE, n.loc = 211, col.lab = 1, col.pop = 0, col.others = NULL, pop = popFac, row.marknames = 0)
indNames(obj_155_pop) 
ploidy(obj_155_pop) # should return 2 since we gave it 2 alleles for each marker.
pop(obj_155_pop)

#add popfac
popFac <- read.csv("pop.csv", header = FALSE)
popFac <- factor(unlist(popFac))
pop(obj_155_pop) <- popFac
pop(obj_155_pop)


genind_155 <- genind2genpop(obj_155_pop, pop = popFac, quiet = FALSE, process.other = FALSE,
              other.action = mean)

pop_dist <- dist.genpop(genind_155, method = 1, diag = TRUE, upper = TRUE)
pop_mat <- as.matrix(pop_dist)

heatmap(as.matrix(pop_dist), symm=TRUE, keep.dendro = FALSE)

########### IN GGPLOT2? #############
(p <- ggplot(nba.m, aes(variable, Name)) + geom_tile(aes(fill = rescale),
                                                       +     colour = "white") + scale_fill_gradient(low = "white",
                                                                                                     +     high = "steelblue"))






