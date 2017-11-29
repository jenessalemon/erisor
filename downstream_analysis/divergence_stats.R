## This R code calculates Jost's D using the package "mmod" (David Winter).
## Load required packages
library("ape")
library("genetics")    
library("seqinr")
library("ggplot2")
library("adegenet")
library("mmod")
library("swfscMisc")

## Convert structure file to genind
obj_155_d <- read.structure('/Users/jimblotter/Desktop/Grad_School/Data_Analysis/erisor/downstream_analysis/mmod/lemon.str', n.ind = 155, onerowperind= FALSE, n.loc = 211, col.lab = 1, col.pop = 2, col.others = NULL, row.marknames = 0)
obj_155_d
obj_155_d@pop

## Remove a locus with no data
genotyped <- obj_155_d[,locFac(obj_155_d) != "L011"]

## Pairwise D (This is just the same as D because I only have one popualtion.)
pairwise_D(genotyped, linearized = FALSE, hsht_mean = "arithmetic") #0.01178946

## Obtain useful statistics
diff_stats(genotyped)

## Manual calculation of D:
## [(Ht - Hs)/(1-Hs)]*[(n/(n-1))]             #Eq. 11 from Jost 2008. [(n/(n-1))] simplifies to 2, because I just have 2 populations.
((0.05167627-0.04605300)/(1-0.04605300))*2    #0.01178948 

## A final (third) way to calculate global D
mean(diff_stats(genotyped)$per.locus[,5]) #0.01588839
