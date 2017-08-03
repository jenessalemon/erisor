#order species and populations together, use ggplot for heat map
##First, create a small .str file for practice:
setwd('/Users/jimblotter/Desktop/Grad_School/Data_Analysis/QC/') #remember to put the .str file there!
#install.packages("genetics")
library("ape")
library("genetics")        
library("pegas")
library("seqinr")
library("ggplot2")
library("adegenet")
obj1 <- read.structure("erisor_reps_9sout.str")
indNames(obj1) 
ploidy(obj1)
D <- dist(tab(obj1))
D
M <- as(D, "matrix")               #Makes a normal matrix of the distance values.
M                                  
str(M)
#I don't really find this helpful...but I'm willing to try it because I don't know what to do!

plot(M)
