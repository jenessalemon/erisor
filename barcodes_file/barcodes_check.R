setwd('/Users/jimblotter/Desktop/Grad_School/Data_Analysis/erisor/barcodes_file/')
install.packages(compare)
library(compare)
before = read.csv("before.csv", header = FALSE)
after = read.csv("after.csv", header = FALSE)
before
after

for(i in before[,1]){
  if(after)
}


barcodes_after[,1] == barcodes_before[,1]
