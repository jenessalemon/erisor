setwd('/Users/jimblotter/Desktop/Grad_School/Data_Analysis/erisor/QC/')
gomp_layout <- read.csv("gompert_layout_p1.csv", header = TRUE, stringsAsFactors = FALSE)
well_sample <- read.csv("well_sample.csv", header = TRUE, stringsAsFactors = FALSE)
gomp_layout
well_sample

barcodes <- list()
for(i in well_sample$Well){
  for(j in gomp_layout$well){
    if (well_sample[i,1] == gomp_layout[j,1]){
      barcodes <- c(barcodes, gomp_layout[j,2])
    }
    j = j+1
  }
  i = i+1
}
barcodes

  
#Small tests
well_sample[1,1] == gomp_layout[1,1]
#=TRUE
 
hey <- c(1, 2, 3, 4) 
hey <- c(hey, gomp_layout[1,1])
hey
#works


hi <- list()
for(i in 1:11){
  for(j in 11:11){
    if (i == j){
      hi <- c(hi, gomp_layout[j,2]) 
      j = j+1
      i = i+1
    }
  }
}
hi
#
barcodes <- list()
well_sample_well <- well_sample$Well
gomp_layout_well <- gomp_layout$well
for(i in well_sample_well){
  for(j in gomp_layout_well){
    if ( == j){
      barcodes <- c(barcodes, gomp_layout[j,2])
    }
    j = j+1
  }
  i = i+1
}
barcodes

#what I can't do, is access the value behind the ith iteration 

  
gomp_layout[1,2] #works, so why does it not work when j is there in place of 1?
  