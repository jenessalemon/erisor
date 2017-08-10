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
for(i in 1:10){
  for(j in 11:20){
    if (i == 13){
      hi <- c(hi, j) #if I put i, it prints 1:10, if I put j, it prints 13 ten times.
      j = j+1
      i = i+1
    }
  }
}
hi
#
barcodes <- list()
for(i in well_sample$Well){
  for(j in gomp_layout$well){
    if (i == j){
      barcodes <- c(barcodes, gomp_layout[j,2]) #if I put i, it prints 1:10, if I put j, it prints 13 ten times.
    }
    j = as.numeric(j)+1
  }
  i = as.numeric(i)+1
}
barcodes

  
gomp_layout[1,2] #works, so why does it not work when j is there in place of 1?
  