setwd('/Users/jimblotter/Desktop/Grad_School/Data_Analysis/erisor/QC/')
gomp_layout <- read.csv("gompert_layout_p1.csv", header = TRUE, stringsAsFactors = FALSE)
well_sample <- read.csv("well_sample.csv", header = TRUE, stringsAsFactors = FALSE)
str(gomp_layout)
well_sample

barcodes <- list()
for(i in well_sample$Well){
  for(j in gomp_layout$well){
    if (well_sample[i,1] == gomp_layout[j,1]){
      barcodes <- c(barcodes, gomp_layout$[j,2])
      j = j+1
      i = i+1
    }
  }
}
barcodes

  
  
  
  
  
  
  