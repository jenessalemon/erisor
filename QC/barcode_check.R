setwd('/Users/jimblotter/Desktop/Grad_School/Data_Analysis/erisor/QC/')
gomp_layout <- read.csv("gompert_layout_p1.csv", header = TRUE, stringsAsFactors = FALSE)
well_sample <- read.csv("well_sample.csv", header = TRUE, stringsAsFactors = FALSE)
str(gomp_layout)
well_sample

'''df <- data.frame(matrix(ncol = 2, nrow = 283))
#for each sample in well_sample, find the well, find it in the other dataframe, append the barcode from the other data frame
for(i in well_sample){
  print(well_sample[i,1])
}



for(i in well_sample$Well){
  if (well_sample[i, 1] == gomp_layout$well){
    well_sample$Barcode <- gomp_layout[i,2]
  }
}
well_sample'''

for(i in well_sample$Well){
  for(j in gomp_layout$well){
    if (i == j){
    well_sample$Barcode <- gomp_layout[i,2]
    }
  }
}
well_sample

  
  
  
  
  
  
  