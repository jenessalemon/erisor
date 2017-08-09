setwd('/Users/jimblotter/Desktop/Grad_School/Data_Analysis/erisor/QC/')
gomp_layout <- read.csv("gompert_layout_p1.csv", header = TRUE)
well_sample <- read.csv("well_sample_p1.csv", header = TRUE)
str(gomp_layout)
well_sample

#df <- data.frame(matrix(ncol = 2, nrow = 283))
#for each sample in well_sample, find the well, find it in the other dataframe, append the barcode from the other data frame
for(i in well_sample[,2]){
  print(i)
}
