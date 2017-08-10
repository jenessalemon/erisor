## The objective here is to solve the "columns vs wells" problem.
## The two input files are well_sample, which has two columns, "Well" and "Sample",
## And gomp_layout which also has two columns, "well" and "barcode." What I try to do is:
## Loop through the list of wells in well_sample. For each of the wells, I want to
## compare to the wells in gomp_layout. When there is a match, that barcode needs to be
## added to an emtpy list. Then I will add the list of barcodes in (now in appropriate order)
## to well_sample as a new column, which will give me three columns, "Well","Sample", "Barcode."

setwd('/Users/jimblotter/Desktop/Grad_School/Data_Analysis/erisor/QC/')
gomp_layout <- read.csv("gompert_layout_p1.csv", header = TRUE, stringsAsFactors = FALSE)
well_sample <- read.csv("well_sample.csv", header = TRUE, stringsAsFactors = FALSE)
gomp_layout #has well and barcode FROM JUST ONE PLATE
well_sample #has well and sample FROM JUST ONE PLATE

barcodes <- list()                                 #initialize empty list
for(i in well_sample$Well){                        #for each well in well_sample
  for(j in gomp_layout$well){                      #search each well gomp_layout
    if (well_sample[i,1] == gomp_layout[j,1]){     #when you find the well in the barcodes file
      barcodes <- c(barcodes, gomp_layout[j,2])    #append the barcode to a list (then I plan to add the list to the dataframe as a new column)
    }
    j = j+1
  }
  i = i+1
}
barcodes

well_sample$Barcode <- barcodes                    #this will append the list of barcodes
  
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
#works

#what I can't do, is access the value behind the ith iteration 
  