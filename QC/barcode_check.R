## The objective here is to solve the "columns vs wells" problem. Barcodes are pulled
## from a different dataframe, and assigned to the appropriate spot in the dataframe
## containing samples and wells. Note that this must be done one plate at a time!
## The two input files are well_sample, which has two columns, "Well" and "Sample",
## And gomp_layout which also has two columns, "well" and "barcode."

setwd('/Users/jimblotter/Desktop/Grad_School/Data_Analysis/erisor/QC/')

diff_barcode_check <- function(gomp_layouts, well_samples){
  gomp_layout <- read.csv(gomp_layouts, header = TRUE, stringsAsFactors = FALSE)
  well_sample <- read.csv(well_samples, header = TRUE, stringsAsFactors = FALSE)
  A1A2A3 <- c(gomp_layout$well)
  for(i in A1A2A3){
    samp <- well_sample[which(well_sample$Well==i),2]        #bcode is a holder variable, which gets, from gompert_layout, the variable in the second column from the row where the well = value at ith iteration in the list.
    gomp_layout[which(gomp_layout$well==i),2] <- samp        #assign that variable (a barcode), to the (previously empty) cell in Barcode column of well_sample, on the row where that same variable is found.
  }
  #gomp_layout$well <- NULL                                   #drop the well column
  return(gomp_layout)
}

#Now, let's run it on all three
product1 <- diff_barcode_check("gompert_layout_p1.csv", "well_samp_p1.csv")                  #save output of the function to a variable
product1
write_tsv(product1, path = '/Users/jimblotter/Desktop/Grad_School/Data_Analysis/erisor/QC/Ripy_barcodes1')                              #tab delineated

product2 <- diff_barcode_check("gompert_layout_p2.csv", "well_samp_p2.csv")                  #save output of the function to a variable
product2
write_tsv(product2, path = '/Users/jimblotter/Desktop/Grad_School/Data_Analysis/erisor/QC/Ripy_barcodes2')                              #tab delineated

product3 <- diff_barcode_check("gompert_layout_p3.csv", "well_samp_p3.csv")                  #save output of the function to a variable
product3
write_tsv(product3, path = '/Users/jimblotter/Desktop/Grad_School/Data_Analysis/erisor/QC/Ripy_barcodes3')                              #tab delineated

#Let's do all three:

