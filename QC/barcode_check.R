## The objective here is to solve the "columns vs wells" problem. Barcodes are pulled
## from a different dataframe, and assigned to the appropriate spot in the dataframe
## containing samples and wells. Note that this must be done one plate at a time!
## The two input files are well_sample, which has two columns, "Well" and "Sample",
## And gomp_layout which also has two columns, "well" and "barcode."

setwd('/Users/jimblotter/Desktop/Grad_School/Data_Analysis/erisor/QC/')
gomp_layout <- read.csv("gompert_layout_p1.csv", header = TRUE, stringsAsFactors = FALSE)
well_sample <- read.csv("well_sample.csv", header = TRUE, stringsAsFactors = FALSE)
gomp_layout #has well and barcode FROM JUST ONE PLATE
well_sample #has well and sample FROM JUST ONE PLATE
my_well_list <- c(gomp_layout$well) #a list of wells from the gompert layout.

barcode_check <- function(well_list){
  gomp_layout <- read.csv("gompert_layout_p1.csv", header = TRUE, stringsAsFactors = FALSE)
  well_sample <- read.csv("well_sample.csv", header = TRUE, stringsAsFactors = FALSE)
  for(i in well_list){
    bcode <- gomp_layout[which(gomp_layout$well==i),2]             #bcode is a holder variable, which gets, from gompert_layout, the variable in the second column from the row where the well = value at ith iteration in the list.
    well_sample[which(well_sample$Well==i),3] <- bcode       #assign that variable (a barcode), to the (previously empty) cell in Barcode column of well_sample, on the row where that same variable is found.
  }
  return(well_sample)
}
barcode_check(my_well_list)

#Here I need to erase the well column, and export. I also need to run for the other 2 plates.
#Note that the way my brain decided to think of this, I can't run diff to find differences.
#If I wanted to do that they would need to be in the gompert order, not my origional order.
#This should be somewhat easy to change in necessary. I would just pull samples from well_samples
#instead of pulling barcodes from gomp_layout. In fact, let's try:

diff_barcode_check <- function(well_list){
  gomp_layout <- read.csv("gompert_layout_p1.csv", header = TRUE, stringsAsFactors = FALSE)
  well_sample <- read.csv("well_sample.csv", header = TRUE, stringsAsFactors = FALSE)
  for(i in well_list){
    samp <- well_sample[which(well_sample$Well==i),2]        #bcode is a holder variable, which gets, from gompert_layout, the variable in the second column from the row where the well = value at ith iteration in the list.
    gomp_layout[which(gomp_layout$well==i),2] <- samp        #assign that variable (a barcode), to the (previously empty) cell in Barcode column of well_sample, on the row where that same variable is found.
  }
  return(gomp_layout)
}
product <- diff_barcode_check(my_well_list)                  #save output of the function to a variable
product$well <- NULL                                         #drop the well column
path <- '/Users/jimblotter/Desktop/Grad_School/Data_Analysis/erisor/QC/diff_barcodes' #probably not needed
write_tsv(product, path = path)                              #tab delineated

######################### Practice that led to the function below ################################
which(gomp_layout$Well=='B1')
#works
well_sample <- read.csv("well_sample.csv", header = TRUE, stringsAsFactors = FALSE)
#re-read it in to start fresh
bcode <- gomp_layout[which(gomp_layout$well=='D1'),2]             #bcode is a holder variable, which gets, from gompert_layout, the variable in the second column from the row where the well = "D1".
well_sample$Barcode[which(well_sample$Well=='D1')] <- bcode       #assign that variable (a barcode), to the (previously empty) cell in well_sample$Barcode, on the row where that same variable is found.
well_sample                                                       #check to see if it worked
## Works for one. Can I get it to go A1 through H12? Can't see a way to do this wihout looping...
#First, I will wrap into a function:
barcode_check <- function(well){
  well_sample <- read.csv("well_sample.csv", header = TRUE, stringsAsFactors = FALSE)
  bcode <- gomp_layout[which(gomp_layout$well==well),2]             #bcode is a holder variable, which gets, from gompert_layout, the variable in the second column from the row where the well = "D1".
  well_sample$Barcode[which(well_sample$Well==well)] <- bcode       #assign that variable (a barcode), to the (previously empty) cell in well_sample$Barcode, on the row where that same variable is found.
  well_sample 
  return(well_sample)
}
barcode_check("C1")

bcodes_list <- c(gomp_layout$well)
bcodes_list
#now, back to looping
barcode_check <- function(bcodes_list){
  gomp_layout <- read.csv("gompert_layout_p1.csv", header = TRUE, stringsAsFactors = FALSE)
  well_sample <- read.csv("well_sample.csv", header = TRUE, stringsAsFactors = FALSE)
  for(i in bcodes_list){
    bcode <- gomp_layout[which(gomp_layout$well==i),2]             #bcode is a holder variable, which gets, from gompert_layout, the variable in the second column from the row where the well = "D1".
    well_sample[which(well_sample$Well==i),3] <- bcode       #assign that variable (a barcode), to the (previously empty) cell in well_sample$Barcode, on the row where that same variable is found.
  }
  return(well_sample)
}
barcode_check(bcodes_list)


  