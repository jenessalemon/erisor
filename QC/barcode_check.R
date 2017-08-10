## The objective here is to solve the "columns vs wells" problem.
## The two input files are well_sample, which has two columns, "Well" and "Sample",
## And gomp_layout which also has two columns, "well" and "barcode." What I try to do is:

setwd('/Users/jimblotter/Desktop/Grad_School/Data_Analysis/erisor/QC/')
gomp_layout <- read.csv("gompert_layout_p1.csv", header = TRUE, stringsAsFactors = FALSE)
well_sample <- read.csv("well_sample.csv", header = TRUE, stringsAsFactors = FALSE)
gomp_layout #has well and barcode FROM JUST ONE PLATE
well_sample #has well and sample FROM JUST ONE PLATE

##AN ATTEMPT WITHOUT LOOPING##
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


  