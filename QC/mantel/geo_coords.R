setwd('/Users/jimblotter/Desktop/Grad_School/Data_Analysis/erisor/QC/mantel/str_analysis_geoinput') #remember to put the .str file there!
order <- read.csv("geo_shock177.csv", header = TRUE, stringsAsFactors = FALSE)
order$Sample
coords <- read.csv("geo_for_mantel_names.csv", header = TRUE, stringsAsFactors = FALSE)
coords

coords_matrix <- function(order, coords){
  orders <- c(order$Sample)
  for(i in orders){
    latitude <- coords[which(coords$sample==i),3]        #bcode is a holder variable, which gets, from gompert_layout, the variable in the second column from the row where the well = value at ith iteration in the list.
    order[which(order$Sample==i),3] <- latitude        #assign that variable (a barcode), to the (previously empty) cell in Barcode column of well_sample, on the row where that same variable is found.
  }
  for(i in orders){
    longitude <- coords[which(coords$sample==i),2]        #bcode is a holder variable, which gets, from gompert_layout, the variable in the second column from the row where the well = value at ith iteration in the list.
    order[which(order$Sample==i),2] <- longitude        #assign that variable (a barcode), to the (previously empty) cell in Barcode column of well_sample, on the row where that same variable is found.
  }
  order$Sample <- NULL                                   #drop the well column
  coords_mat <- data.frame(order)
  #out <- barcode_mat[]
  return(coords_mat)
}
heythere <- coords_matrix(order, coords)
heythere
write_tsv(heythere, path = '/Users/jimblotter/Desktop/Grad_School/Data_Analysis/erisor/QC/mantel/str_analysis_geoinput/177shock')                              #tab delineated
