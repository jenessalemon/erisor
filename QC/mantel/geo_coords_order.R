setwd('/Users/jimblotter/Desktop/Grad_School/Data_Analysis/erisor/QC/mantel') #remember to put the .str file there!
coords <- read.csv("geo_for_mantel.csv", header = TRUE, sep=",")
p2_order <- read.csv("p2_order.csv", header = TRUE, sep=",")
p3_order <- read.csv("p3_order.csv", header = TRUE, sep=",")
p4_order <- read.csv("p4_order.csv", header = TRUE, sep=",")
p5_order <- read.csv("p5_order.csv", header = TRUE, sep=",")

#Generate geographic coordinate files.
order <- p5_order$sample
order
coords[,2]
geo_samps <- function(coords_file, order_file){
  order <- order_file$sample
  for(i in order){
    long <- coords_file[which(coords_file$sample==i),2]       #samp is a holder variable, which gets, from gompert_layout, the variable in the second column from the row where the well = value at ith iteration in the list.
    print(long)
    order_file[i,2] <- long        #assign that variable (a barcode), to the (previously empty) cell in Barcode column of well_sample, on the row where that same variable is found.
  }
  #order_file$sample <- NULL                                   #drop the sample column
  geo_complete <- data.frame(order_file)
  return(geo_complete)
}
#which(order_file$sample==i)
#Now, let's run it on all three
hey <- geo_samps(coords, p5_order)    
hey
