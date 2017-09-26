setwd('/Users/jimblotter/Desktop/Grad_School/Data_Analysis/erisor/QC/mantel') #remember to put the .str file there!
coords <- read.csv("geo_for_mantel.csv", header = TRUE, sep=",")
p2_order <- read.csv("p2_order.csv", header = TRUE, sep=",")
p3_order <- read.csv("p3_order.csv", header = TRUE, sep=",")
p4_order <- read.csv("p4_order.csv", header = TRUE, sep=",")
p5_order <- read.csv("p5_order.csv", header = TRUE, sep=",")

#use barcodes check file to generate geographic coordinate files.
