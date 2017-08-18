setwd('/Users/jimblotter/Desktop/Grad_School/Data_Analysis/erisor/QC/')
library(fields)
library(RColorBrewer)
library(mapplots)
source("http://bioconductor.org/biocLite.R")
biocLite("LEA")
source("http://membres-timc.imag.fr/Olivier.Francois/Conversion.R")
source("http://membres-timc.imag.fr/Olivier.Francois/POPSutilities.R")

#download and convert my data:  (I want the .str file)
erisor.input.file =  "/Users/jimblotter/Desktop/Grad_School/Data_Analysis/erisor/QC/good_data.str"
struct2geno(file = input.file, TESS = TRUE, diploid = TRUE, FORMAT = 2,
            extra.row = 0, extra.col = 1, output = "secondary_contact.geno")  #may have more than one extra.col?
#run pop structure analysis that assumes K=3
library(LEA)
obj.snmf = snmf("secondary_contact.geno", K = 3, alpha = 100, project = "new")
qmatrix = Q(obj.snmf, K = 3)
barplot(t(qmatrix), col = c("orange","violet","lightgreen"), border = NA, space = 0,
        xlab = "Individuals", ylab = "Admixture coefficients")
