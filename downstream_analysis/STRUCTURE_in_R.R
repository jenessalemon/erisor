## This script allows the user to run STRUCTURE, which evaluates population structure
## using a Bayesian framework. The publication can be found here: http://membres-timc.imag.fr/Olivier.Francois/tutoRstructure.pdf

setwd('/Users/jimblotter/Desktop/Grad_School/Data_Analysis/erisor/downstream_analysis/')
library(fields)
library(RColorBrewer)
library(mapplots)
source("http://bioconductor.org/biocLite.R")
biocLite("LEA")
source("http://membres-timc.imag.fr/Olivier.Francois/Conversion.R")
source("http://membres-timc.imag.fr/Olivier.Francois/POPSutilities.R")

#run pop structure analysis that assumes K=3
library(LEA)
obj.snmf = snmf("p1_good.geno", K = 2, alpha = 100, project = "new")
qmatrix = Q(obj.snmf, K = 2)
barplot(t(qmatrix), col = c("orange","violet","lightgreen"), border = NA, space = 0,
        xlab = "Individuals", ylab = "Admixture coefficients", main = "erisor plate 1")
#need to label each individual...ugh

#What to choose for K?
obj.snmf = snmf("p1_good.geno", K = 1:8,  ploidy = 2, entropy = T,
                alpha = 100, project = "new")
plot(obj.snmf, col = "blue4", cex = 1.4, pch = 19)