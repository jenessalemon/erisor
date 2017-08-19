install.packages(c("fields","RColorBrewer","mapplots"))
source("http://bioconductor.org/biocLite.R")
biocLite("LEA")

source("http://membres-timc.imag.fr/Olivier.Francois/Conversion.R")
source("http://membres-timc.imag.fr/Olivier.Francois/POPSutilities.R")

#download and convert example data:
input.file =  "http://membres-timc.imag.fr/Olivier.Francois/secondary_contact.str"
struct2geno(file = input.file, TESS = TRUE, diploid = TRUE, FORMAT = 2,
            extra.row = 0, extra.col = 0, output = "secondary_contact.geno")
#might not need to convert .str to .geno...ipyrad provides .geno.

#run pop structure analysis that assumes K=3
library(LEA)
obj.snmf = snmf("secondary_contact.geno", K = 3, alpha = 100, project = "new")
qmatrix = Q(obj.snmf, K = 3)
barplot(t(qmatrix), col = c("orange","violet","lightgreen"), border = NA, space = 0,
        xlab = "Individuals", ylab = "Admixture coefficients")

#Read in coordinates, create sample identifiers (this will be hard to replicate, my data are not as clean)
coord = read.table("coordinates.coord")
pop = rep(1:60, each = 10)

#Obtain pop estimates of admixture cooeficients:
K=3
Npop = length(unique(pop))
qpop = matrix(NA, ncol = K, nrow = Npop) 
coord.pop = matrix(NA, ncol = 2, nrow = Npop) 
for (i in unique(pop)){
  qpop[i,] = apply(qmatrix[pop == i,], 2, mean)
  coord.pop[i,] = apply(coord[pop == i,], 2, mean)}

#Geographic mapping of admixture proportions:
library(mapplots)
plot(coord, xlab = "Longitude", ylab = "Latitude", type = "n")
map(add = T, col = "grey90", fill = TRUE)
for (i in 1:Npop){
  add.pie(z = qpop[i,], x = coord.pop[i,1], y = coord.pop[i,2], labels = "",
  col = c("orange","violet","lightgreen"))}

#Loading a vector of population tables for my data?
pop = scan("mypop.txt"

#Choosing the number of clusters:
#Perform runs for 8 values of K, choose the value of K for which the cross entropy curve exhibits a plateau
obj.snmf = snmf("secondary_contact.geno", K = 1:8,  ploidy = 2, entropy = T,
                alpha = 100, project = "new")
plot(obj.snmf, col = "blue4", cex = 1.4, pch = 19)