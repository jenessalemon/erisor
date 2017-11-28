## lemon.str is the structure file.
library(adegenet)
library(mmod)

## read in the data
obj <- read.structure('/Users/jimblotter/Desktop/Grad_School/Data_Analysis/erisor/QC/mantel/str_files/lemon.str', n.ind = 155, onerowperind= FALSE, n.loc = 211, col.lab = 1, col.pop = 2, col.others = NULL, row.marknames = 0)
## ! When it asks for additional info, just hit enter !

## use mmod to calculate the following:
diff_stats(obj)   
## I get $per.locus, Hs, Ht, Gst_est, Gprine_st D_het... then 
## NA for D_mean, and the following error:
## Warning message:
## In HsHt(g) : Need at least two population to calculate differentiation

pairwise_D(obj, linearized = FALSE, hsht_mean = "arithmetic") 
## I get the same error as before, but I also get the number: 0.01178946

D_Jost(obj)
## I get $per.locus information, $global.het... then NA for $global.harm_mean
## and the same warning message. 