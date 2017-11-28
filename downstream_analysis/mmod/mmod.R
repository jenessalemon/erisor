## lemon.str is the structure file.
library(adegenet)
library(mmod)

## read in the data
obj <- read.structure('/path/lemon.str', n.ind = 155, onerowperind= FALSE, n.loc = 211, col.lab = 1, col.pop = 2, col.others = NULL, row.marknames = 0)
## ! When it asks for additional info, just hit enter !

## use mmod to calculate the following:
diff_stats(obj)   
## I get $per.locus, Hs, Ht, Gst_est, Gprine_st D_het... then 
## D_mean = NA, and the following error:
## Warning message:
## In HsHt(g) : Need at least two population to calculate differentiation

pairwise_D(obj, linearized = FALSE, hsht_mean = "arithmetic") 
## I get the same error as before, but I also get the number: 0.01178946

D_Jost(obj)
## I get $per.locus, $global.het, then $global.harm_mean = NA
## and the same warning message.

##Attempt at manual calculation:
## [(Ht - Hs)/(1-Hs)]*[(n/(n-1))]             #Eq. 11 from Jost 2008. [(n/(n-1))] simplifies to 2, because I just have 2 populations.
((0.05167627-0.04605300)/(1-0.04605300))*2    #0.01178948 
## So, I have what I need. ASSUMING it's correct. I'd like the get the
## code working to verify.

## To check whether or not missing data was the problem, I created a subet
## only including loci that I could see have data. I got the same results as
## above.
sub <- obj[loc = 22:100]        #random subset that I can see have data
diff_stats(sub)  
pairwise_D(sub, linearized = FALSE, hsht_mean = "arithmetic")
D_Jost(sub)

