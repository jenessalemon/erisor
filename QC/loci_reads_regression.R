setwd('/Users/jimblotter/Desktop/Grad_School/Data_Analysis/erisor/QC/')
stats <- read.table("erisor_reps_stats.txt", header = TRUE, fill = TRUE)
reads_raw <- stats[,2]
loci_in_assembly <- stats[,9]

reads_and_loci <- subset(stats, select = c("reads_raw", "loci_in_assembly"))
summary(reads_and_loci)
plot(reads_and_loci)

fit <- lm(reads_raw ~ loci_in_assembly, data=stats)   #linear regression
summary(fit)
plot(fit)

loci_per_read <- loci_in_assembly/reads_raw
mean(loci_per_read, na.rm=TRUE)               #0.0021
lpr <- loci_per_read[!is.na(loci_per_read)]   #to get rid of the Na that is messing up my funciton
count <- function(nums){
  x=0
  for(i in nums){
    if(i >= .0021){
     x=x+1
    }
  }
  return(x)
}
count(lpr)                                    #gives number of samples that passed the filter
