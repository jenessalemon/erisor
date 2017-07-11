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

##Which samples already have 1,750 loci? (exclude)

##Which samples have less than 400,000 reads? (exclude)

##Which samples have less than 70% expected loci? (exclude)
loci_per_read <- loci_in_assembly/reads_raw              #number of loci per read
av_loci <- mean(loci_per_read, na.rm=TRUE)               #0.0021
#lpr <- loci_per_read[!is.na(loci_per_read)]            #to get rid of the Na that is messing up my funciton
expected_loci <- reads_raw*av_loci
observed_vs_expected <- loci_in_assembly/expected_loci
ove <- observed_vs_expected[!is.na(observed_vs_expected)] #to get rid of the Na that is messing up my funciton

count <- function(nums){
  x=0
  for(i in nums){
    if(i >= .70){
      print(i)
     x=x+1
    }
  }
  return(x)
}
count(ove)                  #gives number of samples that passed the filter
