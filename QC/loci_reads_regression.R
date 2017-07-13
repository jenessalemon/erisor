setwd('/Users/jimblotter/Desktop/Grad_School/Data_Analysis/erisor/QC/')
stats <- read.table("erisor_reps_stats-with_spp.txt", header = TRUE, fill = TRUE)
str(stats)
reads_raw <- stats[,2]
loci_in_assembly <- stats[,9]

reads_and_loci <- subset(stats, select = c("reads_raw", "loci_in_assembly"))
summary(reads_and_loci)
plot(reads_and_loci)

fit <- lm(reads_raw ~ loci_in_assembly, data=stats)   #linear regression
summary(fit)
plot(fit)

##How many samples already have 1,750 loci? (exclude 52 samples)
lia <- loci_in_assembly[!is.na(loci_in_assembly)]
count_lia <- function(nums){
  x=0
  for(i in nums){
    if(i >= 1750){
      #print(i)
      x=x+1
    }
  }
  print("Loci with more than 1750 loci:")
  return(x)
}
count_lia(lia) 

##How many samples have less than 400,000 raw reads? (exclude 109 samples)
rr <- reads_raw[!is.na(reads_raw)]
count_rr <- function(nums){
  x=0
  for(i in nums){
    if(i < 400000){
      #print(i)
      x=x+1
    }
  }
  print("Samples with less than 400,000 raw reads:")
  return(x)
}
count_rr(rr)

##How many samples have less than 70% expected loci?
loci_per_read <- loci_in_assembly/reads_raw              #number of loci per read
av_loci <- mean(loci_per_read, na.rm=TRUE)               #0.0021
lpr <- loci_per_read[!is.na(loci_per_read)]            #to get rid of the Na that is messing up my funciton
expected_loci <- reads_raw*av_loci
observed_vs_expected <- loci_in_assembly/expected_loci
ove <- observed_vs_expected[!is.na(observed_vs_expected)] #to get rid of the Na that is messing up my funciton

count_ove <- function(nums){
  x=0
  for(i in nums){
    if(i < .70){
      print(i)
      x=x+1
    }
  }
  print("Samples with less loci than expected (70%):")
  return(x)
}
count_ove(ove)                  #gives number of samples have less loci than expected

##Now I need to put the three filters together, to see how many/which samples I am left with:
##list the samples with less than 1750 loci (we don't need to rerun the good samples)
L = stats$loci_in_assembly <= 1750 
passed1 = stats[L,]                #list of individuals that passed first filter (basically a new dataframe)

##of those-list the samples with more than 400,000 reads (don't want to rerun samples that won't even reach ~1,000,000 reads when doubled.)
R = passed1$reads_raw >= 400000
passed2 = passed1[R,]               #list of individuals that passed second filter

#of those-list the samples with >70% expected loci (not worth rerunning, we won't get many more loci)
E = passed2$reads_raw >= ?
E
passed3 = stats[R,]
passed3                           #these are the candidates for the second library!
