## The relevant _stats.txt file from ipyrad is the input file for this script.
## In this script I plot reads vs. loci. I evaluate the linear regression. I apply
## filters to the data and evaluate how many individuals from each population will
## be retained if I cluster with these new filters. I decided on good traige parameter
## settings int this script, and looked at the plate distribution of poor reads to
## discover two rows of pipetting error.

setwd('/Users/jimblotter/Desktop/Grad_School/Data_Analysis/erisor/QC/')
stats <- read.table("p1_stats_now.txt", header = TRUE, fill = TRUE)
str(stats)

################################## Plotting #######################################
##Simple plot before regression
reads_and_loci <- subset(stats, select = c("reads_raw", "loci_in_assembly"))
summary(reads_and_loci)
plot(reads_and_loci, main= "p1")
#Plot with coloration by species and populations: first, get the population names from the sample names:
sample_name_elements <- strsplit(row.names(stats), "_")[1:262]
stats$population <- as.factor(sapply(sample_name_elements, function(x) x[[2]]))
attach(stats)
# Plot loci vs reads colored by population
plot(reads_raw, loci_in_assembly, col=population, main="Color by population")
# Plot loci vs reads colored by species
plot(reads_raw, loci_in_assembly, col=species, main="Color by species")

##Linear Model
fit <- lm(reads_raw ~ loci_in_assembly, data=stats)   #linear regression
summary(fit)
plot(fit)

############## How many good samples does each population have? ##########################################
L = stats$loci_in_assembly >= 230 #apply filter, was 200
good_loci = stats[L,]              #list of just the samples with good loci
R = good_loci$reads_raw >= 400000  #apply filter to the list of samples with good loci was 350,000
good_loci_and_reads = good_loci[R,]         #list of just the samples with good loci and reads
good_samples = row.names(good_loci_and_reads)              #these are the good samples I have for my project
nrow(good_loci_and_reads)

#I want a function that will tell me how many good individuals I have from each population.
sample_name_elements <- strsplit(row.names(good_loci_and_reads), "_")[1:nrow(good_loci_and_reads)]  #split the string
inds_mat <- matrix(unlist(sample_name_elements), ncol=3, byrow=TRUE)                                #make it a matrix so you can index it
pops <- inds_mat[,2]                            #second column of the new matrix is the population name (followed by an s but who cares)
table(pops)                                     #thanks stack overflow

############ Filtering to determine strategic re-runs ############
#Now I need to apply three (or 4) filters, to see how many/which candidates for resequencing I am left with:

##list the samples with less than 1750 loci (we don't need to rerun the good samples)
filter1 = stats$loci_in_assembly <= 1200    #apply filter
passed1 = stats[filter1,]                   #list of individuals that passed first filter (basically a new dataframe)

##list the samples more than 600 loci
filter1.5 = passed1$loci_in_assembly >= 600
passed1.5 = passed1[filter1.5,]

##of those-list the samples with more than 400,000 reads (don't want to rerun samples that won't even reach ~1,000,000 reads when doubled.)
filter2 = passed1.5$reads_raw >= 200000
passed2 = passed1.5[filter2,]               #list of individuals that passed second filter

#of those-list the samples with >70% expected loci (not worth rerunning, we won't get many more loci)
#simple setup calculations
loci <- passed2$loci_in_assembly       #get the loci counts
reads <- passed2$reads_raw             #get the read counts
loci_per_read <- loci/reads            #how many reads per loci?
expected_loci <- reads*0.0021          #reads*average # of loci per read gives the expected # of reads (mean was calculated earlier)
observed_vs_expected <- loci/expected_loci                      #observed over expected
#now add the observed vs expected column to the filter matrix, and filter it for the third and last time.
passed2$obs_vs_exp<-observed_vs_expected
filter3 = passed2$obs_vs_exp >= .70
reruns = passed2[filter3,]
rerun_names = row.names(reruns)                       #these are the candidates for the second library! But we must check for repeats with the good samples.
#check for replicates with the good samples, and remove:
exclusive <- function(candidates, good_samps){
  i <- 2
  for(i in candidates){
    while (i<=length(candidates)){
      if (i %in% good_samps) { 
        candidates<-candidates[-i,]
      }
      i <- i+1
    }
  }
  print("These are eligible reruns")
  return(candidates)
}
exclusive(rerun_names, good_samples)
################################## Plate distribution of low reads ########################################
B = stats$reads_raw <= 40000 #create a new dataframe with only individuals 10000 and lower
low_reads = stats[B,]
low_reads             

plate1 <- matrix(nrow = 8, ncol = 12)   #initialize 3 empty matrices resembling the plates
plate2 <- matrix(nrow = 8, ncol = 12)
plate3 <- matrix(nrow = 8, ncol = 12)

plate3[6,6] <- "X"   #Plot low coverage samples. I just did this by hand it didn't take long.
plate3[4,3] <- "X"
plate3[5,6] <- "X"
plate1[1,4] <- "X"
plate3[7,4] <- "X"
plate3[3,6] <- "X"
plate3[1,6] <- "X"
plate1[7,12] <- "X"
plate1[8,10] <- "X"
plate3[2,1] <- "X"
plate3[2,2] <- "X"
plate3[8,7] <- "X"
plate3[2,4] <- "X"
plate3[6,2] <- "X"
plate3[8,2] <- "X"
plate3[8,11] <- "X"
plate3[1,10] <- "X"
plate1[8,5] <- "X"
plate3[5,2] <- "X"
plate3[4,6] <- "X"
plate3[3,2] <- "X"
plate3[7,2] <- "X"
plate3[1,2] <- "X"
plate3[7,6] <- "X"
plate2[8,9] <- "X"
plate3[8,8] <- "X"


plate1
plate2
plate3 

reads = stats$reads_raw
sum(reads)

good_samps <- data.frame(good_samples)
write_tsv(good_samps, path = '/Users/jimblotter/Desktop/Grad_School/Data_Analysis/erisor/QC/hi')


