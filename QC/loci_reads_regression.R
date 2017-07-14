setwd('/Users/jimblotter/Desktop/Grad_School/Data_Analysis/erisor/QC/')
stats <- read.table("erisor_reps_stats-with_spp.txt", header = TRUE, fill = TRUE)
str(stats)

################################## Plotting #######################################
##Simple plot before regression
reads_and_loci <- subset(stats, select = c("reads_raw", "loci_in_assembly"))
summary(reads_and_loci)
plot(reads_and_loci)
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

################################## Filtering #######################################
#Now I need to apply three filters, to see how many/which candidates for resequencing I am left with:

##list the samples with less than 1750 loci (we don't need to rerun the good samples)
filter1 = stats$loci_in_assembly <= 1750    #apply filter
passed1 = stats[filter1,]                   #list of individuals that passed first filter (basically a new dataframe)

##of those-list the samples with more than 400,000 reads (don't want to rerun samples that won't even reach ~1,000,000 reads when doubled.)
filter2 = passed1$reads_raw >= 400000
passed2 = passed1[filter2,]               #list of individuals that passed second filter

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
row.names(reruns)                       #these are the candidates for the second library!

############################# How many good samples does each population have? ##########################################
L = stats$loci_in_assembly >= 1000 #apply filter
good_loci = stats[L,]              #list of just the samples with good loci
R = good_loci$reads_raw >= 500000  #apply filter to the list of samples with good loci
good_loci_and_reads = good_loci[R,]         #list of just the samples with good loci and reads
good_samples = row.names(good_loci_and_reads)              #these are the good samples I have for my project
nrow(good_loci_and_reads)

#I want a function that will tell me how many good individuals I have from each population.
sample_name_elements <- strsplit(row.names(good_loci_and_reads), "_")[1:nrow(good_loci_and_reads)]  #split the string
inds_mat <- matrix(unlist(sample_name_elements), ncol=3, byrow=TRUE)                                #make it a matrix so you can index it
pops <- inds_mat[,2]                            #second column of the new matrix is the population name (followed by an s but who cares)
good_loci_and_reads$population <- pops          #add the new column
good_loci_and_reads                             #print to make sure it worked. it did!
