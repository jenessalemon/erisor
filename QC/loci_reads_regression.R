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

###########################################################################################
#add a column to passed2?
step3 <- passed2[c(1,2,3)] = multiple columns

loci <- passed2$loci_in_assembly
reads <- passed2$reads_raw
loci_per_read <- loci/reads
expected_loci <- reads*0.0021
observed_vs_expected <- loci/expected_loci

seventy <- observed_vs_expected >= .7



