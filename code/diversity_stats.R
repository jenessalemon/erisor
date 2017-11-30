## This script estimates nucleotide diversity, and compares expected levels
## of heterozygosity.
## Input files are .vcf, and the ipyrad stats files.

#################### Overall Nucleotide Diversity ###############################
library(PopGenome)
# Read in the data. For some reason, the data must be in a folder. Just put it in a folder by itself.
GENOME.class <- readData('/Users/jimblotter/Desktop/Grad_School/Data_Analysis/erisor/QC/my_folder', format="VCF") 

# Set population data. Here I just want to compare diversity between my two species, so I treat each species as a "population."
GENOME.class <- set.populations(GENOME.class,     
                                list(c("p_004.5s_011","p_004.5s_05","p_004.5s_07","p_004.5s_08","p_004.5s_12","p_004.5s_15","p_004s_04","p_005s_02","p_005s_04","p_005s_06","p_005s_07","p_005s_11","p_006s_01","p_006s_04","p_006s_06","p_006s_10","p_006s_12","p_006s_14","p_007s_01","p_007s_02","p_007s_04","p_007s_06","p_007s_07","p_007s_09","p_008.5s_02","p_009s_01","p_009s_09","p_009s_12","p_009s_13","p_009s_14","p_009s_15","p_010s_02","p_010s_06","p_010s_11","p_011s_06","p_011s_08","p_011s_10","p_011s_11","p_012s_01","p_012s_06","p_012s_08","p_012s_12","p_012s_14","p_013s_02","p_013s_08","p_013s_12","p_013s_13","p_013s_14","p_014s_01","p_014s_02","p_014s_03","p_014s_11","p_014s_13","p_014s_15","p_015s_09","p_015s_13","p_015s_14","p_016s_01","p_016s_02","p_016s_03","p_016s_04","p_016s_06","p_016s_07","p_016s_14","p_016s_15","p_017s_01","p_017s_03","p_017s_05","p_017s_06","p_017s_07","p_018s_02","p_018s_06","p_018s_09","p_018s_14","p_019s_01","p_019s_03","p_019s_12","p_020s_01","p_020s_05","p_020s_11","p_020s_13","p_021s_01","p_021s_05","p_021s_08","p_021s_09","p_021s_12","p_021s_14","p_023s_03","p_024s_01","p_024s_06","p_024s_11","p_024s_12","p_025s_13","p_026s_04","p_026s_05","p_026s_12","p_026s_14","p_027s_01","p_027s_03","p_027s_05","p_027s_06","p_027s_12","p_030s_01","p_030s_02","p_030s_04","p_030s_05","p_030s_06","p_030s_07","p_030s_11","p_032s_05","p_1027s_02","p_1027s_04","p_1027s_08","p_1027s_09","p_1027s_10","p_1027s_12","p_1027s_18","p_1027s_20"), c("p_001s_02","p_001s_03","p_001s_09","p_001s_12","p_001s_13","p_001s_14","p_001s_16","p_002.5s_01","p_002.5s_04","p_002.5s_05","p_002.5s_07","p_002.5s_09","p_002.5s_10","p_002s_01","p_002s_03","p_002s_05","p_002s_06","p_002s_08","p_002s_09","p_002s_10","p_002s_12","p_002s_13","p_002s_14","p_002s_16","p_029s_02","p_029s_03","p_029s_04","p_029s_05","p_029s_09","p_029s_10","p_031s_02","p_031s_04","p_031s_05","p_031s_07","p_031s_08","p_031s_11","p_031s_14")), diploid=TRUE)
GENOME.class@populations

# Get nucleotide diversity.
GENOME.class <- diversity.stats(GENOME.class)     
get.diversity(GENOME.class)
GENOME.class@nuc.diversity.within                 #returns the within population (shockleyi and soredium are the two, collective "populations," nucleotide diversity.
#pop 1 (E. shockeyi)    pop 2 (E. soredium)
#[1,] 1.822106          0.2132544

################### Estimated Heterozygosity ###################################
# Read in the estimated heterozygosity from the ipyrad stats files. (The two species were assembled separately in iPyRAD.)
setwd('')
sored_stats <- read.table("lowmsl_sored_stats.txt", header = TRUE, fill = TRUE)
shock_stats <- read.table("lowmsl_shock_stats.txt", header = TRUE, fill = TRUE)
sored_hetero <- sored_stats$hetero_est
shock_hetero <- shock_stats$hetero_est

# Now we simply generate the histogram of frequency of individuals with certain 
# levels of expected heterozygosity. The following function adds opacity to the 
# overlay so that we can see E. shockleyi's histogram underneath.
add.alpha <- function(col, alpha=1){
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2, 
        function(x) 
          rgb(x[1], x[2], x[3], alpha=alpha))  
}
teal <- add.alpha("#458B74", alpha = 1)
orange <- add.alpha("#FF7F50", alpha = 0.8)

hist(shock_hetero, breaks=seq(0.004,0.028, by=0.001), col = teal, main = "Estimated Heterozygosity of Individuals", xlab = "Estimated Heterozygosity", ylab = "Frequency of Indiviudals", xlim = c(0.003,0.027), xaxt="n")
hist(sored_hetero, breaks=seq(0.004,0.028, by=0.001), col= orange, xlim = c(0.004,0.027), xaxt="n", add=T)

axis(1, at=c(0.004,0.027), labels=c("",""), lwd.ticks=0)
axis(1, at=seq(0.004,0.027, by=.001), lwd=0, lwd.ticks=1)

##################### EXPECTED HETEROZYGOSITY ########################
obj_155_soredium <- read.structure('/Users/jimblotter/Desktop/Grad_School/Data_Analysis/erisor/QC/medmsl_sored.str', n.ind = 37, onerowperind= FALSE, n.loc = 49201, col.lab = 1, col.pop = 0, col.others = NULL, row.marknames = 0)
obj_155_soredium
obj_155_shockleyi <- read.structure('/Users/jimblotter/Desktop/Grad_School/Data_Analysis/erisor/QC/medmsl_shock.str', n.ind = 118, onerowperind= FALSE, n.loc = 13983, col.lab = 1, col.pop = 0, col.others = NULL, row.marknames = 0)
obj_155_shockleyi

div_sored <- summary(obj_155_soredium)
div_shock <- summary(obj_155_shockleyi)

plot(div_sored$Hobs,div_sored$Hexp, xlab="Hobs", ylab="Hexp", 
     main="Exp Heterozygosity as a function of Obs heterozygosity per locus-E.soredium")

plot(div_shock$Hobs,div_shock$Hexp, xlab="Hobs", ylab="Hexp", 
     main="Exp Heterozygosity as a function of Obs heterozygosity per locus-E.shockleyi")


barplot(div_sored$Hexp-div_sored$Hobs, main="Heterozygosity: expected-observed",
        ylab="Hexp - Hobs")
t.test(div_sored$Hexp,div_sored$Hobs,pair=T,var.equal=TRUE,alter="greater")
"          Paired t-test

data:  div_sored$Hexp and div_sored$Hobs
t = 18.567, df = 941, p-value < 2.2e-16
alternative hypothesis: true difference in means is greater than 0
95 percent confidence interval:
0.05429225        Inf
sample estimates:
mean of the differences 
0.05957536 "

barplot(div_shock$Hexp-div_shock$Hobs, main="Heterozygosity: expected-observed",
        ylab="Hexp - Hobs")
t.test(div_shock$Hexp,div_shock$Hobs,pair=T,var.equal=TRUE,alter="greater")
"	        Paired t-test

data:  div_shock$Hexp and div_shock$Hobs
t = 8.1937, df = 336, p-value = 2.689e-15
alternative hypothesis: true difference in means is greater than 0
95 percent confidence interval:
0.02573512        Inf
sample estimates:
mean of the differences 
0.03222134 "



summary(div_sored$Hobs)
summary(div_shock$Hobs)
summary(div_sored$Hexp)
summary(div_shock$Hexp)

mean(div_sored$Hexp) #0.1190287
mean(div_shock$Hexp) #0.05072056
mean(div_sored$Hexp) - mean(div_sored$Hobs) #0.05957536
mean(div_shock$Hexp) - mean(div_shock$Hobs) #0.03222134

################# Diversity statistics in PopGenome ##################
#Read in the data, which needs to be in its own folder for some reason.
GENOME.class <- readData('/Users/jimblotter/Desktop/Grad_School/Data_Analysis/erisor/QC/my_folder', format="VCF") 

# Set population data. Here I just want to compare diversity between my two species, so I treat each species as a "population."
GENOME.class <- set.populations(GENOME.class,     
                                list(c("p_004.5s_011","p_004.5s_05","p_004.5s_07","p_004.5s_08","p_004.5s_12","p_004.5s_15","p_004s_04","p_005s_02","p_005s_04","p_005s_06","p_005s_07","p_005s_11","p_006s_01","p_006s_04","p_006s_06","p_006s_10","p_006s_12","p_006s_14","p_007s_01","p_007s_02","p_007s_04","p_007s_06","p_007s_07","p_007s_09","p_008.5s_02","p_009s_01","p_009s_09","p_009s_12","p_009s_13","p_009s_14","p_009s_15","p_010s_02","p_010s_06","p_010s_11","p_011s_06","p_011s_08","p_011s_10","p_011s_11","p_012s_01","p_012s_06","p_012s_08","p_012s_12","p_012s_14","p_013s_02","p_013s_08","p_013s_12","p_013s_13","p_013s_14","p_014s_01","p_014s_02","p_014s_03","p_014s_11","p_014s_13","p_014s_15","p_015s_09","p_015s_13","p_015s_14","p_016s_01","p_016s_02","p_016s_03","p_016s_04","p_016s_06","p_016s_07","p_016s_14","p_016s_15","p_017s_01","p_017s_03","p_017s_05","p_017s_06","p_017s_07","p_018s_02","p_018s_06","p_018s_09","p_018s_14","p_019s_01","p_019s_03","p_019s_12","p_020s_01","p_020s_05","p_020s_11","p_020s_13","p_021s_01","p_021s_05","p_021s_08","p_021s_09","p_021s_12","p_021s_14","p_023s_03","p_024s_01","p_024s_06","p_024s_11","p_024s_12","p_025s_13","p_026s_04","p_026s_05","p_026s_12","p_026s_14","p_027s_01","p_027s_03","p_027s_05","p_027s_06","p_027s_12","p_030s_01","p_030s_02","p_030s_04","p_030s_05","p_030s_06","p_030s_07","p_030s_11","p_032s_05","p_1027s_02","p_1027s_04","p_1027s_08","p_1027s_09","p_1027s_10","p_1027s_12","p_1027s_18","p_1027s_20"), c("p_001s_02","p_001s_03","p_001s_09","p_001s_12","p_001s_13","p_001s_14","p_001s_16","p_002.5s_01","p_002.5s_04","p_002.5s_05","p_002.5s_07","p_002.5s_09","p_002.5s_10","p_002s_01","p_002s_03","p_002s_05","p_002s_06","p_002s_08","p_002s_09","p_002s_10","p_002s_12","p_002s_13","p_002s_14","p_002s_16","p_029s_02","p_029s_03","p_029s_04","p_029s_05","p_029s_09","p_029s_10","p_031s_02","p_031s_04","p_031s_05","p_031s_07","p_031s_08","p_031s_11","p_031s_14")), diploid=TRUE)
# Check population assignment:
GENOME.class@populations 

# More divergence statistics:
GENOME.class <- F_ST.stats(GENOME.class)
get.F_ST(GENOME.class, pairwise = TRUE)
GENOME.class@nucleotide.F_ST                #0.07481693
GENOME.class@haplotype.F_ST                 #0.2055677
GENOME.class@Nei.G_ST                       #0.1172353  
GENOME.class@Hudson.G_ST                    #0.06442249

##################### Diversity stats in mmod ##############################
## Convert structure file to genind
obj_155_soredium <- read.structure('/Users/jimblotter/Desktop/Grad_School/Data_Analysis/erisor/QC/medmsl_sored.str', n.ind = 37, onerowperind= FALSE, n.loc = 49201, col.lab = 1, col.pop = 0, col.others = NULL, row.marknames = 0)
obj_155_soredium
obj_155_shockleyi <- read.structure('/Users/jimblotter/Desktop/Grad_School/Data_Analysis/erisor/QC/medmsl_shock.str', n.ind = 118, onerowperind= FALSE, n.loc = 13983, col.lab = 1, col.pop = 0, col.others = NULL, row.marknames = 0)
obj_155_shockleyi

## Remove a locus with no data
genotyped_soredium <- obj_155_soredium[,locFac(obj_155_d) != "L011"]
genotyped_shockleyi <- obj_155_shockleyi[,locFac(obj_155_d) != "L011"]

## Calculate useful statistics:
Gst_Nei(genotyped_soredium) 
Gst_Nei(genotyped_shockleyi) 
Gst_Hedrick(genotyped_soredium) 
Gst_Hedrick(genotyped_shockleyi)
diff_stats(genotyped_soredium)
diff_stats(genotyped_shockleyi)





