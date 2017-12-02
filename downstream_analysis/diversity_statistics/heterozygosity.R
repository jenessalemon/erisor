## This script first graphs estimated hereozygosity from the ipyrad stats files. 
## Then, I wrote my own script to estimate heterozygosity across individuals and across 

## Part one- graphing heterozygosity from the ipyrad outfiles.
setwd('/Users/jimblotter/Desktop/Grad_School/Data_Analysis/erisor/QC/')
sored_stats <- read.table("lowmsl_sored_stats.txt", header = TRUE, fill = TRUE)
sored_stats
shock_stats <- read.table("lowmsl_shock_stats.txt", header = TRUE, fill = TRUE)
shock_stats
#both_stats <- read.table("combined_177_stats.txt", header = TRUE, fill = TRUE)

## First we plot shockelyi
ggplot(data=shock_stats, aes(shock_stats$hetero_est)) + 
  geom_histogram(breaks=seq(0.004,0.028, by=0.001), 
                 col="black", 
                 fill="#458B74", 
                 alpha = 1) + 
  ggtitle("Estimated Heterozygosity in E. shockleyi") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(plot.title = element_text(size=22)) +
  theme(text = element_text(size=18)) +
  labs(x="Frequency of Individuals", y="Estimated Heterozygosity") + 
  xlim(c(0.003,0.027)) + 
  ylim(c(0,17)) +
  #geom_vline(xintercept = mean(shock_stats$hetero_est), col = "black", size = 1.7) +
  geom_vline(xintercept = mean(sored_stats$hetero_est), col = "#FF7F50", size = 1.7)

## Then soredium
sored_het <- ggplot(data=sored_stats, aes(sored_stats$hetero_est)) + 
  geom_histogram(breaks=seq(0.004,0.028, by=0.001), 
                 col="black", 
                 fill="#FF7F50", 
                 alpha = 1) + 
  ggtitle("Estimated Heterozygosity in E. soredium") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(plot.title = element_text(size=22)) +
  theme(text = element_text(size=18)) +
  labs(x="Frequency of Individuals", y="Estimated Heterozygosity") + 
  xlim(c(0.003,0.027)) + 
  ylim(c(0,17)) +
  #geom_vline(xintercept = mean(sored_stats$hetero_est), col = "black", size = 1.7) +
  geom_vline(xintercept = mean(shock_stats$hetero_est), col = "#458B74", size = 1.7)

## Now we do a t-test to determine the difference between the two means
t.test(shock_hetero, sored_hetero)
"   	Welch Two Sample t-test

data:  shock_hetero and sored_hetero
t = 5.3594, df = 64.714, p-value = 1.19e-06
alternative hypothesis: true difference in means is not equal to 0
95 percent confidence interval:
0.001996189 0.004367938
sample estimates:
mean of x  mean of y 
0.01588685 0.01270478 
"

## Part 2- Getting heterozygosity estimates from a genind object.
setwd('/Users/jimblotter/Desktop/Grad_School/Data_Analysis/erisor/QC') #remember to put the .str file there!
library("ape")
library("genetics")    
library("seqinr")
library("ggplot2")
library("adegenet")
#SUBTRACT ONE FROM NUMBER OF LOCI
#read in data
obj_155_shock <- read.structure("10msl_shock.str", n.ind = 118, onerowperind= FALSE, n.loc = 41436, col.lab = 1, col.pop = 0, col.others = NULL, row.marknames = 0)
indNames(obj_155_shock) 
ploidy(obj_155) # should return 2 since we gave it 2 alleles for each marker.
obj_155_soredium <- read.structure("medmsl_sored.str", n.ind = 37, onerowperind= FALSE, n.loc = 49200, col.lab = 1, col.pop = 0, col.others = NULL, row.marknames = 0)
indNames(obj_155_soredium) # I only get 247 because (24) low coverage individuals get filtered out.
ploidy(obj_155_soredium)
#Need to build a dataframe with individuals as rows, loci as columns and status as the values.
df <- as.data.frame(obj_155)
df_soredium <- as.data.frame(obj_155_soredium)

################## For each Individual (also provided by ipyrad) ###################
heterozygotes <- function(df_genind){
  print("Percentage of Heterozygous individuals across loci:")
  df_genind[is.na(df_genind)] <- 9                       #convert NAs to 9s
  heterozygosities <- c()
  for(j in 1:nrow(df_genind)){                           #for each individual
    row <- df_genind[j,]                             #set row = to the ith row of the dataframe
    count1 <- 0                                        #set the counts to 1
    count2 <- 0
      for(i in row){
        if(i == 1){                                      #count the heterozygous alleles
          count1 <- count1 + 1
        }
        if(i == 2){                                      #count the homozygous alleles
          count2 <- count2 + 1
        }
      }
    hetero_loci <- count1/2
    homo_loci <- count2
    ind <- hetero_loci/(hetero_loci + homo_loci)
    heterozygosities[j] <- ind
    #print(rownames(df_genind)[row])
  }
  j <- j+1
  return(heterozygosities)
}

#Call on shockleyi and soredium
shock_heterozygosity <- heterozygotes(df)
shock_heterozygosity
mean(shock_heterozygosity)

sored_heterozygosity <- heterozygotes(df_soredium)
sored_heterozygosity
mean(sored_heterozygosity)

hist(shock_heterozygosity, breaks = 20, col = "blue", ylab = "Frequency of Individuals", main = "shockleyi blue, soredium overlaid in pink, purple is the overlap", xlab = "Heterozygosity", xlim = c(0,0.13))
hist(sored_heterozygosity, breaks = 20, col=rgb(1,0,0,0.5), add=T)

t.test(sored_hetero,shock_hetero)

########################### FOR EACH LOCUS ########################################
#Need to get hetero/(hetero + homo) for each individual
individuals <- function(df_genind){
  print("Percentage of Heterozygous loci across individuals:")
  df_genind[is.na(df_genind)] <- 9                       #convert NAs to 9s
  heterozygosities <- c()
  for(j in 1:ncol(df_genind)){                           #for each individual
      col <- df_genind[,j]                               #set row = to the ith row of the dataframe
      count1 <- 0                                        #set the counts to 0
      count2 <- 0
      for(i in col){
        if(i == 1){                                      #count the heterozygous alleles
          count1 <- count1 + 1
        }
        if(i == 2 | i == 0){                             #count the homozygous alleles
          count2 <- count2 + 1
        }
    hetero_inds <- count1
    homo_inds <- count2
    loc <- hetero_inds/(hetero_inds + homo_inds)
    heterozygosities[j] <- loc
    }
  }
  heterozygosities <- heterozygosities[!is.na(heterozygosities)]
  return(heterozygosities)
}
individuals(df_soredium)

#Call on shockleyi and soredium
shock_hetero_inds <- individuals(df)
shock_hetero_inds
mean(shock_hetero_inds)

sored_hetero_inds <- individuals(df_soredium)
sored_hetero_inds
mean(sored_hetero_inds)

hist(shock_hetero_inds, breaks = 20, col = "blue", ylab = "Frequency of Individuals", main = "shockleyi blue, soredium overlaid in pink, purple is the overlap", xlab = "Heterozygosity", xlim = c(0,0.13))
hist(sored_hetero_inds, breaks = 20, col=rgb(1,0,0,0.5), add=T)

#################### Obs vs Exp Heterozygostity ############################
# Read in the data
obj_155_soredium <- read.structure('/Users/jimblotter/Desktop/Grad_School/Data_Analysis/erisor/QC/c155_sored.str', n.ind = 37, onerowperind= FALSE, n.loc = 942, col.lab = 1, col.pop = 0, col.others = NULL, row.marknames = 0)
obj_155_soredium
obj_155_shockleyi <- read.structure('/Users/jimblotter/Desktop/Grad_School/Data_Analysis/erisor/QC/c155_shock.str', n.ind = 118, onerowperind= FALSE, n.loc = 337, col.lab = 1, col.pop = 0, col.others = NULL, row.marknames = 0)
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







