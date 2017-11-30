## This script first graphs estimated hereozygosity from the ipyrad stats files. 
## Then, I wrote my own script to estimate heterozygosity across individuals and across 

## Part one- graphing heterozygosity from the ipyrad outfiles.
setwd('/Users/jimblotter/Desktop/Grad_School/Data_Analysis/erisor/QC/')
sored_stats <- read.table("lowmsl_sored_stats.txt", header = TRUE, fill = TRUE)
sored_stats
shock_stats <- read.table("lowmsl_shock_stats.txt", header = TRUE, fill = TRUE)
shock_stats
#both_stats <- read.table("combined_177_stats.txt", header = TRUE, fill = TRUE)

sored_hetero <- sored_stats$hetero_est
shock_hetero <- shock_stats$hetero_est

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
axis(2, at=c(0,17), labels=c("",""), lwd.ticks=0)
axis(2, at=seq(0,17, by=5), lwd=0, lwd.ticks=1)






mean(sored_hetero)
sd(sored_hetero)
mean(shock_hetero)
sd(shock_hetero)

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

hist(shock_hetero, breaks = 20, col = "blue", main = "shockleyi blue, soredium overlaid in pink", xlab = "Heterozygosity")
hist(sored_hetero, breaks = 20, col=rgb(1,0,0,0.5), add=T)

t.test(sored_hetero,shock_hetero)

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
div <- summary(obj_155)
plot(div$Hobs, xlab="Loci number", ylab="Observed Heterozygosity", 
     main="Observed Heterozygosity per Locus - E.shockleyi")
plot(div$Hexp, xlab="Loci number", ylab="Expected Heterozygosity", 
     main="Expected Heterozygosity per Locus - E.shockleyi")
plot(div$Hobs,div$Hexp, xlab="Hobs", ylab="Hexp", 
     main="Exp Heterozygosity as a function of Obs heterozygosity per locus-E.shockleyi")

div <- summary(obj_155_soredium)
plot(div$Hobs, xlab="Loci number", ylab="Observed Heterozygosity", 
     main="Observed Heterozygosity per Locus - E.soredium")
plot(div$Hexp, xlab="Loci number", ylab="Expected Heterozygosity", 
     main="Expected Heterozygosity per Locus - E.soredium")
plot(div$Hobs,div$Hexp, xlab="Hobs", ylab="Hexp", 
     main="Exp Heterozygosity as a function of Obs heterozygosity per locus-E.soredium")



both_stats <- read.table("/Users/jimblotter/Desktop/Grad_School/Data_Analysis/erisor/QC/final163_stats.txt", header = TRUE, fill = TRUE)
oops <- both_stats$error_est
mean(oops)
mean(both_stats$reads_consens)



