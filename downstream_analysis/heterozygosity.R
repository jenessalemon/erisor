## This script first graphs estimated hereozygosity from the ipyrad stats files. 
## Then, I wrote my own script to estimate heterozygosity across individuals and across 
## loci using a genind object.

## Part one- graphing heterozygosity from the ipyrad outfiles.
setwd('/Users/jimblotter/Desktop/Grad_School/Data_Analysis/erisor/downstream_analysis/')
sored_stats <- read.table("soredium_c177_stats.txt", header = TRUE, fill = TRUE)
sored_stats
shock_stats <- read.table("shockleyi_c177_stats.txt", header = TRUE, fill = TRUE)
shock_stats
both_stats <- read.table("combined_177_stats.txt", header = TRUE, fill = TRUE)
both_stats

sored_hetero <- sored_stats$hetero_est
shock_hetero <- shock_stats$hetero_est
both_hetero <- both_stats$hetero_est

mean(sored_hetero)
mean(shock_hetero)
mean(both_hetero)


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
obj_155 <- read.structure("c155_shock.str", n.ind = 118, onerowperind= FALSE, n.loc = 330, col.lab = 1, col.pop = 0, col.others = NULL, row.marknames = 0)
indNames(obj_155) 
ploidy(obj_155) # should return 2 since we gave it 2 alleles for each marker.
obj_155_soredium <- read.structure("c155_sored.str", n.ind = 37, onerowperind= FALSE, n.loc = 941, col.lab = 1, col.pop = 0, col.others = NULL, row.marknames = 0)
indNames(obj_155_soredium) # I only get 247 because (24) low coverage individuals get filtered out.
ploidy(obj_155_soredium)
#Need to build a dataframe with individuals as rows, loci as columns and status as the values.
df <- as.data.frame(obj_155)
df_soredium <- as.data.frame(obj_155_soredium)

#Need to get hetero/(hetero + homo) for each individual
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

sored_heterozygosity <- heterozygotes(df_soredium)
sored_heterozygosity

hist(shock_heterozygosity, breaks = 20, col = "blue", ylab = "Frequency of Individuals", main = "shockleyi blue, soredium overlaid in pink, purple is the overlap", xlab = "Heterozygosity", xlim = c(0,0.13))
hist(sored_heterozygosity, breaks = 20, col=rgb(1,0,0,0.5), add=T)

t.test(sored_hetero,shock_hetero)

#######################FOR EACH LOCUS#################################
#Need to get hetero/(hetero + homo) for each individual
individuals <- function(df_genind){
  print("Percentage of Heterozygous loci across individuals:")
  df_genind[is.na(df_genind)] <- 9                       #convert NAs to 9s
  heterozygosities <- c()
  for(j in 1:ncol(df_genind)){                           #for each individual
    if(is.odd(j) == TRUE){                               #I want to single out one allele
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


