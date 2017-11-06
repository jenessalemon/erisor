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

##Need levels of heterzygous individuals not loci
setwd('/Users/jimblotter/Desktop/Grad_School/Data_Analysis/erisor/QC') #remember to put the .str file there!
library("ape")
library("genetics")    
library("seqinr")
library("ggplot2")
library("adegenet")
#SUBTRACT ONE FROM NUMBER OF LOCI

#read in data
obj_155 <- read.structure("c155_shock.str", n.ind = 118, onerowperind= FALSE, n.loc = 330, col.lab = 1, col.pop = 0, col.others = NULL, row.marknames = 0)
indNames(obj_155) # I only get 247 because (24) low coverage individuals get filtered out.
ploidy(obj_155) # should return 2 since we gave it 2 alleles for each marker.
#Need to build a dataframe with individuals as rows, loci as columns and status as the values.
df <- as.data.frame((obj_155))
df$L001.01
 
#Need to get hetero/(hetero + homo) for each individual
heterozygotes <- function(df_genind){
  print("Percentage of Heterozygous loci:")
  for(row in df_genind){      #for each individual
    count1 <- 0               #set the counts to 1
    count2 <- 0
    for(i in row){            #for each allele at each locus
      if(is.na(i)){           #set NAs to 9
        i <- 9
      }
      if(i == 1){             #count the heterozygous alleles
        count1 <- count1 + 1
      }
      if(i == 2){            #count the homozygous alleles
        count2 <- count2 + 1
      }
      #print("Number of Heterozygous loci:")
      #print(count1)
      #print("Number of Homozygous loci:")
      #print(count2)
      #print(rownames(df_genind)[row])
      print(count1/(count1 + count2)) #calculate heterozygosity
    }
  }
}
heterozygotes(df)

