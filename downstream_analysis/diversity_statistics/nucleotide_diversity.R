#################### Overall Nucleotide Diversity ###############################
library(PopGenome)
GENOME.class <- readData("/Users/jimblotter/Desktop/Grad_School/Data_Analysis/erisor/QC/my_folder", format="VCF") #read in the data
GENOME.class <- set.populations(GENOME.class,     #set populations
                                list(c("p_004.5s_011","p_004.5s_05","p_004.5s_07","p_004.5s_08","p_004.5s_12","p_004.5s_15","p_004s_04","p_005s_02","p_005s_04","p_005s_06","p_005s_07","p_005s_11","p_006s_01","p_006s_04","p_006s_06","p_006s_10","p_006s_12","p_006s_14","p_007s_01","p_007s_02","p_007s_04","p_007s_06","p_007s_07","p_007s_09","p_008.5s_02","p_009s_01","p_009s_09","p_009s_12","p_009s_13","p_009s_14","p_009s_15","p_010s_02","p_010s_06","p_010s_11","p_011s_06","p_011s_08","p_011s_10","p_011s_11","p_012s_01","p_012s_06","p_012s_08","p_012s_12","p_012s_14","p_013s_02","p_013s_08","p_013s_12","p_013s_13","p_013s_14","p_014s_01","p_014s_02","p_014s_03","p_014s_11","p_014s_13","p_014s_15","p_015s_09","p_015s_13","p_015s_14","p_016s_01","p_016s_02","p_016s_03","p_016s_04","p_016s_06","p_016s_07","p_016s_14","p_016s_15","p_017s_01","p_017s_03","p_017s_05","p_017s_06","p_017s_07","p_018s_02","p_018s_06","p_018s_09","p_018s_14","p_019s_01","p_019s_03","p_019s_12","p_020s_01","p_020s_05","p_020s_11","p_020s_13","p_021s_01","p_021s_05","p_021s_08","p_021s_09","p_021s_12","p_021s_14","p_023s_03","p_024s_01","p_024s_06","p_024s_11","p_024s_12","p_025s_13","p_026s_04","p_026s_05","p_026s_12","p_026s_14","p_027s_01","p_027s_03","p_027s_05","p_027s_06","p_027s_12","p_030s_01","p_030s_02","p_030s_04","p_030s_05","p_030s_06","p_030s_07","p_030s_11","p_032s_05","p_1027s_02","p_1027s_04","p_1027s_08","p_1027s_09","p_1027s_10","p_1027s_12","p_1027s_18","p_1027s_20"), c("p_001s_02","p_001s_03","p_001s_09","p_001s_12","p_001s_13","p_001s_14","p_001s_16","p_002.5s_01","p_002.5s_04","p_002.5s_05","p_002.5s_07","p_002.5s_09","p_002.5s_10","p_002s_01","p_002s_03","p_002s_05","p_002s_06","p_002s_08","p_002s_09","p_002s_10","p_002s_12","p_002s_13","p_002s_14","p_002s_16","p_029s_02","p_029s_03","p_029s_04","p_029s_05","p_029s_09","p_029s_10","p_031s_02","p_031s_04","p_031s_05","p_031s_07","p_031s_08","p_031s_11","p_031s_14")), diploid=TRUE)
GENOME.class@populations                          #checking to see if it worked.
GENOME.class <- diversity.stats(GENOME.class)     
get.diversity(GENOME.class)
GENOME.class@nuc.diversity.within                 #returns the within population (shockleyi and soredium are the two, collective "populations," nucleotide diversity.)
#pop 1     pop 2
#[1,] 1.822106 0.2132544

PopGplot(GENOME.class@nuc.diversity.within, colours)

#GENOME.class.slide <- sliding.window.transform(GENOME.class, 10,10, type=1)
#genome.pos <- sapply(GENOME.class.slide@region.names, function(x){
#  split <- strsplit(x," ")[[1]][c(1,3)]
#  val   <- mean(as.numeric(split))
#  return(val)
#})
#plot(genome.pos, <slide.statistic.values>)
#GENOME.class@region.stats
