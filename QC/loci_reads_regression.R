setwd('/Users/jimblotter/Desktop/Grad_School/Data_Analysis/QC/')
stats <- read.table("erisor_reps_stats.txt", header = TRUE, fill = TRUE)
reads_raw <- stats[,2]
loci_in_assembly <- stats[,9]

reads_and_loci <- subset(stats, select = c("reads_raw", "loci_in_assembly"))
summary(reads_and_loci)
plot(reads_and_loci)

fit <- lm(reads_raw ~ loci_in_assembly, data=stats)
summary(fit)
plot(fit)