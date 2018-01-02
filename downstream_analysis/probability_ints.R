## Plot STRUCTURE point estimates with 90% confidence intervals.
library("readr")
library("ggplot2")

## prob_intervals is obtained from the "inferred ancestry of individuals" section of (cat structure-test-K-2-rep-0_f) found in structure_out, the Jupyter notebook working directory at the CHPC.
## Read this file in, and manipulate it using R dataframes.
prob_intervals <- read.table("/Users/jimblotter/Desktop/Grad_School/Data_Analysis/erisor/downstream_analysis/probability_intervals.txt", header = TRUE, fill = TRUE)
prob_intervals$Label <- NULL
prob_intervals$Inferred <- NULL
prob_intervals[,7] <- NULL
names(prob_intervals)[1] <- "sample"
names(prob_intervals)[2] <- "percent_missing"
names(prob_intervals)[3] <- "cluster1"
names(prob_intervals)[4] <- "cluster2"
names(prob_intervals)[5] <- "interval_low"
names(prob_intervals)[6] <- "interval_high"

## Now, unfortunately, you have to export to separate the confidence intervals into two separate values.
write_csv(prob_intervals, path = '/Users/jimblotter/Desktop/Grad_School/Data_Analysis/erisor/downstream_analysis/probs_int.csv')

## Now upload this to google sheets "probability intervals", run that, and read it back in.
prob_int <- read.csv("/Users/jimblotter/Desktop/Grad_School/Data_Analysis/erisor/downstream_analysis/probability_intervals_bargraph.csv", header = TRUE, fill = TRUE)
names(prob_int)[3] <- "point_est"
prob_int

## Create a little dataframe.
pidf <- data.frame(
  sample = prob_int$sample,
  point = prob_int$point_est,
  upper = prob_int$interval_high,
  lower = prob_int$interval_low
)

## Plot the dataframe, and add errorbars and points.
p <- ggplot(pidf, aes(sample, point))
p + geom_errorbar(aes(ymin = lower, ymax = upper)) + geom_pointrange(aes(ymin = lower, ymax = upper), size = .13)
  
  