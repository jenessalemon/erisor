## Plot STRUCTURE with a bar graph, including confidence intervals.
library("readr")
library("ggplot2")
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
write_csv(prob_intervals, path = '/Users/jimblotter/Desktop/Grad_School/Data_Analysis/erisor/downstream_analysis/probs_int.csv')

prob_int <- read.csv("/Users/jimblotter/Desktop/Grad_School/Data_Analysis/erisor/downstream_analysis/probability_intervals_bargraph.csv", header = TRUE, fill = TRUE)
names(prob_int)[3] <- "point_est"
prob_int


ggplot(data=prob_int, aes(x=prob_int$sample, y=prob_int$point_est)) +
  geom_col(prob_int$point_est) +
  geom_errorbar(aes(ymin = prob_int$interval_low, ymax = prob_int$interval_high, width=0.1))


pidf <- data.frame(
  sample = prob_int$sample,
  point = prob_int$point_est,
  upper = prob_int$interval_high,
  lower = prob_int$interval_low
)
  
p <- ggplot(pidf, aes(sample, point))
p + geom_pointrange(aes(ymin = lower, ymax = upper))
p + geom_errorbar(aes(ymin = lower, ymax = upper)) + geom_pointrange(aes(ymin = lower, ymax = upper), size = .13)
  
  