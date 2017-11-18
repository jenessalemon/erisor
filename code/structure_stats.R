## This script analysis the Q matrix output (genome composition) from STRUCTURE.
## First, we get the Q matrix. Input is a .geno file.
## We also need a .csv with 5 columns:
## 1. Sample name
## 2. The Q matrix values for the first cluster.
## 3. The Q matrix values for the second cluster.
## 4. The sequencing run that the sample was associated with (optional).
## 5. Morphological identification of the sample.

library(LEA)
obj.snmf = snmf("erisor.geno", K = 2, alpha = 100, project = "new")
qmatrix = Q(obj.snmf, K = 2)

# Plot the structure output as a bar graph.
#barplot(t(qmatrix), col = c("orange","violet","lightgreen","blue"), border = NA, space = 0,
#        xlab = "Individuals", ylab = "Admixture coefficients", main = "msl = 50%") # names.arg = in_order2)

# Convert the Q matrix to a dataframe.
qmat <- as.data.frame(qmatrix)
qmat

# Read in the other columns.
Q <- read.csv("ipy_order.csv", header = TRUE, sep=",")
Q[3] <- qmat[1] 
Q[2] <- qmat[2] 
Q


# This function determines the number of individuals from each species with 
# significant admixture.
Q_admix <- function(df, threshold){
  spec <- c()
  count_sored_admix <- 0
  count_shock_admix <- 0
  count_total_admix <- 0
  count_shock <- 0
  count_sored <- 0
  for(i in 1:nrow(df)){
    if(df[i,5] == "soredium" & df[i,2] < threshold){  #if str says that this individual is soredium, and I specified that it is, call it soredium.
      count_sored_admix <- count_sored_admix + 1
      count_total_admix <- count_total_admix + 1
      spec[i] <- "soredium"
    }
    if(df[i,5] == "shockleyi" & df[i,3] < threshold){                          #if the Q matrix value for soredium is above the threshold...
      count_shock_admix <- count_shock_admix + 1
      count_total_admix <- count_total_admix + 1
    }
    if(df[i,5] == "shockleyi"){                           #if the Q matrix value for shockleyi is above the threshold...
      count_shock <- count_shock + 1
    }
    if(df[i,5] == "soredium"){                           #if the Q matrix value for shockleyi is above the threshold...
      count_sored <- count_sored + 1
    }
  }
  cat("\n")
  cat("Threshold: ", threshold, "\n")
  cat("\n")
  print("number of soredium in the assembly:")
  print(count_sored)
  print("Morphological soredium that show admixture:")
  cat(count_sored_admix, "or", count_sored_admix/count_sored, "\n")
  cat("\n")
  print("number of shockleyi in the assembly:") 
  print(count_shock)
  print("Morphological shockleyi that show admixture:")
  cat(count_shock_admix, "or", count_shock_admix/count_shock, "\n")
  cat("\n")
  print("Total that show admixture:")
  cat(count_total_admix, "or", count_total_admix/nrow(df), "\n")
  cat("\n")
  cat("\n")
}
Q_admix(Q, .89999) #call with .90 threshold.

## This function counts the number of samples with genome composition higher than 90%
Q_90 <- function(df, threshold){
  spec <- c()
  count_sored_thresh <- 0
  count_shock_thresh <- 0
  count_shock <- 0
  count_sored <- 0
  for(i in 1:nrow(df)){
    if(df[i,5] == "soredium" & df[i,2] > threshold){  
      count_sored_thresh <- count_sored_thresh + 1
      spec[i] <- "soredium"
    }
    if(df[i,5] == "shockleyi" & df[i,3] > threshold){                          #if the Q matrix value for soredium is above the threshold...
      count_shock_thresh <- count_shock_thresh + 1
    }
    if(df[i,5] == "shockleyi"){                           #if the Q matrix value for shockleyi is above the threshold...
      count_shock <- count_shock + 1
    }
    if(df[i,5] == "soredium"){                           #if the Q matrix value for shockleyi is above the threshold...
      count_sored <- count_sored + 1
    }
  }
  cat("\n")
  cat("Threshold: ", threshold, "\n")
  cat("\n")
  print("number of soredium in the assembly:")
  print(count_sored)
  print("Samples with genome compositions above threshold:")
  cat(count_sored_thresh, "or", count_sored_thresh/count_sored, "\n")
  cat("\n")
  print("number of shockleyi in the assembly:") 
  print(count_shock)
  print("Samples with genome compositions above threshold:")
  cat(count_shock_thresh, "or", count_shock_thresh/count_shock, "\n")
  cat("\n")
  cat("\n")
}
Q_90(Q, .89999) #call with .90 threshold.


## This function counts the number of discordant individuals.
Q_discord <- function(df){
  count_shock <- 0
  count_sored <- 0
  count_sored_wrong <- 0
  count_shock_wrong <- 0
  for(i in 1:nrow(df)){
    if(df[i,5] == "soredium" & df[i,2] < 0.5){  #if str says that this individual is soredium, and I specified that it is, call it soredium.
      count_sored_wrong <- count_sored_wrong + 1
    }
    if(df[i,5] == "shockleyi" & df[i,3] < 0.5){                          #if the Q matrix value for soredium is above the threshold...
      count_shock_wrong <- count_shock_wrong + 1
    }
    if(df[i,5] == "shockleyi"){                           #if the Q matrix value for shockleyi is above the threshold...
      count_shock <- count_shock + 1
    }
    if(df[i,5] == "soredium"){                           #if the Q matrix value for shockleyi is above the threshold...
      count_sored <- count_sored + 1
    }
  }
  cat("\n")
  print("Number of soredium in the assembly:")
  print(count_sored)
  print("Discordant Samples:")
  cat(count_sored_wrong, "or", count_sored_wrong/count_sored, "\n")
  cat("\n")
  print("Number of shockleyi in the assembly:") 
  print(count_shock)
  print("Discordant Samples:")
  cat(count_shock_wrong, "or", count_shock_wrong/count_shock, "\n")
  cat("\n")
  cat("\n")
}
Q_discord(Q) 
