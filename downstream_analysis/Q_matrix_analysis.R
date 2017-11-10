#Number (and %) of soredium coming out as soredium (and all four combinations)
#Q_threshold - determines the number of matching samples- considers Q matrix identitity true identity
#Q_threshold_morph - determines the number of matching samples- considers morphological identitity true identity
#Q_admix - determines the number and percentage of samples with admixture. Counts any sample with less than .90 genome composition as admixed (even if this is .00001)


setwd('/Users/jimblotter/Desktop/Grad_School/Data_Analysis/erisor/QC')
#Read in the data
Q <- read.csv("Q_mat_155.csv", header = TRUE, sep=",")
Q

#Set threshold at .90
Q_threshold <- function(df, threshold){
  spec <- c()
  count_sodso <- 0
  count_shdsh <- 0
  count_sored <- 0
  count_shock <- 0
  for(i in 1:nrow(df)){
    if(df[i,2] > threshold){                          #if the Q matrix value for soredium is above the threshold...
      count_sored <- count_sored + 1
    }
    if(df[i,2] > threshold & df[i,5] == "soredium"){  #if str says that this individual is soredium, and I specified that it is, call it soredium.
      count_sodso <- count_sodso + 1
      spec[i] <- "soredium"
    }
    if(df[i,3] > threshold){                           #if the Q matrix value for shockleyi is above the threshold...
      count_shock <- count_shock + 1
    }
    if(df[i,3] > threshold & df[i,5] == "shockleyi"){  #if str says that this individual is shockleyi, and I specified that it is, call it shockleyi.
      count_shdsh <- count_shdsh +1
    }
  }
  cat("\n")
  cat("Threshold: ", threshold, "\n")
  cat("\n")
  print("number of soredium passed the threshold:")
  print(count_sored)
  print("soredium that passed the threshold AND were correctly determined:")
  cat(count_sodso, "or", count_sodso/count_sored, "\n")
  cat("\n")
  print("number of shockleyi passed the threshold:") 
  print(count_shock)
  print("shockleyi that passed the threshold AND were correctly determined:")
  cat(count_shdsh, "or", count_shdsh/count_shock, "\n")
  cat("\n")
}

#Call with .90 and .95 as threshold, with and without GH/GHE.
Q_threshold(Q, .89999)
Q_threshold(Q, .94444)

##########################################################################
Q_threshold_morph <- function(df, threshold){
  spec <- c()
  count_sodso <- 0
  count_shdsh <- 0
  count_sored <- 0
  count_shock <- 0
  for(i in 1:nrow(df)){
    if(df[i,5] == "soredium"){  #if str says that this individual is soredium, and I specified that it is, call it soredium.
      count_sored <- count_sored + 1
      spec[i] <- "soredium"
    }
    if(df[i,2] > threshold & df[i,5] == "soredium"){                          #if the Q matrix value for soredium is above the threshold...
      count_sodso <- count_sodso + 1
    }
    if(df[i,5] == "shockleyi"){                           #if the Q matrix value for shockleyi is above the threshold...
      count_shock <- count_shock + 1
    }
    if(df[i,3] > threshold & df[i,5] == "shockleyi"){  #if str says that this individual is shockleyi, and I specified that it is, call it shockleyi.
      count_shdsh <- count_shdsh +1
    }
  }
  cat("\n")
  cat("Threshold: ", threshold, "\n")
  cat("\n")
  print("number of soredium in the assembly:")
  print(count_sored)
  print("Morphological soredium that were correctly determined:")
  cat(count_sodso, "or", count_sodso/count_sored, "\n")
  cat("\n")
  print("number of shockleyi in the assembly:") 
  print(count_shock)
  print("Morphological shockleyi that were correctly determined:")
  cat(count_shdsh, "or", count_shdsh/count_shock, "\n")
  cat("\n")
}

#Call with .90 and .95 as threshold, with and without GH/GHE.
Q_threshold_morph(Q, .89999)
Q_threshold_morph(Q, .94444)

##########################################################################

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
Q_admix(Q, .89999)
Q_admix(Q, .94444)
