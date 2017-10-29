#Number (and %) of soredium coming out as soredium (and all four combinations)
#Counting .90 and above as "shockleyi" and .95 and above as shockleyi
#With and without GH/GHE.

#Set threshold at .90
Q_threshold <- function(df, threshold){
  spec <- c()
  count_sodso <- 0
  count_shdsh <- 0
  count_sored <- 0
  count_shock <- 0
  for(i in 1:nrow(df)){
    if(df[i,2] > threshold){
      count_sored <- count_sored + 1
    }
    if(df[i,2] > threshold & df[i,5] == "soredium"){  #if str says that this individual is soredium, and I specified that it is, call it soredium.
      count_sodso <- count_sodso + 1
      spec[i] <- "soredium"
    }
    if(df[i,2] < 1-threshold){
      count_shock <- count_shock + 1
    }
    if(df[i,2] < 1-threshold & df[i,5] == "shockleyi"){
      count_shdsh <- count_shdsh +1
    }
  }
  cat("\n")
  cat("Threshold: ", threshold, "\n")
  cat("\n")
  print("number of soredium passed the threshold:")
  print(count_sored)
  print("soredium correctly determined:")
  cat(count_sodso, "or", count_sodso/count_sored, "\n")
  cat("\n")
  print("number of shockleyi passed the threshold:") 
  print(count_shock)
  print("shockleyi correctly determined:")
  cat(count_shdsh, "or", count_shdsh/count_shock)
}

##########################################################################
#Read in the data
Q <- read.csv("Q_matrix.csv", header = TRUE, sep=",")
#Q
#Read in without GH and GHE
Q2 <- read.csv("Q_matrix_GHandE_removed.csv", header = TRUE, sep=",")
#Q2

#Call with .90 and .95 as threshold, with and without GH/GHE.
Q_threshold(Q, .89999)
Q_threshold(Q, .94444)

Q_threshold(Q2, .89999)
Q_threshold(Q2, .94444)


