library(readxl)
library(splitstackshape)
library(stringr)
library(caret)
library(igraph)
library(e1071)
library(pROC)

preprocessing <- function(data, proteins) {
  na_values <- apply(data[,proteins],2,function(x) sum(is.na(x)))
  set_good <- names(na_values[na_values == 0])
  set_middle <- names(na_values[na_values <= nrow(data)/10 & na_values>0])
  set_bad <- names(na_values[na_values > nrow(data)/10])
  
  for(p in set_middle){
    data[is.na(data[,p]), p] <- mean(data[,p], na.rm = TRUE)
  }
  if(length(set_bad) != 0) print(str_c("Proteins ",paste0(set_bad, collapse = ', '), " were excluded as persentage of missing values for them more than 10%."))
  if(length(set_middle) != 0) print(str_c("Missing values for ",paste0(set_middle, collapse = ', '), " were changed by mean."))
  if(length(set_good) != 0) print(str_c("Proteins ",paste0(set_good, collapse = ', '), " don't have missing values."))
  
  
  final_set_of_proteins <- c(set_good, set_middle)
  
  return(list(data, final_set_of_proteins))
}

