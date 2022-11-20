A_function <- function(Y, t) {
  deltaY <- c()
  deltat <- c()
  for (i in 1:(length(Y)-1)) {
    deltaY[i] <- Y[i+1] - Y[i]
    deltat[i] <- t[i+1] - t[i]
  }
  return(sum((deltaY*deltat)/2)/(length(t)-1))
}

B_function <- function(Y, t) {
  return(sqrt(sum((Y - mean(Y))**2)/length(t))/mean(Y))
}

C_function <- function(Y, t) {
  return(sum(Y*t)/sum(t))
}

derivative_function <- function(Y, t) {
  w <- c()
  deltaY <- c()
  deltat <- c()
  for (i in 1:(length(Y)-1)) {
    deltaY[i] <- Y[i+1] - Y[i]
    deltat[i] <- t[i+1] - t[i]
    w[i] <- 1/( t[length(t)] - (t[i+1]+t[i])/2)
  }
  return(sum(w*(deltaY/deltat)))
}

index_function <- function(Y, t, name) {
  if (name == "derivative") return(derivative_function(Y, t))
  if (name == "A") return(A_function(Y, t))
  if (name == "B") return(B_function(Y, t))
  if (name == "C") return(C_function(Y, t))
}

new_record_with_indices <- function(new_temp,proteins) {
  rownames(new_temp) <- new_temp$id
  new_id <- str_c(unique(new_temp$main_id),"_",nrow(new_temp))
  new_samples_names <- data.frame(id = new_id,
                                  score = new_temp[nrow(new_temp), ]$score,
                                  label = new_temp[nrow(new_temp), ]$label,
                                  main_id = unique(new_temp$main_id),
                                  number_longit = nrow(new_temp),
                                  range_from_first = new_temp[nrow(new_temp), ]$time -
                                    new_temp[1, ]$time,
                                  range_from_prev = new_temp[nrow(new_temp), ]$time -
                                    new_temp[(nrow(new_temp) -1), ]$time,
                                  time_group = new_temp[nrow(new_temp), ]$time_group)
  p <- proteins[1]
  for(p in proteins) {
    A_index <- (A_function(new_temp[,p], new_temp$time))
    B_index <- (B_function(new_temp[,p], new_temp$time))
    C_index <- (C_function(new_temp[,p], new_temp$time))
    D_index <- (derivative_function(new_temp[,p], new_temp$time))
    temp_index <- data.frame(A_index, B_index,C_index,D_index)
    colnames(temp_index)  <- str_c(p,"_",colnames(temp_index))
    new_samples_names <- cbind(new_samples_names, temp_index)
  }
  return(new_samples_names)
}
