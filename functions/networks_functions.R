
#function calculation_PN(list_of_data, type_PN, type_data, DIR)
# create DIR in 2_NETWORKS which will called type_data_type_PN
# where for each dataframe from list_of_data will be created 2 files: chars and networks

calculation_PN <- function(all_data, type_PN, type_data, DIR) {
  
  subDir <- paste0(DIR,type_data,"_", type_PN)
  if (file.exists(subDir)){
    setwd(subDir)
  } else {
    dir.create(subDir)
    setwd(subDir)
  }
  
  for(i in 1:length(all_data)) {
    print("================================")
    print(i)
    print("================================")
    
    data_object <- all_data[[i]]
    
    data = data_object[[1]]
    proteins = data_object[[2]]
    name_data = data_object[[3]]
    
    list_of_pair <- get_list_of_pair(proteins)
    parametr <- c()
    
    edges_dataframe <- Main_Calculation(parametr, "temp", data, list_of_pair, type_PN) # calculated network
    
    chars <- get_char(data, edges_dataframe, proteins, 0.5) # calculated characteristics of network
    save(chars, file = str_c(type_data,"_",name_data,"_chars.RData")) #save characteristics of network
    
    lst <- list(edges_dataframe, data[,c("score", "id", "label")])
    save(lst, file = str_c(type_data,"_",name_data,"_networks.RData")) #save network
  }
}




#### Main function:
## Return networks
## -- If type_of_PN == ZSCORE - wieght between 2 proteins - zscore from linear regression on controls (Zanin)
## -- If type_of_PN == KDE - wieght between 2 proteins - "distance" from kde on controls (Harry Whitwell, OC)
## -- If type_of_PN == SVM/GLM - wieght between 2 proteins - probability (obtained by the SVM/GLM model on 2 proteins) of belonging to a class of cases
## -- If type_of_PN == LONG_xgbTree(/glmnet/nnet)" - wieght between 2 proteins - probability (obtained by the xgbTree/glmnet/nnet model on 8 parameters: A, B, C, D indices for 2 proteins) of belonging to a class of cases

Main_Calculation <- function(parametr, filename,data1, list_of_pair, type_of_PN) {
  if(type_of_PN == "ZSCORE") {
    st <- Sys.time()
    edges_dataframe <- lapply(list_of_pair, function(x) 
      results_calculation_zscore(x[1], x[2], data1,c(parametr)))
    print(Sys.time() - st)
  }
  
  if(type_of_PN == "KDE") {
    st <- Sys.time()
    edges_dataframe <- lapply(list_of_pair, function(x) 
      results_calculation_kde(x[1], x[2], data1,c(parametr)))
    print(Sys.time() - st)
  }
  
  if(type_of_PN == "SVM") {
    st <- Sys.time()
    edges_dataframe <- lapply(list_of_pair, function(x) 
      results_calculation_svm(x[1], x[2], data1,c(parametr)))
    print(Sys.time() - st)
  }
  
  if(type_of_PN == "GLM") {
    st <- Sys.time()
    edges_dataframe <- lapply(list_of_pair, function(x) 
      results_calculation_glm(x[1], x[2], data1,c(parametr)))
    print(Sys.time() - st)
  }
  
  if(type_of_PN == "LONG_xgbTree") {
    st <- Sys.time()
    edges_dataframe <- lapply(list_of_pair, function(x) 
      results_calculation_long_xgbTree(x[1], x[2], data1,c(parametr)))
    print(Sys.time() - st)
  }
  if(type_of_PN == "LONG_glmnet") {
    st <- Sys.time()
    edges_dataframe <- lapply(list_of_pair, function(x) 
      results_calculation_long_glmnet(x[1], x[2], data1,c(parametr)))
    print(Sys.time() - st)
  }
  if(type_of_PN == "LONG_nnet") {
    st <- Sys.time()
    edges_dataframe <- lapply(list_of_pair, function(x) 
      results_calculation_long_nnet(x[1], x[2], data1,c(parametr)))
    print(Sys.time() - st)
  }

  edges_dataframe <- do.call(rbind, edges_dataframe)
  edges_dataframe <- data.frame(edges_dataframe)
  colnames(edges_dataframe) <- c("p1","p2",data1$id)
  edges_dataframe[, 3:(ncol(edges_dataframe))] <- 
    sapply(edges_dataframe[, 3:ncol(edges_dataframe)], as.numeric)
  return(edges_dataframe)
}




get_list_of_pair <- function(proteins) {
  list_of_pair <- list()
  k=0
  for(i in 1:(length(proteins)-1)) {
    for(j in (i+1):length(proteins)) {
      k = k + 1
      list_of_pair[[k]] <- c(proteins[i], proteins[j])
    }
  }
  return(list_of_pair)
}



########## Parenclitic ZSCORE
results_calculation_zscore <- function(p1,p2, data,dop_param) {
  st <- Sys.time()
  df <- data[c(p1, p2,dop_param, "score", "label")]
  rownames(df) <- data$id
  answ <- kernel_parenclitic_zscore(df)
  return(c(p1, p2,answ))
}

kernel_parenclitic_zscore <- function(df) {
  set.seed(1)
  df_train_0 <- subset(df,df$score == 0 & df$label == "TRAIN")
  df_train_1 <- subset(df,df$score == 1 & df$label == "TRAIN")
  
  colnames(df_train_0) <- c("p1", "p2", "score", "label")
  m <- lm(p1~p2, data = df_train_0)
  colnames(df) <- c("p1", "p2", "score", "label")
  df$prediction <- predict(m, df)
  df$resud <- (df$prediction-df$p1)
  
  df_train_0 <- subset(df,df$score == 0 & df$label == "TRAIN")
  Mean <- mean(df_train_0$resud)
  SD <- sd(df_train_0$resud)
  df$resud2 <- abs(df$resud - Mean)/SD
  return(df$resud2)
}

########## Parenclitic KDE
results_calculation_kde <- function(p1,p2, data,dop_param) {
  st <- Sys.time()
  df <- data[c(p1, p2,dop_param, "score", "label")]
  rownames(df) <- data$id
  #df_train <- subset(df,df$label == "TRAIN")[1:(ncol(df)-1)]
  answ <- kernel_parenclitic_kde(df)
  return(c(p1, p2,answ))
}


kernel_parenclitic_kde <- function(df) {
  set.seed(1)
  df_train_0 <- subset(df,df$score == 0 & df$label == "TRAIN")
  dens <- try(MASS::kde2d(df_train_0[,1],df_train_0[,2]))
  if(class(dens)  != "try-error") {
    gr <- data.frame(with(dens, expand.grid(x,y)), as.vector(dens$z))
    names(gr) <- c("xgr", "ygr", "zgr")
    mod <- loess(zgr~xgr*ygr, data=gr)
    
    df$pointdens <- predict(mod, newdata=data.frame(xgr=df[,1], ygr=df[,2]))
    gr$preds <- predict(mod, newdata=data.frame(xgr=gr[,1], ygr=gr[,2]))
    df$pointdens <- df$pointdens/sum(gr$preds)
    gr$preds <- gr$preds/sum(gr$preds)
    
    center <- subset(gr, gr$preds == max(gr$preds,na.rm=T))
    gr$dist_to_center <- sqrt((gr$xgr - as.numeric(center[1,1]))**2 + (gr$ygr - as.numeric(center[1,2]))**2)
    
    
    gr_countour <- subset(gr,gr$xgr == max(gr$xgr) | gr$ygr == max(gr$ygr) | gr$xgr == min(gr$xgr) | gr$ygr == min(gr$ygr))
    v <-  gr$preds
    df$pointdens3 <- unlist(lapply(df$pointdens, function(x) 1-sum(v[v<=x], na.rm=T)))
    
    df$pointdens2 <- apply(df,1,function(x) ifelse(is.na(x[6]) == T, distance_extension(x, gr_countour,center), as.numeric(x[6])))
  } else {
    df$pointdens2 = 0
  }
  return(df$pointdens2)
}

distance_extension <- function(x, gr_countour,center) {
  return(1)
}

########## Synolitic SVM
results_calculation_svm <- function(p1,p2, data,dop_param) {
  st <- Sys.time()
  df <- data[c(p1, p2,dop_param, "score", "label")]
  rownames(df) <- data$id
  df_train <- subset(df,df$label == "TRAIN")[1:(ncol(df)-1)]
  svmfit <- kernel_parenclitic_prob_cv_fast(df_train)
  print(Sys.time() - st)
  answ <- as.numeric(attr(predict(svmfit, df[1:(ncol(df)-1)], probability = T), 'probabilities')[,2])
  return(c(p1, p2,answ))
}

kernel_parenclitic_prob_cv_fast <- function(df_train) {
  set.seed(123)
  #print(df_train$score)
  df_train$score <- factor(df_train$score, levels=c(0,1))
  svmFit <- svm(score ~ ., data = df_train, probability = T, scale=T)
  return(svmFit)
}

########## Synolitic GLM
results_calculation_glm <- function(p1,p2, data,dop_param) {
  st <- Sys.time()
  df <- data[c(p1, p2,dop_param, "score", "label")]
  rownames(df) <- data$id
  df_train <- subset(df,df$label == "TRAIN")[1:(ncol(df)-1)]
  glmfit <- kernel_parenclitic_prob_glm(df_train)
  print(Sys.time() - st)
  answ <- predict(glmfit, df[1:(ncol(df)-1)], type="response")
  return(c(p1, p2,answ))
}

kernel_parenclitic_prob_glm <- function(df_train) {
  set.seed(123)
  #df_train$score <- as.factor(df_train$score)
  glmfit <- glm(score ~ ., data = df_train, family = binomial)
  return(glmfit)
}

predicted_models <- function(df_train, type) {
  set.seed(123)
  df_train$score <- as.factor(df_train$score)
  levels(df_train$score) <- c("first_class","second_class")
  modelfit <- try(train(score ~ .,
                        data = df_train,
                        method = type,
                        preProc = c("center", "scale"),
                        metric = "ROC",
                        trControl = trainControl(method = "cv",number=5,
                                                 classProbs = TRUE, 
                                                 summaryFunction = twoClassSummary,
                                                 allowParallel = T)))
  print(.Last.value)
  return(modelfit)
}


getting_predictions_on_test_for_raw_data <- function(all_data, type_of_data,loocv_range, num_of_main, models) {
  results <- data.frame()
  for(i in loocv_range) {
    data <- all_data[[i]][[1]]
    proteins <- all_data[[i]][[2]]
    temp <- getting_prediction_on_raw_data_for_test(data, proteins, type_of_data, models)
    temp$run <- i
    results <- rbind(results, temp)
  }
  data <- all_data[[num_of_main]][[1]]
  proteins <- all_data[[num_of_main]][[2]]
  temp <- getting_prediction_on_raw_data_for_test(data, proteins, type_of_data, models)
  temp$run <- "main"
  results <- rbind(results, temp)
  
  return(results)
}

getting_prediction_on_raw_data_for_test <- function(data, proteins, type_of_data, models) {
  results = data.frame()
  
  if (type_of_data == "Long") {
    proteins <-  c(str_c(proteins,"_A_index"), str_c(proteins,"_B_index"),str_c(proteins,"_C_index"),str_c(proteins,"_D_index"))
  }
  
  df <- data[c(proteins, "score", "label")]
  is.na(df)<-sapply(df, is.infinite)
  df[is.na(df)]<-0
  
  rownames(df) <- data$id
  df_train <- subset(df,df$label == "TRAIN")[1:(ncol(df)-1)]
  
  results <- data.frame()
  for(m in models) {
    M <- predicted_models(df_train, m)
    df$prediction <- predict(M, df, type="prob")[,2]
    temp_test <- subset(df, df$label == "TEST")
    
    train_roc <- pROC::roc(subset(df, df$label == "TRAIN")$score, subset(df, df$label == "TRAIN")$prediction )
    
    results <- rbind(results, data.frame(train_directions = train_roc$direction, prediction = temp_test$prediction, 
                                         score = temp_test$score, id = rownames(temp_test),
                                         characteristic = m, type_of_approach = str_c(type_of_data,"_", m)))
  }
return(results)
}
########## Synolitic Longitudinal 

##### Synolitic Longitudinal xgbTree
results_calculation_long_xgbTree <- function(p1,p2, data,dop_param) {
  st <- Sys.time()
  params <- c(str_c(p1,"_A_index"), str_c(p1,"_B_index"),str_c(p1,"_C_index"),str_c(p1,"_D_index"),
              str_c(p2,"_A_index"), str_c(p2,"_B_index"),str_c(p2,"_C_index"),str_c(p2,"_D_index"))
  df <- data[c(params,dop_param, "score", "label")]
  is.na(df)<-sapply(df, is.infinite)
  df[is.na(df)]<-0
  
  rownames(df) <- data$id
  df_train <- subset(df,df$label == "TRAIN")[1:(ncol(df)-1)]
  
  M <- predicted_models(df_train, "xgbTree")
  answ <- predict(M, df, type="prob")[,2]
  return(c(p1, p2,answ))
}

##### Synolitic Longitudinal glmnet
results_calculation_long_glmnet <- function(p1,p2, data,dop_param) {
  st <- Sys.time()
  params <- c(str_c(p1,"_A_index"), str_c(p1,"_B_index"),str_c(p1,"_C_index"),str_c(p1,"_D_index"),
              str_c(p2,"_A_index"), str_c(p2,"_B_index"),str_c(p2,"_C_index"),str_c(p2,"_D_index"))
  df <- data[c(params,dop_param, "score", "label")]
  is.na(df)<-sapply(df, is.infinite)
  df[is.na(df)]<-0
  
  rownames(df) <- data$id
  df_train <- subset(df,df$label == "TRAIN")[1:(ncol(df)-1)]
  
  M <- predicted_models(df_train, "glmnet")
  answ <- predict(M, df, type="prob")[,2]
  return(c(p1, p2,answ))
}

##### Synolitic Longitudinal nnet
results_calculation_long_nnet <- function(p1,p2, data,dop_param) {
  st <- Sys.time()
  params <- c(str_c(p1,"_A_index"), str_c(p1,"_B_index"),str_c(p1,"_C_index"),str_c(p1,"_D_index"),
              str_c(p2,"_A_index"), str_c(p2,"_B_index"),str_c(p2,"_C_index"),str_c(p2,"_D_index"))
  df <- data[c(params,dop_param, "score", "label")]
  is.na(df)<-sapply(df, is.infinite)
  df[is.na(df)]<-0
  
  rownames(df) <- data$id
  df_train <- subset(df,df$label == "TRAIN")[1:(ncol(df)-1)]
  
  M <- predicted_models(df_train, "nnet")
  answ <- predict(M, df, type="prob")[,2]
  return(c(p1, p2,answ))
}

#### Characteristics
get_char <- function(data1, edges_dataframe, proteins,TR) {
  graph_objects <- lapply(as.character(data1$id), function(x) graph_objects_getting(x,edges_dataframe,proteins))
  data_covid <- data1
  important_nodes <- c()
  st <- Sys.time()
  res_characteristics <- data.frame()
  for(j in 1: length(graph_objects)) {
    gr <- graph_objects[[j]]
    res_characteristics <- rbind(res_characteristics, 
                                 data.frame(t(getting_characterictics(gr,important_nodes))))
  }
  print(Sys.time() - st)
  
  res_characteristics$score <- data_covid$score
  res_characteristics$id <- data_covid$id
  res_characteristics$label <- data_covid$label
  
  return(res_characteristics)
}

getting_characterictics<- function(gr, important_nodes) {
  closeness = as.vector(closeness(gr,weights=E(gr)$width))
  closeness_char <- main_stat_char(closeness)
  names(closeness_char) <- str_c("closeness_", names(closeness_char))
  
  betweenness = as.vector(betweenness(gr,weights=E(gr)$width))
  betweenness_char <- main_stat_char(betweenness)
  names(betweenness_char) <- str_c("betweenness_", names(betweenness_char))
  
  edge_betweenness = as.vector(edge_betweenness(gr,weights=E(gr)$width))
  edge_betweenness_char <- main_stat_char(edge_betweenness)
  names(edge_betweenness_char) <- str_c("edge_betweenness_", names(edge_betweenness_char))
  
  page.rank = as.vector(page.rank(gr,weights=E(gr)$width)$vector)
  page.rank_char <- main_stat_char(page.rank)
  names(page.rank_char) <- str_c("page.rank_", names(page.rank_char))
  
  eigen_centrality = as.vector(eigen_centrality(gr,weights=E(gr)$width)$vector)
  eigen_centrality_char <- main_stat_char(eigen_centrality)
  names(eigen_centrality_char) <- str_c("eigen_centrality_", names(eigen_centrality_char))
  
  authority_score = as.vector(authority_score(gr,weights=E(gr)$width)$vector)
  authority_score_char <- main_stat_char(authority_score)
  names(authority_score_char) <- str_c("authority_score_", names(authority_score_char))
  
  s <- as.vector(strength(gr,weights=E(gr)$width))
  s_char <- main_stat_char(s)
  names(s_char) <- str_c("strength_", names(s_char))
  
  Eweghts <- E(gr)$width
  Eweghts_char <- main_stat_char(Eweghts)
  names(Eweghts_char) <- str_c("Eweights_", names(Eweghts_char))
  
  return(c(closeness_char, betweenness_char, edge_betweenness_char, 
           page.rank_char, eigen_centrality_char, authority_score_char,s_char, Eweghts_char))
  
}

main_stat_char <- function(x) {
  v <- c(sum(x==0),min(x), max(x), mean(x), sd(x))
  coefvar <- ifelse(mean(x) != 0, sd(x)/mean(x), 0)
  v <- c(v, coefvar)
  names(v) <- c("zeros","min", "max", "mean", "sd", "coefvar")
  return(v)
}
graph_objects_getting <- function(id, very_good, proteins) {
  df <- data.frame(p1=very_good$p1, p2=very_good$p2,id=very_good[,id])
  df$id2 <- ifelse(is.na(df$id), 0, df$id)
  df_good <- df
  from1 <- as.vector(df_good$p1)
  to1 <- as.vector(df_good$p2)
  relations <- data.frame(from=from1,
                          to=to1)
  g <- graph.data.frame(relations, directed=F, vertices=proteins)
  
  
  E(g)$width = ifelse(df_good$id2 <= 0.00001, 0.00001,df_good$id2)
  return(g)
}

######### getting predictions on characteristics

get_chars_predictions_for_tests <- function(chars_set, name_of_characteristics) {

  results <- data.frame()
  is.na(chars_set)<-sapply(chars_set, is.infinite)
  chars_set[is.na(chars_set)]<-0
  
  chars_set_train_data <- subset(chars_set,chars_set$label == "TRAIN")
  chars_set_test_data <- subset(chars_set,chars_set$label == "TEST")
  
  for(i in 1:length(name_of_characteristics)) {
    temp_train <- chars_set_train_data[,c(name_of_characteristics[i],"score")]
    temp_test <- chars_set_test_data[,c(name_of_characteristics[i],"score")]
    m <- try(glm(factor(score)~.,data = temp_train, family = "binomial"))
    
    #if(class(m)  != "try-error") {
      ress_train <- predict(m, temp_train, type="response")
      ress_test <- predict(m, temp_test, type="response")
      
      train_roc <- pROC::roc(temp_train$score, ress_train )
      results <- rbind(results, data.frame(train_directions = train_roc$direction, prediction = ress_test, 
                                   score = temp_test$score, id = chars_set_test_data$id,
                                   characteristic = name_of_characteristics[i]))
   # }
  }
  return(results)
}

collect_predictions_on_edge_CA125_HE4 <- function(all_data, type_PN, type_data, DIR, loocv_range, num_of_main, characteristics) {
  subDir <- paste0(DIR,type_data,"_", type_PN,"/")
  all_predicitons_loocv = data.frame()
  for(i in loocv_range) {
    data <- all_data[[i]][[1]]
    name_of_set <- all_data[[i]][[3]]
    load(paste0(subDir,type_data,"_",name_of_set,"_networks.RData")) 
    data$predictions <- as.vector(t(lst[[1]][6,data$id]))
    all_predicitons_loocv <- rbind(all_predicitons_loocv, subset(data,data$label == "TEST")[,c("score", "predictions")])
  }
  
  data_main <- all_data[[num_of_main]][[1]]
  name_of_set <- all_data[[num_of_main]][[3]]
  load(paste0(subDir,type_data,"_",name_of_set,"_networks.RData")) 
  data_main$predictions <- as.vector(t(lst[[1]][6,data_main$id]))
  all_predicitons_main <- subset(data_main,data_main$label == "TEST")[,c("score", "predictions")]
  
  auc_main <- pROC::roc(all_predicitons_main$score, all_predicitons_main$predictions)$auc
  auc_loocv <- pROC::roc(all_predicitons_loocv$score, all_predicitons_loocv$predictions)$auc
  return(str_c("loocv: ", auc_loocv,"; main: ",auc_main))
}

collect_predictions_on_test <- function(all_data, type_PN, type_data, DIR, loocv_range, num_of_main, characteristics) {
  subDir <- paste0(DIR,type_data,"_", type_PN,"/")
  all_predicitons = data.frame()
  for(i in loocv_range) {
    name_of_set <- all_data[[i]][[3]]
    load(paste0(subDir,type_data,"_",name_of_set,"_chars.RData"))    
    predictions <- get_chars_predictions_for_tests(chars, characteristics)
    predictions$run = i
    all_predicitons <- rbind(all_predicitons, predictions)
  }

  name_of_main_set <- all_data[[num_of_main]][[3]]
  load(paste0(subDir,type_data,"_",name_of_main_set,"_chars.RData"))    
  predictions_main <- get_chars_predictions_for_tests(chars, characteristics)
  predictions_main$run = "main"
  all_predicitons <- rbind(all_predicitons, predictions_main)
  all_predicitons$type_of_approach <- str_c(type_data,"_", type_PN)
  return(all_predicitons)
}

auc_calculation<- function(set_of_predictions, name_of_main_run, characteristics) {
  auc_results <- data.frame()
  main <- subset(set_of_predictions,set_of_predictions$run == name_of_main_run)
  loocv <- subset(set_of_predictions,set_of_predictions$run != name_of_main_run)
  for(ch in characteristics) {
    temp_main <- subset(main,main$characteristic == ch)
    temp_loocv <- subset(loocv,loocv$characteristic == ch)
    
    auc_main <- pROC::roc(temp_main$score, temp_main$prediction)$auc
    auc_loocv <- pROC::roc(temp_loocv$score, temp_loocv$prediction)$auc
    auc_results <- rbind(auc_results, data.frame(type_of_approach = temp_main[1,]$type_of_approach, 
                              characteristic = ch,
                              auc = auc_main,
                              type_of_experement = "main"),
                         data.frame(type_of_approach = temp_main[1,]$type_of_approach, 
                                    characteristic = ch,
                                    auc = auc_loocv,
                                    type_of_experement = "loocv"))
  }
  return(auc_results)
}
  

