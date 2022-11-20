PATH = "/Users/Shared/Previously Relocated Items/Security/OC_synolitic_and_longit/" #PATH to main folder

source(paste0(PATH, "functions/general_functions.R")) #add functions
source(paste0(PATH, "functions/networks_functions.R")) #add functions of networks

# load data
load(paste0(PATH,"1_DATA/all_data.RData")) 
load(paste0(PATH,"1_DATA/all_data_long.RData"))

############## getting predictions "xgbTree", "nnet", "glmnet" applying to raw data or raw longit data (to indices of proteins)
models <- c("xgbTree", "nnet", "glmnet")
results_on_raw_data <- data.frame()

raw_data <- getting_predictions_on_test_for_raw_data(all_data, "Usual",c(1:18), 19, models)
raw_data_long <- getting_predictions_on_test_for_raw_data(all_data_long, "Long",c(1:18), 19, models)

############## getting predictions for best characteristics of networks


DIR_input <- paste0(PATH, "2_NETWORKS/")

best_characteristics <- c("closeness_mean","closeness_max","closeness_min", 
         "Eweights_max", "Eweights_mean",
         "strength_mean","strength_min","strength_max")

glm <- collect_predictions_on_test(all_data, "GLM", "Usual",DIR_input,c(1:18), 19, best_characteristics)
svm <- collect_predictions_on_test(all_data, "SVM", "Usual",DIR_input,c(1:18), 19, best_characteristics)

long_xgbtree <- collect_predictions_on_test(all_data_long, "LONG_xgbTree", "Long",DIR_input,c(1:18), 19, best_characteristics)
long_glmnet <- collect_predictions_on_test(all_data_long, "LONG_glmnet", "Long",DIR_input,c(1:18), 19, best_characteristics)
long_nnet <- collect_predictions_on_test(all_data_long, "LONG_nnet", "Long",DIR_input,c(1:18), 19, best_characteristics)  

##############
# Non-longitudinal algorithms (unlike longitudinal ones) contain predictions for samples at the first time point. 
# To make the results comparable, I removed these samples.
non_long_data <- all_data[[19]][[1]]

unique_patients <- unique(non_long_data$main_id) 
exclude_samples <- c()
p <- unique_patients[1]
for (p in unique_patients) {
  temp <- subset(non_long_data,non_long_data$main_id == p) #get patient
  temp <-  temp[order(temp$Age.at.sample.taken..years),]  #order by age
  exclude_samples <- c(exclude_samples, temp[1,]$id)
}

glm$excluded <- unlist(lapply(glm$id, function(x)  length(intersect(exclude_samples,x))))
svm$excluded <- unlist(lapply(svm$id, function(x)  length(intersect(exclude_samples,x))))
raw_data$excluded <- unlist(lapply(raw_data$id, function(x)  length(intersect(exclude_samples,x))))
raw_data_long$excluded <- unlist(lapply(raw_data_long$id, function(x)  length(intersect(exclude_samples,x))))

glm <- subset(glm,glm$excluded == 0)
svm <- subset(svm,svm$excluded == 0)
raw_data <- subset(raw_data,raw_data$excluded == 0)
raw_data_long <- subset(raw_data_long,raw_data_long$excluded == 0)

############## calculate AUC for each method, for each characteristic, for loocv and for main
raw_data_aucs <- auc_calculation(raw_data, "main", models)
raw_data_long_aucs <- auc_calculation(raw_data_long, "main", models)
raw_data_aucs$type_of_approach <- str_c("Raw_Usual_", raw_data_aucs$characteristic)
raw_data_long_aucs$type_of_approach <- str_c("Raw_Long_", raw_data_long_aucs$characteristic)
glm_aucs <- auc_calculation(glm, "main", best_characteristics)
svm_aucs <- auc_calculation(svm, "main", best_characteristics)
long_xgbtree_aucs <- auc_calculation(long_xgbtree, "main", best_characteristics)
long_glmnet_aucs <- auc_calculation(long_glmnet, "main", best_characteristics)
long_nnet_aucs <- auc_calculation(long_nnet, "main", best_characteristics)

all_aucs <- rbind(raw_data_aucs, raw_data_long_aucs, glm_aucs, svm_aucs, long_xgbtree_aucs, long_glmnet_aucs, long_nnet_aucs)
unique(all_aucs$type_of_approach)
fig1 <- ggplot(all_aucs, aes(factor(type_of_approach, levels=unique(all_aucs$type_of_approach)), auc)) + geom_boxplot() + facet_grid(.~type_of_experement) + 
  theme(axis.text.x = element_text(angle=90, hjust=1)) + labs(x="approach")
 
loocv <- subset(all_aucs,all_aucs$type_of_experement == "loocv")
##### best result on loocv
best_loocv <- subset(loocv, loocv$auc == max(loocv$auc))
##### result on test
subset(all_aucs, all_aucs$characteristic == best_loocv$characteristic & all_aucs$type_of_approach == best_loocv$type_of_approach & 
         all_aucs$type_of_experement == "main")
##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### #####

a<- collect_predictions_on_edge_CA125_HE4(all_data_long, "LONG_xgbTree", "Long",DIR_input,c(1:18), 19, best_characteristics)

print("only edge of CA125 and HE4, Long, XgbTree:")
print(a)
print("Best Network result, Long, XgbTree:")
print(str_c("loocv: ",best_loocv$auc, "; main: ",subset(all_aucs, all_aucs$characteristic == best_loocv$characteristic & all_aucs$type_of_approach == best_loocv$type_of_approach & 
                                                          all_aucs$type_of_experement == "main")$auc ))
