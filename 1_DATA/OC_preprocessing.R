PATH = "/Users/Shared/Previously Relocated Items/Security/OC_synolitic_and_longit/" #PATH to main folder

source(paste0(PATH, "functions/general_functions.R")) #add functions
source(paste0(PATH, "functions/index_functions.R")) #add functions on indices

data <- read_excel(paste0(PATH,"1_DATA/AlldataTypeIIOC2.xlsx"), na = c('', 'NA', 'NA')) #read dataset
data <- data.frame(data)
data <- subset(data, data$QC.warning == "Pass") #exluded samples which didn't pass quality control


##################################################
##### PROTEINS
# excluded <- c("PROZ", "SHBG",	"CalProt",	"FN",	"GRP78",	"LCAT",	"CRP",	"IGFBP2") - proteins, which we initially excluded, but maybe we can use them now. Please ask Oleg
# inhouse<- colnames(data)[27:36] - "inhouse" proteins
# OLink <- colnames(data)[46:137] - "OLink" proteins
# CA125 (form "inhouse") have a twin in "OLink": MUC16
# HE4 (form "inhouse") have a twin in "OLink": WFDC2
# I've tried different combinations of these protein sets, but the best result I've got for "inhouse" set. I'll use it here for example

proteins <- colnames(data)[27:36]

preproc <- preprocessing(data, proteins) #you can find this function in functions/general_functions.R. Excluded proteins with missing values > 10%; missing values for others changed by mean
data <- preproc[[1]]
proteins <- preproc[[2]]
# Please not LRG1, PAEP were excluded!

data[,proteins] <- scale(data[,proteins]) #scale proteins of interest

##################################################
##### STRTIFICATION
# The idea was to divide the patients into TRAIN and TEST (without taking into account the number of measurements they have) by age (in the LAST time point), BMI and cancer status.

unique_patients <- unique(data$Rand.Vol.ID) 
strata_data <- data.frame()
for (id in unique_patients) {
  temp <- subset(data,data$Rand.Vol.ID == id) #get 1 patient
  temp <-  temp[order(temp$Age.at.sample.taken..years,decreasing = T),]  #order by age
  strata_data <- rbind(strata_data, temp[1,c("Rand.Vol.ID","Age.at.sample.taken..years", "BMI", "Cancer.Control.Status.1")]) #get age, bmi and status from LAST record (LAST time point)
}
rownames(strata_data) <- strata_data$Rand.Vol.ID

strata_data$q_age <- ifelse(strata_data$Age.at.sample.taken..years < median(strata_data$Age.at.sample.taken..years),0,1) #2 groups by age (through median)
strata_data$q_bmi<- ifelse(strata_data$BMI < median(strata_data$BMI),0,1) #2 groups by BMI (through median)

SEED = 2 #hust seed, can be changed
set.seed(SEED)

out <- stratified(strata_data, c("Cancer.Control.Status.1", "q_age","q_bmi"), size =0.6,bothSets=T) #size 0.6 allows you to get equal number of cases and controls in both subsets
train <- data.frame(out$SAMP1)
test <- data.frame(out$SAMP2)

train$label <- "TRAIN"
test$label <- "TEST"

strata_data_new <- rbind(train, test)
data <- merge(data, strata_data_new[,c("Rand.Vol.ID", "label")], by="Rand.Vol.ID") #add TRIAN/TEST labels to the main dataset

##################################################
#For the uniformity of the application of algorithms, the dataset should contain the parameters of the 
# - main_id: patient 
# - id: sample id 
# - label: TRAIN/TEST 
# - score: cases/controls (in OC dataset Cancer.Control.Status.1)
# - time: for longitudinal algorithms mainly (in OC dataset, this is the age at the time of sampling: Age.at.sample.taken..years)
# - time_group: time to diagnosis for cases and separate label for controls (in OC dataset Time.to.Dx.group..years)
data$main_id <- data$Rand.Vol.ID
data$id <- data$PROMISE.ID
data$label <- data$label
data$score <- data$Cancer.Control.Status.1
data$time <- data$Age.at.sample.taken..years
data$time_group <- data$Time.to.Dx.group..years
##################################################
##### Longitudinal prepаration
# I've tried different approaches, but the best was this one:
# - for each patient, for each of it's protein calculate A, B, C, D indices (please look at ptx file, slide 1 and functions/index_functions.R). 
#   Please note, we considering here all of possible lengths. If patient has a mesurment of protein p1 in time: p1_1, p1_2, p1_3, p1_4, when we will get three A (for example) indices:
#   A(p1_1,p1_2), A(p1_1,p1_2,p1_3), A(p1_1,p1_2,p1_3,p1_4)
# - In the kernel of synolitic approach instead of building svm or glm model on pair of proteins p1 and p2 we will use all 4 indices (A, B, C, D) for both proteins
#   Since the number of parameters in the kernel of the longitudinal approach is 8, and not 2 (as in the case of the usual sinolytic approach where we are using glm or svm), 
#   more complex models will be used inside this kernel, allowing automatic feature selection

# Here a dataset is created, in which indices are calculated for the original proteins.

data_long  <- data.frame()

patient_id <- unique_patients[1]
times <- 2
for(patient_id in unique_patients) {
  patient <- subset(data, data$main_id == patient_id)
  patient <-  patient[order(patient$time),]
  if(nrow(patient) >=2) {
    for(times in 1:(nrow(patient)-1)) { #choosing all possible lengths of longitudinal vectors
      temp <- patient[1:(times+1),]
      new_record <- new_record_with_indices(temp, proteins) #original dataset should have main_id, id, label, score, time, time_group
      data_long <- rbind(data_long, new_record) 
    }
  }
}

##################################################
##### LOOCV preparation
# Since the dataset is small, the loocv (namely, l2ocv) procedure on the TRIAN set was used to determine the best model.
# From the TRIAN set (with 36 patients, 18 cases and 18 controls), 18 datasets were made, in each of which a random pair of 1 case and 1 control was named as a TEST. 
# For each dataset all models were trained on 34 others patients and then the result was applied to the excluded pair. 
# At the end all these results were combined and some kind of metric was calculated (for example, AUC).
#
# For convenience of calculation by such a procedure here were created
# - a list of random pairs on the TRIAN set is created
# - for usual and for logitudinal approaches separatly: 
# --- a list of 18 datasets is created iteratively
# --- at the end (as 19 dataset) the main dataset is added (where labels are the main labels of the TRAIN and TEST)
#
# Later this 2 lists will be used in all procedures.


# getting id for cases and controls in TRAIN set
TRAIN_samples_1 <- subset(strata_data_new, strata_data_new$Cancer.Control.Status.1 ==1 & strata_data_new$label =="TRAIN")$Rand.Vol.ID
TRAIN_samples_0 <- subset(strata_data_new, strata_data_new$Cancer.Control.Status.1 ==0 & strata_data_new$label =="TRAIN")$Rand.Vol.ID

pairs_for_l2ocv <- data.frame()

# Collecting randomly choosen pairs
while(length(TRAIN_samples_1) >1) {
  print(length(TRAIN_samples_1))
  set.seed(2)
  case <- sample(TRAIN_samples_1, 1)
  set.seed(2)
  control <- sample(TRAIN_samples_0, 1)
  
  pairs_for_l2ocv <- rbind(pairs_for_l2ocv, data.frame(case, control))
  TRAIN_samples_0 <- setdiff(TRAIN_samples_0, control)
  TRAIN_samples_1 <- setdiff(TRAIN_samples_1, case)
}
pairs_for_l2ocv <- rbind(pairs_for_l2ocv, data.frame(case=TRAIN_samples_1, control=TRAIN_samples_0))

# creating 2 list (for usual and longit approaches) of datasets for l2ocv procedure
all_data <- list()
all_data_long <- list()

TRAIN_data <- subset(data, data$label == "TRAIN")
TRAIN_data_long <- subset(data_long, data_long$label == "TRAIN")
i=1
for(i in 1:nrow(pairs_for_l2ocv)) {
  print(i)
  
  samps <- as.vector(t(pairs_for_l2ocv[i,]))
  TRAIN_data$label <- "TRAIN"
  temp <- subset(TRAIN_data,TRAIN_data$main_id == samps[1] | TRAIN_data$main_id == samps[2])
  TRAIN_data[rownames(temp),]$label = "TEST"
  
  all_data[[i]] <- list()
  all_data[[i]][[1]]  <- TRAIN_data # dataset with right labels
  all_data[[i]][[2]]  <- proteins # set of proteins which will be used
  all_data[[i]][[3]]  <- str_c("inhouse_loocv_",i) #name of dataset
  
  TRAIN_data_long$label <- "TRAIN"
  temp_long <- subset(TRAIN_data_long,TRAIN_data_long$main_id == samps[1] | TRAIN_data_long$main_id == samps[2])
  TRAIN_data_long[rownames(temp_long),]$label = "TEST"
  
  all_data_long[[i]] <- list()
  all_data_long[[i]][[1]]  <- TRAIN_data_long # dataset with right labels
  all_data_long[[i]][[2]]  <- proteins # set of proteins which will be used
  all_data_long[[i]][[3]]  <- str_c("inhouse_long_loocv_",i) #name of dataset
}

i=i+1
all_data[[i]] <- list()
all_data[[i]][[1]]  <- data
all_data[[i]][[2]]  <- proteins
all_data[[i]][[3]]  <- str_c("inhouse_main_",i)

all_data_long[[i]] <- list()
all_data_long[[i]][[1]]  <- data_long
all_data_long[[i]][[2]]  <- proteins
all_data_long[[i]][[3]]  <- str_c("inhouse_long_main_",i)
##################################################

save(all_data, file = paste0(PATH, "1_DATA/all_data.RData"))
save(all_data_long, file = paste0(PATH, "1_DATA/all_data_long.RData"))
