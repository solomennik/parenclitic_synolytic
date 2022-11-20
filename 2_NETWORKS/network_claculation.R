PATH = "/Users/Shared/Previously Relocated Items/Security/OC_synolitic_and_longit/" #PATH to main folder

source(paste0(PATH, "functions/general_functions.R")) #add functions
source(paste0(PATH, "functions/networks_functions.R")) #add functions of networks

# load data
load(paste0(PATH,"1_DATA/all_data.RData")) 
load(paste0(PATH,"1_DATA/all_data_long.RData"))

###################################
# resulting directory
DIR <- paste0(PATH, "2_NETWORKS/")
#function calculation_PN(list_of_data, type_PN, type_data, DIR)
# create DIR in 2_NETWORKS which will called type_data_type_PN
# where for each dataframe from list_of_data will be created 2 files: chars and networks


calculation_PN(all_data[c(1:19)], "GLM", "Usual",DIR)
calculation_PN(all_data[c(1:19)], "SVM", "Usual",DIR)

calculation_PN(all_data_long[c(1:19)], "LONG_xgbTree", "Long",DIR)
calculation_PN(all_data_long[c(1:19)], "LONG_glmnet", "Long",DIR)
calculation_PN(all_data_long[c(1:19)], "LONG_nnet", "Long",DIR)
