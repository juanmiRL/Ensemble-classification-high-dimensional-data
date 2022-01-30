#===============================================================================
#=========================== READ Ensemble =====================================
#===============================================================================


library(tidyverse)
library(MLSeq)
library(S4Vectors)
library(DESeq2)
library(kernlab)
library(e1071)
library(gbm)
library(adabag)
library(RSNNS)
library(glmnet)
library(HDclassif)
library(gpls)
library(TCGAbiolinks)
library(limma)
library(edgeR)
library(caret)
library(factoextra)
library(tidymodels)
library(skimr)
library(DataExplorer)
library(ggpubr)
library(univariateML)
library(GGally)
library(parallel)
library(doParallel)
library(dplyr)
library(purrr)
library(doMC)
library(fastAdaboost)
library(xgboost)
library(caretEnsemble)
library(Rcpp)
library(h2o)


READ_results_anova_pvalue <- readRDS("C:/Users/juanm/Desktop/Ensemble-classification-high-dimensional-data/Models/READ/READ_results_anova_pvalue.rds")

filtrado_anova_pvalue_500 <-  READ_results_anova_pvalue %>% pull(gen) %>% head(500)
filtrado_anova_pvalue_100 <- READ_results_anova_pvalue %>% pull(gen) %>% head(100)
filtrado_anova_pvalue_50  <- READ_results_anova_pvalue %>% pull(gen) %>% head(50)


# load data with the next sentence
tcga_data <- readRDS("C:/Users/juanm/Desktop/Ensemble-classification-high-dimensional-data/Data/tcga_data_READ.RDS")


# Preprocessing 

clinical_data = colData(tcga_data)
group = factor(clinical_data$definition)
group = relevel(group, ref="Solid Tissue Normal")
design = model.matrix(~group)
head(design)
dge = DGEList( # creating a DGEList object
  counts=assay(tcga_data),
  samples=colData(tcga_data),
  genes=as.data.frame(rowData(tcga_data)))

# filtering
keep = filterByExpr(dge,design) # defining which genes to keep
dge = dge[keep,,keep.lib.sizes=FALSE] # filtering the dge object
rm(keep) #  use rm() to remove objects from memory if you don't need them anymore

dge = calcNormFactors(dge,method="TMM")
v = voom(dge,design,plot=TRUE)

fit = lmFit(v, design)
fit = eBayes(fit)

topGenes = topTable(fit, coef=1, sort.by="p")
print(topGenes)


# Split data (80% train 20% test)

# Transpose and make it into a matrix object
d_mat = as.matrix(t(v$E))
dim(d_mat)

# As before, we want this to be a factor
d_resp = as.factor(v$targets$definition)
levels(d_resp)

df <- cbind(d_mat,d_resp)
df <- as.data.frame(df)
colnames(df)[length(df)] # check class of the tumor 
names(df)[length(df)] <- "Class"
df$Class <- as.factor(df$Class)


# Train/ test split 
set.seed(86)
train <- createDataPartition(y = df$Class, p = 0.7, list = FALSE, times = 1)
x_train <- df[train, ]
x_test  <- df[-train, ]
y_train <- df$Class[train]
y_test  <- df$Class[-train]

rm(list = c("df","clinical_data", "d_mat","d_resp","dge", "train",
            "design","fit","group","tcga_data","topGenes","v"))


#===============================================================================
#=========================== Majority Vote READ Ensemble =======================
#===============================================================================


#-------------------#
# ANOVA P-VALUE 50
#-------------------#

library(statip)

# load the models 
READ_rf_pvalue_50 <- readRDS("C:/Users/juanm/Desktop/Ensemble-classification-high-dimensional-data/Models/READ/READ_rf_pvalue_50.rds")
READ_svmrad_pvalue_50 <- readRDS("C:/Users/juanm/Desktop/Ensemble-classification-high-dimensional-data/Models/READ/READ_svmrad_pvalue_50.rds")
READ_nnet_pvalue_50 <- readRDS("C:/Users/juanm/Desktop/Ensemble-classification-high-dimensional-data/Models/READ/READ_nnet_pvalue_50.rds")
READ_gbm_pvalue_50 <- readRDS("C:/Users/juanm/Desktop/Ensemble-classification-high-dimensional-data/Models/READ/READ_gbm_pvalue_50.rds")
READ_xgbm_pvalue_50 <- readRDS("C:/Users/juanm/Desktop/Ensemble-classification-high-dimensional-data/Models/READ/READ_xgbm_pvalue_50.rds")
READ_glmnet_pvalue_50 <- readRDS("C:/Users/juanm/Desktop/Ensemble-classification-high-dimensional-data/Models/READ/READ_glmnet_pvalue_50.rds")
READ_hdda_pvalue_50 <- readRDS("C:/Users/juanm/Desktop/Ensemble-classification-high-dimensional-data/Models/READ/READ_hdda_pvalue_50.rds")
READ_logitBoost_pvalue_50 <- readRDS("C:/Users/juanm/Desktop/Ensemble-classification-high-dimensional-data/Models/READ/READ_logitBoost_pvalue_50.rds")


# Predicting the classes

# library(Rcpp)
# READ_majority_rf_50 <- predict(object = READ_rf_pvalue_50, newdata = x_test)

READ_majority_svm_50 <- predict(object = READ_svmrad_pvalue_50, newdata = x_test)
READ_majority_nnet_50 <- predict(object = READ_nnet_pvalue_50, newdata = x_test)
READ_majority_gbm_50 <- predict(object = READ_gbm_pvalue_50, newdata = x_test)
READ_majority_xgbm_50 <- predict(object = READ_xgbm_pvalue_50, newdata = x_test)
READ_majority_glmnet_50 <- predict(object = READ_glmnet_pvalue_50, newdata = x_test)
READ_majority_hdda_50 <- predict(object = READ_hdda_pvalue_50, newdata = x_test)
READ_majority_logitBoost_50 <- predict(object = READ_logitBoost_pvalue_50, newdata = x_test)


READ_majority_vote_50 <- as.data.frame(READ_majority_svm_50)


READ_majority_vote_50 <- cbind(READ_majority_vote_50 ,READ_majority_nnet_50, READ_majority_gbm_50,
                               READ_majority_xgbm_50, READ_majority_glmnet_50, READ_majority_hdda_50,  
                               READ_majority_logitBoost_50)


# apply mode function for all predictions 

READ_majority_vote_50$Majority_vote <- apply(READ_majority_vote_50, MARGIN = 1, FUN = mfv)

READ_majority_vote_50$Majority_vote <- as.factor(READ_majority_vote_50$Majority_vote)

caret::confusionMatrix(READ_majority_vote_50$Majority_vote, y_test)

# print all confusion matrix 

for (i in 1:ncol(READ_majority_vote_50)) {
  cat(colnames(READ_majority_vote_50[i]))
  cat("\n", "\n")
  print(caret::confusionMatrix(READ_majority_vote_50[,i], y_test))
}

#-------------------#
# ANOVA P-VALUE 100
#-------------------#


# load the models 
READ_rf_pvalue_100 <- readRDS("C:/Users/juanm/Desktop/Ensemble-classification-high-dimensional-data/Models/READ/READ_rf_pvalue_100.rds")
READ_svmrad_pvalue_100 <- readRDS("C:/Users/juanm/Desktop/Ensemble-classification-high-dimensional-data/Models/READ/READ_svmrad_pvalue_100.rds")
READ_nnet_pvalue_100 <- readRDS("C:/Users/juanm/Desktop/Ensemble-classification-high-dimensional-data/Models/READ/READ_nnet_pvalue_100.rds")
READ_gbm_pvalue_100 <- readRDS("C:/Users/juanm/Desktop/Ensemble-classification-high-dimensional-data/Models/READ/READ_gbm_pvalue_100.rds")
READ_xgbm_pvalue_100 <- readRDS("C:/Users/juanm/Desktop/Ensemble-classification-high-dimensional-data/Models/READ/READ_xgbm_pvalue_100.rds")
READ_glmnet_pvalue_100 <- readRDS("C:/Users/juanm/Desktop/Ensemble-classification-high-dimensional-data/Models/READ/READ_glmnet_pvalue_100.rds")
READ_hdda_pvalue_100 <- readRDS("C:/Users/juanm/Desktop/Ensemble-classification-high-dimensional-data/Models/READ/READ_hdda_pvalue_100.rds")
READ_logitBoost_pvalue_100 <- readRDS("C:/Users/juanm/Desktop/Ensemble-classification-high-dimensional-data/Models/READ/READ_logitBoost_pvalue_100.rds")


# Predicting the classes

READ_majority_svm_100 <- predict(object = READ_svmrad_pvalue_100, newdata = x_test)
READ_majority_nnet_100 <- predict(object = READ_nnet_pvalue_100, newdata = x_test)
READ_majority_gbm_100 <- predict(object = READ_gbm_pvalue_100, newdata = x_test)
READ_majority_xgbm_100 <- predict(object = READ_xgbm_pvalue_100, newdata = x_test)
READ_majority_glmnet_100 <- predict(object = READ_glmnet_pvalue_100, newdata = x_test)
READ_majority_hdda_100 <- predict(object = READ_hdda_pvalue_100, newdata = x_test)
READ_majority_logitBoost_100 <- predict(object = READ_logitBoost_pvalue_100, newdata = x_test)


READ_majority_vote_100 <- as.data.frame(READ_majority_svm_100)


READ_majority_vote_100 <- cbind(READ_majority_vote_100 ,READ_majority_nnet_100, READ_majority_gbm_100,
                                READ_majority_xgbm_100, READ_majority_glmnet_100, READ_majority_hdda_100,  
                                READ_majority_logitBoost_100)


# apply mode function for all predictions 

READ_majority_vote_100$Majority_vote <- apply(READ_majority_vote_100, MARGIN = 1, FUN = mfv)

READ_majority_vote_100$Majority_vote <- as.factor(READ_majority_vote_100$Majority_vote)

caret::confusionMatrix(READ_majority_vote_100$Majority_vote, y_test)

# print all confusion matrix 

for (i in 1:ncol(READ_majority_vote_100)) {
  cat(colnames(READ_majority_vote_100[i]))
  cat("\n", "\n")
  print(caret::confusionMatrix(READ_majority_vote_100[,i], y_test))
}


#-------------------#
# ANOVA P-VALUE 500
#-------------------#


# load the models 
READ_rf_pvalue_500 <- readRDS("C:/Users/juanm/Desktop/Ensemble-classification-high-dimensional-data/READ/READ_rf_pvalue_500.rds")
READ_svmrad_pvalue_500 <- readRDS("C:/Users/juanm/Desktop/Ensemble-classification-high-dimensional-data/READ/READ_svmrad_pvalue_500.rds")
READ_nnet_pvalue_500 <- readRDS("C:/Users/juanm/Desktop/Ensemble-classification-high-dimensional-data/READ/READ_nnet_pvalue_500.rds")
READ_gbm_pvalue_500 <- readRDS("C:/Users/juanm/Desktop/Ensemble-classification-high-dimensional-data/READ/READ_gbm_pvalue_500.rds")
READ_xgbm_pvalue_500 <- readRDS("C:/Users/juanm/Desktop/Ensemble-classification-high-dimensional-data/READ/READ_xgbm_pvalue_500.rds")
READ_glmnet_pvalue_500 <- readRDS("C:/Users/juanm/Desktop/Ensemble-classification-high-dimensional-data/READ/READ_glmnet_pvalue_500.rds")
READ_hdda_pvalue_500 <- readRDS("C:/Users/juanm/Desktop/Ensemble-classification-high-dimensional-data/READ/READ_hdda_pvalue_500.rds")
READ_logitBoost_pvalue_500 <- readRDS("C:/Users/juanm/Desktop/Ensemble-classification-high-dimensional-data/READ/READ_logitBoost_pvalue_500.rds")


# Predicting the classes

READ_majority_svm_500 <- predict(object = READ_svmrad_pvalue_500, newdata = x_test)
READ_majority_nnet_500 <- predict(object = READ_nnet_pvalue_500, newdata = x_test)
READ_majority_gbm_500 <- predict(object = READ_gbm_pvalue_500, newdata = x_test)
READ_majority_xgbm_500 <- predict(object = READ_xgbm_pvalue_500, newdata = x_test)
READ_majority_glmnet_500 <- predict(object = READ_glmnet_pvalue_500, newdata = x_test)
READ_majority_hdda_500 <- predict(object = READ_hdda_pvalue_500, newdata = x_test)
READ_majority_logitBoost_500 <- predict(object = READ_logitBoost_pvalue_500, newdata = x_test)


READ_majority_vote_500 <- as.data.frame(READ_majority_svm_500)


READ_majority_vote_500 <- cbind(READ_majority_vote_500 ,READ_majority_nnet_500, READ_majority_gbm_500,
                                READ_majority_xgbm_500, READ_majority_glmnet_500, READ_majority_hdda_500,  
                                READ_majority_logitBoost_500)


# apply mode function for all predictions 

READ_majority_vote_500$Majority_vote <- apply(READ_majority_vote_500, MARGIN = 1, FUN = mfv)

READ_majority_vote_500$Majority_vote <- as.factor(READ_majority_vote_500$Majority_vote)

caret::confusionMatrix(READ_majority_vote_500$Majority_vote, y_test)

# print all confusion matrix 

for (i in 1:ncol(READ_majority_vote_500)) {
  cat(colnames(READ_majority_vote_500[i]))
  cat("\n", "\n")
  print(caret::confusionMatrix(READ_majority_vote_500[,i], y_test))
}


#-------------------#
#      RFE 50
#-------------------#

# load the models 
READ_rf_RFE_50 <- readRDS("C:/Users/juanm/Desktop/Ensemble-classification-high-dimensional-data/READ/READ_rf_RFE_50.rds")
READ_svmrad_RFE_50 <- readRDS("C:/Users/juanm/Desktop/Ensemble-classification-high-dimensional-data/READ/READ_svmrad_RFE_50.rds")
READ_nnet_RFE_50 <- readRDS("C:/Users/juanm/Desktop/Ensemble-classification-high-dimensional-data/READ/READ_nnet_RFE_50.rds")
READ_gbm_RFE_50 <- readRDS("C:/Users/juanm/Desktop/Ensemble-classification-high-dimensional-data/READ/READ_gbm_RFE_50.rds")
READ_xgbm_RFE_50 <- readRDS("C:/Users/juanm/Desktop/Ensemble-classification-high-dimensional-data/READ/READ_xgbm_RFE_50.rds")
READ_glmnet_RFE_50 <- readRDS("C:/Users/juanm/Desktop/Ensemble-classification-high-dimensional-data/READ/READ_glmnet_RFE_50.rds")
READ_hdda_RFE_50 <- readRDS("C:/Users/juanm/Desktop/Ensemble-classification-high-dimensional-data/READ/READ_hdda_RFE_50.rds")
READ_logitBoost_RFE_50 <- readRDS("C:/Users/juanm/Desktop/Ensemble-classification-high-dimensional-data/READ/READ_logitBoost_RFE_50.rds")


# load RFE selected Variables 
RFE_results <- readRDS("C:/Users/juanm/Desktop/Ensemble-classification-high-dimensional-data/READ/READ_RFE_results.rds")
selected_vars <- RFE_results$variables
best_vars_50 <- RFE_results$control$functions$selectVar(selected_vars, 50)

x_train_rfe_50 <- x_train[, c("Class", best_vars_50)] 
x_test_rfe_50 <- x_test[, c("Class", best_vars_50)] 



# Predicting the classes


READ_majority_svm_RFE_50 <- predict(object = READ_svmrad_RFE_50, newdata = x_test_rfe_50)
READ_majority_nnet_RFE_50 <- predict(object = READ_nnet_RFE_50, newdata = x_test_rfe_50)
READ_majority_gbm_RFE_50 <- predict(object = READ_gbm_RFE_50, newdata = x_test_rfe_50)
READ_majority_xgbm_RFE_50 <- predict(object = READ_xgbm_RFE_50, newdata = x_test_rfe_50)
READ_majority_glmnet_RFE_50 <- predict(object = READ_glmnet_RFE_50, newdata = x_test_rfe_50)
READ_majority_hdda_RFE_50 <- predict(object = READ_hdda_RFE_50, newdata = x_test_rfe_50)
READ_majority_logitBoost_RFE_50 <- predict(object = READ_logitBoost_RFE_50, newdata = x_test_rfe_50)


READ_majority_vote_RFE_50 <- as.data.frame(READ_majority_svm_RFE_50)


READ_majority_vote_RFE_50 <- cbind(READ_majority_vote_RFE_50 ,READ_majority_nnet_RFE_50, READ_majority_gbm_RFE_50,
                                   READ_majority_xgbm_RFE_50, READ_majority_glmnet_RFE_50, READ_majority_hdda_RFE_50,  
                                   READ_majority_logitBoost_RFE_50)


# apply mode function for all predictions 

READ_majority_vote_RFE_50$Majority_vote_RFE_50 <- apply(READ_majority_vote_RFE_50, MARGIN = 1, FUN = mfv)

READ_majority_vote_RFE_50$Majority_vote_RFE_50 <- as.factor(READ_majority_vote_RFE_50$Majority_vote_RFE_50)

caret::confusionMatrix(READ_majority_vote_RFE_50$Majority_vote_RFE_50, y_test)

# print all confusion matrix 

for (i in 1:ncol(READ_majority_vote_RFE_50)) {
  cat(colnames(READ_majority_vote_RFE_50[i]))
  cat("\n", "\n")
  print(caret::confusionMatrix(READ_majority_vote_RFE_50[,i], y_test))
}


#-------------------#
#      RFE 100
#-------------------#

# load the models 
READ_rf_RFE_100 <- readRDS("C:/Users/juanm/Desktop/Ensemble-classification-high-dimensional-data/READ/READ_rf_RFE_100.rds")
READ_svmrad_RFE_100 <- readRDS("C:/Users/juanm/Desktop/Ensemble-classification-high-dimensional-data/READ/READ_svmrad_RFE_100.rds")
READ_nnet_RFE_100 <- readRDS("C:/Users/juanm/Desktop/Ensemble-classification-high-dimensional-data/READ/READ_nnet_RFE_100.rds")
READ_gbm_RFE_100 <- readRDS("C:/Users/juanm/Desktop/Ensemble-classification-high-dimensional-data/READ/READ_gbm_RFE_100.rds")
READ_xgbm_RFE_100 <- readRDS("C:/Users/juanm/Desktop/Ensemble-classification-high-dimensional-data/READ/READ_xgbm_RFE_100.rds")
READ_glmnet_RFE_100 <- readRDS("C:/Users/juanm/Desktop/Ensemble-classification-high-dimensional-data/READ/READ_glmnet_RFE_100.rds")
READ_hdda_RFE_100 <- readRDS("C:/Users/juanm/Desktop/Ensemble-classification-high-dimensional-data/READ/READ_hdda_RFE_100.rds")
READ_logitBoost_RFE_100 <- readRDS("C:/Users/juanm/Desktop/Ensemble-classification-high-dimensional-data/READ/READ_logitBoost_RFE_100.rds")


# load RFE selected Variables 
RFE_results <- readRDS("C:/Users/juanm/Desktop/Ensemble-classification-high-dimensional-data/READ/READ_RFE_results.rds")
selected_vars <- RFE_results$variables
best_vars_100 <- RFE_results$control$functions$selectVar(selected_vars, 100)

x_train_rfe_100 <- x_train[, c("Class", best_vars_100)] 
x_test_rfe_100 <- x_test[, c("Class", best_vars_100)] 

# Predicting the classes

READ_majority_svm_RFE_100 <- predict(object = READ_svmrad_RFE_100, newdata = x_test_rfe_100)
READ_majority_nnet_RFE_100 <- predict(object = READ_nnet_RFE_100, newdata = x_test_rfe_100)
READ_majority_gbm_RFE_100 <- predict(object = READ_gbm_RFE_100, newdata = x_test_rfe_100)
READ_majority_xgbm_RFE_100 <- predict(object = READ_xgbm_RFE_100, newdata = x_test_rfe_100)
READ_majority_glmnet_RFE_100 <- predict(object = READ_glmnet_RFE_100, newdata = x_test_rfe_100)
READ_majority_hdda_RFE_100 <- predict(object = READ_hdda_RFE_100, newdata = x_test_rfe_100)
READ_majority_logitBoost_RFE_100 <- predict(object = READ_logitBoost_RFE_100, newdata = x_test_rfe_100)


READ_majority_vote_RFE_100 <- as.data.frame(READ_majority_svm_RFE_100)


READ_majority_vote_RFE_100 <- cbind(READ_majority_vote_RFE_100 ,READ_majority_nnet_RFE_100, READ_majority_gbm_RFE_100,
                                    READ_majority_xgbm_RFE_100, READ_majority_glmnet_RFE_100, READ_majority_hdda_RFE_100,  
                                    READ_majority_logitBoost_RFE_100)


# apply mode function for all predictions 

READ_majority_vote_RFE_100$Majority_vote_RFE_100 <- apply(READ_majority_vote_RFE_100, MARGIN = 1, FUN = mfv)

READ_majority_vote_RFE_100$Majority_vote_RFE_100 <- as.factor(READ_majority_vote_RFE_100$Majority_vote_RFE_100)

caret::confusionMatrix(READ_majority_vote_RFE_100$Majority_vote_RFE_100, y_test)

# print all confusion matrix 

for (i in 1:ncol(READ_majority_vote_RFE_100)) {
  cat(colnames(READ_majority_vote_RFE_100[i]))
  cat("\n", "\n")
  print(caret::confusionMatrix(READ_majority_vote_RFE_100[,i], y_test))
}


#===============================================================================
#=========================== Averaging  READ Ensemble ==========================
#===============================================================================


#-------------------#
#  ANOVA PVALUE 50
#-------------------#


READ_averaging_nnet_50 <- predict(object = READ_nnet_pvalue_50, newdata = x_test, type = "prob")
READ_averaging_gbm_50 <- predict(object = READ_gbm_pvalue_50, newdata = x_test, type = "prob")
READ_averaging_xgbm_50 <- predict(object = READ_xgbm_pvalue_50, newdata = x_test, type = "prob")
READ_averaging_logitBoost_50 <- predict(object = READ_logitBoost_pvalue_50, newdata = x_test, type = "prob")

# Taking average of predictions

READ_averaging_pvalue_50 <-(READ_averaging_nnet_50$`1` + READ_averaging_gbm_50$`1` + 
                              READ_averaging_xgbm_50$`1`  + READ_averaging_logitBoost_50$`1`) / 4

# Splitting into binary classes at 0.5

READ_averaging_pvalue_50 <- as.factor(if_else(READ_averaging_pvalue_50 > 0.5, 1, 3))

caret::confusionMatrix(READ_averaging_pvalue_50, y_test)


#-------------------#
#  ANOVA PVALUE 100
#-------------------#


READ_averaging_nnet_100 <- predict(object = READ_nnet_pvalue_100, newdata = x_test, type = "prob")
READ_averaging_gbm_100 <- predict(object = READ_gbm_pvalue_100, newdata = x_test, type = "prob")
READ_averaging_xgbm_100 <- predict(object = READ_xgbm_pvalue_100, newdata = x_test, type = "prob")
READ_averaging_logitBoost_100 <- predict(object = READ_logitBoost_pvalue_100, newdata = x_test, type = "prob")

# Taking average of predictions

READ_averaging_pvalue_100 <-(READ_averaging_nnet_100$`1` + READ_averaging_gbm_100$`1` + 
                               READ_averaging_xgbm_100$`1` + + READ_averaging_logitBoost_100$`1`) / 4

# Splitting into binary classes at 0.5

READ_averaging_pvalue_100<- as.factor(if_else(READ_averaging_pvalue_100 > 0.5, 1, 3))

caret::confusionMatrix(READ_averaging_pvalue_100, y_test)


#-------------------#
#  ANOVA PVALUE 500
#-------------------#


READ_averaging_nnet_500 <- predict(object = READ_nnet_pvalue_500, newdata = x_test, type = "prob")
READ_averaging_gbm_500 <- predict(object = READ_gbm_pvalue_500, newdata = x_test, type = "prob")
READ_averaging_xgbm_500 <- predict(object = READ_xgbm_pvalue_500, newdata = x_test, type = "prob")
READ_averaging_logitBoost_500 <- predict(object = READ_logitBoost_pvalue_500, newdata = x_test, type = "prob")

# Taking average of predictions

READ_averaging_pvalue_500 <-(READ_averaging_nnet_500$`1` + READ_averaging_gbm_500$`1` + 
                               READ_averaging_xgbm_500$`1` + + READ_averaging_logitBoost_500$`1`) / 4

# Splitting into binary classes at 0.5

READ_averaging_pvalue_500 <- as.factor(if_else(READ_averaging_pvalue_500 > 0.5, 1, 3))

caret::confusionMatrix(READ_averaging_pvalue_500, y_test)


#-------------------#
#      RFE 50
#-------------------#


READ_averaging_nnet_RFE_50 <- predict(object = READ_nnet_RFE_50, newdata = x_test_rfe_50, type = "prob")
READ_averaging_gbm_RFE_50 <- predict(object = READ_gbm_RFE_50, newdata = x_test_rfe_50, type = "prob")
READ_averaging_xgbm_RFE_50 <- predict(object = READ_xgbm_RFE_50, newdata = x_test_rfe_50, type = "prob")
READ_averaging_logitBoost_RFE_50 <- predict(object = READ_logitBoost_RFE_50, newdata = x_test_rfe_50, type = "prob")

# Taking average of predictions

READ_averaging_RFE_50 <- (READ_averaging_nnet_RFE_50$`1` + READ_averaging_gbm_RFE_50$`1` + 
                            READ_averaging_xgbm_RFE_50$`1` + + READ_averaging_logitBoost_RFE_50$`1`) / 4

# Splitting into binary classes at 0.5

READ_averaging_RFE_50 <- as.factor(if_else(READ_averaging_RFE_50 > 0.5, 1, 3))

caret::confusionMatrix(READ_averaging_RFE_50, y_test)


#-------------------#
#      RFE 100
#-------------------#


READ_averaging_nnet_RFE_100 <- predict(object = READ_nnet_RFE_100, newdata = x_test_rfe_100, type = "prob")
READ_averaging_gbm_RFE_100 <- predict(object = READ_gbm_RFE_100, newdata = x_test_rfe_100, type = "prob")
READ_averaging_xgbm_RFE_100 <- predict(object = READ_xgbm_RFE_100, newdata = x_test_rfe_100, type = "prob")
READ_averaging_logitBoost_RFE_100 <- predict(object = READ_logitBoost_RFE_100, newdata = x_test_rfe_100, type = "prob")

# Taking average of predictions

READ_averaging_RFE_100 <-(READ_averaging_nnet_RFE_100$`1` + READ_averaging_gbm_RFE_100$`1` + 
                            READ_averaging_xgbm_RFE_100$`1` + + READ_averaging_logitBoost_RFE_100$`1`) / 4

# Splitting into binary classes at 0.5

READ_averaging_RFE_100 <- as.factor(if_else(READ_averaging_RFE_100 > 0.5, 1, 3))

caret::confusionMatrix(READ_averaging_RFE_100, y_test)


#===============================================================================
#=========================== Weighted Average  READ Ensemble ===================
#===============================================================================

#-------------------#
#  ANOVA PVALUE 50
#-------------------#

# Taking Weighted average of predictions

READ_weigthed_average_pvalue_50 <- (0.15 * READ_averaging_nnet_50$`1` + 0.15 * READ_averaging_xgbm_50$`1` +
                                      + 0.15 * READ_averaging_gbm_50$`1` + 0.45 * READ_averaging_logitBoost_50$`1`) 

# Splitting into binary classes at 0.5

READ_weigthed_average_pvalue_50 <- as.factor(if_else(READ_weigthed_average_pvalue_50 > 0.5, 1, 3))

caret::confusionMatrix(READ_weigthed_average_pvalue_50, y_test)

#-------------------#
#  ANOVA PVALUE 100
#-------------------#


# Taking Weighted average of predictions

READ_weigthed_average_pvalue_100 <-(0.45 * READ_averaging_nnet_100$`1` + 0.15 * READ_averaging_xgbm_100$`1` +
                                      + 0.15 * READ_averaging_gbm_100$`1` + 0.15 * READ_averaging_logitBoost_100$`1`) 

# Splitting into binary classes at 0.5

READ_weigthed_average_pvalue_100 <- as.factor(if_else(READ_weigthed_average_pvalue_100 > 0.5, 1, 3))

caret::confusionMatrix(READ_weigthed_average_pvalue_100, y_test)

#-------------------#
#  ANOVA PVALUE 500
#-------------------#


# Taking Weighted average of predictions

READ_weigthed_average_pvalue_500 <-(0.45 * READ_averaging_nnet_500$`1` + 0.15 * READ_averaging_xgbm_500$`1` +
                                      + 0.15 * READ_averaging_gbm_500$`1` + 0.15 * READ_averaging_logitBoost_500$`1`) 

# Splitting into binary classes at 0.5

READ_weigthed_average_pvalue_500 <- as.factor(if_else(READ_weigthed_average_pvalue_500 > 0.5, 1, 3))

caret::confusionMatrix(READ_weigthed_average_pvalue_500, y_test)


#-------------------#
#      RFE 50
#-------------------#


# Taking Weighted average of predictions

READ_weigthed_average_RFE_50 <- (0.45 * READ_averaging_nnet_RFE_50$`1` + 0.15 * READ_averaging_gbm_RFE_50$`1` + 
                                   0.15 * READ_averaging_xgbm_RFE_50$`1`  + 0.15 * READ_averaging_logitBoost_RFE_50$`1`) 

# Splitting into binary classes at 0.5

READ_weigthed_average_RFE_50 <- as.factor(if_else(READ_weigthed_average_RFE_50 > 0.5, 1, 3))

caret::confusionMatrix(READ_weigthed_average_RFE_50, y_test)

#-------------------#
#      RFE 100
#-------------------#


# Taking Weighted average of predictions

READ_weigthed_average_RFE_100 <- (0.45 * READ_averaging_nnet_RFE_100$`1` + 0.15 * READ_averaging_gbm_RFE_100$`1` + 
                                    0.15 * READ_averaging_xgbm_RFE_100$`1`  + 0.15 * READ_averaging_logitBoost_RFE_100$`1`) 

# Splitting into binary classes at 0.5

READ_weigthed_average_RFE_100 <- as.factor(if_else(READ_weigthed_average_RFE_100 > 0.5, 1, 3))

caret::confusionMatrix(READ_weigthed_average_RFE_100, y_test)


#-------------------------------------------------------------------------------------------------#
#---------------------- Stacked Ensembles Using H20 ----------------------------------------------# 
#-------------------------------------------------------------------------------------------------#


#-------------------#
#  ANOVA PVALUE 50
#-------------------#

library(h2o)

h2o.init()

READ_results_anova_pvalue <- readRDS("C:/Users/juanm/Desktop/Ensemble-classification-high-dimensional-data/Models/READ/READ_results_anova_pvalue.rds")

filtrado_anova_pvalue_500 <-  READ_results_anova_pvalue %>% pull(gen) %>% head(500)
filtrado_anova_pvalue_100 <- READ_results_anova_pvalue %>% pull(gen) %>% head(100)
filtrado_anova_pvalue_50  <- READ_results_anova_pvalue %>% pull(gen) %>% head(50)


# load data with the next sentence
tcga_data <- readRDS("C:/Users/juanm/Desktop/Ensemble-classification-high-dimensional-data/Data/tcga_data_READ.RDS")


# Preprocessing 

clinical_data = colData(tcga_data)
group = factor(clinical_data$definition)
group = relevel(group, ref="Solid Tissue Normal")
design = model.matrix(~group)
head(design)
dge = DGEList( # creating a DGEList object
  counts=assay(tcga_data),
  samples=colData(tcga_data),
  genes=as.data.frame(rowData(tcga_data)))

# filtering
keep = filterByExpr(dge,design) # defining which genes to keep
dge = dge[keep,,keep.lib.sizes=FALSE] # filtering the dge object
rm(keep) #  use rm() to remove objects from memory if you don't need them anymore

dge = calcNormFactors(dge,method="TMM")
v = voom(dge,design,plot=TRUE)

fit = lmFit(v, design)
fit = eBayes(fit)

topGenes = topTable(fit, coef=1, sort.by="p")
print(topGenes)


# Split data (80% train 20% test)

# Transpose and make it into a matrix object
d_mat = as.matrix(t(v$E))
dim(d_mat)

# As before, we want this to be a factor
d_resp = as.factor(v$targets$definition)
levels(d_resp)

df <- cbind(d_mat,d_resp)
df <- as.data.frame(df)
colnames(df)[length(df)] # check class of the tumor 
names(df)[length(df)] <- "Class"
df$Class <- as.factor(df$Class)


# Train/ test split 
set.seed(86)
train <- createDataPartition(y = df$Class, p = 0.7, list = FALSE, times = 1)
x_train <- df[train, ]
x_test  <- df[-train, ]
y_train <- df$Class[train]
y_test  <- df$Class[-train]


x_train <- x_train[c("Class", filtrado_anova_pvalue_50)]
x_test <- x_test[c("Class", filtrado_anova_pvalue_50)]

rm(list = c("df","clinical_data", "d_mat","d_resp","dge", "train",
            "design","fit","group","tcga_data","topGenes","v"))

# create an h20 object 
x_train <- as.h2o(x_train)
x_test <- as.h2o(x_test)
y_train <- as.h2o(y_train)
y_test <- as.h2o(y_test)

# Identify predictors and response
y = "Class"
x <- setdiff(names(x_train), y)


# AutoML: Automatic Machine Learning
aml <- h2o.automl(x = x, y = y,
                  training_frame = x_train,
                  max_models = 10,
                  seed = 86)

# View the AutoML Leaderboard
lb <- aml@leaderboard
print(lb, n = nrow(lb))  # Print all rows instead of default (6 rows)

# The leader model is stored here
aml@leader


# To generate predictions on a test set, you can make predictions
# directly on the `H2OAutoML` object or on the leader model
# object directly
pred <- h2o.predict(aml, x_test)  # predict(aml, test) also works

# Get leaderboard with all possible columns
lb <- h2o.get_leaderboard(object = aml, extra_columns = "ALL")
lb

# # Get the best model using the metric
# m <- aml@leader
# # this is equivalent to
# m <- h2o.get_best_model(aml)
# 
# # Get the best model using a non-default metric
# m <- h2o.get_best_model(aml, criterion = "logloss")
# 
# # Get the best XGBoost model using default sort metric
# xgb <- h2o.get_best_model(aml, algorithm = "xgboost")
# 
# # Get the best XGBoost model, ranked by logloss
# xgb <- h2o.get_best_model(aml, algorithm = "xgboost", criterion = "logloss")


# Auto ML Stacked Ensemble


# Get a specific model by model ID
m <- h2o.getModel("StackedEnsemble_AllModels_AutoML_20220130_080356")

m_pred <- h2o.predict(aml, x_test)
m_pred <- m_pred[, 1]
m_pred <- as.data.frame(m_pred)
m_pred <- as.factor(m_pred$predict)
y_test <- as.data.frame(y_test)
y_test <- as.factor(y_test$`x`)
levels(m_pred) <- levels(y_test)
caret::confusionMatrix(m_pred, y_test)

#########################################################################################
# --------------------------- Stacked Ensemble -----------------------------------------#
#########################################################################################



# Number of CV folds (to generate level-one data for stacking)
nfolds <- 5

# There are a few ways to assemble a list of models to stack together:
# 1. Train individual models and put them in a list
# 2. Train a grid of models
# 3. Train several grids of models
# Note: All base models must have the same cross-validation folds and
# the cross-validated predicted values must be kept.


# 1. Generate a 2-model ensemble (GBM + RF)

# Train & Cross-validate a GBM
my_gbm <- h2o.gbm(x = x,
                  y = y,
                  training_frame = x_train,
                  distribution = "multinomial",
                  ntrees = 500,
                  max_depth = 3,
                  min_rows = 2,
                  learn_rate = 0.01,
                  nfolds = nfolds,
                  keep_cross_validation_predictions = TRUE,
                  seed = 86)

h2o.performance(my_gbm, newdata = x_test)



# Train & Cross-validate a RF
my_rf <- h2o.randomForest(x = x,
                          y = y,
                          training_frame = x_train,
                          ntrees = 500,
                          nfolds = nfolds,
                          keep_cross_validation_predictions = TRUE,
                          seed = 86)

h2o.performance(my_rf, newdata = x_test)


my_glm <- h2o.glm(
  y = y,
  x = x,
  training_frame = x_train,
  family = "multinomial",
  link   = "family_default",
  standardize  = TRUE,
  balance_classes   = FALSE,
  ignore_const_cols = TRUE,
  missing_values_handling = "Skip",
  lambda_search = TRUE,
  solver = "AUTO",
  alpha  = 0.95,
  seed = 86,
  nfolds = nfolds,
  keep_cross_validation_predictions = TRUE,
  model_id = "my_glm"
)

h2o.performance(my_glm, newdata = x_test)


# Train a stacked ensemble using the GBM and RF above
ensemble <- h2o.stackedEnsemble(x = x,
                                y = y,
                                training_frame = x_train,
                                base_models = list(my_gbm, my_rf, my_glm))

# Eval ensemble performance on a test set
perf <- h2o.performance(ensemble, newdata = x_test)

# Compare to base learner performance on the test set
perf_gbm_test <- h2o.performance(my_gbm, newdata = x_test)
perf_rf_test <- h2o.performance(my_rf, newdata = x_test)
baselearner_best_auc_test <- max(h2o.auc(perf_gbm_test), h2o.auc(perf_rf_test))
ensemble_auc_test <- h2o.auc(perf)
print(sprintf("Best Base-learner Test AUC:  %s", baselearner_best_auc_test))
print(sprintf("Ensemble Test AUC:  %s", ensemble_auc_test))
# [1] "Best Base-learner Test AUC:  1"
# [1] "Ensemble Test AUC:  1"

# Generate predictions on a test set 
pred <- h2o.predict(ensemble, newdata = x_test)
pred <- pred[, 1]
pred <- as.data.frame(pred)
pred <- as.factor(pred$predict)
levels(pred) <- levels(y_test)
caret::confusionMatrix(pred, y_test)


my_gbm_pred <- h2o.predict(my_gbm, newdata = x_test)
my_gbm_pred <- my_gbm_pred[, 1]
my_gbm_pred <- as.data.frame(my_gbm_pred)
my_gbm_pred <- as.factor(my_gbm_pred$predict)
levels(my_gbm_pred) <- levels(y_test)
caret::confusionMatrix(my_gbm_pred, y_test)

my_rf_pred <- h2o.predict(my_rf, newdata = x_test)
my_rf_pred <- my_rf_pred[, 1]
my_rf_pred <- as.data.frame(my_rf_pred)
my_rf_pred <- as.factor(my_rf_pred$predict)
levels(my_rf_pred) <- levels(y_test)
caret::confusionMatrix(my_rf_pred, y_test)


my_glm_pred <- h2o.predict(my_glm, newdata = x_test)
my_glm_pred <- my_glm_pred[, 1]
my_glm_pred <- as.data.frame(my_glm_pred)
my_glm_pred <- as.factor(my_glm_pred$predict)
levels(my_glm_pred) <- levels(y_test)
caret::confusionMatrix(my_glm_pred, y_test)


#-------------------#
#  ANOVA PVALUE 100
#-------------------#

library(h2o)

h2o.init()

READ_results_anova_pvalue <- readRDS("C:/Users/juanm/Desktop/Ensemble-classification-high-dimensional-data/Models/READ/READ_results_anova_pvalue.rds")

filtrado_anova_pvalue_500 <-  READ_results_anova_pvalue %>% pull(gen) %>% head(500)
filtrado_anova_pvalue_100 <- READ_results_anova_pvalue %>% pull(gen) %>% head(100)
filtrado_anova_pvalue_50  <- READ_results_anova_pvalue %>% pull(gen) %>% head(50)


# load data with the next sentence
tcga_data <- readRDS("C:/Users/juanm/Desktop/Ensemble-classification-high-dimensional-data/Data/tcga_data_READ.RDS")


# Preprocessing 

clinical_data = colData(tcga_data)
group = factor(clinical_data$definition)
group = relevel(group, ref="Solid Tissue Normal")
design = model.matrix(~group)
head(design)
dge = DGEList( # creating a DGEList object
  counts=assay(tcga_data),
  samples=colData(tcga_data),
  genes=as.data.frame(rowData(tcga_data)))

# filtering
keep = filterByExpr(dge,design) # defining which genes to keep
dge = dge[keep,,keep.lib.sizes=FALSE] # filtering the dge object
rm(keep) #  use rm() to remove objects from memory if you don't need them anymore

dge = calcNormFactors(dge,method="TMM")
v = voom(dge,design,plot=TRUE)

fit = lmFit(v, design)
fit = eBayes(fit)

topGenes = topTable(fit, coef=1, sort.by="p")
print(topGenes)


# Split data (80% train 20% test)

# Transpose and make it into a matrix object
d_mat = as.matrix(t(v$E))
dim(d_mat)

# As before, we want this to be a factor
d_resp = as.factor(v$targets$definition)
levels(d_resp)

df <- cbind(d_mat,d_resp)
df <- as.data.frame(df)
colnames(df)[length(df)] # check class of the tumor 
names(df)[length(df)] <- "Class"
df$Class <- as.factor(df$Class)

# Train/ test split 
set.seed(86)
train <- createDataPartition(y = df$Class, p = 0.7, list = FALSE, times = 1)
x_train <- df[train, ]
x_test  <- df[-train, ]
y_train <- df$Class[train]
y_test  <- df$Class[-train]


x_train <- x_train[c("Class", filtrado_anova_pvalue_100)]
x_test <- x_test[c("Class", filtrado_anova_pvalue_100)]

rm(list = c("df","clinical_data", "d_mat","d_resp","dge", "train",
            "design","fit","group","tcga_data","topGenes","v"))

# create an h20 object 
x_train <- as.h2o(x_train)
x_test <- as.h2o(x_test)
y_train <- as.h2o(y_train)
y_test <- as.h2o(y_test)

# Identify predictors and response
y = "Class"
x <- setdiff(names(x_train), y)


# AutoML: Automatic Machine Learning
aml <- h2o.automl(x = x, y = y,
                  training_frame = x_train,
                  max_models = 10,
                  seed = 86)

# View the AutoML Leaderboard
lb <- aml@leaderboard
print(lb, n = nrow(lb))  # Print all rows instead of default (6 rows)

# The leader model is stored here
aml@leader


# To generate predictions on a test set, you can make predictions
# directly on the `H2OAutoML` object or on the leader model
# object directly
pred <- h2o.predict(aml, x_test)  # predict(aml, test) also works


# Get leaderboard with all possible columns
lb <- h2o.get_leaderboard(object = aml, extra_columns = "ALL")
lb

# # Get the best model using the metric
# m <- aml@leader
# # this is equivalent to
# m <- h2o.get_best_model(aml)
# 
# # Get the best model using a non-default metric
# m <- h2o.get_best_model(aml, criterion = "logloss")
# 
# # Get the best XGBoost model using default sort metric
# xgb <- h2o.get_best_model(aml, algorithm = "xgboost")
# 
# # Get the best XGBoost model, ranked by logloss
# xgb <- h2o.get_best_model(aml, algorithm = "xgboost", criterion = "logloss")


# Auto ML Stacked Ensemble


# Get a specific model by model ID
m <- h2o.getModel("StackedEnsemble_AllModels_AutoML_20220130_081131")

m_pred <- h2o.predict(aml, x_test)
m_pred <- m_pred[, 1]
m_pred <- as.data.frame(m_pred)
m_pred <- as.factor(m_pred$predict)
y_test <- as.data.frame(y_test)
y_test <- as.factor(y_test$`x`)
levels(m_pred) <- levels(y_test)
caret::confusionMatrix(m_pred, y_test)

#########################################################################################
# --------------------------- Stacked Ensemble -----------------------------------------#
#########################################################################################



# Number of CV folds (to generate level-one data for stacking)
nfolds <- 5

# There are a few ways to assemble a list of models to stack together:
# 1. Train individual models and put them in a list
# 2. Train a grid of models
# 3. Train several grids of models
# Note: All base models must have the same cross-validation folds and
# the cross-validated predicted values must be kept.


# 1. Generate a 2-model ensemble (GBM + RF)

# Train & Cross-validate a GBM
my_gbm <- h2o.gbm(x = x,
                  y = y,
                  training_frame = x_train,
                  distribution = "multinomial",
                  ntrees = 500,
                  max_depth = 3,
                  min_rows = 2,
                  learn_rate = 0.01,
                  nfolds = nfolds,
                  keep_cross_validation_predictions = TRUE,
                  seed = 86)

h2o.performance(my_gbm, newdata = x_test)



# Train & Cross-validate a RF
my_rf <- h2o.randomForest(x = x,
                          y = y,
                          training_frame = x_train,
                          ntrees = 500,
                          nfolds = nfolds,
                          keep_cross_validation_predictions = TRUE,
                          seed = 86)

h2o.performance(my_rf, newdata = x_test)


my_glm <- h2o.glm(
  y = y,
  x = x,
  training_frame = x_train,
  family = "multinomial",
  link   = "family_default",
  standardize  = TRUE,
  balance_classes   = FALSE,
  ignore_const_cols = TRUE,
  missing_values_handling = "Skip",
  lambda_search = TRUE,
  solver = "AUTO",
  alpha  = 0.95,
  seed = 86,
  nfolds = nfolds,
  keep_cross_validation_predictions = TRUE,
  model_id = "my_glm"
)

h2o.performance(my_glm, newdata = x_test)


# Train a stacked ensemble using the GBM and RF above
ensemble <- h2o.stackedEnsemble(x = x,
                                y = y,
                                training_frame = x_train,
                                base_models = list(my_gbm, my_rf, my_glm))

# Eval ensemble performance on a test set
perf <- h2o.performance(ensemble, newdata = x_test)

# Compare to base learner performance on the test set
perf_gbm_test <- h2o.performance(my_gbm, newdata = x_test)
perf_rf_test <- h2o.performance(my_rf, newdata = x_test)
baselearner_best_auc_test <- max(h2o.auc(perf_gbm_test), h2o.auc(perf_rf_test))
ensemble_auc_test <- h2o.auc(perf)
print(sprintf("Best Base-learner Test AUC:  %s", baselearner_best_auc_test))
print(sprintf("Ensemble Test AUC:  %s", ensemble_auc_test))
# [1] "Best Base-learner Test AUC:  1"
# [1] "Ensemble Test AUC:  1"

# Generate predictions on a test set 
pred <- h2o.predict(ensemble, newdata = x_test)
pred <- pred[, 1]
pred <- as.data.frame(pred)
pred <- as.factor(pred$predict)
levels(pred) <- levels(y_test)
caret::confusionMatrix(pred, y_test)


my_gbm_pred <- h2o.predict(my_gbm, newdata = x_test)
my_gbm_pred <- my_gbm_pred[, 1]
my_gbm_pred <- as.data.frame(my_gbm_pred)
my_gbm_pred <- as.factor(my_gbm_pred$predict)
levels(my_gbm_pred) <- levels(y_test)
caret::confusionMatrix(my_gbm_pred, y_test)

my_rf_pred <- h2o.predict(my_rf, newdata = x_test)
my_rf_pred <- my_rf_pred[, 1]
my_rf_pred <- as.data.frame(my_rf_pred)
my_rf_pred <- as.factor(my_rf_pred$predict)
levels(my_rf_pred) <- levels(y_test)
caret::confusionMatrix(my_rf_pred, y_test)

my_glm_pred <- h2o.predict(my_glm, newdata = x_test)
my_glm_pred <- my_glm_pred[, 1]
my_glm_pred <- as.data.frame(my_glm_pred)
my_glm_pred <- as.factor(my_glm_pred$predict)
levels(my_glm_pred) <- levels(y_test)
caret::confusionMatrix(my_glm_pred, y_test)


#-------------------#
#  ANOVA PVALUE 500
#-------------------#

library(h2o)

h2o.init()

READ_results_anova_pvalue <- readRDS("C:/Users/juanm/Desktop/Ensemble-classification-high-dimensional-data/Models/READ/READ_results_anova_pvalue.rds")

filtrado_anova_pvalue_500 <-  READ_results_anova_pvalue %>% pull(gen) %>% head(500)
filtrado_anova_pvalue_100 <- READ_results_anova_pvalue %>% pull(gen) %>% head(100)
filtrado_anova_pvalue_50  <- READ_results_anova_pvalue %>% pull(gen) %>% head(50)


# load data with the next sentence
tcga_data <- readRDS("C:/Users/juanm/Desktop/Ensemble-classification-high-dimensional-data/Data/tcga_data_READ.RDS")


# Preprocessing 

clinical_data = colData(tcga_data)
group = factor(clinical_data$definition)
group = relevel(group, ref="Solid Tissue Normal")
design = model.matrix(~group)
head(design)
dge = DGEList( # creating a DGEList object
  counts=assay(tcga_data),
  samples=colData(tcga_data),
  genes=as.data.frame(rowData(tcga_data)))

# filtering
keep = filterByExpr(dge,design) # defining which genes to keep
dge = dge[keep,,keep.lib.sizes=FALSE] # filtering the dge object
rm(keep) #  use rm() to remove objects from memory if you don't need them anymore

dge = calcNormFactors(dge,method="TMM")
v = voom(dge,design,plot=TRUE)

fit = lmFit(v, design)
fit = eBayes(fit)

topGenes = topTable(fit, coef=1, sort.by="p")
print(topGenes)


# Split data (80% train 20% test)

# Transpose and make it into a matrix object
d_mat = as.matrix(t(v$E))
dim(d_mat)

# As before, we want this to be a factor
d_resp = as.factor(v$targets$definition)
levels(d_resp)

df <- cbind(d_mat,d_resp)
df <- as.data.frame(df)
colnames(df)[length(df)] # check class of the tumor 
names(df)[length(df)] <- "Class"
df$Class <- as.factor(df$Class)

# Train/ test split 
set.seed(86)
train <- createDataPartition(y = df$Class, p = 0.7, list = FALSE, times = 1)
x_train <- df[train, ]
x_test  <- df[-train, ]
y_train <- df$Class[train]
y_test  <- df$Class[-train]


x_train <- x_train[c("Class", filtrado_anova_pvalue_500)]
x_test <- x_test[c("Class", filtrado_anova_pvalue_500)]

rm(list = c("df","clinical_data", "d_mat","d_resp","dge", "train",
            "design","fit","group","tcga_data","topGenes","v"))

# create an h20 object 
x_train <- as.h2o(x_train)
x_test <- as.h2o(x_test)
y_train <- as.h2o(y_train)
y_test <- as.h2o(y_test)

# Identify predictors and response
y = "Class"
x <- setdiff(names(x_train), y)


# AutoML: Automatic Machine Learning
aml <- h2o.automl(x = x, y = y,
                  training_frame = x_train,
                  max_models = 10,
                  seed = 86)

# View the AutoML Leaderboard
lb <- aml@leaderboard
print(lb, n = nrow(lb))  # Print all rows instead of default (6 rows)

# The leader model is stored here
aml@leader


# To generate predictions on a test set, you can make predictions
# directly on the `H2OAutoML` object or on the leader model
# object directly
pred <- h2o.predict(aml, x_test)  # predict(aml, test) also works


# Get leaderboard with all possible columns
lb <- h2o.get_leaderboard(object = aml, extra_columns = "ALL")
lb

# # Get the best model using the metric
# m <- aml@leader
# # this is equivalent to
# m <- h2o.get_best_model(aml)
# 
# # Get the best model using a non-default metric
# m <- h2o.get_best_model(aml, criterion = "logloss")
# 
# # Get the best XGBoost model using default sort metric
# xgb <- h2o.get_best_model(aml, algorithm = "xgboost")
# 
# # Get the best XGBoost model, ranked by logloss
# xgb <- h2o.get_best_model(aml, algorithm = "xgboost", criterion = "logloss")


# Auto ML Stacked Ensemble

# 
# # Get a specific model by model ID
m <- h2o.getModel("StackedEnsemble_AllModels_AutoML_20220130_081627")

m_pred <- h2o.predict(aml, x_test)
m_pred <- m_pred[, 1]
m_pred <- as.data.frame(m_pred)
m_pred <- as.factor(m_pred$predict)
y_test <- as.data.frame(y_test)
y_test <- as.factor(y_test$`x`)
levels(m_pred) <- levels(y_test)
caret::confusionMatrix(m_pred, y_test)

#########################################################################################
# --------------------------- Stacked Ensemble -----------------------------------------#
#########################################################################################



# Number of CV folds (to generate level-one data for stacking)
nfolds <- 5

# There are a few ways to assemble a list of models to stack together:
# 1. Train individual models and put them in a list
# 2. Train a grid of models
# 3. Train several grids of models
# Note: All base models must have the same cross-validation folds and
# the cross-validated predicted values must be kept.


# 1. Generate a 2-model ensemble (GBM + RF)

# Train & Cross-validate a GBM
my_gbm <- h2o.gbm(x = x,
                  y = y,
                  training_frame = x_train,
                  distribution = "multinomial",
                  ntrees = 500,
                  max_depth = 3,
                  min_rows = 2,
                  learn_rate = 0.01,
                  nfolds = nfolds,
                  keep_cross_validation_predictions = TRUE,
                  seed = 86)

h2o.performance(my_gbm, newdata = x_test)



# Train & Cross-validate a RF
my_rf <- h2o.randomForest(x = x,
                          y = y,
                          training_frame = x_train,
                          ntrees = 500,
                          nfolds = nfolds,
                          keep_cross_validation_predictions = TRUE,
                          seed = 86)

h2o.performance(my_rf, newdata = x_test)


my_glm <- h2o.glm(
  y = y,
  x = x,
  training_frame = x_train,
  family = "multinomial",
  link   = "family_default",
  standardize  = TRUE,
  balance_classes   = FALSE,
  ignore_const_cols = TRUE,
  missing_values_handling = "Skip",
  lambda_search = TRUE,
  solver = "AUTO",
  alpha  = 0.95,
  seed = 86,
  nfolds = nfolds,
  keep_cross_validation_predictions = TRUE,
  model_id = "my_glm"
)

h2o.performance(my_glm, newdata = x_test)


# Train a stacked ensemble using the GBM and RF above
ensemble <- h2o.stackedEnsemble(x = x,
                                y = y,
                                training_frame = x_train,
                                base_models = list(my_gbm, my_rf, my_glm))

# Eval ensemble performance on a test set
perf <- h2o.performance(ensemble, newdata = x_test)

# Compare to base learner performance on the test set
perf_gbm_test <- h2o.performance(my_gbm, newdata = x_test)
perf_rf_test <- h2o.performance(my_rf, newdata = x_test)
baselearner_best_auc_test <- max(h2o.auc(perf_gbm_test), h2o.auc(perf_rf_test))
ensemble_auc_test <- h2o.auc(perf)
print(sprintf("Best Base-learner Test AUC:  %s", baselearner_best_auc_test))
print(sprintf("Ensemble Test AUC:  %s", ensemble_auc_test))
# [1] "Best Base-learner Test AUC:  1"
# [1] "Ensemble Test AUC:  1"

# Generate predictions on a test set 
pred <- h2o.predict(ensemble, newdata = x_test)
pred <- pred[, 1]
pred <- as.data.frame(pred)
pred <- as.factor(pred$predict)
levels(pred) <- levels(y_test)
caret::confusionMatrix(pred, y_test)


my_gbm_pred <- h2o.predict(my_gbm, newdata = x_test)
my_gbm_pred <- my_gbm_pred[, 1]
my_gbm_pred <- as.data.frame(my_gbm_pred)
my_gbm_pred <- as.factor(my_gbm_pred$predict)
levels(my_gbm_pred) <- levels(y_test)
caret::confusionMatrix(my_gbm_pred, y_test)

my_rf_pred <- h2o.predict(my_rf, newdata = x_test)
my_rf_pred <- my_rf_pred[, 1]
my_rf_pred <- as.data.frame(my_rf_pred)
my_rf_pred <- as.factor(my_rf_pred$predict)
levels(my_rf_pred) <- levels(y_test)
caret::confusionMatrix(my_rf_pred, y_test)

my_glm_pred <- h2o.predict(my_glm, newdata = x_test)
my_glm_pred <- my_glm_pred[, 1]
my_glm_pred <- as.data.frame(my_glm_pred)
my_glm_pred <- as.factor(my_glm_pred$predict)
levels(my_glm_pred) <- levels(y_test)
caret::confusionMatrix(my_glm_pred, y_test)

#-------------------#
#  RFE 50
#-------------------#

library(h2o)

h2o.init()


# load data with the next sentence
tcga_data <- readRDS("C:/Users/juanm/Desktop/Ensemble-classification-high-dimensional-data/Data/tcga_data_READ.RDS")


# Preprocessing 

clinical_data = colData(tcga_data)
group = factor(clinical_data$definition)
group = relevel(group, ref="Solid Tissue Normal")
design = model.matrix(~group)
head(design)
dge = DGEList( # creating a DGEList object
  counts=assay(tcga_data),
  samples=colData(tcga_data),
  genes=as.data.frame(rowData(tcga_data)))

# filtering
keep = filterByExpr(dge,design) # defining which genes to keep
dge = dge[keep,,keep.lib.sizes=FALSE] # filtering the dge object
rm(keep) #  use rm() to remove objects from memory if you don't need them anymore

dge = calcNormFactors(dge,method="TMM")
v = voom(dge,design,plot=TRUE)

fit = lmFit(v, design)
fit = eBayes(fit)

topGenes = topTable(fit, coef=1, sort.by="p")
print(topGenes)


# Split data (80% train 20% test)

# Transpose and make it into a matrix object
d_mat = as.matrix(t(v$E))
dim(d_mat)

# As before, we want this to be a factor
d_resp = as.factor(v$targets$definition)
levels(d_resp)

df <- cbind(d_mat,d_resp)
df <- as.data.frame(df)
colnames(df)[length(df)] # check class of the tumor 
names(df)[length(df)] <- "Class"
df$Class <- as.factor(df$Class)


# Train/ test split 
set.seed(86)
train <- createDataPartition(y = df$Class, p = 0.7, list = FALSE, times = 1)
x_train <- df[train, ]
x_test  <- df[-train, ]
y_train <- df$Class[train]
y_test  <- df$Class[-train]



# load RFE selected Variables 
RFE_results <- readRDS("C:/Users/juanm/Desktop/Ensemble-classification-high-dimensional-data/READ/READ_RFE_results.rds")
selected_vars <- RFE_results$variables
best_vars_50 <- RFE_results$control$functions$selectVar(selected_vars, 50)

x_train <- x_train[, c("Class", best_vars_50)] 
x_test <- x_test[, c("Class", best_vars_50)] 


rm(list = c("df","clinical_data", "d_mat","d_resp","dge", "train",
            "design","fit","group","tcga_data","topGenes","v"))

# create an h20 object 
x_train <- as.h2o(x_train)
x_test <- as.h2o(x_test)
y_train <- as.h2o(y_train)
y_test <- as.h2o(y_test)

# Identify predictors and response
y = "Class"
x <- setdiff(names(x_train), y)


# AutoML: Automatic Machine Learning
aml <- h2o.automl(x = x, y = y,
                  training_frame = x_train,
                  max_models = 10,
                  seed = 86)

# View the AutoML Leaderboard
lb <- aml@leaderboard
print(lb, n = nrow(lb))  # Print all rows instead of default (6 rows)

# The leader model is stored here
aml@leader


# To generate predictions on a test set, you can make predictions
# directly on the `H2OAutoML` object or on the leader model
# object directly
pred <- h2o.predict(aml, x_test)  # predict(aml, test) also works


# Get leaderboard with all possible columns
lb <- h2o.get_leaderboard(object = aml, extra_columns = "ALL")
lb

# # Get the best model using the metric
# m <- aml@leader
# # this is equivalent to
# m <- h2o.get_best_model(aml)
# 
# # Get the best model using a non-default metric
# m <- h2o.get_best_model(aml, criterion = "logloss")
# 
# # Get the best XGBoost model using default sort metric
# xgb <- h2o.get_best_model(aml, algorithm = "xgboost")
# 
# # Get the best XGBoost model, ranked by logloss
# xgb <- h2o.get_best_model(aml, algorithm = "xgboost", criterion = "logloss")


# Auto ML Stacked Ensemble

# 
# # Get a specific model by model ID
m <- h2o.getModel("StackedEnsemble_AllModels_AutoML_20220130_082441")

m_pred <- h2o.predict(aml, x_test)
m_pred <- m_pred[, 1]
m_pred <- as.data.frame(m_pred)
m_pred <- as.factor(m_pred$predict)
y_test <- as.data.frame(y_test)
y_test <- as.factor(y_test$`x`)
levels(m_pred) <- levels(y_test)
caret::confusionMatrix(m_pred, y_test)



#########################################################################################
# --------------------------- Stacked Ensemble -----------------------------------------#
#########################################################################################



# Number of CV folds (to generate level-one data for stacking)
nfolds <- 5

# There are a few ways to assemble a list of models to stack together:
# 1. Train individual models and put them in a list
# 2. Train a grid of models
# 3. Train several grids of models
# Note: All base models must have the same cross-validation folds and
# the cross-validated predicted values must be kept.


# 1. Generate a 2-model ensemble (GBM + RF)

# Train & Cross-validate a GBM
my_gbm <- h2o.gbm(x = x,
                  y = y,
                  training_frame = x_train,
                  distribution = "multinomial",
                  ntrees = 500,
                  max_depth = 3,
                  min_rows = 2,
                  learn_rate = 0.01,
                  nfolds = nfolds,
                  keep_cross_validation_predictions = TRUE,
                  seed = 86)

h2o.performance(my_gbm, newdata = x_test)



# Train & Cross-validate a RF
my_rf <- h2o.randomForest(x = x,
                          y = y,
                          training_frame = x_train,
                          ntrees = 500,
                          nfolds = nfolds,
                          keep_cross_validation_predictions = TRUE,
                          seed = 86)

h2o.performance(my_rf, newdata = x_test)


my_glm <- h2o.glm(
  y = y,
  x = x,
  training_frame = x_train,
  family = "multinomial",
  link   = "family_default",
  standardize  = TRUE,
  balance_classes   = FALSE,
  ignore_const_cols = TRUE,
  missing_values_handling = "Skip",
  lambda_search = TRUE,
  solver = "AUTO",
  alpha  = 0.95,
  seed = 86,
  nfolds = nfolds,
  keep_cross_validation_predictions = TRUE,
  model_id = "my_glm"
)

h2o.performance(my_glm, newdata = x_test)


# Train a stacked ensemble using the GBM and RF above
ensemble <- h2o.stackedEnsemble(x = x,
                                y = y,
                                training_frame = x_train,
                                base_models = list(my_gbm, my_rf, my_glm))

# Eval ensemble performance on a test set
perf <- h2o.performance(ensemble, newdata = x_test)

# Compare to base learner performance on the test set
perf_gbm_test <- h2o.performance(my_gbm, newdata = x_test)
perf_rf_test <- h2o.performance(my_rf, newdata = x_test)
baselearner_best_auc_test <- max(h2o.auc(perf_gbm_test), h2o.auc(perf_rf_test))
ensemble_auc_test <- h2o.auc(perf)
print(sprintf("Best Base-learner Test AUC:  %s", baselearner_best_auc_test))
print(sprintf("Ensemble Test AUC:  %s", ensemble_auc_test))
# [1] "Best Base-learner Test AUC:  1"
# [1] "Ensemble Test AUC:  1"

# Generate predictions on a test set 
pred <- h2o.predict(ensemble, newdata = x_test)
pred <- pred[, 1]
pred <- as.data.frame(pred)
pred <- as.factor(pred$predict)
levels(pred) <- levels(y_test)
caret::confusionMatrix(pred, y_test)

my_gbm_pred <- h2o.predict(my_gbm, newdata = x_test)
my_gbm_pred <- my_gbm_pred[, 1]
my_gbm_pred <- as.data.frame(my_gbm_pred)
my_gbm_pred <- as.factor(my_gbm_pred$predict)
levels(my_gbm_pred) <- levels(y_test)
caret::confusionMatrix(my_gbm_pred, y_test)

my_rf_pred <- h2o.predict(my_rf, newdata = x_test)
my_rf_pred <- my_rf_pred[, 1]
my_rf_pred <- as.data.frame(my_rf_pred)
my_rf_pred <- as.factor(my_rf_pred$predict)
levels(my_rf_pred) <- levels(y_test)
caret::confusionMatrix(my_rf_pred, y_test)

my_glm_pred <- h2o.predict(my_glm, newdata = x_test)
my_glm_pred <- my_glm_pred[, 1]
my_glm_pred <- as.data.frame(my_glm_pred)
my_glm_pred <- as.factor(my_glm_pred$predict)
levels(my_glm_pred) <- levels(y_test)
caret::confusionMatrix(my_glm_pred, y_test)


#-------------------#
#  RFE 100
#-------------------#

library(h2o)

h2o.init()

# load data with the next sentence
tcga_data <- readRDS("C:/Users/juanm/Desktop/Ensemble-classification-high-dimensional-data/Data/tcga_data_READ.RDS")


# Preprocessing 

clinical_data = colData(tcga_data)
group = factor(clinical_data$definition)
group = relevel(group, ref="Solid Tissue Normal")
design = model.matrix(~group)
head(design)
dge = DGEList( # creating a DGEList object
  counts=assay(tcga_data),
  samples=colData(tcga_data),
  genes=as.data.frame(rowData(tcga_data)))

# filtering
keep = filterByExpr(dge,design) # defining which genes to keep
dge = dge[keep,,keep.lib.sizes=FALSE] # filtering the dge object
rm(keep) #  use rm() to remove objects from memory if you don't need them anymore

dge = calcNormFactors(dge,method="TMM")
v = voom(dge,design,plot=TRUE)

fit = lmFit(v, design)
fit = eBayes(fit)

topGenes = topTable(fit, coef=1, sort.by="p")
print(topGenes)


# Split data (80% train 20% test)

# Transpose and make it into a matrix object
d_mat = as.matrix(t(v$E))
dim(d_mat)

# As before, we want this to be a factor
d_resp = as.factor(v$targets$definition)
levels(d_resp)

df <- cbind(d_mat,d_resp)
df <- as.data.frame(df)
colnames(df)[length(df)] # check class of the tumor 
names(df)[length(df)] <- "Class"
df$Class <- as.factor(df$Class)


# Train/ test split 
set.seed(86)
train <- createDataPartition(y = df$Class, p = 0.7, list = FALSE, times = 1)
x_train <- df[train, ]
x_test  <- df[-train, ]
y_train <- df$Class[train]
y_test  <- df$Class[-train]


# load RFE selected Variables 
RFE_results <- readRDS("C:/Users/juanm/Desktop/Ensemble-classification-high-dimensional-data/READ/READ_RFE_results.rds")
selected_vars <- RFE_results$variables
best_vars_100 <- RFE_results$control$functions$selectVar(selected_vars, 100)

x_train <- x_train[, c("Class", best_vars_100)] 
x_test <- x_test[, c("Class", best_vars_100)] 


rm(list = c("df","clinical_data", "d_mat","d_resp","dge", "train",
            "design","fit","group","tcga_data","topGenes","v"))

# create an h20 object 
x_train <- as.h2o(x_train)
x_test <- as.h2o(x_test)
y_train <- as.h2o(y_train)
y_test <- as.h2o(y_test)

# Identify predictors and response
y = "Class"
x <- setdiff(names(x_train), y)


# AutoML: Automatic Machine Learning
aml <- h2o.automl(x = x, y = y,
                  training_frame = x_train,
                  max_models = 10,
                  seed = 86)

# View the AutoML Leaderboard
lb <- aml@leaderboard
print(lb, n = nrow(lb))  # Print all rows instead of default (6 rows)

# The leader model is stored here
aml@leader


# To generate predictions on a test set, you can make predictions
# directly on the `H2OAutoML` object or on the leader model
# object directly
pred <- h2o.predict(aml, x_test)  # predict(aml, test) also works


# Get leaderboard with all possible columns
lb <- h2o.get_leaderboard(object = aml, extra_columns = "ALL")
lb

# # Get the best model using the metric
# m <- aml@leader
# # this is equivalent to
# m <- h2o.get_best_model(aml)
# 
# # Get the best model using a non-default metric
# m <- h2o.get_best_model(aml, criterion = "logloss")
# 
# # Get the best XGBoost model using default sort metric
# xgb <- h2o.get_best_model(aml, algorithm = "xgboost")
# 
# # Get the best XGBoost model, ranked by logloss
# xgb <- h2o.get_best_model(aml, algorithm = "xgboost", criterion = "logloss")


# Auto ML Stacked Ensemble

# 
# # Get a specific model by model ID
m <- h2o.getModel("StackedEnsemble_AllModels_AutoML_20220130_082822")

m_pred <- h2o.predict(aml, x_test)
m_pred <- m_pred[, 1]
m_pred <- as.data.frame(m_pred)
m_pred <- as.factor(m_pred$predict)
y_test <- as.data.frame(y_test)
y_test <- as.factor(y_test$`x`)
levels(m_pred) <- levels(y_test)
caret::confusionMatrix(m_pred, y_test)

#########################################################################################
# --------------------------- Stacked Ensemble -----------------------------------------#
#########################################################################################



# Number of CV folds (to generate level-one data for stacking)

nfolds <- 5


# Train & Cross-validate a GBM
my_gbm <- h2o.gbm(x = x,
                  y = y,
                  training_frame = x_train,
                  distribution = "multinomial",
                  ntrees = 500,
                  max_depth = 3,
                  min_rows = 2,
                  learn_rate = 0.01,
                  nfolds = nfolds,
                  keep_cross_validation_predictions = TRUE,
                  seed = 86)

h2o.performance(my_gbm, newdata = x_test)



# Train & Cross-validate a RF
my_rf <- h2o.randomForest(x = x,
                          y = y,
                          training_frame = x_train,
                          ntrees = 500,
                          nfolds = nfolds,
                          keep_cross_validation_predictions = TRUE,
                          seed = 86)

h2o.performance(my_rf, newdata = x_test)


my_glm <- h2o.glm(
  y = y,
  x = x,
  training_frame = x_train,
  family = "multinomial",
  link   = "family_default",
  standardize  = TRUE,
  balance_classes   = FALSE,
  ignore_const_cols = TRUE,
  missing_values_handling = "Skip",
  lambda_search = TRUE,
  solver = "AUTO",
  alpha  = 0.95,
  seed = 86,
  nfolds = nfolds,
  keep_cross_validation_predictions = TRUE,
  model_id = "my_glm"
)

h2o.performance(my_glm, newdata = x_test)



# Train a stacked ensemble using the GBM and RF above
ensemble <- h2o.stackedEnsemble(x = x,
                                y = y,
                                training_frame = x_train,
                                base_models = list(my_gbm, my_rf, my_glm))

# Eval ensemble performance on a test set
perf <- h2o.performance(ensemble, newdata = x_test)

# Compare to base learner performance on the test set
perf_gbm_test <- h2o.performance(my_gbm, newdata = x_test)
perf_rf_test <- h2o.performance(my_rf, newdata = x_test)
baselearner_best_auc_test <- max(h2o.auc(perf_gbm_test), h2o.auc(perf_rf_test))
ensemble_auc_test <- h2o.auc(perf)
print(sprintf("Best Base-learner Test AUC:  %s", baselearner_best_auc_test))
print(sprintf("Ensemble Test AUC:  %s", ensemble_auc_test))
# [1] "Best Base-learner Test AUC:  1"
# [1] "Ensemble Test AUC:  1"

# Generate predictions on a test set 
pred <- h2o.predict(ensemble, newdata = x_test)
pred <- pred[, 1]
pred <- as.data.frame(pred)
pred <- as.factor(pred$predict)
levels(pred) <- levels(y_test)
caret::confusionMatrix(pred, y_test)

my_gbm_pred <- h2o.predict(my_gbm, newdata = x_test)
my_gbm_pred <- my_gbm_pred[, 1]
my_gbm_pred <- as.data.frame(my_gbm_pred)
my_gbm_pred <- as.factor(my_gbm_pred$predict)
levels(my_gbm_pred) <- levels(y_test)
caret::confusionMatrix(my_gbm_pred, y_test)

my_rf_pred <- h2o.predict(my_rf, newdata = x_test)
my_rf_pred <- my_rf_pred[, 1]
my_rf_pred <- as.data.frame(my_rf_pred)
my_rf_pred <- as.factor(my_rf_pred$predict)
levels(my_rf_pred) <- levels(y_test)
caret::confusionMatrix(my_rf_pred, y_test)

my_glm_pred <- h2o.predict(my_glm, newdata = x_test)
my_glm_pred <- my_glm_pred[, 1]
my_glm_pred <- as.data.frame(my_glm_pred)
my_glm_pred <- as.factor(my_glm_pred$predict)
levels(my_glm_pred) <- levels(y_test)
caret::confusionMatrix(my_glm_pred, y_test)

# Close the cluster 

h2o.shutdown()
