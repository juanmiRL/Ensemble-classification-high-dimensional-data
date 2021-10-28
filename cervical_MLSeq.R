if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.12")

BiocManager::install("MLSeq")
BiocManager::install("S4Vectors")
BiocManager::install("gpls")
BiocManager::install("TCGAbiolinks")

library(tidyverse)
library(MLSeq)
library(S4Vectors)
library(DESeq2)

########################################################################
####################### Cervical Data ##################################
########################################################################


filepath <- system.file("extdata/cervical.txt", package = "MLSeq")
cervical <- read.table(filepath, header=TRUE)
head(cervical[ ,1:10])
class <- DataFrame(condition = factor(rep(c("N","T"), c(29, 29))))


# We do not perform a differential expression analysis to select differentially
# expressed genes. However, in practice, DE analysis might be performed before
# fitting classifiers. Here, we selected top 100 features having the highest
# gene-wise variances in order to decrease computational cost.

set.seed(86)
vars <- sort(apply(cervical, 1, var, na.rm = TRUE), decreasing = TRUE)
data <- cervical[names(vars)[1:100], ]
nTest <- ceiling(ncol(data) * 0.3)
ind <- sample(ncol(data), nTest, FALSE)

# Minimum count is set to 1 in order to prevent 0 division problem within classification models.

data.train <- as.matrix(data[ ,-ind] + 1)
data.test <- as.matrix(data[ ,ind] + 1)
classtr <- DataFrame(condition = class[-ind, ])
classts <- DataFrame(condition = class[ind, ])

data.trainS4 = DESeqDataSetFromMatrix(countData = data.train, colData = classtr,
                                      design = formula(~condition))
data.testS4 = DESeqDataSetFromMatrix(countData = data.test, colData = classts,
                                     design = formula(~condition))

printAvailableMethods()


# The possible normalization-transformation combinations are:
# - deseq-vst: Normalization is applied with deseq median ratio method. Variance stabilizing transformation is applied to the normalized data
# - deseq-rlog: Normalization is applied with deseq median ratio method. Regularized logarithmic
# transformation is applied to the normalized data
# - deseq-logcpm: Normalization is applied with deseq median ratio method. Log of counts-per-million
# transformation is applied to the normalized data
# - tmm-logcpm: Normalization is applied with trimmed mean of M values (TMM) method. Log of
# counts-per-million transformation is applied to the normalized data.


# Support Vector Machines with Radial Kernel
svm_radial_fit <- classify(data = data.trainS4, method = "svmRadial",
                           preProcessing = "deseq-rlog", ref = "T", tuneLength = 10,
                           control = trainControl(method = "repeatedcv", number = 5,
                                                  repeats = 10, classProbs = TRUE))
show(svm_radial_fit)
trained(svm_radial_fit)
plot(svm_radial_fit)

pred_radial_svm <- predictClassify(svm_radial_fit, data.testS4)
pred_radial_svm
pred_radial_svm <- relevel(pred_radial_svm, ref = "T")
actual <- relevel(classts$condition, ref = "T")
tbl <- table(Predicted = pred_radial_svm, Actual = actual)
confusionMatrix(tbl, positive = "T")

# Support Vector Machines with Linear Kernel
svm_linear_fit <- classify(data = data.trainS4, method = "svmLinear",
                           preProcessing = "deseq-rlog", ref = "T", tuneLength = 10,
                           control = trainControl(method = "repeatedcv", number = 5,
                                                  repeats = 10, classProbs = TRUE))
show(svm_linear_fit)
trained(svm_linear_fit)

pred_linear_svm <- predict(svm_linear_fit, data.testS4)
pred_linear_svm
pred_linear_svm <- relevel(pred_linear_svm, ref = "T")
actual <- relevel(classts$condition, ref = "T")
tbl <- table(Predicted = pred_linear_svm, Actual = actual)
confusionMatrix(tbl, positive = "T")

# Support Vector Machines with Polynomial Kernel
svm_poly_fit <- classify(data = data.trainS4, method = "svmPoly",
                         preProcessing = "deseq-rlog", ref = "T", tuneLength = 10,
                         control = trainControl(method = "repeatedcv", number = 5,
                                                repeats = 10, classProbs = TRUE))
show(svm_poly_fit)
trained(svm_poly_fit)
plot(svm_poly_fit)

pred_poly_svm <- predictClassify(svm_poly_fit, data.testS4)
pred_poly_svm
pred_poly_svm <- relevel(pred_poly_svm, ref = "T")
actual <- relevel(classts$condition, ref = "T")
tbl <- table(Predicted = pred_poly_svm, Actual = actual)
confusionMatrix(tbl, positive = "T")

# Random Forest 
rf_fit <- classify(data = data.trainS4, method = "rf",
                   preProcessing = "deseq-vst", ref = "T", tuneLength = 10,
                   control = trainControl(method = "repeatedcv", number = 5,
                                          repeats = 10, classProbs = TRUE))
show(rf_fit)
trained(rf_fit)
plot(rf_fit)

pred_rf <- predictClassify(rf_fit, data.testS4)
pred_rf
pred_rf <- relevel(pred_rf, ref = "T")
actual <- relevel(classts$condition, ref = "T")
tbl <- table(Predicted = pred_rf, Actual = actual)
confusionMatrix(tbl, positive = "T")

# Stochastic Gradient Boosting
gbm_fit <- classify(data = data.trainS4, method = "gbm",
                    preProcessing = "deseq-vst", ref = "T", tuneLength = 10,
                    control = trainControl(method = "repeatedcv", number = 5,
                                           repeats = 10, classProbs = TRUE))

pred_gbm <- predictClassify(gbm_fit, data.testS4)
pred_gbm
pred_gbm <- relevel(pred_gbm, ref = "T")
actual <- relevel(classts$condition, ref = "T")
tbl <- table(Predicted = pred_gbm, Actual = actual)
confusionMatrix(tbl, positive = "T")


# Neural Network
nnet_fit <- classify(data = data.trainS4, method = "nnet",
                     preProcessing = "deseq-vst", ref = "T", tuneLength = 10,
                     control = trainControl(method = "repeatedcv", number = 5,
                                            repeats = 10, classProbs = TRUE))

show(nnet_fit)
trained(nnet_fit)
plot(nnet_fit)

pred_nnet <- predictClassify(nnet_fit, data.testS4)
pred_nnet
pred_nnet <- relevel(pred_nnet, ref = "T")
actual <- relevel(classts$condition, ref = "T")
tbl <- table(Predicted = pred_nnet, Actual = actual)
confusionMatrix(tbl, positive = "T")


# glmnet
glmnet_fit <- classify(data = data.trainS4, method = "glmnet",
                       preProcessing = "deseq-vst", ref = "T", tuneLength = 10,
                       control = trainControl(method = "repeatedcv", number = 5,
                                              repeats = 10, classProbs = TRUE))

show(glmnet_fit)
trained(glmnet_fit)
plot(glmnet_fit)

pred_glmnet <- predictClassify(glmnet_fit, data.testS4)
pred_glmnet
pred_glmnet <- relevel(pred_glmnet, ref = "T")
actual <- relevel(classts$condition, ref = "T")
tbl <- table(Predicted = pred_glmnet, Actual = actual)
confusionMatrix(tbl, positive = "T")


# High Dimensional Discriminant Analysis
hdda_fit <- classify(data = data.trainS4, method = "hdda",
                     preProcessing = "deseq-vst", ref = "T", tuneLength = 10,
                     control = trainControl(method = "repeatedcv", number = 5,
                                            repeats = 10, classProbs = TRUE))

show(hdda_fit)
trained(hdda_fit)
plot(hdda_fit)

pred_hdda <- predictClassify(hdda_fit, data.testS4)
pred_hdda
pred_hdda <- relevel(pred_hdda, ref = "T")
actual <- relevel(classts$condition, ref = "T")
tbl <- table(Predicted = pred_hdda, Actual = actual)
confusionMatrix(tbl, positive = "T")


# Generalized Partial Least Squares

gpls_fit <- classify(data = data.trainS4, method = "gpls",
                     preProcessing = "deseq-vst", ref = "T", tuneLength = 10,
                     control = trainControl(method = "repeatedcv", number = 5,
                                            repeats = 10, classProbs = TRUE))

show(gpls_fit)
trained(gpls_fit)
plot(gpls_fit)

pred_gpls <- predictClassify(gpls_fit, data.testS4)
pred_gpls
pred_gpls <- relevel(pred_gpls, ref = "T")
actual <- relevel(classts$condition, ref = "T")
tbl <- table(Predicted = pred_gpls, Actual = actual)
confusionMatrix(tbl, positive = "T")


# Boosted Logistic Regression

LogitBoost_fit <- classify(data = data.trainS4, method = "LogitBoost",
                           preProcessing = "deseq-vst", ref = "T", tuneLength = 10,
                           control = trainControl(method = "repeatedcv", number = 5,
                                                  repeats = 10, classProbs = TRUE))

show(LogitBoost_fit)
trained(LogitBoost_fit)
plot(LogitBoost_fit)

pred_LogitBoost <- predictClassify(LogitBoost_fit, data.testS4)
pred_LogitBoost
pred_LogitBoost <- relevel(pred_LogitBoost, ref = "T")
actual <- relevel(classts$condition, ref = "T")
tbl <- table(Predicted = pred_LogitBoost, Actual = actual)
confusionMatrix(tbl, positive = "T")
