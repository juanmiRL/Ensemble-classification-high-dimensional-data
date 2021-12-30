#___________________________________________________________
#                  Rectum Adenocarcinoma 
#___________________________________________________________

library(tidyverse)
library(MLSeq)
library(S4Vectors)
library(DESeq2)
library(kernlab)
library(e1071)
library(READ)
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
library(ROCit)
library(precrec)
library(pROC)
library(Epi)

# load data with the next sentence
tcga_data = readRDS(file = "tcga_data_READ.RDS")


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

# drop level with a single value 
f <- which(df$Class == 2)
f <- df[-156,]
df <- droplevels(f)
levels(df$Class)

# write.csv(df, "df_STAD.csv")

# Set (random-number-generator) seed so that results are consistent between runs
set.seed(86)
train <- createDataPartition(y = df$Class, p = 0.7, list = FALSE, times = 1)
x_train <- df[train, ]
x_test  <- df[-train, ]
y_train = df$Class[train]
y_test  = df$Class[-train]



# The dataframes df, clinical_data, d_mat, d_res are
# eliminated so as not to occupy so much memory
rm(list = c("df","clinical_data", "d_mat","d_resp","dge", "train",
            "design","fit","group","tcga_data","topGenes","v"))


prop.table(table(x_train$Class)) %>% round(3)
prop.table(table(x_test$Class)) %>% round(3)

summary(x_train$ENSG00000000003)
head(x_train[,1:6])


# Representation of the expression of 100 randomly selected genes.
set.seed(86)
x_train %>% select_at(sample(4:ncol(x_train), 100)) %>%
  gather(key = "gen", value = "expresion") %>%
  ggplot(aes(x = gen, y = expresion)) +
  geom_boxplot(outlier.size = 0.3, fill = "gray70") +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

custom_anova <- function(x,y){
  anova <- summary(aov(x ~ as.factor(y)))
  return(unlist(anova)["Pr(>F)1"])
}

p_values <- x_train %>%
  dplyr::select(-Class) %>%
  map_dbl(.f = custom_anova, y = x_train$Class) %>%
  sort() 

p_values %>% head(10)


#===============================================================================
# Parallel version of bootstrapping for ANOVA filter
#===============================================================================

n_cores <- parallel::detectCores() - 1
registerDoParallel(makeCluster(n_cores))
getDoParWorkers()


# number of iterations bootstrapping
n_boot <- 100

# seeds for reproducibility
set.seed(86)
seeds = sample.int(1000, size = n_boot)


results_anova_pvalue <- foreach(i = 1:n_boot) %dopar% {
  library(dplyr)
  library(purrr)
  # bootstrapping sample is created
  set.seed(seeds[i])
  index <- sample(1:nrow(x_train), size = nrow(x_train), replace = TRUE)
  pseudo_sample <- x_train[index, ]
  
  # p-values for the new sample are computed
  p_values <- dplyr::select(pseudo_sample,-Class) %>%
    purrr::map_dbl(.f = custom_anova, y = pseudo_sample$Class) 
  
  # return p-values
  p_values
}

options(cores = 1)

# The results saved in a list are transformed to data frame
names(results_anova_pvalue) <-  paste("resample", 1:n_boot, sep = "_") 
results_anova_pvalue <- data.frame(results_anova_pvalue)
results_anova_pvalue <- results_anova_pvalue %>% rownames_to_column(var = "gen")
results_anova_pvalue <- results_anova_pvalue %>%
  mutate(pvalue_mean = rowMeans(results_anova_pvalue[, -1])) %>%
  arrange(pvalue_mean)

# The created object is saved to disk so as not to have to repeat the entire computing
saveRDS(object = results_anova_pvalue, file = "READ_results_anova_pvalue.rds")

results_anova_pvalue %>% dplyr::select(1,2,3,4) %>% head()

# The 50, 100 and 500 genes identified as most relevant are filtered by anova

READ_filter_anova_pvalue_50 <-  results_anova_pvalue %>% pull(gen) %>% head(50)
READ_filter_anova_pvalue_100 <- results_anova_pvalue %>% pull(gen) %>% head(100)
READ_filter_anova_pvalue_500 <- results_anova_pvalue %>% pull(gen) %>% head(500)


#*******************************************************************************
#*************************** ANOVA P-VALUE 50 features *************************
#*******************************************************************************


#===============================================================================
#=========================== RANDOM FOREST ===================================== 
#===============================================================================


# PARALLEL PROCESS
#===============================================================================

#install.packages("doMC", repos="http://R-Forge.R-project.org")

registerDoMC(cores = 11)

# HYPERPARAMETERS, NUMBER OF REPETITIONS AND SEEDS FOR EACH REPEAT
#===============================================================================
repetitions_boot <- 50

# Hyperparameters
hyperparameters <- expand.grid(mtry = c(2, 5, 10, 50),
                               min.node.size = c(2, 3, 4, 5, 10),
                               splitrule = "gini")

set.seed(86)
seeds <- vector(mode = "list", length = repetitions_boot + 1)
for (i in 1:repetitions_boot) {
  seeds[[i]] <- sample.int(1000, nrow(hyperparameters))
}
seeds[[repetitions_boot + 1]] <- sample.int(1000, 1)

# DEFINITION OF TRAINING
#===============================================================================
control_train <- trainControl(method = "boot", number = repetitions_boot,
                              seeds = seeds, returnResamp = "final",
                              verboseIter = TRUE, allowParallel = TRUE)

# FIT MODEL 
# ==============================================================================
x_train$Class <- as.factor(x_train$Class)

set.seed(86)
READ_rf_pvalue_50 <- caret::train(
  form = Class ~ .,
  data = x_train[c("Class", READ_filter_anova_pvalue_50)],
  method = "ranger",
  tuneGrid = hyperparameters,
  metric = "Accuracy",
  trControl = control_train,
  num.trees = 500)

saveRDS(object = READ_rf_pvalue_50, file = "READ_rf_pvalue_50.rds")
registerDoMC(cores = 1)


READ_rf_pvalue_50

# GRAPHIC REPRESENTATION
# ==============================================================================
ggplot(READ_rf_pvalue_50, highlight = TRUE) +
  labs(title = "Evolution of the accuracy of the Random Forest model") +
  guides(color = guide_legend(title = "mtry"),
         shape = guide_legend(title = "mtry")) +
  theme_bw()

# TEST PREDICTIONS
# ==============================================================================
predic_READ_rf_pvalue_50 <- predict(object = READ_rf_pvalue_50, newdata = x_test)
predic_READ_rf_pvalue_50
y_test <- as.factor(y_test)
caret::confusionMatrix(predic_READ_rf_pvalue_50, y_test)


#===============================================================================
#=========================== SVM Kernel Radial =================================
#===============================================================================


# PARALLEL PROCESS
#===============================================================================

#install.packages("doMC", repos="http://R-Forge.R-project.org")

registerDoMC(cores = 11)

# HYPERPARAMETERS, NUMBER OF REPETITIONS AND SEEDS FOR EACH REPEAT
#===============================================================================
repetitions_boot <- 50

# Hyperparameters
hyperparameters <- expand.grid(sigma = c(0.0001, 0.001, 0.01),
                               C = c(1, 10, 50, 100, 250, 500, 700, 1000))

set.seed(86)
seeds <- vector(mode = "list", length = repetitions_boot + 1)
for (i in 1:repetitions_boot) {
  seeds[[i]] <- sample.int(1000, nrow(hyperparameters))
}
seeds[[repetitions_boot + 1]] <- sample.int(1000, 1)

# DEFINITION OF TRAINING
#===============================================================================
control_train <- trainControl(method = "boot", number = repetitions_boot,
                              seeds = seeds, returnResamp = "final",
                              verboseIter = TRUE, allowParallel = TRUE)

# FIT MODEL 
# ==============================================================================

set.seed(86)
READ_svmrad_pvalue_50 <- caret::train(
  form = Class ~ .,
  data = x_train[c("Class", READ_filter_anova_pvalue_50)],
  method = "svmRadial",
  tuneGrid = hyperparameters,
  metric = "Accuracy",
  trControl = control_train
)

saveRDS(object = READ_svmrad_pvalue_50, file = "READ_svmrad_pvalue_50.rds")
registerDoMC(cores = 1)


READ_svmrad_pvalue_50

# GRAPHIC REPRESENTATION
# ==============================================================================
ggplot(READ_svmrad_pvalue_50, highlight = TRUE) +
  labs(title = "Evolution of the accuracy of the SVM model") +
  guides(color = guide_legend(title = "Sigma"),
         shape = guide_legend(title = "Sigma")) +
  theme_bw()

# TEST PREDICTIONS
# ==============================================================================
predic_READ_svm_pvalue_50 <- predict(object = READ_svmrad_pvalue_50, newdata = x_test)
predic_READ_svm_pvalue_50
y_test <- as.factor(y_test)
caret::confusionMatrix(predic_READ_svm_pvalue_50, y_test)


#===============================================================================
#=========================== Neural Network ====================================
#===============================================================================


# PARALLEL PROCESS
#===============================================================================

#install.packages("doMC", repos="http://R-Forge.R-project.org")

registerDoMC(cores = 11)

# HYPERPARAMETERS, NUMBER OF REPETITIONS AND SEEDS FOR EACH REPEAT
#===============================================================================
repetitions_boot <- 50

# Hyperparameters
hyperparameters <- expand.grid(size = c(2, 4, 8, 10, 16),
                               decay = c(0.01, 0.1))

set.seed(86)
seeds <- vector(mode = "list", length = repetitions_boot + 1)
for (i in 1:repetitions_boot) {
  seeds[[i]] <- sample.int(1000, nrow(hyperparameters))
}
seeds[[repetitions_boot + 1]] <- sample.int(1000, 1)

# DEFINITION OF TRAINING
#===============================================================================
control_train <- trainControl(method = "boot", number = repetitions_boot,
                              seeds = seeds, returnResamp = "final",
                              verboseIter = TRUE, allowParallel = TRUE)

# FIT MODEL 
# ==============================================================================

set.seed(86)
READ_nnet_pvalue_50 <- caret::train(
  form = Class ~ .,
  data = x_train[c("Class", READ_filter_anova_pvalue_50)],
  method = "nnet",
  tuneGrid = hyperparameters,
  metric = "Accuracy",
  trControl = control_train,
  rang = c(-0.7, 0.7),
  MaxNWts = 21000,
  trace = FALSE
)

saveRDS(object = READ_nnet_pvalue_50, file = "READ_nnet_pvalue_50.rds")
registerDoMC(cores = 1)

READ_nnet_pvalue_50

# GRAPHIC REPRESENTATION
# ==============================================================================
ggplot(READ_nnet_pvalue_50, highlight = TRUE) +
  labs(title = "Evolution of the accuracy of the Neural Net model") +
  guides(color = guide_legend(title = "Decay"),
         shape = guide_legend(title = "Decay")) +
  theme_bw()

# TEST PREDICTIONS
# ==============================================================================
predic_READ_nnet_pvalue_50 <- predict(object = READ_nnet_pvalue_50, newdata = x_test)
predic_READ_nnet_pvalue_50
y_test <- as.factor(y_test)
caret::confusionMatrix(predic_READ_nnet_pvalue_50, y_test)



#===============================================================================
#=========================== Gradient boosting =================================
#===============================================================================


# PARALLEL PROCESS
#===============================================================================

#install.packages("doMC", repos="http://R-Forge.R-project.org")

registerDoMC(cores = 11)

# HYPERPARAMETERS, NUMBER OF REPETITIONS AND SEEDS FOR EACH REPEAT
#===============================================================================
repetitions_boot <- 50

# Hyperparameters
hyperparameters <- expand.grid(
  shrinkage = c(0.01, 0.1, 0.3),
  interaction.depth = c(1, 3, 5),
  n.minobsinnode = c(5, 10, 15),
  n.trees = c(500, 1000)
  # bag.fraction = c(0.65, 0.8, 1), 
  # optimal_trees = 0,               
  # min_RMSE = 0,                    
  # min_cor = 0
)

set.seed(86)
seeds <- vector(mode = "list", length = repetitions_boot + 1)
for (i in 1:repetitions_boot) {
  seeds[[i]] <- sample.int(1000, nrow(hyperparameters))
}
seeds[[repetitions_boot + 1]] <- sample.int(1000, 1)

# DEFINITION OF TRAINING
#===============================================================================
control_train <- trainControl(method = "boot", number = repetitions_boot,
                              seeds = seeds, returnResamp = "final",
                              verboseIter = TRUE, allowParallel = TRUE)

# FIT MODEL 
# ==============================================================================

set.seed(86)
READ_gbm_pvalue_50 <- caret::train(
  form = Class ~ .,
  data = x_train[c("Class", READ_filter_anova_pvalue_50)],
  method = "gbm",
  tuneGrid = hyperparameters,
  metric = "Accuracy",
  preProcess = c("center", "scale"),
  trControl = control_train
)

saveRDS(object = READ_gbm_pvalue_50, file = "READ_gbm_pvalue_50.rds")
registerDoMC(cores = 1)

READ_READ_pvalue_50

# GRAPHIC REPRESENTATION
# ==============================================================================
ggplot(READ_gbm_pvalue_50, highlight = TRUE) +
  labs(title = "Evolution of the accuracy of the READ model") +
  guides(color = guide_legend(title = "Size"),
         shape = guide_legend(title = "Size")) +
  theme_bw()

# TEST PREDICTIONS
# ==============================================================================
predic_READ_gbm_pvalue_50 <- predict(object = READ_gbm_pvalue_50, newdata = x_test)
predic_READ_gbm_pvalue_50
y_test <- as.factor(y_test)
caret::confusionMatrix(predic_READ_gbm_pvalue_50, y_test)


#===============================================================================
#=========================== XGBM ==============================================
#===============================================================================


# PARALLEL PROCESS
#===============================================================================

#install.packages("doMC", repos="http://R-Forge.R-project.org")

registerDoMC(cores = 11)

# HYPERPARAMETERS, NUMBER OF REPETITIONS AND SEEDS FOR EACH REPEAT
#===============================================================================
repetitions_boot <- 50

# Hyperparameters
hyperparameters <- expand.grid(
  nrounds = c(500, 1000),
  eta = c(0.01, 0.001), 
  max_depth = c(2, 4, 6),
  gamma = 1,
  colsample_bytree = c(0.2, 0.4),
  min_child_weight = c(1, 5),
  subsample = 1
)


set.seed(86)
seeds <- vector(mode = "list", length = repetitions_boot + 1)
for (i in 1:repetitions_boot) {
  seeds[[i]] <- sample.int(1000, nrow(hyperparameters))
}
seeds[[repetitions_boot + 1]] <- sample.int(1000, 1)

# DEFINITION OF TRAINING
#===============================================================================
control_train <- trainControl(method = "boot", number = repetitions_boot,
                              seeds = seeds, returnResamp = "final",
                              verboseIter = TRUE, allowParallel = TRUE)

# FIT MODEL 
# ==============================================================================

set.seed(86)
READ_xgbm_pvalue_50 <- caret::train(
  form = Class ~ .,
  data = x_train[c("Class", READ_filter_anova_pvalue_50)],
  method = "xgbTree",
  tuneGrid = hyperparameters,
  metric = "Accuracy",
  preProcess = c("center", "scale"),
  trControl = control_train
)

saveRDS(object = READ_xgbm_pvalue_50, file = "READ_xgbm_pvalue_50.rds")
registerDoMC(cores = 1)

READ_xgbm_pvalue_50

# GRAPHIC REPRESENTATION
# ==============================================================================
plot(READ_xgbm_pvalue_50, highlight = TRUE)

# TEST PREDICTIONS
# ==============================================================================
predic_READ_xgbm_pvalue_50 <- predict(object = READ_xgbm_pvalue_50, newdata = x_test)
predic_READ_xgbm_pvalue_50
y_test <- as.factor(y_test)
caret::confusionMatrix(predic_READ_xgbm_pvalue_50, y_test)


#===============================================================================
#=========================== glmnet ============================================
#===============================================================================


# PARALLEL PROCESS
#===============================================================================

#install.packages("doMC", repos="http://R-Forge.R-project.org")

registerDoMC(cores = 11)

# HYPERPARAMETERS, NUMBER OF REPETITIONS AND SEEDS FOR EACH REPEAT
#===============================================================================
repetitions_boot <- 50

# Hyperparameters
hyperparameters <- expand.grid(
  lambda = c(0, 1, 10, 100),
  alpha = c (0.1, 0.01, 0.001)
)

set.seed(86)
seeds <- vector(mode = "list", length = repetitions_boot + 1)
for (i in 1:repetitions_boot) {
  seeds[[i]] <- sample.int(1000, nrow(hyperparameters))
}
seeds[[repetitions_boot + 1]] <- sample.int(1000, 1)

# DEFINITION OF TRAINING
#===============================================================================
control_train <- trainControl(method = "boot", number = repetitions_boot,
                              seeds = seeds, returnResamp = "final",
                              verboseIter = TRUE, allowParallel = TRUE)

# FIT MODEL 
# ==============================================================================

set.seed(86)
READ_glmnet_pvalue_50 <- caret::train(
  form = Class ~ .,
  data = x_train[c("Class", READ_filter_anova_pvalue_50)],
  method = "glmnet",
  family = "multinomial",
  tuneGrid = hyperparameters,
  metric = "Accuracy",
  trControl = control_train
)

saveRDS(object = READ_glmnet_pvalue_50, file = "READ_glmnet_pvalue_50.rds")
registerDoMC(cores = 1)

READ_glmnet_pvalue_50

# GRAPHIC REPRESENTATION
# ==============================================================================
ggplot(READ_glmnet_pvalue_50, highlight = TRUE) +
  labs(title = "Evolution of the accuracy of the GLMnet model") +
  guides(color = guide_legend(title = "Alpha"),
         shape = guide_legend(title = "Alpha")) +
  theme_bw()

# TEST PREDICTIONS
# ==============================================================================
predic_READ_glmnet_pvalue_50 <- predict(object = READ_glmnet_pvalue_50, newdata = x_test)
predic_READ_glmnet_pvalue_50
y_test <- as.factor(y_test)
caret::confusionMatrix(predic_READ_glmnet_pvalue_50, y_test)


#===============================================================================
#=========================== HDDA ==============================================
#===============================================================================


# PARALLEL PROCESS
#===============================================================================

#install.packages("doMC", repos="http://R-Forge.R-project.org")

registerDoMC(cores = 11)

# HYPERPARAMETERS, NUMBER OF REPETITIONS AND SEEDS FOR EACH REPEAT
#===============================================================================
repetitions_boot <- 50

# Hyperparameters
hyperparameters <- expand.grid(
  threshold = seq(0.1,0.9,by=0.1),
  model = "ALL")

set.seed(86)
seeds <- vector(mode = "list", length = repetitions_boot + 1)
for (i in 1:repetitions_boot) {
  seeds[[i]] <- sample.int(1000, nrow(hyperparameters))
}
seeds[[repetitions_boot + 1]] <- sample.int(1000, 1)

# DEFINITION OF TRAINING
#===============================================================================
control_train <- trainControl(method = "boot", number = repetitions_boot,
                              seeds = seeds, returnResamp = "final",
                              verboseIter = TRUE, allowParallel = TRUE)

# FIT MODEL 
# ==============================================================================

set.seed(86)
READ_hdda_pvalue_50 <- caret::train(
  form = Class ~ .,
  data = x_train[c("Class", READ_filter_anova_pvalue_50)],
  method = "hdda",
  tuneGrid = hyperparameters,
  metric = "Accuracy",
  trControl = control_train
)

saveRDS(object = READ_hdda_pvalue_50, file = "READ_hdda_pvalue_50.rds")
registerDoMC(cores = 1)

READ_hdda_pvalue_50

# GRAPHIC REPRESENTATION
# ==============================================================================
ggplot(READ_hdda_pvalue_50, highlight = TRUE) +
  labs(title = "Evolution of the accuracy of the HDDA model") +
  guides(color = guide_legend(title = "threshold"),
         shape = guide_legend(title = "threshold")) +
  theme_bw()

# TEST PREDICTIONS
# ==============================================================================
predic_READ_hdda_pvalue_50 <- predict(object = READ_hdda_pvalue_50, newdata = x_test)
predic_READ_hdda_pvalue_50
y_test <- as.factor(y_test)
caret::confusionMatrix(predic_READ_hdda_pvalue_50, y_test)


#===============================================================================
#=========================== LogitBoost ========================================
#===============================================================================


# PARALLEL PROCESS
#===============================================================================

#install.packages("doMC", repos="http://R-Forge.R-project.org")

registerDoMC(cores = 11)

# HYPERPARAMETERS, NUMBER OF REPETITIONS AND SEEDS FOR EACH REPEAT
#===============================================================================
repetitions_boot <- 50

# Hyperparameters
hyperparameters <- expand.grid(
  nIter = c(25, 50, 100, 250)
)

set.seed(86)
seeds <- vector(mode = "list", length = repetitions_boot + 1)
for (i in 1:repetitions_boot) {
  seeds[[i]] <- sample.int(1000, nrow(hyperparameters))
}
seeds[[repetitions_boot + 1]] <- sample.int(1000, 1)

# DEFINITION OF TRAINING
#===============================================================================
control_train <- trainControl(method = "boot", number = repetitions_boot,
                              seeds = seeds, returnResamp = "final",
                              verboseIter = TRUE, allowParallel = TRUE)

# FIT MODEL 
# ==============================================================================

set.seed(86)
READ_logitBoost_pvalue_50 <- caret::train(
  form = Class ~ .,
  data = x_train[c("Class", READ_filter_anova_pvalue_50)],
  method = "LogitBoost",
  tuneGrid = hyperparameters,
  metric = "Accuracy",
  trControl = control_train
)

saveRDS(object = READ_logitBoost_pvalue_50, file = "READ_logitBoost_pvalue_50.rds")
registerDoMC(cores = 1)

READ_logitBoost_pvalue_50

# GRAPHIC REPRESENTATION
# ==============================================================================
ggplot(READ_logitBoost_pvalue_50, highlight = TRUE) +
  labs(title = "Evolution of the accuracy of the Logitboost model") +
  guides(color = guide_legend(title = "nIter"),
         shape = guide_legend(title = "nIter")) +
  theme_bw()

# TEST PREDICTIONS
# ==============================================================================
predic_READ_logitBoost_pvalue_50 <- predict(object = READ_logitBoost_pvalue_50, newdata = x_test)
predic_READ_logitBoost_pvalue_50
y_test <- as.factor(y_test)
caret::confusionMatrix(predic_READ_logitBoost_pvalue_50, y_test)


#===============================================================================
#=========================== Stacking Ensemble  ================================
#===============================================================================

# DEFINICIÓN DEL ENTRENAMIENTO
#===============================================================================

stacking_df <- cbind(x_train[, READ_filter_anova_pvalue_50], x_train$Class)
names(stacking_df) [51] <- "Class"
levels(stacking_df$Class)
levels(stacking_df$Class) <- make.names(levels(stacking_df$Class))


trainControl <- trainControl(method = "repeatedcv", number = 5, repeats = 3,
                             savePredictions = "all", classProbs = TRUE,
                             index = createFolds(stacking_df$Class, 5),verboseIter = TRUE)


# Run multiplealgorithms in one call.
algorithm <- c("nnet","rf","svmRadial","xgbTree","hdda","LogitBoost")

set.seed(86)
models <- caretList(Class ~ ., data = stacking_df, 
                    trControl = trainControl,
                    methodList = algorithm)

results <- resamples(models)
summary(results)

# Box plots to compare models
scales <- list(x = list(relation = "free"), y = list(relation = "free"))
bwplot(results, scales = scales)

# Ensemble the predictions of 'models' to form a new combined prediction based on glm
set.seed(86)
stack_glm <- caretStack(models, method = "glm")
pred_stack_glm <- predict(stack_glm, newdata = x_test)
levels(pred_stack_glm) <- c("1","3")
caret::confusionMatrix(pred_stack_glm ,y_test)



#*******************************************************************************
#*************************** ANOVA P-VALUE 100 features *************************
#*******************************************************************************


#===============================================================================
#=========================== RANDOM FOREST ===================================== 
#===============================================================================


# PARALLEL PROCESS
#===============================================================================

#install.packages("doMC", repos="http://R-Forge.R-project.org")

registerDoMC(cores = 11)

# HYPERPARAMETERS, NUMBER OF REPETITIONS AND SEEDS FOR EACH REPEAT
#===============================================================================
repetitions_boot <- 50

# Hyperparameters
hyperparameters <- expand.grid(mtry = c(2, 5, 10, 50),
                               min.node.size = c(2, 3, 4, 5, 10),
                               splitrule = "gini")

set.seed(86)
seeds <- vector(mode = "list", length = repetitions_boot + 1)
for (i in 1:repetitions_boot) {
  seeds[[i]] <- sample.int(1000, nrow(hyperparameters))
}
seeds[[repetitions_boot + 1]] <- sample.int(1000, 1)

# DEFINITION OF TRAINING
#===============================================================================
control_train <- trainControl(method = "boot", number = repetitions_boot,
                              seeds = seeds, returnResamp = "final",
                              verboseIter = TRUE, allowParallel = TRUE)

# FIT MODEL 
# ==============================================================================
x_train$Class <- as.factor(x_train$Class)

set.seed(86)
READ_rf_pvalue_100 <- caret::train(
  form = Class ~ .,
  data = x_train[c("Class", READ_filter_anova_pvalue_100)],
  method = "ranger",
  tuneGrid = hyperparameters,
  metric = "Accuracy",
  trControl = control_train,
  num.trees = 500)

saveRDS(object = READ_rf_pvalue_100, file = "READ_rf_pvalue_100.rds")
registerDoMC(cores = 1)


READ_rf_pvalue_100

# GRAPHIC REPRESENTATION
# ==============================================================================
ggplot(READ_rf_pvalue_100, highlight = TRUE) +
  labs(title = "Evolution of the accuracy of the Random Forest model") +
  guides(color = guide_legend(title = "mtry"),
         shape = guide_legend(title = "mtry")) +
  theme_bw()

# TEST PREDICTIONS
# ==============================================================================
predic_READ_rf_pvalue_100 <- predict(object = READ_rf_pvalue_100, newdata = x_test)
predic_READ_rf_pvalue_100
y_test <- as.factor(y_test)
caret::confusionMatrix(predic_READ_rf_pvalue_100, y_test)


#===============================================================================
#=========================== SVM Kernel Radial =================================
#===============================================================================


# PARALLEL PROCESS
#===============================================================================

#install.packages("doMC", repos="http://R-Forge.R-project.org")

registerDoMC(cores = 11)

# HYPERPARAMETERS, NUMBER OF REPETITIONS AND SEEDS FOR EACH REPEAT
#===============================================================================
repetitions_boot <- 50

# Hyperparameters
hyperparameters <- expand.grid(sigma = c(0.0001, 0.001, 0.01),
                               C = c(1, 10, 50, 100, 250, 500, 700, 1000))

set.seed(86)
seeds <- vector(mode = "list", length = repetitions_boot + 1)
for (i in 1:repetitions_boot) {
  seeds[[i]] <- sample.int(1000, nrow(hyperparameters))
}
seeds[[repetitions_boot + 1]] <- sample.int(1000, 1)

# DEFINITION OF TRAINING
#===============================================================================
control_train <- trainControl(method = "boot", number = repetitions_boot,
                              seeds = seeds, returnResamp = "final",
                              verboseIter = TRUE, allowParallel = TRUE)

# FIT MODEL 
# ==============================================================================

set.seed(86)
READ_svmrad_pvalue_100 <- caret::train(
  form = Class ~ .,
  data = x_train[c("Class", READ_filter_anova_pvalue_100)],
  method = "svmRadial",
  tuneGrid = hyperparameters,
  metric = "Accuracy",
  trControl = control_train
)

saveRDS(object = READ_svmrad_pvalue_100, file = "READ_svmrad_pvalue_100.rds")
registerDoMC(cores = 1)


READ_svmrad_pvalue_100

# GRAPHIC REPRESENTATION
# ==============================================================================
ggplot(READ_svmrad_pvalue_100, highlight = TRUE) +
  labs(title = "Evolution of the accuracy of the SVM model") +
  guides(color = guide_legend(title = "Sigma"),
         shape = guide_legend(title = "Sigma")) +
  theme_bw()

# TEST PREDICTIONS
# ==============================================================================
predic_READ_svmrad_pvalue_100 <- predict(object = READ_svmrad_pvalue_100, newdata = x_test)
predic_READ_svmrad_pvalue_100
y_test <- as.factor(y_test)
caret::confusionMatrix(predic_READ_svmrad_pvalue_100, y_test)


#===============================================================================
#=========================== Neural Network ====================================
#===============================================================================


# PARALLEL PROCESS
#===============================================================================

#install.packages("doMC", repos="http://R-Forge.R-project.org")

registerDoMC(cores = 11)

# HYPERPARAMETERS, NUMBER OF REPETITIONS AND SEEDS FOR EACH REPEAT
#===============================================================================
repetitions_boot <- 50

# Hyperparameters
hyperparameters <- expand.grid(size = c(4, 8, 10, 16),
                               decay = c(0.01, 0.1))

set.seed(86)
seeds <- vector(mode = "list", length = repetitions_boot + 1)
for (i in 1:repetitions_boot) {
  seeds[[i]] <- sample.int(1000, nrow(hyperparameters))
}
seeds[[repetitions_boot + 1]] <- sample.int(1000, 1)

# DEFINITION OF TRAINING
#===============================================================================
control_train <- trainControl(method = "boot", number = repetitions_boot,
                              seeds = seeds, returnResamp = "final",
                              verboseIter = TRUE, allowParallel = TRUE)

# FIT MODEL 
# ==============================================================================

set.seed(86)
READ_nnet_pvalue_100 <- caret::train(
  form = Class ~ .,
  data = x_train[c("Class", READ_filter_anova_pvalue_100)],
  method = "nnet",
  tuneGrid = hyperparameters,
  metric = "Accuracy",
  trControl = control_train,
  rang = c(-0.7, 0.7),
  MaxNWts = 21000,
  trace = FALSE
)

saveRDS(object = READ_nnet_pvalue_100, file = "READ_nnet_pvalue_100.rds")
registerDoMC(cores = 1)

READ_nnet_pvalue_100

# GRAPHIC REPRESENTATION
# ==============================================================================
ggplot(READ_nnet_pvalue_100, highlight = TRUE) +
  labs(title = "Evolution of the accuracy of the Neural Net model") +
  guides(color = guide_legend(title = "Decay"),
         shape = guide_legend(title = "Decay")) +
  theme_bw()

# TEST PREDICTIONS
# ==============================================================================
predic_READ_nnet_pvalue_100 <- predict(object = READ_nnet_pvalue_100, newdata = x_test)
predic_READ_nnet_pvalue_100
y_test <- as.factor(y_test)
caret::confusionMatrix(predic_READ_nnet_pvalue_100, y_test)



#===============================================================================
#=========================== Gradient boosting =================================
#===============================================================================


# PARALLEL PROCESS
#===============================================================================

#install.packages("doMC", repos="http://R-Forge.R-project.org")

registerDoMC(cores = 11)

# HYPERPARAMETERS, NUMBER OF REPETITIONS AND SEEDS FOR EACH REPEAT
#===============================================================================
repetitions_boot <- 50

# Hyperparameters
hyperparameters <- expand.grid(
  shrinkage = c(0.01, 0.1, 0.3),
  interaction.depth = c(1, 3, 5),
  n.minobsinnode = c(5, 10, 15),
  n.trees = c(500, 1000)
  # bag.fraction = c(0.65, 0.8, 1), 
  # optimal_trees = 0,               
  # min_RMSE = 0,                    
  # min_cor = 0
)

set.seed(86)
seeds <- vector(mode = "list", length = repetitions_boot + 1)
for (i in 1:repetitions_boot) {
  seeds[[i]] <- sample.int(1000, nrow(hyperparameters))
}
seeds[[repetitions_boot + 1]] <- sample.int(1000, 1)

# DEFINITION OF TRAINING
#===============================================================================
control_train <- trainControl(method = "boot", number = repetitions_boot,
                              seeds = seeds, returnResamp = "final",
                              verboseIter = TRUE, allowParallel = TRUE)

# FIT MODEL 
# ==============================================================================

set.seed(86)
READ_gbm_pvalue_100 <- caret::train(
  form = Class ~ .,
  data = x_train[c("Class", READ_filter_anova_pvalue_100)],
  method = "gbm",
  tuneGrid = hyperparameters,
  metric = "Accuracy",
  preProcess = c("center", "scale"),
  trControl = control_train
)

saveRDS(object = READ_gbm_pvalue_100, file = "READ_gbm_pvalue_100.rds")
registerDoMC(cores = 1)

READ_gbm_pvalue_100

# GRAPHIC REPRESENTATION
# ==============================================================================
ggplot(READ_gbm_pvalue_100, highlight = TRUE) +
  labs(title = "Evolution of the accuracy of the READ model") +
  guides(color = guide_legend(title = "Size"),
         shape = guide_legend(title = "Size")) +
  theme_bw()

# TEST PREDICTIONS
# ==============================================================================
predic_READ_gbm_pvalue_100 <- predict(object = READ_gbm_pvalue_100, newdata = x_test)
predic_READ_gbm_pvalue_100
y_test <- as.factor(y_test)
caret::confusionMatrix(predic_READ_gbm_pvalue_100, y_test)


#===============================================================================
#=========================== XGBM ==============================================
#===============================================================================


# PARALLEL PROCESS
#===============================================================================

#install.packages("doMC", repos="http://R-Forge.R-project.org")

registerDoMC(cores = 11)

# HYPERPARAMETERS, NUMBER OF REPETITIONS AND SEEDS FOR EACH REPEAT
#===============================================================================
repetitions_boot <- 50

# Hyperparameters
hyperparameters <- expand.grid(
  nrounds = c(500, 1000),
  eta = c(0.01, 0.001), 
  max_depth = c(2, 4, 6),
  gamma = 1,
  colsample_bytree = c(0.2, 0.4),
  min_child_weight = c(1, 5),
  subsample = 1
)


set.seed(86)
seeds <- vector(mode = "list", length = repetitions_boot + 1)
for (i in 1:repetitions_boot) {
  seeds[[i]] <- sample.int(1000, nrow(hyperparameters))
}
seeds[[repetitions_boot + 1]] <- sample.int(1000, 1)

# DEFINITION OF TRAINING
#===============================================================================
control_train <- trainControl(method = "boot", number = repetitions_boot,
                              seeds = seeds, returnResamp = "final",
                              verboseIter = TRUE, allowParallel = TRUE)

# FIT MODEL 
# ==============================================================================

set.seed(86)
READ_xgbm_pvalue_100 <- caret::train(
  form = Class ~ .,
  data = x_train[c("Class", READ_filter_anova_pvalue_100)],
  method = "xgbTree",
  tuneGrid = hyperparameters,
  metric = "Accuracy",
  preProcess = c("center", "scale"),
  trControl = control_train
)

saveRDS(object = READ_xgbm_pvalue_100, file = "READ_xgbm_pvalue_100.rds")
registerDoMC(cores = 1)

READ_xgbm_pvalue_100

# GRAPHIC REPRESENTATION
# ==============================================================================
plot(READ_xgbm_pvalue_100, highlight = TRUE)

# TEST PREDICTIONS
# ==============================================================================
predic_READ_xgbm_pvalue_100 <- predict(object = READ_xgbm_pvalue_100, newdata = x_test)
predic_READ_xgbm_pvalue_100
y_test <- as.factor(y_test)
caret::confusionMatrix(predic_READ_xgbm_pvalue_100, y_test)


#===============================================================================
#=========================== glmnet ============================================
#===============================================================================


# PARALLEL PROCESS
#===============================================================================

#install.packages("doMC", repos="http://R-Forge.R-project.org")

registerDoMC(cores = 11)

# HYPERPARAMETERS, NUMBER OF REPETITIONS AND SEEDS FOR EACH REPEAT
#===============================================================================
repetitions_boot <- 50

# Hyperparameters
hyperparameters <- expand.grid(
  lambda = c(0, 1, 10, 100),
  alpha = c (0.1, 0.01, 0.001)
)

set.seed(86)
seeds <- vector(mode = "list", length = repetitions_boot + 1)
for (i in 1:repetitions_boot) {
  seeds[[i]] <- sample.int(1000, nrow(hyperparameters))
}
seeds[[repetitions_boot + 1]] <- sample.int(1000, 1)

# DEFINITION OF TRAINING
#===============================================================================
control_train <- trainControl(method = "boot", number = repetitions_boot,
                              seeds = seeds, returnResamp = "final",
                              verboseIter = TRUE, allowParallel = TRUE)

# FIT MODEL 
# ==============================================================================

set.seed(86)
READ_glmnet_pvalue_100 <- caret::train(
  form = Class ~ .,
  data = x_train[c("Class", READ_filter_anova_pvalue_100)],
  method = "glmnet",
  family = "multinomial",
  tuneGrid = hyperparameters,
  metric = "Accuracy",
  trControl = control_train
)

saveRDS(object = READ_glmnet_pvalue_100, file = "READ_glmnet_pvalue_100.rds")
registerDoMC(cores = 1)

READ_glmnet_pvalue_100

# GRAPHIC REPRESENTATION
# ==============================================================================
ggplot(READ_glmnet_pvalue_100, highlight = TRUE) +
  labs(title = "Evolution of the accuracy of the GLMnet model") +
  guides(color = guide_legend(title = "Alpha"),
         shape = guide_legend(title = "Alpha")) +
  theme_bw()

# TEST PREDICTIONS
# ==============================================================================
predic_READ_glmnet_pvalue_100 <- predict(object = READ_glmnet_pvalue_100, newdata = x_test)
predic_READ_glmnet_pvalue_100
y_test <- as.factor(y_test)
caret::confusionMatrix(predic_READ_glmnet_pvalue_100, y_test)


#===============================================================================
#=========================== HDDA ==============================================
#===============================================================================


# PARALLEL PROCESS
#===============================================================================

#install.packages("doMC", repos="http://R-Forge.R-project.org")

registerDoMC(cores = 11)

# HYPERPARAMETERS, NUMBER OF REPETITIONS AND SEEDS FOR EACH REPEAT
#===============================================================================
repetitions_boot <- 50

# Hyperparameters
hyperparameters <- expand.grid(
  threshold = seq(0.1,0.9,by=0.1),
  model = "ALL")

set.seed(86)
seeds <- vector(mode = "list", length = repetitions_boot + 1)
for (i in 1:repetitions_boot) {
  seeds[[i]] <- sample.int(1000, nrow(hyperparameters))
}
seeds[[repetitions_boot + 1]] <- sample.int(1000, 1)

# DEFINITION OF TRAINING
#===============================================================================
control_train <- trainControl(method = "boot", number = repetitions_boot,
                              seeds = seeds, returnResamp = "final",
                              verboseIter = TRUE, allowParallel = TRUE)

# FIT MODEL 
# ==============================================================================

set.seed(86)
READ_hdda_pvalue_100 <- caret::train(
  form = Class ~ .,
  data = x_train[c("Class", READ_filter_anova_pvalue_100)],
  method = "hdda",
  tuneGrid = hyperparameters,
  metric = "Accuracy",
  trControl = control_train
)

saveRDS(object = READ_hdda_pvalue_100, file = "READ_hdda_pvalue_100.rds")
registerDoMC(cores = 1)

READ_hdda_pvalue_100

# GRAPHIC REPRESENTATION
# ==============================================================================
ggplot(READ_hdda_pvalue_100, highlight = TRUE) +
  labs(title = "Evolution of the accuracy of the HDDA model") +
  guides(color = guide_legend(title = "threshold"),
         shape = guide_legend(title = "threshold")) +
  theme_bw()

# TEST PREDICTIONS
# ==============================================================================
predic_READ_hdda_pvalue_100 <- predict(object = READ_hdda_pvalue_100, newdata = x_test)
predic_READ_hdda_pvalue_100
y_test <- as.factor(y_test)
caret::confusionMatrix(predic_READ_hdda_pvalue_100, y_test)


#===============================================================================
#=========================== LogitBoost ========================================
#===============================================================================


# PARALLEL PROCESS
#===============================================================================

#install.packages("doMC", repos="http://R-Forge.R-project.org")

registerDoMC(cores = 11)

# HYPERPARAMETERS, NUMBER OF REPETITIONS AND SEEDS FOR EACH REPEAT
#===============================================================================
repetitions_boot <- 50

# Hyperparameters
hyperparameters <- expand.grid(
  nIter = c(25, 50, 100, 250)
)

set.seed(86)
seeds <- vector(mode = "list", length = repetitions_boot + 1)
for (i in 1:repetitions_boot) {
  seeds[[i]] <- sample.int(1000, nrow(hyperparameters))
}
seeds[[repetitions_boot + 1]] <- sample.int(1000, 1)

# DEFINITION OF TRAINING
#===============================================================================
control_train <- trainControl(method = "boot", number = repetitions_boot,
                              seeds = seeds, returnResamp = "final",
                              verboseIter = TRUE, allowParallel = TRUE)

# FIT MODEL 
# ==============================================================================

set.seed(86)
READ_logitBoost_pvalue_100 <- caret::train(
  form = Class ~ .,
  data = x_train[c("Class", READ_filter_anova_pvalue_100)],
  method = "LogitBoost",
  tuneGrid = hyperparameters,
  metric = "Accuracy",
  trControl = control_train
)

saveRDS(object = READ_logitBoost_pvalue_100, file = "READ_logitBoost_pvalue_100.rds")
registerDoMC(cores = 1)

READ_logitBoost_pvalue_100

# GRAPHIC REPRESENTATION
# ==============================================================================
ggplot(READ_logitBoost_pvalue_100, highlight = TRUE) +
  labs(title = "Evolution of the accuracy of the Logitboost model") +
  guides(color = guide_legend(title = "nIter"),
         shape = guide_legend(title = "nIter")) +
  theme_bw()

# TEST PREDICTIONS
# ==============================================================================
predic_READ_logitBoost_pvalue_100 <- predict(object = READ_logitBoost_pvalue_100, newdata = x_test)
predic_READ_logitBoost_pvalue_100
y_test <- as.factor(y_test)
caret::confusionMatrix(predic_READ_logitBoost_pvalue_100, y_test)


#===============================================================================
#=========================== Stacking Ensemble  ================================
#===============================================================================

# DEFINICIÓN DEL ENTRENAMIENTO
#===============================================================================

stacking_df <- cbind(x_train[,READ_filter_anova_pvalue_100],x_train$Class)
names(stacking_df) [101] <- "Class"
levels(stacking_df$Class)
levels(stacking_df$Class) <- make.names(levels(stacking_df$Class))


trainControl <- trainControl(method = "repeatedcv", number = 5, repeats = 3,
                             savePredictions = "all", classProbs = TRUE,
                             index = createFolds(stacking_df$Class, 5),verboseIter = TRUE)


# Run multiplealgorithms in one call.
algorithm <- c("nnet","rf","svmRadial","READ","xgbTree","hdda","LogitBoost")

set.seed(86)
models <- caretList(Class ~ ., data = stacking_df, 
                    trControl = trainControl,
                    methodList = algorithm)

results <- resamples(models)
summary(results)

# Box plots to compare models
scales <- list(x = list(relation = "free"), y = list(relation = "free"))
bwplot(results, scales = scales)

# Ensemble the predictions of 'models' to form a new combined prediction based on glm
set.seed(86)
stack_glm <- caretStack(models, method = "glm")
pred_stack_glm <- predict(stack_glm, newdata = x_test)
levels(pred_stack_glm) <- c("1","2", "3")
caret::confusionMatrix(pred_stack_glm ,y_test)



#*******************************************************************************
#*************************** ANOVA P-VALUE 500 features *************************
#*******************************************************************************


#===============================================================================
#=========================== RANDOM FOREST ===================================== 
#===============================================================================


# PARALLEL PROCESS
#===============================================================================

#install.packages("doMC", repos="http://R-Forge.R-project.org")

registerDoMC(cores = 11)

# HYPERPARAMETERS, NUMBER OF REPETITIONS AND SEEDS FOR EACH REPEAT
#===============================================================================
repetitions_boot <- 50

# Hyperparameters
hyperparameters <- expand.grid(mtry = c(2, 5, 10, 50),
                               min.node.size = c(2, 3, 4, 5, 10),
                               splitrule = "gini")

set.seed(86)
seeds <- vector(mode = "list", length = repetitions_boot + 1)
for (i in 1:repetitions_boot) {
  seeds[[i]] <- sample.int(1000, nrow(hyperparameters))
}
seeds[[repetitions_boot + 1]] <- sample.int(1000, 1)

# DEFINITION OF TRAINING
#===============================================================================
control_train <- trainControl(method = "boot", number = repetitions_boot,
                              seeds = seeds, returnResamp = "final",
                              verboseIter = TRUE, allowParallel = TRUE)

# FIT MODEL 
# ==============================================================================
x_train$Class <- as.factor(x_train$Class)

set.seed(86)
READ_rf_pvalue_500 <- caret::train(
  form = Class ~ .,
  data = x_train[c("Class", READ_filter_anova_pvalue_500)],
  method = "ranger",
  tuneGrid = hyperparameters,
  metric = "Accuracy",
  trControl = control_train,
  num.trees = 500)

saveRDS(object = READ_rf_pvalue_500, file = "READ_rf_pvalue_500.rds")
registerDoMC(cores = 1)


READ_rf_pvalue_500

# GRAPHIC REPRESENTATION
# ==============================================================================
ggplot(READ_rf_pvalue_500, highlight = TRUE) +
  labs(title = "Evolution of the accuracy of the Random Forest model") +
  guides(color = guide_legend(title = "mtry"),
         shape = guide_legend(title = "mtry")) +
  theme_bw()

# TEST PREDICTIONS
# ==============================================================================
predic_READ_rf_pvalue_500 <- predict(object = READ_rf_pvalue_500, newdata = x_test)
predic_READ_rf_pvalue_500
y_test <- as.factor(y_test)
caret::confusionMatrix(predic_READ_rf_pvalue_500, y_test)


#===============================================================================
#=========================== SVM Kernel Radial =================================
#===============================================================================


# PARALLEL PROCESS
#===============================================================================

#install.packages("doMC", repos="http://R-Forge.R-project.org")

registerDoMC(cores = 11)

# HYPERPARAMETERS, NUMBER OF REPETITIONS AND SEEDS FOR EACH REPEAT
#===============================================================================
repetitions_boot <- 50

# Hyperparameters
hyperparameters <- expand.grid(sigma = c(0.0001, 0.001, 0.01),
                               C = c(1, 10, 50, 100, 250, 500, 700, 1000))

set.seed(86)
seeds <- vector(mode = "list", length = repetitions_boot + 1)
for (i in 1:repetitions_boot) {
  seeds[[i]] <- sample.int(1000, nrow(hyperparameters))
}
seeds[[repetitions_boot + 1]] <- sample.int(1000, 1)

# DEFINITION OF TRAINING
#===============================================================================
control_train <- trainControl(method = "boot", number = repetitions_boot,
                              seeds = seeds, returnResamp = "final",
                              verboseIter = TRUE, allowParallel = TRUE)

# FIT MODEL 
# ==============================================================================

set.seed(86)
READ_svmrad_pvalue_500 <- caret::train(
  form = Class ~ .,
  data = x_train[c("Class", READ_filter_anova_pvalue_500)],
  method = "svmRadial",
  tuneGrid = hyperparameters,
  metric = "Accuracy",
  trControl = control_train
)

saveRDS(object = READ_svmrad_pvalue_500, file = "READ_svmrad_pvalue_500.rds")
registerDoMC(cores = 1)


READ_svmrad_pvalue_500

# GRAPHIC REPRESENTATION
# ==============================================================================
ggplot(READ_svmrad_pvalue_500, highlight = TRUE) +
  labs(title = "Evolution of the accuracy of the SVM model") +
  guides(color = guide_legend(title = "Sigma"),
         shape = guide_legend(title = "Sigma")) +
  theme_bw()

# TEST PREDICTIONS
# ==============================================================================
predic_READ_svmrad_pvalue_500 <- predict(object = READ_svmrad_pvalue_500, newdata = x_test)
predic_READ_svmrad_pvalue_500
y_test <- as.factor(y_test)
caret::confusionMatrix(predic_READ_svmrad_pvalue_500, y_test)


#===============================================================================
#=========================== Neural Network ====================================
#===============================================================================


# PARALLEL PROCESS
#===============================================================================

#install.packages("doMC", repos="http://R-Forge.R-project.org")

registerDoMC(cores = 11)

# HYPERPARAMETERS, NUMBER OF REPETITIONS AND SEEDS FOR EACH REPEAT
#===============================================================================
repetitions_boot <- 50

# Hyperparameters
hyperparameters <- expand.grid(size = c(2, 4, 6, 8, 10),
                               decay = c(0.01, 0.1))

set.seed(86)
seeds <- vector(mode = "list", length = repetitions_boot + 1)
for (i in 1:repetitions_boot) {
  seeds[[i]] <- sample.int(1000, nrow(hyperparameters))
}
seeds[[repetitions_boot + 1]] <- sample.int(1000, 1)

# DEFINITION OF TRAINING
#===============================================================================
control_train <- trainControl(method = "boot", number = repetitions_boot,
                              seeds = seeds, returnResamp = "final",
                              verboseIter = TRUE, allowParallel = TRUE)

# FIT MODEL 
# ==============================================================================

set.seed(86)
READ_nnet_pvalue_500 <- caret::train(
  form = Class ~ .,
  data = x_train[c("Class", READ_filter_anova_pvalue_500)],
  method = "nnet",
  tuneGrid = hyperparameters,
  metric = "Accuracy",
  trControl = control_train,
  rang = c(-0.7, 0.7),
  MaxNWts = 21000,
  trace = FALSE
)

saveRDS(object = READ_nnet_pvalue_500, file = "READ_nnet_pvalue_500.rds")
registerDoMC(cores = 1)

READ_nnet_pvalue_500

# GRAPHIC REPRESENTATION
# ==============================================================================
ggplot(READ_nnet_pvalue_500, highlight = TRUE) +
  labs(title = "Evolution of the accuracy of the Neural Net model") +
  guides(color = guide_legend(title = "Decay"),
         shape = guide_legend(title = "Decay")) +
  theme_bw()

# TEST PREDICTIONS
# ==============================================================================
predic_READ_nnet_pvalue_500 <- predict(object = READ_nnet_pvalue_500, newdata = x_test)
predic_READ_nnet_pvalue_500
y_test <- as.factor(y_test)
caret::confusionMatrix(predic_READ_nnet_pvalue_500, y_test)



#===============================================================================
#=========================== Gradient boosting =================================
#===============================================================================


# PARALLEL PROCESS
#===============================================================================

#install.packages("doMC", repos="http://R-Forge.R-project.org")

registerDoMC(cores = 11)

# HYPERPARAMETERS, NUMBER OF REPETITIONS AND SEEDS FOR EACH REPEAT
#===============================================================================
repetitions_boot <- 50

# Hyperparameters
hyperparameters <- expand.grid(
  shrinkage = c(0.01, 0.1, 0.3),
  interaction.depth = c(1, 3, 5),
  n.minobsinnode = c(5, 10, 15),
  n.trees = c(500, 1000)
  # bag.fraction = c(0.65, 0.8, 1), 
  # optimal_trees = 0,               
  # min_RMSE = 0,                    
  # min_cor = 0
)

set.seed(86)
seeds <- vector(mode = "list", length = repetitions_boot + 1)
for (i in 1:repetitions_boot) {
  seeds[[i]] <- sample.int(1000, nrow(hyperparameters))
}
seeds[[repetitions_boot + 1]] <- sample.int(1000, 1)

# DEFINITION OF TRAINING
#===============================================================================
control_train <- trainControl(method = "boot", number = repetitions_boot,
                              seeds = seeds, returnResamp = "final",
                              verboseIter = TRUE, allowParallel = TRUE)

# FIT MODEL 
# ==============================================================================

set.seed(86)
READ_gbm_pvalue_500 <- caret::train(
  form = Class ~ .,
  data = x_train[c("Class", READ_filter_anova_pvalue_500)],
  method = "gbm",
  tuneGrid = hyperparameters,
  metric = "Accuracy",
  preProcess = c("center", "scale"),
  trControl = control_train
)

saveRDS(object = READ_gbm_pvalue_500, file = "READ_gbm_pvalue_500.rds")
registerDoMC(cores = 1)

READ_gbm_pvalue_500

# GRAPHIC REPRESENTATION
# ==============================================================================
ggplot(READ_gbm_pvalue_500, highlight = TRUE) +
  labs(title = "Evolution of the accuracy of the READ model") +
  guides(color = guide_legend(title = "Size"),
         shape = guide_legend(title = "Size")) +
  theme_bw()

# TEST PREDICTIONS
# ==============================================================================
predic_READ_gbm_pvalue_500 <- predict(object = READ_gbm_pvalue_500, newdata = x_test)
predic_READ_gbm_pvalue_500
y_test <- as.factor(y_test)
caret::confusionMatrix(predic_READ_gbm_pvalue_500, y_test)


#===============================================================================
#=========================== XGBM ==============================================
#===============================================================================


# PARALLEL PROCESS
#===============================================================================

#install.packages("doMC", repos="http://R-Forge.R-project.org")

registerDoMC(cores = 11)

# HYPERPARAMETERS, NUMBER OF REPETITIONS AND SEEDS FOR EACH REPEAT
#===============================================================================
repetitions_boot <- 50

# Hyperparameters
hyperparameters <- expand.grid(
  nrounds = c(500, 1000),
  eta = c(0.01, 0.001), 
  max_depth = c(2, 4, 6),
  gamma = 1,
  colsample_bytree = c(0.2, 0.4),
  min_child_weight = c(1, 5),
  subsample = 1
)


set.seed(86)
seeds <- vector(mode = "list", length = repetitions_boot + 1)
for (i in 1:repetitions_boot) {
  seeds[[i]] <- sample.int(1000, nrow(hyperparameters))
}
seeds[[repetitions_boot + 1]] <- sample.int(1000, 1)

# DEFINITION OF TRAINING
#===============================================================================
control_train <- trainControl(method = "boot", number = repetitions_boot,
                              seeds = seeds, returnResamp = "final",
                              verboseIter = TRUE, allowParallel = TRUE)

# FIT MODEL 
# ==============================================================================

set.seed(86)
READ_xgbm_pvalue_500 <- caret::train(
  form = Class ~ .,
  data = x_train[c("Class", READ_filter_anova_pvalue_500)],
  method = "xgbTree",
  tuneGrid = hyperparameters,
  metric = "Accuracy",
  preProcess = c("center", "scale"),
  trControl = control_train
)

saveRDS(object = READ_xgbm_pvalue_500, file = "READ_xgbm_pvalue_500.rds")
registerDoMC(cores = 1)

READ_xgbm_pvalue_500

# GRAPHIC REPRESENTATION
# ==============================================================================
plot(READ_xgbm_pvalue_500, highlight = TRUE)

# TEST PREDICTIONS
# ==============================================================================
predic_READ_xgbm_pvalue_500 <- predict(object = READ_xgbm_pvalue_500, newdata = x_test)
predic_READ_xgbm_pvalue_500
y_test <- as.factor(y_test)
caret::confusionMatrix(predic_READ_xgbm_pvalue_500, y_test)


#===============================================================================
#=========================== glmnet ============================================
#===============================================================================


# PARALLEL PROCESS
#===============================================================================

#install.packages("doMC", repos="http://R-Forge.R-project.org")

registerDoMC(cores = 11)

# HYPERPARAMETERS, NUMBER OF REPETITIONS AND SEEDS FOR EACH REPEAT
#===============================================================================
repetitions_boot <- 50

# Hyperparameters
hyperparameters <- expand.grid(
  lambda = c(0, 1, 10, 100),
  alpha = c (0.1, 0.01, 0.001)
)

set.seed(86)
seeds <- vector(mode = "list", length = repetitions_boot + 1)
for (i in 1:repetitions_boot) {
  seeds[[i]] <- sample.int(1000, nrow(hyperparameters))
}
seeds[[repetitions_boot + 1]] <- sample.int(1000, 1)

# DEFINITION OF TRAINING
#===============================================================================
control_train <- trainControl(method = "boot", number = repetitions_boot,
                              seeds = seeds, returnResamp = "final",
                              verboseIter = TRUE, allowParallel = TRUE)

# FIT MODEL 
# ==============================================================================

set.seed(86)
READ_glmnet_pvalue_500 <- caret::train(
  form = Class ~ .,
  data = x_train[c("Class", READ_filter_anova_pvalue_500)],
  method = "glmnet",
  family = "multinomial",
  tuneGrid = hyperparameters,
  metric = "Accuracy",
  trControl = control_train
)

saveRDS(object = READ_glmnet_pvalue_500, file = "READ_glmnet_pvalue_500.rds")
registerDoMC(cores = 1)

READ_glmnet_pvalue_500

# GRAPHIC REPRESENTATION
# ==============================================================================
ggplot(READ_glmnet_pvalue_500, highlight = TRUE) +
  labs(title = "Evolution of the accuracy of the GLMnet model") +
  guides(color = guide_legend(title = "Alpha"),
         shape = guide_legend(title = "Alpha")) +
  theme_bw()

# TEST PREDICTIONS
# ==============================================================================
predic_READ_glmnet_pvalue_500 <- predict(object = READ_glmnet_pvalue_500, newdata = x_test)
predic_READ_glmnet_pvalue_500
y_test <- as.factor(y_test)
caret::confusionMatrix(predic_READ_glmnet_pvalue_500, y_test)


#===============================================================================
#=========================== HDDA ==============================================
#===============================================================================


# PARALLEL PROCESS
#===============================================================================

#install.packages("doMC", repos="http://R-Forge.R-project.org")

registerDoMC(cores = 11)

# HYPERPARAMETERS, NUMBER OF REPETITIONS AND SEEDS FOR EACH REPEAT
#===============================================================================
repetitions_boot <- 50

# Hyperparameters
hyperparameters <- expand.grid(
  threshold = seq(0.1,0.9,by=0.1),
  model = "ALL")

set.seed(86)
seeds <- vector(mode = "list", length = repetitions_boot + 1)
for (i in 1:repetitions_boot) {
  seeds[[i]] <- sample.int(1000, nrow(hyperparameters))
}
seeds[[repetitions_boot + 1]] <- sample.int(1000, 1)

# DEFINITION OF TRAINING
#===============================================================================
control_train <- trainControl(method = "boot", number = repetitions_boot,
                              seeds = seeds, returnResamp = "final",
                              verboseIter = TRUE, allowParallel = TRUE)

# FIT MODEL 
# ==============================================================================

set.seed(86)
READ_hdda_pvalue_500 <- caret::train(
  form = Class ~ .,
  data = x_train[c("Class", READ_filter_anova_pvalue_500)],
  method = "hdda",
  tuneGrid = hyperparameters,
  metric = "Accuracy",
  trControl = control_train
)

saveRDS(object = READ_hdda_pvalue_500, file = "READ_hdda_pvalue_500.rds")
registerDoMC(cores = 1)

READ_hdda_pvalue_500

# GRAPHIC REPRESENTATION
# ==============================================================================
ggplot(READ_hdda_pvalue_500, highlight = TRUE) +
  labs(title = "Evolution of the accuracy of the HDDA model") +
  guides(color = guide_legend(title = "threshold"),
         shape = guide_legend(title = "threshold")) +
  theme_bw()

# TEST PREDICTIONS
# ==============================================================================
predic_READ_hdda_pvalue_500 <- predict(object = READ_hdda_pvalue_500, newdata = x_test)
predic_READ_hdda_pvalue_500
y_test <- as.factor(y_test)
caret::confusionMatrix(predic_READ_hdda_pvalue_500, y_test)


#===============================================================================
#=========================== LogitBoost ========================================
#===============================================================================


# PARALLEL PROCESS
#===============================================================================

#install.packages("doMC", repos="http://R-Forge.R-project.org")

registerDoMC(cores = 11)

# HYPERPARAMETERS, NUMBER OF REPETITIONS AND SEEDS FOR EACH REPEAT
#===============================================================================
repetitions_boot <- 50

# Hyperparameters
hyperparameters <- expand.grid(
  nIter = c(25, 50, 100, 250)
)

set.seed(86)
seeds <- vector(mode = "list", length = repetitions_boot + 1)
for (i in 1:repetitions_boot) {
  seeds[[i]] <- sample.int(1000, nrow(hyperparameters))
}
seeds[[repetitions_boot + 1]] <- sample.int(1000, 1)

# DEFINITION OF TRAINING
#===============================================================================
control_train <- trainControl(method = "boot", number = repetitions_boot,
                              seeds = seeds, returnResamp = "final",
                              verboseIter = TRUE, allowParallel = TRUE)

# FIT MODEL 
# ==============================================================================

set.seed(86)
READ_logitBoost_pvalue_500 <- caret::train(
  form = Class ~ .,
  data = x_train[c("Class", READ_filter_anova_pvalue_500)],
  method = "LogitBoost",
  tuneGrid = hyperparameters,
  metric = "Accuracy",
  trControl = control_train
)

saveRDS(object = READ_logitBoost_pvalue_500, file = "READ_logitBoost_pvalue_500.rds")
registerDoMC(cores = 1)

READ_logitBoost_pvalue_500

# GRAPHIC REPRESENTATION
# ==============================================================================
ggplot(READ_logitBoost_pvalue_500, highlight = TRUE) +
  labs(title = "Evolution of the accuracy of the Logitboost model") +
  guides(color = guide_legend(title = "nIter"),
         shape = guide_legend(title = "nIter")) +
  theme_bw()

# TEST PREDICTIONS
# ==============================================================================
predic_READ_logitBoost_pvalue_500 <- predict(object = READ_logitBoost_pvalue_500, newdata = x_test)
predic_READ_logitBoost_pvalue_500
y_test <- as.factor(y_test)
caret::confusionMatrix(predic_READ_logitBoost_pvalue_500, y_test)


#===============================================================================
#=========================== Stacking Ensemble  ================================
#===============================================================================

# DEFINICIÓN DEL ENTRENAMIENTO
#===============================================================================

stacking_df <- cbind(x_train[,READ_filter_anova_pvalue_500],x_train$Class)
names(stacking_df) [501] <- "Class"
levels(stacking_df$Class)
levels(stacking_df$Class) <- make.names(levels(stacking_df$Class))


trainControl <- trainControl(method = "repeatedcv", number = 5, repeats = 3,
                             savePredictions = "all", classProbs = TRUE,
                             index = createFolds(stacking_df$Class, 5),verboseIter = TRUE)


# Run multiplealgorithms in one call.
algorithm <- c("nnet","rf","svmRadial","READ","xgbTree","hdda","LogitBoost")

set.seed(86)
models <- caretList(Class ~ ., data = stacking_df, 
                    trControl = trainControl,
                    methodList = algorithm)

results <- resamples(models)
summary(results)

# Box plots to compare models
scales <- list(x = list(relation = "free"), y = list(relation = "free"))
bwplot(results, scales = scales)

# Ensemble the predictions of 'models' to form a new combined prediction based on glm
set.seed(86)
stack_glm <- caretStack(models, method = "glm")
pred_stack_glm <- predict(stack_glm, newdata = x_test)
levels(pred_stack_glm) <- c("1","2")
caret::confusionMatrix(pred_stack_glm ,y_test)



#*******************************************************************************
#************************ Recursive Feature Elimination ************************
#*******************************************************************************


# sets a limit on the number of nested expressions 
options(expressions = 500000)


zero_var <- nearZeroVar(x_train)
class_index <- which(colnames(x_train) == "Class")
p <- ncol(x_train)/2
x <- findCorrelation(cor(x_train[, - class_index]), .8, verbose = TRUE)
x_train_corr <- x_train[, x]
x_train_corr <- cbind(x_train_corr, x_train$Class) 
colnames(x_train_corr) [length(x_train_corr)] <- "Class"

# define the control using a random forest selection function
set.seed(86)
rfe_control <- rfeControl(functions = rfFuncs, 
                          method = "cv", 
                          number = 10,
                          repeats = 3,
                          allowParallel = TRUE)

# run the RFE algorithm
RFE_results <- rfe(form = Class ~ . ,
                   data = x_train_corr, 
                   sizes = c(1:124), 
                   rfeControl = rfe_control)

# summarize the results
print(RFE_results)

# list the chosen features
predictors(RFE_results)

# plot the results
plot(RFE_results, type=c("g", "o"))

selected_vars <- RFE_results$variables
write.csv(selected_vars, "READ_selected_vars_RFE.csv")
saveRDS(RFE_results, "READ_RFE_results.rds")
best_vars_50 <- RFE_results$control$functions$selectVar(selected_vars, 50)

x_train_rfe_50 <- x_train[, c("Class", best_vars_50)] 
x_test_rfe_50 <- x_test[, c("Class", best_vars_50)] 


#===============================================================================
#=========================== RANDOM FOREST ===================================== 
#===============================================================================


# PARALLEL PROCESS
#===============================================================================

#install.packages("doMC", repos="http://R-Forge.R-project.org")

registerDoMC(cores = 11)

# HYPERPARAMETERS, NUMBER OF REPETITIONS AND SEEDS FOR EACH REPEAT
#===============================================================================
repetitions_boot <- 50

# Hyperparameters
hyperparameters <- expand.grid(mtry = c(2, 5, 10, 50),
                               min.node.size = c(2, 3, 4, 5, 10),
                               splitrule = "gini")

set.seed(86)
seeds <- vector(mode = "list", length = repetitions_boot + 1)
for (i in 1:repetitions_boot) {
  seeds[[i]] <- sample.int(1000, nrow(hyperparameters))
}
seeds[[repetitions_boot + 1]] <- sample.int(1000, 1)

# DEFINITION OF TRAINING
#===============================================================================
control_train <- trainControl(method = "boot", number = repetitions_boot,
                              seeds = seeds, returnResamp = "final",
                              verboseIter = TRUE, allowParallel = TRUE)

# FIT MODEL 
# ==============================================================================
x_train_rfe_50$Class <- as.factor(x_train_rfe_50$Class)

set.seed(86)
READ_rf_RFE_50 <- caret::train(
  form = Class ~ .,
  data = x_train_rfe_50,
  method = "ranger",
  tuneGrid = hyperparameters,
  metric = "Accuracy",
  trControl = control_train,
  num.trees = 500)

saveRDS(object = READ_rf_RFE_50, file = "READ_rf_RFE_50.rds")
registerDoMC(cores = 1)


READ_rf_RFE_50

# GRAPHIC REPRESENTATION
# ==============================================================================
ggplot(READ_rf_RFE_50, highlight = TRUE) +
  labs(title = "Evolution of the accuracy of the Random Forest model") +
  guides(color = guide_legend(title = "mtry"),
         shape = guide_legend(title = "mtry")) +
  theme_bw()

# TEST PREDICTIONS
# ==============================================================================
predic_READ_rf_RFE_50 <- predict(object = READ_rf_RFE_50, newdata = x_test_rfe_50)
predic_READ_rf_RFE_50
y_test <- as.factor(y_test)
caret::confusionMatrix(predic_READ_rf_RFE_50, y_test)


#===============================================================================
#=========================== SVM Kernel Radial =================================
#===============================================================================


# PARALLEL PROCESS
#===============================================================================

#install.packages("doMC", repos="http://R-Forge.R-project.org")

registerDoMC(cores = 11)

# HYPERPARAMETERS, NUMBER OF REPETITIONS AND SEEDS FOR EACH REPEAT
#===============================================================================
repetitions_boot <- 50

# Hyperparameters
hyperparameters <- expand.grid(sigma = c(0.0001, 0.001, 0.01),
                               C = c(1, 10, 50, 100, 250, 500, 700, 1000))

set.seed(86)
seeds <- vector(mode = "list", length = repetitions_boot + 1)
for (i in 1:repetitions_boot) {
  seeds[[i]] <- sample.int(1000, nrow(hyperparameters))
}
seeds[[repetitions_boot + 1]] <- sample.int(1000, 1)

# DEFINITION OF TRAINING
#===============================================================================
control_train <- trainControl(method = "boot", number = repetitions_boot,
                              seeds = seeds, returnResamp = "final",
                              verboseIter = TRUE, allowParallel = TRUE)

# FIT MODEL 
# ==============================================================================

set.seed(86)
READ_svmrad_RFE_50 <- caret::train(
  form = Class ~ .,
  data = x_train_rfe_50,
  method = "svmRadial",
  tuneGrid = hyperparameters,
  metric = "Accuracy",
  trControl = control_train
)

saveRDS(object = READ_svmrad_RFE_50, file = "READ_svmrad_RFE_50.rds")
registerDoMC(cores = 1)


READ_svmrad_RFE_50

# GRAPHIC REPRESENTATION
# ==============================================================================
ggplot(READ_svmrad_RFE_50, highlight = TRUE) +
  labs(title = "Evolution of the accuracy of the SVM model") +
  guides(color = guide_legend(title = "Sigma"),
         shape = guide_legend(title = "Sigma")) +
  theme_bw()

# TEST PREDICTIONS
# ==============================================================================
predic_READ_svmrad_RFE_50 <- predict(object = READ_svmrad_RFE_50, newdata = x_test_rfe_50)
predic_READ_svmrad_RFE_50
y_test <- as.factor(y_test)
caret::confusionMatrix(predic_READ_svmrad_RFE_50, y_test)


#===============================================================================
#=========================== Neural Network ====================================
#===============================================================================


# PARALLEL PROCESS
#===============================================================================

#install.packages("doMC", repos="http://R-Forge.R-project.org")

registerDoMC(cores = 11)

# HYPERPARAMETERS, NUMBER OF REPETITIONS AND SEEDS FOR EACH REPEAT
#===============================================================================
repetitions_boot <- 50

# Hyperparameters
hyperparameters <- expand.grid(size = c(2, 4, 8, 10, 16),
                               decay = c(0.01, 0.1))

set.seed(86)
seeds <- vector(mode = "list", length = repetitions_boot + 1)
for (i in 1:repetitions_boot) {
  seeds[[i]] <- sample.int(1000, nrow(hyperparameters))
}
seeds[[repetitions_boot + 1]] <- sample.int(1000, 1)

# DEFINITION OF TRAINING
#===============================================================================
control_train <- trainControl(method = "boot", number = repetitions_boot,
                              seeds = seeds, returnResamp = "final",
                              verboseIter = TRUE, allowParallel = TRUE)

# FIT MODEL 
# ==============================================================================

set.seed(86)
READ_nnet_RFE_50 <- caret::train(
  form = Class ~ .,
  data = x_train_rfe_50,
  method = "nnet",
  tuneGrid = hyperparameters,
  metric = "Accuracy",
  trControl = control_train,
  rang = c(-0.7, 0.7),
  MaxNWts = 21000,
  trace = FALSE
)

saveRDS(object = READ_nnet_RFE_50, file = "READ_nnet_RFE_50.rds")
registerDoMC(cores = 1)

READ_nnet_RFE_50

# GRAPHIC REPRESENTATION
# ==============================================================================
ggplot(READ_nnet_RFE_50, highlight = TRUE) +
  labs(title = "Evolution of the accuracy of the Neural Net model") +
  guides(color = guide_legend(title = "Decay"),
         shape = guide_legend(title = "Decay")) +
  theme_bw()

# TEST PREDICTIONS
# ==============================================================================
predic_READ_nnet_RFE_50 <- predict(object = READ_nnet_RFE_50, newdata = x_test_rfe_50)
predic_READ_nnet_RFE_50
y_test <- as.factor(y_test)
caret::confusionMatrix(predic_READ_nnet_RFE_50, y_test)



#===============================================================================
#=========================== Gradient boosting =================================
#===============================================================================


# PARALLEL PROCESS
#===============================================================================

#install.packages("doMC", repos="http://R-Forge.R-project.org")

registerDoMC(cores = 11)

# HYPERPARAMETERS, NUMBER OF REPETITIONS AND SEEDS FOR EACH REPEAT
#===============================================================================
repetitions_boot <- 50

# Hyperparameters
hyperparameters <- expand.grid(
  shrinkage = c(0.01, 0.1, 0.3),
  interaction.depth = c(1, 3, 5),
  n.minobsinnode = c(5, 10, 15),
  n.trees = c(500, 1000)
  # bag.fraction = c(0.65, 0.8, 1), 
  # optimal_trees = 0,               
  # min_RMSE = 0,                    
  # min_cor = 0
)

set.seed(86)
seeds <- vector(mode = "list", length = repetitions_boot + 1)
for (i in 1:repetitions_boot) {
  seeds[[i]] <- sample.int(1000, nrow(hyperparameters))
}
seeds[[repetitions_boot + 1]] <- sample.int(1000, 1)

# DEFINITION OF TRAINING
#===============================================================================
control_train <- trainControl(method = "boot", number = repetitions_boot,
                              seeds = seeds, returnResamp = "final",
                              verboseIter = TRUE, allowParallel = TRUE)

# FIT MODEL 
# ==============================================================================

set.seed(86)
READ_gbm_RFE_50 <- caret::train(
  form = Class ~ .,
  data = x_train_rfe_50,
  method = "gbm",
  tuneGrid = hyperparameters,
  metric = "Accuracy",
  preProcess = c("center", "scale"),
  trControl = control_train
)

saveRDS(object = READ_gbm_RFE_50, file = "READ_gbm_RFE_50.rds")
registerDoMC(cores = 1)

READ_gbm_RFE_50

# GRAPHIC REPRESENTATION
# ==============================================================================
ggplot(READ_gbm_RFE_50, highlight = TRUE) +
  labs(title = "Evolution of the accuracy of the GBM model") +
  guides(color = guide_legend(title = "Shrinkage"),
         shape = guide_legend(title = "Shrinkage")) +
  theme_bw()

# TEST PREDICTIONS
# ==============================================================================
predic_READ_gbm_RFE_50 <- predict(object = READ_gbm_RFE_50, newdata = x_test_rfe_50)
predic_READ_gbm_RFE_50
y_test <- as.factor(y_test)
caret::confusionMatrix(predic_READ_gbm_RFE_50, y_test)


#===============================================================================
#=========================== XGBM ==============================================
#===============================================================================


# PARALLEL PROCESS
#===============================================================================

#install.packages("doMC", repos="http://R-Forge.R-project.org")

registerDoMC(cores = 11)

# HYPERPARAMETERS, NUMBER OF REPETITIONS AND SEEDS FOR EACH REPEAT
#===============================================================================
repetitions_boot <- 50

# Hyperparameters
hyperparameters <- expand.grid(
  nrounds = c(500, 1000),
  eta = c(0.01, 0.001), 
  max_depth = c(2, 4, 6),
  gamma = 1,
  colsample_bytree = c(0.2, 0.4),
  min_child_weight = c(1, 5),
  subsample = 1
)


set.seed(86)
seeds <- vector(mode = "list", length = repetitions_boot + 1)
for (i in 1:repetitions_boot) {
  seeds[[i]] <- sample.int(1000, nrow(hyperparameters))
}
seeds[[repetitions_boot + 1]] <- sample.int(1000, 1)

# DEFINITION OF TRAINING
#===============================================================================
control_train <- trainControl(method = "boot", number = repetitions_boot,
                              seeds = seeds, returnResamp = "final",
                              verboseIter = TRUE, allowParallel = TRUE)

# FIT MODEL 
# ==============================================================================

set.seed(86)
READ_xgbm_RFE_50 <- caret::train(
  form = Class ~ .,
  data = x_train_rfe_50,
  method = "xgbTree",
  tuneGrid = hyperparameters,
  metric = "Accuracy",
  preProcess = c("center", "scale"),
  trControl = control_train
)

saveRDS(object = READ_xgbm_RFE_50, file = "READ_xgbm_RFE_50.rds")
registerDoMC(cores = 1)

READ_xgbm_RFE_50

# GRAPHIC REPRESENTATION
# ==============================================================================
plot(READ_xgbm_RFE_50, highlight = TRUE)

# TEST PREDICTIONS
# ==============================================================================
predic_READ_xgbm_RFE_50 <- predict(object = READ_xgbm_RFE_50, newdata = x_test_rfe_50)
predic_READ_xgbm_RFE_50
y_test <- as.factor(y_test)
caret::confusionMatrix(predic_READ_xgbm_RFE_50, y_test)


#===============================================================================
#=========================== glmnet ============================================
#===============================================================================


# PARALLEL PROCESS
#===============================================================================

#install.packages("doMC", repos="http://R-Forge.R-project.org")

registerDoMC(cores = 11)

# HYPERPARAMETERS, NUMBER OF REPETITIONS AND SEEDS FOR EACH REPEAT
#===============================================================================
repetitions_boot <- 50

# Hyperparameters
hyperparameters <- expand.grid(
  lambda = c(0, 1, 10, 100),
  alpha = c (0.1, 0.01, 0.001)
)

set.seed(86)
seeds <- vector(mode = "list", length = repetitions_boot + 1)
for (i in 1:repetitions_boot) {
  seeds[[i]] <- sample.int(1000, nrow(hyperparameters))
}
seeds[[repetitions_boot + 1]] <- sample.int(1000, 1)

# DEFINITION OF TRAINING
#===============================================================================
control_train <- trainControl(method = "boot", number = repetitions_boot,
                              seeds = seeds, returnResamp = "final",
                              verboseIter = TRUE, allowParallel = TRUE)

# FIT MODEL 
# ==============================================================================

set.seed(86)
READ_glmnet_RFE_50 <- caret::train(
  form = Class ~ .,
  data = x_train_rfe_50,
  method = "glmnet",
  family = "multinomial",
  tuneGrid = hyperparameters,
  metric = "Accuracy",
  trControl = control_train
)

saveRDS(object = READ_glmnet_RFE_50, file = "READ_glmnet_RFE_50.rds")
registerDoMC(cores = 1)

READ_glmnet_RFE_50

# GRAPHIC REPRESENTATION
# ==============================================================================
ggplot(READ_glmnet_RFE_50, highlight = TRUE) +
  labs(title = "Evolution of the accuracy of the GLMnet model") +
  guides(color = guide_legend(title = "Alpha"),
         shape = guide_legend(title = "Alpha")) +
  theme_bw()

# TEST PREDICTIONS
# ==============================================================================
predic_READ_glmnet_RFE_50 <- predict(object = READ_glmnet_RFE_50, newdata = x_test_rfe_50)
predic_READ_glmnet_RFE_50
y_test <- as.factor(y_test)
caret::confusionMatrix(predic_READ_glmnet_RFE_50, y_test)


#===============================================================================
#=========================== HDDA ==============================================
#===============================================================================


# PARALLEL PROCESS
#===============================================================================

#install.packages("doMC", repos="http://R-Forge.R-project.org")

registerDoMC(cores = 11)

# HYPERPARAMETERS, NUMBER OF REPETITIONS AND SEEDS FOR EACH REPEAT
#===============================================================================
repetitions_boot <- 50

# Hyperparameters
hyperparameters <- expand.grid(
  threshold = seq(0.1,0.9,by=0.1),
  model = "ALL")

set.seed(86)
seeds <- vector(mode = "list", length = repetitions_boot + 1)
for (i in 1:repetitions_boot) {
  seeds[[i]] <- sample.int(1000, nrow(hyperparameters))
}
seeds[[repetitions_boot + 1]] <- sample.int(1000, 1)

# DEFINITION OF TRAINING
#===============================================================================
control_train <- trainControl(method = "boot", number = repetitions_boot,
                              seeds = seeds, returnResamp = "final",
                              verboseIter = TRUE, allowParallel = TRUE)

# FIT MODEL 
# ==============================================================================

set.seed(86)
READ_hdda_RFE_50 <- caret::train(
  form = Class ~ .,
  data = x_train_rfe_50,
  method = "hdda",
  tuneGrid = hyperparameters,
  metric = "Accuracy",
  trControl = control_train
)

saveRDS(object = READ_hdda_RFE_50, file = "READ_hdda_RFE_50.rds")
registerDoMC(cores = 1)

READ_hdda_RFE_50

# GRAPHIC REPRESENTATION
# ==============================================================================
ggplot(READ_hdda_RFE_50, highlight = TRUE) +
  labs(title = "Evolution of the accuracy of the HDDA model") +
  guides(color = guide_legend(title = "threshold"),
         shape = guide_legend(title = "threshold")) +
  theme_bw()

# TEST PREDICTIONS
# ==============================================================================
predic_READ_hdda_RFE_50 <- predict(object = READ_hdda_RFE_50, newdata = x_test_rfe_50)
predic_READ_hdda_RFE_50
y_test <- as.factor(y_test)
caret::confusionMatrix(predic_READ_hdda_RFE_50, y_test)


#===============================================================================
#=========================== LogitBoost ========================================
#===============================================================================


# PARALLEL PROCESS
#===============================================================================

#install.packages("doMC", repos="http://R-Forge.R-project.org")

registerDoMC(cores = 11)

# HYPERPARAMETERS, NUMBER OF REPETITIONS AND SEEDS FOR EACH REPEAT
#===============================================================================
repetitions_boot <- 50

# Hyperparameters
hyperparameters <- expand.grid(
  nIter = c(25, 50, 100, 250)
)

set.seed(86)
seeds <- vector(mode = "list", length = repetitions_boot + 1)
for (i in 1:repetitions_boot) {
  seeds[[i]] <- sample.int(1000, nrow(hyperparameters))
}
seeds[[repetitions_boot + 1]] <- sample.int(1000, 1)

# DEFINITION OF TRAINING
#===============================================================================
control_train <- trainControl(method = "boot", number = repetitions_boot,
                              seeds = seeds, returnResamp = "final",
                              verboseIter = TRUE, allowParallel = TRUE)

# FIT MODEL 
# ==============================================================================

set.seed(86)
READ_logitBoost_RFE_50 <- caret::train(
  form = Class ~ .,
  data = x_train_rfe_50,
  method = "LogitBoost",
  tuneGrid = hyperparameters,
  metric = "Accuracy",
  trControl = control_train
)

saveRDS(object = READ_logitBoost_RFE_50, file = "READ_logitBoost_RFE_50.rds")
registerDoMC(cores = 1)

READ_logitBoost_RFE_50

# GRAPHIC REPRESENTATION
# ==============================================================================
ggplot(READ_logitBoost_RFE_50, highlight = TRUE) +
  labs(title = "Evolution of the accuracy of the Logitboost model") +
  guides(color = guide_legend(title = "nIter"),
         shape = guide_legend(title = "nIter")) +
  theme_bw()

# TEST PREDICTIONS
# ==============================================================================
predic_READ_logitBoost_RFE_50 <- predict(object = READ_logitBoost_RFE_50, newdata = x_test_rfe_50)
predic_READ_logitBoost_RFE_50
y_test <- as.factor(y_test)
caret::confusionMatrix(predic_READ_logitBoost_RFE_50, y_test)


#*******************************************************************************
#*************************** RFE 100 features **********************************
#*******************************************************************************


best_vars_100 <- RFE_results$control$functions$selectVar(selected_vars, 100)
x_train_rfe_100 <- x_train[, c("Class", best_vars_100)] 
x_test_rfe_100 <- x_test[, c("Class", best_vars_100)] 



#===============================================================================
#=========================== RANDOM FOREST ===================================== 
#===============================================================================


# PARALLEL PROCESS
#===============================================================================

#install.packages("doMC", repos="http://R-Forge.R-project.org")

registerDoMC(cores = 11)

# HYPERPARAMETERS, NUMBER OF REPETITIONS AND SEEDS FOR EACH REPEAT
#===============================================================================
repetitions_boot <- 50

# Hyperparameters
hyperparameters <- expand.grid(mtry = c(2, 5, 10, 50),
                               min.node.size = c(2, 3, 4, 5, 10),
                               splitrule = "gini")

set.seed(86)
seeds <- vector(mode = "list", length = repetitions_boot + 1)
for (i in 1:repetitions_boot) {
  seeds[[i]] <- sample.int(1000, nrow(hyperparameters))
}
seeds[[repetitions_boot + 1]] <- sample.int(1000, 1)

# DEFINITION OF TRAINING
#===============================================================================
control_train <- trainControl(method = "boot", number = repetitions_boot,
                              seeds = seeds, returnResamp = "final",
                              verboseIter = TRUE, allowParallel = TRUE)

# FIT MODEL 
# ==============================================================================
x_train$Class <- as.factor(x_train$Class)

set.seed(86)
READ_rf_RFE_100 <- caret::train(
  form = Class ~ .,
  data = x_train_rfe_100,
  method = "ranger",
  tuneGrid = hyperparameters,
  metric = "Accuracy",
  trControl = control_train,
  num.trees = 500)

saveRDS(object = READ_rf_RFE_100, file = "READ_rf_RFE_100.rds")
registerDoMC(cores = 1)


READ_rf_RFE_100

# GRAPHIC REPRESENTATION
# ==============================================================================
ggplot(READ_rf_RFE_100, highlight = TRUE) +
  labs(title = "Evolution of the accuracy of the Random Forest model") +
  guides(color = guide_legend(title = "mtry"),
         shape = guide_legend(title = "mtry")) +
  theme_bw()

# TEST PREDICTIONS
# ==============================================================================
predic_READ_rf_RFE_100 <- predict(object = READ_rf_RFE_100, newdata = x_test_rfe_100)
predic_READ_rf_RFE_100
y_test <- as.factor(y_test)
caret::confusionMatrix(predic_READ_rf_RFE_100, y_test)


#===============================================================================
#=========================== SVM Kernel Radial =================================
#===============================================================================


# PARALLEL PROCESS
#===============================================================================

#install.packages("doMC", repos="http://R-Forge.R-project.org")

registerDoMC(cores = 11)

# HYPERPARAMETERS, NUMBER OF REPETITIONS AND SEEDS FOR EACH REPEAT
#===============================================================================
repetitions_boot <- 50

# Hyperparameters
hyperparameters <- expand.grid(sigma = c(0.0001, 0.001, 0.01),
                               C = c(1, 10, 50, 100, 250, 500, 700, 1000))

set.seed(86)
seeds <- vector(mode = "list", length = repetitions_boot + 1)
for (i in 1:repetitions_boot) {
  seeds[[i]] <- sample.int(1000, nrow(hyperparameters))
}
seeds[[repetitions_boot + 1]] <- sample.int(1000, 1)

# DEFINITION OF TRAINING
#===============================================================================
control_train <- trainControl(method = "boot", number = repetitions_boot,
                              seeds = seeds, returnResamp = "final",
                              verboseIter = TRUE, allowParallel = TRUE)

# FIT MODEL 
# ==============================================================================

set.seed(86)
READ_svmrad_RFE_100 <- caret::train(
  form = Class ~ .,
  data = x_train_rfe_100,
  method = "svmRadial",
  tuneGrid = hyperparameters,
  metric = "Accuracy",
  trControl = control_train
)

saveRDS(object = READ_svmrad_RFE_100, file = "READ_svmrad_RFE_100.rds")
registerDoMC(cores = 1)


READ_svmrad_RFE_100

# GRAPHIC REPRESENTATION
# ==============================================================================
ggplot(READ_svmrad_RFE_100, highlight = TRUE) +
  labs(title = "Evolution of the accuracy of the SVM model") +
  guides(color = guide_legend(title = "Sigma"),
         shape = guide_legend(title = "Sigma")) +
  theme_bw()

# TEST PREDICTIONS
# ==============================================================================
predic_READ_svmrad_RFE_100 <- predict(object = READ_svmrad_RFE_100, newdata = x_test_rfe_100)
predic_READ_svmrad_RFE_100
y_test <- as.factor(y_test)
caret::confusionMatrix(predic_READ_svmrad_RFE_100, y_test)


#===============================================================================
#=========================== Neural Network ====================================
#===============================================================================


# PARALLEL PROCESS
#===============================================================================

#install.packages("doMC", repos="http://R-Forge.R-project.org")

registerDoMC(cores = 11)

# HYPERPARAMETERS, NUMBER OF REPETITIONS AND SEEDS FOR EACH REPEAT
#===============================================================================
repetitions_boot <- 50

# Hyperparameters
hyperparameters <- expand.grid(size = c(2, 4, 8, 10, 16, 20),
                               decay = c(0.01, 0.1))

set.seed(86)
seeds <- vector(mode = "list", length = repetitions_boot + 1)
for (i in 1:repetitions_boot) {
  seeds[[i]] <- sample.int(1000, nrow(hyperparameters))
}
seeds[[repetitions_boot + 1]] <- sample.int(1000, 1)

# DEFINITION OF TRAINING
#===============================================================================
control_train <- trainControl(method = "boot", number = repetitions_boot,
                              seeds = seeds, returnResamp = "final",
                              verboseIter = TRUE, allowParallel = TRUE)

# FIT MODEL 
# ==============================================================================

set.seed(86)
READ_nnet_RFE_100 <- caret::train(
  form = Class ~ .,
  data = x_train_rfe_100,
  method = "nnet",
  tuneGrid = hyperparameters,
  metric = "Accuracy",
  trControl = control_train,
  rang = c(-0.7, 0.7),
  MaxNWts = 21000,
  trace = FALSE
)

saveRDS(object = READ_nnet_RFE_100, file = "READ_nnet_RFE_100.rds")
registerDoMC(cores = 1)

READ_nnet_RFE_100

# GRAPHIC REPRESENTATION
# ==============================================================================
ggplot(READ_nnet_RFE_100, highlight = TRUE) +
  labs(title = "Evolution of the accuracy of the Neural Net model") +
  guides(color = guide_legend(title = "Decay"),
         shape = guide_legend(title = "Decay")) +
  theme_bw()

# TEST PREDICTIONS
# ==============================================================================
predic_READ_nnet_RFE_100 <- predict(object = READ_nnet_RFE_100, newdata = x_test_rfe_100)
predic_READ_nnet_RFE_100
y_test <- as.factor(y_test)
caret::confusionMatrix(predic_READ_nnet_RFE_100, y_test)



#===============================================================================
#=========================== Gradient boosting =================================
#===============================================================================


# PARALLEL PROCESS
#===============================================================================

#install.packages("doMC", repos="http://R-Forge.R-project.org")

registerDoMC(cores = 11)

# HYPERPARAMETERS, NUMBER OF REPETITIONS AND SEEDS FOR EACH REPEAT
#===============================================================================
repetitions_boot <- 50

# Hyperparameters
hyperparameters <- expand.grid(
  shrinkage = c(0.01, 0.1, 0.3),
  interaction.depth = c(1, 3, 5),
  n.minobsinnode = c(5, 10, 15),
  n.trees = c(500, 1000)
  # bag.fraction = c(0.65, 0.8, 1), 
  # optimal_trees = 0,               
  # min_RMSE = 0,                    
  # min_cor = 0
)

set.seed(86)
seeds <- vector(mode = "list", length = repetitions_boot + 1)
for (i in 1:repetitions_boot) {
  seeds[[i]] <- sample.int(1000, nrow(hyperparameters))
}
seeds[[repetitions_boot + 1]] <- sample.int(1000, 1)

# DEFINITION OF TRAINING
#===============================================================================
control_train <- trainControl(method = "boot", number = repetitions_boot,
                              seeds = seeds, returnResamp = "final",
                              verboseIter = TRUE, allowParallel = TRUE)

# FIT MODEL 
# ==============================================================================

set.seed(86)
READ_gbm_RFE_100 <- caret::train(
  form = Class ~ .,
  data = x_train_rfe_100,
  method = "gbm",
  tuneGrid = hyperparameters,
  metric = "Accuracy",
  preProcess = c("center", "scale"),
  trControl = control_train
)

saveRDS(object = READ_gbm_RFE_100, file = "READ_gbm_RFE_100.rds")
registerDoMC(cores = 1)

READ_gbm_RFE_100

# GRAPHIC REPRESENTATION
# ==============================================================================
ggplot(READ_gbm_RFE_100, highlight = TRUE) +
  labs(title = "Evolution of the accuracy of the GBM model") +
  guides(color = guide_legend(title = "Shrinkage"),
         shape = guide_legend(title = "Shrinkage")) +
  theme_bw()

# TEST PREDICTIONS
# ==============================================================================
predic_READ_gbm_RFE_100 <- predict(object = READ_gbm_RFE_100, newdata = x_test_rfe_100)
predic_READ_gbm_RFE_100
y_test <- as.factor(y_test)
caret::confusionMatrix(predic_READ_gbm_RFE_100, y_test)


#===============================================================================
#=========================== XGBM ==============================================
#===============================================================================


# PARALLEL PROCESS
#===============================================================================

#install.packages("doMC", repos="http://R-Forge.R-project.org")

registerDoMC(cores = 11)

# HYPERPARAMETERS, NUMBER OF REPETITIONS AND SEEDS FOR EACH REPEAT
#===============================================================================
repetitions_boot <- 50

# Hyperparameters
hyperparameters <- expand.grid(
  nrounds = c(500, 1000),
  eta = c(0.01, 0.001), 
  max_depth = c(2, 4, 6),
  gamma = 1,
  colsample_bytree = c(0.2, 0.4),
  min_child_weight = c(1, 5),
  subsample = 1
)


set.seed(86)
seeds <- vector(mode = "list", length = repetitions_boot + 1)
for (i in 1:repetitions_boot) {
  seeds[[i]] <- sample.int(1000, nrow(hyperparameters))
}
seeds[[repetitions_boot + 1]] <- sample.int(1000, 1)

# DEFINITION OF TRAINING
#===============================================================================
control_train <- trainControl(method = "boot", number = repetitions_boot,
                              seeds = seeds, returnResamp = "final",
                              verboseIter = TRUE, allowParallel = TRUE)

# FIT MODEL 
# ==============================================================================

set.seed(86)
READ_xgbm_RFE_100 <- caret::train(
  form = Class ~ .,
  data = x_train_rfe_100,
  method = "xgbTree",
  tuneGrid = hyperparameters,
  metric = "Accuracy",
  preProcess = c("center", "scale"),
  trControl = control_train
)

saveRDS(object = READ_xgbm_RFE_100, file = "READ_xgbm_RFE_100.rds")
registerDoMC(cores = 1)

READ_xgbm_RFE_100

# GRAPHIC REPRESENTATION
# ==============================================================================
plot(READ_xgbm_RFE_100, highlight = TRUE)

# TEST PREDICTIONS
# ==============================================================================
predic_READ_xgbm_RFE_100 <- predict(object = READ_xgbm_RFE_100, newdata = x_test_rfe_100)
predic_READ_xgbm_RFE_100
y_test <- as.factor(y_test)
caret::confusionMatrix(predic_READ_xgbm_RFE_100, y_test)


#===============================================================================
#=========================== glmnet ============================================
#===============================================================================


# PARALLEL PROCESS
#===============================================================================

#install.packages("doMC", repos="http://R-Forge.R-project.org")

registerDoMC(cores = 11)

# HYPERPARAMETERS, NUMBER OF REPETITIONS AND SEEDS FOR EACH REPEAT
#===============================================================================
repetitions_boot <- 50

# Hyperparameters
hyperparameters <- expand.grid(
  lambda = c(0, 1, 10, 100),
  alpha = c (0.1, 0.01, 0.001)
)

set.seed(86)
seeds <- vector(mode = "list", length = repetitions_boot + 1)
for (i in 1:repetitions_boot) {
  seeds[[i]] <- sample.int(1000, nrow(hyperparameters))
}
seeds[[repetitions_boot + 1]] <- sample.int(1000, 1)

# DEFINITION OF TRAINING
#===============================================================================
control_train <- trainControl(method = "boot", number = repetitions_boot,
                              seeds = seeds, returnResamp = "final",
                              verboseIter = TRUE, allowParallel = TRUE)

# FIT MODEL 
# ==============================================================================

set.seed(86)
READ_glmnet_RFE_100 <- caret::train(
  form = Class ~ .,
  data = x_train_rfe_100,
  method = "glmnet",
  family = "multinomial",
  tuneGrid = hyperparameters,
  metric = "Accuracy",
  trControl = control_train
)

saveRDS(object = READ_glmnet_RFE_100, file = "READ_glmnet_RFE_100.rds")
registerDoMC(cores = 1)

READ_glmnet_RFE_100

# GRAPHIC REPRESENTATION
# ==============================================================================
ggplot(READ_glmnet_RFE_100, highlight = TRUE) +
  labs(title = "Evolution of the accuracy of the GLMnet model") +
  guides(color = guide_legend(title = "Alpha"),
         shape = guide_legend(title = "Alpha")) +
  theme_bw()

# TEST PREDICTIONS
# ==============================================================================
predic_READ_glmnet_RFE_100 <- predict(object = READ_glmnet_RFE_100, newdata = x_test_rfe_100)
predic_READ_glmnet_RFE_100
y_test <- as.factor(y_test)
caret::confusionMatrix(predic_READ_glmnet_RFE_100, y_test)


#===============================================================================
#=========================== HDDA ==============================================
#===============================================================================


# PARALLEL PROCESS
#===============================================================================

#install.packages("doMC", repos="http://R-Forge.R-project.org")

registerDoMC(cores = 11)

# HYPERPARAMETERS, NUMBER OF REPETITIONS AND SEEDS FOR EACH REPEAT
#===============================================================================
repetitions_boot <- 50

# Hyperparameters
hyperparameters <- expand.grid(
  threshold = seq(0.1,0.9,by=0.1),
  model = "ALL")

set.seed(86)
seeds <- vector(mode = "list", length = repetitions_boot + 1)
for (i in 1:repetitions_boot) {
  seeds[[i]] <- sample.int(1000, nrow(hyperparameters))
}
seeds[[repetitions_boot + 1]] <- sample.int(1000, 1)

# DEFINITION OF TRAINING
#===============================================================================
control_train <- trainControl(method = "boot", number = repetitions_boot,
                              seeds = seeds, returnResamp = "final",
                              verboseIter = TRUE, allowParallel = TRUE)

# FIT MODEL 
# ==============================================================================

set.seed(86)
READ_hdda_RFE_100 <- caret::train(
  form = Class ~ .,
  data = x_train_rfe_100,
  method = "hdda",
  tuneGrid = hyperparameters,
  metric = "Accuracy",
  trControl = control_train
)

saveRDS(object = READ_hdda_RFE_100, file = "READ_hdda_RFE_100.rds")
registerDoMC(cores = 1)

READ_hdda_RFE_100

# GRAPHIC REPRESENTATION
# ==============================================================================
ggplot(READ_hdda_RFE_100, highlight = TRUE) +
  labs(title = "Evolution of the accuracy of the HDDA model") +
  guides(color = guide_legend(title = "threshold"),
         shape = guide_legend(title = "threshold")) +
  theme_bw()

# TEST PREDICTIONS
# ==============================================================================
predic_READ_hdda_RFE_100 <- predict(object = READ_hdda_RFE_100, newdata = x_test_rfe_100)
predic_READ_hdda_RFE_100
y_test <- as.factor(y_test)
caret::confusionMatrix(predic_READ_hdda_RFE_100, y_test)


#===============================================================================
#=========================== LogitBoost ========================================
#===============================================================================


# PARALLEL PROCESS
#===============================================================================

#install.packages("doMC", repos="http://R-Forge.R-project.org")

registerDoMC(cores = 11)

# HYPERPARAMETERS, NUMBER OF REPETITIONS AND SEEDS FOR EACH REPEAT
#===============================================================================
repetitions_boot <- 50

# Hyperparameters
hyperparameters <- expand.grid(
  nIter = c(25, 50, 100, 250)
)

set.seed(86)
seeds <- vector(mode = "list", length = repetitions_boot + 1)
for (i in 1:repetitions_boot) {
  seeds[[i]] <- sample.int(1000, nrow(hyperparameters))
}
seeds[[repetitions_boot + 1]] <- sample.int(1000, 1)

# DEFINITION OF TRAINING
#===============================================================================
control_train <- trainControl(method = "boot", number = repetitions_boot,
                              seeds = seeds, returnResamp = "final",
                              verboseIter = TRUE, allowParallel = TRUE)

# FIT MODEL 
# ==============================================================================

set.seed(86)
READ_logitBoost_RFE_100 <- caret::train(
  form = Class ~ .,
  data = x_train_rfe_100,
  method = "LogitBoost",
  tuneGrid = hyperparameters,
  metric = "Accuracy",
  trControl = control_train
)

saveRDS(object = READ_logitBoost_RFE_100, file = "READ_logitBoost_RFE_100.rds")
registerDoMC(cores = 1)

READ_logitBoost_RFE_100

# GRAPHIC REPRESENTATION
# ==============================================================================
ggplot(READ_logitBoost_RFE_100, highlight = TRUE) +
  labs(title = "Evolution of the accuracy of the Logitboost model") +
  guides(color = guide_legend(title = "nIter"),
         shape = guide_legend(title = "nIter")) +
  theme_bw()

# TEST PREDICTIONS
# ==============================================================================
predic_READ_logitBoost_RFE_100 <- predict(object = READ_logitBoost_RFE_100, newdata = x_test_rfe_100)
predic_READ_logitBoost_RFE_100
y_test <- as.factor(y_test)
caret::confusionMatrix(predic_READ_logitBoost_RFE_100, y_test)






