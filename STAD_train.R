
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
library(ROCit)
library(precrec)
library(pROC)
library(Epi)

########################################################################
####################### Cervical Data ##################################
########################################################################


filepath <- system.file("extdata/cervical.txt", package = "MLSeq")
cervical <- STAD.table(filepath, header=TRUE)
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
# show(gbm_fit)
# trained(gbm_fit)
# plot(gbm_fit)

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


########################################################################
####################### TCGA Data ######################################
########################################################################

GDCprojects = getGDCprojects()
names_projects <- GDCprojects[c("project_id", "name")]

# Stomach Adenocarcinoma
TCGAbiolinks:::getProjectSummary("TCGA-STAD")

query_TCGA = GDCquery(
  project = "TCGA-STAD",
  data.category = "Transcriptome Profiling", # parameter enforced by GDCquery
  experimental.strategy = "RNA-Seq",
  workflow.type = "HTSeq - Counts")

STAD_res = getResults(query_TCGA) # make results as table
head(STAD_res) # data of the first 6 patients.
colnames(STAD_res) # columns present in the table

head(STAD_res$sample_type) # first 6 types of tissue.
table(STAD_res$sample_type) # summary of distinct tissues types present in this study

GDCdownload(query = query_TCGA) # download the data into my current workdirectory 

tcga_data = GDCprepare(query_TCGA) # load the data 
dim(tcga_data)

# colData () allows us to access the clinical data associated with our samples. 
# The colnames () and rownames () functions can be used to extract the column and row names of a given table, respectively.

colnames(colData(tcga_data))

# tables 

table(tcga_data@colData$gender)
table(tcga_data@colData$vital_status)
table(tcga_data@colData$tissue_or_organ_of_origin)
table(tcga_data@colData$race)
addmargins(round(prop.table(table(tcga_data@colData$vital_status,tcga_data@colData$gender))*100,3))
addmargins(round(prop.table(table(tcga_data@colData$vital_status,tcga_data@colData$race))*100,3))
addmargins(round(prop.table(table(tcga_data@colData$tissue_or_organ_of_origin,tcga_data@colData$gender))*100,3))
addmargins(round(prop.table(table(tcga_data@colData$tissue_or_organ_of_origin,tcga_data@colData$race))*100,3))
addmargins(round(prop.table(table(tcga_data@colData$vital_status,tcga_data@colData$tissue_or_organ_of_origin))*100,3))

# assay function 
head(assay(tcga_data)[,1:10]) # expression of first 6 genes and first 10 samples
head(rowData(tcga_data))     # ensembl id and gene id of the first 6 genes.

# Save the data as a file, if you need it later, you can just load this file
# instead of having to run the whole pipeline again

# saveRDS(object = tcga_data,
#         file = "tcga_data_STAD.RDS",
#         compress = FALSE)

# load data with the next sentence
tcga_data = STADRDS(file = "tcga_data_STAD.RDS")


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

# x_train$Class <- NULL
# x_test$Class <- NULL

prop.table(table(x_train$Class)) %>% round(3)
prop.table(table(x_test$Class)) %>% round(3)

summary(x_train$ENSG00000000003)
head(x_train[,1:6])

# ## Genes with variance close to zero
# 
# # We proceed to eliminate those genes whose maximum expression 
# # does not exceed 5 times the minimum expression (max/min) < 5 
# # and whose absolute difference between maximum and minimum does not exceed 500 units (max − min < 500)
# 
# 
# filter_variance <- function(x){
#   # This function returns TRUE for all non-numeric columns and, in the case
#   # of the numerical ones, those that exceed the minimum variance conditions.
#   if(is.numeric(x)){
#     maximum <- max(x)
#     minimum <- min(x)
#     ratio <- maximum / minimum
#     range <- maximum - minimum
#     return(ratio >= 5 & range >= 500)
#   }
#   
#   else{
#     return(TRUE)
#   }
# }
# 
# # The columns that meet the condition are identified
# genes_variance <- map_lgl(.x = x_train, .f = filter_variance)
# 
# # Number of columns (genes) excluded
# sum(genes_variance == FALSE)

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
# VERSIÓN PARALELIZADA DE BOOTSTRAPPING PARA FILTRADO POR ANOVA
#===============================================================================

n_cores <- parallel::detectCores() - 1
registerDoParallel(makeCluster(n_cores))
getDoParWorkers()


# Número de iteraciones bootstrapping
n_boot <- 100

# Semillas para que los muestreos sean reproducibles
set.seed(86)
seeds = sample.int(1000, size = n_boot)


resultados_anova_pvalue <- foreach(i = 1:n_boot) %dopar% {
  library(dplyr)
  library(purrr)
  # Se crea una muestra por bootstrapping
  set.seed(seeds[i])
  indices <- sample(1:nrow(x_train), size = nrow(x_train), replace = TRUE)
  pseudo_muestra <- x_train[indices, ]
  
  # Se calculan los p-values para la nueva muestra 
  p_values <- dplyr::select(pseudo_muestra,-Class) %>%
    purrr::map_dbl(.f = custom_anova, y = pseudo_muestra$Class) 
  
  # Se devuelven los p-value
  p_values
}

options(cores = 1)


# Los resultados almacenados en forma de lista se convierten en dataframe
names(resultados_anova_pvalue) <-  paste("resample", 1:n_boot, sep = "_") 
resultados_anova_pvalue <- data.frame(resultados_anova_pvalue)
resultados_anova_pvalue <- resultados_anova_pvalue %>% rownames_to_column(var = "gen")
resultados_anova_pvalue <- resultados_anova_pvalue %>%
  mutate(pvalue_medio = rowMeans(resultados_anova_pvalue[, -1])) %>%
  arrange(pvalue_medio)

# Se guarda en disco el objeto cSTADo para no tener que repetir de nuevo toda la
# computación.
saveRDS(object = resultados_anova_pvalue, file = "resultados_anova_pvalue.rds")

resultados_anova_pvalue %>% dplyr::select(1,2,3,4) %>% head()

# Se filtran los 100, 50 y 500 genes identificados como más relevantes mediante anova
filtrado_anova_pvalue_500 <-  resultados_anova_pvalue %>% pull(gen) %>% head(500)
filtrado_anova_pvalue_100 <- resultados_anova_pvalue %>% pull(gen) %>% head(100)
filtrado_anova_pvalue_50  <- resultados_anova_pvalue %>% pull(gen) %>% head(50)



#*******************************************************************************
#*************************** ANOVA PVALUE 100 feautures ************************
#*******************************************************************************


#===============================================================================
#=========================== RANDOM FOREST ===================================== 
#===============================================================================


# PARALELIZACIÓN DE PROCESO
#===============================================================================

#install.packages("doMC", repos="http://R-Forge.R-project.org")

registerDoMC(cores = 11)

# HIPERPARÁMETROS, NÚMERO DE REPETICIONES Y SEMILLAS PARA CADA REPETICIÓN
#===============================================================================
repeticiones_boot <- 50

# Hiperparámetros
hiperparametros <- expand.grid(mtry = c(2, 5, 10, 50),
                               min.node.size = c(2, 3, 4, 5, 10),
                               splitrule = "gini")

set.seed(86)
seeds <- vector(mode = "list", length = repeticiones_boot + 1)
for (i in 1:repeticiones_boot) {
  seeds[[i]] <- sample.int(1000, nrow(hiperparametros))
}
seeds[[repeticiones_boot + 1]] <- sample.int(1000, 1)

# DEFINICIÓN DEL ENTRENAMIENTO
#===============================================================================
control_train <- trainControl(method = "boot", number = repeticiones_boot,
                              seeds = seeds, returnResamp = "final",
                              verboseIter = TRUE, allowParallel = TRUE)

# AJUSTE DEL MODELO
# ==============================================================================
x_train$Class <- as.factor(x_train$Class)

set.seed(86)
rf_pvalue_100 <- caret::train(
  form = Class ~ .,
  data = x_train[c("Class", filtrado_anova_pvalue_100)],
  method = "ranger",
  tuneGrid = hiperparametros,
  metric = "Accuracy",
  trControl = control_train,
  # Número de árboles ajustados
  num.trees = 500)

saveRDS(object = rf_pvalue_100, file = "rf_pvalue_100.rds")
registerDoMC(cores = 1)


rf_pvalue_100

# REPRESENTACIÓN GRÁFICA
# ==============================================================================
ggplot(rf_pvalue_100, highlight = TRUE) +
  labs(title = "Evolución del accuracy del modelo Random Forest") +
  guides(color = guide_legend(title = "mtry"),
         shape = guide_legend(title = "mtry")) +
  theme_bw()

# PREDICCIONES SOBRE TEST 
# ==============================================================================
predic_rd_pvalue_100 <- predict(object = rf_pvalue_100, newdata = x_test)
predic_rd_pvalue_100
y_test <- as.factor(y_test)
caret::confusionMatrix(predic_rd_pvalue_100, y_test)



#===============================================================================
#=========================== SVM Kernel Radial ================================= 
#===============================================================================


# PARALELIZACIÓN DE PROCESO
#===============================================================================

registerDoMC(cores = 11)

# HIPERPARÁMETROS, NÚMERO DE REPETICIONES Y SEMILLAS PARA CADA REPETICIÓN
#===============================================================================
repeticiones_boot <- 50

# Hiperparámetros
hiperparametros <- expand.grid(sigma = c(0.0001, 0.001, 0.01),
                               C = c(1, 10, 50, 100, 250, 500, 700, 1000))
set.seed(123)
seeds <- vector(mode = "list", length = repeticiones_boot + 1)
for (i in 1:repeticiones_boot) {
  seeds[[i]] <- sample.int(1000, nrow(hiperparametros))
}
seeds[[repeticiones_boot + 1]] <- sample.int(1000, 1)

# DEFINICIÓN DEL ENTRENAMIENTO
#===============================================================================
control_train <- trainControl(method = "boot", number = repeticiones_boot,
                              seeds = seeds, returnResamp = "final",
                              verboseIter = TRUE, allowParallel = TRUE)
# AJUSTE DEL MODELO
# ==============================================================================
set.seed(86)
svmrad_pvalue_100 <- caret::train(
  form = Class ~ .,
  data = x_train[c("Class", filtrado_anova_pvalue_100)],
  method = "svmRadial",
  tuneGrid = hiperparametros,
  metric = "Accuracy",
  trControl = control_train
)
registerDoMC(cores = 1)
saveRDS(object = svmrad_pvalue_100, file = "svmrad_pvalue_100.rds")

svmrad_pvalue_100

# REPRESENTACIÓN GRÁFICA
# ==============================================================================
ggplot(svmrad_pvalue_100, highlight = TRUE) +
  labs(title = "Evolución del accuracy del modelo SVM Radial") +
  theme_bw()

# PREDICCIONES SOBRE TEST 
# ==============================================================================
predic_SVM_pvalue_100 <- predict(object = svmrad_pvalue_100, newdata = x_test)
caret::confusionMatrix(predic_SVM_pvalue_100, y_test)



#===============================================================================
#=========================== Neural Network  =================================== 
#===============================================================================

# PARALELIZACIÓN DE PROCESO
#===============================================================================

registerDoMC(cores = 11)

# HIPERPARÁMETROS, NÚMERO DE REPETICIONES Y SEMILLAS PARA CADA REPETICIÓN
#===============================================================================
repeticiones_boot <- 50

# Hiperparámetros
hiperparametros <- expand.grid(size = c(5, 10, 15, 20, 40),
                               decay = c(0.01, 0.1))

set.seed(86)
seeds <- vector(mode = "list", length = repeticiones_boot + 1)
for (i in 1:repeticiones_boot) {
  seeds[[i]] <- sample.int(1000, nrow(hiperparametros))
}
seeds[[repeticiones_boot + 1]] <- sample.int(1000, 1)

# DEFINICIÓN DEL ENTRENAMIENTO
#===============================================================================
control_train <- trainControl(method = "boot", number = repeticiones_boot,
                              seeds = seeds, returnResamp = "final",
                              verboseIter = TRUE, allowParallel = TRUE)

# AJUSTE DEL MODELO
# ==============================================================================
set.seed(86)
nnet_pvalue_100 <- caret::train(
  form = Class ~ .,
  data = x_train[c("Class", filtrado_anova_pvalue_100)],
  method = "nnet",
  tuneGrid = hiperparametros,
  metric = "Accuracy",
  trControl = control_train,
  # Rango de inicialización de los pesos
  rang = c(-0.7, 0.7),
  # Número máximo de pesos
  MaxNWts = 10000,
  # Para que no se muestre cada iteración por pantalla
  trace = FALSE
)

saveRDS(object = nnet_pvalue_100, file = "nnet_pvalue_100.rds")
registerDoMC(cores = 1)

nnet_pvalue_100

# REPRESENTACIÓN GRÁFICA
# ==============================================================================
ggplot(nnet_pvalue_100, highlight = TRUE) +
  labs(title = "Evolución del accuracy del modelo NNET") +
  theme_bw()

# PREDICCIONES SOBRE TEST 
# ==============================================================================
predic_nnet_pvalue_100 <- predict(object = nnet_pvalue_100, newdata = x_test)
caret::confusionMatrix(predic_nnet_pvalue_100, y_test)



#===============================================================================
#=========================== GBM  ============================================== 
#===============================================================================

# PARALELIZACIÓN DE PROCESO
#===============================================================================

registerDoMC(cores = 11)

# HIPERPARÁMETROS, NÚMERO DE REPETICIONES Y SEMILLAS PARA CADA REPETICIÓN
#===============================================================================
repeticiones_boot <- 50

# Hiperparámetros
hiperparametros <- expand.grid(
  shrinkage = c(0.01, 0.1, 0.3),
  interaction.depth = c(1, 3, 5),
  n.minobsinnode = c(5, 10, 15),
  n.trees = c(100,250,500,1000)
  # bag.fraction = c(0.65, 0.8, 1), 
  # optimal_trees = 0,               
  # min_RMSE = 0,                    
  # min_cor = 0
)

set.seed(86)
seeds <- vector(mode = "list", length = repeticiones_boot + 1)
for (i in 1:repeticiones_boot) {
  seeds[[i]] <- sample.int(1000, nrow(hiperparametros))
}
seeds[[repeticiones_boot + 1]] <- sample.int(1000, 1)

# DEFINICIÓN DEL ENTRENAMIENTO
#===============================================================================
control_train <- trainControl(method = "boot", number = repeticiones_boot,
                              seeds = seeds, returnResamp = "final",
                              verboseIter = TRUE, allowParallel = TRUE)

# AJUSTE DEL MODELO
# ==============================================================================
set.seed(86)
gbm_pvalue_100 <- caret::train(
  form = Class ~ .,
  data = x_train[c("Class", filtrado_anova_pvalue_100)],
  method = "gbm",
  tuneGrid = hiperparametros,
  metric = "Accuracy",
  trControl = control_train
)

saveRDS(object = gbm_pvalue_100, file = "gbm_pvalue_100.rds")
registerDoMC(cores = 1)

gbm_pvalue_100

# REPRESENTACIÓN GRÁFICA
# ==============================================================================
ggplot(gbm_pvalue_100, highlight = TRUE) +
  labs(title = "Evolución del accuracy del modelo GBM") +
  theme_bw()

# PREDICCIONES SOBRE TEST 
# ==============================================================================
predic_gbm_pvalue_100 <- predict(object = gbm_pvalue_100, newdata = x_test)
caret::confusionMatrix(predic_gbm_pvalue_100, y_test)



#===============================================================================
#=========================== XGBM  ============================================= 
#===============================================================================

# PARALELIZACIÓN DE PROCESO
#===============================================================================

registerDoMC(cores = 11)

# HIPERPARÁMETROS, NÚMERO DE REPETICIONES Y SEMILLAS PARA CADA REPETICIÓN
#===============================================================================
repeticiones_boot <- 50

hiperparametros <- expand.grid(
  nrounds = c(500, 1000),
  eta = c(0.01, 0.001), 
  max_depth = c(2, 4, 6),
  gamma = 1,
  colsample_bytree = c(0.2, 0.4),
  min_child_weight = c(1, 5),
  subsample = 1
)



set.seed(86)
seeds <- vector(mode = "list", length = repeticiones_boot + 1)
for (i in 1:repeticiones_boot) {
  seeds[[i]] <- sample.int(1000, nrow(hiperparametros))
}
seeds[[repeticiones_boot + 1]] <- sample.int(1000, 1)

# DEFINICIÓN DEL ENTRENAMIENTO
#===============================================================================
control_train <- trainControl(method = "boot", number = repeticiones_boot,
                              seeds = seeds, returnResamp = "final",
                              verboseIter = TRUE, allowParallel = TRUE)

# AJUSTE DEL MODELO
# ==============================================================================
set.seed(86)
xgbm_pvalue_100 <- caret::train(
  form = Class ~ .,
  data = x_train[c("Class", filtrado_anova_pvalue_100)],
  method = "xgbTree",
  tuneGrid = hiperparametros,
  metric = "Accuracy",
  preProcess = c("center", "scale"),
  trControl = control_train
)


saveRDS(object = xgbm_pvalue_100, file = "xgbm_pvalue_100.rds")
registerDoMC(cores = 1)

xgbm_pvalue_100

# REPRESENTACIÓN GRÁFICA
# ==============================================================================
ggplot(xgbm_pvalue_100, highlight = TRUE) +
  labs(title = "Evolución del accuracy del modelo XGBM") +
  theme_bw()

# PREDICCIONES SOBRE TEST 
# ==============================================================================
predic_xgbm_pvalue_100 <- predict(object = xgbm_pvalue_100, newdata = x_test)
caret::confusionMatrix(predic_xgbm_pvalue_100, y_test)


#===============================================================================
#=========================== glmnet  ===========================================
#===============================================================================

# PARALELIZACIÓN DE PROCESO
#===============================================================================

registerDoMC(cores = 11)

# HIPERPARÁMETROS, NÚMERO DE REPETICIONES Y SEMILLAS PARA CADA REPETICIÓN
#===============================================================================
repeticiones_boot <- 50

hiperparametros <- expand.grid(
  lambda = c(0, 1, 10, 100),
  alpha = c (0.1, 0.01, 0.001)
)



set.seed(86)
seeds <- vector(mode = "list", length = repeticiones_boot + 1)
for (i in 1:repeticiones_boot) {
  seeds[[i]] <- sample.int(1000, nrow(hiperparametros))
}
seeds[[repeticiones_boot + 1]] <- sample.int(1000, 1)

# DEFINICIÓN DEL ENTRENAMIENTO
#===============================================================================
control_train <- trainControl(method = "boot", number = repeticiones_boot,
                              seeds = seeds, returnResamp = "final",
                              verboseIter = TRUE, allowParallel = TRUE)

# AJUSTE DEL MODELO
# ==============================================================================
set.seed(86)
glmnet_pvalue_100 <- caret::train(
  form = Class ~ .,
  data = x_train[c("Class", filtrado_anova_pvalue_100)],
  method = "glmnet",
  family = "binomial",
  tuneGrid = hiperparametros,
  metric = "Accuracy",
  trControl = control_train
)


saveRDS(object = glmnet_pvalue_100, file = "glmnet_pvalue_100.rds")
registerDoMC(cores = 1)

glmnet_pvalue_100

# REPRESENTACIÓN GRÁFICA
# ==============================================================================
ggplot(glmnet_pvalue_100, highlight = TRUE) +
  labs(title = "Evolución del accuracy del modelo glmnet") +
  theme_bw()

# PREDICCIONES SOBRE TEST 
# ==============================================================================
predic_glmnet_pvalue_100 <- predict(object = glmnet_pvalue_100, newdata = x_test)
caret::confusionMatrix(predic_glmnet_pvalue_100, y_test)


#===============================================================================
#=========================== hdda  =============================================
#===============================================================================

# PARALELIZACIÓN DE PROCESO
#===============================================================================

registerDoMC(cores = 11)

# HIPERPARÁMETROS, NÚMERO DE REPETICIONES Y SEMILLAS PARA CADA REPETICIÓN
#===============================================================================
repeticiones_boot <- 50

hiperparametros <- expand.grid(
  threshold = seq(0.1,0.9,by=0.1),
  model = "ALL"
)



set.seed(86)
seeds <- vector(mode = "list", length = repeticiones_boot + 1)
for (i in 1:repeticiones_boot) {
  seeds[[i]] <- sample.int(1000, nrow(hiperparametros))
}
seeds[[repeticiones_boot + 1]] <- sample.int(1000, 1)

# DEFINICIÓN DEL ENTRENAMIENTO
#===============================================================================
control_train <- trainControl(method = "boot", number = repeticiones_boot,
                              seeds = seeds, returnResamp = "final",
                              verboseIter = TRUE, allowParallel = TRUE)

# AJUSTE DEL MODELO
# ==============================================================================
set.seed(86)
hdda_pvalue_100 <- caret::train(
  form = Class ~ .,
  data = x_train[c("Class", filtrado_anova_pvalue_100)],
  method = "hdda",
  tuneGrid = hiperparametros,
  metric = "Accuracy",
  trControl = control_train
)


saveRDS(object = hdda_pvalue_100, file = "hdda_pvalue_100.rds")
registerDoMC(cores = 1)

hdda_pvalue_100

# REPRESENTACIÓN GRÁFICA
# ==============================================================================
ggplot(hdda_pvalue_100, highlight = TRUE) +
  labs(title = "Evolución del accuracy del modelo hdda") +
  theme_bw()

# PREDICCIONES SOBRE TEST 
# ==============================================================================
predic_hdda_pvalue_100 <- predict(object = hdda_pvalue_100, newdata = x_test)
caret::confusionMatrix(predic_hdda_pvalue_100, y_test)


#===============================================================================
#=========================== LogitBoost  =======================================
#===============================================================================

# PARALELIZACIÓN DE PROCESO
#===============================================================================

registerDoMC(cores = 11)

# HIPERPARÁMETROS, NÚMERO DE REPETICIONES Y SEMILLAS PARA CADA REPETICIÓN
#===============================================================================
repeticiones_boot <- 50

hiperparametros <- expand.grid(
  nIter = c(25, 50, 100,250)
)



set.seed(86)
seeds <- vector(mode = "list", length = repeticiones_boot + 1)
for (i in 1:repeticiones_boot) {
  seeds[[i]] <- sample.int(1000, nrow(hiperparametros))
}
seeds[[repeticiones_boot + 1]] <- sample.int(1000, 1)

# DEFINICIÓN DEL ENTRENAMIENTO
#===============================================================================
control_train <- trainControl(method = "boot", number = repeticiones_boot,
                              seeds = seeds, returnResamp = "final",
                              verboseIter = TRUE, allowParallel = TRUE)

# AJUSTE DEL MODELO
# ==============================================================================
set.seed(86)
LogitBoost_pvalue_100 <- caret::train(
  form = Class ~ .,
  data = x_train[c("Class", filtrado_anova_pvalue_100)],
  method = "LogitBoost",
  tuneGrid = hiperparametros,
  metric = "Accuracy",
  trControl = control_train
)


saveRDS(object = LogitBoost_pvalue_100, file = "LogitBoost_pvalue_100.rds")
registerDoMC(cores = 1)

LogitBoost_pvalue_100

# REPRESENTACIÓN GRÁFICA
# ==============================================================================
ggplot(LogitBoost_pvalue_100, highlight = TRUE) +
  labs(title = "Evolución del accuracy del modelo LogitBoost") +
  theme_bw()

# PREDICCIONES SOBRE TEST 
# ==============================================================================
predic_LogitBoost_pvalue_100 <- predict(object = LogitBoost_pvalue_100, newdata = x_test)
caret::confusionMatrix(predic_LogitBoost_pvalue_100, y_test)


#===============================================================================
#=========================== Stacking Ensemble  ================================
#===============================================================================

# DEFINICIÓN DEL ENTRENAMIENTO
#===============================================================================


trainControl <- trainControl(method = "repeatedcv", number = 5, repeats = 3,
                             savePredictions = "all", classProbs = TRUE,
                             index = createFolds(stacking_df$Class, 5),verboseIter = TRUE)


# Run multiplealgorithms in one call.
algorithm <- c("nnet","rf","svmRadial","gbm","xgbTree","hdda","LogitBoost")


stacking_df <- cbind(x_train[,filtrado_anova_pvalue_100],x_train$Class)
names(stacking_df) [101] <- "Class"
levels(stacking_df$Class)
levels(stacking_df$Class) <- make.names(levels(stacking_df$Class))

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
#*************************** ANOVA PVALUE 50 feautures ************************
#*******************************************************************************


#===============================================================================
#=========================== RANDOM FOREST ===================================== 
#===============================================================================


# PARALELIZACIÓN DE PROCESO
#===============================================================================

#install.packages("doMC", repos="http://R-Forge.R-project.org")

registerDoMC(cores = 11)

# HIPERPARÁMETROS, NÚMERO DE REPETICIONES Y SEMILLAS PARA CADA REPETICIÓN
#===============================================================================
repeticiones_boot <- 50

# Hiperparámetros
hiperparametros <- expand.grid(mtry = c(2, 5, 10, 50),
                               min.node.size = c(2, 3, 4, 5, 10),
                               splitrule = "gini")

set.seed(86)
seeds <- vector(mode = "list", length = repeticiones_boot + 1)
for (i in 1:repeticiones_boot) {
  seeds[[i]] <- sample.int(1000, nrow(hiperparametros))
}
seeds[[repeticiones_boot + 1]] <- sample.int(1000, 1)

# DEFINICIÓN DEL ENTRENAMIENTO
#===============================================================================
control_train <- trainControl(method = "boot", number = repeticiones_boot,
                              seeds = seeds, returnResamp = "final",
                              verboseIter = TRUE, allowParallel = TRUE)

# AJUSTE DEL MODELO
# ==============================================================================
x_train$Class <- as.factor(x_train$Class)

set.seed(86)
rf_pvalue_50 <- caret::train(
  form = Class ~ .,
  data = x_train[c("Class", filtrado_anova_pvalue_50)],
  method = "ranger",
  tuneGrid = hiperparametros,
  metric = "Accuracy",
  trControl = control_train,
  # Número de árboles ajustados
  num.trees = 500)

saveRDS(object = rf_pvalue_50, file = "rf_pvalue_50.rds")
registerDoMC(cores = 1)


rf_pvalue_50

# REPRESENTACIÓN GRÁFICA
# ==============================================================================
ggplot(rf_pvalue_50, highlight = TRUE) +
  labs(title = "Evolución del accuracy del modelo Random Forest") +
  guides(color = guide_legend(title = "mtry"),
         shape = guide_legend(title = "mtry")) +
  theme_bw()

# PREDICCIONES SOBRE TEST 
# ==============================================================================
predic_rd_pvalue_50 <- predict(object = rf_pvalue_50, newdata = x_test)
predic_rd_pvalue_50
y_test <- as.factor(y_test)
caret::confusionMatrix(predic_rd_pvalue_50, y_test)



#===============================================================================
#=========================== SVM Kernel Radial ================================= 
#===============================================================================


# PARALELIZACIÓN DE PROCESO
#===============================================================================

registerDoMC(cores = 11)

# HIPERPARÁMETROS, NÚMERO DE REPETICIONES Y SEMILLAS PARA CADA REPETICIÓN
#===============================================================================
repeticiones_boot <- 50

# Hiperparámetros
hiperparametros <- expand.grid(sigma = c(0.0001, 0.001, 0.01),
                               C = c(1, 10, 50, 100, 250, 500, 700, 1000))
set.seed(123)
seeds <- vector(mode = "list", length = repeticiones_boot + 1)
for (i in 1:repeticiones_boot) {
  seeds[[i]] <- sample.int(1000, nrow(hiperparametros))
}
seeds[[repeticiones_boot + 1]] <- sample.int(1000, 1)

# DEFINICIÓN DEL ENTRENAMIENTO
#===============================================================================
control_train <- trainControl(method = "boot", number = repeticiones_boot,
                              seeds = seeds, returnResamp = "final",
                              verboseIter = TRUE, allowParallel = TRUE)
# AJUSTE DEL MODELO
# ==============================================================================
set.seed(86)
svmrad_pvalue_50 <- caret::train(
  form = Class ~ .,
  data = x_train[c("Class", filtrado_anova_pvalue_50)],
  method = "svmRadial",
  tuneGrid = hiperparametros,
  metric = "Accuracy",
  trControl = control_train
)
registerDoMC(cores = 1)
saveRDS(object = svmrad_pvalue_50, file = "svmrad_pvalue_50.rds")

svmrad_pvalue_50

# REPRESENTACIÓN GRÁFICA
# ==============================================================================
ggplot(svmrad_pvalue_50, highlight = TRUE) +
  labs(title = "Evolución del accuracy del modelo SVM Radial") +
  theme_bw()

# PREDICCIONES SOBRE TEST 
# ==============================================================================
predic_SVM_pvalue_50 <- predict(object = svmrad_pvalue_50, newdata = x_test)
caret::confusionMatrix(predic_SVM_pvalue_50, y_test)



#===============================================================================
#=========================== Neural Network  =================================== 
#===============================================================================

# PARALELIZACIÓN DE PROCESO
#===============================================================================

registerDoMC(cores = 11)

# HIPERPARÁMETROS, NÚMERO DE REPETICIONES Y SEMILLAS PARA CADA REPETICIÓN
#===============================================================================
repeticiones_boot <- 50

# Hiperparámetros
hiperparametros <- expand.grid(size = c(5, 10, 15, 20, 40),
                               decay = c(0.01, 0.1))

set.seed(86)
seeds <- vector(mode = "list", length = repeticiones_boot + 1)
for (i in 1:repeticiones_boot) {
  seeds[[i]] <- sample.int(1000, nrow(hiperparametros))
}
seeds[[repeticiones_boot + 1]] <- sample.int(1000, 1)

# DEFINICIÓN DEL ENTRENAMIENTO
#===============================================================================
control_train <- trainControl(method = "boot", number = repeticiones_boot,
                              seeds = seeds, returnResamp = "final",
                              verboseIter = TRUE, allowParallel = TRUE)

# AJUSTE DEL MODELO
# ==============================================================================
set.seed(86)
nnet_pvalue_50 <- caret::train(
  form = Class ~ .,
  data = x_train[c("Class", filtrado_anova_pvalue_50)],
  method = "nnet",
  tuneGrid = hiperparametros,
  metric = "Accuracy",
  trControl = control_train,
  # Rango de inicialización de los pesos
  rang = c(-0.7, 0.7),
  # Número máximo de pesos
  MaxNWts = 10000,
  # Para que no se muestre cada iteración por pantalla
  trace = FALSE
)

saveRDS(object = nnet_pvalue_50, file = "nnet_pvalue_50.rds")
registerDoMC(cores = 1)

nnet_pvalue_50

# REPRESENTACIÓN GRÁFICA
# ==============================================================================
ggplot(nnet_pvalue_50, highlight = TRUE) +
  labs(title = "Evolución del accuracy del modelo NNET") +
  theme_bw()

# PREDICCIONES SOBRE TEST 
# ==============================================================================
predic_nnet_pvalue_50 <- predict(object = nnet_pvalue_50, newdata = x_test)
caret::confusionMatrix(predic_nnet_pvalue_50, y_test)



#===============================================================================
#=========================== GBM  ============================================== 
#===============================================================================

# PARALELIZACIÓN DE PROCESO
#===============================================================================

registerDoMC(cores = 11)

# HIPERPARÁMETROS, NÚMERO DE REPETICIONES Y SEMILLAS PARA CADA REPETICIÓN
#===============================================================================
repeticiones_boot <- 50

# Hiperparámetros
hiperparametros <- expand.grid(
  shrinkage = c(0.01, 0.1, 0.3),
  interaction.depth = c(1, 3, 5),
  n.minobsinnode = c(5, 10, 15),
  n.trees = c(100,250,500,1000)
  # bag.fraction = c(0.65, 0.8, 1), 
  # optimal_trees = 0,               
  # min_RMSE = 0,                    
  # min_cor = 0
)

set.seed(86)
seeds <- vector(mode = "list", length = repeticiones_boot + 1)
for (i in 1:repeticiones_boot) {
  seeds[[i]] <- sample.int(1000, nrow(hiperparametros))
}
seeds[[repeticiones_boot + 1]] <- sample.int(1000, 1)

# DEFINICIÓN DEL ENTRENAMIENTO
#===============================================================================
control_train <- trainControl(method = "boot", number = repeticiones_boot,
                              seeds = seeds, returnResamp = "final",
                              verboseIter = TRUE, allowParallel = TRUE)

# AJUSTE DEL MODELO
# ==============================================================================
set.seed(86)
gbm_pvalue_50 <- caret::train(
  form = Class ~ .,
  data = x_train[c("Class", filtrado_anova_pvalue_50)],
  method = "gbm",
  tuneGrid = hiperparametros,
  metric = "Accuracy",
  trControl = control_train
)

saveRDS(object = gbm_pvalue_50, file = "gbm_pvalue_50.rds")
registerDoMC(cores = 1)

gbm_pvalue_50

# REPRESENTACIÓN GRÁFICA
# ==============================================================================
ggplot(gbm_pvalue_50, highlight = TRUE) +
  labs(title = "Evolución del accuracy del modelo GBM") +
  theme_bw()

# PREDICCIONES SOBRE TEST 
# ==============================================================================
predic_gbm_pvalue_50 <- predict(object = gbm_pvalue_50, newdata = x_test)
caret::confusionMatrix(predic_gbm_pvalue_50, y_test)



#===============================================================================
#=========================== XGBM  ============================================= 
#===============================================================================

# PARALELIZACIÓN DE PROCESO
#===============================================================================

registerDoMC(cores = 11)

# HIPERPARÁMETROS, NÚMERO DE REPETICIONES Y SEMILLAS PARA CADA REPETICIÓN
#===============================================================================
repeticiones_boot <- 50

hiperparametros <- expand.grid(
  nrounds = c(500, 1000),
  eta = c(0.01, 0.001), 
  max_depth = c(2, 4, 6),
  gamma = 1,
  colsample_bytree = c(0.2, 0.4),
  min_child_weight = c(1, 5),
  subsample = 1
)



set.seed(86)
seeds <- vector(mode = "list", length = repeticiones_boot + 1)
for (i in 1:repeticiones_boot) {
  seeds[[i]] <- sample.int(1000, nrow(hiperparametros))
}
seeds[[repeticiones_boot + 1]] <- sample.int(1000, 1)

# DEFINICIÓN DEL ENTRENAMIENTO
#===============================================================================
control_train <- trainControl(method = "boot", number = repeticiones_boot,
                              seeds = seeds, returnResamp = "final",
                              verboseIter = TRUE, allowParallel = TRUE)

# AJUSTE DEL MODELO
# ==============================================================================
set.seed(86)
xgbm_pvalue_50 <- caret::train(
  form = Class ~ .,
  data = x_train[c("Class", filtrado_anova_pvalue_50)],
  method = "xgbTree",
  tuneGrid = hiperparametros,
  metric = "Accuracy",
  preProcess = c("center", "scale"),
  trControl = control_train
)


saveRDS(object = xgbm_pvalue_50, file = "xgbm_pvalue_50.rds")
registerDoMC(cores = 1)

xgbm_pvalue_50

# REPRESENTACIÓN GRÁFICA
# ==============================================================================
ggplot(xgbm_pvalue_50, highlight = TRUE) +
  labs(title = "Evolución del accuracy del modelo XGBM") +
  theme_bw()

# PREDICCIONES SOBRE TEST 
# ==============================================================================
predic_xgbm_pvalue_50 <- predict(object = xgbm_pvalue_50, newdata = x_test)
caret::confusionMatrix(predic_xgbm_pvalue_50, y_test)


#===============================================================================
#=========================== glmnet  ===========================================
#===============================================================================

# PARALELIZACIÓN DE PROCESO
#===============================================================================

registerDoMC(cores = 11)

# HIPERPARÁMETROS, NÚMERO DE REPETICIONES Y SEMILLAS PARA CADA REPETICIÓN
#===============================================================================
repeticiones_boot <- 50

hiperparametros <- expand.grid(
  lambda = c(0, 1, 10, 100),
  alpha = c (0.1, 0.01, 0.001)
)



set.seed(86)
seeds <- vector(mode = "list", length = repeticiones_boot + 1)
for (i in 1:repeticiones_boot) {
  seeds[[i]] <- sample.int(1000, nrow(hiperparametros))
}
seeds[[repeticiones_boot + 1]] <- sample.int(1000, 1)

# DEFINICIÓN DEL ENTRENAMIENTO
#===============================================================================
control_train <- trainControl(method = "boot", number = repeticiones_boot,
                              seeds = seeds, returnResamp = "final",
                              verboseIter = TRUE, allowParallel = TRUE)

# AJUSTE DEL MODELO
# ==============================================================================
set.seed(86)
glmnet_pvalue_50 <- caret::train(
  form = Class ~ .,
  data = x_train[c("Class", filtrado_anova_pvalue_50)],
  method = "glmnet",
  family = "binomial",
  tuneGrid = hiperparametros,
  metric = "Accuracy",
  trControl = control_train
)


saveRDS(object = glmnet_pvalue_50, file = "glmnet_pvalue_50.rds")
registerDoMC(cores = 1)

glmnet_pvalue_50

# REPRESENTACIÓN GRÁFICA
# ==============================================================================
ggplot(glmnet_pvalue_50, highlight = TRUE) +
  labs(title = "Evolución del accuracy del modelo glmnet") +
  theme_bw()

# PREDICCIONES SOBRE TEST 
# ==============================================================================
predic_glmnet_pvalue_50 <- predict(object = glmnet_pvalue_50, newdata = x_test)
caret::confusionMatrix(predic_glmnet_pvalue_50, y_test)


#===============================================================================
#=========================== hdda  =============================================
#===============================================================================

# PARALELIZACIÓN DE PROCESO
#===============================================================================

registerDoMC(cores = 11)

# HIPERPARÁMETROS, NÚMERO DE REPETICIONES Y SEMILLAS PARA CADA REPETICIÓN
#===============================================================================
repeticiones_boot <- 50

hiperparametros <- expand.grid(
  threshold = seq(0.1,0.9,by=0.1),
  model = "ALL"
)



set.seed(86)
seeds <- vector(mode = "list", length = repeticiones_boot + 1)
for (i in 1:repeticiones_boot) {
  seeds[[i]] <- sample.int(1000, nrow(hiperparametros))
}
seeds[[repeticiones_boot + 1]] <- sample.int(1000, 1)

# DEFINICIÓN DEL ENTRENAMIENTO
#===============================================================================
control_train <- trainControl(method = "boot", number = repeticiones_boot,
                              seeds = seeds, returnResamp = "final",
                              verboseIter = TRUE, allowParallel = TRUE)

# AJUSTE DEL MODELO
# ==============================================================================
set.seed(86)
hdda_pvalue_50 <- caret::train(
  form = Class ~ .,
  data = x_train[c("Class", filtrado_anova_pvalue_50)],
  method = "hdda",
  tuneGrid = hiperparametros,
  metric = "Accuracy",
  trControl = control_train
)


saveRDS(object = hdda_pvalue_50, file = "hdda_pvalue_50.rds")
registerDoMC(cores = 1)

hdda_pvalue_50

# REPRESENTACIÓN GRÁFICA
# ==============================================================================
ggplot(hdda_pvalue_50, highlight = TRUE) +
  labs(title = "Evolución del accuracy del modelo hdda") +
  theme_bw()

# PREDICCIONES SOBRE TEST 
# ==============================================================================
predic_hdda_pvalue_50 <- predict(object = hdda_pvalue_50, newdata = x_test)
caret::confusionMatrix(predic_hdda_pvalue_50, y_test)


#===============================================================================
#=========================== LogitBoost  =======================================
#===============================================================================

# PARALELIZACIÓN DE PROCESO
#===============================================================================

registerDoMC(cores = 11)

# HIPERPARÁMETROS, NÚMERO DE REPETICIONES Y SEMILLAS PARA CADA REPETICIÓN
#===============================================================================
repeticiones_boot <- 50

hiperparametros <- expand.grid(
  nIter = c(25, 50, 100,250)
)



set.seed(86)
seeds <- vector(mode = "list", length = repeticiones_boot + 1)
for (i in 1:repeticiones_boot) {
  seeds[[i]] <- sample.int(1000, nrow(hiperparametros))
}
seeds[[repeticiones_boot + 1]] <- sample.int(1000, 1)

# DEFINICIÓN DEL ENTRENAMIENTO
#===============================================================================
control_train <- trainControl(method = "boot", number = repeticiones_boot,
                              seeds = seeds, returnResamp = "final",
                              verboseIter = TRUE, allowParallel = TRUE)

# AJUSTE DEL MODELO
# ==============================================================================
set.seed(86)
LogitBoost_pvalue_50 <- caret::train(
  form = Class ~ .,
  data = x_train[c("Class", filtrado_anova_pvalue_50)],
  method = "LogitBoost",
  tuneGrid = hiperparametros,
  metric = "Accuracy",
  trControl = control_train
)


saveRDS(object = LogitBoost_pvalue_50, file = "LogitBoost_pvalue_50.rds")
registerDoMC(cores = 1)

LogitBoost_pvalue_50

# REPRESENTACIÓN GRÁFICA
# ==============================================================================
ggplot(LogitBoost_pvalue_50, highlight = TRUE) +
  labs(title = "Evolución del accuracy del modelo LogitBoost") +
  theme_bw()

# PREDICCIONES SOBRE TEST 
# ==============================================================================
predic_LogitBoost_pvalue_50 <- predict(object = LogitBoost_pvalue_50, newdata = x_test)
caret::confusionMatrix(predic_LogitBoost_pvalue_50, y_test)


#===============================================================================
#=========================== Stacking Ensemble  ================================
#===============================================================================

# DEFINICIÓN DEL ENTRENAMIENTO
#===============================================================================


trainControl <- trainControl(method = "repeatedcv", number = 5, repeats = 3,
                             savePredictions = "all", classProbs = TRUE,
                             index = createFolds(stacking_df$Class, 5),verboseIter = TRUE)


# Run multiplealgorithms in one call.
algorithm <- c("nnet","rf","svmRadial","gbm","xgbTree","hdda","LogitBoost")


stacking_df <- cbind(x_train[,filtrado_anova_pvalue_50],x_train$Class)
names(stacking_df) [51] <- "Class"
levels(stacking_df$Class)
levels(stacking_df$Class) <- make.names(levels(stacking_df$Class))
# typeof(caretList)

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
#*************************** ANOVA PVALUE 500 feautures ************************
#*******************************************************************************


#===============================================================================
#=========================== RANDOM FOREST ===================================== 
#===============================================================================


# PARALELIZACIÓN DE PROCESO
#===============================================================================

#install.packages("doMC", repos="http://R-Forge.R-project.org")

registerDoMC(cores = 11)

# HIPERPARÁMETROS, NÚMERO DE REPETICIONES Y SEMILLAS PARA CADA REPETICIÓN
#===============================================================================
repeticiones_boot <- 50

# Hiperparámetros
hiperparametros <- expand.grid(mtry = c(2, 5, 10, 50),
                               min.node.size = c(2, 3, 4, 5, 10),
                               splitrule = "gini")

set.seed(86)
seeds <- vector(mode = "list", length = repeticiones_boot + 1)
for (i in 1:repeticiones_boot) {
  seeds[[i]] <- sample.int(1000, nrow(hiperparametros))
}
seeds[[repeticiones_boot + 1]] <- sample.int(1000, 1)

# DEFINICIÓN DEL ENTRENAMIENTO
#===============================================================================
control_train <- trainControl(method = "boot", number = repeticiones_boot,
                              seeds = seeds, returnResamp = "final",
                              verboseIter = TRUE, allowParallel = TRUE)

# AJUSTE DEL MODELO
# ==============================================================================
x_train$Class <- as.factor(x_train$Class)

set.seed(86)
rf_pvalue_500 <- caret::train(
  form = Class ~ .,
  data = x_train[c("Class", filtrado_anova_pvalue_500)],
  method = "ranger",
  tuneGrid = hiperparametros,
  metric = "Accuracy",
  trControl = control_train,
  # Número de árboles ajustados
  num.trees = 500)

saveRDS(object = rf_pvalue_500, file = "rf_pvalue_500.rds")
registerDoMC(cores = 1)


rf_pvalue_500

# REPRESENTACIÓN GRÁFICA
# ==============================================================================
ggplot(rf_pvalue_500, highlight = TRUE) +
  labs(title = "Evolución del accuracy del modelo Random Forest") +
  guides(color = guide_legend(title = "mtry"),
         shape = guide_legend(title = "mtry")) +
  theme_bw()

# PREDICCIONES SOBRE TEST 
# ==============================================================================
predic_rd_pvalue_500 <- predict(object = rf_pvalue_500, newdata = x_test)
predic_rd_pvalue_500
y_test <- as.factor(y_test)
caret::confusionMatrix(predic_rd_pvalue_500, y_test)



#===============================================================================
#=========================== SVM Kernel Radial ================================= 
#===============================================================================


# PARALELIZACIÓN DE PROCESO
#===============================================================================

registerDoMC(cores = 11)

# HIPERPARÁMETROS, NÚMERO DE REPETICIONES Y SEMILLAS PARA CADA REPETICIÓN
#===============================================================================
repeticiones_boot <- 50

# Hiperparámetros
hiperparametros <- expand.grid(sigma = c(0.0001, 0.001, 0.01),
                               C = c(1, 10, 50, 100, 250, 500, 700, 1000))
set.seed(123)
seeds <- vector(mode = "list", length = repeticiones_boot + 1)
for (i in 1:repeticiones_boot) {
  seeds[[i]] <- sample.int(1000, nrow(hiperparametros))
}
seeds[[repeticiones_boot + 1]] <- sample.int(1000, 1)

# DEFINICIÓN DEL ENTRENAMIENTO
#===============================================================================
control_train <- trainControl(method = "boot", number = repeticiones_boot,
                              seeds = seeds, returnResamp = "final",
                              verboseIter = TRUE, allowParallel = TRUE)
# AJUSTE DEL MODELO
# ==============================================================================
set.seed(86)
svmrad_pvalue_500 <- caret::train(
  form = Class ~ .,
  data = x_train[c("Class", filtrado_anova_pvalue_500)],
  method = "svmRadial",
  tuneGrid = hiperparametros,
  metric = "Accuracy",
  trControl = control_train
)
registerDoMC(cores = 1)
saveRDS(object = svmrad_pvalue_500, file = "svmrad_pvalue_500.rds")

svmrad_pvalue_500

# REPRESENTACIÓN GRÁFICA
# ==============================================================================
ggplot(svmrad_pvalue_500, highlight = TRUE) +
  labs(title = "Evolución del accuracy del modelo SVM Radial") +
  theme_bw()

# PREDICCIONES SOBRE TEST 
# ==============================================================================
predic_SVM_pvalue_500 <- predict(object = svmrad_pvalue_500, newdata = x_test)
caret::confusionMatrix(predic_SVM_pvalue_500, y_test)



#===============================================================================
#=========================== Neural Network  =================================== 
#===============================================================================

# PARALELIZACIÓN DE PROCESO
#===============================================================================

registerDoMC(cores = 11)

# HIPERPARÁMETROS, NÚMERO DE REPETICIONES Y SEMILLAS PARA CADA REPETICIÓN
#===============================================================================
repeticiones_boot <- 50

# Hiperparámetros
hiperparametros <- expand.grid(size = c(5, 10, 15, 20, 40),
                               decay = c(0.01, 0.1))

set.seed(86)
seeds <- vector(mode = "list", length = repeticiones_boot + 1)
for (i in 1:repeticiones_boot) {
  seeds[[i]] <- sample.int(1000, nrow(hiperparametros))
}
seeds[[repeticiones_boot + 1]] <- sample.int(1000, 1)

# DEFINICIÓN DEL ENTRENAMIENTO
#===============================================================================
control_train <- trainControl(method = "boot", number = repeticiones_boot,
                              seeds = seeds, returnResamp = "final",
                              verboseIter = TRUE, allowParallel = TRUE)

# AJUSTE DEL MODELO
# ==============================================================================
set.seed(86)
nnet_pvalue_500 <- caret::train(
  form = Class ~ .,
  data = x_train[c("Class", filtrado_anova_pvalue_500)],
  method = "nnet",
  tuneGrid = hiperparametros,
  metric = "Accuracy",
  trControl = control_train,
  # Rango de inicialización de los pesos
  rang = c(-0.7, 0.7),
  # Número máximo de pesos
  MaxNWts = 10000,
  # Para que no se muestre cada iteración por pantalla
  trace = FALSE
)

saveRDS(object = nnet_pvalue_500, file = "nnet_pvalue_100.rds")
registerDoMC(cores = 1)

nnet_pvalue_500

# REPRESENTACIÓN GRÁFICA
# ==============================================================================
ggplot(nnet_pvalue_500, highlight = TRUE) +
  labs(title = "Evolución del accuracy del modelo NNET") +
  theme_bw()

# PREDICCIONES SOBRE TEST 
# ==============================================================================
predic_nnet_pvalue_500 <- predict(object = nnet_pvalue_500, newdata = x_test)
caret::confusionMatrix(predic_nnet_pvalue_500, y_test)



#===============================================================================
#=========================== GBM  ============================================== 
#===============================================================================

# PARALELIZACIÓN DE PROCESO
#===============================================================================

registerDoMC(cores = 11)

# HIPERPARÁMETROS, NÚMERO DE REPETICIONES Y SEMILLAS PARA CADA REPETICIÓN
#===============================================================================
repeticiones_boot <- 50

# Hiperparámetros
hiperparametros <- expand.grid(
  shrinkage = c(0.01, 0.1, 0.3),
  interaction.depth = c(1, 3, 5),
  n.minobsinnode = c(5, 10, 15),
  n.trees = c(100,250,500,1000)
  # bag.fraction = c(0.65, 0.8, 1), 
  # optimal_trees = 0,               
  # min_RMSE = 0,                    
  # min_cor = 0
)

set.seed(86)
seeds <- vector(mode = "list", length = repeticiones_boot + 1)
for (i in 1:repeticiones_boot) {
  seeds[[i]] <- sample.int(1000, nrow(hiperparametros))
}
seeds[[repeticiones_boot + 1]] <- sample.int(1000, 1)

# DEFINICIÓN DEL ENTRENAMIENTO
#===============================================================================
control_train <- trainControl(method = "boot", number = repeticiones_boot,
                              seeds = seeds, returnResamp = "final",
                              verboseIter = TRUE, allowParallel = TRUE)

# AJUSTE DEL MODELO
# ==============================================================================
set.seed(86)
gbm_pvalue_500 <- caret::train(
  form = Class ~ .,
  data = x_train[c("Class", filtrado_anova_pvalue_500)],
  method = "gbm",
  tuneGrid = hiperparametros,
  metric = "Accuracy",
  trControl = control_train
)

saveRDS(object = gbm_pvalue_500, file = "gbm_pvalue_500.rds")
registerDoMC(cores = 1)

gbm_pvalue_500

# REPRESENTACIÓN GRÁFICA
# ==============================================================================
ggplot(gbm_pvalue_500, highlight = TRUE) +
  labs(title = "Evolución del accuracy del modelo GBM") +
  theme_bw()

# PREDICCIONES SOBRE TEST 
# ==============================================================================
predic_gbm_pvalue_500 <- predict(object = gbm_pvalue_500, newdata = x_test)
caret::confusionMatrix(predic_gbm_pvalue_500, y_test)



#===============================================================================
#=========================== XGBM  ============================================= 
#===============================================================================

# PARALELIZACIÓN DE PROCESO
#===============================================================================

registerDoMC(cores = 11)

# HIPERPARÁMETROS, NÚMERO DE REPETICIONES Y SEMILLAS PARA CADA REPETICIÓN
#===============================================================================
repeticiones_boot <- 50

hiperparametros <- expand.grid(
  nrounds = c(500, 1000),
  eta = c(0.01, 0.001), 
  max_depth = c(2, 4, 6),
  gamma = 1,
  colsample_bytree = c(0.2, 0.4),
  min_child_weight = c(1, 5),
  subsample = 1
)



set.seed(86)
seeds <- vector(mode = "list", length = repeticiones_boot + 1)
for (i in 1:repeticiones_boot) {
  seeds[[i]] <- sample.int(1000, nrow(hiperparametros))
}
seeds[[repeticiones_boot + 1]] <- sample.int(1000, 1)

# DEFINICIÓN DEL ENTRENAMIENTO
#===============================================================================
control_train <- trainControl(method = "boot", number = repeticiones_boot,
                              seeds = seeds, returnResamp = "final",
                              verboseIter = TRUE, allowParallel = TRUE)

# AJUSTE DEL MODELO
# ==============================================================================
set.seed(86)
xgbm_pvalue_500 <- caret::train(
  form = Class ~ .,
  data = x_train[c("Class", filtrado_anova_pvalue_500)],
  method = "xgbTree",
  tuneGrid = hiperparametros,
  metric = "Accuracy",
  preProcess = c("center", "scale"),
  trControl = control_train
)


saveRDS(object = xgbm_pvalue_500, file = "xgbm_pvalue_500.rds")
registerDoMC(cores = 1)

xgbm_pvalue_100

# REPRESENTACIÓN GRÁFICA
# ==============================================================================
ggplot(xgbm_pvalue_100, highlight = TRUE) +
  labs(title = "Evolución del accuracy del modelo XGBM") +
  theme_bw()

# PREDICCIONES SOBRE TEST 
# ==============================================================================
predic_xgbm_pvalue_100 <- predict(object = xgbm_pvalue_100, newdata = x_test)
caret::confusionMatrix(predic_xgbm_pvalue_100, y_test)


#===============================================================================
#=========================== glmnet  ===========================================
#===============================================================================

# PARALELIZACIÓN DE PROCESO
#===============================================================================

registerDoMC(cores = 11)

# HIPERPARÁMETROS, NÚMERO DE REPETICIONES Y SEMILLAS PARA CADA REPETICIÓN
#===============================================================================
repeticiones_boot <- 50

hiperparametros <- expand.grid(
  lambda = c(0, 1, 10, 100),
  alpha = c (0.1, 0.01, 0.001)
)



set.seed(86)
seeds <- vector(mode = "list", length = repeticiones_boot + 1)
for (i in 1:repeticiones_boot) {
  seeds[[i]] <- sample.int(1000, nrow(hiperparametros))
}
seeds[[repeticiones_boot + 1]] <- sample.int(1000, 1)

# DEFINICIÓN DEL ENTRENAMIENTO
#===============================================================================
control_train <- trainControl(method = "boot", number = repeticiones_boot,
                              seeds = seeds, returnResamp = "final",
                              verboseIter = TRUE, allowParallel = TRUE)

# AJUSTE DEL MODELO
# ==============================================================================
set.seed(86)
glmnet_pvalue_500 <- caret::train(
  form = Class ~ .,
  data = x_train[c("Class", filtrado_anova_pvalue_500)],
  method = "glmnet",
  family = "binomial",
  tuneGrid = hiperparametros,
  metric = "Accuracy",
  trControl = control_train
)


saveRDS(object = glmnet_pvalue_500, file = "glmnet_pvalue_500.rds")
registerDoMC(cores = 1)

glmnet_pvalue_500

# REPRESENTACIÓN GRÁFICA
# ==============================================================================
ggplot(glmnet_pvalue_500, highlight = TRUE) +
  labs(title = "Evolución del accuracy del modelo glmnet") +
  theme_bw()

# PREDICCIONES SOBRE TEST 
# ==============================================================================
predic_glmnet_pvalue_500 <- predict(object = glmnet_pvalue_500, newdata = x_test)
caret::confusionMatrix(predic_glmnet_pvalue_500, y_test)


#===============================================================================
#=========================== hdda  =============================================
#===============================================================================

# PARALELIZACIÓN DE PROCESO
#===============================================================================

registerDoMC(cores = 11)

# HIPERPARÁMETROS, NÚMERO DE REPETICIONES Y SEMILLAS PARA CADA REPETICIÓN
#===============================================================================
repeticiones_boot <- 50

hiperparametros <- expand.grid(
  threshold = seq(0.1,0.9,by=0.1),
  model = "ALL"
)



set.seed(86)
seeds <- vector(mode = "list", length = repeticiones_boot + 1)
for (i in 1:repeticiones_boot) {
  seeds[[i]] <- sample.int(1000, nrow(hiperparametros))
}
seeds[[repeticiones_boot + 1]] <- sample.int(1000, 1)

# DEFINICIÓN DEL ENTRENAMIENTO
#===============================================================================
control_train <- trainControl(method = "boot", number = repeticiones_boot,
                              seeds = seeds, returnResamp = "final",
                              verboseIter = TRUE, allowParallel = TRUE)

# AJUSTE DEL MODELO
# ==============================================================================
set.seed(86)
hdda_pvalue_100 <- caret::train(
  form = Class ~ .,
  data = x_train[c("Class", filtrado_anova_pvalue_100)],
  method = "hdda",
  tuneGrid = hiperparametros,
  metric = "Accuracy",
  trControl = control_train
)


saveRDS(object = hdda_pvalue_100, file = "hdda_pvalue_100.rds")
registerDoMC(cores = 1)

hdda_pvalue_100

# REPRESENTACIÓN GRÁFICA
# ==============================================================================
ggplot(hdda_pvalue_100, highlight = TRUE) +
  labs(title = "Evolución del accuracy del modelo hdda") +
  theme_bw()

# PREDICCIONES SOBRE TEST 
# ==============================================================================
predic_hdda_pvalue_100 <- predict(object = hdda_pvalue_100, newdata = x_test)
caret::confusionMatrix(predic_hdda_pvalue_100, y_test)


#===============================================================================
#=========================== LogitBoost  =======================================
#===============================================================================

# PARALELIZACIÓN DE PROCESO
#===============================================================================

registerDoMC(cores = 11)

# HIPERPARÁMETROS, NÚMERO DE REPETICIONES Y SEMILLAS PARA CADA REPETICIÓN
#===============================================================================
repeticiones_boot <- 50

hiperparametros <- expand.grid(
  nIter = c(25, 50, 100,250)
)



set.seed(86)
seeds <- vector(mode = "list", length = repeticiones_boot + 1)
for (i in 1:repeticiones_boot) {
  seeds[[i]] <- sample.int(1000, nrow(hiperparametros))
}
seeds[[repeticiones_boot + 1]] <- sample.int(1000, 1)

# DEFINICIÓN DEL ENTRENAMIENTO
#===============================================================================
control_train <- trainControl(method = "boot", number = repeticiones_boot,
                              seeds = seeds, returnResamp = "final",
                              verboseIter = TRUE, allowParallel = TRUE)

# AJUSTE DEL MODELO
# ==============================================================================
set.seed(86)
LogitBoost_pvalue_100 <- caret::train(
  form = Class ~ .,
  data = x_train[c("Class", filtrado_anova_pvalue_100)],
  method = "LogitBoost",
  tuneGrid = hiperparametros,
  metric = "Accuracy",
  trControl = control_train
)


saveRDS(object = LogitBoost_pvalue_100, file = "LogitBoost_pvalue_100.rds")
registerDoMC(cores = 1)

LogitBoost_pvalue_100

# REPRESENTACIÓN GRÁFICA
# ==============================================================================
ggplot(LogitBoost_pvalue_100, highlight = TRUE) +
  labs(title = "Evolución del accuracy del modelo LogitBoost") +
  theme_bw()

# PREDICCIONES SOBRE TEST 
# ==============================================================================
predic_LogitBoost_pvalue_100 <- predict(object = LogitBoost_pvalue_100, newdata = x_test)
caret::confusionMatrix(predic_LogitBoost_pvalue_100, y_test)


#===============================================================================
#=========================== Stacking Ensemble  ================================
#===============================================================================

# DEFINICIÓN DEL ENTRENAMIENTO
#===============================================================================


trainControl <- trainControl(method = "repeatedcv", number = 5, repeats = 3,
                             savePredictions = "all", classProbs = TRUE,
                             index = createFolds(stacking_df$Class, 5),verboseIter = TRUE)


# Run multiplealgorithms in one call.
algorithm <- c("nnet","rf","svmRadial","gbm","xgbTree","hdda","LogitBoost")


stacking_df <- cbind(x_train[,filtrado_anova_pvalue_100],x_train$Class)
names(stacking_df) [101] <- "Class"
levels(stacking_df$Class)
levels(stacking_df$Class) <- make.names(levels(stacking_df$Class))

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



=======
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

########################################################################
####################### Cervical Data ##################################
########################################################################


filepath <- system.file("extdata/cervical.txt", package = "MLSeq")
cervical <- STAD.table(filepath, header=TRUE)
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
# show(gbm_fit)
# trained(gbm_fit)
# plot(gbm_fit)

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


########################################################################
####################### TCGA Data ######################################
########################################################################

GDCprojects = getGDCprojects()
names_projects <- GDCprojects[c("project_id", "name")]

# Stomach Adenocarcinoma
TCGAbiolinks:::getProjectSummary("TCGA-STAD")

query_TCGA = GDCquery(
  project = "TCGA-STAD",
  data.category = "Transcriptome Profiling", # parameter enforced by GDCquery
  experimental.strategy = "RNA-Seq",
  workflow.type = "HTSeq - Counts")

STAD_res = getResults(query_TCGA) # make results as table
head(STAD_res) # data of the first 6 patients.
colnames(STAD_res) # columns present in the table

head(STAD_res$sample_type) # first 6 types of tissue.
table(STAD_res$sample_type) # summary of distinct tissues types present in this study

GDCdownload(query = query_TCGA) # download the data into my current workdirectory 

tcga_data = GDCprepare(query_TCGA) # load the data 
dim(tcga_data)

# colData () allows us to access the clinical data associated with our samples. 
# The colnames () and rownames () functions can be used to extract the column and row names of a given table, respectively.

colnames(colData(tcga_data))

# tables 

table(tcga_data@colData$gender)
table(tcga_data@colData$vital_status)
table(tcga_data@colData$tissue_or_organ_of_origin)
table(tcga_data@colData$race)
addmargins(round(prop.table(table(tcga_data@colData$vital_status,tcga_data@colData$gender))*100,3))
addmargins(round(prop.table(table(tcga_data@colData$vital_status,tcga_data@colData$race))*100,3))
addmargins(round(prop.table(table(tcga_data@colData$tissue_or_organ_of_origin,tcga_data@colData$gender))*100,3))
addmargins(round(prop.table(table(tcga_data@colData$tissue_or_organ_of_origin,tcga_data@colData$race))*100,3))
addmargins(round(prop.table(table(tcga_data@colData$vital_status,tcga_data@colData$tissue_or_organ_of_origin))*100,3))

# assay function 
head(assay(tcga_data)[,1:10]) # expression of first 6 genes and first 10 samples
head(rowData(tcga_data))     # ensembl id and gene id of the first 6 genes.

# Save the data as a file, if you need it later, you can just load this file
# instead of having to run the whole pipeline again

# saveRDS(object = tcga_data,
#         file = "tcga_data_STAD.RDS",
#         compress = FALSE)

# load data with the next sentence
tcga_data = STADRDS(file = "tcga_data_STAD.RDS")


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

# x_train$Class <- NULL
# x_test$Class <- NULL

prop.table(table(x_train$Class)) %>% round(3)
prop.table(table(x_test$Class)) %>% round(3)

summary(x_train$ENSG00000000003)
head(x_train[,1:6])

# ## Genes with variance close to zero
# 
# # We proceed to eliminate those genes whose maximum expression 
# # does not exceed 5 times the minimum expression (max/min) < 5 
# # and whose absolute difference between maximum and minimum does not exceed 500 units (max − min < 500)
# 
# 
# filter_variance <- function(x){
#   # This function returns TRUE for all non-numeric columns and, in the case
#   # of the numerical ones, those that exceed the minimum variance conditions.
#   if(is.numeric(x)){
#     maximum <- max(x)
#     minimum <- min(x)
#     ratio <- maximum / minimum
#     range <- maximum - minimum
#     return(ratio >= 5 & range >= 500)
#   }
#   
#   else{
#     return(TRUE)
#   }
# }
# 
# # The columns that meet the condition are identified
# genes_variance <- map_lgl(.x = x_train, .f = filter_variance)
# 
# # Number of columns (genes) excluded
# sum(genes_variance == FALSE)

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
# VERSIÓN PARALELIZADA DE BOOTSTRAPPING PARA FILTRADO POR ANOVA
#===============================================================================

n_cores <- parallel::detectCores() - 1
registerDoParallel(makeCluster(n_cores))
getDoParWorkers()


# Número de iteraciones bootstrapping
n_boot <- 100

# Semillas para que los muestreos sean reproducibles
set.seed(86)
seeds = sample.int(1000, size = n_boot)


resultados_anova_pvalue <- foreach(i = 1:n_boot) %dopar% {
  library(dplyr)
  library(purrr)
  # Se crea una muestra por bootstrapping
  set.seed(seeds[i])
  indices <- sample(1:nrow(x_train), size = nrow(x_train), replace = TRUE)
  pseudo_muestra <- x_train[indices, ]
  
  # Se calculan los p-values para la nueva muestra 
  p_values <- dplyr::select(pseudo_muestra,-Class) %>%
    purrr::map_dbl(.f = custom_anova, y = pseudo_muestra$Class) 
  
  # Se devuelven los p-value
  p_values
}

options(cores = 1)


# Los resultados almacenados en forma de lista se convierten en dataframe
names(resultados_anova_pvalue) <-  paste("resample", 1:n_boot, sep = "_") 
resultados_anova_pvalue <- data.frame(resultados_anova_pvalue)
resultados_anova_pvalue <- resultados_anova_pvalue %>% rownames_to_column(var = "gen")
resultados_anova_pvalue <- resultados_anova_pvalue %>%
  mutate(pvalue_medio = rowMeans(resultados_anova_pvalue[, -1])) %>%
  arrange(pvalue_medio)

# Se guarda en disco el objeto cSTADo para no tener que repetir de nuevo toda la
# computación.
saveRDS(object = resultados_anova_pvalue, file = "resultados_anova_pvalue.rds")

resultados_anova_pvalue %>% dplyr::select(1,2,3,4) %>% head()

# Se filtran los 100, 50 y 500 genes identificados como más relevantes mediante anova
filtrado_anova_pvalue_500 <-  resultados_anova_pvalue %>% pull(gen) %>% head(500)
filtrado_anova_pvalue_100 <- resultados_anova_pvalue %>% pull(gen) %>% head(100)
filtrado_anova_pvalue_50  <- resultados_anova_pvalue %>% pull(gen) %>% head(50)



#*******************************************************************************
#*************************** ANOVA PVALUE 100 feautures ************************
#*******************************************************************************


#===============================================================================
#=========================== RANDOM FOREST ===================================== 
#===============================================================================


# PARALELIZACIÓN DE PROCESO
#===============================================================================

#install.packages("doMC", repos="http://R-Forge.R-project.org")

registerDoMC(cores = 11)

# HIPERPARÁMETROS, NÚMERO DE REPETICIONES Y SEMILLAS PARA CADA REPETICIÓN
#===============================================================================
repeticiones_boot <- 50

# Hiperparámetros
hiperparametros <- expand.grid(mtry = c(2, 5, 10, 50),
                               min.node.size = c(2, 3, 4, 5, 10),
                               splitrule = "gini")

set.seed(86)
seeds <- vector(mode = "list", length = repeticiones_boot + 1)
for (i in 1:repeticiones_boot) {
  seeds[[i]] <- sample.int(1000, nrow(hiperparametros))
}
seeds[[repeticiones_boot + 1]] <- sample.int(1000, 1)

# DEFINICIÓN DEL ENTRENAMIENTO
#===============================================================================
control_train <- trainControl(method = "boot", number = repeticiones_boot,
                              seeds = seeds, returnResamp = "final",
                              verboseIter = TRUE, allowParallel = TRUE)

# AJUSTE DEL MODELO
# ==============================================================================
x_train$Class <- as.factor(x_train$Class)

set.seed(86)
rf_pvalue_100 <- caret::train(
  form = Class ~ .,
  data = x_train[c("Class", filtrado_anova_pvalue_100)],
  method = "ranger",
  tuneGrid = hiperparametros,
  metric = "Accuracy",
  trControl = control_train,
  # Número de árboles ajustados
  num.trees = 500)

saveRDS(object = rf_pvalue_100, file = "rf_pvalue_100.rds")
registerDoMC(cores = 1)


rf_pvalue_100

# REPRESENTACIÓN GRÁFICA
# ==============================================================================
ggplot(rf_pvalue_100, highlight = TRUE) +
  labs(title = "Evolución del accuracy del modelo Random Forest") +
  guides(color = guide_legend(title = "mtry"),
         shape = guide_legend(title = "mtry")) +
  theme_bw()

# PREDICCIONES SOBRE TEST 
# ==============================================================================
predic_rd_pvalue_100 <- predict(object = rf_pvalue_100, newdata = x_test)
predic_rd_pvalue_100
y_test <- as.factor(y_test)
caret::confusionMatrix(predic_rd_pvalue_100, y_test)



#===============================================================================
#=========================== SVM Kernel Radial ================================= 
#===============================================================================


# PARALELIZACIÓN DE PROCESO
#===============================================================================

registerDoMC(cores = 11)

# HIPERPARÁMETROS, NÚMERO DE REPETICIONES Y SEMILLAS PARA CADA REPETICIÓN
#===============================================================================
repeticiones_boot <- 50

# Hiperparámetros
hiperparametros <- expand.grid(sigma = c(0.0001, 0.001, 0.01),
                               C = c(1, 10, 50, 100, 250, 500, 700, 1000))
set.seed(123)
seeds <- vector(mode = "list", length = repeticiones_boot + 1)
for (i in 1:repeticiones_boot) {
  seeds[[i]] <- sample.int(1000, nrow(hiperparametros))
}
seeds[[repeticiones_boot + 1]] <- sample.int(1000, 1)

# DEFINICIÓN DEL ENTRENAMIENTO
#===============================================================================
control_train <- trainControl(method = "boot", number = repeticiones_boot,
                              seeds = seeds, returnResamp = "final",
                              verboseIter = TRUE, allowParallel = TRUE)
# AJUSTE DEL MODELO
# ==============================================================================
set.seed(86)
svmrad_pvalue_100 <- caret::train(
  form = Class ~ .,
  data = x_train[c("Class", filtrado_anova_pvalue_100)],
  method = "svmRadial",
  tuneGrid = hiperparametros,
  metric = "Accuracy",
  trControl = control_train
)
registerDoMC(cores = 1)
saveRDS(object = svmrad_pvalue_100, file = "svmrad_pvalue_100.rds")

svmrad_pvalue_100

# REPRESENTACIÓN GRÁFICA
# ==============================================================================
ggplot(svmrad_pvalue_100, highlight = TRUE) +
  labs(title = "Evolución del accuracy del modelo SVM Radial") +
  theme_bw()

# PREDICCIONES SOBRE TEST 
# ==============================================================================
predic_SVM_pvalue_100 <- predict(object = svmrad_pvalue_100, newdata = x_test)
caret::confusionMatrix(predic_SVM_pvalue_100, y_test)



#===============================================================================
#=========================== Neural Network  =================================== 
#===============================================================================

# PARALELIZACIÓN DE PROCESO
#===============================================================================

registerDoMC(cores = 11)

# HIPERPARÁMETROS, NÚMERO DE REPETICIONES Y SEMILLAS PARA CADA REPETICIÓN
#===============================================================================
repeticiones_boot <- 50

# Hiperparámetros
hiperparametros <- expand.grid(size = c(5, 10, 15, 20, 40),
                               decay = c(0.01, 0.1))

set.seed(86)
seeds <- vector(mode = "list", length = repeticiones_boot + 1)
for (i in 1:repeticiones_boot) {
  seeds[[i]] <- sample.int(1000, nrow(hiperparametros))
}
seeds[[repeticiones_boot + 1]] <- sample.int(1000, 1)

# DEFINICIÓN DEL ENTRENAMIENTO
#===============================================================================
control_train <- trainControl(method = "boot", number = repeticiones_boot,
                              seeds = seeds, returnResamp = "final",
                              verboseIter = TRUE, allowParallel = TRUE)

# AJUSTE DEL MODELO
# ==============================================================================
set.seed(86)
nnet_pvalue_100 <- caret::train(
  form = Class ~ .,
  data = x_train[c("Class", filtrado_anova_pvalue_100)],
  method = "nnet",
  tuneGrid = hiperparametros,
  metric = "Accuracy",
  trControl = control_train,
  # Rango de inicialización de los pesos
  rang = c(-0.7, 0.7),
  # Número máximo de pesos
  MaxNWts = 10000,
  # Para que no se muestre cada iteración por pantalla
  trace = FALSE
)

saveRDS(object = nnet_pvalue_100, file = "nnet_pvalue_100.rds")
registerDoMC(cores = 1)

nnet_pvalue_100

# REPRESENTACIÓN GRÁFICA
# ==============================================================================
ggplot(nnet_pvalue_100, highlight = TRUE) +
  labs(title = "Evolución del accuracy del modelo NNET") +
  theme_bw()

# PREDICCIONES SOBRE TEST 
# ==============================================================================
predic_nnet_pvalue_100 <- predict(object = nnet_pvalue_100, newdata = x_test)
caret::confusionMatrix(predic_nnet_pvalue_100, y_test)



#===============================================================================
#=========================== GBM  ============================================== 
#===============================================================================

# PARALELIZACIÓN DE PROCESO
#===============================================================================

registerDoMC(cores = 11)

# HIPERPARÁMETROS, NÚMERO DE REPETICIONES Y SEMILLAS PARA CADA REPETICIÓN
#===============================================================================
repeticiones_boot <- 50

# Hiperparámetros
hiperparametros <- expand.grid(
  shrinkage = c(0.01, 0.1, 0.3),
  interaction.depth = c(1, 3, 5),
  n.minobsinnode = c(5, 10, 15),
  n.trees = c(100,250,500,1000)
  # bag.fraction = c(0.65, 0.8, 1), 
  # optimal_trees = 0,               
  # min_RMSE = 0,                    
  # min_cor = 0
)

set.seed(86)
seeds <- vector(mode = "list", length = repeticiones_boot + 1)
for (i in 1:repeticiones_boot) {
  seeds[[i]] <- sample.int(1000, nrow(hiperparametros))
}
seeds[[repeticiones_boot + 1]] <- sample.int(1000, 1)

# DEFINICIÓN DEL ENTRENAMIENTO
#===============================================================================
control_train <- trainControl(method = "boot", number = repeticiones_boot,
                              seeds = seeds, returnResamp = "final",
                              verboseIter = TRUE, allowParallel = TRUE)

# AJUSTE DEL MODELO
# ==============================================================================
set.seed(86)
gbm_pvalue_100 <- caret::train(
  form = Class ~ .,
  data = x_train[c("Class", filtrado_anova_pvalue_100)],
  method = "gbm",
  tuneGrid = hiperparametros,
  metric = "Accuracy",
  trControl = control_train
)

saveRDS(object = gbm_pvalue_100, file = "gbm_pvalue_100.rds")
registerDoMC(cores = 1)

gbm_pvalue_100

# REPRESENTACIÓN GRÁFICA
# ==============================================================================
ggplot(gbm_pvalue_100, highlight = TRUE) +
  labs(title = "Evolución del accuracy del modelo GBM") +
  theme_bw()

# PREDICCIONES SOBRE TEST 
# ==============================================================================
predic_gbm_pvalue_100 <- predict(object = gbm_pvalue_100, newdata = x_test)
caret::confusionMatrix(predic_gbm_pvalue_100, y_test)



#===============================================================================
#=========================== XGBM  ============================================= 
#===============================================================================

# PARALELIZACIÓN DE PROCESO
#===============================================================================

registerDoMC(cores = 11)

# HIPERPARÁMETROS, NÚMERO DE REPETICIONES Y SEMILLAS PARA CADA REPETICIÓN
#===============================================================================
repeticiones_boot <- 50

hiperparametros <- expand.grid(
  nrounds = c(500, 1000),
  eta = c(0.01, 0.001), 
  max_depth = c(2, 4, 6),
  gamma = 1,
  colsample_bytree = c(0.2, 0.4),
  min_child_weight = c(1, 5),
  subsample = 1
)



set.seed(86)
seeds <- vector(mode = "list", length = repeticiones_boot + 1)
for (i in 1:repeticiones_boot) {
  seeds[[i]] <- sample.int(1000, nrow(hiperparametros))
}
seeds[[repeticiones_boot + 1]] <- sample.int(1000, 1)

# DEFINICIÓN DEL ENTRENAMIENTO
#===============================================================================
control_train <- trainControl(method = "boot", number = repeticiones_boot,
                              seeds = seeds, returnResamp = "final",
                              verboseIter = TRUE, allowParallel = TRUE)

# AJUSTE DEL MODELO
# ==============================================================================
set.seed(86)
xgbm_pvalue_100 <- caret::train(
  form = Class ~ .,
  data = x_train[c("Class", filtrado_anova_pvalue_100)],
  method = "xgbTree",
  tuneGrid = hiperparametros,
  metric = "Accuracy",
  preProcess = c("center", "scale"),
  trControl = control_train
)


saveRDS(object = xgbm_pvalue_100, file = "xgbm_pvalue_100.rds")
registerDoMC(cores = 1)

xgbm_pvalue_100

# REPRESENTACIÓN GRÁFICA
# ==============================================================================
ggplot(xgbm_pvalue_100, highlight = TRUE) +
  labs(title = "Evolución del accuracy del modelo XGBM") +
  theme_bw()

# PREDICCIONES SOBRE TEST 
# ==============================================================================
predic_xgbm_pvalue_100 <- predict(object = xgbm_pvalue_100, newdata = x_test)
caret::confusionMatrix(predic_xgbm_pvalue_100, y_test)


#===============================================================================
#=========================== glmnet  ===========================================
#===============================================================================

# PARALELIZACIÓN DE PROCESO
#===============================================================================

registerDoMC(cores = 11)

# HIPERPARÁMETROS, NÚMERO DE REPETICIONES Y SEMILLAS PARA CADA REPETICIÓN
#===============================================================================
repeticiones_boot <- 50

hiperparametros <- expand.grid(
  lambda = c(0, 1, 10, 100),
  alpha = c (0.1, 0.01, 0.001)
)



set.seed(86)
seeds <- vector(mode = "list", length = repeticiones_boot + 1)
for (i in 1:repeticiones_boot) {
  seeds[[i]] <- sample.int(1000, nrow(hiperparametros))
}
seeds[[repeticiones_boot + 1]] <- sample.int(1000, 1)

# DEFINICIÓN DEL ENTRENAMIENTO
#===============================================================================
control_train <- trainControl(method = "boot", number = repeticiones_boot,
                              seeds = seeds, returnResamp = "final",
                              verboseIter = TRUE, allowParallel = TRUE)

# AJUSTE DEL MODELO
# ==============================================================================
set.seed(86)
glmnet_pvalue_100 <- caret::train(
  form = Class ~ .,
  data = x_train[c("Class", filtrado_anova_pvalue_100)],
  method = "glmnet",
  family = "binomial",
  tuneGrid = hiperparametros,
  metric = "Accuracy",
  trControl = control_train
)


saveRDS(object = glmnet_pvalue_100, file = "glmnet_pvalue_100.rds")
registerDoMC(cores = 1)

glmnet_pvalue_100

# REPRESENTACIÓN GRÁFICA
# ==============================================================================
ggplot(glmnet_pvalue_100, highlight = TRUE) +
  labs(title = "Evolución del accuracy del modelo glmnet") +
  theme_bw()

# PREDICCIONES SOBRE TEST 
# ==============================================================================
predic_glmnet_pvalue_100 <- predict(object = glmnet_pvalue_100, newdata = x_test)
caret::confusionMatrix(predic_glmnet_pvalue_100, y_test)


#===============================================================================
#=========================== hdda  =============================================
#===============================================================================

# PARALELIZACIÓN DE PROCESO
#===============================================================================

registerDoMC(cores = 11)

# HIPERPARÁMETROS, NÚMERO DE REPETICIONES Y SEMILLAS PARA CADA REPETICIÓN
#===============================================================================
repeticiones_boot <- 50

hiperparametros <- expand.grid(
  threshold = seq(0.1,0.9,by=0.1),
  model = "ALL"
)



set.seed(86)
seeds <- vector(mode = "list", length = repeticiones_boot + 1)
for (i in 1:repeticiones_boot) {
  seeds[[i]] <- sample.int(1000, nrow(hiperparametros))
}
seeds[[repeticiones_boot + 1]] <- sample.int(1000, 1)

# DEFINICIÓN DEL ENTRENAMIENTO
#===============================================================================
control_train <- trainControl(method = "boot", number = repeticiones_boot,
                              seeds = seeds, returnResamp = "final",
                              verboseIter = TRUE, allowParallel = TRUE)

# AJUSTE DEL MODELO
# ==============================================================================
set.seed(86)
hdda_pvalue_100 <- caret::train(
  form = Class ~ .,
  data = x_train[c("Class", filtrado_anova_pvalue_100)],
  method = "hdda",
  tuneGrid = hiperparametros,
  metric = "Accuracy",
  trControl = control_train
)


saveRDS(object = hdda_pvalue_100, file = "hdda_pvalue_100.rds")
registerDoMC(cores = 1)

hdda_pvalue_100

# REPRESENTACIÓN GRÁFICA
# ==============================================================================
ggplot(hdda_pvalue_100, highlight = TRUE) +
  labs(title = "Evolución del accuracy del modelo hdda") +
  theme_bw()

# PREDICCIONES SOBRE TEST 
# ==============================================================================
predic_hdda_pvalue_100 <- predict(object = hdda_pvalue_100, newdata = x_test)
caret::confusionMatrix(predic_hdda_pvalue_100, y_test)


#===============================================================================
#=========================== LogitBoost  =======================================
#===============================================================================

# PARALELIZACIÓN DE PROCESO
#===============================================================================

registerDoMC(cores = 11)

# HIPERPARÁMETROS, NÚMERO DE REPETICIONES Y SEMILLAS PARA CADA REPETICIÓN
#===============================================================================
repeticiones_boot <- 50

hiperparametros <- expand.grid(
  nIter = c(25, 50, 100,250)
)



set.seed(86)
seeds <- vector(mode = "list", length = repeticiones_boot + 1)
for (i in 1:repeticiones_boot) {
  seeds[[i]] <- sample.int(1000, nrow(hiperparametros))
}
seeds[[repeticiones_boot + 1]] <- sample.int(1000, 1)

# DEFINICIÓN DEL ENTRENAMIENTO
#===============================================================================
control_train <- trainControl(method = "boot", number = repeticiones_boot,
                              seeds = seeds, returnResamp = "final",
                              verboseIter = TRUE, allowParallel = TRUE)

# AJUSTE DEL MODELO
# ==============================================================================
set.seed(86)
LogitBoost_pvalue_100 <- caret::train(
  form = Class ~ .,
  data = x_train[c("Class", filtrado_anova_pvalue_100)],
  method = "LogitBoost",
  tuneGrid = hiperparametros,
  metric = "Accuracy",
  trControl = control_train
)


saveRDS(object = LogitBoost_pvalue_100, file = "LogitBoost_pvalue_100.rds")
registerDoMC(cores = 1)

LogitBoost_pvalue_100

# REPRESENTACIÓN GRÁFICA
# ==============================================================================
ggplot(LogitBoost_pvalue_100, highlight = TRUE) +
  labs(title = "Evolución del accuracy del modelo LogitBoost") +
  theme_bw()

# PREDICCIONES SOBRE TEST 
# ==============================================================================
predic_LogitBoost_pvalue_100 <- predict(object = LogitBoost_pvalue_100, newdata = x_test)
caret::confusionMatrix(predic_LogitBoost_pvalue_100, y_test)


#===============================================================================
#=========================== Stacking Ensemble  ================================
#===============================================================================

# DEFINICIÓN DEL ENTRENAMIENTO
#===============================================================================


trainControl <- trainControl(method = "repeatedcv", number = 5, repeats = 3,
                             savePredictions = "all", classProbs = TRUE,
                             index = createFolds(stacking_df$Class, 5),verboseIter = TRUE)


# Run multiplealgorithms in one call.
algorithm <- c("nnet","rf","svmRadial","gbm","xgbTree","hdda","LogitBoost")


stacking_df <- cbind(x_train[,filtrado_anova_pvalue_100],x_train$Class)
names(stacking_df) [101] <- "Class"
levels(stacking_df$Class)
levels(stacking_df$Class) <- make.names(levels(stacking_df$Class))

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
#*************************** ANOVA PVALUE 50 feautures ************************
#*******************************************************************************


#===============================================================================
#=========================== RANDOM FOREST ===================================== 
#===============================================================================


# PARALELIZACIÓN DE PROCESO
#===============================================================================

#install.packages("doMC", repos="http://R-Forge.R-project.org")

registerDoMC(cores = 11)

# HIPERPARÁMETROS, NÚMERO DE REPETICIONES Y SEMILLAS PARA CADA REPETICIÓN
#===============================================================================
repeticiones_boot <- 50

# Hiperparámetros
hiperparametros <- expand.grid(mtry = c(2, 5, 10, 50),
                               min.node.size = c(2, 3, 4, 5, 10),
                               splitrule = "gini")

set.seed(86)
seeds <- vector(mode = "list", length = repeticiones_boot + 1)
for (i in 1:repeticiones_boot) {
  seeds[[i]] <- sample.int(1000, nrow(hiperparametros))
}
seeds[[repeticiones_boot + 1]] <- sample.int(1000, 1)

# DEFINICIÓN DEL ENTRENAMIENTO
#===============================================================================
control_train <- trainControl(method = "boot", number = repeticiones_boot,
                              seeds = seeds, returnResamp = "final",
                              verboseIter = TRUE, allowParallel = TRUE)

# AJUSTE DEL MODELO
# ==============================================================================
x_train$Class <- as.factor(x_train$Class)

set.seed(86)
rf_pvalue_50 <- caret::train(
  form = Class ~ .,
  data = x_train[c("Class", filtrado_anova_pvalue_50)],
  method = "ranger",
  tuneGrid = hiperparametros,
  metric = "Accuracy",
  trControl = control_train,
  # Número de árboles ajustados
  num.trees = 500)

saveRDS(object = rf_pvalue_50, file = "rf_pvalue_50.rds")
registerDoMC(cores = 1)


rf_pvalue_50

# REPRESENTACIÓN GRÁFICA
# ==============================================================================
ggplot(rf_pvalue_50, highlight = TRUE) +
  labs(title = "Evolución del accuracy del modelo Random Forest") +
  guides(color = guide_legend(title = "mtry"),
         shape = guide_legend(title = "mtry")) +
  theme_bw()

# PREDICCIONES SOBRE TEST 
# ==============================================================================
predic_rd_pvalue_50 <- predict(object = rf_pvalue_50, newdata = x_test)
predic_rd_pvalue_50
y_test <- as.factor(y_test)
caret::confusionMatrix(predic_rd_pvalue_50, y_test)



#===============================================================================
#=========================== SVM Kernel Radial ================================= 
#===============================================================================


# PARALELIZACIÓN DE PROCESO
#===============================================================================

registerDoMC(cores = 11)

# HIPERPARÁMETROS, NÚMERO DE REPETICIONES Y SEMILLAS PARA CADA REPETICIÓN
#===============================================================================
repeticiones_boot <- 50

# Hiperparámetros
hiperparametros <- expand.grid(sigma = c(0.0001, 0.001, 0.01),
                               C = c(1, 10, 50, 100, 250, 500, 700, 1000))
set.seed(123)
seeds <- vector(mode = "list", length = repeticiones_boot + 1)
for (i in 1:repeticiones_boot) {
  seeds[[i]] <- sample.int(1000, nrow(hiperparametros))
}
seeds[[repeticiones_boot + 1]] <- sample.int(1000, 1)

# DEFINICIÓN DEL ENTRENAMIENTO
#===============================================================================
control_train <- trainControl(method = "boot", number = repeticiones_boot,
                              seeds = seeds, returnResamp = "final",
                              verboseIter = TRUE, allowParallel = TRUE)
# AJUSTE DEL MODELO
# ==============================================================================
set.seed(86)
svmrad_pvalue_50 <- caret::train(
  form = Class ~ .,
  data = x_train[c("Class", filtrado_anova_pvalue_50)],
  method = "svmRadial",
  tuneGrid = hiperparametros,
  metric = "Accuracy",
  trControl = control_train
)
registerDoMC(cores = 1)
saveRDS(object = svmrad_pvalue_50, file = "svmrad_pvalue_50.rds")

svmrad_pvalue_50

# REPRESENTACIÓN GRÁFICA
# ==============================================================================
ggplot(svmrad_pvalue_50, highlight = TRUE) +
  labs(title = "Evolución del accuracy del modelo SVM Radial") +
  theme_bw()

# PREDICCIONES SOBRE TEST 
# ==============================================================================
predic_SVM_pvalue_50 <- predict(object = svmrad_pvalue_50, newdata = x_test)
caret::confusionMatrix(predic_SVM_pvalue_50, y_test)



#===============================================================================
#=========================== Neural Network  =================================== 
#===============================================================================

# PARALELIZACIÓN DE PROCESO
#===============================================================================

registerDoMC(cores = 11)

# HIPERPARÁMETROS, NÚMERO DE REPETICIONES Y SEMILLAS PARA CADA REPETICIÓN
#===============================================================================
repeticiones_boot <- 50

# Hiperparámetros
hiperparametros <- expand.grid(size = c(5, 10, 15, 20, 40),
                               decay = c(0.01, 0.1))

set.seed(86)
seeds <- vector(mode = "list", length = repeticiones_boot + 1)
for (i in 1:repeticiones_boot) {
  seeds[[i]] <- sample.int(1000, nrow(hiperparametros))
}
seeds[[repeticiones_boot + 1]] <- sample.int(1000, 1)

# DEFINICIÓN DEL ENTRENAMIENTO
#===============================================================================
control_train <- trainControl(method = "boot", number = repeticiones_boot,
                              seeds = seeds, returnResamp = "final",
                              verboseIter = TRUE, allowParallel = TRUE)

# AJUSTE DEL MODELO
# ==============================================================================
set.seed(86)
nnet_pvalue_50 <- caret::train(
  form = Class ~ .,
  data = x_train[c("Class", filtrado_anova_pvalue_50)],
  method = "nnet",
  tuneGrid = hiperparametros,
  metric = "Accuracy",
  trControl = control_train,
  # Rango de inicialización de los pesos
  rang = c(-0.7, 0.7),
  # Número máximo de pesos
  MaxNWts = 10000,
  # Para que no se muestre cada iteración por pantalla
  trace = FALSE
)

saveRDS(object = nnet_pvalue_50, file = "nnet_pvalue_50.rds")
registerDoMC(cores = 1)

nnet_pvalue_50

# REPRESENTACIÓN GRÁFICA
# ==============================================================================
ggplot(nnet_pvalue_50, highlight = TRUE) +
  labs(title = "Evolución del accuracy del modelo NNET") +
  theme_bw()

# PREDICCIONES SOBRE TEST 
# ==============================================================================
predic_nnet_pvalue_50 <- predict(object = nnet_pvalue_50, newdata = x_test)
caret::confusionMatrix(predic_nnet_pvalue_50, y_test)



#===============================================================================
#=========================== GBM  ============================================== 
#===============================================================================

# PARALELIZACIÓN DE PROCESO
#===============================================================================

registerDoMC(cores = 11)

# HIPERPARÁMETROS, NÚMERO DE REPETICIONES Y SEMILLAS PARA CADA REPETICIÓN
#===============================================================================
repeticiones_boot <- 50

# Hiperparámetros
hiperparametros <- expand.grid(
  shrinkage = c(0.01, 0.1, 0.3),
  interaction.depth = c(1, 3, 5),
  n.minobsinnode = c(5, 10, 15),
  n.trees = c(100,250,500,1000)
  # bag.fraction = c(0.65, 0.8, 1), 
  # optimal_trees = 0,               
  # min_RMSE = 0,                    
  # min_cor = 0
)

set.seed(86)
seeds <- vector(mode = "list", length = repeticiones_boot + 1)
for (i in 1:repeticiones_boot) {
  seeds[[i]] <- sample.int(1000, nrow(hiperparametros))
}
seeds[[repeticiones_boot + 1]] <- sample.int(1000, 1)

# DEFINICIÓN DEL ENTRENAMIENTO
#===============================================================================
control_train <- trainControl(method = "boot", number = repeticiones_boot,
                              seeds = seeds, returnResamp = "final",
                              verboseIter = TRUE, allowParallel = TRUE)

# AJUSTE DEL MODELO
# ==============================================================================
set.seed(86)
gbm_pvalue_50 <- caret::train(
  form = Class ~ .,
  data = x_train[c("Class", filtrado_anova_pvalue_50)],
  method = "gbm",
  tuneGrid = hiperparametros,
  metric = "Accuracy",
  trControl = control_train
)

saveRDS(object = gbm_pvalue_50, file = "gbm_pvalue_50.rds")
registerDoMC(cores = 1)

gbm_pvalue_50

# REPRESENTACIÓN GRÁFICA
# ==============================================================================
ggplot(gbm_pvalue_50, highlight = TRUE) +
  labs(title = "Evolución del accuracy del modelo GBM") +
  theme_bw()

# PREDICCIONES SOBRE TEST 
# ==============================================================================
predic_gbm_pvalue_50 <- predict(object = gbm_pvalue_50, newdata = x_test)
caret::confusionMatrix(predic_gbm_pvalue_50, y_test)



#===============================================================================
#=========================== XGBM  ============================================= 
#===============================================================================

# PARALELIZACIÓN DE PROCESO
#===============================================================================

registerDoMC(cores = 11)

# HIPERPARÁMETROS, NÚMERO DE REPETICIONES Y SEMILLAS PARA CADA REPETICIÓN
#===============================================================================
repeticiones_boot <- 50

hiperparametros <- expand.grid(
  nrounds = c(500, 1000),
  eta = c(0.01, 0.001), 
  max_depth = c(2, 4, 6),
  gamma = 1,
  colsample_bytree = c(0.2, 0.4),
  min_child_weight = c(1, 5),
  subsample = 1
)



set.seed(86)
seeds <- vector(mode = "list", length = repeticiones_boot + 1)
for (i in 1:repeticiones_boot) {
  seeds[[i]] <- sample.int(1000, nrow(hiperparametros))
}
seeds[[repeticiones_boot + 1]] <- sample.int(1000, 1)

# DEFINICIÓN DEL ENTRENAMIENTO
#===============================================================================
control_train <- trainControl(method = "boot", number = repeticiones_boot,
                              seeds = seeds, returnResamp = "final",
                              verboseIter = TRUE, allowParallel = TRUE)

# AJUSTE DEL MODELO
# ==============================================================================
set.seed(86)
xgbm_pvalue_50 <- caret::train(
  form = Class ~ .,
  data = x_train[c("Class", filtrado_anova_pvalue_50)],
  method = "xgbTree",
  tuneGrid = hiperparametros,
  metric = "Accuracy",
  preProcess = c("center", "scale"),
  trControl = control_train
)


saveRDS(object = xgbm_pvalue_50, file = "xgbm_pvalue_50.rds")
registerDoMC(cores = 1)

xgbm_pvalue_50

# REPRESENTACIÓN GRÁFICA
# ==============================================================================
ggplot(xgbm_pvalue_50, highlight = TRUE) +
  labs(title = "Evolución del accuracy del modelo XGBM") +
  theme_bw()

# PREDICCIONES SOBRE TEST 
# ==============================================================================
predic_xgbm_pvalue_50 <- predict(object = xgbm_pvalue_50, newdata = x_test)
caret::confusionMatrix(predic_xgbm_pvalue_50, y_test)


#===============================================================================
#=========================== glmnet  ===========================================
#===============================================================================

# PARALELIZACIÓN DE PROCESO
#===============================================================================

registerDoMC(cores = 11)

# HIPERPARÁMETROS, NÚMERO DE REPETICIONES Y SEMILLAS PARA CADA REPETICIÓN
#===============================================================================
repeticiones_boot <- 50

hiperparametros <- expand.grid(
  lambda = c(0, 1, 10, 100),
  alpha = c (0.1, 0.01, 0.001)
)



set.seed(86)
seeds <- vector(mode = "list", length = repeticiones_boot + 1)
for (i in 1:repeticiones_boot) {
  seeds[[i]] <- sample.int(1000, nrow(hiperparametros))
}
seeds[[repeticiones_boot + 1]] <- sample.int(1000, 1)

# DEFINICIÓN DEL ENTRENAMIENTO
#===============================================================================
control_train <- trainControl(method = "boot", number = repeticiones_boot,
                              seeds = seeds, returnResamp = "final",
                              verboseIter = TRUE, allowParallel = TRUE)

# AJUSTE DEL MODELO
# ==============================================================================
set.seed(86)
glmnet_pvalue_50 <- caret::train(
  form = Class ~ .,
  data = x_train[c("Class", filtrado_anova_pvalue_50)],
  method = "glmnet",
  family = "binomial",
  tuneGrid = hiperparametros,
  metric = "Accuracy",
  trControl = control_train
)


saveRDS(object = glmnet_pvalue_50, file = "glmnet_pvalue_50.rds")
registerDoMC(cores = 1)

glmnet_pvalue_50

# REPRESENTACIÓN GRÁFICA
# ==============================================================================
ggplot(glmnet_pvalue_50, highlight = TRUE) +
  labs(title = "Evolución del accuracy del modelo glmnet") +
  theme_bw()

# PREDICCIONES SOBRE TEST 
# ==============================================================================
predic_glmnet_pvalue_50 <- predict(object = glmnet_pvalue_50, newdata = x_test)
caret::confusionMatrix(predic_glmnet_pvalue_50, y_test)


#===============================================================================
#=========================== hdda  =============================================
#===============================================================================

# PARALELIZACIÓN DE PROCESO
#===============================================================================

registerDoMC(cores = 11)

# HIPERPARÁMETROS, NÚMERO DE REPETICIONES Y SEMILLAS PARA CADA REPETICIÓN
#===============================================================================
repeticiones_boot <- 50

hiperparametros <- expand.grid(
  threshold = seq(0.1,0.9,by=0.1),
  model = "ALL"
)



set.seed(86)
seeds <- vector(mode = "list", length = repeticiones_boot + 1)
for (i in 1:repeticiones_boot) {
  seeds[[i]] <- sample.int(1000, nrow(hiperparametros))
}
seeds[[repeticiones_boot + 1]] <- sample.int(1000, 1)

# DEFINICIÓN DEL ENTRENAMIENTO
#===============================================================================
control_train <- trainControl(method = "boot", number = repeticiones_boot,
                              seeds = seeds, returnResamp = "final",
                              verboseIter = TRUE, allowParallel = TRUE)

# AJUSTE DEL MODELO
# ==============================================================================
set.seed(86)
hdda_pvalue_50 <- caret::train(
  form = Class ~ .,
  data = x_train[c("Class", filtrado_anova_pvalue_50)],
  method = "hdda",
  tuneGrid = hiperparametros,
  metric = "Accuracy",
  trControl = control_train
)


saveRDS(object = hdda_pvalue_50, file = "hdda_pvalue_100.rds")
registerDoMC(cores = 1)

hdda_pvalue_50

# REPRESENTACIÓN GRÁFICA
# ==============================================================================
ggplot(hdda_pvalue_50, highlight = TRUE) +
  labs(title = "Evolución del accuracy del modelo hdda") +
  theme_bw()

# PREDICCIONES SOBRE TEST 
# ==============================================================================
predic_hdda_pvalue_50 <- predict(object = hdda_pvalue_50, newdata = x_test)
caret::confusionMatrix(predic_hdda_pvalue_50, y_test)


#===============================================================================
#=========================== LogitBoost  =======================================
#===============================================================================

# PARALELIZACIÓN DE PROCESO
#===============================================================================

registerDoMC(cores = 11)

# HIPERPARÁMETROS, NÚMERO DE REPETICIONES Y SEMILLAS PARA CADA REPETICIÓN
#===============================================================================
repeticiones_boot <- 50

hiperparametros <- expand.grid(
  nIter = c(25, 50, 100,250)
)



set.seed(86)
seeds <- vector(mode = "list", length = repeticiones_boot + 1)
for (i in 1:repeticiones_boot) {
  seeds[[i]] <- sample.int(1000, nrow(hiperparametros))
}
seeds[[repeticiones_boot + 1]] <- sample.int(1000, 1)

# DEFINICIÓN DEL ENTRENAMIENTO
#===============================================================================
control_train <- trainControl(method = "boot", number = repeticiones_boot,
                              seeds = seeds, returnResamp = "final",
                              verboseIter = TRUE, allowParallel = TRUE)

# AJUSTE DEL MODELO
# ==============================================================================
set.seed(86)
LogitBoost_pvalue_50 <- caret::train(
  form = Class ~ .,
  data = x_train[c("Class", filtrado_anova_pvalue_50)],
  method = "LogitBoost",
  tuneGrid = hiperparametros,
  metric = "Accuracy",
  trControl = control_train
)


saveRDS(object = LogitBoost_pvalue_50, file = "LogitBoost_pvalue_50.rds")
registerDoMC(cores = 1)

LogitBoost_pvalue_50

# REPRESENTACIÓN GRÁFICA
# ==============================================================================
ggplot(LogitBoost_pvalue_50, highlight = TRUE) +
  labs(title = "Evolución del accuracy del modelo LogitBoost") +
  theme_bw()

# PREDICCIONES SOBRE TEST 
# ==============================================================================
predic_LogitBoost_pvalue_50 <- predict(object = LogitBoost_pvalue_50, newdata = x_test)
caret::confusionMatrix(predic_LogitBoost_pvalue_50, y_test)


#===============================================================================
#=========================== Stacking Ensemble  ================================
#===============================================================================

# DEFINICIÓN DEL ENTRENAMIENTO
#===============================================================================


trainControl <- trainControl(method = "repeatedcv", number = 5, repeats = 3,
                             savePredictions = "all", classProbs = TRUE,
                             index = createFolds(stacking_df$Class, 5),verboseIter = TRUE)


# Run multiplealgorithms in one call.
algorithm <- c("nnet","rf","svmRadial","gbm","xgbTree","hdda","LogitBoost")


stacking_df <- cbind(x_train[,filtrado_anova_pvalue_50],x_train$Class)
names(stacking_df) [51] <- "Class"
levels(stacking_df$Class)
levels(stacking_df$Class) <- make.names(levels(stacking_df$Class))
# typeof(caretList)

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
#*************************** ANOVA PVALUE 500 features ************************
#*******************************************************************************


#===============================================================================
#=========================== RANDOM FOREST ===================================== 
#===============================================================================


# PARALELIZACIÓN DE PROCESO
#===============================================================================

#install.packages("doMC", repos="http://R-Forge.R-project.org")

registerDoMC(cores = 11)

# HIPERPARÁMETROS, NÚMERO DE REPETICIONES Y SEMILLAS PARA CADA REPETICIÓN
#===============================================================================
repeticiones_boot <- 50

# Hiperparámetros
hiperparametros <- expand.grid(mtry = c(2, 5, 10, 50),
                               min.node.size = c(2, 3, 4, 5, 10),
                               splitrule = "gini")

set.seed(86)
seeds <- vector(mode = "list", length = repeticiones_boot + 1)
for (i in 1:repeticiones_boot) {
  seeds[[i]] <- sample.int(1000, nrow(hiperparametros))
}
seeds[[repeticiones_boot + 1]] <- sample.int(1000, 1)

# DEFINICIÓN DEL ENTRENAMIENTO
#===============================================================================
control_train <- trainControl(method = "boot", number = repeticiones_boot,
                              seeds = seeds, returnResamp = "final",
                              verboseIter = TRUE, allowParallel = TRUE)

# AJUSTE DEL MODELO
# ==============================================================================
x_train$Class <- as.factor(x_train$Class)

set.seed(86)
rf_pvalue_500 <- caret::train(
  form = Class ~ .,
  data = x_train[c("Class", filtrado_anova_pvalue_500)],
  method = "ranger",
  tuneGrid = hiperparametros,
  metric = "Accuracy",
  trControl = control_train,
  # Número de árboles ajustados
  num.trees = 500)

saveRDS(object = rf_pvalue_500, file = "rf_pvalue_500.rds")
registerDoMC(cores = 1)


rf_pvalue_500

# REPRESENTACIÓN GRÁFICA
# ==============================================================================
ggplot(rf_pvalue_500, highlight = TRUE) +
  labs(title = "Evolución del accuracy del modelo Random Forest") +
  guides(color = guide_legend(title = "mtry"),
         shape = guide_legend(title = "mtry")) +
  theme_bw()

# PREDICCIONES SOBRE TEST 
# ==============================================================================
predic_rd_pvalue_500 <- predict(object = rf_pvalue_500, newdata = x_test)
predic_rd_pvalue_500
y_test <- as.factor(y_test)
caret::confusionMatrix(predic_rd_pvalue_500, y_test)



#===============================================================================
#=========================== SVM Kernel Radial ================================= 
#===============================================================================


# PARALELIZACIÓN DE PROCESO
#===============================================================================

registerDoMC(cores = 11)

# HIPERPARÁMETROS, NÚMERO DE REPETICIONES Y SEMILLAS PARA CADA REPETICIÓN
#===============================================================================
repeticiones_boot <- 50

# Hiperparámetros
hiperparametros <- expand.grid(sigma = c(0.0001, 0.001, 0.01),
                               C = c(1, 10, 50, 100, 250, 500, 700, 1000))
set.seed(123)
seeds <- vector(mode = "list", length = repeticiones_boot + 1)
for (i in 1:repeticiones_boot) {
  seeds[[i]] <- sample.int(1000, nrow(hiperparametros))
}
seeds[[repeticiones_boot + 1]] <- sample.int(1000, 1)

# DEFINICIÓN DEL ENTRENAMIENTO
#===============================================================================
control_train <- trainControl(method = "boot", number = repeticiones_boot,
                              seeds = seeds, returnResamp = "final",
                              verboseIter = TRUE, allowParallel = TRUE)
# AJUSTE DEL MODELO
# ==============================================================================
set.seed(86)
svmrad_pvalue_500 <- caret::train(
  form = Class ~ .,
  data = x_train[c("Class", filtrado_anova_pvalue_500)],
  method = "svmRadial",
  tuneGrid = hiperparametros,
  metric = "Accuracy",
  trControl = control_train
)
registerDoMC(cores = 1)
saveRDS(object = svmrad_pvalue_500, file = "svmrad_pvalue_500.rds")

svmrad_pvalue_500

# REPRESENTACIÓN GRÁFICA
# ==============================================================================
ggplot(svmrad_pvalue_500, highlight = TRUE) +
  labs(title = "Evolución del accuracy del modelo SVM Radial") +
  theme_bw()

# PREDICCIONES SOBRE TEST 
# ==============================================================================
predic_SVM_pvalue_500 <- predict(object = svmrad_pvalue_500, newdata = x_test)
caret::confusionMatrix(predic_SVM_pvalue_500, y_test)



#===============================================================================
#=========================== Neural Network  =================================== 
#===============================================================================

# PARALELIZACIÓN DE PROCESO
#===============================================================================

registerDoMC(cores = 11)

# HIPERPARÁMETROS, NÚMERO DE REPETICIONES Y SEMILLAS PARA CADA REPETICIÓN
#===============================================================================
repeticiones_boot <- 50

# Hiperparámetros
hiperparametros <- expand.grid(size = c(5, 10, 15, 20),
                               decay = c(0.01, 0.1))

set.seed(86)
seeds <- vector(mode = "list", length = repeticiones_boot + 1)
for (i in 1:repeticiones_boot) {
  seeds[[i]] <- sample.int(1000, nrow(hiperparametros))
}
seeds[[repeticiones_boot + 1]] <- sample.int(1000, 1)

# DEFINICIÓN DEL ENTRENAMIENTO
#===============================================================================
control_train <- trainControl(method = "boot", number = repeticiones_boot,
                              seeds = seeds, returnResamp = "final",
                              verboseIter = TRUE, allowParallel = TRUE)

# AJUSTE DEL MODELO
# ==============================================================================
set.seed(86)
nnet_pvalue_500 <- caret::train(
  form = Class ~ .,
  data = x_train[c("Class", filtrado_anova_pvalue_500)],
  method = "nnet",
  tuneGrid = hiperparametros,
  metric = "Accuracy",
  trControl = control_train,
  # Rango de inicialización de los pesos
  rang = c(-0.7, 0.7),
  # Número máximo de pesos
  MaxNWts = 21000,
  # Para que se muestre cada iteración por pantalla
  trace = TRUE
)

saveRDS(object = nnet_pvalue_500, file = "nnet_pvalue_500.rds")
registerDoMC(cores = 1)

nnet_pvalue_500

# REPRESENTACIÓN GRÁFICA
# ==============================================================================
ggplot(nnet_pvalue_500, highlight = TRUE) +
  labs(title = "Evolución del accuracy del modelo NNET") +
  theme_bw()

# PREDICCIONES SOBRE TEST 
# ==============================================================================
predic_nnet_pvalue_500 <- predict(object = nnet_pvalue_500, newdata = x_test)
caret::confusionMatrix(predic_nnet_pvalue_500, y_test)



#===============================================================================
#=========================== GBM  ============================================== 
#===============================================================================

# PARALELIZACIÓN DE PROCESO
#===============================================================================

registerDoMC(cores = 11)

# HIPERPARÁMETROS, NÚMERO DE REPETICIONES Y SEMILLAS PARA CADA REPETICIÓN
#===============================================================================
repeticiones_boot <- 50

# Hiperparámetros
hiperparametros <- expand.grid(
  shrinkage = c(0.01, 0.1, 0.3),
  interaction.depth = c(1, 3, 5),
  n.minobsinnode = c(5, 10, 15),
  n.trees = c(100,250,500,1000)
  # bag.fraction = c(0.65, 0.8, 1), 
  # optimal_trees = 0,               
  # min_RMSE = 0,                    
  # min_cor = 0
)

set.seed(86)
seeds <- vector(mode = "list", length = repeticiones_boot + 1)
for (i in 1:repeticiones_boot) {
  seeds[[i]] <- sample.int(1000, nrow(hiperparametros))
}
seeds[[repeticiones_boot + 1]] <- sample.int(1000, 1)

# DEFINICIÓN DEL ENTRENAMIENTO
#===============================================================================
control_train <- trainControl(method = "boot", number = repeticiones_boot,
                              seeds = seeds, returnResamp = "final",
                              verboseIter = TRUE, allowParallel = TRUE)

# AJUSTE DEL MODELO
# ==============================================================================
set.seed(86)
gbm_pvalue_500 <- caret::train(
  form = Class ~ .,
  data = x_train[c("Class", filtrado_anova_pvalue_500)],
  method = "gbm",
  tuneGrid = hiperparametros,
  metric = "Accuracy",
  trControl = control_train
)

saveRDS(object = gbm_pvalue_500, file = "gbm_pvalue_500.rds")
registerDoMC(cores = 1)

gbm_pvalue_500

# REPRESENTACIÓN GRÁFICA
# ==============================================================================
ggplot(gbm_pvalue_500, highlight = TRUE) +
  labs(title = "Evolución del accuracy del modelo GBM") +
  theme_bw()

# PREDICCIONES SOBRE TEST 
# ==============================================================================
predic_gbm_pvalue_500 <- predict(object = gbm_pvalue_500, newdata = x_test)
caret::confusionMatrix(predic_gbm_pvalue_500, y_test)



#===============================================================================
#=========================== XGBM  ============================================= 
#===============================================================================

# PARALELIZACIÓN DE PROCESO
#===============================================================================

registerDoMC(cores = 11)

# HIPERPARÁMETROS, NÚMERO DE REPETICIONES Y SEMILLAS PARA CADA REPETICIÓN
#===============================================================================
repeticiones_boot <- 50

hiperparametros <- expand.grid(
  nrounds = c(500, 1000),
  eta = c(0.01, 0.001), 
  max_depth = c(2, 4, 6),
  gamma = 1,
  colsample_bytree = c(0.2, 0.4),
  min_child_weight = c(1, 5),
  subsample = 1
)



set.seed(86)
seeds <- vector(mode = "list", length = repeticiones_boot + 1)
for (i in 1:repeticiones_boot) {
  seeds[[i]] <- sample.int(1000, nrow(hiperparametros))
}
seeds[[repeticiones_boot + 1]] <- sample.int(1000, 1)

# DEFINICIÓN DEL ENTRENAMIENTO
#===============================================================================
control_train <- trainControl(method = "boot", number = repeticiones_boot,
                              seeds = seeds, returnResamp = "final",
                              verboseIter = TRUE, allowParallel = TRUE)

# AJUSTE DEL MODELO
# ==============================================================================
set.seed(86)
xgbm_pvalue_500 <- caret::train(
  form = Class ~ .,
  data = x_train[c("Class", filtrado_anova_pvalue_500)],
  method = "xgbTree",
  tuneGrid = hiperparametros,
  metric = "Accuracy",
  preProcess = c("center", "scale"),
  trControl = control_train
)


saveRDS(object = xgbm_pvalue_500, file = "xgbm_pvalue_500.rds")
registerDoMC(cores = 1)

xgbm_pvalue_500

# REPRESENTACIÓN GRÁFICA
# ==============================================================================
ggplot(xgbm_pvalue_500, highlight = TRUE) +
  labs(title = "Evolución del accuracy del modelo XGBM") +
  theme_bw()

# PREDICCIONES SOBRE TEST 
# ==============================================================================
predic_xgbm_pvalue_500 <- predict(object = xgbm_pvalue_500, newdata = x_test)
caret::confusionMatrix(predic_xgbm_pvalue_500, y_test)


#===============================================================================
#=========================== glmnet  ===========================================
#===============================================================================

# PARALELIZACIÓN DE PROCESO
#===============================================================================

registerDoMC(cores = 11)

# HIPERPARÁMETROS, NÚMERO DE REPETICIONES Y SEMILLAS PARA CADA REPETICIÓN
#===============================================================================
repeticiones_boot <- 50

hiperparametros <- expand.grid(
  lambda = c(0, 1, 10, 100),
  alpha = c (0.1, 0.01, 0.001)
)



set.seed(86)
seeds <- vector(mode = "list", length = repeticiones_boot + 1)
for (i in 1:repeticiones_boot) {
  seeds[[i]] <- sample.int(1000, nrow(hiperparametros))
}
seeds[[repeticiones_boot + 1]] <- sample.int(1000, 1)

# DEFINICIÓN DEL ENTRENAMIENTO
#===============================================================================
control_train <- trainControl(method = "boot", number = repeticiones_boot,
                              seeds = seeds, returnResamp = "final",
                              verboseIter = TRUE, allowParallel = TRUE)

# AJUSTE DEL MODELO
# ==============================================================================
set.seed(86)
glmnet_pvalue_500 <- caret::train(
  form = Class ~ .,
  data = x_train[c("Class", filtrado_anova_pvalue_500)],
  method = "glmnet",
  family = "binomial",
  tuneGrid = hiperparametros,
  metric = "Accuracy",
  trControl = control_train
)


saveRDS(object = glmnet_pvalue_500, file = "glmnet_pvalue_500.rds")
registerDoMC(cores = 1)

glmnet_pvalue_500

# REPRESENTACIÓN GRÁFICA
# ==============================================================================
ggplot(glmnet_pvalue_500, highlight = TRUE) +
  labs(title = "Evolución del accuracy del modelo glmnet") +
  theme_bw()

# PREDICCIONES SOBRE TEST 
# ==============================================================================
predic_glmnet_pvalue_500 <- predict(object = glmnet_pvalue_500, newdata = x_test)
caret::confusionMatrix(predic_glmnet_pvalue_500, y_test)


#===============================================================================
#=========================== hdda  =============================================
#===============================================================================

# PARALELIZACIÓN DE PROCESO
#===============================================================================

registerDoMC(cores = 11)

# HIPERPARÁMETROS, NÚMERO DE REPETICIONES Y SEMILLAS PARA CADA REPETICIÓN
#===============================================================================
repeticiones_boot <- 50

hiperparametros <- expand.grid(
  threshold = seq(0.1,0.9,by=0.1),
  model = "ALL"
)



set.seed(86)
seeds <- vector(mode = "list", length = repeticiones_boot + 1)
for (i in 1:repeticiones_boot) {
  seeds[[i]] <- sample.int(1000, nrow(hiperparametros))
}
seeds[[repeticiones_boot + 1]] <- sample.int(1000, 1)

# DEFINICIÓN DEL ENTRENAMIENTO
#===============================================================================
control_train <- trainControl(method = "boot", number = repeticiones_boot,
                              seeds = seeds, returnResamp = "final",
                              verboseIter = TRUE, allowParallel = TRUE)

# AJUSTE DEL MODELO
# ==============================================================================
set.seed(86)
hdda_pvalue_500 <- caret::train(
  form = Class ~ .,
  data = x_train[c("Class", filtrado_anova_pvalue_500)],
  method = "hdda",
  tuneGrid = hiperparametros,
  metric = "Accuracy",
  trControl = control_train
)


saveRDS(object = hdda_pvalue_500, file = "hdda_pvalue_500.rds")
registerDoMC(cores = 1)

hdda_pvalue_500

# REPRESENTACIÓN GRÁFICA
# ==============================================================================
ggplot(hdda_pvalue_500, highlight = TRUE) +
  labs(title = "Evolución del accuracy del modelo hdda") +
  theme_bw()

# PREDICCIONES SOBRE TEST 
# ==============================================================================
predic_hdda_pvalue_500 <- predict(object = hdda_pvalue_500, newdata = x_test)
caret::confusionMatrix(predic_hdda_pvalue_500, y_test)


#===============================================================================
#=========================== LogitBoost  =======================================
#===============================================================================

# PARALELIZACIÓN DE PROCESO
#===============================================================================

registerDoMC(cores = 11)

# HIPERPARÁMETROS, NÚMERO DE REPETICIONES Y SEMILLAS PARA CADA REPETICIÓN
#===============================================================================
repeticiones_boot <- 50

hiperparametros <- expand.grid(
  nIter = c(25, 50, 100,250)
)



set.seed(86)
seeds <- vector(mode = "list", length = repeticiones_boot + 1)
for (i in 1:repeticiones_boot) {
  seeds[[i]] <- sample.int(1000, nrow(hiperparametros))
}
seeds[[repeticiones_boot + 1]] <- sample.int(1000, 1)

# DEFINICIÓN DEL ENTRENAMIENTO
#===============================================================================
control_train <- trainControl(method = "boot", number = repeticiones_boot,
                              seeds = seeds, returnResamp = "final",
                              verboseIter = TRUE, allowParallel = TRUE)

# AJUSTE DEL MODELO
# ==============================================================================
set.seed(86)
LogitBoost_pvalue_500 <- caret::train(
  form = Class ~ .,
  data = x_train[c("Class", filtrado_anova_pvalue_500)],
  method = "LogitBoost",
  tuneGrid = hiperparametros,
  metric = "Accuracy",
  trControl = control_train
)


saveRDS(object = LogitBoost_pvalue_500, file = "LogitBoost_pvalue_500.rds")
registerDoMC(cores = 1)

LogitBoost_pvalue_500

# REPRESENTACIÓN GRÁFICA
# ==============================================================================
ggplot(LogitBoost_pvalue_500, highlight = TRUE) +
  labs(title = "Evolución del accuracy del modelo LogitBoost") +
  theme_bw()

# PREDICCIONES SOBRE TEST 
# ==============================================================================
predic_LogitBoost_pvalue_500 <- predict(object = LogitBoost_pvalue_500, newdata = x_test)
caret::confusionMatrix(predic_LogitBoost_pvalue_500, y_test)


#===============================================================================
#=========================== Stacking Ensemble  ================================
#===============================================================================

# DEFINICIÓN DEL ENTRENAMIENTO
#===============================================================================


trainControl <- trainControl(method = "repeatedcv", number = 5, repeats = 3,
                             savePredictions = "all", classProbs = TRUE,
                             index = createFolds(stacking_df$Class, 5),verboseIter = TRUE)


# Run multiplealgorithms in one call.
algorithm <- c("nnet","rf","svmRadial","gbm","xgbTree","hdda","LogitBoost")


stacking_df <- cbind(x_train[,filtrado_anova_pvalue_500],x_train$Class)
names(stacking_df) [501] <- "Class"
levels(stacking_df$Class)
levels(stacking_df$Class) <- make.names(levels(stacking_df$Class))

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
                   sizes = c(1:286), 
                   rfeControl = rfe_control)

# summarize the results
print(RFE_results)

# list the chosen features
predictors(RFE_results)

# plot the results
plot(RFE_results, type=c("g", "o"))

selected_vars <- RFE_results$variables
write.csv(selected_vars, "STAD_selected_vars_RFE.csv")
saveRDS(RFE_results, "STAD_RFE_results.rds")
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
STAD_rf_RFE_50 <- caret::train(
  form = Class ~ .,
  data = x_train_rfe_50,
  method = "ranger",
  tuneGrid = hyperparameters,
  metric = "Accuracy",
  trControl = control_train,
  num.trees = 500)

saveRDS(object = STAD_rf_RFE_50, file = "STAD_rf_RFE_50.rds")
registerDoMC(cores = 1)


STAD_rf_RFE_50

# GRAPHIC REPRESENTATION
# ==============================================================================
ggplot(STAD_rf_RFE_50, highlight = TRUE) +
  labs(title = "Evolution of the accuracy of the Random Forest model") +
  guides(color = guide_legend(title = "mtry"),
         shape = guide_legend(title = "mtry")) +
  theme_bw()

# TEST PREDICTIONS
# ==============================================================================
predic_STAD_rf_RFE_50 <- predict(object = STAD_rf_RFE_50, newdata = x_test_rfe_50)
predic_STAD_rf_RFE_50
y_test <- as.factor(y_test)
caret::confusionMatrix(predic_STAD_rf_RFE_50, y_test)


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
STAD_svmrad_RFE_50 <- caret::train(
  form = Class ~ .,
  data = x_train_rfe_50,
  method = "svmRadial",
  tuneGrid = hyperparameters,
  metric = "Accuracy",
  trControl = control_train
)

saveRDS(object = STAD_svmrad_RFE_50, file = "STAD_svmrad_RFE_50.rds")
registerDoMC(cores = 1)


STAD_svmrad_RFE_50

# GRAPHIC REPRESENTATION
# ==============================================================================
ggplot(STAD_svmrad_RFE_50, highlight = TRUE) +
  labs(title = "Evolution of the accuracy of the SVM model") +
  guides(color = guide_legend(title = "Sigma"),
         shape = guide_legend(title = "Sigma")) +
  theme_bw()

# TEST PREDICTIONS
# ==============================================================================
predic_STAD_svmrad_RFE_50 <- predict(object = STAD_svmrad_RFE_50, newdata = x_test_rfe_50)
predic_STAD_svmrad_RFE_50
y_test <- as.factor(y_test)
caret::confusionMatrix(predic_STAD_svmrad_RFE_50, y_test)


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
STAD_nnet_RFE_50 <- caret::train(
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

saveRDS(object = STAD_nnet_RFE_50, file = "STAD_nnet_RFE_50.rds")
registerDoMC(cores = 1)

STAD_nnet_RFE_50

# GRAPHIC REPRESENTATION
# ==============================================================================
ggplot(STAD_nnet_RFE_50, highlight = TRUE) +
  labs(title = "Evolution of the accuracy of the Neural Net model") +
  guides(color = guide_legend(title = "Decay"),
         shape = guide_legend(title = "Decay")) +
  theme_bw()

# TEST PREDICTIONS
# ==============================================================================
predic_STAD_nnet_RFE_50 <- predict(object = STAD_nnet_RFE_50, newdata = x_test_rfe_50)
predic_STAD_nnet_RFE_50
y_test <- as.factor(y_test)
caret::confusionMatrix(predic_STAD_nnet_RFE_50, y_test)



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
STAD_gbm_RFE_50 <- caret::train(
  form = Class ~ .,
  data = x_train_rfe_50,
  method = "gbm",
  tuneGrid = hyperparameters,
  metric = "Accuracy",
  preProcess = c("center", "scale"),
  trControl = control_train
)

saveRDS(object = STAD_gbm_RFE_50, file = "STAD_gbm_RFE_50.rds")
registerDoMC(cores = 1)

STAD_gbm_RFE_50

# GRAPHIC REPRESENTATION
# ==============================================================================
ggplot(STAD_gbm_RFE_50, highlight = TRUE) +
  labs(title = "Evolution of the accuracy of the GBM model") +
  guides(color = guide_legend(title = "Shrinkage"),
         shape = guide_legend(title = "Shrinkage")) +
  theme_bw()

# TEST PREDICTIONS
# ==============================================================================
predic_STAD_gbm_RFE_50 <- predict(object = STAD_gbm_RFE_50, newdata = x_test_rfe_50)
predic_STAD_gbm_RFE_50
y_test <- as.factor(y_test)
caret::confusionMatrix(predic_STAD_gbm_RFE_50, y_test)


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
STAD_xgbm_RFE_50 <- caret::train(
  form = Class ~ .,
  data = x_train_rfe_50,
  method = "xgbTree",
  tuneGrid = hyperparameters,
  metric = "Accuracy",
  preProcess = c("center", "scale"),
  trControl = control_train
)

saveRDS(object = STAD_xgbm_RFE_50, file = "STAD_xgbm_RFE_50.rds")
registerDoMC(cores = 1)

STAD_xgbm_RFE_50

# GRAPHIC REPRESENTATION
# ==============================================================================
plot(STAD_xgbm_RFE_50, highlight = TRUE)

# TEST PREDICTIONS
# ==============================================================================
predic_STAD_xgbm_RFE_50 <- predict(object = STAD_xgbm_RFE_50, newdata = x_test_rfe_50)
predic_STAD_xgbm_RFE_50
y_test <- as.factor(y_test)
caret::confusionMatrix(predic_STAD_xgbm_RFE_50, y_test)


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
STAD_glmnet_RFE_50 <- caret::train(
  form = Class ~ .,
  data = x_train_rfe_50,
  method = "glmnet",
  family = "multinomial",
  tuneGrid = hyperparameters,
  metric = "Accuracy",
  trControl = control_train
)

saveRDS(object = STAD_glmnet_RFE_50, file = "STAD_glmnet_RFE_50.rds")
registerDoMC(cores = 1)

STAD_glmnet_RFE_50

# GRAPHIC REPRESENTATION
# ==============================================================================
ggplot(STAD_glmnet_RFE_50, highlight = TRUE) +
  labs(title = "Evolution of the accuracy of the GLMnet model") +
  guides(color = guide_legend(title = "Alpha"),
         shape = guide_legend(title = "Alpha")) +
  theme_bw()

# TEST PREDICTIONS
# ==============================================================================
predic_STAD_glmnet_RFE_50 <- predict(object = STAD_glmnet_RFE_50, newdata = x_test_rfe_50)
predic_STAD_glmnet_RFE_50
y_test <- as.factor(y_test)
caret::confusionMatrix(predic_STAD_glmnet_RFE_50, y_test)


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
STAD_hdda_RFE_50 <- caret::train(
  form = Class ~ .,
  data = x_train_rfe_50,
  method = "hdda",
  tuneGrid = hyperparameters,
  metric = "Accuracy",
  trControl = control_train
)

saveRDS(object = STAD_hdda_RFE_50, file = "STAD_hdda_RFE_50.rds")
registerDoMC(cores = 1)

STAD_hdda_RFE_50

# GRAPHIC REPRESENTATION
# ==============================================================================
ggplot(STAD_hdda_RFE_50, highlight = TRUE) +
  labs(title = "Evolution of the accuracy of the HDDA model") +
  guides(color = guide_legend(title = "threshold"),
         shape = guide_legend(title = "threshold")) +
  theme_bw()

# TEST PREDICTIONS
# ==============================================================================
predic_STAD_hdda_RFE_50 <- predict(object = STAD_hdda_RFE_50, newdata = x_test_rfe_50)
predic_STAD_hdda_RFE_50
y_test <- as.factor(y_test)
caret::confusionMatrix(predic_STAD_hdda_RFE_50, y_test)


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
STAD_logitBoost_RFE_50 <- caret::train(
  form = Class ~ .,
  data = x_train_rfe_50,
  method = "LogitBoost",
  tuneGrid = hyperparameters,
  metric = "Accuracy",
  trControl = control_train
)

saveRDS(object = STAD_logitBoost_RFE_50, file = "STAD_logitBoost_RFE_50.rds")
registerDoMC(cores = 1)

STAD_logitBoost_RFE_50

# GRAPHIC REPRESENTATION
# ==============================================================================
ggplot(STAD_logitBoost_RFE_50, highlight = TRUE) +
  labs(title = "Evolution of the accuracy of the Logitboost model") +
  guides(color = guide_legend(title = "nIter"),
         shape = guide_legend(title = "nIter")) +
  theme_bw()

# TEST PREDICTIONS
# ==============================================================================
predic_STAD_logitBoost_RFE_50 <- predict(object = STAD_logitBoost_RFE_50, newdata = x_test_rfe_50)
predic_STAD_logitBoost_RFE_50
y_test <- as.factor(y_test)
caret::confusionMatrix(predic_STAD_logitBoost_RFE_50, y_test)


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
STAD_rf_RFE_100 <- caret::train(
  form = Class ~ .,
  data = x_train_rfe_100,
  method = "ranger",
  tuneGrid = hyperparameters,
  metric = "Accuracy",
  trControl = control_train,
  num.trees = 500)

saveRDS(object = STAD_rf_RFE_100, file = "STAD_rf_RFE_100.rds")
registerDoMC(cores = 1)


STAD_rf_RFE_100

# GRAPHIC REPRESENTATION
# ==============================================================================
ggplot(STAD_rf_RFE_100, highlight = TRUE) +
  labs(title = "Evolution of the accuracy of the Random Forest model") +
  guides(color = guide_legend(title = "mtry"),
         shape = guide_legend(title = "mtry")) +
  theme_bw()

# TEST PREDICTIONS
# ==============================================================================
predic_STAD_rf_RFE_100 <- predict(object = STAD_rf_RFE_100, newdata = x_test_rfe_100)
predic_STAD_rf_RFE_100
y_test <- as.factor(y_test)
caret::confusionMatrix(predic_STAD_rf_RFE_100, y_test)


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
STAD_svmrad_RFE_100 <- caret::train(
  form = Class ~ .,
  data = x_train_rfe_100,
  method = "svmRadial",
  tuneGrid = hyperparameters,
  metric = "Accuracy",
  trControl = control_train
)

saveRDS(object = STAD_svmrad_RFE_100, file = "STAD_svmrad_RFE_100.rds")
registerDoMC(cores = 1)


STAD_svmrad_RFE_100

# GRAPHIC REPRESENTATION
# ==============================================================================
ggplot(STAD_svmrad_RFE_100, highlight = TRUE) +
  labs(title = "Evolution of the accuracy of the SVM model") +
  guides(color = guide_legend(title = "Sigma"),
         shape = guide_legend(title = "Sigma")) +
  theme_bw()

# TEST PREDICTIONS
# ==============================================================================
predic_STAD_svmrad_RFE_100 <- predict(object = STAD_svmrad_RFE_100, newdata = x_test_rfe_100)
predic_STAD_svmrad_RFE_100
y_test <- as.factor(y_test)
caret::confusionMatrix(predic_STAD_svmrad_RFE_100, y_test)


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
STAD_nnet_RFE_100 <- caret::train(
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

saveRDS(object = STAD_nnet_RFE_100, file = "STAD_nnet_RFE_100.rds")
registerDoMC(cores = 1)

STAD_nnet_RFE_100

# GRAPHIC REPRESENTATION
# ==============================================================================
ggplot(STAD_nnet_RFE_100, highlight = TRUE) +
  labs(title = "Evolution of the accuracy of the Neural Net model") +
  guides(color = guide_legend(title = "Decay"),
         shape = guide_legend(title = "Decay")) +
  theme_bw()

# TEST PREDICTIONS
# ==============================================================================
predic_STAD_nnet_RFE_100 <- predict(object = STAD_nnet_RFE_100, newdata = x_test_rfe_100)
predic_STAD_nnet_RFE_100
y_test <- as.factor(y_test)
caret::confusionMatrix(predic_STAD_nnet_RFE_100, y_test)



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
STAD_gbm_RFE_100 <- caret::train(
  form = Class ~ .,
  data = x_train_rfe_100,
  method = "gbm",
  tuneGrid = hyperparameters,
  metric = "Accuracy",
  preProcess = c("center", "scale"),
  trControl = control_train
)

saveRDS(object = STAD_gbm_RFE_100, file = "STAD_gbm_RFE_100.rds")
registerDoMC(cores = 1)

STAD_gbm_RFE_100

# GRAPHIC REPRESENTATION
# ==============================================================================
ggplot(STAD_gbm_RFE_100, highlight = TRUE) +
  labs(title = "Evolution of the accuracy of the GBM model") +
  guides(color = guide_legend(title = "Shrinkage"),
         shape = guide_legend(title = "Shrinkage")) +
  theme_bw()

# TEST PREDICTIONS
# ==============================================================================
predic_STAD_gbm_RFE_100 <- predict(object = STAD_gbm_RFE_100, newdata = x_test_rfe_100)
predic_STAD_gbm_RFE_100
y_test <- as.factor(y_test)
caret::confusionMatrix(predic_STAD_gbm_RFE_100, y_test)


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
STAD_xgbm_RFE_100 <- caret::train(
  form = Class ~ .,
  data = x_train_rfe_100,
  method = "xgbTree",
  tuneGrid = hyperparameters,
  metric = "Accuracy",
  preProcess = c("center", "scale"),
  trControl = control_train
)

saveRDS(object = STAD_xgbm_RFE_100, file = "STAD_xgbm_RFE_100.rds")
registerDoMC(cores = 1)

STAD_xgbm_RFE_100

# GRAPHIC REPRESENTATION
# ==============================================================================
plot(STAD_xgbm_RFE_100, highlight = TRUE)

# TEST PREDICTIONS
# ==============================================================================
predic_STAD_xgbm_RFE_100 <- predict(object = STAD_xgbm_RFE_100, newdata = x_test_rfe_100)
predic_STAD_xgbm_RFE_100
y_test <- as.factor(y_test)
caret::confusionMatrix(predic_STAD_xgbm_RFE_100, y_test)


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
STAD_glmnet_RFE_100 <- caret::train(
  form = Class ~ .,
  data = x_train_rfe_100,
  method = "glmnet",
  family = "multinomial",
  tuneGrid = hyperparameters,
  metric = "Accuracy",
  trControl = control_train
)

saveRDS(object = STAD_glmnet_RFE_100, file = "STAD_glmnet_RFE_100.rds")
registerDoMC(cores = 1)

STAD_glmnet_RFE_100

# GRAPHIC REPRESENTATION
# ==============================================================================
ggplot(STAD_glmnet_RFE_100, highlight = TRUE) +
  labs(title = "Evolution of the accuracy of the GLMnet model") +
  guides(color = guide_legend(title = "Alpha"),
         shape = guide_legend(title = "Alpha")) +
  theme_bw()

# TEST PREDICTIONS
# ==============================================================================
predic_STAD_glmnet_RFE_100 <- predict(object = STAD_glmnet_RFE_100, newdata = x_test_rfe_100)
predic_STAD_glmnet_RFE_100
y_test <- as.factor(y_test)
caret::confusionMatrix(predic_STAD_glmnet_RFE_100, y_test)


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
STAD_hdda_RFE_100 <- caret::train(
  form = Class ~ .,
  data = x_train_rfe_100,
  method = "hdda",
  tuneGrid = hyperparameters,
  metric = "Accuracy",
  trControl = control_train
)

saveRDS(object = STAD_hdda_RFE_100, file = "STAD_hdda_RFE_100.rds")
registerDoMC(cores = 1)

STAD_hdda_RFE_100

# GRAPHIC REPRESENTATION
# ==============================================================================
ggplot(STAD_hdda_RFE_100, highlight = TRUE) +
  labs(title = "Evolution of the accuracy of the HDDA model") +
  guides(color = guide_legend(title = "threshold"),
         shape = guide_legend(title = "threshold")) +
  theme_bw()

# TEST PREDICTIONS
# ==============================================================================
predic_STAD_hdda_RFE_100 <- predict(object = STAD_hdda_RFE_100, newdata = x_test_rfe_100)
predic_STAD_hdda_RFE_100
y_test <- as.factor(y_test)
caret::confusionMatrix(predic_STAD_hdda_RFE_100, y_test)


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
STAD_logitBoost_RFE_100 <- caret::train(
  form = Class ~ .,
  data = x_train_rfe_100,
  method = "LogitBoost",
  tuneGrid = hyperparameters,
  metric = "Accuracy",
  trControl = control_train
)

saveRDS(object = STAD_logitBoost_RFE_100, file = "STAD_logitBoost_RFE_100.rds")
registerDoMC(cores = 1)

STAD_logitBoost_RFE_100

# GRAPHIC REPRESENTATION
# ==============================================================================
ggplot(STAD_logitBoost_RFE_100, highlight = TRUE) +
  labs(title = "Evolution of the accuracy of the Logitboost model") +
  guides(color = guide_legend(title = "nIter"),
         shape = guide_legend(title = "nIter")) +
  theme_bw()

# TEST PREDICTIONS
# ==============================================================================
predic_STAD_logitBoost_RFE_100 <- predict(object = STAD_logitBoost_RFE_100, newdata = x_test_rfe_100)
predic_STAD_logitBoost_RFE_100
y_test <- as.factor(y_test)
caret::confusionMatrix(predic_STAD_logitBoost_RFE_100, y_test)


#===============================================================================
#=========================== ROC curves ========================================
#===============================================================================


nnet_predict_obj <- mmdata(as.numeric(predic_nnet_pvalue_500),x_test$Class)
nnet_performance <- evalmod(mdat = nnet_predict_obj) 

hdda_predict_obj <- mmdata(as.numeric(predic_hdda_pvalue_50),x_test$Class)
hdda_performance <- evalmod(mdat = hdda_predict_obj) 


nnet_df <- fortify(nnet_performance)
hdda_df <- fortify(hdda_performance)

nnet_df$classifier <- "nnet"
hdda_df$classifier <- "hdda"

performance_df <- rbind(nnet_df, hdda_df)

roc <- performance_df[performance_df$curvetype == "ROC",]

ggplot(roc, aes(x=x, y=y, group = classifier)) + 
  geom_line(aes(color = classifier)) +
  xlab("1 - Specifity") +
  ylab("Sensitivity") +
  ggtitle("ROC curve") +
  theme_minimal()

# Another way to plot 
#===================================================================


nnet_predict_obj <- mmdata(as.numeric(predic_nnet_pvalue_500),x_test$Class)
nnet_performance <- evalmod(mdat = nnet_predict_obj) 
plot(nnet_performance)

precrec_obj2 <- evalmod(scores = as.numeric(predic_nnet_pvalue_500), 
                        labels = x_test$Class, mode="basic")
plot(precrec_obj2)   


# Another way to plot 
#===================================================================

ROCit_obj <- rocit(score = as.numeric(predic_nnet_pvalue_500), 
                   class = x_test$Class)
plot(ROCit_obj)


ROCit_obj <- rocit(score = as.numeric(predic_hdda_pvalue_50), 
                   class = x_test$Class)
plot(ROCit_obj)
ksplot(ROCit_obj)

# Another way to plot 
#===================================================================

rocplot <- ggplot(df, aes(m = as.numeric(predic_hdda_pvalue_50), d = labels))+ geom_roc(n.cuts=20,labels=FALSE)
rocplot + style_roc(theme = theme_grey) + geom_rocci(fill="pink") 


# Another way to plot 
#===================================================================

ROC(test=predic_hdda_pvalue_50, stat=x_test$Class, plot="ROC", AUC=T, main="hdda")


# Another way to plot 
#===================================================================

#' Generate an ROC curve plot with error bars showing 95 percent
#' confidence intervals
#'
#' This code builds off of code written by Vincent Guillemot found here:
#' \url{https://rpubs.com/vguillem/465086}.
#'
#' @param df The df as a data.frame.
#' @param outcome A character string containing the name of the column
#'   containing the outcomes (expressed as 0/1s).
#' @param prediction A character string containing the name of the column
#'   containing the predictions.
#' @param ci Show confidence interval ribbon.
#'   Defaults to FALSE.
#' @param plot_title A character string containing the title for the resulting
#'   plot.
#' @return A ggplot containing the calibration plot
#' @examples
#' data(single_model_dataset)
#' roc_plot(single_model_dataset, outcome = 'outcomes', prediction = 'predictions', ci = TRUE)
#' @export

roc_plot <- function(df, outcome, prediction, ci = FALSE, plot_title = '') {
  obj <- pROC::roc_(df, response = outcome, predictor = prediction, ci = ci, plot=FALSE)
  ciobj <- pROC::ci.se(obj, specificities = seq(0, 1, l = 25))
  dat.ci <- data.frame(x = as.numeric(rownames(ciobj)),
                       lower = ciobj[, 1],
                       upper = ciobj[, 3])
  
  g1 = pROC::ggroc(obj) +
    ggplot2::theme_minimal() +
    ggplot2::geom_abline(
      slope = 1,
      intercept = 1,
      linetype = "dashed",
      alpha = 0.7,
      color = "grey"
    ) +
    ggplot2::coord_equal() +
    ggplot2::ggtitle(plot_title) + 
    ggplot2::xlab("1 - Specificity") + 
    ggplot2::ylab("Sensitivity")
  
  if(ci){
    g2 = g1 + ggplot2::geom_ribbon(
      data = dat.ci,
      ggplot2::aes(x = x, ymin = lower, ymax = upper),
      # fill = "steelblue",
      alpha = 0.2
    )
  } else g2 = g1
  g2
}

#' Generate an ROC curve plot with error bars showing 95 percent
#' confidence intervals
#' @param df The df as a data.frame.
#' @param outcome A character string containing the name of the column
#'   containing the outcomes (expressed as 0/1s).
#' @param prediction A character string containing the name of the column
#'   containing the predictions.
#' @param model A character string containing the name of the column
#'   containing the model label.
#' @param ci Show confidence interval ribbon.
#'   Defaults to FALSE.
#' @param plot_title A character string containing the title for the resulting
#'   plot.
#' @return A ggplot containing the ROC plot
#' @examples
#' data(multi_model_dataset)
#' roc_plot_multi(multi_model_dataset, outcome = 'outcomes', prediction = 'predictions', model = 'model_name', ci = TRUE)
#' @export
roc_plot_multi <- function(df, outcome, prediction, model, ci = FALSE, plot_title = '') {
  
  how_many_models = df[[model]] %>% unique() %>% length()
  
  ci_data = df %>%
    dplyr::group_by(!!rlang::parse_expr(model)) %>%
    dplyr::group_nest() %>%
    dplyr::mutate(roc = purrr::map(data, ~ pROC::roc_(data = ., response = outcome, predictor = prediction, ci = ci, plot=FALSE)),
                  roc = roc %>% setNames(!!rlang::parse_expr(model)),
                  ci_spec = purrr::map(roc, ~ pROC::ci.se(., specificities = seq(0, 1, l = 25))),
                  ci_ribbon = purrr::map(ci_spec, ~ build_ci_data(.)))
  
  # build ROC curves
  g1 = pROC::ggroc(ci_data$roc) +
    ggplot2::theme_minimal() +
    ggplot2::geom_abline(
      slope = 1,
      intercept = 1,
      linetype = "dashed",
      alpha = 0.7,
      color = "grey"
    ) +
    ggplot2::coord_equal() +
    ggplot2::scale_color_brewer(name = 'Models', palette = 'Set1') +
    ggplot2::scale_fill_brewer(name = 'Models', palette = 'Set1') +
    ggplot2::ggtitle(plot_title) + 
    ggplot2::xlab("Specificity") + 
    ggplot2::ylab("Sensitivity")
  
  # build CI intervals
  if(ci){
    ribbon = ci_data %>%
      dplyr::select(!!rlang::parse_expr(model), ci_ribbon) %>%
      dplyr::rename(name = !!rlang::parse_expr(model)) %>%
      tidyr::unnest_wider(ci_ribbon) %>%
      tidyr::unnest(cols = c(x, lower, upper))
    
    # ci_values = ci_data %>%
    #   dplyr::select(model_name, roc) %>%
    #   dplyr::mutate(ci = purrr::map(roc, ~ .$ci)) %>%
    #   dplyr::select(-roc)
    
    ci_values = ci_data %>%
      dplyr::select(!!rlang::parse_expr(model), roc) %>%
      dplyr::mutate(ci = purrr::map(roc, ~ pROC::ci(.)))
    
    g2 = g1 + ggplot2::geom_ribbon(data = ribbon,
                                   ggplot2::aes(fill = name, x = x, ymin = lower, ymax = upper),
                                   alpha = 1/how_many_models,
                                   inherit.aes = FALSE)
  } else g2 = g1
  g2
}

build_ci_data = function(obj){
  data.frame(x = as.numeric(rownames(obj)),
             lower = obj[, 1],
             upper = obj[, 3])
}


df <- data.frame(predictions = as.numeric(predic_hdda_pvalue_50), outcomes = x_test$Class)
roc_plot(df, outcome = 'outcomes', prediction = 'predictions', ci = TRUE, plot_title = "ROC Curve")

