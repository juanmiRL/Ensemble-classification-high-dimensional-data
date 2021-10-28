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

#===============================================================================
#=========================== STDA ============================================== 
#===============================================================================

# load data with the next sentence
tcga_data = readRDS(file = "tcga_data_STAD.RDS")


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

# Se guarda en disco el objeto creado para no tener que repetir de nuevo toda la
# computación.
saveRDS(object = resultados_anova_pvalue, file = "resultados_anova_pvalue.rds")

resultados_anova_pvalue %>% dplyr::select(1,2,3,4) %>% head()

# Se filtran los 100, 50 y 500 genes identificados como más relevantes mediante anova
filtrado_anova_pvalue_500 <-  resultados_anova_pvalue %>% pull(gen) %>% head(500)
filtrado_anova_pvalue_100 <- resultados_anova_pvalue %>% pull(gen) %>% head(100)
filtrado_anova_pvalue_50  <- resultados_anova_pvalue %>% pull(gen) %>% head(50)


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
  threshold = seq(0,1,by=0.1),
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
