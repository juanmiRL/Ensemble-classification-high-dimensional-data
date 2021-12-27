########################################################################
####################### TCGA Data ######################################
########################################################################


library(TCGAbiolinks)
library(tidyverse)
library(MLSeq)
library(S4Vectors)
library(DESeq2)


GDCprojects = getGDCprojects()
names_projects <- GDCprojects[c("project_id", "name")]

# Stomach Adenocarcinoma
#--------------------------------------------------------------------------

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

saveRDS(object = tcga_data,
        file = "tcga_data_STAD.RDS",
        compress = FALSE)



# Glioblastoma Multiforme
#--------------------------------------------------------------------------

TCGAbiolinks:::getProjectSummary("TCGA-GBM")

query_TCGA = GDCquery(
  project = "TCGA-GBM",
  data.category = "Transcriptome Profiling", # parameter enforced by GDCquery
  experimental.strategy = "RNA-Seq",
  workflow.type = "HTSeq - Counts")

GBM_res = getResults(query_TCGA) # make results as table
head(GBM_res) # data of the first 6 patients.
colnames(GBM_res) # columns present in the table

head(GBM_res$sample_type) # first 6 types of tissue.
table(GBM_res$sample_type) # summary of distinct tissues types present in this study

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


# assay function 
head(assay(tcga_data)[,1:10]) # expression of first 6 genes and first 10 samples
head(rowData(tcga_data))     # ensembl id and gene id of the first 6 genes.

# Save the data as a file, if you need it later, you can just load this file
# instead of having to run the whole pipeline again

saveRDS(object = tcga_data,
        file = "tcga_data_STAD.RDS",
        compress = FALSE)
