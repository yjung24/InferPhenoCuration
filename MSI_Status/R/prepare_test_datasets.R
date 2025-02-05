### From https://github.com/bpheng/GSS-Britney-Fieldwork-2024/blob/main/Random%20Forest%20Classification%20Model.qmd 


## Using Test Data with the RF Model -----------------
labels <- factorTb_testdata[[target_attr]]
nonNALabels <- which(!is.na(labels))
data <- sampleScore_sub_testdata[,validated_RAVs]

test_data <- data[nonNALabels,]
test_labels <- labels[nonNALabels]

test_data$MSI_Status <- test_labels


# Test RF model with individual TCGA data cancer types: COAD, UCEC, STAD ----------
## COAD -----------------
coad_rna_cancer <- readRDS(file="~/Projects/InferPhenoCuration/MSI_Status/data/coad_rna_cancer.rds")
target_attr <- "patient.microsatellite_instability_test_results.microsatellite_instability_test_result.mononucleotide_and_dinucleotide_marker_panel_analysis_status"

# coad_rna_cancer <- coad_rna_cancer[,coad_rna_cancer@colData@listData$patient.microsatellite_instability_test_results.microsatellite_instability_test_result.mononucleotide_and_dinucleotide_marker_panel_analysis_status %in% c('msi-h', 'mss')]
coad_rna_cancer <- coad_rna_cancer[,coad_rna_cancer[[target_attr]] %in% c("msi-h", "mss")]
coad.meta <- colData(coad_rna_cancer)
coad.sampleScore <- calculateScore(coad_rna_cancer, RAVmodel) %>% as.data.frame()

var_type <- coad.meta@listData
unique(sapply(var_type, type))
charcTb <- coad.meta[, sapply(var_type, class) == 'character']
factorTb <- charcTb
factorTb[sapply(factorTb, is.character)] <- lapply(factorTb[sapply(factorTb, is.character)], factor)

levels(factorTb[[target_attr]])[levels(factorTb[[target_attr]]) == "msi-h"] <- "msi"
labels <- factorTb[[target_attr]]
nonNALabels <- which(!is.na(labels))
coad_data <- coad.sampleScore[,validated_RAVs]

new_coad_data <- coad_data[nonNALabels,]
new_labels <- labels[nonNALabels]
new_coad_data$msi <- new_labels

## UCEC -------------
ucec_rna_cancer <- readRDS(file="~/Projects/InferPhenoCuration/MSI_Status/data/ucec_rna_cancer.rds")

# ucec_rna_cancer <- ucec_rna_cancer[,ucec_rna_cancer@colData@listData$patient.microsatellite_instability_test_results.microsatellite_instability_test_result.mononucleotide_and_dinucleotide_marker_panel_analysis_status %in% c('msi-h', 'mss')]
ucec_rna_cancer <- ucec_rna_cancer[,ucec_rna_cancer[[target_attr]] %in% c("msi-h", "mss")]
ucec.meta <- colData(ucec_rna_cancer)
ucec.sampleScore <- calculateScore(ucec_rna_cancer, RAVmodel) %>% as.data.frame()

var_type <- ucec.meta@listData
unique(sapply(var_type, type))
charcTb <- ucec.meta[, sapply(var_type, class) == 'character']
factorTb <- charcTb
factorTb[sapply(factorTb, is.character)] <- lapply(factorTb[sapply(factorTb, is.character)], factor)

levels(factorTb[[target_attr]])[levels(factorTb[[target_attr]]) == "msi-h"] <- "msi"
labels <- factorTb[[target_attr]]
nonNALabels <- which(!is.na(labels))
ucec_data <- ucec.sampleScore[,validated_RAVs]

new_ucec_data <- ucec_data[nonNALabels,]
new_labels <- labels[nonNALabels]
new_ucec_data$msi <- new_labels


## STAD ---------------
stad_rna_cancer <- readRDS(file="~/Projects/InferPhenoCuration/MSI_Status/data/stad_rna_cancer.rds")
stad_rna_cancer <- stad_rna_cancer[,stad_rna_cancer[[target_attr]] %in% c("msi-h", "mss")]
stad.meta <- colData(stad_rna_cancer)
stad.sampleScore <- calculateScore(stad_rna_cancer, RAVmodel) %>% as.data.frame()

var_type <- stad.meta@listData
unique(sapply(var_type, type))
charcTb <- stad.meta[, sapply(var_type, class) == 'character']
factorTb <- charcTb
factorTb[sapply(factorTb, is.character)] <- lapply(factorTb[sapply(factorTb, is.character)], factor)

levels(factorTb[[target_attr]])[levels(factorTb[[target_attr]]) == "msi-h"] <- "msi"
labels <- factorTb[[target_attr]]
nonNALabels <- which(!is.na(labels))
stad_data <- stad.sampleScore[,validated_RAVs]

new_stad_data <- stad_data[nonNALabels,]
new_labels <- labels[nonNALabels]
new_stad_data$msi <- new_labels


# # Test with curatedCRCData
# library(curatedCRCData)
# data(package="curatedCRCData")
# 
# ## curatedCRC Dataset 1: GSE13067_eset --------
# data(GSE13067_eset)
# library(SummarizedExperiment)
# 
# mySummarizedExperiment <- makeSummarizedExperimentFromExpressionSet(GSE13067_eset)
# assay(mySummarizedExperiment) <- log2(assay(mySummarizedExperiment) + 1)
# 
# mySummarizedExperiment@colData@listData[["msi"]] <- gsub("MSS", "mss", mySummarizedExperiment@colData@listData[["msi"]] )
# mySummarizedExperiment@colData@listData[["msi"]] <- gsub("MSI", "msi", mySummarizedExperiment@colData@listData[["msi"]] )
# 
# GSE13067.meta <- colData(mySummarizedExperiment)
# GSE13067.sampleScore <- calculateScore(mySummarizedExperiment, RAVmodel) %>% as.data.frame()
# 
# target_attr <- "msi"
# 
# var_type <- GSE13067.meta@listData
# unique(sapply(var_type, type))
# charcTb <- GSE13067.meta[, sapply(var_type, class) == 'character']
# factorTb <- charcTb
# factorTb[sapply(factorTb, is.character)] <- lapply(factorTb[sapply(factorTb, is.character)], factor)
# 
# labels <- factorTb[[target_attr]]
# nonNALabels <- which(!is.na(labels))
# data <- GSE13067.sampleScore[,validated_RAVs]
# 
# new_data <- data[nonNALabels,]
# new_labels <- labels[nonNALabels]
# new_data$msi <- new_labels
# 
# 
# 
# 
# ## curatedCRC Dataset 2: GSE13294_eset
# 
# ```{r}
# data(GSE13294_eset)
# GSE13294_eset <- makeSummarizedExperimentFromExpressionSet(GSE13294_eset)
# assay(GSE13294_eset) <- log2(assay(GSE13294_eset) + 1)
# 
# 
# GSE13294_eset@colData@listData[["msi"]] <- gsub("MSS", "mss", GSE13294_eset@colData@listData[["msi"]] )
# GSE13294_eset@colData@listData[["msi"]] <- gsub("MSI", "msi", GSE13294_eset@colData@listData[["msi"]] )
# 
# GSE13294.meta <- colData(GSE13294_eset)
# GSE13294.sampleScore <- calculateScore(GSE13294_eset, RAVmodel) %>% as.data.frame()
# 
# ```
# 
# ```{r include=FALSE}
# target_attr <- "msi"
# 
# var_type <- GSE13294.meta@listData
# unique(sapply(var_type, type))
# charcTb <- GSE13294.meta[, sapply(var_type, class) == 'character']
# 
# factorTb <- charcTb
# factorTb[sapply(factorTb, is.character)] <- lapply(factorTb[sapply(factorTb, is.character)], factor)
# 
# labels <- factorTb[[target_attr]]
# nonNALabels <- which(!is.na(labels))
# data <- GSE13294.sampleScore[,validated_RAVs]
# 
# new_data <- data[nonNALabels,]
# new_labels <- labels[nonNALabels]
# new_data$msi <- new_labels
# ```
# 
# ```{r}
# GSE13294.prediction <- predict(rf, new_data)
# confusionMatrix(GSE13294.prediction, new_data$msi)
# ```
