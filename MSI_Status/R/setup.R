### From https://github.com/bpheng/GSS-Britney-Fieldwork-2024/blob/main/Random%20Forest%20Classification%20Model.qmd 


## Packages -------------
suppressPackageStartupMessages({
    # BiocManager
    library(GenomicSuperSignature)
    library(curatedTCGAData)
    library(MultiAssayExperiment)
    library(TCGAutils)
    library(ComplexHeatmap)
    
    # CRAN
    library(tidyverse) # includes dplyr, ggplot2, magrittr, tidyr
    library(magick)
    library(wordcloud)
    library(ztable)
    library(metafolio)
    library(randomForest)
    library(caret)
})


## Load RAVmodel -------
RAVmodel <- getModel('C2', load=TRUE)


## Load combined_data (COAD, STAD, UCEC cancer sample data) ----------
combined_data <- readRDS("~/Projects/InferPhenoCuration/MSI_Status/data/combined_data.rds")

tcga_dat <- combined_data[,combined_data@colData@listData$patient.microsatellite_instability_test_results.microsatellite_instability_test_result.mononucleotide_and_dinucleotide_marker_panel_analysis_status %in% c('msi-h', 'mss')]

combined_data_meta <- colData(tcga_dat)


## Split combined_data for training ---------
# Randomly select 70% of samples for training
set.seed(1)
num_sample <- ncol(tcga_dat)
train_sample_ind <- sample(seq_len(num_sample), round(num_sample*0.7))

meta_train <- combined_data_meta[train_sample_ind,] 
data_train <- tcga_dat[, train_sample_ind] # 20,501 genes x 612 samples x 3205 attributes

# Remove batch effect variables from metadata table
batch_var <- "analyte|portion|procurement|aliquot|uuid|barcode"
batch_ind <- grep(batch_var, colnames(meta_train))

meta_train <- meta_train[,-batch_ind] # 612 samples x 1064 metadata attributes

meta_test <- combined_data_meta[-train_sample_ind,-batch_ind]
data_test <- tcga_dat[, -train_sample_ind]


## Split Variable Types for Train Data ------------
## Check for data types in listData
var_type <- meta_train@listData
unique(sapply(var_type, type))

## Separate training data's metadata into two subsets: 
## character variables (~ categorical) and numeric variables (~ continuous)
charcTb <- meta_train[, sapply(var_type, class) == 'character']
numTb <- meta_train[, sapply(var_type, class) %in% c('numeric', 'integer')]


## Split Variable Types for Test Data -------------
## Check for data types in listData
var_type_testdata <- meta_test@listData

## Separate test data's metadata into two subsets: 
## character variables (~ categorical) and numeric variables (~ continuous)
charcTb_testdata <- meta_test[, sapply(var_type_testdata, class) == 'character']
numTb_testdata <- meta_test[, sapply(var_type_testdata, class) %in% c('numeric', 'integer')]


## Sample Scores ---------------
## Calculate validation scores
sampleScore <- calculateScore(data_train, RAVmodel)
rownames(sampleScore) <- gsub("\\.", "-", rownames(sampleScore))

## Test data: calculate validation scores
sampleScore_testdata <- calculateScore(data_test, RAVmodel)
rownames(sampleScore_testdata) <- gsub("\\.", "-", rownames(sampleScore_testdata))


## Sample scores for all RAVs -----------
## Training Data
validate_data <- validate(data_train, RAVmodel)
#heatmapTable(validate_data, RAVmodel)

validated_ind <- validatedSignatures(validate_data,
                                     RAVmodel,
                                     num.out = 4764, # We want to validate all RAVs so we can select which ones to include in the RF model
                                     #scoreCutoff = 0.45,
                                     indexOnly = TRUE)
# saveRDS(validated_ind, file="/Users/bpheng/GSS-Britney-Fieldwork-2024/traindata_validatedRAVs.rds")
# validated_ind <- readRDS("~/Teaching/GSS-Britney-Fieldwork-2024/traindata_validatedRAVs.rds")
#Already saved as an RDS file, called above

## Subset sampleScore
sampleScore_sub <- sampleScore[,validated_ind] %>% as.data.frame()

## Test Data
validate_testdata <- validate(data_test, RAVmodel)
# heatmapTable(validate_testdata, RAVmodel, num.out = 15)
# 
validated_ind_testdata <- validatedSignatures(validate_testdata,
                                     RAVmodel,
                                     num.out = 4764,
                                     #scoreCutoff = 0.45,
                                     indexOnly = TRUE)
# 
# saveRDS(validated_ind_testdata, file="/Users/bpheng/GSS-Britney-Fieldwork-2024/testdata_validatedRAVs.rds")

# validated_ind_testdata <- readRDS("/Users/bpheng/GSS-Britney-Fieldwork-2024/testdata_validatedRAVs.rds")

## Subset sampleScore
sampleScore_sub_testdata <- sampleScore_testdata[,validated_ind_testdata] %>% as.data.frame()


## Clean up factor level name for microsatellite instability - high to “msi” -----------
## Convert character variables into the factor data type
factorTb <- charcTb
factorTb[sapply(factorTb, is.character)] <- lapply(factorTb[sapply(factorTb, is.character)], factor)

levels(factorTb@listData$patient.microsatellite_instability_test_results.microsatellite_instability_test_result.mononucleotide_and_dinucleotide_marker_panel_analysis_status)[levels(factorTb@listData$patient.microsatellite_instability_test_results.microsatellite_instability_test_result.mononucleotide_and_dinucleotide_marker_panel_analysis_status) == "msi-h"] <- "msi"

## Convert character variables into the factor data type
factorTb_testdata <- charcTb_testdata
factorTb_testdata[sapply(factorTb_testdata, is.character)] <- lapply(factorTb_testdata[sapply(factorTb_testdata, is.character)], factor)

levels(factorTb_testdata@listData$patient.microsatellite_instability_test_results.microsatellite_instability_test_result.mononucleotide_and_dinucleotide_marker_panel_analysis_status)[levels(factorTb_testdata@listData$patient.microsatellite_instability_test_results.microsatellite_instability_test_result.mononucleotide_and_dinucleotide_marker_panel_analysis_status) == "msi-h"] <- "msi"


## Use top validated RAVs for TCGA-COAD,STAD, UCEC ---------------
RAVs_combinedTCGA <- c(517, 220, 2109, 1303, 324, 438, 868, #RAVs that have statistically significant pairwise wilcoxon p-values of mss vs msi-h
                       834, 190, 1166, #RAVs with significant KW test statistic (p-value < 0.05) for COAD
                       2344, #significant KW test value for STAD, includes 324, 868, 517 above
                       357) #UCEC KW test value (p-value = 0.056)

validated_RAVs <- paste("RAV", RAVs_combinedTCGA, sep="")
#validated_RAVs <- c("RAV357", "RAV27", "RAV834", "RAV190", "RAV1166", "RAV517",
#"RAV2344", "RAV324", "RAV438", "RAV220", "RAV868", "RAV1008", "RAV625")
target_attr <- "patient.microsatellite_instability_test_results.microsatellite_instability_test_result.mononucleotide_and_dinucleotide_marker_panel_analysis_status"

labels <- factorTb[[target_attr]]
nonNALabels <- which(!is.na(labels))
data <- sampleScore_sub[,validated_RAVs]

train_data <- data[nonNALabels,]
train_labels <- labels[nonNALabels]
train_data$MSI_Status <- train_labels
