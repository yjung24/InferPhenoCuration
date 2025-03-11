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
combined_data <- readRDS("~/InferPhenoCuration/MSI_Status/data/combined_data.rds")

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
heatmapTable(validate_data, RAVmodel)

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

# Numeric Variables
## Calculate R-Squared Value
## R squared value function
calculateRsq <- function (x, y) stats::cor(x, y, use = 'na.or.complete') ^ 2

## Calculate r-squared for numeric attributes
rsq_numAttr <- as.data.frame(matrix(nrow = ncol(numTb),
                                    ncol = ncol(sampleScore_sub)))

colnames(rsq_numAttr) <- colnames(sampleScore_sub)
rownames(rsq_numAttr) <- colnames(numTb)

for (i in seq_len(ncol(numTb))) {
  for (j in seq_len(ncol(sampleScore_sub))) {
    rsq <- calculateRsq(numTb[[i]], sampleScore_sub[, j])
    rsq_numAttr[i, j] <- rsq
  }
}

rsq_numAttr <- na.omit(rsq_numAttr)

# Plot the high R-squred
cutoff <- 0.3 # cutoff of the minimum r-sq value

max_rav <- apply(rsq_numAttr, 1, max)
max_attr <- which(max_rav > cutoff) # the maximum r-sq of the sample score should be above this cutoff
target_rsq <- rsq_numAttr[max_attr,]
target_rsq <- as.data.frame(t(target_rsq))

target_rsq %>% filter(patient.number_of_abnormal_loci > 0.3)

# RAVs 172 and 22 also have R-squared values > 0.4 for number of abnormal loci but not included

## Use top validated RAVs for TCGA-COAD,STAD, UCEC ---------------
RAVs_combinedTCGA <- c(517, 220, 2109, 1303, 324, 438, 868, #RAVs that have statistically significant pairwise wilcoxon p-values of mss vs msi-h
                       834, 190, 1166, #RAVs with significant KW test statistic (p-value < 0.05) for COAD
                       2344, #significant KW test value for STAD, includes 324, 868, 517 above
                       357) #UCEC KW test value (p-value = 0.056)

validated_RAVs <- paste("RAV", RAVs_combinedTCGA, sep="")
#validated_RAVs <- c("RAV357", "RAV27", "RAV834", "RAV190", "RAV1166", "RAV517",
#"RAV2344", "RAV324", "RAV438", "RAV220", "RAV868", "RAV1008", "RAV625")
target_attr <- "patient.number_of_abnormal_loci"

combined_data_msi <- combined_data[,!is.na(combined_data[[target_attr]])]

#Creating new data frames for train and test data with only the validated ravs for samples
#containing nonNA values of the target attribute (number of abnormal loci)

###access the target_attr data from the hierarchical summarized experiment dataset 
labels <- numTb[[target_attr]]

###vector/list of the positions where values are not NA
nonNALabels <- which(!is.na(labels))

###subset sample score data for validate RAVs only
data <- sampleScore_sub[,validated_RAVs]

###subset sample scores values using the position labels of nonNA number of abnormal loci values
train_data <- data[nonNALabels,]

###abnormal number of loci values 
train_labels <- labels[nonNALabels]
train_data$MSI_numeric <- train_labels

###repeat for the test data
data_test <- sampleScore_sub_testdata[,validated_RAVs]
test_data <- data_test[nonNALabels,]
test_labels <- labels[nonNALabels]
test_data$MSI_numeric <- test_labels

# train model
rf_model <- randomForest(MSI_numeric ~ ., data = train_data, ntree = 500, 
                         keep.forest = TRUE, importance = TRUE)
lm <- lm(MSI_numeric ~ ., data = train_data)

#pred model
rf_pred <- predict(rf_model, test_data, predict.all = TRUE)
lm_pred <- predict(lm, test_data)

