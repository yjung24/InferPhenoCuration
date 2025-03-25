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
  library(pROC)
})

RAVmodel <- getModel('C2', load=TRUE)


#Identifying sex predictor RAVs from sex-specific disease data sets 

tcga_brca <- curatedTCGAData(diseaseCode = 'BRCA',
                             assays = 'RNA*',
                             version = '2.0.1',
                             dry.run = FALSE)

tcga_ov <- curatedTCGAData(diseaseCode = 'OV',
                           assays = 'RNA*',
                           version = '2.0.1',
                           dry.run = FALSE)

tcga_ucec <- curatedTCGAData(diseaseCode = 'UCEC',
                             assays = 'RNA*',
                             version = '2.0.1',
                             dry.run = FALSE)

tcga_prad <- curatedTCGAData(diseaseCode = 'PRAD',
                             assays = 'RNA*',
                             version = '2.0.1',
                             dry.run = FALSE)

## Parse out cancer vs normal samples
tcga_brca <- TCGAsplitAssays(tcga_brca, c("01", "11"))
sampleTables(tcga_brca)

tcga_ov <- TCGAsplitAssays(tcga_ov, c("01"))
sampleTables(tcga_ov)

tcga_ucec <- TCGAsplitAssays(tcga_ucec, c("01"))
sampleTables(tcga_ucec)

tcga_prad <- TCGAsplitAssays(tcga_prad, c("01", "11"))
sampleTables(tcga_prad)

tcga_brca <- getWithColData(tcga_brca,
                            '01_BRCA_RNASeq2Gene-20160128')
assay(tcga_brca) <- log2(assay(tcga_brca) + 1)
meta_tcga_brca <- colData(tcga_brca)
tcga_brca_gender <- tcga_brca@colData@listData$patient.gender

tcga_ov <- getWithColData(tcga_ov,
                          '01_OV_RNASeq2Gene-20160128')
assay(tcga_ov) <- log2(assay(tcga_ov) + 1)
meta_tcga_ov <- colData(tcga_ov)
tcga_ov_gender <- tcga_ov@colData@listData$patient.gender

tcga_ucec <- getWithColData(tcga_ucec,
                            '01_UCEC_RNASeq2Gene-20160128')
assay(tcga_ucec) <- log2(assay(tcga_ucec) + 1)
meta_tcga_ucec <- colData(tcga_ucec)
tcga_ucec_gender <- tcga_ucec@colData@listData$patient.gender

tcga_prad <- getWithColData(tcga_prad,
                            '01_PRAD_RNASeq2Gene-20160128')
assay(tcga_prad) <- log2(assay(tcga_prad) + 1)
meta_tcga_prad <- colData(tcga_prad)
tcga_prad_gender <- tcga_prad@colData@listData$patient.gender

#validation & heatmap of top 15
validated_tcga_brca <- validate(tcga_brca, RAVmodel)
heatmapTable(validated_tcga_brca, RAVmodel, num.out=15, column_title = "BRCA Cancer")

validated_tcga_ov <- validate(tcga_ov, RAVmodel)
heatmapTable(validated_tcga_ov, RAVmodel, num.out=15, column_title = "OV Cancer")

validated_tcga_ucec <- validate(tcga_ucec, RAVmodel)
heatmapTable(validated_tcga_ucec, RAVmodel, num.out=15, column_title = "UCEC Cancer")

validated_tcga_prad <- validate(tcga_prad, RAVmodel)
heatmapTable(validated_tcga_prad, RAVmodel, num.out=15, column_title = "PRAD Cancer")

#overlapping RAVs between at least 3 sex-specific diseases
gender_RAV <- c("RAV99", "RAV119", "RAV135", "RAV222", "RAV312", 
                "RAV468", "RAV503", "RAV504", "RAV683", "RAV711", 
                "RAV868", "RAV1016", "RAV1303", "RAV1575", "RAV2515")

#BRCA sample scores for RF model test data
tcga_brca_sampleScores <- calculateScore(tcga_brca, RAVmodel) %>% as.data.frame()
tcga_brca_sampleScores$gender <- tcga_brca_gender

#building sex predictor RF model
tcga <- curatedTCGAData(diseaseCode = c('COAD', 'STAD', 'LUAD'),
                        assays = 'RNASeq2Gene',
                        version = '2.1.1',
                        dry.run = FALSE)
sampleTables(tcga)
tcga <- TCGAsplitAssays(tcga, c("01", "11"))

coad_cancer <- getWithColData(tcga,
                              '01_COAD_RNASeq2Gene-20160128')
stad_cancer <- getWithColData(tcga,
                              '01_STAD_RNASeq2Gene-20160128')
luad_cancer <- getWithColData(tcga,
                              '01_LUAD_RNASeq2Gene-20160128')

combined_data <- cbind(coad_cancer, stad_cancer, luad_cancer)
assay(combined_data) <- log2(assay(combined_data) + 1)
table(combined_data@colData@listData$patient.gender)

## coad, ucec, stad combined dataset
#combined_data <- readRDS("~/InferPhenoCuration/MSI_Status/data/combined_data.rds")
#table(combined_data@colData@listData$patient.gender)

combined_data_meta <- colData(combined_data)
combined_gender <- combined_data@colData@listData$patient.gender

#see what the top 10 validated RAVs are for this combined cancer dataset
combined_validate <- validate(combined_data, RAVmodel)
heatmapTable(combined_validate, RAVmodel, num.out=10)

topRAVs <- combined_validate %>% filter(score > 0.4) %>% rownames()

## Calculate sample scores
sampleScore <- calculateScore(combined_data, RAVmodel)

#validated_ind <- validatedSignatures(combined_validate, RAVmodel, num.out = 15, scoreCutoff = 0.45, indexOnly = TRUE) #Using Pearson Coefficient

## Subset sampleScore to sex predictor RAVs
sampleScore_gender <- sampleScore[, gender_RAV] %>% as.data.frame()
sampleScore_gender$gender <- combined_gender

## gender = binary variable - wilcoxon test

for (i in seq_len(ncol(sampleScore_gender)-1)) {
    print(colnames(sampleScore_gender)[i])
    ## wilcoxon test
    wilcox_test <- wilcox.test(sampleScore_gender[, i] ~ sampleScore_gender$gender, alternative="two.sided")
    print(wilcox_test$p.value)
}

# RAVs significant for patient gender based on wilcoxon p-values < 0.05
combined_genderRAV <- c("RAV504", "RAV868", "RAV1303", "RAV1575", "RAV2515")

sampleScore_gender_sub <- sampleScore_gender[,combined_genderRAV]
sampleScore_gender_sub$gender <- combined_gender
sampleScore_gender_sub$gender <- factor(sampleScore_gender_sub$gender)

## Split data between training data and test data
## Randomly select 70% of samples for training
set.seed(3)
num_sample <- nrow(sampleScore_gender_sub)

## 70%
train_sample_index <- sample(seq_len(num_sample), round(num_sample*0.7))
train_gender_data <- sampleScore_gender_sub[train_sample_index,]

## 30%
test_gender_data <- sampleScore_gender_sub[-train_sample_index,]

##RF model
rf_gender_model <- randomForest(gender ~ ., train_gender_data, proximity = TRUE)

##predict
rf_gender_pred <- predict(rf_gender_model, test_gender_data, type = "prob")

##ROC curves
roc_gender_rf <- roc(test_gender_data$gender, rf_gender_pred[, "female"])
## Plot ROC curves
plot(roc_gender_rf, col = "blue")

##feature importance
importance(rf_gender_model)
varImpPlot(rf_gender_model)

##confusion matrix
confusionMatrix(predict(rf_gender_model, test_gender_data), test_gender_data$gender)

#another test data - BRCA
brca_sampleScore_sub <- tcga_brca_sampleScores[,combined_genderRAV]
brca_sampleScore_sub$gender <- tcga_brca_gender
brca_sampleScore_sub$gender <- factor(brca_sampleScore_sub$gender)

rf_gender_pred2 <- predict(rf_gender_model, brca_sampleScore_sub, type = "prob")
confusionMatrix(predict(rf_gender_model, brca_sampleScore_sub), brca_sampleScore_sub$gender)
