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
#sampleTables(tcga_brca)

tcga_ov <- TCGAsplitAssays(tcga_ov, c("01"))
#sampleTables(tcga_ov)

tcga_ucec <- TCGAsplitAssays(tcga_ucec, c("01"))
#sampleTables(tcga_ucec)

tcga_prad <- TCGAsplitAssays(tcga_prad, c("01", "11"))
#sampleTables(tcga_prad)

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
tcga_rf <- curatedTCGAData(diseaseCode = c('COAD', 'STAD', 'LUAD'),
                        assays = 'RNASeq2Gene',
                        version = '2.1.1',
                        dry.run = FALSE)
sampleTables(tcga_rf)
tcga_rf <- TCGAsplitAssays(tcga_rf, c("01", "11"))

coad_cancer <- getWithColData(tcga_rf,
                              '01_COAD_RNASeq2Gene-20160128')
stad_cancer <- getWithColData(tcga_rf,
                              '01_STAD_RNASeq2Gene-20160128')
luad_cancer <- getWithColData(tcga_rf,
                              '01_LUAD_RNASeq2Gene-20160128')

combined_data_rf <- cbind(coad_cancer, stad_cancer, luad_cancer)
assay(combined_data_rf) <- log2(assay(combined_data_rf) + 1)
table(combined_data_rf@colData@listData$patient.gender)

combined_data_meta <- colData(combined_data_rf)
combined_gender <- combined_data_rf@colData@listData$patient.gender

#see what the top 10 validated RAVs are for this combined cancer dataset
combined_validate <- validate(combined_data_rf, RAVmodel)
heatmapTable(combined_validate, RAVmodel, num.out=10)

topRAVs <- combined_validate %>% filter(score > 0.4) %>% rownames()

## Calculate sample scores
sampleScore <- calculateScore(combined_data_rf, RAVmodel)

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

#repeat with RAVs validated using combined 4 sex-specific cancer TCGA data sets 
tcga <- curatedTCGAData(diseaseCode = c('BRCA', 'UCEC', 'OV', 'PRAD'),
                        assays = 'RNASeq2Gene',
                        version = '2.1.1',
                        dry.run = FALSE)
sampleTables(tcga)
tcga <- TCGAsplitAssays(tcga, c("01", "11"))
br_cancer <- getWithColData(tcga,
                              '01_BRCA_RNASeq2Gene-20160128')
ut_cancer <- getWithColData(tcga,
                              '01_UCEC_RNASeq2Gene-20160128')
ov_cancer <- getWithColData(tcga,
                              '01_OV_RNASeq2Gene-20160128')
pr_cancer <- getWithColData(tcga,
                            '01_PRAD_RNASeq2Gene-20160128')

sex_tcga <- cbind(br_cancer, ut_cancer, ov_cancer, pr_cancer)
assay(sex_tcga) <- log2(assay(sex_tcga) + 1)
table(sex_tcga@colData@listData$patient.gender)

sex_tcga_meta <- colData(sex_tcga)
sex_tcga_gender <- sex_tcga@colData@listData$patient.gender

##combining validation data sets
combined_data <- cbind(coad_cancer, stad_cancer, luad_cancer)
assay(combined_data) <- log2(assay(combined_data) + 1)
table(combined_data@colData@listData$patient.gender)

sex_tcga_val <- validate(sex_tcga, RAVmodel)
heatmapTable(sex_tcga_val, RAVmodel, num.out=15, column_title = "Sex Specific Cancer TCGA")
sex_tcga_samplescore <- calculateScore(sex_tcga, RAVmodel)
sex_tcga_ravs <- c("RAV99", "RAV220", "RAV221", "RAV222", "RAV468", "RAV503",  
                   "RAV504", "RAV532", "RAV625", "RAV683", "RAV868", "RAV1241",
                   "RAV1303", "RAV1575")
sex_tcga_rav_subset <- sex_tcga_samplescore[,sex_tcga_ravs] %>% as.data.frame()
sex_tcga_rav_subset$gender <- sex_tcga_gender %>% as.factor()

for (i in seq_len(ncol(sex_tcga_rav_subset)-1)) {
  print(colnames(sex_tcga_rav_subset)[i])
  ## wilcoxon test
  wilcox_test <- wilcox.test(sex_tcga_rav_subset[, i] ~ sex_tcga_rav_subset$gender, 
                             alternative="two.sided")
  print(wilcox_test$p.value)
}

## excluding RAVs with Wilcoxon p-values > 0.05 and pearson corr. coeff. > 0.50
sex_tcga_ravs <- c("RAV99", "RAV220", "RAV221", "RAV468", "RAV503",  
                   "RAV504", "RAV532", "RAV683", "RAV1241")

## scatter plot for a pair of RAVs
df_x <- data.frame(RAV99 = sex_tcga_rav_subset$RAV99, 
                   RAV220 = sex_tcga_rav_subset$RAV220,
                   RAV221 = sex_tcga_rav_subset$RAV221,
                   RAV468 = sex_tcga_rav_subset$RAV468,
                   RAV503 = sex_tcga_rav_subset$RAV503,
                   RAV504 = sex_tcga_rav_subset$RAV504,
                   RAV532 = sex_tcga_rav_subset$RAV532,
                   RAV683 = sex_tcga_rav_subset$RAV683,
                   RAV1241 = sex_tcga_rav_subset$RAV1241)

df_x <- data.frame(RAV99 = sex_tcga_rav_subset$RAV99, 
                   RAV220 = sex_tcga_rav_subset$RAV220,
                   RAV221 = sex_tcga_rav_subset$RAV221,
                   RAV468 = sex_tcga_rav_subset$RAV468,
                   RAV503 = sex_tcga_rav_subset$RAV503,
                   RAV504 = sex_tcga_rav_subset$RAV504,
                   RAV532 = sex_tcga_rav_subset$RAV532,
                   RAV683 = sex_tcga_rav_subset$RAV683,
                   RAV1241 = sex_tcga_rav_subset$RAV1241)

combos <- expand.grid(xvar = names(df_x), yvar = names(df_y), stringsAsFactors = FALSE)
plot_data <- pmap_dfr(combos, function(xvar, yvar) {
  tibble(
    x = df_x[[xvar]],
    y = df_y[[yvar]],
    group = sex_tcga_rav_subset$gender,
    xvar = xvar,
    yvar = yvar
  )
})
ggplot(plot_data, aes(x = x, y = y, color = group)) +
  geom_point(size = 1) +
  facet_grid(rows = vars(yvar), cols = vars(xvar), scales = "free") +
  theme_minimal() +
  labs(title = "Scatter Plot of all Permutations of RAVs Colored by Gender")

## rf model using non-sex specific TCGA datasets predicted by sex-specific RAVs
## Subset sampleScore to sex predictor RAVs
sampleScore_sex_tcga <- sampleScore[, sex_tcga_ravs] %>% as.data.frame()
sampleScore_sex_tcga$gender <- combined_gender %>% as.factor()

for (i in seq_len(ncol(sampleScore_sex_tcga)-1)) {
  print(colnames(sampleScore_sex_tcga)[i])
  ## wilcoxon test
  wilcox_test <- wilcox.test(sampleScore_sex_tcga[, i] ~ sampleScore_sex_tcga$gender, 
                             alternative="two.sided")
  print(wilcox_test$p.value)
}

## only RAV220 and RAV504 have Wilcoxon p-values < 0.05 using test dataset

## 1st RF model will be with samples scores of all sex_tcga RAVs
## Randomly select 70% of samples for training
set.seed(5)
num_sample <- nrow(sampleScore_sex_tcga)

## 70%
train_sample_index <- sample(seq_len(num_sample), round(num_sample*0.7))
train_sex_tcga <- sampleScore_sex_tcga[train_sample_index,]

## 30%
test_sex_tcga <- sampleScore_sex_tcga[-train_sample_index,]

##RF model
rf_sex_tcga_model <- randomForest(gender ~ ., train_sex_tcga, proximity = TRUE)

##predict
rf_sex_tcga_pred <- predict(rf_sex_tcga_model, test_sex_tcga, type = "prob")

##ROC curves
roc_rf_sex_tcga <- roc(test_sex_tcga$gender, rf_sex_tcga_pred[, "female"])
## Plot ROC curves
plot(roc_rf_sex_tcga, col = "blue")

##feature importance
importance(rf_sex_tcga_model)
varImpPlot(rf_sex_tcga_model)

##confusion matrix
confusionMatrix(predict(rf_sex_tcga_model, test_sex_tcga), test_sex_tcga$gender)

## 2nd RF model will be with samples scores of RAVs 220 and 504

##RF model
rf_sex_tcga_model2 <- randomForest(gender ~ RAV220 + RAV504, train_sex_tcga, proximity = TRUE)

##predict
rf_sex_tcga_pred2 <- predict(rf_sex_tcga_model2, test_sex_tcga, type = "prob")

##ROC curves
roc_rf_sex_tcga2 <- roc(test_sex_tcga$gender, rf_sex_tcga_pred2[, "female"])
## Plot ROC curves
plot(roc_rf_sex_tcga2, col = "blue")

##feature importance
importance(rf_sex_tcga_model2)
varImpPlot(rf_sex_tcga_model2)

##confusion matrix
confusionMatrix(predict(rf_sex_tcga_model2, test_sex_tcga), test_sex_tcga$gender)


