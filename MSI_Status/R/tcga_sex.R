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

RAVmodel <- getModel('C2', load=TRUE)

## load data all studies - took too long did not finish running
## tcga_data <- curatedTCGAData(assays = 'RNASeq2Gene',
##                        version = '2.1.1',
##                        dry.run = FALSE)

## coad, ucec, stad combined dataset
combined_data <- readRDS("~/InferPhenoCuration/MSI_Status/data/combined_data.rds")
table(combined_data@colData@listData$patient.gender)

combined_data_meta <- colData(combined_data)

combined_validate <- validate(combined_data, RAVmodel)
heatmapTable(combined_validate, RAVmodel, num.out=10)

## Select columns with >10% completeness
keep_attribute_ind <- which(colSums(!is.na(combined_data_meta)) > round(nrow(combined_data_meta)/10))
dat <- combined_data_meta[keep_attribute_ind]
dat <- subset(combined_data_meta, select= -patientID)

## Calculate sample scores
sampleScore <- calculateScore(combined_data, RAVmodel)

validated_ind <- validatedSignatures(combined_validate, RAVmodel, num.out = 15, scoreCutoff = 0.45, indexOnly = TRUE) #Using Pearson Coefficient

## Subset sampleScore to join with MCPcounter
sampleScore_sub <- sampleScore[, validated_ind] %>% as.data.frame()

## convert character type variables to factors
factorTb <- dat[, sapply(dat, class) == 'character']
factorTb[sapply(factorTb, is.character)] <- lapply(factorTb[sapply(factorTb, is.character)], factor, exclude = NULL)

## gender = binary variable - wilcoxon test
for (i in seq_len(ncol(sampleScore_sub))) {
  print(colnames(sampleScore_sub)[i])
  xy <- pairwise.wilcox.test(sampleScore_sub[, i], factorTb[,"patient.gender"], p.adjust.method = "bonferroni")
  print(xy)
}

wilcox_test_res <- as.data.frame(matrix(nrow = ncol(factorTb),
                                        ncol = ncol(sampleScore_sub)))

rownames(wilcox_test_res) <- colnames(factorTb)
colnames(wilcox_test_res) <- colnames(sampleScore_sub)

wtest_coad_wvalue <- wilcox_test_res
wtest_coad_pvalue <- wilcox_test_res

for (i in seq_len(ncol(sampleScore_sub))) {
    ## wilcoxon test
    wilcox_test <- wilcox.test(sampleScore_sub[, i] ~ factorTb[,"patient.gender"], alternative="two.sided")
    
    ## W value
    wval <- wilcox_test$statistic
    wtest_coad_wvalue["patient.gender", i] <- wval
    
    ## p-value
    pval <- wilcox_test$p.value
    wtest_coad_pvalue["patient.gender", i] <- pval
    
}
wtest_coad_wvalue <- wtest_coad_wvalue["patient.gender",]
wtest_coad_pvalue <- wtest_coad_pvalue["patient.gender",]

# RAVs significant for patient gender based on wilcoxon p-values < 0.05
genderRAV <- c("RAV504", "RAV468", "RAV220", "RAV1575", "RAV1303", "RAV119", "RAV324", "RAV1243",
               "RAV312", "RAV99", "RAV438")
