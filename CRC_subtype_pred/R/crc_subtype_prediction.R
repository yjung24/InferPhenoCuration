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

## Load data (Results/CRC/data/eSets/setNames.RData)
load("~/InferPhenoCuration/setNames.RData")   # 18 CRC datasets
names(validated_ind_all) <- setNames

## validated_RAVs 834 and 833 - sample scores and CMS labels already available in following dataset
test <- readRDS("~/InferPhenoCuration/CRC_subtype_pred/data/SummaryForFig4.rds")

## convert CMS labels into factor
test<- as.factor(test$cms_label_crc)

# Data preparation for building RF model

## Split data between training data and test data
## Randomly select 70% of samples for training
target_attr <- "cms_label_crc"
set.seed(2)
num_sample <- nrow(test)

## random sampling of row numbers
train_sample_index <- sample(seq_len(num_sample), round(num_sample*0.7))
train_cms_data <- test[train_sample_index,]

## remaining 30% of data will be for validation of RF model
test_cms_data <- test[-train_sample_index,]

## random forest model
rf_cms_model <- randomForest(cms_label_RF ~ RAV834 + RAV833, train_cms_data, ntree = 500, 
                       keep.forest = TRUE, importance = TRUE)

## RF prediction
rf_cms_pred <- predict(rf_cms_model, test_cms_data)

conMatrix <- confusionMatrix(rf_cms_pred, test_cms_data$cms_label_RF)
conMatrix





