## Packages -------------
suppressPackageStartupMessages({
  # Classifiers
  library(CMSclassifier)
  library(Biobase)
  library(limma)
  library(CMScaller)
  
  # Gene Annotations
  library(org.Hs.eg.db)
  
  # CRAN
  library(tidyverse) # includes dplyr, ggplot2, magrittr, tidyr
  library(magick)
  library(ztable)
  library(metafolio)
  library(randomForest)
  library(readr)
  library(dplyr)
  library(caret)
})

# CMSclassifier on msk_2022 - originally used CMScaller

## Load data
#msk_2022 <- read_tsv("~/InferPhenoCuration/CRC_subtype_pred/data/rectal_msk_2022_clinical_data.tsv")
#msk_2022_exp <- read_delim("InferPhenoCuration/CRC_subtype_pred/data/msk_2022_data_mrna_seq_expression.txt", 
#                           delim = "\t", escape_double = FALSE, 
#                           trim_ws = TRUE)

## export as RDS objects to make this code reproducible

saveRDS(msk_2022, file = "~/InferPhenoCuration/CRC_subtype_pred/data/msk_2022.rds")
saveRDS(msk_2022_exp, file = "~/InferPhenoCuration/CRC_subtype_pred/data/msk_2022_exp.rds")

## Account for duplicate Hugo Symbols
msk_2022_gep <- msk_2022_exp %>% 
  distinct(Hugo_Symbol, .keep_all = TRUE) %>%
  as.data.frame()

### only 113 samples have expression data available
msk_2022_gep_samples <- colnames(msk_2022_gep)

## convert hugo symbols to entrez id probes
entrezids <- mapIds(org.Hs.eg.db,
                                keys = msk_2022_gep$Hugo_Symbol,
                                column = "ENTREZID",
                                keytype = "SYMBOL",
                                multiVals = "first")
## add probes as rownames
msk_2022_gep$probe <- entrezids
msk_2022_gep_filtered <- msk_2022_gep[!is.na(msk_2022_gep$probe), ]
rownames(msk_2022_gep_filtered) <- msk_2022_gep_filtered$probe

## remove hugo symbol and probe columns
msk_2022_gep_filtered$Hugo_Symbol <- NULL
msk_2022_gep_filtered$probe <- NULL

## conduct log transformation
#msk_2022_gep_log <- log2(msk_2022_gep_filtered + 1) <- determined exprs data was already log-transformed

## call CMSclassifier
predictCMS_msk_2022 <- classifyCMS(msk_2022_gep_filtered,"RF")

## call CMScaller
predictCMS_msk_2022_CMScaller <- CMScaller(msk_2022_gep_filtered)

## select samples from meta data that have exp data available
msk_2022_meta_sub <- msk_2022[msk_2022$`Sample ID` %in% msk_2022_gep_samples,]
benchmark <- dplyr::bind_cols(CMSclassifier = predictCMS_msk_2022$predictedCMS$RF,
                              original_CMS_preds = msk_2022_meta_sub$CMS)
benchmark_CMScaller <- dplyr::bind_cols(CMScaller = predictCMS_msk_2022_CMScaller$prediction,
                                        original_CMS_preds = msk_2022_meta_sub$CMS)
table(benchmark_CMScaller$CMScaller, benchmark_CMScaller$original_CMS_preds, useNA = "always")

## recode CMS variable to harmonize with validation data and filter CMS == NA 
benchmark <- benchmark %>% 
  transform(original_CMS_preds = 
              recode_factor(original_CMS_preds,"CMS1" = "CMS1", 
                            "CMS2" = "CMS2",
                            "CMS3" = "CMS3",
                            "CMS4" = "CMS4",
                            "1" = "NA")) %>%
  filter(original_CMS_preds != "NA")

benchmark$CMSclassifier <- as.factor(benchmark$CMSclassifier)
confusionMatrix(benchmark$CMSclassifier, benchmark$original_CMS_preds)

benchmark_CMScaller2 <- benchmark_CMScaller %>% 
  transform(original_CMS_preds = 
              recode_factor(original_CMS_preds,"CMS1" = "CMS1", 
                            "CMS2" = "CMS2",
                            "CMS3" = "CMS3",
                            "CMS4" = "CMS4",
                            "1" = "NA")) %>%
  filter(original_CMS_preds != "NA")

benchmark_CMScaller2$CMScaller <- as.factor(benchmark_CMScaller2$CMScaller)
confusionMatrix(benchmark_CMScaller2$CMScaller, benchmark_CMScaller2$original_CMS_preds)

#CMScaller on silu_2022

silu_2022 <- read_tsv("~/InferPhenoCuration/CRC_subtype_pred/data/coad_silu_2022_clinical_data.tsv")
silu_2022_exp <- read_delim("InferPhenoCuration/CRC_subtype_pred/data/silu_2022_data_mrna_seq_expression.txt", 
                            delim = "\t", escape_double = FALSE, 
                            trim_ws = TRUE)

saveRDS(silu_2022, file = "~/InferPhenoCuration/CRC_subtype_pred/data/silu_2022.rds")
saveRDS(silu_2022_exp, file = "~/InferPhenoCuration/CRC_subtype_pred/data/silu_2022_exp.rds")


silu_2022_gep <- silu_2022_exp %>% 
  distinct(Hugo_Symbol, .keep_all = TRUE) %>%
  as.data.frame()

entrezids2 <- mapIds(org.Hs.eg.db,
                    keys = silu_2022_gep$Hugo_Symbol,
                    column = "ENTREZID",
                    keytype = "SYMBOL",
                    multiVals = "first")

silu_2022_gep$probe <- entrezids2
silu_2022_gep_filtered <- silu_2022_gep[!is.na(silu_2022_gep$probe), ]
rownames(silu_2022_gep_filtered) <- silu_2022_gep_filtered$probe

silu_2022_gep_filtered$Hugo_Symbol <- NULL
silu_2022_gep_filtered$probe <- NULL

# silu_2022_gep_log <- log2(silu_2022_gep_filtered + 1) <- determined that data was already log-transformed

#table(rownames(silu_2022_gep_log) %in% templates.CMS$probe, useNA="always")
      
predictCMS_silu_2022 <- CMScaller(silu_2022_gep_filtered)

silu_2022_meta_sub <- silu_2022[silu_2022$`Sample ID` %in% colnames(silu_2022_gep_filtered),]

benchmark2 <- dplyr::bind_cols(CMScaller = predictCMS_silu_2022$prediction,
                              original_CMS_preds = silu_2022_meta_sub$CMS)
table(benchmark2$CMScaller, benchmark2$original_CMS_preds, useNA = "always")
benchmark2 <- benchmark2 %>% 
  transform(original_CMS_preds = 
              recode_factor(original_CMS_preds,"CMS1" = "CMS1", 
                            "CMS2" = "CMS2",
                            "CMS3" = "CMS3",
                            "CMS4" = "CMS4",
                            "mixed" = "NA")) %>%
  filter(original_CMS_preds != "NA") 

benchmark2$original_CMS_preds <- droplevels(benchmark2$original_CMS_preds)

benchmark2$CMScaller <- as.factor(benchmark2$CMScaller)

confusionMatrix(benchmark2$CMScaller, benchmark2$original_CMS_preds)


