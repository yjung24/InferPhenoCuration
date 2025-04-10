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
  library(readr)
  library(dplyr)
  library(caret)
})


## Load RAVmodel -------
RAVmodel <- getModel('C2', load=TRUE)

load("~/InferPhenoCuration/CRC_subtype_pred/data/eSets/setNames.RData")
setNames

## Load validation samples
for (set in setNames) {
  load(paste0("~/InferPhenoCuration/CRC_subtype_pred/data/eSets/", set, '.RData'))
}

## final dataframe = phenotype/meta data and sample scores
## phenotype tables combined
pdata_df <- setNames %>% lapply(function(set) {
  eSet <- get(set)
  pdata <- pData(eSet)
  eSet_tmp <- eSet[, pdata$sample_type %in% "tumor"]
  pdata_tmp <- pData(eSet_tmp) #tumor phenotype
  
  ind_rm <- grep("CRIS_", colnames(pdata_tmp))
  if (length(ind_rm) != 0) {pdata_tmp <- pdata_tmp[,-ind_rm]} #remove "CRIS" column 
  pdata_tmp$study <- set  # add 'study' column
  
  return(pdata_tmp) #returns phenotype data for each CRC dataset
}) %>% Reduce('rbind', .)

## common genes between all validation datasets
all_genes <- list()
for (set in setNames) {
  eSet <- get(set)
  exprs <- exprs(eSet) %>% rmNaInf #Remove rows with missing and Inf values from a matrix
  all_genes[[set]] <- rownames(exprs)
}
cg <- Reduce(intersect, all_genes)

## expression matrix combined <- dataframe for applying sample scores to
exprs_df <- setNames %>% lapply(function(set) {
  eSet <- get(set)
  pdata <- pData(eSet)
  eSet_tmp <- eSet[cg, pdata$sample_type %in% "tumor"]
  exprs_tmp <- exprs(eSet_tmp) %>% rmNaInf
  exprs_tmp <- apply(exprs_tmp, 1, function(x) x - mean(x)) %>% t
  return(exprs_tmp)  
}) %>% Reduce('cbind', .) 

# calculate sample scores
validate_crc <- validate(exprs_df, RAVmodel)
heatmapTable(validate_crc, RAVmodel, num.out = 10)
validated_crc_RAVs <- c("RAV188", "RAV1575", "RAV834", "RAV832", "RAV833", 
                        "RAV438", "RAV324", "RAV192", "RAV981", "RAV220")
sampleScore <- calculateScore(exprs_df, RAVmodel)
sampleScore_sub <- sampleScore[,validated_crc_RAVs]
data_all <- cbind(sampleScore_sub, pdata_df)

## CMS = categorical variable - kruskal-wallis

kruskal_p.value <- list()

for (i in seq_len(ncol(sampleScore_sub))) {
  RAV <- colnames(sampleScore_sub)[i]
  ## kruskal-wallis
  kruskal <- kruskal.test(sampleScore_sub[, i] ~ data_all$cms_label_crc)
  p <- kruskal$p.value
  kruskal_p.value[[RAV]] <- p
}

kruskal_p.value

## convert CMS labels into factor
data_all$cms_label_crc <- as.factor(data_all$cms_label_crc)
## transform cms_label_crc
data_all_clean <- data_all %>% 
  transform(cms_label_crc = 
              recode_factor(cms_label_crc, 
                            "CMS1" = "CMS1",
                            "CMS2" = "CMS2",
                            "CMS3" = "CMS3",
                            "CMS4" = "CMS4",
                            "NOLBL" = "NA")) %>%
  mutate(cms_label = cms_label_crc) %>%
  filter(cms_label != "NA") 

data_all_clean$cms_label <- factor(data_all_clean$cms_label)

# Data preparation for building RF model

## Split data between training data and test data
## Randomly select 70% of samples for training
target_attr <- "cms_label"
set.seed(2)
num_sample <- nrow(data_all_clean)

## random sampling of row numbers
train_sample_index <- sample(seq_len(num_sample), round(num_sample*0.7))
train_cms_data <- data_all_clean[train_sample_index,]

## remaining 30% of data will be for validation of RF model
test_cms_data <- data_all_clean[-train_sample_index,]

## random forest model
rf_cms_model <- randomForest(cms_label ~ RAV188 + RAV1575 + RAV834 + RAV832 + 
                              RAV833 + RAV438 + RAV324 + RAV192 + RAV981 + 
                              RAV220, train_cms_data, ntree = 500, 
                              keep.forest = TRUE, importance = TRUE)

## RF prediction using remainder of validation data
rf_cms_pred <- predict(rf_cms_model, test_cms_data)

conMatrix <- confusionMatrix(rf_cms_pred, test_cms_data$cms_label)
conMatrix

# Test model using silu_2022 data 
## load clinical data and gene set 
silu_2022 <- read_tsv("~/InferPhenoCuration/CRC_subtype_pred/data/coad_silu_2022_clinical_data.tsv")
silu_2022_exp <- read.delim("~/InferPhenoCuration/CRC_subtype_pred/data/data_mrna_seq_expression.txt")

silu_2022$cms <- factor(silu_2022$CMS)


## harmonize CMS label between silu_2022 and validation dataset
silu_2022_tf <- silu_2022 %>% 
  transform(CMS = recode_factor(CMS,"CMS1" = "CMS1", 
                                    "CMS2" = "CMS2",
                                    "CMS3" = "CMS3",
                                    "CMS4" = "CMS4",
                                    "mixed" = "NA")) %>%
  mutate(cms_label = CMS)

## convert expression df to matrix
## filtering for numeric columns/removing gene symbol column before converting to matrix
silu_numeric <- silu_2022_exp[, sapply(silu_2022_exp, is.numeric)]

## storing gene names 
gene_names <- silu_2022_exp$Hugo_Symbol

## converting to matrix
silu_2022_exp2 <- as.matrix(silu_numeric)

## adding gene names as rownames
rownames(silu_2022_exp2) <- gene_names
validate_silu <- validate(silu_2022_exp2, RAVmodel)
heatmapTable(validate_silu, RAVmodel, num.out = 10)

## calculate sample scores
silu_samplescore <- calculateScore(silu_2022_exp2, RAVmodel)
## subset sample scores for predictor RAVs only as a df
silu_rav_subset <- as.data.frame(silu_samplescore[,validated_crc_RAVs])
## combine sample scores and meta data
silu_combined <- cbind(silu_rav_subset, cms_label = silu_2022_tf$cms_label)

## drop na values in cms_label column
silu_combined <- silu_combined %>% 
  filter(cms_label != "NA") 
silu_combined$cms_label <- factor(silu_combined$cms_label)

## predict cms subtype label for silu_2022 dataset
rf_cms_pred2 <-predict(rf_cms_model, silu_combined)
confusionMatrix(rf_cms_pred2, silu_combined$cms_label)

## GSEA
RAVmodel2 <- RAVmodel[,validated_crc_RAVs]

## currently hard coded to run with RAVmodel2 with C5 gene set; need to convert to function
gsea_script <- "~/InferPhenoCuration/CRC_subtype_pred/R/msigdb c5.R"
source(gsea_script)  

gsea_dir <- file.path("~/InferPhenoCuration/CRC_subtype_pred", paste0("gsea"))  # GSEA C5 DB is saved here

## from source script in GenomicSuperSignature 
## code edit: "qvalues" -> "qvalue"
searchPathways_edit <- function(RAVmodel, gsea.dir) {
  
  ## If you want to select only a subset of RAVs with the specific cluster size
  # ind <- which(metadata(RAVmodel)$size > 3)
  # gsea_all <- vector(mode = "list", length = length(ind))
  # names(gsea_all) <- colnames(RAVmodel)[ind]
  
  gsea_all <- vector(mode = "list", length = ncol(RAVmodel))
  names(gsea_all) <- colnames(RAVmodel)
  gsea.dir <- gsea.dir
  
  for (i in seq_len(ncol(RAVmodel))) {
    pathToRes <- file.path(gsea.dir, paste0("gsea_", i, ".rds"))
    res <- readRDS(pathToRes)
    
    # If there is no enriched pathways
    if (nrow(res) == 0) {
      resName <- paste0("RAV", i)
      gsea_all[[resName]] <- res[, c("Description", "NES", "pvalue", "qvalue"), drop = FALSE]
      print(paste("RAV", i, "has no enriched pathways."))
      next
    }
    
    res <- res[which(res$qvalue == min(res$qvalue)), c("Description", "NES", "pvalue", "qvalue"), drop = FALSE]
    resName <- colnames(RAVmodel)[i]#paste0("RAV", i)
    gsea_all[[resName]] <- res
    print(paste("RAV", i, "is added."))
  }
  
  return(gsea_all)
}

gsea_all <- searchPathways_edit(RAVmodel2, gsea_dir)  
gsea(RAVmodel2) <- gsea_all

annotateRAV(RAVmodel2, 192)
