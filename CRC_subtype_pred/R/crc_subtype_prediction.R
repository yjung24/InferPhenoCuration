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
silu_2022_exp <- read_delim("InferPhenoCuration/CRC_subtype_pred/data/silu_2022_data_mrna_seq_expression.txt", 
                            delim = "\t", escape_double = FALSE, 
                            trim_ws = TRUE)

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
silu_numeric <- apply(silu_numeric, 1, function(x) x - mean(x))
silu_numeric <- t(silu_numeric)

## storing gene names 
gene_names <- silu_2022_exp$Hugo_Symbol

common_genes <- intersect(rownames(exprs_df), gene_names)

## converting to matrix
silu_2022_exp2 <- as.matrix(silu_numeric)

## adding gene names as rownames
rownames(silu_2022_exp2) <- gene_names
silu_2022_exp3 <- silu_2022_exp2[common_genes,]
validate_silu <- validate(silu_2022_exp3, RAVmodel)
heatmapTable(validate_silu, RAVmodel, num.out = 10)

## calculate sample scores
silu_samplescore <- calculateScore(silu_2022_exp3, RAVmodel)
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

# Test model using rectal_msk_2022 data from cBioPortal (study used CMScaller)
## load clinical data and gene set 
msk_2022 <- read_tsv("~/InferPhenoCuration/CRC_subtype_pred/data/rectal_msk_2022_clinical_data.tsv")
msk_2022_exp <- read_delim("InferPhenoCuration/CRC_subtype_pred/data/msk_2022_data_mrna_seq_expression.txt", 
                           delim = "\t", escape_double = FALSE, 
                           trim_ws = TRUE)
table(msk_2022$CMS)

## recode CMS variable to harmonize with validation data
msk_2022$CMS <- as.factor(msk_2022$CMS)
msk_2022 <- msk_2022 %>% transform(CMS = recode_factor(CMS,"CMS1" = "CMS1", 
                                 "CMS2" = "CMS2",
                                 "CMS3" = "CMS3",
                                 "CMS4" = "CMS4",
                                 "1" = "NA")) %>%
                        mutate(cms_label = CMS)

## fix duplicate gene/row names  
msk_2022_gene_names <- make.unique(msk_2022_exp$Hugo_Symbol)

## only some samples had expression data available: store these sample names
msk_2022_sample_id <- intersect(msk_2022$`Sample ID`, colnames(msk_2022_exp))

## subset samples in both meta and exp datasets and center exp data
msk_2022_exp_sub <- msk_2022_exp[,msk_2022_sample_id]
msk_2022_exp_sub <- as.matrix(msk_2022_exp_sub)
msk_2022_exp_sub <- apply(msk_2022_exp_sub, 1, function(x) x - mean(x))
msk_2022_exp_sub <- t(msk_2022_exp_sub)

## assign rownames to exp data subset
rownames(msk_2022_exp_sub) <- msk_2022_gene_names

## validate using RAVmodel
validate_msk <- validate(msk_2022_exp_sub, RAVmodel)
heatmapTable(validate_msk, RAVmodel, num.out = 10)

## obtain sample scores
msk_samplescore <- calculateScore(msk_2022_exp_sub, RAVmodel)
msk_rav_subset <- as.data.frame(msk_samplescore[,validated_crc_RAVs])


msk_rav_subset <- msk_rav_subset %>% mutate('Sample ID' = rownames(msk_rav_subset))
## subset meta data for samples with exp data
msk_2022_2 <- as.data.frame(msk_2022[msk_2022$'Sample ID' %in% msk_2022_sample_id,])
## combine meta data and exp data using merge on Sample ID column
msk_combined <- merge(msk_2022_2, msk_rav_subset, by = 'Sample ID')
## remove NAs for CMS
msk_combined <- msk_combined %>% 
  filter(cms_label != "NA") 

## drop NA factor level
msk_combined$cms_label <- factor(msk_combined$cms_label)

## prediction model for rectal_msk_2022
rf_cms_pred3 <- predict(rf_cms_model, msk_combined)
confusionMatrix(rf_cms_pred3, msk_combined$cms_label)

## GSEA
RAVmodel2 <- RAVmodel[,validated_crc_RAVs]

## currently hard coded to run with RAVmodel2 with C_ gene set; need to convert to function
gsea_script <- "~/InferPhenoCuration/CRC_subtype_pred/R/msigdb_gmt.R"
source(gsea_script)  

gsea_dir <- file.path("~/InferPhenoCuration/CRC_subtype_pred", paste0("gsea"))  # GSEA C_ DB is saved here

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

val_RAV_numeric <- c("188","1575","834","832","833","438","324","192","981","220")
  for (i in seq_len(val_RAV_numeric)) {
print(annotateRAV(RAVmodel2, val_RAV_numeric[i]))
  }
