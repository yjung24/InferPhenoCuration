## Packages -------------
suppressPackageStartupMessages({
  library(GenomicSuperSignature)
  library(clusterProfiler)
})

# MSigDB C5

dat_dir <- "~/InferPhenoCuration/CRC_subtype_pred"

term2gene <- clusterProfiler::read.gmt("~/InferPhenoCuration/CRC_subtype_pred/data/c5.all.v2024.1.Hs.symbols.gmt")
colnames(term2gene) <- c("gs_name", "entrez_gene")

for (i in seq_len(ncol(RAVmodel2))) {
  fname <- paste0("gsea_", i, ".rds")
  fpath <- file.path(dat_dir, "gsea", fname)
  
  geneList <- RAVindex(RAVmodel2)[,i]
  geneList <- sort(geneList, decreasing = TRUE)
  res <- clusterProfiler::GSEA(geneList, TERM2GENE = term2gene,
                               pvalueCutoff = 0.05, seed = TRUE)
  saveRDS(res, fpath)
}