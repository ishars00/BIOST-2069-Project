######## Libraries ########
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# for GO pathway analysis
# BiocManager::install("org.Hs.eg.db")

library(clusterProfiler)
library("org.Hs.eg.db", character.only = TRUE)
library(tidyverse)
library(here)




### Read in full DESeq2 analysis for human data
hmn_res = read_csv(here('Results', 'human-deseq2-results.csv'))

## Top DE genes were determined by pvalue < 0.05
hmn_de = hmn_res |> filter(pvalue < 0.05)

## Create a GeneList named vec structure for GO GSE analysis
genes = hmn_res |>
  arrange(desc(log2FoldChange))

geneList = genes |> pull(log2FoldChange)
names(geneList) = genes |> pull(gene)

# load db 
ego <- gseGO(
  geneList = geneList,
  OrgDb = org.Hs.eg.db,
  ont = "ALL",
  keyType = "SYMBOL",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05, 
  verbose = TRUE,
  eps = 0,
)

# plot results...?
dotplot(ego, showCategory=7, split=".sign") + facet_grid(.~.sign)
