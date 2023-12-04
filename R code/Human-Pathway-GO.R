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
gse <- gseGO(
  geneList = geneList,
  OrgDb = org.Hs.eg.db,
  ont = "ALL",
  keyType = "SYMBOL",
  pAdjustMethod = "none", # or BH?
  pvalueCutoff = 0.05, 
  verbose = TRUE,
  eps = 0,
)

# plot results...?
dotplot(gse, showCategory=7, split=".sign") + facet_grid(.~.sign)

# do basic enrichment
ego = enrichGO(
  gene = hmn_de |> pull(gene),
  universe = hmn_res |> pull(gene),
  OrgDb = org.Hs.eg.db,
  ont = "ALL",
  keyType = "SYMBOL",
  pvalueCutoff  = 0.05,
  pAdjustMethod = 'BH',
)

dotplot(ego)

# save out enriched pathways for REVIGO analysis
ego@result |>
  as_tibble() |>
  select(ID, pvalue) |>
  write_delim(here('Results', 'Human-Enriched-Pathways.txt'), delim='\t')
