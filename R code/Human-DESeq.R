
######## Libraries ########
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
# BiocManager::install(version = "3.18")

## For fetching data from GEO
BiocManager::install("GEOquery")


library(GEOquery)
library(tidyverse)
library(here)
library(DESeq2)
library(EnhancedVolcano)

# helpers

to_named <- function(df, converter, rowcol) {
  res = converter(df |> select(!{{rowcol}}))
  rownames(res) = df |> pull({{rowcol}})
  return (res)
}

# First, get the data from GEO using BioConductor
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE167536
# geo_accension = 'GSE167536'
# gds = getGEO(geo_accension, GSEMatrix=FALSE)


hmn_expr = read_csv(here('GSE167536', 'GSE167536_humanstress_rawcounts.txt')) |>
  to_named(as.matrix, 'gene')

hmn_samples = read_csv(here('GSE167536', 'BioSample', 'human_sample_info.csv')) |>
  mutate(group = relevel(as.factor(str_replace(group, ' ', '_')), 'low_stress')) |>
  dplyr::rename(provider = `biomaterial provider`) |>
  to_named(as.data.frame, 'SampleName') |>
  (function(df, o) df[o, ])(colnames(hmn_expr)) |>
  AnnotatedDataFrame()

hmn_eset = ExpressionSet(hmn_expr, hmn_samples)
hmn_se = makeSummarizedExperimentFromExpressionSet(hmn_eset)
hmn_dss = DESeqDataSet(hmn_se, ~ group)

# do not filter the genes?

# levels set when reading data
hmn_dss = DESeq(hmn_dss)

# they used pvalue < 0.05 for human DE genes
hmn_res = results(hmn_dss, contrast = c('group', 'low_stress', 'high_stress'), alpha=0.05)
hmn_res = hmn_res[complete.cases(hmn_res),]
EnhancedVolcano(hmn_res, lab=rownames(hmn_res), x='log2FoldChange', y='pvalue',
                pCutoff = 0.05)

# results to match with the paper's Sup Fig 7
paper_res = hmn_res[hmn_res$pvalue < 0.05,]



