
######## Libraries ########
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

library(tidyverse)
library(here)
library(DESeq2)
library(EnhancedVolcano)
library(pheatmap)
library(RColorBrewer)

# helpers
to_named <- function(df, converter, rowcol) {
  res = converter(df |> select(!{{rowcol}}))
  rownames(res) = df |> pull({{rowcol}})
  return (res)
}


hmn_expr = read_csv(here('GSE167536', 'GSE167536_humanstress_rawcounts.txt')) |>
  to_named(as.matrix, 'gene')

hmn_mi = read_csv(here('GSE167536', 'BioSample', 'human_mi_status.csv')) |>
  dplyr::rename(SampleName = sample_name) |>
  select(SampleName, MI_status) |>
  mutate(MI_status = relevel(as.factor(MI_status), 'control'))

hmn_samples = read_csv(here('GSE167536', 'BioSample', 'human_sample_info.csv')) |>
  mutate(group = relevel(as.factor(str_replace(group, ' ', '_')), 'low_stress')) |>
  dplyr::rename(provider = `biomaterial provider`) |>
  left_join(hmn_mi, by = join_by(SampleName)) |>
  to_named(as.data.frame, 'SampleName') |>
  (function(df, o) df[o, ])(colnames(hmn_expr)) |>
  AnnotatedDataFrame()

hmn_eset = ExpressionSet(hmn_expr, hmn_samples)
hmn_se = makeSummarizedExperimentFromExpressionSet(hmn_eset)
# from https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html:
#   In order to benefit from the default settings of the package, 
#   you should put the variable of interest at the end of the formula 
#   and make sure the control level is the first level.
# this is different than what they had in the paper
# hmn_dss = DESeqDataSet(hmn_se, ~ group + MI_status)
hmn_dss = DESeqDataSet(hmn_se, ~ MI_status + group)

# do not filter the genes?

# proper reference levels were set when reading in the data
# and the design is set up to compare for the stress group
# But, maybe better if this was explicit in the call?
hmn_dss = DESeq(hmn_dss)

# they used pvalue < 0.05 for human DE genes
#hmn_contrast = c('group', 'low_stress', 'high_stress')
hmn_contrast = c('group', 'high_stress', 'low_stress')
hmn_res = results(hmn_dss, contrast = hmn_contrast)
hmn_res = hmn_res[complete.cases(hmn_res),]
EnhancedVolcano(hmn_res, 
                lab=rownames(hmn_res), 
                x='log2FoldChange', 
                y='pvalue',
                pCutoff = 0.05)

# shrink lfc for plotting
hmn.norm = lfcShrink(hmn_dss, contrast = hmn_contrast, type="normal")
hmn.norm = hmn.norm[complete.cases(hmn.norm),]
EnhancedVolcano(hmn.norm, 
                lab=rownames(hmn_res), 
                x='log2FoldChange', 
                y='pvalue',
                pCutoff = 0.05)


# results to match with the paper's Sup Fig 7
paper_res = hmn_res[hmn_res$pvalue < 0.05,]
nrow(paper_res) # should be 1959

# create pheatmap to match Figure 5B
# top 100 DE genes by significance on Z-scaled counts `rows`
top_100 = paper_res|> 
  as_tibble(rownames="gene") |>
  arrange(pvalue) |>
  head(100) |>
  pull(gene)

pheatmap(
  exprs(hmn_eset[top_100]), 
  color = colorRampPalette(rev(brewer.pal(n = 5, name ="RdBu")))(50),
  scale='row', 
  clustering_method = 'ward.D2', 
  annotation_col = pData(hmn_eset) |> select(group, MI_status),
  annotation_colors = list(
    group = c('low_stress' = 'blue', 'high_stress' = 'red'),
    MI_status = c('control' = 'limegreen', 'MI'='gold')
    ),
  show_rownames = FALSE,
  show_colnames = FALSE
  )

### This doesn't seem to match their paper's results...
# They did provide a list of their top human DE genes (along with pvalues)
# so, will that recreate their heatmap...? 

# save out all of the DESeq2 results for further analysis (pathway)
hmn_res |> 
  as_tibble(rownames="gene") |>
  arrange(gene) |>
  write_csv(here('Results','human-deseq2-results.csv'))

#### for Enrichr analysis
# which needs newline separated ENTREZ gene symbols
de_out = paper_res |> 
  as_tibble(rownames='gene')

de_out |>
  filter(log2FoldChange > 0) |>
  select(gene) |>
  write_delim(
    here('Results', 'human-DE-pathway-up-genes.txt'), 
    delim='\n', 
    col_names=FALSE
  )

de_out |>
  filter(log2FoldChange < 0) |>
  select(gene) |>
  write_delim(
    here('Results', 'human-DE-pathway-down-genes.txt'), 
    delim='\n', 
    col_names=FALSE
  )

