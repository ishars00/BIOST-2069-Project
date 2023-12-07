library(clusterProfiler)
library("org.Hs.eg.db", character.only = TRUE)
library(tidyverse)
library(here)

rrho_upreg = read_csv(
  here('Results', 'RRHO-lfc', 'RRHO_GO_MostUpregulatedHuman DE_VS_Mouse DE.csv'),
  col_names = c('gene'))

rrho_downreg = read_csv(
  here('Results', 'RRHO-lfc', 'RRHO_GO_MostDownregulatedHuman DE_VS_Mouse DE.csv'),
  col_names = c('gene'))

hmn_res = read_csv(here('Results', 'human-deseq2-results.csv'))

down_ego = enrichGO(
  gene = rrho_downreg |> pull(gene),
  universe = hmn_res |> pull(gene),
  OrgDb = org.Hs.eg.db,
  ont = "ALL",
  keyType = "SYMBOL",
  pvalueCutoff  = 0.05,
  pAdjustMethod = 'BH',
)

dotplot(down_ego)

down_ego@result |> 
  as_tibble() |>
  select(ID, qvalue) |>
  write_delim(here('Results', 'RRHO-lfc', 'downreg-go.tsv'), delim='\t')

up_ego = enrichGO(
  gene = rrho_upreg |> pull(gene),
  universe = hmn_res |> pull(gene),
  OrgDb = org.Hs.eg.db,
  ont = "ALL",
  keyType = "SYMBOL",
  pvalueCutoff  = 0.05,
  pAdjustMethod = 'BH',
)

dotplot(up_ego)

up_ego@result |> 
  as_tibble() |>
  select(ID, qvalue) |>
  write_delim(here('Results', 'RRHO-lfc', 'upreg-go.tsv'), delim='\t')
