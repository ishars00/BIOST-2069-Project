if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# BiocManager::install("RRHO")

library(RRHO)
library(here)
library(tidyverse)


# get the proper gene sets to convert mouse genes to human orthologs
# data retrieved directly from BioMart
# https://www.ensembl.info/2009/01/21/how-to-get-all-the-orthologous-genes-between-two-species/
species_map = read_csv(here('Results', 'mus-to-sapiens-genemap.csv')) |>
  filter(!is.na(`Gene name`))

# load in human DE results and filter based on the paper's criteria
#   p-value < 0.05, baseMeans ∼4 or larger ==> 1959 genes
human_de = read_csv(here('Results', 'human-deseq2-results.csv')) |>
  filter(pvalue < 0.05 & baseMean >= 4) |>
  arrange(pvalue)
