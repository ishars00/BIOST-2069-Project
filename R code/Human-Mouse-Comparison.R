if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# BiocManager::install("RRHO")

library("AnnotationDbi")
library("org.Mm.eg.db")
library("org.Hs.eg.db", character.only = TRUE)
library(clusterProfiler)

library(RRHO)
library(here)
library(tidyverse)
library(ggplot2)


# get the proper gene sets to convert mouse genes to human orthologs
# data retrieved directly from BioMart
# https://www.ensembl.info/2009/01/21/how-to-get-all-the-orthologous-genes-between-two-species/
# Also, have to map between the mouse RefSeq gene symbols and the Ensembl gene IDs 
# from BioMart
species_map = read_csv(here('Results', 'mus-to-sapiens-genemap.csv')) |>
  filter(!is.na(`Gene name`)) |>
  dplyr::rename(hum_gene = `Human gene name`, mus_gene = `Gene name`) |>
  mutate(mus_refseq = mapIds(
    org.Mm.eg.db,
    keys = `Gene stable ID`,
    keytype = "ENSEMBL",
    column = "REFSEQ",
    multiVals = first
  )) |>
  filter(!is.na(mus_refseq))

# load in human DE results and filter based on the paper's criteria
### p-value < 0.05, baseMeans ∼4 or larger ==> 1959 genes
human_de = read_csv(here('Results', 'human-deseq2-results.csv')) |>
  arrange(pvalue)

human_filt = human_de |>
  filter(pvalue < 0.05)

# load in mouse DE individual RNA-Seq results
### adjusted p-value < 0.05 ==> 5926 genes
mouse_de = read_csv(here('Results', 'mus-individual-deseq2-results.csv')) |>
  arrange(padj)

mouse_de |>
  filter(padj < 0.05) |>
  nrow()

# map mouse genename to human orthology genename
### 5926 mus genes ==> 3676 human ortholog genes
# we got 6990 mouse DE ==> 3318 homologues
mouse_mapped = mouse_de |> 
  inner_join(species_map, by = join_by(gene == mus_refseq))

mouse_mapped_filt = mouse_mapped |>
  filter(padj < 0.05)
  
# this filters down to 3318 mouse to human mapped genes

# create tibble with mouse and human => log2FoldChange, pvalue, padj
# this is pre-filtered, so re-apply the filters to get the filtered dataset
# they used in the paper

comparison = human_de |>
  dplyr::select(gene, log2FoldChange, pvalue, padj) |>
  rename_with(\(x) paste0("hum_",x)) |>
  inner_join(
    mouse_mapped |>
      dplyr::select(hum_gene, mus_gene, log2FoldChange, pvalue, padj) |>
      rename_with(\(x) paste0('mus_',x), .cols=c(log2FoldChange, pvalue, padj)),
    by = join_by(hum_gene)) |>
  dplyr::rename(hum_l2fc = hum_log2FoldChange, mus_l2fc = mus_log2FoldChange) |>
  mutate(category = case_when(
    hum_l2fc > 0 & mus_l2fc > 0 ~ "Up in Both",
    hum_l2fc < 0 & mus_l2fc < 0 ~ "Down in Both",
    hum_l2fc > 0 & mus_l2fc < 0 ~ "Up Human, Down Mouse",
    hum_l2fc < 0 & mus_l2fc > 0 ~ "Down Human, Up Mouse",
  )) |>
  mutate(category = as.factor(category))

# This has 339 observations, as compared to 421 in the paper
comparison_filt = comparison |>
  filter(hum_pvalue < 0.05 & mus_padj < 0.05)

# make double volcano plot
comparison_filt |>
  ggplot() +
  geom_point(aes(mus_l2fc, hum_l2fc, color=category), alpha=0.75) +
  scale_x_continuous('Mouse Log2 Fold Change', limits=c(-2,2)) +
  scale_y_continuous('Human Log2 Fold Change', limits=c(-2,2)) +
  scale_color_manual(values = c(
    "Up in Both" = "red",
    "Down in Both" = "dodgerblue",
    "Up Human, Down Mouse" = "limegreen",
    "Down Human, Up Mouse" = "gold"
  ))

comparison_filt |>
  group_by(category) |>
  summarise(n = n()) |>
  arrange(category)

# make RRHO plot
# this uses the full DE results from the mouse and human datasets
df_rrho = comparison |>
  dplyr::select(hum_gene, hum_pvalue, hum_l2fc, mus_pvalue, mus_l2fc) |>
  mutate(
    hum_rrho = -log10(hum_pvalue) * sign(hum_l2fc),
    mus_rrho = -log10(mus_pvalue) * sign(mus_l2fc),
  ) |>
  group_by(hum_gene) |>
  summarise(across(everything(), first))

res_rrho = RRHO(
  df_rrho |> select(hum_gene, hum_rrho) |> as.data.frame(),
  df_rrho |> select(hum_gene, mus_rrho) |> as.data.frame(),
  labels = c('Human DE', 'Mouse DE'),
  alternative = 'enrichment',
  plots = TRUE,
  outputdir = here('Results', 'RRHO-pvalue')
)

lattice::levelplot(res_rrho$hypermat)

# try running with the log fold change only
res_lfc = RRHO(
  df_rrho |> select(hum_gene, hum_l2fc) |> as.data.frame(),
  df_rrho |> select(hum_gene, mus_l2fc) |> as.data.frame(),
  alternative = 'enrichment',
  labels = c('Human DE', 'Mouse DE'),
  plots = TRUE,
  outputdir = here('Results', 'RRHO-lfc')
)

lattice::levelplot(res_lfc$hypermat)

# run RRHO on just the DE genes from human and mouse

df_filt_rrho = comparison_filt |>
  dplyr::select(hum_gene, hum_pvalue, hum_l2fc, mus_pvalue, mus_l2fc) |>
  mutate(
    hum_rrho = -log10(hum_pvalue) * sign(hum_l2fc),
    mus_rrho = -log10(mus_pvalue) * sign(mus_l2fc),
  ) |>
  group_by(hum_gene) |>
  summarise(across(everything(), first))

res_filtered_rrho = RRHO(
  df_filt_rrho |> select(hum_gene, hum_rrho) |> as.data.frame(),
  df_filt_rrho |> select(hum_gene, mus_rrho) |> as.data.frame(),
  labels = c('Human DE', 'Mouse DE'),
  alternative = 'enrichment',
)

lattice::levelplot(res_filtered_rrho$hypermat)

# Try to make a rank list for the RRHO website
# https://systems.crump.ucla.edu/rankrank/rankranksimple.php#example


# pathway analysis?
# get the genes that are concurrently expressed in both human and mouse
gene_up_hum_mus = comparison_filt |>
  filter(category == "Up in Both" | category == "Down in Both") |>
  pull(hum_gene)

ego = enrichGO(
  gene = gene_up_hum_mus,
  universe = comparison |> pull(hum_gene),
  OrgDb = org.Hs.eg.db,
  ont = "ALL",
  keyType = "SYMBOL",
  pvalueCutoff  = 0.05,
  pAdjustMethod = 'BH',
  readable = TRUE
)

dotplot(ego)
# save out enriched pathways for REVIGO analysis
ego@result |>
  as_tibble() |>
  select(ID, qvalue) |>
  write_delim(here('Results', 'Human-Mouse-Enriched-Pathways.txt'), delim='\t')
