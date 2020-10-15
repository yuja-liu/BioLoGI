library(dplyr)
library(ggplot2)
library(magrittr)
library(ggplot2)
library(biomaRt)
library(tidyr)
library(BioCircos)
library(latex2exp)

# load logic gates
setwd("~/Documents/graduation_project/YABMSS")
detected_gates <- read.csv("data/GSE75748/detected_logi.tsv", sep = "\t")

# get gene symbol
detected_gates <- detected_gates %>%
  separate(col = name_inducer, into = c("id_inducer", "ver_inducer")) %>%
  separate(col = name_reporter, into = c("id_reporter", "ver_reporter")) %>%
  glimpse()
id_2_symbol <- read.csv("data/GSE75748/hg38_ensembl_id_2_symbol.tsv", sep = "\t") %>%
  rename(gene_id = Gene.stable.ID) %>%
  rename(gene_symbol = Gene.name) %>%
  glimpse()
detected_gates_w_symbol <- detected_gates %>%
  left_join(id_2_symbol, by = c(id_inducer = "gene_id")) %>%
  rename(symbol_inducer = gene_symbol) %>%
  left_join(id_2_symbol, by = c(id_reporter = "gene_id")) %>%
  rename(symbol_reporter = gene_symbol) %>%
  # fill NA in symbol by id
  mutate(symbol_inducer = coalesce(symbol_inducer, id_inducer),
         symbol_reporter = coalesce(symbol_reporter, id_reporter))

# filter Delta AIC
detected_gates_w_symbol <- detected_gates_w_symbol %>%
  filter((Delta_AIC < -20) & (AIC < -1000))

# summarize inducer genes
detected_gates_w_symbol <- detected_gates_w_symbol %>%
  group_by(id_inducer, symbol_inducer) %>%
  mutate(n_inducer = n()) %>%
  ungroup() %>%
  arrange(desc(n_inducer)) %>%
  glimpse()

# get gene coordinates
mart = useEnsembl(biomart="ensembl", 
                  dataset="hsapiens_gene_ensembl")
coordinates_inducer <- getBM(attributes = c("ensembl_gene_id", "chromosome_name", 
                                            "start_position", "end_position", "gene_biotype"),
                              filters = "ensembl_gene_id",
                              values = detected_gates %>% distinct(id_inducer) %>% pull(id_inducer),
                              mart = mart)
coordinates_reporter <- getBM(attributes = c("ensembl_gene_id", "chromosome_name", 
                                            "start_position", "end_position", "gene_biotype"),
                             filters = "ensembl_gene_id",
                             values = detected_gates %>% distinct(id_reporter) %>% pull(id_reporter),
                             mart = mart)
detected_gates_w_coordinates <- detected_gates_w_symbol %>%
  left_join(coordinates_inducer, by = c("id_inducer" = "ensembl_gene_id")) %>%
  rename('chr_inducer' = "chromosome_name",
         "start_inducer" = "start_position",
         "end_inducer" = "end_position",
         "type_inducer" = "gene_biotype") %>%
  left_join(coordinates_reporter, by = c("id_reporter" = "ensembl_gene_id")) %>%
  rename('chr_reporter' = "chromosome_name",
         "start_reporter" = "start_position",
         "end_reporter" = "end_position",
         "type_reporter" = "gene_biotype") %>%
  glimpse()

# Biotype conforming
biotype_conform <- function(biotype){
  if(is.na(biotype))
    return("others")
  if (biotype %in% c("lncRNA", "misc_RNA",
                     "snRNA", "snoRNA","scaRNA", "miRNA")) 
    return(biotype)
  if (grepl("pseudogene", biotype))
    return("pseudogene")
  if (biotype == "protein_coding" |
      grepl("TR|IG", biotype))
    return("protein_coding")
  return("others")
}
detected_gates_w_coordinates <- detected_gates_w_coordinates %>%
  mutate(type_inducer = sapply(type_inducer, biotype_conform)) %>%
  mutate(type_reporter = sapply(type_reporter, biotype_conform))
# save to file
detected_gates_w_coordinates %>% write.table(
  "data/GSE75748/detected_logi_w_coor.tsv",
  sep = '\t', row.names = FALSE)
# join to detected_gates_n
detected_gates_n <- detected_gates_w_coordinates %>%
  group_by(id_inducer, symbol_inducer, n_inducer, type_inducer) %>%
  distinct(id_inducer) %>%
  glimpse()

# plot
detected_gates_n %>%
  head(50) %>%
  ggplot(aes(x = reorder(id_inducer, desc(n_inducer)), y = n_inducer,
             fill = type_inducer)) +
  geom_bar(stat = "identity") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        axis.title.x = element_blank()) +
  scale_x_discrete(labels = detected_gates_n$symbol_inducer) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) +
  labs(y = "Number of ogic gates as an input", fill = "input RNA type")

# also violin plot
detected_gates_w_coordinates %>%
  filter(type_inducer != 'others') %>%
  ggplot(aes(x = type_inducer, y = abs(Delta_AIC), fill = type_inducer)) +
  geom_violin() +
  theme_classic() +
  labs(y = TeX("$\\Delta$AIC")) +
  theme(axis.title.x = element_blank(),
        legend.position = 'none')
  

# Circos Plot
selected_RNA_itrxn <- detected_gates_w_coordinates %>%
  dplyr::filter(id_inducer == "ENSG00000141431") %>%
  drop_na()
track_list = BioCircosLinkTrack("RNA-RNA interaction",
                                # add a zero-length link to label the start node
                                gene1Chromosomes = c(selected_RNA_itrxn$chr_inducer, 
                                                     selected_RNA_itrxn$chr_inducer[1]),
                                gene1Starts = c(selected_RNA_itrxn$start_inducer,
                                                selected_RNA_itrxn$start_inducer[1]),
                                gene1Ends = c(selected_RNA_itrxn$end_inducer,
                                              selected_RNA_itrxn$end_inducer[1]),
                                gene2Chromosomes = c(selected_RNA_itrxn$chr_reporter,
                                                     selected_RNA_itrxn$chr_inducer[1]),
                                gene2Starts = c(selected_RNA_itrxn$start_reporter,
                                                selected_RNA_itrxn$start_inducer[1]),
                                gene2Ends = c(selected_RNA_itrxn$end_reporter,
                                              selected_RNA_itrxn$end_inducer[1]),
                                maxRadius = 0.9,
                                labels = c(selected_RNA_itrxn$symbol_reporter,
                                           selected_RNA_itrxn$symbol_inducer[1]),
                                width = "0.05em",
                                labelSize = "0.5em",
                                labelPadding = 70)
track_list = track_list + BioCircosBarTrack("Delta AIC",
                                chromosomes = selected_RNA_itrxn$chr_reporter,
                                starts = selected_RNA_itrxn$start_reporter,
                                # using fake ends to achieve constant width
                                ends = selected_RNA_itrxn$start_reporter + 2000000,
                                values = abs(selected_RNA_itrxn$Delta_AIC),
                                maxRadius = 0.9,
                                minRadius = 1.0,
                                range = c(0, 200),
                                color = 'red')
BioCircos(track_list,
          genomeFillColor = "PuOr",
          chrPad = 0.02,
          displayGenomeBorder = FALSE,
          genomeTicksDisplay = FALSE,
          yChr = FALSE)
