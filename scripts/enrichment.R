library(dplyr)
library(ggplot2)
library(magrittr)
library(ggplot2)
library(biomaRt)
library(tidyr)

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

# summarize inducer genes
detected_gates_n <- detected_gates_w_symbol %>%
  group_by(id_inducer, symbol_inducer) %>%
  summarise(n_inducer = n()) %>%
  arrange(desc(n_inducer)) %>%
  glimpse()

# plot
detected_gates_n %>%
  head(20) %>%
  ggplot(aes(x = reorder(symbol_inducer, desc(n_inducer)), y = n_inducer)) +
  geom_bar(stat = "identity") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        axis.title.x = element_blank()) +
  labs(y = "Number of involved logic gates as an input")

# Also sort the original table to see related reporters
detected_gates_w_symbol <- detected_gates_w_symbol %>%
  left_join(detected_gates_n, by = "id_inducer") %>%
  arrange(desc(n_inducer))

# biomRt gene symbol conversion. Removed due to temporary shutdown of ensembl
#
# mart = useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl")
# gene_symbols_inducer <- getBM(attributes = c("hgnc_symbol", "ensembl_gene_id"),
#                       filters = "ensembl_gene_id",
#                       values = detected_gates$id_inducer %>% distinct(),
#                       mart = mart)
# gene_symbols_reporter <- getBM(attributes = c("hgnc_symbol", "ensembl_gene_id"),
#                               filters = "ensembl_gene_id",
#                               values = detected_gates$id_reporter %>% distinct(),
#                               mart = mart)
# detected_gates_w_symbol <- detected_gates %>%
#   left_join(gene_symbols_inducer, by = c("id_inducer" = "ensembl_gene_id")) %>%
#   rename(symbol_inducer = hgnc_symbol) %>%
#   mutate(symbol_inducer = replace(symbol_inducer, symbol_inducer == "", NA)) %>%
#   left_join(gene_symbols_reporter, by = c("id_reporter" = "ensembl_gene_id")) %>%
#   rename(symbol_reporter = hgnc_symbol) %>%
#   mutate(symbol_reporter = replace(symbol_reporter, symbol_reporter == "", NA)) %>%
#   glimpse()
