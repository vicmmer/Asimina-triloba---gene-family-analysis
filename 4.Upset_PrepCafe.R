#!/usr/bin/env Rscripte


library(UpSetR)

# ==== 1. Paths ====


orthogroups_file <- file.path("orthofinder/Results_Dec03/Orthogroups/Orthogroups.GeneCount.tsv")
include_file     <- file.path("interproscan_input/output/include_orthogroups.txt")
output_pdf       <- file.path("upset_plot.pdf")

cat("Reading Orthogroups table...\n")
orthogroups_df <- read.table(orthogroups_file, header = TRUE, stringsAsFactors = FALSE, sep="\t")

cat("Reading include list...\n")
include_ogs <- read.table(include_file, header = FALSE, stringsAsFactors = FALSE)[,1]

# ==== 2. Filter ====
cat("Filtering orthogroups...\n")
filtered_df <- orthogroups_df[orthogroups_df$Orthogroup %in% include_ogs, ]

# Make binary presence/absence matrix (convert >0 → 1)
filtered_df[ filtered_df > 0 ] <- 1

# Select all species columns (excluding Orthogroup + Total column)
selected_species <- colnames(orthogroups_df)[2:(ncol(orthogroups_df)-1)]

# ==== 3. Plot ====
cat("Generating UpSet plot...\n")
pdf(output_pdf, width = 10, height = 8)

upset(filtered_df,
      nsets           = length(selected_species),
      sets            = rev(selected_species),
      keep.order      = TRUE,
      order.by        = "freq",
      number.angles   = 30,
      empty.intersections = "on"
)

dev.off()

cat("Done! Output written to: ", output_pdf, "\n")

################## Now prepare Input for CAFE ########################
library(dplyr)
library(tidyr)
library(readr)
library(ape)
library(phytools)

## ---- Paths (relative to ~/genefamily) ----
orthogroups_file <- "orthofinder/Results_Dec03/Orthogroups/Orthogroups.GeneCount.tsv"
go_file          <- "interproscan_input/output/all_go_unique.tsv"
tree_file        <- "orthofinder/Results_Dec03/Species_Tree/SpeciesTree_rooted.txt"

## ---- 1) Read Orthogroups.GeneCount.tsv ----
orthogroups_df <- read_tsv(orthogroups_file, show_col_types = FALSE)

counts <- orthogroups_df %>%
  rename(FamilyID = Orthogroup) %>%
  select(-Total) %>%
  rename_with(~ sub("_annotated$", "", .x), -FamilyID)

## ---- 2) Read GO table (OG \t GO stuff) ----
go_df <- read.delim(go_file, header = FALSE, stringsAsFactors = FALSE)
colnames(go_df) <- c("FamilyID", "Annotation")

go_collapsed <- go_df %>%
  group_by(FamilyID) %>%
  summarise(Desc = paste(unique(Annotation), collapse = "|"), .groups = "drop")

## ---- 3) Join + format for CAFE ----
cafe_input <- counts %>%
  left_join(go_collapsed, by = "FamilyID") %>%
  mutate(Desc = tidyr::replace_na(Desc, "(null)")) %>%
  relocate(Desc, .before = FamilyID) %>%
  rename(`Family ID` = FamilyID)

write_tsv(cafe_input, "cafe_gene_families.tsv")

## ---- 4) Make ultrametric tree for CAFE ----
tree <- read.tree(tree_file)

# Make tip labels match the species names in cafe_gene_families.tsv
# (strip any trailing "_annotated")
tree$tip.label <- sub("_annotated$", "", tree$tip.label)

ult.tree <- force.ultrametric(tree, method = "extend")
write.tree(ult.tree, file = "SpeciesTree_ultrametric.txt")
