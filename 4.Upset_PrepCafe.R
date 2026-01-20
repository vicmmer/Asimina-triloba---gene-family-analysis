#!/usr/bin/env Rscript

# Purpose:
#   (1) Make an UpSet plot for a subset of orthogroups
#       (based on include_orthogroups.txt)
#   (2) Create CAFE input table (gene family counts + GO desc)
#       and an ultrametric species tree for CAFE.
#
# Notes:
#   - Auto-detects newest orthofinder/Results_* folder
#   - Strips trailing "_annotated" from species names
# ------------------------------------------------------------

suppressPackageStartupMessages({
  library(UpSetR)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ape)
  library(phytools)
})

# ===== 0) Auto-detect newest OrthoFinder results =====
results_dir <- list.dirs("orthofinder", recursive = FALSE, full.names = TRUE)
results_dir <- results_dir[grepl("Results_", basename(results_dir))]
if (length(results_dir) == 0) stop("No orthofinder/Results_* directory found")

results_dir <- results_dir[order(file.info(results_dir)$mtime, decreasing = TRUE)]
OF_DIR <- results_dir[1]
cat("Using OrthoFinder results directory:", OF_DIR, "\n")

# ===== 1) Paths =====
orthogroups_file <- file.path(OF_DIR, "Orthogroups", "Orthogroups.GeneCount.tsv")
tree_file        <- file.path(OF_DIR, "Species_Tree", "SpeciesTree_rooted.txt")

include_file <- file.path("interproscan_input", "output", "include_orthogroups.txt")
go_file      <- file.path("interproscan_input", "output", "all_go_unique.tsv")

output_pdf   <- "upset_plot.pdf"
cafe_out     <- "cafe_gene_families.tsv"
tree_out     <- "SpeciesTree_ultrametric.txt"

# ===== 2) Read Orthogroups counts (once) =====
cat("Reading Orthogroups table...\n")
orthogroups_df <- read_tsv(orthogroups_file, show_col_types = FALSE)

# Clean species names (strip trailing "_annotated")
orthogroups_df <- orthogroups_df %>%
  rename_with(~ sub("_annotated$", "", .x))

# Identify species columns (everything except Orthogroup + Total)
species_cols <- setdiff(names(orthogroups_df), c("Orthogroup", "Total"))

# ===== 3) UpSet plot for included orthogroups =====
cat("Reading include list...\n")
include_ogs <- read_table(include_file, col_names = FALSE, show_col_types = FALSE) %>% pull(1)

cat("Filtering orthogroups...\n")
filtered_df <- orthogroups_df %>%
  filter(Orthogroup %in% include_ogs) %>%
  select(Orthogroup, all_of(species_cols))

# Convert counts -> presence/absence (only species numeric columns)
filtered_mat <- filtered_df
filtered_mat[species_cols] <- lapply(filtered_mat[species_cols], function(x) as.integer(x > 0))

cat("Generating UpSet plot...\n")
pdf(output_pdf, width = 10, height = 8)
upset(
  filtered_mat,
  nsets = length(species_cols),
  sets = rev(species_cols),
  keep.order = TRUE,
  order.by = "freq",
  number.angles = 30,
  empty.intersections = "on"
)
dev.off()
cat("Done! Output written to:", output_pdf, "\n\n")

# ===== 4) Prepare CAFE input =====
cat("Preparing CAFE input...\n")

counts <- orthogroups_df %>%
  rename(`Family ID` = Orthogroup) %>%
  select(`Family ID`, all_of(species_cols))  # drop Total automatically

cat("Reading GO table...\n")
go_df <- read_tsv(go_file, col_names = FALSE, show_col_types = FALSE)
colnames(go_df) <- c("FamilyID", "Annotation")

go_collapsed <- go_df %>%
  group_by(FamilyID) %>%
  summarise(Desc = paste(unique(Annotation), collapse = "|"), .groups = "drop")

cafe_input <- counts %>%
  mutate(FamilyID = `Family ID`) %>%
  select(-`Family ID`) %>%
  left_join(go_collapsed, by = "FamilyID") %>%
  mutate(Desc = replace_na(Desc, "(null)")) %>%
  select(Desc, FamilyID, everything()) %>%
  rename(`Family ID` = FamilyID)

write_tsv(cafe_input, cafe_out)
cat("Wrote CAFE table:", cafe_out, "\n\n")

# ===== 5) Make ultrametric tree for CAFE =====
cat("Reading species tree...\n")
tree <- read.tree(tree_file)

# Make tip labels match counts table (strip "_annotated")
tree$tip.label <- sub("_annotated$", "", tree$tip.label)

cat("Forcing ultrametric tree...\n")
ult.tree <- force.ultrametric(tree, method = "extend")
write.tree(ult.tree, file = tree_out)

cat("Wrote ultrametric tree:", tree_out, "\n")
cat("All done.\n")
