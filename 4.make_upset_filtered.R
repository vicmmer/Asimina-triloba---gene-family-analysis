#!/usr/bin/env Rscript

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
