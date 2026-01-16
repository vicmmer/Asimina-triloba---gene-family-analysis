# Trying to figure out TOPGO analysis based on Suzy's code with some modifications
## tutorial she linked:
# http://avrilomics.blogspot.com/2015/07/using-topgo-to-test-for-go-term.html

# ---- Load packages ----
setwd("C:/Users/vicme/OneDrive/Desktop/Manuscript_Pawpaw") 
library(topGO)       # GO enrichment functions 
library(Rgraphviz)   # topGO needs this to draw the little GO network PDFs
library(readr)
library(dplyr)
library(tidyr)
library(stringr)


# Editing input to collapse entries 
all_tsv <- read.delim("all_go_unique_clean.tsv",
                      sep = "\t", 
                      header = FALSE)

colnames(all_tsv) <- c("Orthogroup", "GO")

#2) Collapse to one row per Orthogroup and remove duplicate GO terms
collapsed <- all_tsv %>%
  separate_rows(GO, sep = ";") %>%
  mutate(GO = str_trim(GO)) %>%
  filter(!is.na(GO), GO != "") %>%
  distinct(Orthogroup, GO) %>%               # drop duplicates within each OG
  arrange(Orthogroup, GO) %>%                # optional: sort terms
  group_by(Orthogroup) %>%
  summarise(GO = paste(GO, collapse = ";"), .groups = "drop")
write_tsv(collapsed, "all_go_collapsed.tsv") #THIS WILL BE INPUT FOR NEXT STEPS !

# ---- Define the gene universe (mapping from gene/OG -> GO terms) ----

geneID2GO <- readMappings(file = "all_go_collapsed.tsv")  # read the mapping (each line: gene <tab> GO1,GO2,...)
geneUniverse <- names(geneID2GO)  #all the orthogroup ids 

read_interest <- function(path) {
  x <- read.table(path, header = FALSE, stringsAsFactors = FALSE)  # no header, keep as characters
  as.character(x[[1]])                                             # return just the first column as a character vector
}

#### Import files from cafe #### 
chg <- read_tsv("Gamma_change.tab", show_col_types = FALSE)
chg <- chg %>%
  rename(Orthogroup = FamilyID) 
names(chg) <- names(chg) %>%
  str_replace("<[0-9]+>", "") 
chg <- chg[, !(is.na(names(chg)) | names(chg) == "")]

  
#Import significance file from cafe 
fam_sig <- read_tsv(
  "Gamma_family_results.txt",
  comment = "#",
  col_names = c("Orthogroup", "pvalue", "Significant")
)

chg2 <- chg %>%
  left_join(fam_sig, by = "Orthogroup")

table(chg2$Significant, useNA = "ifany")

pawpaw_any <- chg2 %>%
  filter(Asimina_triloba > 0) %>%
  pull(Orthogroup) %>%
  unique()
pawpaw_sig <- chg2 %>%
  filter(
    Asimina_triloba > 0,
    Significant == "y"
  ) %>%
  pull(Orthogroup) %>%
  unique()

annon_any <- chg2 %>%
  filter(
    Asimina_triloba > 0 |
      Annona_cherimola > 0 |
      Annona_montana > 0 |
      Annona_muricata > 0
  ) %>%
  pull(Orthogroup) %>%
  unique()

annon_sig <- chg2 %>%
  filter(
    Significant == "y",
    Asimina_triloba > 0 |
      Annona_cherimola > 0 |
      Annona_montana > 0 |
      Annona_muricata > 0
  ) %>%
  pull(Orthogroup) %>%
  unique()

# =========================
# ===== NEW CODE STARTS HERE
# =========================

# Pick which sets you want to run (swap _any vs _sig)
pawpaw_expansions <- pawpaw_sig   # or pawpaw_any
annon_ids         <- annon_sig    # or annon_any

# build the 0/1 selection vectors over the full universe
pawpaw_flag <- factor(as.integer(geneUniverse %in% pawpaw_expansions))
names(pawpaw_flag) <- geneUniverse

annon_flag <- factor(as.integer(geneUniverse %in% annon_ids))
names(annon_flag) <- geneUniverse

# ---- 1. Pawpaw figure (BP and MF) ----

# Pawpaw  - BP (Biological Processes)
myGOdataBP <- new("topGOdata",
                  description = "IC unique OGs - Pawpaw BP",
                  ontology    = "BP",
                  allGenes    = pawpaw_flag,
                  annot       = annFUN.gene2GO,
                  gene2GO     = geneID2GO)

resultFisher <- runTest(myGOdataBP, algorithm = "classic", statistic = "fisher")

bp_scores <- score(resultFisher)
bp_n_all  <- sum(!is.na(bp_scores))

if (bp_n_all > 0) {
  pawpaw_BP <- GenTable(myGOdataBP,
                        classicFisher = resultFisher,
                        orderBy       = "classicFisher",
                        ranksOf       = "classicFisher",
                        topNodes      = bp_n_all)
  
  pawpaw_BP$P_Value <- bp_scores[pawpaw_BP$GO.ID]
  pawpaw_BP$FDR_BH  <- p.adjust(bp_scores, method = "BH")[pawpaw_BP$GO.ID]
  
  names(pawpaw_BP)[names(pawpaw_BP) == "GO.ID"]       <- "GO_ID"
  names(pawpaw_BP)[names(pawpaw_BP) == "Term"]        <- "GO_Term"
  names(pawpaw_BP)[names(pawpaw_BP) == "Annotated"]   <- "N_Annotated_Universe"
  names(pawpaw_BP)[names(pawpaw_BP) == "Significant"] <- "N_Significant_Selected"
  names(pawpaw_BP)[names(pawpaw_BP) == "Expected"]    <- "Expected_Significant"
  
  pawpaw_BP$classicFisher <- NULL
  pawpaw_BP$Label    <- "pawpaw"
  pawpaw_BP$Category <- "BP"
  
  pawpaw_BP <- pawpaw_BP[, c("Label","Category","GO_ID","GO_Term",
                             "N_Annotated_Universe","N_Significant_Selected","Expected_Significant",
                             "P_Value","FDR_BH"), drop = FALSE]
  
  write.csv(pawpaw_BP, "pawpaw_BP_topgo.csv", row.names = FALSE)
} else {
  write.csv(data.frame(Label=character(),Category=character(),GO_ID=character(),GO_Term=character(),
                       N_Annotated_Universe=integer(),N_Significant_Selected=integer(),Expected_Significant=numeric(),
                       P_Value=numeric(),FDR_BH=numeric()),
            "pawpaw_BP_topgo.csv", row.names = FALSE)
}

showSigOfNodes(myGOdataBP, bp_scores, firstSigNodes = min(10, max(1, bp_n_all)), useInfo = "all")
printGraph(myGOdataBP, resultFisher, firstSigNodes = min(10, max(1, bp_n_all)),
           fn.prefix = "tGO_pawpaw_BP", useInfo = "all", pdfSW = TRUE)

# Pawpaw  - MF (molecular function)
myGOdataMF <- new("topGOdata",
                  description = "IC unique OGs - Pawpaw MF",
                  ontology    = "MF",
                  allGenes    = pawpaw_flag,
                  annot       = annFUN.gene2GO,
                  gene2GO     = geneID2GO)

resultFisher <- runTest(myGOdataMF, algorithm = "classic", statistic = "fisher")

mf_scores <- score(resultFisher)
mf_n_all  <- sum(!is.na(mf_scores))

if (mf_n_all > 0) {
  pawpaw_MF <- GenTable(myGOdataMF,
                        classicFisher = resultFisher,
                        orderBy       = "classicFisher",
                        ranksOf       = "classicFisher",
                        topNodes      = mf_n_all)
  
  pawpaw_MF$P_Value <- mf_scores[pawpaw_MF$GO.ID]
  pawpaw_MF$FDR_BH  <- p.adjust(mf_scores, method = "BH")[pawpaw_MF$GO.ID]
  
  names(pawpaw_MF)[names(pawpaw_MF) == "GO.ID"]       <- "GO_ID"
  names(pawpaw_MF)[names(pawpaw_MF) == "Term"]        <- "GO_Term"
  names(pawpaw_MF)[names(pawpaw_MF) == "Annotated"]   <- "N_Annotated_Universe"
  names(pawpaw_MF)[names(pawpaw_MF) == "Significant"] <- "N_Significant_Selected"
  names(pawpaw_MF)[names(pawpaw_MF) == "Expected"]    <- "Expected_Significant"
  
  pawpaw_MF$classicFisher <- NULL
  pawpaw_MF$Label    <- "pawpaw"
  pawpaw_MF$Category <- "MF"
  
  pawpaw_MF <- pawpaw_MF[, c("Label","Category","GO_ID","GO_Term",
                             "N_Annotated_Universe","N_Significant_Selected","Expected_Significant",
                             "P_Value","FDR_BH"), drop = FALSE]
  
  write.csv(pawpaw_MF, "pawpaw_MF_topgo.csv", row.names = FALSE)
} else {
  write.csv(data.frame(Label=character(),Category=character(),GO_ID=character(),GO_Term=character(),
                       N_Annotated_Universe=integer(),N_Significant_Selected=integer(),Expected_Significant=numeric(),
                       P_Value=numeric(),FDR_BH=numeric()),
            "pawpaw_MF_topgo.csv", row.names = FALSE)
}

showSigOfNodes(myGOdataMF, mf_scores, firstSigNodes = min(10, max(1, mf_n_all)), useInfo = "all")
printGraph(myGOdataMF, resultFisher, firstSigNodes = min(10, max(1, mf_n_all)),
           fn.prefix = "tGO_pawpaw_MF", useInfo = "all", pdfSW = TRUE)

#---- 2. Annonaceae figure (BP, CC and MF) ----

# Annonaceae - BP
myGOdataBP <- new("topGOdata",
                  description = "IC unique OGs - Annonaceae BP",
                  ontology    = "BP",
                  allGenes    = annon_flag,
                  annot       = annFUN.gene2GO,
                  gene2GO     = geneID2GO)

resultFisher <- runTest(myGOdataBP, algorithm = "classic", statistic = "fisher")

bp_scores <- score(resultFisher)
bp_n_all  <- sum(!is.na(bp_scores))

if (bp_n_all > 0) {
  annon_BP <- GenTable(myGOdataBP,
                       classicFisher = resultFisher,
                       orderBy       = "classicFisher",
                       ranksOf       = "classicFisher",
                       topNodes      = bp_n_all)
  
  annon_BP$P_Value <- bp_scores[annon_BP$GO.ID]
  annon_BP$FDR_BH  <- p.adjust(bp_scores, method = "BH")[annon_BP$GO.ID]
  
  names(annon_BP)[names(annon_BP) == "GO.ID"]       <- "GO_ID"
  names(annon_BP)[names(annon_BP) == "Term"]        <- "GO_Term"
  names(annon_BP)[names(annon_BP) == "Annotated"]   <- "N_Annotated_Universe"
  names(annon_BP)[names(annon_BP) == "Significant"] <- "N_Significant_Selected"
  names(annon_BP)[names(annon_BP) == "Expected"]    <- "Expected_Significant"
  
  annon_BP$classicFisher <- NULL
  annon_BP$Label    <- "annonaceae"
  annon_BP$Category <- "BP"
  
  annon_BP <- annon_BP[, c("Label","Category","GO_ID","GO_Term",
                           "N_Annotated_Universe","N_Significant_Selected","Expected_Significant",
                           "P_Value","FDR_BH"), drop = FALSE]
  
  write.csv(annon_BP, "annonaceae_BP_topgo.csv", row.names = FALSE)
} else {
  write.csv(data.frame(Label=character(),Category=character(),GO_ID=character(),GO_Term=character(),
                       N_Annotated_Universe=integer(),N_Significant_Selected=integer(),Expected_Significant=numeric(),
                       P_Value=numeric(),FDR_BH=numeric()),
            "annonaceae_BP_topgo.csv", row.names = FALSE)
}

showSigOfNodes(myGOdataBP, bp_scores, firstSigNodes = min(10, max(1, bp_n_all)), useInfo = "all")
printGraph(myGOdataBP, resultFisher, firstSigNodes = min(10, max(1, bp_n_all)),
           fn.prefix = "tGO_annonaceae_BP", useInfo = "all", pdfSW = TRUE)

# Annonaceae - MF
myGOdataMF <- new("topGOdata",
                  description = "IC unique OGs - Annonaceae MF",
                  ontology    = "MF",
                  allGenes    = annon_flag,
                  annot       = annFUN.gene2GO,
                  gene2GO     = geneID2GO)

resultFisher <- runTest(myGOdataMF, algorithm = "classic", statistic = "fisher")

mf_scores <- score(resultFisher)
mf_n_all  <- sum(!is.na(mf_scores))

if (mf_n_all > 0) {
  annon_MF <- GenTable(myGOdataMF,
                       classicFisher = resultFisher,
                       orderBy       = "classicFisher",
                       ranksOf       = "classicFisher",
                       topNodes      = mf_n_all)
  
  annon_MF$P_Value <- mf_scores[annon_MF$GO.ID]
  annon_MF$FDR_BH  <- p.adjust(mf_scores, method = "BH")[annon_MF$GO.ID]
  
  names(annon_MF)[names(annon_MF) == "GO.ID"]       <- "GO_ID"
  names(annon_MF)[names(annon_MF) == "Term"]        <- "GO_Term"
  names(annon_MF)[names(annon_MF) == "Annotated"]   <- "N_Annotated_Universe"
  names(annon_MF)[names(annon_MF) == "Significant"] <- "N_Significant_Selected"
  names(annon_MF)[names(annon_MF) == "Expected"]    <- "Expected_Significant"
  
  annon_MF$classicFisher <- NULL
  annon_MF$Label    <- "annonaceae"
  annon_MF$Category <- "MF"
  
  annon_MF <- annon_MF[, c("Label","Category","GO_ID","GO_Term",
                           "N_Annotated_Universe","N_Significant_Selected","Expected_Significant",
                           "P_Value","FDR_BH"), drop = FALSE]
  
  write.csv(annon_MF, "annonaceae_MF_topgo.csv", row.names = FALSE)
} else {
  write.csv(data.frame(Label=character(),Category=character(),GO_ID=character(),GO_Term=character(),
                       N_Annotated_Universe=integer(),N_Significant_Selected=integer(),Expected_Significant=numeric(),
                       P_Value=numeric(),FDR_BH=numeric()),
            "annonaceae_MF_topgo.csv", row.names = FALSE)
}

showSigOfNodes(myGOdataMF, mf_scores, firstSigNodes = min(10, max(1, mf_n_all)), useInfo = "all")
printGraph(myGOdataMF, resultFisher, firstSigNodes = min(10, max(1, mf_n_all)),
           fn.prefix = "tGO_annonaceae_MF", useInfo = "all", pdfSW = TRUE)

# Annonaceae - CC
myGOdataCC <- new("topGOdata",
                  description = "IC unique OGs - Annonaceae CC",
                  ontology    = "CC",
                  allGenes    = annon_flag,
                  annot       = annFUN.gene2GO,
                  gene2GO     = geneID2GO)

resultFisher <- runTest(myGOdataCC, algorithm = "classic", statistic = "fisher")

cc_scores <- score(resultFisher)
cc_n_all  <- sum(!is.na(cc_scores))

if (cc_n_all > 0) {
  annon_CC <- GenTable(myGOdataCC,
                       classicFisher = resultFisher,
                       orderBy       = "classicFisher",
                       ranksOf       = "classicFisher",
                       topNodes      = cc_n_all)
  
  annon_CC$P_Value <- cc_scores[annon_CC$GO.ID]
  annon_CC$FDR_BH  <- p.adjust(cc_scores, method = "BH")[annon_CC$GO.ID]
  
  names(annon_CC)[names(annon_CC) == "GO.ID"]       <- "GO_ID"
  names(annon_CC)[names(annon_CC) == "Term"]        <- "GO_Term"
  names(annon_CC)[names(annon_CC) == "Annotated"]   <- "N_Annotated_Universe"
  names(annon_CC)[names(annon_CC) == "Significant"] <- "N_Significant_Selected"
  names(annon_CC)[names(annon_CC) == "Expected"]    <- "Expected_Significant"
  
  annon_CC$classicFisher <- NULL
  annon_CC$Label    <- "annonaceae"
  annon_CC$Category <- "CC"
  
  annon_CC <- annon_CC[, c("Label","Category","GO_ID","GO_Term",
                           "N_Annotated_Universe","N_Significant_Selected","Expected_Significant",
                           "P_Value","FDR_BH"), drop = FALSE]
  
  write.csv(annon_CC, "annonaceae_CC_topgo.csv", row.names = FALSE)
} else {
  write.csv(data.frame(Label=character(),Category=character(),GO_ID=character(),GO_Term=character(),
                       N_Annotated_Universe=integer(),N_Significant_Selected=integer(),Expected_Significant=numeric(),
                       P_Value=numeric(),FDR_BH=numeric()),
            "annonaceae_CC_topgo.csv", row.names = FALSE)
}

showSigOfNodes(myGOdataCC, cc_scores, firstSigNodes = min(10, max(1, cc_n_all)), useInfo = "all")
printGraph(myGOdataCC, resultFisher, firstSigNodes = min(10, max(1, cc_n_all)),
           fn.prefix = "tGO_annonaceae_CC", useInfo = "all", pdfSW = TRUE)

# ---- 3. Combine per-ontology CSVS into one tsv per label  ----
combine_topgo <- function(label) {
  files <- list.files(pattern = paste0("^", label, "_(BP|MF|CC)_topgo\\.csv$"),
                      ignore.case = TRUE)
  if (!length(files)) { message("No files for ", label); return(invisible(NULL)) }
  
  dfs <- lapply(files, function(f) {
    d <- read.csv(f, check.names = FALSE)
    if (nrow(d) == 0) return(NULL)
    d
  })
  keep <- !vapply(dfs, is.null, logical(1))
  if (!any(keep)) { message("All files empty for ", label); return(invisible(NULL)) }
  dfs <- dfs[keep]
  
  all_cols <- Reduce(union, lapply(dfs, names))
  dfs <- lapply(dfs, function(d) {
    miss <- setdiff(all_cols, names(d))
    if (length(miss)) for (m in miss) d[[m]] <- NA
    d[all_cols]
  })
  out <- do.call(rbind, dfs)
  
  front <- c("Label","Category","GO_ID","GO_Term",
             "N_Annotated_Universe","N_Significant_Selected","Expected_Significant",
             "P_Value","FDR_BH")
  front <- front[front %in% names(out)]
  out <- out[, c(front, setdiff(names(out), front)), drop = FALSE]
  
  write.table(out, paste0(label, "_ALL_topgo.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  message("Wrote ", label, "_ALL_topgo.tsv with ", nrow(out), " rows (", length(dfs), " files)")
}

combine_topgo("pawpaw")
combine_topgo("annonaceae")

#  ---- 4. Load the combined TSVs back in and look at them  ----
annonaceae <- read.delim("annonaceae_ALL_topgo.tsv", sep = "\t",
                         stringsAsFactors = FALSE, check.names = FALSE)
View(annonaceae)

pawpaw <- read.delim("pawpaw_ALL_topgo.tsv", sep = "\t",
                     stringsAsFactors = FALSE, check.names = FALSE)
View(pawpaw)

###acetogenin: 
acetogenin_mf_keywords <- c("oxidoreductase","dehydrogenase","reductase","monooxygenase","NAD","FAD")

pawpaw_aceto_MF_anyP <- pawpaw %>%
  filter(Category == "MF") %>%
  filter(grepl(paste(acetogenin_mf_keywords, collapse="|"), GO_Term, ignore.case=TRUE))

nrow(pawpaw_aceto_MF_anyP)
head(pawpaw_aceto_MF_anyP, 20)

View(pawpaw_aceto_MF)



