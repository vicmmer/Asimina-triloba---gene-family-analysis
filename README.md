# Gene Family Evolution Pipeline (Pawpaw-Focused Comparative Genomics)

This repository contains scripts and workflows used to analyze **gene family evolution** across *Annonaceae* and related magnoliid taxa, with a **primary focus on pawpaw (*Asimina triloba*)–specific gene family changes**.

The workflow integrates **genome completeness assessment, orthogroup inference, protein functional annotation, gene family filtering, phylogenetic dating, and GO enrichment analyses**, enabling functional interpretation of pawpaw-specific evolutionary patterns.

---

## Repository Overview

This repository includes scripts that perform the following functions:
- Genome and proteome quality assessment (BUSCO)
- Orthogroup inference across Annonaceae and outgroups (OrthoFinder)
- Functional annotation of orthogroups (InterProScan)
- Filtering of transposable-element–associated annotations
- Identification of pawpaw-specific and shared gene families
- Gene family intersection analyses (UpSet plots)
- Phylogenetic dating for contextualizing gene family evolution (treePL)
- GO enrichment analysis of pawpaw-focused gene sets (topGO)

All scripts are intended to be run sequentially unless otherwise noted or can also be run with an automated script (stil working on it)

---

## Species and Data Sources

| Species | File Name | Source Link | Source Type | Year |
|------|----------|------------|------------|------|
| *Annona cherimola* | `Annona_cherimola.fa` | https://ihsmsubtropicals.uma.es/downloads/Annona%20cherimola | CNCB / UMA (annotated proteins) | 2023 |
| *Annona montana* | `Annona_montana.fa` | https://ngdc.cncb.ac.cn/gwh/Genome/154848/show | CNCB–GWH | 2023 |
| *Annona muricata* | `Annona_muricata_annotated.fasta` | DOI:10.1111/1755-0998.13353 | From authors (local copy) | 2021 |
| *Asimina triloba* | `Asimina_triloba_annotated.fasta` | Local assembly | Our lab’s assembly | 2023 |
| *Cinnamomum micranthum* | `Cinnamomum_micranthum.fa` | https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_003546025.1/ | NCBI Datasets | 2019 |
| *Lindera megaphylla* | `Lindera_megaphylla.fa` | https://ngdc.cncb.ac.cn/gwh/Assembly/26279/show | CNCB–GWH | 2023 |
| *Magnolia kwangsiensis* | `Magnolia_kwangsiensis.fa` | https://ngdc.cncb.ac.cn/gwh/Assembly/98328/show | CNCB–GWH | 2025 |
| *Persea americana* | `Persea_americana.fa` | https://genomevolution.org/coge/api/v1/genomes/29302/sequence | Science Data Bank | 2025 |

---

## Important Notes & Caveats

- **OrthoFinder output must be manually cleaned** prior to functional annotation.
- `interproscan.sh` **must be run inside the directory containing the cleaned OrthoFinder output**.
- After InterProScan completes, `filter_interpro_output.sh` **must be run** to:
  - Remove empty annotation files  
  - Remove annotations associated with transposable elements  

These filtering steps are critical for accurately identifying **pawpaw-specific gene family signals**.

---

##  Workflow Overview

Scripts are numbered in the recommended execution order.

- ### 0. Data preparation

0.download_accessions.sh

- ### 1. Genome Completeness
1.busco.sh
1b.busco_summaries.sh
  -Runs BUSCO an generates completeness summaries across species

- ### 2. Orthogroup inference
2.orthofinder.sh
Runs Orthofinder to identify orthougroups in taxa
2b.prepare_interproinput.sh
Prepares and puts Orthofinder output in ready-to-go form for functional annotation

- ### 3. Functional Annotation
3. interproscan.sh 

