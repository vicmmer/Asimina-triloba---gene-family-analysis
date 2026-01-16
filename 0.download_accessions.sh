#!/bin/bash
set -e

# ===========================
#   SETUP
# ===========================
mkdir -p downloads
FINAL_DIR="."
echo "=== Downloading genomes and proteins ==="

# ===========================
#   Annona cherimola (GFF only)
# ===========================
echo "=== Annona cherimola ==="
wget -O downloads/Annona_cherimola.faa.gz \
  "https://ihsmsubtropicals.uma.es/downloads/Annona%20cherimola/Sequences/anche102_proteins_annot.fasta.gz"

gunzip -c downloads/Annona_cherimola.faa.gz > Annona_cherimola.fa
# ===========================
#   Annona montana – from CNCB (Protein)
# ===========================
echo "=== Annona montana ==="
wget -O downloads/Annona_montana.faa.gz \
  "https://download.cncb.ac.cn/gwh/Plants/Annona_montana_Am_v1.0_GWHDQZG00000000/GWHDQZG00000000.Protein.faa.gz"

gunzip -c downloads/Annona_montana.faa.gz > Annona_montana.fa

# ===========================
#   Cinnamomum micranthum – NCBI
# ===========================
echo "=== Cinnamomum micranthum ==="
datasets download genome accession GCA_003546025.1 \
  --include genome,protein,gff3 \
  --filename downloads/Cinnamomum_micranthum.zip

unzip downloads/Cinnamomum_micranthum.zip -d downloads/Cinnamomum_micranthum_unzip

cp downloads/Cinnamomum_micranthum_unzip/ncbi_dataset/data/*/protein.faa \
  Cinnamomum_micranthum.fa

# ===========================
#   Lindera megaphylla – CNCB
# ===========================
echo "=== Lindera megaphylla ==="
wget -O downloads/Lindera_megaphylla.faa.gz \
  "https://download.cncb.ac.cn/gwh/Plants/Lindera_megaphylla_LMv1_GWHBKHA00000000/GWHBKHA00000000.Protein.faa.gz"

gunzip -c downloads/Lindera_megaphylla.faa.gz > Lindera_megaphylla.fa

# ===========================
#   Magnolia kwangsiensis – CNCB
# ===========================
echo "=== Magnolia kwangsiensis ==="
wget -O downloads/Magnolia_kwangsiensis.faa.gz \
  "https://download.cncb.ac.cn/gwh/Plants/Magnolia_kwangsiensis_Mkmh2.hap2_GWHGEUJ00000000.1/GWHGEUJ00000000.1.Protein.faa.gz"

gunzip -c downloads/Magnolia_kwangsiensis.faa.gz > Magnolia_kwangsiensis.fa

# ===========================
#   Persea americana
# ===========================
echo "=== Downloaded directly from link provided in publication onto computer>


# ===========================
#   Annona muricata & Asimina triloba (local)
# ===========================
echo "=== Already had Annona muricata & Asimina triloba; Asimina triloba is our own assembly and Annona muricata was obtained by contacting authors of  doi: 10.1111/1755-0998.13353 ==="


# ===========================
#   Done!
# ===========================
echo "=== All final protein files ready in: $FINAL_DIR ==="
ls -lh *.fa
