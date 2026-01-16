# paw-paw-visualization
| Species                        | File Name                         | Source Link                                                                                                                    | Source Type                   | Year |
| ------------------------------ | --------------------------------- | ------------------------------------------------------------------------------------------------------------------------------ | ----------------------------- | ---- |
| *Annona cherimola*             | `Annona_cherimola.fa`             | [https://ihsmsubtropicals.uma.es/downloads/Annona%20cherimola](https://ihsmsubtropicals.uma.es/downloads/Annona%20cherimola)   | CNCB/UMA (annotated proteins) | 2023 |
| *Annona montana*               | `Annona_montana.fa`               | [https://ngdc.cncb.ac.cn/gwh/Genome/154848/show](https://ngdc.cncb.ac.cn/gwh/Genome/154848/show)                               | CNCB–GWH                      | 2023 |
| *Annona muricata*              | `Annona_muricata_annotated.fasta` | Provided by authors of DOI:10.1111/1755-0998.13353                                                                             | From authors (local)          | 2021 |
| *Asimina triloba*              | `Asimina_triloba_annotated.fasta` | Local assembly                                                                                                                 | Our lab’s assembly            | 2023 |
| *Cinnamomum micranthum*        | `Cinnamomum_micranthum.fa`        | [https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_003546025.1/](https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_003546025.1/) | NCBI Datasets                 | 2019 |
| *Lindera megaphylla*           | `Lindera_megaphylla.fa`           | [https://ngdc.cncb.ac.cn/gwh/Assembly/26279/show](https://ngdc.cncb.ac.cn/gwh/Assembly/26279/show)                             | CNCB–GWH                      | 2023 |
| *Magnolia kwangsiensis*        | `Magnolia_kwangsiensis.fa`        | [https://ngdc.cncb.ac.cn/gwh/Assembly/98328/show](https://ngdc.cncb.ac.cn/gwh/Assembly/98328/show)                             | CNCB–GWH                      | 2025 |
| *Persea americana*  | `Persea_americana.fa`         | [[https://genomevolution.org/coge/api/v1/genomes/29302/sequence](https://genomevolution.org/coge/api/v1/genomes/29302/sequence](https://academic.oup.com/hr/article/11/7/uhae119/7656377#471718547)) | Science Data Bank                          | 2025 |

Note: Orthofinder output had to be cleaned, and the interpro script had to be run within a folder with the cleaned orthofinder output, only way it worked 
After running interpro.sh and it finished running. in the output of that I need to run the filter_interpro_output.sh which will filter out empty files and files that contains transposonn elements 

current workflow of scripts: 



0.download_accessions.sh 
1.busco.sh
1b.busco_summaries.sh 
2.orthofinder.sh 
2b.prepare_interproinput.sh
3.interproscan.sh
3b.filter_interpro_output.sh
4.make_upset_filtered.R  
5.config.sh
6: run treepl with treePL config.cfg
7a.topgoprep.sh 
7b.topgo.R

#then I have the topgo analysis R script in "C:/Users/vicme/OneDrive/Desktop/Manuscript_Pawpaw" 
