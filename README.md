# NarendranLab

Narendran Laboratory, Arnie Charbonneau Cancer Institute, Alberta Children's Hospital and University of Calgary, Calgary, Alberta, Canada

Corresponding author email: anarendr@ucalgary.ca

This repository contains R code for the curation of RNAseq data from TARGET pediatric tumour tissue samples provided by the TCGAbiolinks Bioconductor/R package, and GTEX adult normal tissue samples provided by the recount2 project through the TCGAbiolinks Bioconductor/R package. Differential expression analysis (DEA) is performed using the limma-voom method with cut-offs at fold change > 2 or < -2 and p-value < 10e-16. Gene ontology (GO) enrichment analysis was performed using the topGO Bioconductor/R package. Project codes matching corresponding tumour and normal tissues are provided in Project_Codes_New.xlsx. PID-related genes were curated from the Fulgent Genetics 'Comprehensive Primary Immunodeficiency NGS Panel' and is provided in Panel_PIDGenes.xlsx. 

This repository also contains R code for the survival analysis of the experimentally identified differentially expressed genes in pediatric cancers from the TARGET database provided by the TCGAbiolinks Bioconductor/R package. Cox regression for survival analysis was performed and Kaplan-Meier curves were generated. The genes used are outlined in Pediatric_DEGs.csv.

Standing S, Tran S, Murguia-Favela L, Kovalchuk O, Bose P, Narendran A. Identification of Altered Primary Immunodeficiency-Associated Genes and Their Implications in Pediatric Cancers. Cancers (Basel). 2022 Nov 30;14(23):5942. doi: 10.3390/cancers14235942. PMID: 36497424; PMCID: PMC9741011.
