# iMotifEvolution
Code and core analysis data for the study: **i-Motifs are shaped by positive selection and drive human-specific regulatory innovation and disease risk**

## Overview
This repository contains the reproducible analysis code and processed small-scale data for the evolutionary and functional analysis of i-Motifs across the tree of life, mammals, and great apes.

## Repository Structure
- **01_Macro_Evolution**: Random GC simulation, robust Z-score calculation, phylogenetic tree visualization (Figure 1)
- **02_Mammal_Conservation**: i-Motif conservation classification (C1/C2/C3), stability scoring, and GO/KEGG enrichment (Figure 2)
- **03_Great_Ape_Comparative**: Great ape genome comparison, chromosome-level density heatmap (Figure 3)
- **04_Functional_Enrichment**: UpSet plot, promoter enrichment, and functional/diabetes pathway analysis (Figure 4)
- **05_TSS_Expression**: Strand-specific TSS profiling and gene expression (TPM) analysis (Figure 5)

## Software Versions
- R 4.5.2
- Python 3.11.10

## Data Availability
The core i-Motif coordinates, classification annotations, and enrichment results generated in this study are available in this repository.
The number of prokaryotic and eukaryotic genomes used, as well as the version information of great ape T2T genome assemblies, are provided in Supplementary Tables 1 and 2.
RNA-seq data were obtained from GEO under accession GSE105160.
Public genome assemblies are available from RefSeq and NCBI Assembly.
Zoonomia multiple genome alignments are available from the UCSC Genome Browser.

## License
MIT License
