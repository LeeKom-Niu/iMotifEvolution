# 07_Functional_Validation

Functional validation of predicted i-Motifs: overlap with ATAC-seq open chromatin and
expression of iM-associated genes.

Companion to Supplementary Figures S21–S22.

## Overview

Two lines of evidence are assembled here:

1. **Chromatin accessibility** — i-Motif peaks are intersected with ATAC-seq peaks
   within TSS-proximal promoters, and the resulting overlaps are summarised as
   two-set Venn diagrams.
2. **Expression** — genes whose promoters carry a unique i-Motif are compared with
   genes whose promoters carry a unique G4 and with all expressed genes, using
   HEK293T TPM values.

## Files

| File | Purpose | Input | Output |
|---|---|---|---|
| `run_analysis_2.sh` | T2T pipeline: build TSS ± 1 kb promoters, filter iM peaks at score >= 0.5, intersect with ATAC peaks, write summary statistics and draw a Venn diagram | `Homo_sapiens_all.bed`, `promoter_regions.bed`, ATAC peak file | `results_v2/` (`im_promoter.bed`, `im_promoter_atac.bed`, `statistics.txt`, `plot_venn.py`, `venn_im_promoter_vs_atac_T2T_v2.pdf`, `.png`) |
| `venn_diagram.py` | Single-panel Venn diagram of iM promoter peaks vs ATAC peaks; counts are read at run time with `wc -l` | `im_promoter.bed`, `GSE260781_ATAC.bed3`, `im_promoter_atac.bed` | `venn_im_promoter_vs_atac.pdf`, `.png` |
| `venn_combine.py` | Two-panel combined Venn figure: panel A, iMab CUT&Tag vs ATAC (HEK293T); panel B, predicted iM vs ATAC (T2T, score >= 0.5). Counts are fixed in the script | — | `venn_combined.pdf`, `.png` |
| `plot_promoter_expression_v2.py` | Violin plot with jittered scatter comparing TPM distributions of iM-associated, G4-associated and all genes. Outliers are removed by the 1.5 x IQR rule; significance is computed by a two-sided Mann-Whitney U test | `data/HEK293T_gene_TPM.tsv`, `data/iM_unique_promoter_overlap.bed`, `data/G4_unique_promoter_overlap.bed` | `promoter_v2_unique.pdf`, `.png` |

## Data sources

- iMab CUT&Tag peaks: Zanin et al. 2023
- ATAC-seq: GSE260781, GSE274771
- Gene expression (TPM): HEK293T, three replicates

## Dependencies

- Python 3 with `pandas`, `numpy`, `matplotlib`, `scipy`
- `bedtools`

## Notes

`plot_promoter_expression_v2.py` opens its input files by bare filename and therefore
expects to be run with `data/` as the working directory. `run_analysis_2.sh` and
`venn_diagram.py` reference inputs that were produced on the original compute cluster
and are not all distributed with this repository.
