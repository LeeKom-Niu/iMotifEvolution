# 06_Selection_Analysis

Selection analysis of i-Motif stem and loop regions using the HKA test and PhyloFit
branch-length estimation.

Companion to Supplementary Figures S13–S14.

## Overview

Two independent pipelines run on promoter regions for three great-ape groups
(hominid, homininae, hominini):

1. **HKA arm** — partitions promoter i-Motifs into stem, loop and 4-fold degenerate
   (4FD) site classes, counts polymorphic sites and fixed divergence, and applies
   Fisher's exact test to each pairwise comparison.
2. **PhyloFit arm** — builds MAFs for the same partitions, fits a phylogenetic model
   of sequence evolution per partition, and compares the human branch length of stems
   against that of loops.

## Files

| File | Purpose | Input | Output |
|---|---|---|---|
| `run_HKA_promoter.sh` | Orchestrator for the HKA arm; runs 01 → 02 → 03 | — | — |
| `01_prepare_hka_counts.sh` | Extract promoter G4s and their stem/loop partitions; count SNPs and divergence; define 4FD sites | `input/stems.bed`, `input/loops.bed`, `input/all_divergence_genomewide.bed`, `input/hprc.snps.t2t.bed`, `input/4fd_all.bed`, `input/groups_0based/human.*G4s_0based.bed`, `promoter_regions.bed.gz` | `output/counts_promoter.tsv`, `output/4fd_counts.txt`, intermediate BEDs |
| `02_fisher_test_hka.R` | Fisher's exact test for stem vs loop, stem vs 4FD and loop vs 4FD | `output/counts_promoter.tsv` | `output/fisher_results_promoter.tsv` |
| `03_plot_hka_combined.R` | Three-panel odds-ratio figure with confidence intervals | `output/fisher_results_promoter.tsv` | `output/Figure_HKA_combined.pdf`, `.png` |
| `run_PhyloFit_promoter.sh` | Orchestrator for the PhyloFit arm; runs 04 → 05 | — | — |
| `04_run_phylofit_models.sh` | Generate MAFs with hal2maf, filter them, fit PhyloFit models, extract the hs1 branch length | `hal_data/*.hal`, `hal_input/{group}_{region}_{part}.bed`, `filter_maf_first.py` | `maf_output/*.maf`, `phylofit_models/*.mod`, `branch_lengths.tsv` |
| `05_plot_phylofit_promoter.py` | Stem/loop branch-length ratio per group, promoter region | `branch_lengths.tsv` | `Figure_phylofit_promoter.pdf`, `.png` |
| `filter_maf_first.py` | Helper invoked by 04: deduplicate MAF blocks and simplify sequence names | Raw MAF | Filtered MAF |

## Run order

```
run_HKA_promoter.sh          # HKA arm:      01 -> 02 -> 03
run_PhyloFit_promoter.sh     # PhyloFit arm: 04 -> 05
```

The two arms are independent and may be run in either order. On a cluster,
`run_PhyloFit_promoter.sh` submits `04_run_phylofit_models.sh` via `sbatch --wait`
so that the plotting step runs only after the models and `branch_lengths.tsv` exist.

## Dependencies

- R with `data.table`, `ggplot2`, `patchwork`
- Python 3 with `pandas`, `matplotlib`, `numpy`
- `phast` (for `phyloFit`)
- `hal2maf` (Cactus)
- `bedtools`

## Data

`data/` contains the small result tables produced by these scripts
(`branch_lengths.tsv`, `fisher_results_promoter.tsv`). The large inputs
(HAL alignment, SNP and divergence BEDs, promoter regions) are not distributed here;
see the main README for data availability.
