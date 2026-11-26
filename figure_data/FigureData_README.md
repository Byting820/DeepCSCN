# Processed Data for Manuscript Figures

This directory contains processed data used to generate the data-driven panels in the main and supplementary figures of the DeepCSCN manuscript. The tables below map each figure or panel to its corresponding data file(s).

## File Formats and Organization

- All `.tsv` files are UTF-8-encoded, tab-delimited tables with a header row.
- For files covering multiple panels, use the `panel` column where applicable, or consult the panel descriptions below.
- These are processed plotting data, not raw sequencing data. Original data sources and accession numbers are described in the manuscript's Data Availability section and Supplementary Table S2.

**Note:** The figure numbers and panel labels correspond to the revised manuscript.

## Main Figures

| Figure | Processed data | Relationship to the figure |
| --- | --- | --- |
| Fig. 2b | [`Fig2/Fig2b.tsv`](Fig2/Fig2b.tsv) | Gene-level values used for the hESC module plot. |
| Fig. 2c–f | [`four_module_validation_metrics_plot_data.tsv`](Fig2/Fig2c-f/four_module_validation_metrics_plot_data.tsv) | Method-level values plotted in panels c–f: GO-significant modules, state-elevated modules, COSG marker-enriched modules, and marker proportion. |
| Fig. 2c–f, supporting data | [`three_primary_metrics_by_method_time.tsv`](Fig2/Fig2c-f/three_primary_metrics_by_method_time.tsv), [`three_primary_metrics_summary.tsv`](Fig2/Fig2c-f/three_primary_metrics_summary.tsv), and [`three_primary_metrics_metadata.json`](Fig2/Fig2c-f/three_primary_metrics_metadata.json) | Per-time-point values, overall summaries, and metric definitions underlying panels c–f. |
| Fig. 3a | [`Fig3/Fig3a.tsv`](Fig3/Fig3a.tsv) | Mean and SD of AUROC and AUPRC ratio for the model-comparison panel. |
| Fig. 3b | [`Fig3/Fig3b.tsv`](Fig3/Fig3b.tsv) | Mean and SD of AUROC and AUPRC ratio for the cell-order permutation experiment. |
| Fig. 3c | [`Fig3/Fig3c.tsv`](Fig3/Fig3c.tsv) | Module gene retention and marker-gene proportions by PBMC cell type under different cell-ordering conditions. |
| Fig. 4b–e | [`Fig4b.tsv`](Fig4/Fig4b.tsv), [`Fig4c.tsv`](Fig4/Fig4c.tsv), [`Fig4d.tsv`](Fig4/Fig4d.tsv), [`Fig4e.tsv`](Fig4/Fig4e.tsv) | Data for the PCA plot, expression heatmap, gene violin plots, and CDH2/PRTG expression plots, respectively. |
| Fig. 5a–e | [`Fig5a.tsv`](Fig5/Fig5a.tsv), [`Fig5b.tsv`](Fig5/Fig5b.tsv), [`Fig5c.tsv`](Fig5/Fig5c.tsv), [`Fig5d.tsv`](Fig5/Fig5d.tsv), [`Fig5e.tsv`](Fig5/Fig5e.tsv) | Data for the PBMC UMAP, top-gene heatmap, GO heatmap, B-cell network nodes/layout, and CD4+ T-cell network nodes/layout, respectively. |
| Fig. 6a–e | [`Fig6a.tsv`](Fig6/Fig6a.tsv), [`Fig6b.tsv`](Fig6/Fig6b.tsv), [`Fig6c.tsv`](Fig6/Fig6c.tsv), [`Fig6d.tsv`](Fig6/Fig6d.tsv), [`Fig6e.tsv`](Fig6/Fig6e.tsv) | Data for the spatial cell-type plot, expression heatmap, GO plot, cancer-network nodes/layout, and normal-network nodes/layout, respectively. |

## Supplementary Figures

| Figure | Processed data | Relationship to the figure |
| --- | --- | --- |
| Supplementary Fig. S3 | [`FigS3.tsv`](Supplementary_Figures/FigS3.tsv) | Component-ablation results (Supplementary Fig. S3b). |
| Supplementary Fig. S4 | [`FigS4.tsv`](Supplementary_Figures/FigS4.tsv) | Feature-depth and rare-cell sensitivity results for all nine panels. |
| Supplementary Fig. S5 | [`FigS5.tsv`](Supplementary_Figures/FigS5.tsv) | AUROC and AUPRC-ratio heatmap values for ten methods at four sparsity levels. |
| Supplementary Fig. S6a–f | [`FigS6.tsv`](Supplementary_Figures/FigS6.tsv) | Module-specific data for hHep, mDC, mESC, mHSC_E, mHSC_GM, and mHSC_L. |
| Supplementary Fig. S7 | [`FigS7.tsv`](Supplementary_Figures/FigS7.tsv) | Module-function benchmark values across six datasets. |
| Supplementary Fig. S8 | [`FigS8.tsv`](Supplementary_Figures/FigS8.tsv) | AUROC, AUPRC ratio, EPR, partial AUPR, and precision@k heatmap values. |
| Supplementary Fig. S9a–f | [`FigS9.tsv`](Supplementary_Figures/FigS9.tsv) | Network node and module data for hESC at 00, 12, 24, 36, 72, and 96 h. |
| Supplementary Fig. S10 | [`FigS10.tsv`](Supplementary_Figures/FigS10.tsv) | B-cell GO-enrichment results. |
| Supplementary Fig. S11a–c | [`FigS11.tsv`](Supplementary_Figures/FigS11.tsv) | Network node and layout data for Cytotoxic T cells, CD14+ monocytes, and CD16+ monocytes. |
| Supplementary Fig. S12 | [`FigS12.tsv`](Supplementary_Figures/FigS12.tsv) | Spatial spot and panel data for cancer epithelial module 1. |
| Supplementary Fig. S13 | [`FigS13.tsv`](Supplementary_Figures/FigS13.tsv) | Spatial spot and panel data for cancer epithelial module 2. |
| Supplementary Fig. S14 | [`FigS14.tsv`](Supplementary_Figures/FigS14.tsv) | Spatial spot and panel data for cancer epithelial module 3. |
| Supplementary Fig. S15a–e | [`FigS15.tsv`](Supplementary_Figures/FigS15.tsv) | Runtime and memory-scaling values. |
| Supplementary Fig. S16d–e | [`FigS16.tsv`](Supplementary_Figures/FigS16.tsv) | Training-loss and NMI values by epoch. |
| Supplementary Fig. S17 | [`FigS17.tsv`](Supplementary_Figures/FigS17.tsv) | AUROC and AUPRC ratio values for ten methods evaluated on the PBMC dataset. |


## Figures Without Separate Processed Data

Some figures are schematic illustrations or model architecture diagrams rather than data-driven plots. No separate processed numerical data are required for these figures:

| Figure | Description | Processed data |
| --- | --- | --- |
| Fig. 2a | Schematic illustration of the DeepCSCN framework. | Not applicable. |
| Fig. 4a | Schematic illustration of the analysis workflow. | Not applicable. |

The same principle applies to supplementary panels that consist solely of schematic illustrations or model architecture diagrams. Processed data for all data-driven supplementary panels included in this directory are mapped above.

## Using the Figure Data

To locate the data for a figure, find its number in the corresponding table and follow the relative file link. Files contain the processed values or plotting inputs associated with the indicated panel. The original datasets and accession numbers are documented in the manuscript's Data Availability section and Supplementary Table S2.

For quantitative panels based on replicate analyses, summary statistics in these files should be interpreted according to the corresponding figure legends and Supplementary Methods. Reproducing individual plotted replicate points requires the underlying run-level values when those are shown.
