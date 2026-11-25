# Processed Data for Manuscript Figures

This directory contains the processed data used to generate the
data-driven panels in the main and supplementary figures of the DeepCSCN
manuscript. The tables below map each figure or figure panel to the
corresponding processed data file(s) and briefly describe their
contents.

## File formats and organization

-   All `.tsv` files are UTF-8 encoded, tab-delimited text files with a
    header row.
-   For a file used to generate multiple panels, the relevant panel or
    data subset is specified in the description or identified by a
    `panel` column, where applicable.
-   This directory contains processed figure-generation data, not raw
    sequencing files. Original dataset sources and accession numbers are
    provided in the manuscript's Data Availability section and
    Supplementary Table S2.

## Main figures

  --------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  Figure                  Processed data file(s)                                                                                                     Relationship to the
                                                                                                                                                     figure
  ----------------------- -------------------------------------------------------------------------------------------------------------------------- -----------------------
  Fig. 2b                 [`Fig2/Fig2b.tsv`](Fig2/Fig2b.tsv)                                                                                         Gene-level values used
                                                                                                                                                     to generate the hESC
                                                                                                                                                     module plot.

  Fig. 2c--f              [`Fig2/Fig2c-f/four_module_validation_metrics_plot_data.tsv`](Fig2/Fig2c-f/four_module_validation_metrics_plot_data.tsv)   Method-level values
                                                                                                                                                     plotted in panels c--f,
                                                                                                                                                     including
                                                                                                                                                     GO-significant modules,
                                                                                                                                                     state-elevated modules,
                                                                                                                                                     marker-enriched
                                                                                                                                                     modules, and
                                                                                                                                                     marker-gene
                                                                                                                                                     proportions.

  Fig. 2c--f (supporting  [`Fig2/Fig2c-f/three_primary_metrics_by_method_time.tsv`](Fig2/Fig2c-f/three_primary_metrics_by_method_time.tsv),          Per-time-point values,
  data)                   [`Fig2/Fig2c-f/three_primary_metrics_summary.tsv`](Fig2/Fig2c-f/three_primary_metrics_summary.tsv), and                    overall summaries, and
                          [`Fig2/Fig2c-f/three_primary_metrics_metadata.json`](Fig2/Fig2c-f/three_primary_metrics_metadata.json)                     metric definitions
                                                                                                                                                     underlying panels c--f.

  Fig. 3a                 [`Fig3/Fig3a.tsv`](Fig3/Fig3a.tsv)                                                                                         Mean and standard
                                                                                                                                                     deviation of AUROC and
                                                                                                                                                     AUPRC ratio for the
                                                                                                                                                     architecture/model
                                                                                                                                                     comparison.

  Fig. 3b                 [`Fig3/Fig3b.tsv`](Fig3/Fig3b.tsv)                                                                                         Mean and standard
                                                                                                                                                     deviation of AUROC and
                                                                                                                                                     AUPRC ratio under the
                                                                                                                                                     original,
                                                                                                                                                     within-cell-type, and
                                                                                                                                                     global cell-ordering
                                                                                                                                                     conditions.

  Fig. 3c                 [`Fig3/Fig3c.tsv`](Fig3/Fig3c.tsv)                                                                                         Module-gene retention
                                                                                                                                                     and marker-gene
                                                                                                                                                     proportions by PBMC
                                                                                                                                                     cell type under the
                                                                                                                                                     cell-ordering
                                                                                                                                                     permutation conditions.

  Fig. 4b--e              [`Fig4b.tsv`](Fig4b.tsv), [`Fig4c.tsv`](Fig4c.tsv), [`Fig4d.tsv`](Fig4d.tsv), and [`Fig4e.tsv`](Fig4e.tsv)                 Data used for the PCA
                                                                                                                                                     plot, expression
                                                                                                                                                     heatmap, gene violin
                                                                                                                                                     plots, and CDH2/PRTG
                                                                                                                                                     expression plots,
                                                                                                                                                     respectively.

  Fig. 5a--e              [`Fig5a.tsv`](Fig5a.tsv), [`Fig5b.tsv`](Fig5b.tsv), [`Fig5c.tsv`](Fig5c.tsv), [`Fig5d.tsv`](Fig5d.tsv), and                Data used for the PBMC
                          [`Fig5e.tsv`](Fig5e.tsv)                                                                                                   UMAP, top-gene heatmap,
                                                                                                                                                     GO-enrichment heatmap,
                                                                                                                                                     B-cell network, and
                                                                                                                                                     CD4+ T-cell network,
                                                                                                                                                     respectively.

  Fig. 6a--e              [`Fig6a.tsv`](Fig6a.tsv), [`Fig6b.tsv`](Fig6b.tsv), [`Fig6c.tsv`](Fig6c.tsv), [`Fig6d.tsv`](Fig6d.tsv), and                Data used for the
                          [`Fig6e.tsv`](Fig6e.tsv)                                                                                                   spatial cell-type plot,
                                                                                                                                                     expression heatmap,
                                                                                                                                                     GO-enrichment plot,
                                                                                                                                                     cancer epithelial
                                                                                                                                                     network, and normal
                                                                                                                                                     epithelial network,
                                                                                                                                                     respectively.
  --------------------------------------------------------------------------------------------------------------------------------------------------------------------------

## Supplementary figures

  ------------------------------------------------------------------------------------------------------------------------
  Figure                  Processed data file(s)                                                   Relationship to the
                                                                                                   figure
  ----------------------- ------------------------------------------------------------------------ -----------------------
  Supplementary Fig. S3   [`Supplementary_Figures/FigS3.tsv`](Supplementary_Figures/FigS3.tsv)     Component-ablation
                                                                                                   results shown in panel
                                                                                                   b.

  Supplementary Fig. S4   [`Supplementary_Figures/FigS4.tsv`](Supplementary_Figures/FigS4.tsv)     Feature-depth and
                                                                                                   rare-cell sensitivity
                                                                                                   results for all nine
                                                                                                   panels.

  Supplementary Fig. S5   [`Supplementary_Figures/FigS5.tsv`](Supplementary_Figures/FigS5.tsv)     AUROC and AUPRC-ratio
                                                                                                   heatmap values for ten
                                                                                                   methods across four
                                                                                                   sparsity levels.

  Supplementary Fig.      [`Supplementary_Figures/FigS6.tsv`](Supplementary_Figures/FigS6.tsv)     Module-specific data
  S6a--f                                                                                           for hHep, mDC, mESC,
                                                                                                   mHSC_E, mHSC_GM, and
                                                                                                   mHSC_L.

  Supplementary Fig. S7   [`Supplementary_Figures/FigS7.tsv`](Supplementary_Figures/FigS7.tsv)     Module-function
                                                                                                   benchmark values across
                                                                                                   six datasets.

  Supplementary Fig. S8   [`Supplementary_Figures/FigS8.tsv`](Supplementary_Figures/FigS8.tsv)     AUROC, AUPRC ratio,
                                                                                                   EPR, partial AUPR, and
                                                                                                   Precision@k values
                                                                                                   across benchmark
                                                                                                   datasets and methods.

  Supplementary Fig.      [`Supplementary_Figures/FigS9.tsv`](Supplementary_Figures/FigS9.tsv)     Network-node and module
  S9a--f                                                                                           data for hESC at 00,
                                                                                                   12, 24, 36, 72, and 96
                                                                                                   h.

  Supplementary Fig. S10  [`Supplementary_Figures/FigS10.tsv`](Supplementary_Figures/FigS10.tsv)   B-cell GO-enrichment
                                                                                                   results.

  Supplementary Fig.      [`Supplementary_Figures/FigS11.tsv`](Supplementary_Figures/FigS11.tsv)   Network-node and layout
  S11a--c                                                                                          data for cytotoxic T
                                                                                                   cells, CD14+ monocytes,
                                                                                                   and CD16+ monocytes.

  Supplementary Fig. S12  [`Supplementary_Figures/FigS12.tsv`](Supplementary_Figures/FigS12.tsv)   Spatial-spot and panel
                                                                                                   data for cancer
                                                                                                   epithelial module 1.

  Supplementary Fig. S13  [`Supplementary_Figures/FigS13.tsv`](Supplementary_Figures/FigS13.tsv)   Spatial-spot and panel
                                                                                                   data for cancer
                                                                                                   epithelial module 2.

  Supplementary Fig. S14  [`Supplementary_Figures/FigS14.tsv`](Supplementary_Figures/FigS14.tsv)   Spatial-spot and panel
                                                                                                   data for cancer
                                                                                                   epithelial module 3.

  Supplementary Fig.      [`Supplementary_Figures/FigS15.tsv`](Supplementary_Figures/FigS15.tsv)   Runtime and
  S15a--e                                                                                          memory-scaling values.

  Supplementary Fig.      [`Supplementary_Figures/FigS16.tsv`](Supplementary_Figures/FigS16.tsv)   Training-loss and NMI
  S16d--e                                                                                          values across training
                                                                                                   epochs.

  Supplementary Fig. S17  [`Supplementary_Figures/FigS17.tsv`](Supplementary_Figures/FigS17.tsv)   AUROC and AUPRC-ratio
                                                                                                   benchmark values for
                                                                                                   ten methods evaluated
                                                                                                   on the PBMC dataset.
  ------------------------------------------------------------------------------------------------------------------------

## Figures without separate processed data files

Some figures are schematic illustrations or model/framework diagrams
rather than plots generated from processed numerical data. Therefore, no
separate processed data files are applicable to these figures:

-   **Fig. 2a:** schematic illustration of the DeepCSCN framework.
-   **Fig. 4a:** schematic illustration of the analysis workflow.
-   **Supplementary figures consisting solely of schematic illustrations
    or model-architecture diagrams:** no separate processed data files
    are required.

For supplementary figures containing data-driven panels, the
corresponding processed data files are listed in the table above.

## How to use these files

Each processed data file is linked to its corresponding figure or panel
in the tables above. The files contain the numerical values or
intermediate plotting data used to generate the corresponding
visualizations. To locate the data for a particular figure, identify the
figure number and panel in the relevant table and follow the listed
relative file path.

For information on the original datasets, accession numbers, and data
sources, please refer to the Data Availability section and Supplementary
Table S2 of the manuscript.
