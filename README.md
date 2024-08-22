## Summary

This repository contains the scripts yielding to the processed scRNAseq data of NK cells in Multiple Myeloma patients.
The study utilizing this dataset is published in   the results in our [article](https://doi.org/10.1182/blood.2023023529).

The processed data of 8 NK cell clusters determined in the MM patients and healthy donors were stored in Zenodo as an accessible Seurat object [link](prospective_link_here)

Below structure of scripts are the step by step data curation and then investigation of NK cells in the downstream analysis.

## Structural layout
```{r eval=FALSE, include=TRUE}
├── src
│   └── 01_individually_processing_per_emulsion.R  # Reading Cellranger / Velocyto outputs, then annotating cell types per emulsion to ref PBMC in a loop
│   └── 02_merging_processed_emulsions.R #Merging each emulsions' data into one data object
│   └── 03a_nuc_frac_calcul.R #Calculating nuclear RNA fraction of cells
│   └── 03b_nuclear_RNA_fraction_metric.py # Helper python script of nuclear RNA fraction
│   └── 04_QC_check_and_trim.R # QC plots and trimming based-on QC
│   └── 05a_preprocessing.R # Normalization,ScaleData,Dimensional Reduction,PC-correction etc.
│   └── 05b_plots_leading_to_nk_subsetting.R # Intermediate step of assessing inferred cell types
│   └── 06a_nk_subsetting.R # Separating NK cells from the rest of the cell types 
│   └── 06b_removal_of_prolif_nk_portion.R # Keeping non-proliferating, more stable NK cells for downstream analysis
│   └── multi_criteria_filt_GSEA.R # Multiple parameter assessment oriented GSEA script
│   └── combin_presto_function.R # Helper function of combinatorial Presto running to rank all genes for comparison groups (ClusterXvsY)
│   └── combin_gsea_function.R # Helper function to run fGSEA for comparison pairs and keep them in list of lists format
│   └── # Placeholder for "Slingshot-custom-script"
│   └── # Placeholder for "scComp" script
│   └── # Placeholder for "MELD" script
```
## Analysis diagram

[AnalysisDiagram.pdf](https://github.com/user-attachments/files/16709177/AnalysisDiagram.pdf)

## Analysis medium

All of the analyses of this study except *MELD* were performed with the local machine (workstation) containing Intel® Core™ i9-9900K CPU and 128 GB RAM, running on Linux (Ubuntu 22.04.2 LTS) operation system.
