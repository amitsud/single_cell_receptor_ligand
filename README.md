# Tabula Sapiens & Multiple Myeloma Microenvironment Analysis

This repository contains R code for analyzing cell-cell interactions and ligand-target predictions in the context of plasma cell dyscrasias, as detailed in Went et al. Proteome-wide Mendelian Randomisation Identifies Protein Associations and Therapeutic Targets for B-cell Malignancy. Blood Neoplasia. 2026.

## Analysis Overview

The analysis is performed in two primary stages:

1.  **Tabula Sapiens Processing & CellChat Analysis:**
    * Loads and processes single-cell RNA-seq data from the Tabula Sapiens atlas.
    * Filters for immune-relevant tissues (Blood, Bone Marrow, Lymph Node, Spleen, Thymus).
    * Integrates external results of proteome-wide Mendelian randomisation of B-cell tumors to filter the CellChat interaction database.
    * Performs CellChat analysis to infer cell-cell communication networks, focusing on specific immune interactions.

2.  **NicheNet Analysis:**
    * Focuses on Plasma cells as the receiver population.
    * Investigates specific ligands of interest (`TNFSF13`/APRIL and `TNFSF13B`/BAFF).
    * Utilizes NicheNet to predict the regulatory potential of these ligands on a Multiple Myeloma/MGUS gene signature derived from differential expression analysis.

## Prerequisites

This code relies on the following R packages:

* **Seurat**
* **CellChat**
* **nichenetr**
* **tidyverse** (dplyr, ggplot2, vroom)
* **future** (for parallel processing)

You can install the main dependencies via CRAN or GitHub (for CellChat/NicheNet/SeuratWrappers).

## Data Availability

The analysis requires the following inputs:

1.  **Tabula Sapiens Object:** A Seurat object containing the Tabula Sapiens dataset.
    * *Source:* [Tabula Sapiens on Figshare]([https://figshare.com/articles/dataset/Tabula_Sapiens_v2/27921984](https://figshare.com/articles/dataset/Tabula_Sapiens_release_1_0/14267219/1))
2.  **NicheNet Networks:** The code automatically downloads the necessary ligand-receptor networks and weighted matrices from Zenodo.
3.  **Differential Expression Data:** A text file containing differentially expressed genes (e.g., `mm_scrnaseq_deg.txt`). This file should include columns for `disease_stage`, `gene`, and adjusted p-values.

## Usage

1.  Clone this repository.
2.  Open the R script.
3.  **Update File Paths:** specific lines marked with `PLACEHOLDER` in the code must be updated to point to your local data directories:
4.  Run the script interactively or via command line.
