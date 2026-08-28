IRI_Niche_2025

# Code supporting:

Rietjens RGJ, Manzato B, Liao Y, Sol WMPJ, Mahfouz A, Giera M, van den Berg BM, Dumas SJ, Wang G, Rabelink TJ.
Impaired Metabolic Remodeling Characterizes Injured Niches during Post-ischemic Kidney Repair.
Journal of the American Society of Nephrology. 2026.
https://doi.org/10.1681/ASN.0000001216

<img width="1259" height="772" alt="Screenshot 2026-08-28 at 08 54 15" src="https://github.com/user-attachments/assets/5d5f3771-b07a-4741-900d-adc5442fd031" />

# Overview

This repository contains the analysis code supporting the multimodal spatial omics analyses described in the manuscript.

The study integrates spatial metabolomics using quantitative mass spectrometry imaging (qMSI) with high-resolution spatial transcriptomics (Stereo-seq) to identify metabolically distinct tissue niches following ischemia-reperfusion injury and to characterize metabolic alterations in proximal tubule cells residing within injured microenvironments.

The computational workflow includes spatial domain identification using BANKSY, spatial metabolite analysis, projection of metabolomics-derived niches onto Stereo-seq data, and cross-modal validation of cell-type composition.

## Analysis workflow

<img width="458" height="405" alt="Screenshot 2026-08-28 at 08 50 34" src="https://github.com/user-attachments/assets/cff2bd13-7d67-4445-b4cc-ebee603f9151" />

## Repository structure

**spatial_metabolomics_qmsi_analysis/**
Analysis of the spatial metabolomics data.
* sqMSI_imputation.ipynb — spatial imputation of missing qMSI measurements.
* banksy-multisample-qmsi.R — multi-sample BANKSY analysis of qMSI data and identification of spatial domains.
* scils_like_visualization.R — visualization of metabolite distributions within defined niches.
* volcano3D_metabolomics.R — statistical comparison and visualization of metabolite differences between proximal tubule populations across niches.

**spatial_transcriptomics_stereoseq_analysis/**
Processing and integration of Stereo-seq data.
* tissue_segmentation_stseq.ipynb — segmentation of individual kidneys from Stereo-seq datasets.
* transcriptomics_overlay.R — projection of qMSI/BANKSY-derived spatial domains onto Stereo-seq coordinates.

**validation/**
Cross-modal validation of spatial niche assignments.
* knn_matching.ipynb — k-nearest-neighbour matching of qMSI spatial spots to Stereo-seq annotations to assess cell-type composition of metabolic spatial domains.
* overlaying_tissues.ipynb — tissue-level spatial overlay/registration analyses.

The code is provided to document and support the analyses reported in the manuscript.

## Citation
If you use this code, please cite:

Rietjens RGJ, Manzato B, Liao Y, Sol WMPJ, Mahfouz A, Giera M, van den Berg BM, Dumas SJ, Wang G, Rabelink TJ. Impaired Metabolic Remodeling Characterizes Injured Niches during Post-ischemic Kidney Repair. Journal of the American Society of Nephrology. 2026. doi:10.1681/ASN.0000001216


