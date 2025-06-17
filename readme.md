
# LAMPS
This repository contains code for <u>L</u>U<u>A</u>D <u>M</u>ulti-<u>P</u>latform <u>S</u>equencing data analysis. Scripts are organized according to figure order and example data was supplied in input directory.

## Input files
Input files will be available upon publication.

## Scripts
Scripts|Data|Description
-----------|------------|-------------------------------------
01_fig1_cell_qc&annotation.r|LAMPS single-cell & LSA single-cell|Quality control and cell type annotation for scRNA-seq data 
02_fig1_malignant_cell_identification.r|LAMPS single-cell & LSA single-cell|Identify malignant cell from epithelial cell cluster according to marker expression, GII and cluster distribution 
03_fig1_mp_discovery.r|LAMPS single-cell & LSA single-cell|Meta-program identification from scRNA-seq data 
04_fig1_scoring&qc.r|LAMPS single-cell & LSA single-cell|Score malignant cells with meta-programs and re-cluster for quality control  
05_fig1_trajectory(individual).r|LAMPS single-cell & LSA single-cell| Construct malignant cell trajectory for each patient separately  
06_fig2_trajectory(unified).r|LAMPS single-cell & LSA single-cell| Construct the consensus trajectory  
07_fig2_bulk_assignment.r|LSA bulk| Assign patients with only bulk-seq data to a branch 
08_fig2_multi-region_transition_inference.r|LAMPS multi-region| Branch transition direction inference 
09_fig3_TME_clustering.r|LSA bulk| Cluster patients according to TME cell composition for TME subtype discovery 
10_fig3_EM_discovery.r|LSA bulk & LSA single-cell| Eco-module discovery 
11_fig4_ST.r|LAMPS ST| Analysis for spatial transcriptomic data including spot deconvolution, cell colocalization and cell interaction 
12_fig5_cohort_comparision.r|LSA bulk| Cohort comparsion between EAS and EUR in malignant cell and TME cell 
13_fig6_prognosis.r|LSA bulk| COX model construction 
14_fig6_treatment.r|LSA bulk & SU2C-MARK| Treatment prediction 
  
## System

```         
LSB Version:    :core-4.1-amd64:core-4.1-noarch:cxx-4.1-amd64:cxx-4.1-noarch:desktop-4.1-amd64:desktop-4.1-noarch:languages-4.1-amd64:languages-4.1-noarch:printing-4.1-amd64:printing-4.1-noarch
Distributor ID: CentOS
Description:    CentOS Linux release 7.5.1804 (Core) 
Release:        7.5.1804
Codename:       Core
```
