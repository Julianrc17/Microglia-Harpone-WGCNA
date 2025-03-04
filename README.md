# README

## Overview
This project performs Weighted Gene Co-expression Network Analysis (WGCNA) on microglia gene expression data. The script processes expression and trait data, constructs gene co-expression networks, identifies modules, and correlates them with genotype traits.

## Prerequisites
Ensure you have the following dependencies installed in R before running the script:

```r
install.packages("WGCNA")
install.packages("gplots")
```

Additionally, set the working directory to the folder containing the data files:

```r
workingDir = "."
setwd(workingDir)
```

## Data Files
- `combatExpr.txt`: Expression data after batch correction.
- `joint_samplesdf.txt`: Sample metadata containing genotype and cell type.
- `GO_screen_annot.txt`: Annotation file for gene symbols.
- `Microglia-01-dataInput.RData`: Processed expression data.
- `Microglia-02-networkConstruction-auto.RData`: Network construction results.

## Workflow
### 1. Load and Filter Data
- Read expression data (`combatExpr.txt`).
- Read sample metadata (`joint_samplesdf.txt`) and filter for relevant genotypes (`E2WT`, `E3WT`, `E4WT`, `KO`).
- Extract microglia-specific samples and subset expression data.
- Merge expression data with annotation.

### 2. Preprocess Data for WGCNA
- Convert data to numeric format.
- Remove genes and samples with excessive missing values.
- Check for outliers using hierarchical clustering.

### 3. Construct Network and Identify Modules
- Choose soft-thresholding power.
- Construct network and identify gene modules.
- Merge similar modules.

### 4. Correlate Modules with Traits
- Correlate module eigengenes with genotype traits (`KO`, `E2WT`, `E3WT`, `E4WT`).
- Identify significant modules.
- Generate heatmaps and scatter plots.

### 5. Save and Export Results
- Save module correlations (`module_correlations_KO.txt`, `module_correlations_E2WT.txt`, etc.).
- Export gene significance values for each genotype.

## Outputs
- `module_correlations_KO.txt`: Correlation values for KO genotype.
- `module_correlations_E2WT.txt`: Correlation values for E2WT genotype.
- `module_correlations_E3WT.txt`: Correlation values for E3WT genotype.
- `module_correlations_E4WT.txt`: Correlation values for E4WT genotype.
- `Microglia-02-networkConstruction-auto.RData`: Network construction results.
- `Sample clustering dendrogram`: Used to identify outliers.
- `Module-trait relationships heatmap`: Visualizes module associations with traits.

## Notes
- Adjust the cut height (`h = 80`) if necessary when detecting outliers.
- Modify the list of modules in `modules` to visualize specific ones.
- Ensure sample names in `samplesdf_filtered` match row names in `datExpr`.

## Contact
For any issues or inquiries, please reach out to the project contributors.

