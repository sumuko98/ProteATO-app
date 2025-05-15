# ProteATO app
 R Shiny App for the Analysis of Proteomics Data

This repository contains the code for a Shiny application designed for reproducible analysis of proteomics data. The app is built using R and leverages the DEP package to streamline common steps in a proteomics workflow, from data preprocessing to statistical analysis of the differentially expressed data.

## Installation

The app could be installed locally via [GitHub](https://github.com/sumuko98/ProteATO-app) using the following commands (please note that `devtools` package should be installed before installing the app):


1. Check if devtools is installed, and install if not

`if (!requireNamespace("devtools", quietly = TRUE)) { install.packages("devtools") }`

2. Use `devtools` to install the ShinyApp

`devtools::install_github("sumuko98/ProteATO-app")`


3. Finally, run the app

`library(ProteATO)`

`ProteATO::runApp()`

## Features

**Data Input**: Upload of the proteomics data and experimental design files.

**Preprocessing**: Removes contaminants, handles duplicate gene names, and provides options for missing value filtering.

**Normalization**: Offers Variance Stabilizing Normalization (VSN) to correct systematic biases.

**Imputation**: Supports multiple methods for handling missing values, including "bpca", "knn", "MinProb", and "QRILC".

**Exploratory Analysis**: Includes a Venn diagram, overlap plots, protein counts visualization, missing value heatmaps, and intensity distribution plots.

**Differential Analysis**: Conducts differential enrichment analysis with customizable contrasts and visualizations like volcano plots and p-value histograms.

**Functional Enrichment**: Performs Gene Ontology (GO) and KEGG pathway enrichment analyses.

**PCA and ANOVA**: Provides PCA plots and supports ANOVA for selected conditions.
