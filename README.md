# Selection Shadow across Metazoa

## Description 

This repository contains code used for the analyses in the manuscript entitled:

[**"Under the shadow: Old-biased genes are subject to weak purifying selection at both the tissue and cell type-specific levels"**](https://www.biorxiv.org/)

Any finding reported on the manuscript can be reproduced using the contents of the repository. 

## Directory Structure

### data

[*data*](data) directory includes all of the gene expression matrices used in the study (except for the gene expression matrices of the Tabula Muris Senis dataset). The main *data* directory contains subdirectories for each of the species used in the study. How expression matrices were generated are explained in the Methods section of the manuscript. 

For the Tabula Muris Senis dataset, information on how to download single cell expression values are given in the README file under

```
mus_musculus/tabula_muris_senis/
```

and processed expression values can be found under the directory

```
mus_musculus/tabula_muris_senis/R/results/
```
### analysis scripts

The main directory contains 5 different directories named after each of the species used in the study.

R scripts used for the generation of statistical analysis results can be found under these directories:

```
<organism>/<accession>/R/<accesion>_Analysis.R
```

### figure scipts

[*figure_scripts*](figure_scripts) directory contains all the scripts used to generate figures used in the study. Running these scripts outputs the figures under the [*results_graphs*](results_graphs) directory. 

### supplements

Any supplmentary information such as tables, generated using code, can be found under [*supplements*](supplements) directory. 

## Running order

To reproduce the concents of this reprository, analysis results first must be run individually, which generates all the necessary files for the generation of plots and tables.
