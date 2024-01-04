# CRISPR Screen Processing Scripts

R scripts for processing guide-Seq counts. Input expected to be raw counts from
`mageck count`. Scripts should be run in order.

### Scripts

#### `PrepareScreenCounts.R`

This script should be run in an interactive RStudio environment.
Functions for normalizing counts, checking positive controls, plotting PCA and
plotting heatmaps.

Input: `mageck count` data table; guide library tsv
Output: normalized count (DESeq2), MAUDE compatible; bin statistics matrix, MAUDE compatible

##### Dependencies
```
library(DESeq2)
library(ggpubr)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(tidyr)
library(reshape)
```

#### `MAUDEAnalysis.R`

Perform single-sample MAUDE analysis using output from `PrepareScreenCounts.R`.
All ".tsv" files in "/Counts/" subdirectories are included. `binStats.txt` path
is hard coded currently.

Input: */Counts/<previous output>.tsv; binStats.txt
Output: MAUDE output at the element `_guideLevelStats.tsv` and gene level `_geneLevelStats.tsv`

##### Example

`./scripts/MAUDEAnalysis.R`

##### Dependencies
```
library(MAUDE)
```

#### `PlotRankedGenes.R`

Plot ranked gene significanceZ scores for all genes from MAUDE analysis. Genes with
sigZ > 3 or < -3 are in red. A comma separated list of genes can be added for additional
highlighting.

Input: directory to search; tsv suffix pattern; gene names for highlighting
Output: ranked gene .png

##### Example

`./scripts/PlotRankedGenes.R <top_level_directory> <tsv_suffix> <optional_genes_to_highlight>`
`./scripts/PlotRankedGenes.R ./MAUDE/ _geneLevelStats.tsv Tlr4,Myd88`


##### Dependencies
```
library(ggplot2)
library(dplyr)
library(ggrepel)
```

#### `CorrelationPlotting.R`

Generates scatterplots by plotting gene significanceZ score pairwise amongst all
input files. Only the top "n" (second argument) most significant are plotted and
used for correlation calculation. Some genes are hard coded to be included whether
or not they are in the top "n" genes.

Input: comma separated list of MAUDE output tsvs; number of genes to plot; output file name
Output: output png with pairwise scatterplots of all inputs and pearson correlations

##### Example

`./scripts/CorrelationPlotting.R <input1.tsv,input2.tsv,...,inputN.tsv>
<number_of_genes> <output>`
`./scripts/CorrelationPlotting.R ./MAUDE/<subdir>/xxxx_geneLevelStats.tsv,./MAUDE/<subdir>/yyyy_geneLevelStats.tsv
25 ./outdir/xxxx_yyyyy_correlations.png`

##### Dependencies
```
library(ggplot2)
library(ggrepel)
library(dplyr)
library(purrr)
library(tidyr)
```
