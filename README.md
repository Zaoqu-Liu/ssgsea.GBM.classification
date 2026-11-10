# ssgsea.GBM.classification

R package that classifies glioblastoma (GBM) expression profiles into the Proneural, Classical and Mesenchymal transcriptional subtypes with single-sample GSEA (ssGSEA). Each sample is scored against a 50-gene up-regulated signature per subtype, and the scores are converted to empirical p values with a resampling null, so a sample can be assigned to one subtype or to several.

## Installation

```R
devtools::install_github('Zaoqu-Liu/ssgsea.GBM.classification')
```

The scoring code attaches `gtools`, `verification`, `ROCR`, `MASS` and `MCMCpack`; install them first if they are missing.

## Usage

The input is a tab-delimited GCT file (version 1.2): two header lines (`#1.2`, then the number of rows and columns), a column header line, then one row per gene with the gene symbol in the first column (`NAME`), a description in the second column (`Description`; the gene symbol can be repeated), and one column per sample. Use uncentered, log-transformed expression such as microarray intensities or RPKM; ssGSEA ranks genes within each sample.

```R
library(ssgsea.GBM.classification)

# expr.gct: GCT file in the current directory
# 1000 permutations build the null distribution for the p values
runSsGSEAwithPermutation('expr.gct', 1000)
```

The function writes `p_result_<input file name>.txt` next to the input file (here `p_result_expr.gct.txt`) and removes its temporary files. The table has one row per sample and six columns:

- `Proneural`, `Classical`, `Mesenchymal`: raw ssGSEA scores
- `Proneural_pval`, `Classical_pval`, `Mesenchymal_pval`: empirical p values, computed as (number of random profiles scoring at least as high + 1) / (permutations + 1)

A subtype is called activated when its p value is small (the original study used p < 0.05); samples with significant p values for more than one subtype activate several subtypes. The smallest attainable p value is 1 / (permutations + 1), so use enough permutations for the threshold you need.

## Contents

- `R/runSsGSEAwithPermutationR3.R`: `runSsGSEAwithPermutation()`, the entry point
- `R/msig.library.12.R`: the MSIG/OPAM function library it relies on (`OPAM.apply.model.2` computes the ssGSEA scores)
- `data/Proneural.mod`, `data/Classical.mod`, `data/Mesenchymal.mod`: model files, each holding the 50 up-regulated signature genes of one subtype

## Reference

The approach, a ssGSEA-based classification with an equivalent-distribution resampling test on 50-gene subtype signatures, is described in: Wang Q, Hu B, Hu X, Kim H, et al. Tumor evolution of glioma-intrinsic gene expression subtypes associates with immunological changes in the microenvironment. *Cancer Cell* 2017;32(1):42-56.e6. doi: [10.1016/j.ccell.2017.06.003](https://doi.org/10.1016/j.ccell.2017.06.003). The `DESCRIPTION` file lists Qianghu Wang, first author of that paper, as the package author.
