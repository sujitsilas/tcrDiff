# tcRdiff
TcRdiff is a R package for differential motif analysis. It uses results from GLIPH2 (Grouping of Lymphocyte Interactions by Paratope Hotspots 2), a computational tool that aids in grouping TCR sequences based on shared patterns in the complementarity-determining region 3 (CDR3), facilitating understanding of TCR diversity and antigen recognition patterns in immune responses. TcRdiff helps with an in-depth characterization of T-cell repertoires. It uses the output from GLIPH2 to distinguishing clonally expanded motifs. Clonally expanded and enriched motifs are identified by comparing the summed contribution scores of samples from each cluster (summed template frequency) using exact poisson statistics. Enriched motifs' log2-fold change and p-values can be visualized as volcano plots to identify expanded TCRs.

# Steps
### Step 1: Add a "sample_name" column to the GLIPH output file wtih only sample names. 
### Step 2: Use functions from the tcRdiff package to perform differential motif analysis and make volcano plots.


# Installing package
```R
require(devtools)
install_github("sujitsilas/tcRdiff")
```
tcRdiff can be successfully installed on Mac OS X, Linux, and Windows, using the devtools package directly from GitHub

