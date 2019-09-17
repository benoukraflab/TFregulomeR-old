## Current Functionalities v1.2.0 (Development version)
### What's new compared to last stable version v1.1.0 
1) Changed the function name _TFBSBrowser_ to _dataBrowser_.
2) More detailed information about motifs, such as consistency with exisiting database and correlation between different motif callers, added in the output of _dataBrowser_.
3) Fixed some bugs in function _cofactorReport_.
4) Added y axis and number of peaks with motif in the output of _plotLogo_.



## Previous Stable Version
### Functionalities version 1.1.0 - 31 July 2019

1) Browse the TFregulomeR data warehouse (TFBSBrowser)
2) Load TF peaks (loadPeaks)
3) Search motif matrix and DNA methylation score matrix (searchMotif)
4) Plot motif or MethMotif logo (plotLogo)
5) Export motif matrix and DNA methylation score matrix (exportMMPFM)
6) Get context-independent peaks along with DNA methylation profiles (commonPeaks & commonPeakResult)
7) Get context-dependent peaks along with DNA methylation profiles (exclusivePeaks & exclusivePeakResult)
8) Form a intersected matrix between two lists of peak sets along with DNA methylation profiles and read enrichments, for interactome and co-binding partner studies (intersectPeakMatrix & intersectPeakMatrixResult)
9) Automatically generate a PDF report for TF co-factors along with motif sequences, DNA methylation (within motif and in 200bp regions) and read enrichments (cofactorReport).
10) Plot the TFBS distribution in a given list of peak sets (motifDistrib & plotDistrib)
11) Annotate peak genomic locations (genomeAnnotate)
12) Annotate ontologies of target genes by a peak set (greatAnnotate)
13) Convert a motif matrix to a PFMatrix calss object for *TFBSTools* package (toTFBSTools)

### Functionalities version 1.0.0 - 7 May 2019
1) Browse the TFregulomeR data warehouse (TFBSBrowser)
2) Load TF peaks (loadPeaks)
3) Search motif matrix and DNA methylation score matrix (searchMotif)
4) Plot motif or MethMotif logo (plotLogo)
5) Export motif matrix and DNA methylation score matrix (exportMMPFM)
6) Get context-independent peaks along with DNA methylation profiles (commonPeaks & commonPeakResult)
7) Get context-dependent peaks along with DNA methylation profiles (exclusivePeaks & exclusivePeakResult)
8) Form a intersected matrix between two lists of peak sets along with DNA methylation profiles, for interactome and co-binding partner studies (intersectPeakMatrix & intersectPeakMatrixResult)
9) Plot the TFBS distribution in a given list of peak sets (motifDistrib & plotDistrib)
10) Annotate peak genomic locations (genomeAnnotate)
11) Annotate ontologies of target genes by a peak set (greatAnnotate)
12) Convert a motif matrix to a PFMatrix calss object for *TFBSTools* package (toTFBSTools)
