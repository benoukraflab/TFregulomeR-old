## Current Functionalities v2.0.0 - 5 March 2020
### What's new compared to last version v1.2.2 
1) Added new TFregulomeR compendium server in Canada. Users now can choose servers between Singapore and Canada.
2) linked to mouse TFregulomeR compendium in Canada server, coming soon in Singapore server.


## Previous  Version
### Functionalities version 1.2.2 - 10 January 2020
1) Modified the functions _intersectPeakMatrix_ and _intersectPeakMatrixResult_ to allow profiling of users' input external signals during pair-wise comparison analysis.
2) Added a new function _interactome3D_ which could directly generate a dynamic 3D interface report showing TF interactome coupled with CpG methylation and users' input external signal values.

### Functionalities version 1.2.1 - 25 October 2019
1) Improved the function _greatAnnotate_, which now allows to use peak regions in hg38 directly with rGREAT >= 1.16.1 (no need liftover package).

### Functionalities version 1.2.0 - 6 September 2019
1) Changed the function name _TFBSBrowser_ to _dataBrowser_.
2) More detailed information about motifs, such as consistency with exisiting database and correlation between different motif callers, added in the output of _dataBrowser_.
3) Fixed some bugs in function _cofactorReport_.
4) Added y axis and number of peaks with motif in the output of _plotLogo_.

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
