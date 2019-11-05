<div align="center">
<a name="logo"/>
<img src="./inst/TFregulomeR_logo.png" alt="TFregulomeR Logo" ></img>
</a>
</div>



TFregulomeR

Version: 1.2.0

# Introduction
*TFregulomeR* comprises of a comprehensive compendium of transcription factor binding sites (TFBSs) derived from the MethMotif and GTRD, as well as the ready-to-use functionality in R language facilitating data access, integration and analysis. The binding motifs predicted in-silico from MethMotif and GTRD describe cell specific transcription factor (TF) binding propensities, while the DNA methylation profiles from MethMotif portray a second epigenetic dimension in TF binding events. The whole toolbox allows a better understanding of the TF binding propensities in a cell-specific manner. 

-------

## Release notes
 This repository is TFregulomeR stable release 

 1) Current TFregulomeR stable version: 1.2.0 (Updated on 6 September 2019).


### For development release, please visit [TFregulomeR-dev](https://github.com/linquynus/TFregulomeR-dev) 

1) Current TFregulomeR developement version 1.2.1 (Updated on 25 October 2019).
2) For linking to human and mouse TF data compendium, please use TFregulomeR test version 1.99.1 @ brunch _with_mouse_version_ (Updated on 1 November 2019).

-------


## Documentation
You can check detailed package instructions in [Vignettes](https://bioinfo-csi.nus.edu.sg/methmotif/API_TFregulomeR/TFregulomeR-Vignettes.html)

-------

## Current Functionalities v1.2.0
### Click [here](./inst/update_notes/functionality_update.md) for functionality update notes

1) Browse the TFregulomeR data warehouse (dataBrowser)
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

-------

## Current TFBSs in TFregulomeR compendium

 ### Click [here](./inst/update_notes/compendium_update.md) for TFregulomeR compendium update notes

TFregulomeR data compendium version: 1.2.0

| Item     | Count |
| :---------:|:------:|
| PWM     | 1468   |
| ChIP-seq experiments    | 3891   |
| Unique TF     | 415   |
| PWM with DNA methylation records    | 655   |
| Species     | human (hg38)  |
| Organ   | stem_cell, blood_and_lymph, connective_tissue, colorectum, brain, bone, stomach, prostate, breast, pancreas, skin, kidney, lung, eye, esophagus, heart, muscle, uterus, spleen, cervix, testis, liver, adrenal_gland, neck_and_mouth, pleura, ovary, thymus, fallopian, vagina   |
| Sample type | primary_cells, cell_line, tissue
| Cell or tissue | 414 |
| Disease state | normal, tumor, Simpson_Golabi_Behmel_syndrome, progeria, metaplasia, unknown, immortalized, premetastatic|
| Source | GTRD, MethMotif | 

-------

## Citation

A manuscript describing *TFregulomeR* has been submitted for peer review. If you currently use *TFregulomeR*, please cite us as follows:


Quy Xiao Xuan Lin, Denis Thieffry, Sudhakar Jha, Touati Benoukraf.  _TFregulomeR reveals transcription factors’ context-specific features and functions._ 2019. https://github.com/benoukraflab/TFregulomeR

-------

## Case studies

The scripts of case studies used in our manuscript are available as below.

1. [Case study of CEBPB](./inst/case_study/case_study_of_CEBPB.R)
2. [Case study of MAFF](./inst/case_study/case_study_of_MAFF.R)
3. [Case study of ATF3](./inst/case_study/case_study_of_ATF3.R)


-------

## Installation

#### Prerequisite packages

1) Required packages: the packages below are the basic prerequisite packages for *TFregulomeR* functionalities  (Will be installed automatically)

    - [jsonlite](https://cran.r-project.org/web/packages/jsonlite/index.html) (>= 1.5)
    - [ggplot2](https://cran.r-project.org/web/packages/ggplot2/index.html) (>= 3.0.0)
    - [ggseqlogo](https://cran.r-project.org/web/packages/ggseqlogo/index.html) (>= 0.1)
    - [gridExtra](https://cran.r-project.org/web/packages/gridExtra/index.html) (>= 2.3)
    - [GenomicRanges](https://bioconductor.org/packages/release/bioc/html/GenomicRanges.html) (>= 1.32.7)
    - [curl](https://cran.r-project.org/web/packages/curl/index.html) (>= 3.2)

2) Optional packages: the packages below are optional since they are required only in some functions or some options in a function

    - [rGREAT](https://bioconductor.org/packages/release/bioc/html/rGREAT.html) (>= 1.14.0): only requried in `greatAnnotate()`
    - [liftOver](https://bioconductor.org/packages/release/workflows/html/liftOver.html) (>= 1.4.0): only required when users use hg38 peaks in `greatAnnotate()`. Since GREAT analysis doesn't support hg38, hg38 peaks will be converted to hg19 using liftOver.
    - [rbokeh](https://cran.r-project.org/web/packages/rbokeh/index.html) (>= 0.5.0): only required when users opt to export an intuitive HTML report in `greatAnnotate()`
    - [TxDb.Hsapiens.UCSC.hg38.knownGene](https://bioconductor.org/packages/release/data/annotation/html/TxDb.Hsapiens.UCSC.hg38.knownGene.html) (>= 3.4.0): only required when users opt to annotate hg38 peak locations in `genomeAnnotate()`
    - [TxDb.Hsapiens.UCSC.hg19.knownGene](https://bioconductor.org/packages/release/data/annotation/html/TxDb.Hsapiens.UCSC.hg19.knownGene.html) (>= 3.2.2): only required when users opt to annotate hg19 peak locations in `genomeAnnotate()`
    - [TFBSTools](http://bioconductor.org/packages/release/bioc/html/TFBSTools.html) (>= 1.20.0): only required in `toTFBSTools()`

#### Install

In R console,

```r
# if you have not installed "devtools" package
install.packages("devtools")
devtools::install_github("benoukraflab/TFregulomeR")
```
The step above will automatically install the required packages. However, you still need to install optional packages if you opt to use the functions such as `greatAnnotate()`, `genomeAnnotate()` and `toTFBSTools()`.

-------

## License

Artistic-2.0
