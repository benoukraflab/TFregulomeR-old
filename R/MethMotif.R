MMmotif <- setClass(
  # Set the name for the class
  "MMmotif",
  # define the slots
  slots = c(
    motif_format = "character",
    version = "numeric",
    alphabet = "character",
    strand = "character",
    background = "vector",
    id = "character",
    alternate_name = "character",
    width = "integer",
    nsites = "integer",
    nPeaks = "integer",
    evalue = "numeric",
    motif_matrix = "matrix"
  ),
  # Set the default values for the slots.
  prototype = list(
    motif_format = "MEME",
    version = 4,
    alphabet = "ACGT",
    strand = "+ -",
    background = c("A"=0, "C"=0, "G"=0, "T"=0),
    id = "",
    alternate_name = "",
    width = 0L,
    nsites = 0L,
    nPeaks = 0L,
    evalue = 0,
    motif_matrix = matrix()
  )
)
setGeneric(name = "updateMMmotif",
           def = function(theObject, motif_format, version, alphabet, strand,
                          background, id, alternate_name, width, nsites, nPeaks,
                          evalue, motif_matrix)
           {
             standardGeneric("updateMMmotif")
           }
)
setMethod(f = "updateMMmotif",
          signature(theObject = "MMmotif"),
          definition = function(theObject, motif_format, version, alphabet, strand,
                                background, id, alternate_name, width, nsites, nPeaks,
                                evalue, motif_matrix)
          {
            if (missing(theObject))
            {
              stop("No MMmotif object.Please use 'theObject = '")
            }
            if (!missing(motif_format))
            {
              if (!is.character(motif_format))
              {
                stop("'motif_format' should be character!")
              }
              theObject@motif_format <- motif_format
            }
            if (!missing(version))
            {
              if (!is.numeric(version))
              {
                stop("'version' should be numeric!")
              }
              theObject@version <- version
            }
            if (!missing(alphabet))
            {
              if (!is.character(alphabet))
              {
                stop("'alphabet' should be character!")
              }
              theObject@alphabet <- alphabet
            }
            if (!missing(strand))
            {
              if (!is.character(strand))
              {
                stop("'strand' should be character!")
              }
              theObject@strand <- strand
            }
            if (!missing(background))
            {
              if (!is.vector(background))
              {
                stop("'background' should be named vector!")
              }
              theObject@background <- background
            }
            if (!missing(id))
            {
              if (!is.character(id))
              {
                stop("'id' should be character!")
              }
              theObject@id <- id
            }
            if (!missing(alternate_name))
            {
              if (!is.character(alternate_name))
              {
                stop("'alternate_name' should be character!")
              }
              theObject@alternate_name <- alternate_name
            }
            if (!missing(width))
            {
              if (!is.integer(width))
              {
                stop("'width' should be integer!")
              }
              theObject@width <- width
            }
            if (!missing(nsites))
            {
              if (!is.integer(nsites))
              {
                stop("'nsites' should be integer!")
              }
              theObject@nsites <- nsites
            }
            if (!missing(nPeaks))
            {
              if (!is.integer(nPeaks))
              {
                stop("'nPeaks' should be integer!")
              }
              theObject@nPeaks <- nPeaks
            }
            if (!missing(evalue))
            {
              if (!is.numeric(evalue))
              {
                stop("'evalue' should be numeric!")
              }
              theObject@evalue <- evalue
            }
            if (!missing(motif_matrix))
            {
              if (!is.matrix(motif_matrix))
              {
                stop("'motif_matrix' should be matrix!")
              }
              theObject@motif_matrix <- motif_matrix
            }
            return(theObject)
          }
)

# MethMotif object
MethMotif <- setClass(
  "MethMotif",
  slots = c(MMBetaScore = "matrix", MMmotif = "MMmotif"),
  prototype = list(MMBetaScore = matrix(), MMmotif = MMmotif())
)
setGeneric(name = "updateMethMotif",
           def = function(theObject, MMBetaScore, MMmotif)
           {
             standardGeneric("updateMethMotif")
           }
)
setMethod(f = "updateMethMotif",
          signature(theObject = "MethMotif"),
          definition = function(theObject, MMBetaScore, MMmotif)
          {
            if (missing(theObject))
            {
              stop("No MethMotif object.Please use 'theObject = '")
            }
            if (!missing(MMBetaScore))
            {
              if (!is.matrix(MMBetaScore))
              {
                stop("'MMBetaScore' should be matrix!")
              }
              theObject@MMBetaScore <- MMBetaScore
            }
            if (!missing(MMmotif))
            {
              if (class(MMmotif)[1]!="MMmotif")
              {
                stop("'MMmotif' should be MMmotif class!")
              }
              theObject@MMmotif <- MMmotif
            }
            return(theObject)
          }
)

# IntersectPeakMatrix object

IntersectPeakMatrix <- setClass(
  "IntersectPeakMatrix",
  slots = c(id = "character",
            id_x = "character",
            overlap_percentage_x = "numeric",
            isxTFregulomeID = "logical",
            MethMotif_x = "MethMotif",
            methylation_profile_x = "matrix",
            external_signal_x = "vector",
            tag_density_x = "vector",
            id_y = "character",
            overlap_percentage_y = "numeric",
            isyTFregulomeID = "logical",
            MethMotif_y = "MethMotif",
            methylation_profile_y = "matrix",
            external_signal_y = "vector",
            tag_density_y = "vector"),
  prototype = list(id = "",
                   id_x = "",
                   overlap_percentage_x = 0,
                   isxTFregulomeID = FALSE,
                   MethMotif_x = new('MethMotif'),
                   methylation_profile_x = matrix(),
                   external_signal_x = c(NA),
                   tag_density_x = c("peak_number"=NA, "mean"=NA,"SD"=NA,"median"=NA,
                                     "quartile_25"=NA, "quartile_75"=NA),
                   id_y = "",
                   overlap_percentage_y = 0,
                   isyTFregulomeID = FALSE,
                   MethMotif_y = new('MethMotif'),
                   methylation_profile_y = matrix(),
                   external_signal_y = c(NA),
                   tag_density_y = c("peak_number"=NA, "mean"=NA,"SD"=NA,"median"=NA,
                                     "quartile_25"=NA, "quartile_75"=NA))
)

setGeneric(name = "updateIntersectPeakMatrix",
           def = function(theObject, id, id_x, overlap_percentage_x, isxTFregulomeID,
                          MethMotif_x, methylation_profile_x, external_signal_x,
                          tag_density_x, id_y, overlap_percentage_y,isyTFregulomeID,
                          MethMotif_y, methylation_profile_y, external_signal_y,
                          tag_density_y)
           {
             standardGeneric("updateIntersectPeakMatrix")
           }
)
setMethod(f = "updateIntersectPeakMatrix",
          signature(theObject = "IntersectPeakMatrix"),
          definition = function(theObject, id, id_x, overlap_percentage_x,
                                isxTFregulomeID, MethMotif_x, methylation_profile_x,
                                external_signal_x, tag_density_x, id_y, overlap_percentage_y,
                                isyTFregulomeID, MethMotif_y, methylation_profile_y,
                                external_signal_y, tag_density_y)
          {
            if (missing(theObject))
            {
              stop("No IntersectPeakMatrix object.Please use 'theObject = '")
            }
            if (!missing(id))
            {
              if (!is.character(id))
              {
                stop("'id' should be character!")
              }
              theObject@id <- id
            }
            if (!missing(id_x))
            {
              if (!is.character(id_x))
              {
                stop("'id_x' should be character!")
              }
              theObject@id_x <- id_x
            }
            if (!missing(overlap_percentage_x))
            {
              if (!is.numeric(overlap_percentage_x))
              {
                stop("'overlap_percentage_x' should be numeric!")
              }
              theObject@overlap_percentage_x <- overlap_percentage_x
            }
            if (!missing(isxTFregulomeID))
            {
              if (!is.logical(isxTFregulomeID))
              {
                stop("'isxTFregulomeID' should be logical!")
              }
              theObject@isxTFregulomeID <- isxTFregulomeID
            }
            if (!missing(MethMotif_x))
            {
              if (class(MethMotif_x)[1]!="MethMotif")
              {
                stop("'MethMotif_x' should be class MethMotif!")
              }
              theObject@MethMotif_x <- MethMotif_x
            }
            if (!missing(methylation_profile_x))
            {
              if (class(methylation_profile_x)[1]!="matrix")
              {
                stop("'methylation_profile_x' should be class matrix!")
              }
              theObject@methylation_profile_x <- methylation_profile_x
            }
            if (!missing(external_signal_x))
            {
              if (!is.vector(external_signal_x))
              {
                stop("'external_signal_x' should be class vector!")
              }
              theObject@external_signal_x <- external_signal_x
            }
            if (!missing(tag_density_x))
            {
              if (!is.vector(tag_density_x))
              {
                stop("'tag_density_x' should be class vector!")
              }
              theObject@tag_density_x <- tag_density_x
            }
            if (!missing(id_y))
            {
              if (!is.character(id_y))
              {
                stop("'id_y' should be character!")
              }
              theObject@id_y <- id_y
            }
            if (!missing(overlap_percentage_y))
            {
              if (!is.numeric(overlap_percentage_y))
              {
                stop("'overlap_percentage_y' should be numeric!")
              }
              theObject@overlap_percentage_y <- overlap_percentage_y
            }
            if (!missing(isyTFregulomeID))
            {
              if (!is.logical(isyTFregulomeID))
              {
                stop("'isyTFregulomeID' should be logical!")
              }
              theObject@isyTFregulomeID <- isyTFregulomeID
            }
            if (!missing(MethMotif_y))
            {
              if (class(MethMotif_y)[1]!="MethMotif")
              {
                stop("'MethMotif_y' should be class MethMotif!")
              }
              theObject@MethMotif_y <- MethMotif_y
            }
            if (!missing(methylation_profile_y))
            {
              if (class(methylation_profile_y)[1]!="matrix")
              {
                stop("'methylation_profile_y' should be class matrix!")
              }
              theObject@methylation_profile_y <- methylation_profile_y
            }
            if (!missing(external_signal_y))
            {
              if (!is.vector(external_signal_y))
              {
                stop("'external_signal_y' should be class vector!")
              }
              theObject@external_signal_y <- external_signal_y
            }
            if (!missing(tag_density_y))
            {
              if (!is.vector(tag_density_y))
              {
                stop("'tag_density_y' should be class vector!")
              }
              theObject@tag_density_y <- tag_density_y
            }
            return(theObject)
          }
)

# ExclusivePeaksMM object

ExclusivePeaksMM <- setClass(
  "ExclusivePeaksMM",
  slots = c(id = "character",
            exclusive_percentage = "numeric",
            exclusive_peak = "data.frame",
            isTFregulomeID = "logical",
            MethMotif = "MethMotif",
            methylation_profile = "matrix"),
  prototype = list(id = "",
                   exclusive_percentage = 0,
                   exclusive_peak = data.frame(),
                   isTFregulomeID = FALSE,
                   MethMotif = new('MethMotif'),
                   methylation_profile = matrix())
)
setGeneric(name = "updateExclusivePeaksMM",
           def = function(theObject, id, exclusive_percentage, exclusive_peak,
                          isTFregulomeID, MethMotif, methylation_profile)
           {
             standardGeneric("updateExclusivePeaksMM")
           }
)
setMethod(f = "updateExclusivePeaksMM",
          signature(theObject = "ExclusivePeaksMM"),
          definition = function(theObject, id, exclusive_percentage, exclusive_peak,
                                isTFregulomeID, MethMotif, methylation_profile)
          {
            if (missing(theObject))
            {
              stop("No ExclusivePeaksMM object.Please use 'theObject = '")
            }
            if (!missing(id))
            {
              if (!is.character(id))
              {
                stop("'id' should be character!")
              }
              theObject@id <- id
            }
            if (!missing(exclusive_percentage))
            {
              if (!is.numeric(exclusive_percentage))
              {
                stop("'exclusive_percentage' should be numeric!")
              }
              theObject@exclusive_percentage <- exclusive_percentage
            }
            if (!missing(exclusive_peak))
            {
              if (!is.data.frame(exclusive_peak))
              {
                stop("'exclusive_peak' should be data.frame!")
              }
              theObject@exclusive_peak <- exclusive_peak
            }
            if (!missing(isTFregulomeID))
            {
              if (!is.logical(isTFregulomeID))
              {
                stop("'isTFregulomeID' should be logical!")
              }
              theObject@isTFregulomeID <- isTFregulomeID
            }
            if (!missing(MethMotif))
            {
              if (class(MethMotif)[1]!="MethMotif")
              {
                stop("'MethMotif' should be class MethMotif!")
              }
              theObject@MethMotif <- MethMotif
            }
            if (!missing(methylation_profile))
            {
              if (class(methylation_profile)[1]!="matrix")
              {
                stop("'methylation_profile' should be class matrix!")
              }
              theObject@methylation_profile <- methylation_profile
            }
            return(theObject)
          }
)

# CommonPeaksMM object

CommonPeaksMM <- setClass(
  "CommonPeaksMM",
  slots = c(id = "character",
            common_percentage = "numeric",
            common_peak = "data.frame",
            isTFregulomeID = "logical",
            MethMotif = "MethMotif",
            methylation_profile = "matrix"),
  prototype = list(id = "",
                   common_percentage = 0,
                   common_peak = data.frame(),
                   isTFregulomeID = FALSE,
                   MethMotif = new('MethMotif'),
                   methylation_profile = matrix())
)
setGeneric(name = "updateCommonPeaksMM",
           def = function(theObject, id, common_percentage, common_peak,isTFregulomeID,
                          MethMotif, methylation_profile)
           {
             standardGeneric("updateCommonPeaksMM")
           }
)
setMethod(f = "updateCommonPeaksMM",
          signature(theObject = "CommonPeaksMM"),
          definition = function(theObject, id, common_percentage, common_peak,isTFregulomeID,
                                MethMotif, methylation_profile)
          {
            if (missing(theObject))
            {
              stop("No CommonPeaksMM object.Please use 'theObject = '")
            }
            if (!missing(id))
            {
              if (!is.character(id))
              {
                stop("'id' should be character!")
              }
              theObject@id <- id
            }
            if (!missing(common_percentage))
            {
              if (!is.numeric(common_percentage))
              {
                stop("'common_percentage' should be numeric!")
              }
              theObject@common_percentage <- common_percentage
            }
            if (!missing(common_peak))
            {
              if (!is.data.frame(common_peak))
              {
                stop("'common_peak' should be data.frame!")
              }
              theObject@common_peak <- common_peak
            }
            if (!missing(isTFregulomeID))
            {
              if (!is.logical(isTFregulomeID))
              {
                stop("'isTFregulomeID' should be logical!")
              }
              theObject@isTFregulomeID <- isTFregulomeID
            }
            if (!missing(MethMotif))
            {
              if (class(MethMotif)[1]!="MethMotif")
              {
                stop("'MethMotif' should be class MethMotif!")
              }
              theObject@MethMotif <- MethMotif
            }
            if (!missing(methylation_profile))
            {
              if (class(methylation_profile)[1]!="matrix")
              {
                stop("'methylation_profile' should be class matrix!")
              }
              theObject@methylation_profile <- methylation_profile
            }
            return(theObject)
          }
)
