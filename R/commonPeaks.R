#' commonPeaks
#'
#' This function allows you to obtain the subsets of a list of peak sets common with another list of peak sets (either from MethMotif database or the self-provided).
#' @param target_peak_list Required. List of data.frames, each of which contains bed-format peak regions. They are the peak sets you want to get the common subsets from, and can be loaded from MethMotif database or self-provided.
#' @param target_peak_id Required. Character of vector, each of which is a unique ID corresponding to the element in "target_peak_list". If a peak set is from MethMotif Database, its MethMotif ID should be used here.
#' @param compared_peak_list Required. List of data.frames, each of which contains bed-format peak regions. They are the peak sets you want to compare with the "target_peak_list", and can be loaded from MethMotif database or self-provided.
#' @param compared_peak_id Required. Character of vector, each of which is a unique ID corresponding to the element in "compared_peak_list". If a peak set is from MethMotif Database, its MethMotif ID should be used here.
#' @param motif_format Required. Motif PFM format, either in MEME by default or TRANSFAC.
#' @param TFregulome_url Optional. If the MethMoitf url is NO more "http://bioinfo-csi.nus.edu.sg/methmotif/", please use a new url.
#' @return  matrix of CommonPeaksMM class objects
#' @keywords commonPeaks
#' @export
#' @examples
#' target_peaks <- list(loadPeaks(id = "MM1_HSA_K562_CEBPB"),
#'                      read.delim("my_own_peaks.bed", header = F))
#' target_id <- c("MM1_HSA_K562_CEBPB", "my_own_peaks")
#' compared_peaks <- list(loadPeaks(id = "MM1_HSA_HepG2_CEBPB"),
#'                        read.delim("peaks_to_common_with.bed", header = F))
#' compared_id <- c("MM1_HSA_HepG2_CEBPB", "peaks_to_common_with")
#' commonPeaks_output <- commonPeaks(target_peak_list=target_peaks,
#'                                   target_peak_id=target_id,
#'                                   compared_peak_list=compared_peaks,
#'                                   compared_peak_id=compared_id)

commonPeaks <- function(target_peak_id, motif_only_for_target_peak = F,user_target_peak_list,
                        compared_peak_id, motif_only_for_compared_peak = F, user_compared_peak_list,
                        motif_type = "MEME", TFregulome_url)
{
  # check the input arguments
  if(missing(target_peak_id) && missing(user_target_peak_list))
  {
    stop("No target peak input. Please EITHER input TFregulome peaks using TFregulome ID(s) by 'target_peak_id = ' OR your own peak list using a list of data.frame(s) containing bed-format regions by 'user_target_peak_list = '")
  }
  if(missing(compared_peak_id) && missing(user_compared_peak_list))
  {
    stop("No compared peak input. Please EITHER input TFregulome peaks using TFregulome ID(s) by 'compared_peak_id = ' OR your own peak list using a list of data.frame(s) containing bed-format regions by 'user_compared_peak_list = '")
  }
  if ((!missing(user_target_peak_list) && class(user_target_peak_list) != "list") ||
      (!missing(user_compared_peak_list) && class(user_compared_peak_list) != "list"))
  {
    stop("The class of input 'user_target_peak_list' and 'user_compared_peak_list' should be 'list', a list of bed-like data.frame storing peak regions!")
  }
  if (class(motif_only_for_target_peak) != "logical" || class(motif_only_for_compared_peak) != "logical")
  {
   stop("motif_only_for_target_peak and motif_only_for_compared_peak should be either TRUE or FALSE (default)")
  }
  if (motif_type != "MEME" && motif_type != "TRANSFAC")
  {
    stop("motif_type should be either 'MEME' (default) or 'TRANSFAC'!")
  }

  # make an appropriate API url
  if (missing(TFregulome_url)){
    TFregulome_url <- "http://localhost:8888/api/table_query/"
  } else if (endsWith(TFregulome_url, suffix = "/index.php")==TRUE){
    TFregulome_url <- gsub("index.php", "", TFregulome_url)
    TFregulome_url <- paste0(TFregulome_url, "api/table_query/")
  } else if (endsWith(TFregulome_url, suffix = "/")==TRUE){
    TFregulome_url <- paste0(TFregulome_url, "api/table_query/")
  } else {
    TFregulome_url <- paste0(TFregulome_url, "/api/table_query/")
  }

  message("TFregulomeR::commonPeaks() starting ... ...")

  # loading target peak list
  message("Loading target peak list ... ...")
  target_peak_list_all <- list()
  # loading from TFregulome server
  TFregulome_target_peak_id <- c()
  target_list_count <- 0
  if (!missing(target_peak_id) && length(target_peak_id)>0)
  {
    message(paste0("... You have ", length(target_peak_id)," TFBS(s) requested to be loaded from TFregulome server"))
    if (motif_only_for_target_peak == T)
    {
      message("... You chose to load TF peaks with motif only. Using 'motif_only_for_target_peak' tunes your options")
    }
    else
    {
      message("... You chose to load TF peaks regardless of presence of motif. Using 'motif_only_for_target_peak' tunes your options")
    }
    message("... loading TFBS(s) from TFregulome now")
    for (i in target_peak_id)
    {
      peak_i <- suppressMessages(loadPeaks(id = i, includeMotifOnly = motif_only_for_target_peak))
      if (is.null(peak_i))
      {
        message(paste0("... ... NO peak file for your id '", i,"'."))
      }
      else
      {
        target_list_count <- target_list_count + 1
        target_peak_list_all[[target_list_count]] <- peak_i
        TFregulome_target_peak_id <- c(TFregulome_target_peak_id, i)
        message(paste0(".. ... peak file loaded successfully for id '", i,"'"))
      }
    }
    message("... Done loading TFBS(s) from TFregulome")
  }
  # users' peaks
  user_target_peak_id <- c()
  if (!missing(user_target_peak_list) && length(user_target_peak_list)>0)
  {
    message(paste0("... You have ",length(user_target_peak_list)," customised peak set(s)"))
    if (is.null(names(user_target_peak_list)) ||
        length(unique(names(user_target_peak_list)))!=length(names(user_target_peak_list)))
    {
      message("... ... You didn't provide the name for each customised peak set or your names are not unique. Instead we will use 'user_target_peak1', 'user_target_peak2'..." )
      user_target_peak_id <- paste0("user_target_peak", seq(1,length(user_target_peak_list), 1))
    }
    else
    {
      user_target_peak_id <- names(user_target_peak_list)
    }
    for (i in 1:length(user_target_peak_list))
    {
      peak_i <- user_target_peak_list[[i]]
      peak_i_sub <- peak_i[,1:3]
      colnames(peak_i_sub) <- c("chr","start","end")
      peak_i_sub$id <- paste0(user_target_peak_id[i], "_", as.vector(rownames(peak_i_sub)))
      peak_i_sub <- peak_i_sub[,c("chr","start","end","id")]
      target_list_count <- target_list_count + 1
      target_peak_list_all[[target_list_count]] <- peak_i_sub
    }
  }
  # combine TFregulome ID and user ID
  target_peak_id_all <- c(TFregulome_target_peak_id, user_target_peak_id)

  # loading compared peak list
  message("Loading compared peak list ... ...")
  compared_peak_list_all <- list()
  # loading from TFregulome server
  compared_list_count <- 0
  TFregulome_compared_peak_id <- c()
  if (!missing(compared_peak_id) && length(compared_peak_id)>0)
  {
    message(paste0("... You have ", length(compared_peak_id)," TFBS(s) requested to be loaded from TFregulome server"))
    if (motif_only_for_compared_peak == T)
    {
      message("... You chose to load TF peaks with motif only. Using 'motif_only_for_compared_peak' tunes your options")
    }
    else
    {
      message("... You chose to load TF peaks regardless of presence of motif. Using 'motif_only_for_compared_peak' tunes your options")
    }
    message("... loading TFBS(s) from TFregulome now")
    for (i in compared_peak_id)
    {
      peak_i <- suppressMessages(loadPeaks(id = i, includeMotifOnly = motif_only_for_compared_peak))
      if (is.null(peak_i))
      {
        message(paste0("... ... NO peak file for your id '", i,"'."))
      }
      else
      {
        compared_list_count <- compared_list_count + 1
        compared_peak_list_all[[compared_list_count]] <- peak_i
        TFregulome_compared_peak_id <- c(TFregulome_compared_peak_id, i)
        message(paste0(".. ... peak file loaded successfully for id '", i,"'"))
      }
    }
    message("... Done loading TFBS(s) from TFregulome")
  }
  # users' peaks
  if (!missing(user_compared_peak_list) && length(user_compared_peak_list)>0)
  {
    message(paste0("... You have ",length(user_compared_peak_list)," customised peak set(s)"))
    for (i in 1:length(user_compared_peak_list))
    {
      peak_i <- user_compared_peak_list[[i]]
      peak_i_sub <- peak_i[,1:3]
      colnames(peak_i_sub) <- c("chr","start","end")
      peak_i_sub$id <- paste0("compared_peak_", as.vector(rownames(peak_i_sub)))
      peak_i_sub <- peak_i_sub[,c("chr","start","end","id")]
      compared_list_count <- compared_list_count + 1
      compared_peak_list_all[[compared_list_count]] <- peak_i_sub
    }
  }
  # start analysing
  common_peak_matrix <- list()
  for (i in 1:length(target_peak_list_all))
  {
    target_id_i <- target_peak_id_all[i]
    target_peak_i <- target_peak_list_all[[i]]
    number_of_orignal_target <- nrow(target_peak_i)
    message(paste0("Start analysing: ", target_id_i, "... ..."))

    ## if it is from TFregulome
    if (i <= length(TFregulome_target_peak_id))
    {
      isTFregulome_target <- TRUE
      query_url <- paste0("listTFBS.php?AllTable=F&id=",target_id_i)
      request_content_json <- tryCatch({
        fromJSON(paste0(TFregulome_url,query_url))
      },
      error = function(cond)
      {
        message("There is a warning to connect TFregulome API!")
        message("Advice:")
        message("1) Check internet access;")
        message("2) Check dependent package 'jsonlite';")
        message("3) Current TFregulome server is implemented in MethMotif database, whose homepage is 'http://bioinfo-csi.nus.edu.sg/methmotif/'. If MethMotif homepage url is no more valid, please Google 'MethMotif', and input the valid MethMotif homepage url using 'TFregulome_url = '.")
        message(paste0("warning: ",cond))
        return(NULL)
      })
      request_content_df <- as.data.frame(request_content_json$TFBS_records)
      source_i <- request_content_df[,"source"]
      if (source_i == "MethMotif")
      {
        isMethMotifID_target <- TRUE
        motif_seq_path_target <- request_content_df[1,c("TFBS")]
        meth_file_path_target <- request_content_df[1,c("DNA_methylation_profile")]
        WGBS_replicate_target <- request_content_df[1,c("WGBS_num")]
      }
      else
      {
        isMethMotifID_target <- FALSE
        motif_seq_path_target <- request_content_df[1,c("TFBS")]
      }
    }
    else
    {
      isTFregulome_target <- FALSE
    }
    # comparing with compared peak list
    for (j in 1:length(compared_peak_list_all))
    {
      compared_peak_j <- compared_peak_list_all[[j]]
      if (isTFregulome_target)
      {
        bed_target_i <- with(target_peak_i, GRanges(chr, IRanges(start-99, end+100), id=id))
      }
      else
      {
        bed_target_i <- with(target_peak_i, GRanges(chr, IRanges(start, end), id=id))
      }
      if (j <= length(TFregulome_compared_peak_id))
      {
        bed_compared_j <- with(compared_peak_j, GRanges(chr, IRanges(start-99, end+100), id=id))
      }
      else
      {
        bed_compared_j <- with(compared_peak_j, GRanges(chr, IRanges(start, end), id=id))
      }
      # get target peak intersecting with compared peak
      # subsetOverlaps may mis-think the two sets coming from different references, so suppressWarnings here
      suppressWarnings(bedTarget_with_bedcompared <- subsetByOverlaps(bed_target_i, bed_compared_j))
      peakTarget_with_peakcompared <- unique(as.data.frame(bedTarget_with_bedcompared))
      target_peak_with_peakcompared <- target_peak_i[which(target_peak_i$id %in% peakTarget_with_peakcompared$id), ]
      target_peak_i <- target_peak_with_peakcompared
    }

    MethMotif_target <- new("MethMotif")
    if (isTFregulome_target)
    {
      motif_seq_target <- read.delim(motif_seq_path_target, sep = "\t", header = F)
      if (nrow(target_peak_i)>0)
      {
        #compute motif matrix
        colnames(motif_seq_target) <- c("chr","start","end","strand","weight", "pvalue","qvalue","sequence")
        motif_len_target <- nchar(as.character(motif_seq_target[1,"sequence"]))
        motif_seq_target$id <- paste0(target_id_i,"_motif_sequence_", as.vector(rownames(motif_seq_target)))
        motif_seq_target_grange <- with(motif_seq_target[,c("chr","start","end","id")], GRanges(chr, IRanges(start+1, end), id=id))
        bed_target_done_common <- with(target_peak_i[,c("chr","start","end","id")], GRanges(chr, IRanges(start-99, end+100), id=id))
        suppressWarnings(motif_of_bed_target_done_common <- subsetByOverlaps(motif_seq_target_grange, bed_target_done_common))
        motif_of_peakTarget_done_common <- unique(as.data.frame(motif_of_bed_target_done_common))
        if (nrow(motif_of_peakTarget_done_common)>0)
        {
          motif_of_peakTarget_done_common_allInfo <- motif_seq_target[which(motif_seq_target$id %in% motif_of_peakTarget_done_common$id),]
          motif_matrix_of_peakTarget_done_common <- formMatrixFromSeq(input_sequence = as.vector(motif_of_peakTarget_done_common_allInfo$sequence),
                                                                      motif_format = motif_type)
          # compute beta score matrix
          if (isMethMotifID_target)
          {
            # methylation file can be empty
            meth_level_target <- tryCatch(read.delim(meth_file_path_target, sep = "\t", header = F),
                                          error=function(e) data.frame())
            # methylation file can be empty
            if (nrow(meth_level_target)==0)
            {
              beta_score_matrix_of_peakTarget_done_common <- formBetaScoreFromSeq(input_meth = data.frame(),
                                                                                  WGBS_replicate = WGBS_replicate_target,
                                                                                  motif_len = motif_len_target)
            }
            else
            {
              colnames(meth_level_target) <- c("chr","start","end","meth_score","C_num","T_num","seq_chr","seq_start",
                                               "seq_end","strand","weight","pvalue","qvalue","sequence")
              meth_level_target$id <- paste0(target_id_i,"_motif_with_CG_", as.vector(rownames(meth_level_target)))
              meth_level_target_grange <- with(meth_level_target[,c("seq_chr","seq_start","seq_end","id")],
                                               GRanges(seq_chr, IRanges(seq_start, seq_end), id=id))
              suppressWarnings(meth_level_peakTarget_done_common <- unique(as.data.frame(subsetByOverlaps(meth_level_target_grange,
                                                                         motif_of_bed_target_done_common))))
              meth_level_peakTarget_done_common_allInfo <- meth_level_target[which(meth_level_target$id %in% meth_level_peakTarget_done_common$id),]
              beta_score_matrix_of_peakTarget_done_common <- formBetaScoreFromSeq(input_meth = meth_level_peakTarget_done_common_allInfo,
                                                                                  WGBS_replicate = WGBS_replicate_target,
                                                                                  motif_len = motif_len_target)
            }
          }
          else
          {
            beta_score_matrix_of_peakTarget_done_common <- as.matrix(NA)
          }

          if (motif_type == "TRANSFAC")
          {
            version_target <- 0
          }
          else
          {
            version_target <- 4
          }
          MethMotif_target@MMmotif <- updateMMmotif(MethMotif_target@MMmotif,
                                                    motif_format = motif_type,
                                                    version = version_target,
                                                    background = c("A"=0.25,"C"=0.25,"G"=0.25,"T"=0.25),
                                                    id = paste0(target_id_i,"_common_peaks"),
                                                    alternate_name = target_id_i,
                                                    width = motif_len_target,
                                                    nsites=nrow(motif_of_peakTarget_done_common),
                                                    motif_matrix=motif_matrix_of_peakTarget_done_common)
          MethMotif_target@MMBetaScore <- beta_score_matrix_of_peakTarget_done_common
        }
      }
    }
    new_CommonPeaksMM <- new("CommonPeaksMM")
    new_CommonPeaksMM <- updateCommonPeaksMM(theObject = new_CommonPeaksMM,
                                             id = paste0(target_id_i, "_common_peaks"),
                                             common_percentage = 100*nrow(target_peak_i)/number_of_orignal_target,
                                             common_peak = target_peak_i,
                                             isTFregulomeID = isTFregulome_target,
                                             MethMotif = MethMotif_target)
    common_peak_matrix[[paste0(target_id_i, "_common_peaks")]] <- new_CommonPeaksMM
  }
  message("Done analysing.")
  dim(common_peak_matrix) <- c(length(target_peak_id_all), 1)
  rownames(common_peak_matrix) <- c(target_peak_id_all)
  colnames(common_peak_matrix) <- c("common")
  return(common_peak_matrix)
}



