#' exclusivePeaks
#'
#' This function allows you to obtain a list of exclusive peak subsets along with the DNA methylation profiles.
#' @param target_peak_id Character of vector, each of which is a TFregulome ID. Each of target peak will be compared with all "excluded peaks" to get its exclusive subset.
#' @param motif_only_for_target_peak Either TRUE of FALSE (default). If TRUE, only peaks with motif will be loaded for each TFregulome ID in target_peak_id.
#' @param user_target_peak_list A list of data.frames, each of which contains user's own bed-format target peak regions.
#' @param user_target_peak_id Character of vector, each of which is a unique ID corresponding to each peak set in the list user_target_peak_list. If the IDs are not provided or not unique, the function will automatically generate the IDs of its own. If any of the peak sets is derived from TFregulome database, its TFregulome ID should be used here correspondingly.
#' @param excluded_peak_id Character of vector, each of which is a TFregulome ID.
#' @param motif_only_for_excluded_peak Either TRUE of FALSE (default). If TRUE, only peaks with motif will be loaded for each TFregulome ID in excluded_peak_id.
#' @param user_excluded_peak_list A list of data.frames, each of which contains user's own bed-format excluded peak regions.
#' @param user_excluded_peak_id Character of vector, each of which is a unique ID corresponding to each peak set in the list user_excluded_peak_list. If the IDs are not provided or not unique, the function will automatically generate the IDs of its own. If any of the peak sets is derived from TFregulome database, its TFregulome ID should be used here correspondingly.
#' @param methylation_profile_in_narrow_region Either TRUE (default) of FALSE. If TRUE, methylation states in 200bp window surrounding peak summits for each exclusive peak from target_peak_id and user_target_peak_list (with TFregulome ID).
#' @param motif_type Motif PFM format, either in MEME by default or TRANSFAC.
#' @param TFregulome_url TFregulome server is implemented in MethMotif server. If the MethMoitf url is NO more "http://bioinfo-csi.nus.edu.sg/methmotif/", please use a new url.
#' @return  matrix of ExclusivePeaksMM class objects
#' @keywords exclusivePeaks
#' @export
#' @examples
#' target_id <- "MM1_HSA_K562_CEBPB"
#' excluded_id <- c("MM1_HSA_HepG2_CEBPB", "MM1_HSA_HCT116_CEBPB")
#' excluPeak_output <- exclusivePeaks(target_peak_id=target_id,
#'                                    motif_only_for_target_peak=TRUE,
#'                                    excluded_peak_id=excluded_id,
#'                                    motif_only_for_excluded_peak=TRUE,
#'                                    methylation_profile_in_narrow_region=TRUE)

exclusivePeaks <- function(target_peak_id, motif_only_for_target_peak = F,
                           user_target_peak_list, user_target_peak_id,
                           excluded_peak_id, motif_only_for_excluded_peak = F,
                           user_excluded_peak_list, user_excluded_peak_id,
                           methylation_profile_in_narrow_region = T, motif_type = "MEME",
                           TFregulome_url)
{
  # check the input arguments
  if(missing(target_peak_id) && missing(user_target_peak_list))
  {
    stop("No target peak input. Please input TFregulome peaks using TFregulome ID(s) by 'target_peak_id = ' OR your own peak list using a list of data.frame(s) containing bed-format regions by 'user_target_peak_list = '")
  }
  if(missing(excluded_peak_id) && missing(user_excluded_peak_list))
  {
    stop("No excluded peak input. Please input TFregulome peaks using TFregulome ID(s) by 'excluded_peak_id = ' OR your own peak list using a list of data.frame(s) containing bed-format regions by 'user_excluded_peak_list = '")
  }
  if ((!missing(user_target_peak_list) && class(user_target_peak_list) != "list") ||
      (!missing(user_excluded_peak_list) && class(user_excluded_peak_list) != "list"))
  {
    stop("The class of input 'user_target_peak_list' and 'user_excluded_peak_list' should be 'list', a list of bed-like data.frame storing peak regions!")
  }
  if (class(motif_only_for_target_peak) != "logical" || class(motif_only_for_excluded_peak) != "logical")
  {
    stop("motif_only_for_target_peak and motif_only_for_excluded_peak should be either TRUE or FALSE (default)")
  }
  if (class(methylation_profile_in_narrow_region) != "logical")
  {
    stop("methylation_profile_in_narrow_region should be either TRUE or FALSE (default)")
  }
  if (motif_type != "MEME" && motif_type != "TRANSFAC")
  {
    stop("motif_type should be either 'MEME' (default) or 'TRANSFAC'!")
  }

  # make an appropriate API url
  if (missing(TFregulome_url)){
    TFregulome_url <- "http://bioinfo-csi.nus.edu.sg/methmotif/api/table_query/"
  } else if (endsWith(TFregulome_url, suffix = "/index.php")==TRUE){
    TFregulome_url <- gsub("index.php", "", TFregulome_url)
    TFregulome_url <- paste0(TFregulome_url, "api/table_query/")
  } else if (endsWith(TFregulome_url, suffix = "/")==TRUE){
    TFregulome_url <- paste0(TFregulome_url, "api/table_query/")
  } else {
    TFregulome_url <- paste0(TFregulome_url, "/api/table_query/")
  }

  message("TFregulomeR::exclusivePeaks() starting ... ...")
  if (methylation_profile_in_narrow_region)
  {
    message("You chose to profile the methylation levels in 200bp window around peak summits, if there is any peak loaded from TFregulome")
  }
  else
  {
    message("You chose NOT to profile the methylation levels in 200bp window around peak summits")
  }
  # loading target peak list
  message("Loading target peak list ... ...")
  target_peak_list_all <- list()
  # loading from TFregulome server
  TFregulome_target_peak_id <- c()
  is_taregt_TFregulome <- c()
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
      peak_i <- suppressMessages(loadPeaks(id = i, includeMotifOnly = motif_only_for_target_peak, TFregulome_url = gsub("api/table_query/", "", TFregulome_url)))
      if (is.null(peak_i))
      {
        message(paste0("... ... NO peak file for your id '", i,"'."))
      }
      else
      {
        target_list_count <- target_list_count + 1
        target_peak_list_all[[target_list_count]] <- peak_i
        TFregulome_target_peak_id <- c(TFregulome_target_peak_id, i)
        is_taregt_TFregulome <- c(is_taregt_TFregulome, T)
        message(paste0("... ... peak file loaded successfully for id '", i,"'"))
      }
    }
    message("... Done loading TFBS(s) from TFregulome")
  }
  # users' peaks
  if (!missing(user_target_peak_list) && length(user_target_peak_list)>0)
  {
    message(paste0("... You have ",length(user_target_peak_list)," customised peak set(s)"))
    if (missing(user_target_peak_id) || length(user_target_peak_id)!=length(user_target_peak_list) ||
        length(unique(user_target_peak_id))!=length(user_target_peak_list))
    {
      message("... ... You didn't provide the ID for each customised peak set or your ID number does not uniquely equal to the input user peak number. Instead we will use 'user_target_peak1', 'user_target_peak2'..." )
      user_target_peak_id <- paste0("user_target_peak", seq(1,length(user_target_peak_list), 1))
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
      # test if user input id i match any TFregulome database ID
      motif_matrix_i <- suppressMessages(searchMotif(id = user_target_peak_id[i], TFregulome_url = gsub("api/table_query/", "", TFregulome_url)))
      if (is.null(motif_matrix_i))
      {
        is_taregt_TFregulome <- c(is_taregt_TFregulome, F)
      }
      else
      {
        is_taregt_TFregulome <- c(is_taregt_TFregulome, T)
      }
    }
  }
  else
  {
    user_target_peak_id <- c()
  }
  # combine TFregulome ID and user ID
  target_peak_id_all <- c(TFregulome_target_peak_id, user_target_peak_id)

  # loading excluded peak list
  message("Loading excluded peak list ... ...")
  excluded_peak_list_all <- list()
  is_excluded_TFregulome <- c()
  # loading from TFregulome server
  excluded_list_count <- 0
  TFregulome_excluded_peak_id <- c()
  if (!missing(excluded_peak_id) && length(excluded_peak_id)>0)
  {
    message(paste0("... You have ", length(excluded_peak_id)," TFBS(s) requested to be loaded from TFregulome server"))
    if (motif_only_for_excluded_peak == T)
    {
      message("... You chose to load TF peaks with motif only. Using 'motif_only_for_excluded_peak' tunes your options")
    }
    else
    {
      message("... You chose to load TF peaks regardless of the presence of motif. Using 'motif_only_for_excluded_peak' tunes your options")
    }
    message("... loading TFBS(s) from TFregulome now")
    for (i in excluded_peak_id)
    {
      peak_i <- suppressMessages(loadPeaks(id = i, includeMotifOnly = motif_only_for_excluded_peak, TFregulome_url = gsub("api/table_query/", "", TFregulome_url)))
      if (is.null(peak_i))
      {
        message(paste0("... ... NO peak file for your id '", i,"'."))
      }
      else
      {
        excluded_list_count <- excluded_list_count + 1
        excluded_peak_list_all[[excluded_list_count]] <- peak_i
        is_excluded_TFregulome <- c(is_excluded_TFregulome, T)
        TFregulome_excluded_peak_id <- c(TFregulome_excluded_peak_id, i)
        message(paste0("... ... peak file loaded successfully for id '", i,"'"))
      }
    }
    message("... Done loading TFBS(s) from TFregulome")
  }
  # users' peak
  if (!missing(user_excluded_peak_list) && length(user_excluded_peak_list)>0)
  {
    message(paste0("... You have ",length(user_excluded_peak_list)," customised peak set(s)"))
    if (missing(user_excluded_peak_id) || length(user_excluded_peak_id)!=length(user_excluded_peak_list) ||
        length(unique(user_excluded_peak_id))!=length(user_excluded_peak_list))
    {
      message("... ... You didn't provide the ID for each customised peak set or your ID number does not uniquely equal to the input user peak number. Instead we will use 'user_excluded_peak1', 'user_excluded_peak2'..." )
      user_excluded_peak_id <- paste0("user_excluded_peak", seq(1,length(user_excluded_peak_list), 1))
    }
    for (i in 1:length(user_excluded_peak_list))
    {
      peak_i <- user_excluded_peak_list[[i]]
      peak_i_sub <- peak_i[,1:3]
      colnames(peak_i_sub) <- c("chr","start","end")
      peak_i_sub$id <- paste0("excluded_peak_", as.vector(rownames(peak_i_sub)))
      peak_i_sub <- peak_i_sub[,c("chr","start","end","id")]
      excluded_list_count <- excluded_list_count + 1
      excluded_peak_list_all[[excluded_list_count]] <- peak_i_sub
      # test if user input id i match any TFregulome database ID
      motif_matrix_i <- suppressMessages(searchMotif(id = user_excluded_peak_id[i], TFregulome_url = gsub("api/table_query/", "", TFregulome_url)))
      if (is.null(motif_matrix_i))
      {
        is_excluded_TFregulome <- c(is_excluded_TFregulome, F)
      }
      else
      {
        is_excluded_TFregulome <- c(is_excluded_TFregulome, T)
      }
    }
  }
  # start analysing
  exclusion_matrix <- list()
  for (i in 1:length(target_peak_list_all))
  {
    target_id_i <- target_peak_id_all[i]
    target_peak_i <- target_peak_list_all[[i]]
    number_of_orignal_target <- nrow(target_peak_i)
    message(paste0("Start analysing: ", target_id_i, "... ..."))

    ## if it is from TFregulome
    if (is_taregt_TFregulome[i])
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
        meth_file_200bp_path_target <- request_content_df[1,c("DNA_methylation_profile_200bp")]
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
    # excluding peaks
    for (j in 1:length(excluded_peak_list_all))
    {
      excluded_peak_j <- excluded_peak_list_all[[j]]
      if (isTFregulome_target)
      {
        bed_target_i <- with(target_peak_i, GRanges(chr, IRanges(start-99, end+100), id=id))
      }
      else
      {
        bed_target_i <- with(target_peak_i, GRanges(chr, IRanges(start, end), id=id))
      }
      if (is_excluded_TFregulome[j])
      {
        bed_excluded_j <- with(excluded_peak_j, GRanges(chr, IRanges(start-99, end+100), id=id))
      }
      else
      {
        bed_excluded_j <- with(excluded_peak_j, GRanges(chr, IRanges(start, end), id=id))
      }
      # get target peak intersecting with excluded peak
      # subsetOverlaps may mis-think the two sets coming from different references, so suppressWarnings here
      suppressWarnings(bedTarget_with_bedExcluded <- subsetByOverlaps(bed_target_i, bed_excluded_j))
      peakTarget_with_peakExcluded <- unique(as.data.frame(bedTarget_with_bedExcluded))
      peakTarget_without_peakExcluded <- target_peak_i[which(!(target_peak_i$id %in% peakTarget_with_peakExcluded$id)), ]
      target_peak_i <- peakTarget_without_peakExcluded
    }

    MethMotif_target <- new("MethMotif")
    #initiate methylation profile distribution
    meth_score_distri_target <- matrix()

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
        bed_target_done_exclusive <- with(target_peak_i[,c("chr","start","end","id")], GRanges(chr, IRanges(start-99, end+100), id=id))
        suppressWarnings(motif_of_bed_target_done_exclusive <- subsetByOverlaps(motif_seq_target_grange, bed_target_done_exclusive))
        motif_of_peakTarget_done_exclusive <- unique(as.data.frame(motif_of_bed_target_done_exclusive))
        if (nrow(motif_of_peakTarget_done_exclusive)>0)
        {
          motif_of_peakTarget_done_exclusive_allInfo <- motif_seq_target[which(motif_seq_target$id %in% motif_of_peakTarget_done_exclusive$id),]
          motif_matrix_of_peakTarget_done_exclusive <- formMatrixFromSeq(input_sequence = as.vector(motif_of_peakTarget_done_exclusive_allInfo$sequence),
                                                                      motif_format = motif_type)
          # compute beta score matrix
          if (isMethMotifID_target)
          {
            # methylation file can be empty
            meth_level_target <- tryCatch(read.delim(meth_file_path_target, sep = "\t", header = F),
                                          error=function(e) data.frame())
            if (nrow(meth_level_target)==0)
            {
              beta_score_matrix_of_peakTarget_done_exclusive <- formBetaScoreFromSeq(input_meth = data.frame(),
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
              suppressWarnings(meth_level_peakTarget_done_exclusive <- unique(as.data.frame(subsetByOverlaps(meth_level_target_grange,
                                                                                                             motif_of_bed_target_done_exclusive))))
              meth_level_peakTarget_done_exclusive_allInfo <- meth_level_target[which(meth_level_target$id %in% meth_level_peakTarget_done_exclusive$id),]
              beta_score_matrix_of_peakTarget_done_exclusive <- formBetaScoreFromSeq(input_meth = meth_level_peakTarget_done_exclusive_allInfo,
                                                                                  WGBS_replicate = WGBS_replicate_target,
                                                                                  motif_len = motif_len_target)
            }
          }
          else
          {
            beta_score_matrix_of_peakTarget_done_exclusive <- as.matrix(NA)
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
                                                    id = paste0(target_id_i,"_exclusive_peaks"),
                                                    alternate_name = target_id_i,
                                                    width = motif_len_target,
                                                    nsites=nrow(motif_of_peakTarget_done_exclusive),
                                                    motif_matrix=motif_matrix_of_peakTarget_done_exclusive)
          MethMotif_target@MMBetaScore <- beta_score_matrix_of_peakTarget_done_exclusive
        }

        # profile methylation level in narrow regions
        if (methylation_profile_in_narrow_region)
        {
          ### if in 200bp around peaks
          if (isMethMotifID_target)
          {
            meth_level_200bp_target <- tryCatch(read.delim(meth_file_200bp_path_target, sep = "\t", header = F),
                                                error=function(e) data.frame())
            if (nrow(meth_level_200bp_target) == 0)
            {
              meth_score_distri_target <- formBetaScoreDistri(input_meth = data.frame())
            }
            else
            {
              colnames(meth_level_200bp_target) <- c("chr","start","end",
                                                     "meth_score","C_num","T_num")
              meth_level_200bp_target$id <- paste0(target_id_i,"_200bp_CG_", as.vector(rownames(meth_level_200bp_target)))
              meth_level_200bp_target_grange <- with(meth_level_200bp_target[,c("chr","start","end","id")],
                                                     GRanges(chr, IRanges(start, end), id=id))
              bed_target_i <- with(target_peak_i, GRanges(chr, IRanges(start-99, end+100), id=id))
              suppressWarnings(meth_level_in_exclusive_peaks_200bp <- unique(as.data.frame(subsetByOverlaps(meth_level_200bp_target_grange,
                                                                                                         bed_target_i))))
              meth_level_in_exclusive_peaks_200bp_allInfo <- unique(meth_level_200bp_target[which(meth_level_200bp_target$id
                                                                                               %in% meth_level_in_exclusive_peaks_200bp$id),])
              meth_score_distri_target <- formBetaScoreDistri(input_meth = as.data.frame(meth_level_in_exclusive_peaks_200bp_allInfo$meth_score))
            }
          }
        }
      }
    }
    new_ExclusivePeaksMM <- new("ExclusivePeaksMM")
    new_ExclusivePeaksMM <- updateExclusivePeaksMM(theObject = new_ExclusivePeaksMM,
                                                   id=paste0(target_id_i, "_exclusive_peaks"),
                                                   exclusive_percentage = 100*nrow(target_peak_i)/number_of_orignal_target,
                                                   exclusive_peak = target_peak_i,
                                                   isTFregulomeID = isTFregulome_target,
                                                   MethMotif = MethMotif_target,
                                                   methylation_profile = meth_score_distri_target)
    exclusion_matrix[[paste0(target_id_i, "_exclusive_peaks")]] <- new_ExclusivePeaksMM
  }
  message("Done analysing.")
  dim(exclusion_matrix) <- c(length(target_peak_id_all), 1)
  rownames(exclusion_matrix) <- c(target_peak_id_all)
  colnames(exclusion_matrix) <- c("exclusive")
  return(exclusion_matrix)
}
