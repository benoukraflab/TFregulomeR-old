#' intersectPeakMatrix
#'
#' This function allows you to obtain the pair-wise intersected regions, along with the DNA methylation profiles, between two lists of peak sets, as well as (Meth)Motif logos
#' @param peak_id_x Character of vector, each of which is a unique TFregulome ID.
#' @param motif_only_for_id_x Either TRUE of FALSE (default). If TRUE, only peaks with motif will be loaded for each TFregulome ID in peak_id_x.
#' @param user_peak_list_x A list of data.frames, each of which contains user's own bed-format peak regions for peak list x.
#' @param user_peak_x_id Character of vector, each of which is a unique ID corresponding to each peak set in the list user_peak_list_x. If the IDs are not provided or not unique, the function will automatically generate the IDs of its own. If any of the peak sets is derived from TFregulome database, its TFregulome ID should be used here correspondingly.
#' @param peak_id_y Character of vector, each of which is a unique TFregulome ID.
#' @param motif_only_for_id_y Either TRUE of FALSE (default). If TRUE, only peaks with motif will be loaded for each TFregulome ID in peak_id_y.
#' @param user_peak_list_y A list of data.frames, each of which contains user's own bed-format peak regions for peak list y.
#' @param user_peak_y_id Character of vector, each of which is a unique ID corresponding to each peak set in the list user_peak_list_y. If the IDs are not provided or not unique, the function will automatically generate the IDs of its own. If any of the peak sets is derived from TFregulome database, its TFregulome ID should be used here correspondingly.
#' @param methylation_profile_in_narrow_region Either TRUE (default) of FALSE. If TRUE, methylation states in 200bp window surrounding peak summits for each intersected peak pair from peak_id_x (peak_id_y) and user_peak_list_x (user_peak_list_y) with TFregulome ID.
#' @param motif_type Motif PFM format, either in MEME by default or TRANSFAC.
#' @param TFregulome_url TFregulome server is implemented in MethMotif server. If the MethMoitf url is NO more "http://bioinfo-csi.nus.edu.sg/methmotif/", please use a new url.
#' @return  matrix of IntersectPeakMatrix class objects
#' @keywords intersectPeakMatrix
#' @export
#' @examples
#' peak_id_x <- c("MM1_HSA_K562_CEBPB", "MM1_HSA_HCT116_CEBPB")
#' peak_id_y <- c("MM1_HSA_HepG2_CEBPB", "MM1_HSA_HCT116_CEBPB")
#' intersectPeakMatrix_output <- intersectPeakMatrix(peak_id_x=peak_id_x,
#'                                                   motif_only_for_id_x=T,
#'                                                   peak_id_y=peak_id_y,
#'                                                   motif_only_for_id_y=T,
#'                                                   methylation_profile_in_narrow_region=T)


intersectPeakMatrix <- function(peak_id_x, motif_only_for_id_x = F, user_peak_list_x, user_peak_x_id,
                                peak_id_y, motif_only_for_id_y = F, user_peak_list_y, user_peak_y_id,
                                methylation_profile_in_narrow_region = F,
                                motif_type = "MEME", TFregulome_url)
{
  # check the input argument
  if (missing(peak_id_x) && missing(user_peak_list_x))
  {
    stop("No peak list x input. Please EITHER input TFregulome peaks using TFregulome ID(s) by 'peak_id_x = ' OR your own peak list using a list of data.frame(s) containing bed-format regions by 'user_peak_list_x = '")
  }
  if (missing(peak_id_y) && missing(user_peak_list_y))
  {
    stop("No peak list y input. Please EITHER input TFregulome peaks using TFregulome ID(s) by 'peak_id_y = ' OR your own peak list using a list of data.frame(s) containing bed-format regions by 'user_peak_list_y = '")
  }
  if ((!missing(user_peak_list_x) && class(user_peak_list_x) != "list") ||
      (!missing(user_peak_list_y) && class(user_peak_list_y) != "list"))
  {
    stop("The class of input 'user_peak_list_x' and 'user_peak_list_y' should be 'list', a list of bed-like data.frame storing peak regions!")
  }
  if (class(motif_only_for_id_x) != "logical" || class(motif_only_for_id_y) != "logical")
  {
    stop("motif_only_for_id_x and motif_only_for_id_y should be either TRUE or FALSE (default)")
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
    TFregulome_url <- "http://localhost:8888/api/table_query/"
  } else if (endsWith(TFregulome_url, suffix = "/index.php")==TRUE){
    TFregulome_url <- gsub("index.php", "", TFregulome_url)
    TFregulome_url <- paste0(TFregulome_url, "api/table_query/")
  } else if (endsWith(TFregulome_url, suffix = "/")==TRUE){
    TFregulome_url <- paste0(TFregulome_url, "api/table_query/")
  } else {
    TFregulome_url <- paste0(TFregulome_url, "/api/table_query/")
  }

  message("TFregulomeR::intersectPeakMatrix() starting ... ...")
  if (methylation_profile_in_narrow_region)
  {
    message("You chose to profile the methylation levels in 200bp window around peak summits, if there is any peak loaded from TFregulome. It will make the program slow. Disable it if you want a speedy analysis and do not care about methylation")
  }
  else
  {
    message("You chose NOT to profile the methylation levels in 200bp window around peak summits")
  }
  # loading peak list x
  message("Loading peak list x ... ...")
  peak_list_x_all <- list()
  # loading from TFregulome server
  TFregulome_peak_x_id <- c()
  is_x_TFregulome <- c()
  peak_list_x_count <- 0
  if (!missing(peak_id_x) && length(peak_id_x)>0)
  {
    message(paste0("... You have ", length(peak_id_x)," TFBS(s) requested to be loaded from TFregulome server"))
    if (motif_only_for_id_x == T)
    {
      message("... You chose to load TF peaks with motif only. Using 'motif_only_for_id_x' tunes your options")
    }
    else
    {
      message("... You chose to load TF peaks regardless of presence of motif. Using 'motif_only_for_id_x' tunes your options")
    }
    message("... loading TFBS(s) from TFregulome now")
    for (i in peak_id_x)
    {
      peak_i <- suppressMessages(loadPeaks(id = i, includeMotifOnly = motif_only_for_id_x, TFregulome_url = gsub("api/table_query/", "", TFregulome_url)))
      if (is.null(peak_i))
      {
        message(paste0("... ... NO peak file for your id '", i,"'."))
      }
      else
      {
        peak_list_x_count <- peak_list_x_count + 1
        peak_list_x_all[[peak_list_x_count]] <- peak_i
        is_x_TFregulome <- c(is_x_TFregulome, T)
        TFregulome_peak_x_id <- c(TFregulome_peak_x_id, i)
        message(paste0(".. ... peak file loaded successfully for id '", i,"'"))
      }
    }
    message("... Done loading TFBS(s) from TFregulome")
  }
  # users' peaks
  if (!missing(user_peak_list_x) && length(user_peak_list_x)>0)
  {
    message(paste0("... You have ",length(user_peak_list_x)," customised peak set(s)"))
    if (missing(user_peak_x_id) || length(user_peak_x_id)!=length(user_peak_list_x) ||
        length(unique(user_peak_x_id))!=length(user_peak_list_x))
    {
      message("... ... You didn't provide the ID for each customised peak set or your ID number does not uniquely equal to the input user peak number. Instead we will use 'user_peak_x1', 'user_peak_x2'..." )
      user_peak_x_id <- paste0("user_peak_x", seq(1,length(user_peak_list_x), 1))
    }
    for (i in 1:length(user_peak_list_x))
    {
      peak_i <- user_peak_list_x[[i]]
      peak_i_sub <- peak_i[,1:3]
      colnames(peak_i_sub) <- c("chr","start","end")
      peak_i_sub$id <- paste0(user_peak_x_id[i], "_", as.vector(rownames(peak_i_sub)))
      peak_i_sub <- peak_i_sub[,c("chr","start","end","id")]
      peak_list_x_count <- peak_list_x_count + 1
      peak_list_x_all[[peak_list_x_count]] <- peak_i_sub
      # test if user input id i match any TFregulome database ID
      motif_matrix_i <- suppressMessages(searchMotif(id = user_peak_x_id[i], TFregulome_url = gsub("api/table_query/", "", TFregulome_url)))
      if (is.null(motif_matrix_i))
      {
        is_x_TFregulome <- c(is_x_TFregulome, F)
      }
      else
      {
        is_x_TFregulome <- c(is_x_TFregulome, T)
      }
    }
  }
  else
  {
    user_peak_x_id <- c()
  }
  # combine TFregulome ID and user ID
  peak_id_x_all <- c(TFregulome_peak_x_id, user_peak_x_id)

  # loading compared peak list
  message("Loading peak list y ... ...")
  peak_list_y_all <- list()
  # loading from TFregulome server
  TFregulome_peak_y_id <- c()
  is_y_TFregulome <- c()
  peak_list_y_count <- 0
  if (!missing(peak_id_y) && length(peak_id_y)>0)
  {
    message(paste0("... You have ", length(peak_id_y)," TFBS(s) requested to be loaded from TFregulome server"))
    if (motif_only_for_id_y == T)
    {
      message("... You chose to load TF peaks with motif only. Using 'motif_only_for_id_y' tunes your options")
    }
    else
    {
      message("... You chose to load TF peaks regardless of presence of motif. Using 'motif_only_for_id_y' tunes your options")
    }
    message("... loading TFBS(s) from TFregulome now")
    for (i in peak_id_y)
    {
      peak_i <- suppressMessages(loadPeaks(id = i, includeMotifOnly = motif_only_for_id_y, TFregulome_url = gsub("api/table_query/", "", TFregulome_url)))
      if (is.null(peak_i))
      {
        message(paste0("... ... NO peak file for your id '", i,"'."))
      }
      else
      {
        peak_list_y_count <- peak_list_y_count + 1
        peak_list_y_all[[peak_list_y_count]] <- peak_i
        is_y_TFregulome <- c(is_y_TFregulome, T)
        TFregulome_peak_y_id <- c(TFregulome_peak_y_id, i)
        message(paste0(".. ... peak file loaded successfully for id '", i,"'"))
      }
    }
    message("... Done loading TFBS(s) from TFregulome")
  }
  # users' peaks
  if (!missing(user_peak_list_y) && length(user_peak_list_y)>0)
  {
    message(paste0("... You have ",length(user_peak_list_y)," customised peak set(s)"))
    if (missing(user_peak_y_id) || length(user_peak_y_id)!=length(user_peak_list_y) ||
        length(unique(user_peak_y_id))!=length(user_peak_list_y))
    {
      message("... ... You didn't provide the ID for each customised peak set or your ID number does not uniquely equal to the input user peak number. Instead we will use 'user_peak_y1', 'user_peak_y2'..." )
      user_peak_y_id <- paste0("user_peak_y", seq(1,length(user_peak_list_y), 1))
    }
    for (i in 1:length(user_peak_list_y))
    {
      peak_i <- user_peak_list_y[[i]]
      peak_i_sub <- peak_i[,1:3]
      colnames(peak_i_sub) <- c("chr","start","end")
      peak_i_sub$id <- paste0(user_peak_y_id[i], "_", as.vector(rownames(peak_i_sub)))
      peak_i_sub <- peak_i_sub[,c("chr","start","end","id")]
      peak_list_y_count <- peak_list_y_count + 1
      peak_list_y_all[[peak_list_y_count]] <- peak_i_sub
      # test if user input id i match any TFregulome database ID
      motif_matrix_i <- suppressMessages(searchMotif(id = user_peak_y_id[i], TFregulome_url = gsub("api/table_query/", "", TFregulome_url)))
      if (is.null(motif_matrix_i))
      {
        is_y_TFregulome <- c(is_y_TFregulome, F)
      }
      else
      {
        is_y_TFregulome <- c(is_y_TFregulome, T)
      }
    }
  }
  else
  {
    user_peak_y_id <- c()
  }
  # combine TFregulome ID and user ID
  peak_id_y_all <- c(TFregulome_peak_y_id, user_peak_y_id)


  intersection_matrix <- list()
  for (i in 1:length(peak_list_x_all))
  {
    id_x <- peak_id_x_all[i]
    message(paste0("Start analysing list x:", id_x, "... ..."))
    peak_x <- peak_list_x_all[[i]]
    if (is_x_TFregulome[i])
    {
      isTFregulome_x <- TRUE
      query_url <- paste0("listTFBS.php?AllTable=F&id=",id_x)
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
        isMethMotifID_x <- TRUE
        motif_seq_path_x <- request_content_df[1,c("TFBS")]
        meth_file_path_x <- request_content_df[1,c("DNA_methylation_profile")]
        meth_file_200bp_path_x <- request_content_df[1,c("DNA_methylation_profile_200bp")]
        WGBS_replicate_x <- request_content_df[1,c("WGBS_num")]
      }
      else
      {
        isMethMotifID_x <- FALSE
        motif_seq_path_x <- request_content_df[1,c("TFBS")]
      }
    }
    else
    {
      isTFregulome_x <- FALSE
    }
    # start comparing with peak set y
    for (j in 1:length(peak_list_y_all))
    {
      id_y <- peak_id_y_all[j]
      message(paste0("... ... Start analysing list y:", id_y))
      peak_y <- peak_list_y_all[[j]]
      if (is_y_TFregulome[j])
      {
        isTFregulome_y <- TRUE
        query_url <- paste0("listTFBS.php?AllTable=F&id=",id_y)
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
          isMethMotifID_y <- TRUE
          motif_seq_path_y <- request_content_df[1,c("TFBS")]
          meth_file_path_y <- request_content_df[1,c("DNA_methylation_profile")]
          meth_file_200bp_path_y <- request_content_df[1,c("DNA_methylation_profile_200bp")]
          WGBS_replicate_y <- request_content_df[1,c("WGBS_num")]
        }
        else
        {
          isMethMotifID_y <- FALSE
          motif_seq_path_y <- request_content_df[1,c("TFBS")]
        }
      }
      else
      {
        isTFregulome_y <- FALSE
      }

      # if peak x is from TRregulome database, extend peak regions by 100 bp
      if (isTFregulome_x)
      {
        bed_x <- with(peak_x, GRanges(chr, IRanges(start-99, end+100), id=id))
      }
      else
      {
        bed_x <- with(peak_x, GRanges(chr, IRanges(start, end), id=id))
      }
      # if peak y is from MethMotif database, extend peak regions by 100 bp
      if (isTFregulome_y)
      {
        bed_y <- with(peak_y, GRanges(chr, IRanges(start-99, end+100), id=id))
      }
      else
      {
        bed_y <- with(peak_y, GRanges(chr, IRanges(start, end), id=id))
      }

      # get peak x which intersect with y
      # subsetOverlaps may mis-think the two sets coming from different references, so suppressWarnings here
      suppressWarnings(bedx_with_bedy <- subsetByOverlaps(bed_x, bed_y))
      peakx_with_peaky <- unique(as.data.frame(bedx_with_bedy))
      x_interect_percentage <- 100*nrow(peakx_with_peaky)/nrow(peak_x)
      MethMotif_x <- new('MethMotif')
      # collecting all CpG meth scores in the defined methylation profile area
      meth_score_collection_x <- data.frame()
      # methylation distribution only meaningful if we have WGBS and peaks
      is_methProfile_meaningful_x <- F

      # form MethMotif object if the id is TFregulome id
      if (isTFregulome_x)
      {
        motif_seq_x <- read.delim(motif_seq_path_x, sep = "\t", header = F)
        if (nrow(peakx_with_peaky) > 0)
        {
          #compute motif matrix
          colnames(motif_seq_x) <- c("chr","start","end","strand","weight", "pvalue","qvalue","sequence")
          motif_len_x <- nchar(as.character(motif_seq_x[1,"sequence"]))
          motif_seq_x$id <- paste0(id_x,"_motif_sequence_", as.vector(rownames(motif_seq_x)))
          motif_seq_x_grange <- with(motif_seq_x[,c("chr","start","end","id")], GRanges(chr, IRanges(start+1, end), id=id))
          suppressWarnings(motif_of_peakx_with_peaky_grange <- subsetByOverlaps(motif_seq_x_grange, bedx_with_bedy))
          motif_of_peakx_with_peaky <- unique(as.data.frame(motif_of_peakx_with_peaky_grange))
          if (nrow(motif_of_peakx_with_peaky) > 0)
          {
            motif_of_peakx_with_peaky_allInfo <- motif_seq_x[which(motif_seq_x$id %in% motif_of_peakx_with_peaky$id),]
            motif_matrix_of_peakx_with_y <- formMatrixFromSeq(input_sequence = as.vector(motif_of_peakx_with_peaky_allInfo$sequence),
                                                              motif_format = motif_type)

            # compute beta score matrix,
            if (isMethMotifID_x)
            {
              # methylation file can be empty
              meth_level_x <- tryCatch(read.delim(meth_file_path_x, sep = "\t", header = F),
                                       error=function(e) data.frame())
              # methylation file can be empty
              if (nrow(meth_level_x)==0)
              {
                beta_score_matrix_of_peakx_with_y <- formBetaScoreFromSeq(input_meth = data.frame(),
                                                                          WGBS_replicate = WGBS_replicate_x,
                                                                          motif_len = motif_len_x)
              }
              else
              {
                colnames(meth_level_x) <- c("chr","start","end","meth_score","C_num","T_num","seq_chr","seq_start",
                                            "seq_end","strand","weight","pvalue","qvalue","sequence")
                meth_level_x$id <- paste0(id_x,"_motif_with_CG_", as.vector(rownames(meth_level_x)))
                meth_level_x_grange <- with(meth_level_x[,c("seq_chr","seq_start","seq_end","id")],
                                            GRanges(seq_chr, IRanges(seq_start, seq_end), id=id))
                suppressWarnings(meth_level_x_with_y <- unique(as.data.frame(subsetByOverlaps(meth_level_x_grange, motif_of_peakx_with_peaky_grange))))
                meth_level_x_with_y_allInfo <- meth_level_x[which(meth_level_x$id %in% meth_level_x_with_y$id),]
                beta_score_matrix_of_peakx_with_y <- formBetaScoreFromSeq(input_meth = meth_level_x_with_y_allInfo,
                                                                          WGBS_replicate = WGBS_replicate_x,
                                                                          motif_len = motif_len_x)
              }
            }
            else
            {
              beta_score_matrix_of_peakx_with_y <- as.matrix(NA)
            }

            if (motif_type == "TRANSFAC")
            {
              version_x <- 0
            }
            else
            {
              version_x <- 4
            }
            MethMotif_x@MMmotif <- updateMMmotif(MethMotif_x@MMmotif,
                                                 motif_format = motif_type,
                                                 version = version_x,
                                                 background = c("A"=0.25,"C"=0.25,"G"=0.25,"T"=0.25),
                                                 id = paste0(id_x,"_overlapped_with_", id_y),
                                                 width = motif_len_x,
                                                 nsites=nrow(motif_of_peakx_with_peaky),
                                                 motif_matrix=motif_matrix_of_peakx_with_y)
            MethMotif_x@MMBetaScore <- beta_score_matrix_of_peakx_with_y
          }

          # collecting CpG in x peaks that overlap with y peaks
          if (methylation_profile_in_narrow_region)
          {
            ### if in 200bp around peaks
            if (isMethMotifID_x)
            {
              is_methProfile_meaningful_x <- T
              meth_level_200bp_x <- tryCatch(read.delim(meth_file_200bp_path_x, sep = "\t", header = F),
                                                  error=function(e) data.frame())
              if (nrow(meth_level_200bp_x) > 0)
              {
                colnames(meth_level_200bp_x) <- c("chr","start","end",
                                                       "meth_score","C_num","T_num")
                meth_level_200bp_x$id <- paste0("200bp_CG_", as.vector(rownames(meth_level_200bp_x)))
                meth_level_200bp_x_grange <- with(meth_level_200bp_x[,c("chr","start","end","id")],
                                                       GRanges(chr, IRanges(start, end), id=id))
                suppressWarnings(meth_level_in_peakx_200bp <- unique(as.data.frame(subsetByOverlaps(meth_level_200bp_x_grange,
                                                                                                    bedx_with_bedy))))
                meth_level_in_peakx_200bp_allInfo <- unique(meth_level_200bp_x[which(meth_level_200bp_x$id
                                                                                     %in% meth_level_in_peakx_200bp$id),])
                meth_score_collection_x <- rbind(meth_score_collection_x,
                                               meth_level_in_peakx_200bp_allInfo[,c("chr","start","end","meth_score","C_num","T_num")])
              }
            }
          }
        }
      }
      # form methylation score profile - distribution for peak x
      if (is_methProfile_meaningful_x)
      {
        if (nrow(meth_score_collection_x)>0)
        {
          meth_score_distri_target_x <- formBetaScoreDistri(input_meth = as.data.frame(meth_score_collection_x$meth_score))
        }
        else
        {
          meth_score_distri_target_x <- formBetaScoreDistri(input_meth = data.frame())
        }
      }
      else
      {
        meth_score_distri_target_x <- matrix()
      }

      # get peak y which intersect with x
      suppressWarnings(bedy_with_bedx <- subsetByOverlaps(bed_y, bed_x))
      peaky_with_peakx <- unique(as.data.frame(bedy_with_bedx))
      y_interect_percentage <- 100*nrow(peaky_with_peakx)/nrow(peak_y)
      MethMotif_y <- new('MethMotif')
      # collecting all CpG meth scores in the defined methylation profile area
      meth_score_collection_y <- data.frame()
      # methylation distribution only meaningful if we have WGBS and peaks
      is_methProfile_meaningful_y <- F
      # form MethMotif object if the id is TFregulome id
      if (isTFregulome_y)
      {
        motif_seq_y <- read.delim(motif_seq_path_y, sep = "\t", header = F)
        if (nrow(peaky_with_peakx) > 0)
        {
          #compute motif matrix
          colnames(motif_seq_y) <- c("chr","start","end","strand","weight", "pvalue","qvalue","sequence")
          motif_len_y <- nchar(as.character(motif_seq_y[1,"sequence"]))
          motif_seq_y$id <- paste0(id_y,"_motif_sequence_", as.vector(rownames(motif_seq_y)))
          motif_seq_y_grange <- with(motif_seq_y[,c("chr","start","end","id")], GRanges(chr, IRanges(start+1, end), id=id))
          suppressWarnings(motif_of_peaky_with_peakx_grange <- subsetByOverlaps(motif_seq_y_grange, bedy_with_bedx))
          motif_of_peaky_with_peakx <- unique(as.data.frame(motif_of_peaky_with_peakx_grange))
          if (nrow(motif_of_peaky_with_peakx) > 0)
          {
            motif_of_peaky_with_peakx_allInfo <- motif_seq_y[which(motif_seq_y$id %in% motif_of_peaky_with_peakx$id),]
            motif_matrix_of_peaky_with_x <- formMatrixFromSeq(input_sequence = as.vector(motif_of_peaky_with_peakx_allInfo$sequence),
                                                              motif_format = motif_type)

            # compute beta score matrix,
            if (isMethMotifID_y)
            {
              # methylation file can be empty
              meth_level_y <- tryCatch(read.delim(meth_file_path_y, sep = "\t", header = F),
                                       error=function(e) data.frame())
              # methylation file can be empty
              if (nrow(meth_level_y)==0)
              {
                beta_score_matrix_of_peaky_with_x <- formBetaScoreFromSeq(input_meth = data.frame(),
                                                                          WGBS_replicate = WGBS_replicate_y,
                                                                          motif_len = motif_len_y)
              }
              else
              {
                colnames(meth_level_y) <- c("chr","start","end","meth_score","C_num","T_num","seq_chr","seq_start",
                                            "seq_end","strand","weight","pvalue","qvalue","sequence")
                meth_level_y$id <- paste0(id_y,"_motif_with_CG_", as.vector(rownames(meth_level_y)))
                meth_level_y_grange <- with(meth_level_y[,c("seq_chr","seq_start","seq_end","id")],
                                            GRanges(seq_chr, IRanges(seq_start, seq_end), id=id))
                suppressWarnings(meth_level_y_with_x <- unique(as.data.frame(subsetByOverlaps(meth_level_y_grange, motif_of_peaky_with_peakx_grange))))
                meth_level_y_with_x_allInfo <- meth_level_y[which(meth_level_y$id %in% meth_level_y_with_x$id),]
                beta_score_matrix_of_peaky_with_x <- formBetaScoreFromSeq(input_meth = meth_level_y_with_x_allInfo,
                                                                          WGBS_replicate = WGBS_replicate_y,
                                                                          motif_len = motif_len_y)
              }
            }
            else
            {
              beta_score_matrix_of_peaky_with_x <- as.matrix(NA)
            }

            if (motif_type == "TRANSFAC")
            {
              version_y <- 0
            }
            else
            {
              version_y <- 4
            }
            MethMotif_y@MMmotif <- updateMMmotif(MethMotif_y@MMmotif,
                                                 motif_format = motif_type,
                                                 version = version_y,
                                                 background = c("A"=0.25,"C"=0.25,"G"=0.25,"T"=0.25),
                                                 id = paste0(id_y,"_overlapped_with_", id_x),
                                                 width = motif_len_y,
                                                 nsites=nrow(motif_of_peaky_with_peakx),
                                                 motif_matrix=motif_matrix_of_peaky_with_x)
            MethMotif_y@MMBetaScore <- beta_score_matrix_of_peaky_with_x
          }

          # collecting CpG in y peaks that overlap with x peaks
          if (methylation_profile_in_narrow_region)
          {
            ### if in 200bp around peaks
            if (isMethMotifID_y)
            {
              is_methProfile_meaningful_y <- T
              meth_level_200bp_y <- tryCatch(read.delim(meth_file_200bp_path_y, sep = "\t", header = F),
                                             error=function(e) data.frame())
              if (nrow(meth_level_200bp_y) > 0)
              {
                colnames(meth_level_200bp_y) <- c("chr","start","end",
                                                  "meth_score","C_num","T_num")
                meth_level_200bp_y$id <- paste0("200bp_CG_", as.vector(rownames(meth_level_200bp_y)))
                meth_level_200bp_y_grange <- with(meth_level_200bp_y[,c("chr","start","end","id")],
                                                  GRanges(chr, IRanges(start, end), id=id))
                suppressWarnings(meth_level_in_peaky_200bp <- unique(as.data.frame(subsetByOverlaps(meth_level_200bp_y_grange,
                                                                                                    bedy_with_bedx))))
                meth_level_in_peaky_200bp_allInfo <- unique(meth_level_200bp_y[which(meth_level_200bp_y$id
                                                                                     %in% meth_level_in_peaky_200bp$id),])
                meth_score_collection_y <- rbind(meth_score_collection_y,
                                               meth_level_in_peaky_200bp_allInfo[,c("chr","start","end","meth_score","C_num","T_num")])
              }
            }
          }
        }
      }
      # form methylation score profile - distribution
      if (is_methProfile_meaningful_y)
      {
        if (nrow(meth_score_collection_y)>0)
        {
          meth_score_distri_target_y <- formBetaScoreDistri(input_meth = as.data.frame(meth_score_collection_y$meth_score))
        }
        else
        {
          meth_score_distri_target_y <- formBetaScoreDistri(input_meth = data.frame())
        }
      }
      else
      {
        meth_score_distri_target_y <- matrix()
      }

      #form an IntersectPeakMatrix object
      new_IntersectPeakMatrix <- new("IntersectPeakMatrix")
      new_IntersectPeakMatrix <- updateIntersectPeakMatrix(theObject = new_IntersectPeakMatrix,
                                                           id = paste0(id_x,"_[AND]_",id_y),
                                                           id_x = id_x,
                                                           overlap_percentage_x = x_interect_percentage,
                                                           isxTFregulomeID = isTFregulome_x,
                                                           MethMotif_x = MethMotif_x,
                                                           methylation_profile_x = meth_score_distri_target_x,
                                                           id_y = id_y,
                                                           overlap_percentage_y = y_interect_percentage,
                                                           isyTFregulomeID = isTFregulome_y,
                                                           MethMotif_y = MethMotif_y,
                                                           methylation_profile_y = meth_score_distri_target_y)
      intersection_matrix[[paste0(id_x,"_[AND]_",id_y)]] <- new_IntersectPeakMatrix
    }
  }
  intersection_matrix_matrix <- matrix(intersection_matrix, nrow = length(peak_id_x_all), ncol = length(peak_id_y_all), byrow = TRUE)
  rownames(intersection_matrix_matrix) <- c(peak_id_x_all)
  colnames(intersection_matrix_matrix) <- c(peak_id_y_all)
  return(intersection_matrix_matrix)
}


formMatrixFromSeq <- function(input_sequence, motif_format)
{
  input_sequence <- as.data.frame(input_sequence)
  input_sequence_matrix <- data.frame(do.call('rbind', strsplit(as.character(input_sequence$input_sequence),'',fixed=TRUE)))

  motif_matrix_TRANSFAC <- matrix(rep(-1, 4*ncol(input_sequence_matrix)), ncol = 4)
  alphabet <- c("A","C","G","T")
  colnames(motif_matrix_TRANSFAC) <- alphabet

  for (i in seq(1,4,1)){
    for (j in seq(1,ncol(input_sequence_matrix),1)){
      motif_matrix_TRANSFAC[j, i] <- sum(input_sequence_matrix[,j] == alphabet[i])
    }
  }
  motif_matrix_MEME <- motif_matrix_TRANSFAC/nrow(input_sequence_matrix)

  if (motif_format == "MEME"){
    return(motif_matrix_MEME)
  }
  else{
    return(motif_matrix_TRANSFAC)
  }
}


formBetaScoreFromSeq <- function(input_meth, WGBS_replicate, motif_len)
{
  if (nrow(input_meth) == 0)
  {
    empty_matrix <- TRUE
  }
  else
  {
    if (WGBS_replicate=="2"){
      input_meth <- unique(input_meth[,c("chr","start","end","meth_score","C_num","T_num","seq_chr",
                                         "seq_start","seq_end","strand","sequence")])
      input_meth_d <- input_meth[which(input_meth$strand=="+"),]
      input_meth_r <- input_meth[which(input_meth$strand=="-"),]

      input_meth_d$dis <- input_meth_d$start-input_meth_d$seq_start+1
      input_meth_r$dis <- motif_len-(input_meth_r$start-input_meth_r$seq_start+1)

      if(nrow(input_meth_d)==0 && nrow(input_meth_r)==0){
        empty_matrix <- TRUE
      } else if(nrow(input_meth_d)==0 && nrow(input_meth_r)!=0){
        input_meth_sub <- input_meth_r[,c("dis","meth_score")]
        empty_matrix <- FALSE
      }  else if(nrow(input_meth_d)!=0 && nrow(input_meth_r)==0){
        input_meth_sub <- input_meth_d[,c("dis","meth_score")]
        empty_matrix <- FALSE
      } else if(nrow(input_meth_d)!=0 && nrow(input_meth_r)!=0){
        input_meth_d_sub <- input_meth_d[,c("dis","meth_score")]
        input_meth_r_sub <- input_meth_r[,c("dis","meth_score")]
        input_meth_sub <- rbind(input_meth_d_sub, input_meth_r_sub)
        empty_matrix <- FALSE
      }
    }
    else{
      input_meth <- unique(input_meth[,c("chr","start","end","meth_score","C_num","T_num","seq_chr",
                                         "seq_start","seq_end","strand","sequence")])
      input_meth_d <- input_meth[which(input_meth$strand=="+"),]
      input_meth_r <- input_meth[which(input_meth$strand=="-"),]

      input_meth_d$dis <- input_meth_d$start-input_meth_d$seq_start+1
      input_meth_r$dis <- motif_len-(input_meth_r$start-input_meth_r$seq_start)
      # merge read in both strands
      if(nrow(input_meth_d)>0){
        for (i in 1:nrow(input_meth_d)){
          if (unlist(strsplit(as.character(input_meth_d[i,"sequence"]), split=""))[as.integer(input_meth_d[i,"dis"])]=="G"){
            input_meth_d[i,"dis"] <- input_meth_d[i,"dis"]-1
          }
        }
      }
      # merge read in both strands
      if(nrow(input_meth_r)>0){
        for (i in 1:nrow(input_meth_r)){
          if (unlist(strsplit(as.character(input_meth_r[i,"sequence"]), split=""))[as.integer(input_meth_r[i,"dis"])]=="G"){
            input_meth_r[i,"dis"] <- input_meth_r[i,"dis"]-1
          }
        }
      }
      # calculate overall beta score in both strand
      if(nrow(input_meth_d)>0){
        input_meth_d$id <- paste(input_meth_d$seq_chr,input_meth_d$seq_start,input_meth_d$seq_end,input_meth_d$dis,sep = "")
        input_meth_d_sub <- data.frame()
        input_meth_d_id_uniq <- unique(input_meth_d$id)
        for (i in input_meth_d_id_uniq){
          input_meth_d_temp <- input_meth_d[which(input_meth_d$id==i),c("C_num","T_num","dis")]
          dis_temp <- input_meth_d_temp[1,3]
          meth_temp <- 100*sum(input_meth_d_temp[,1])/sum(input_meth_d_temp[,c(1,2)])
          new_add <- data.frame(i,dis_temp,meth_temp)
          input_meth_d_sub <- rbind(input_meth_d_sub, new_add)
        }
      }
      # calculate overall beta score in both strand
      if(nrow(input_meth_r)>0){
        input_meth_r$id <- paste(input_meth_r$seq_chr,input_meth_r$seq_start,input_meth_r$seq_end,input_meth_r$dis,sep = "")
        input_meth_r_sub <- data.frame()
        input_meth_r_id_uniq <- unique(input_meth_r$id)
        for (i in input_meth_r_id_uniq){
          input_meth_r_temp <- input_meth_r[which(input_meth_r$id==i),c("C_num","T_num","dis")]
          dis_temp <- input_meth_r_temp[1,3]
          meth_temp <- 100*sum(input_meth_r_temp[,1])/sum(input_meth_r_temp[,c(1,2)])
          new_add <- data.frame(i,dis_temp,meth_temp)
          input_meth_r_sub <- rbind(input_meth_r_sub, new_add)
        }
      }

      if(nrow(input_meth_d)==0 && nrow(input_meth_r)==0){
        empty_matrix <- TRUE
      } else if(nrow(input_meth_d)==0 && nrow(input_meth_r)!=0){
        input_meth_sub <- input_meth_r_sub[,2:3]
        colnames(input_meth_sub) <- c("dis","meth_score")
        empty_matrix <- FALSE
      } else if(nrow(input_meth_d)!=0 && nrow(input_meth_r)==0){
        input_meth_sub <- input_meth_d_sub[,2:3]
        colnames(input_meth_sub) <- c("dis","meth_score")
        empty_matrix <- FALSE
      } else if(nrow(input_meth_d)!=0 && nrow(input_meth_r)!=0){
        input_meth_sub <- rbind(input_meth_d_sub[,2:3], input_meth_r_sub[,2:3])
        colnames(input_meth_sub) <- c("dis","meth_score")
        empty_matrix <- FALSE
      }

    }
  }

  if(empty_matrix==TRUE){
    plot_matrix <- data.frame(matrix(0,nrow = 3, ncol = motif_len))
    colnames(plot_matrix) <- c(seq(1,motif_len,1))
  } else{
    plot_matrix <- data.frame(c(1,2,3))
    for (i in seq(1,motif_len,1)){
      input_meth_sub_i <- input_meth_sub[which(input_meth_sub$dis == i),]
      unmeth_count <- nrow(input_meth_sub_i[which(input_meth_sub_i$meth_score<10),])
      meth_count <- nrow(input_meth_sub_i[which(input_meth_sub_i$meth>90),])
      inbetween_count <- nrow(input_meth_sub_i)-unmeth_count-meth_count
      new_column <- c(unmeth_count,inbetween_count,meth_count)
      plot_matrix <- data.frame(plot_matrix,new_column)

    }
    plot_matrix <- plot_matrix[,2:ncol(plot_matrix)]
    colnames(plot_matrix) <- c(seq(1,motif_len,1))
  }

  return(as.matrix(plot_matrix))
}
