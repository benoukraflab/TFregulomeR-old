#' commonPeaks result
#'
#' This function allows you to get the reuslts from the commonPeaks() output, including a list of common peak sets, MethMotif logos, and common peak summary.
#' @param commonPeaks Required. commonPeaks() output, a matrix of CommonPeaksMM class objects.
#' @param return_common_peak_sites Optional. Either TRUE of FALSE (default). If TRUE, a list of data.frames containing common peak sets derived from "target_peak_list" will be returned. If "return_summary" is also TRUE, both will be returned in a list().
#' @param save_MethMotif_logo Optional. Either TRUE of FALSE (default). If TRUE, MethMotif logos for the common peak sets will be saved.
#' @param return_summary Optional. Either TRUE of FALSE (default). If TRUE, a numeric matrix containing the percentage of peaks in "target_peak_list" as common will be returned. If "return_common_peak_sites" is also TRUE, both will be returned in a list().
#' @param logo_type Optional. Logo type for the MethMotif logo to be saved, either "entropy" (default) or "frequency",
#' @param meth_level Optional. Methylation level to be plot for the MethMotif logo, and it should be one of the values, "all" (default), "methylated", and "unmethylated".
#' @return  a list of data.frames, a numeric matrix or MethMotif logo PDF files depending on the options.
#' @keywords commonPeakResult
#' @export
#' @examples
#' target_peaks <- list(loadPeaks(id = "MM1_HSA_K562_CEBPB"),
#'                      read.delim("my_own_peaks.bed", header = F))
#' target_id <- c("MM1_HSA_K562_CEBPB", "my_own_peaks")
#' compared_peaks <- list(loadPeaks(id = "MM1_HSA_HepG2_CEBPB"),
#'                       read.delim("peaks_to_common_with.bed", header = F))
#' compared_id <- c("MM1_HSA_HepG2_CEBPB", "peaks_to_common_with")
#' commonPeaks_output <- commonPeaks(target_peak_list=target_peaks,
#'                                   target_peak_id=target_id,
#'                                   compared_peak_list=compared_peaks,
#'                                   compared_peak_id=compared_id)
#' commonPeaks_result <- commonPeakResult(commonPeaks=commonPeaks_output,
#'                                        return_exclusive_peak_sites=T,
#'                                        save_MethMotif_logo=T,
#'                                        return_summary=T)

commonPeakResult <- function(commonPeaks, return_common_peak_sites = FALSE, save_MethMotif_logo = FALSE,
                                return_summary = FALSE, logo_type="entropy", meth_level="all")
{
  # check input arguments
  if (missing(commonPeaks))
  {
    stop("please provide results from 'commonPeaks()' using 'commonPeaks ='!")
  }
  if (!is.logical(return_common_peak_sites))
  {
    stop("'return_common_peak_sites' should be either TRUE (T) or FALSE (F, default)!")
  }
  if (!is.logical(save_MethMotif_logo))
  {
    stop("'save_MethMotif_logo' should be either TRUE (T) or FALSE (F, default)!")
  }
  if (!is.logical(return_summary))
  {
    stop("'return_summary' should be either TRUE (T) or FALSE (F, default)!")
  }
  if (logo_type != "entropy" && logo_type != "frequency")
  {
    stop("'logo_type' should be either 'entropy' (default) or 'frequency'!")
  }
  if (meth_level != "all" && meth_level != "methylated" && meth_level != "unmethylated")
  {
    stop("'meth_level' should be one of 'all' (default), 'methylated' and 'unmethylated'!")
  }

  # output the arguments
  message("Start getting the results of commonPeakResult ...")
  if (return_common_peak_sites == TRUE)
  {
    message("... ... You chose to return common peak sites;")
    if (return_summary == TRUE)
    {
      message("... ... You chose to return common peak summary;")
      message("... ... ... Both common peak sets and peak summary will be stored in a list, and named with 'common_peak_list' and 'peak_summary' in the list. Use 'names()' in the output for its detials.")
    }
    else
    {
      message("... ... You chose NOT to return common peak summary;")
      message("... ... ... Only list of common peak sets will be returned in a list. Use names() in the output to see the ids of each peak set in the list.")
    }
  }
  else
  {
    message("... ... You chose NOT to return common peak sites;")
    if (return_summary == TRUE)
    {
      message("... ... You chose to return common peak summary;")
      message("... ... ... Only peak summary will be returned in a matrix.")
    }
    else
    {
      message("... ... You chose NOT to return common peak summary;")
    }
  }

  if (save_MethMotif_logo == TRUE)
  {
    message("... ... You chose to save MethMotif logo in PDF if any;")
    if (logo_type == "entropy")
    {
      message("... ... ... You chose entropy logo;")
    }
    else
    {
      message("... ... ... You chose frequency logo;")
    }
    if (meth_level == "all")
    {
      message("... ... ... You chose to show all methylation levels;")
    }
    else if (meth_level == "methylated")
    {
      message("... ... ... You chose to show the methylated only;")
    }
    else
    {
      message("... ... ... You chose to show the unmethylated only;")
    }
  }
  else
  {
    message("... ... You chose NOT to save MethMotif logo in PDF if any;")
  }

  if (return_common_peak_sites == FALSE && save_MethMotif_logo == FALSE && return_summary == FALSE)
  {
    message("... ... You chose no action. EXIT!!")
    return(NULL)
  }
  else
  {
    # if save methmotif logo
    if (save_MethMotif_logo)
    {
      for (i in 1:nrow(commonPeaks))
      {
        logo_id <- commonPeaks[i,1][[1]]@id
        is_TFregulome <- commonPeaks[i,1][[1]]@isTFregulomeID
        nsites <- commonPeaks[i,1][[1]]@MethMotif@MMmotif@nsites
        if (is_TFregulome == TRUE && nsites > 0)
        {
          plotLogo(MM_object = commonPeaks[i,1][[1]]@MethMotif, logo_type = logo_type, meth_level = meth_level)
        }
        else
        {
          message(paste0("... ... ... The input peak set for the results '",logo_id,"' was not orginated from TFregulome or the number of direct binding sites in the common peaks is 0, so no motif logo available."))
        }
      }
    }
    # if return common peaks
    if (return_common_peak_sites)
    {
      common_peak_list <- list()
      for (i in 1:nrow(commonPeaks))
      {
        common_peak_i <- as.data.frame(commonPeaks[i,1][[1]]@common_peak)
        peak_id <- commonPeaks[i,1][[1]]@id
        common_peak_list[[peak_id]] <- common_peak_i
      }
    }
    # if return summary
    if (return_summary)
    {
      summary_matrix <- matrix(nrow = nrow(commonPeaks), ncol = 1)
      id_list <- c()
      for (i in 1:nrow(commonPeaks))
      {
        summary_id <- commonPeaks[i,1][[1]]@id
        id_list <- c(id_list, summary_id)
        summary_matrix[i,1] <- commonPeaks[i,1][[1]]@common_percentage
      }
      rownames(summary_matrix) <- id_list
      colnames(summary_matrix) <- c("percentage_in_original_inputs(%)")
    }
    # return values
    if (return_common_peak_sites == TRUE && return_summary == TRUE)
    {
      return_all <- list()
      return_all[["common_peak_list"]] <- common_peak_list
      return_all[["peak_summary"]] <- summary_matrix
      return(return_all)
    }
    if (return_common_peak_sites == TRUE && return_summary == FALSE)
    {
      return(common_peak_list)
    }
    if (return_common_peak_sites == FALSE && return_summary == TRUE)
    {
      return(summary_matrix)
    }
  }
}
