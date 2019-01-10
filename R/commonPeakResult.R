#' commonPeaks result
#'
#' This function allows you to get the reuslts from the commonPeaks() output, including a list of common peak sets, (Meth)Motif logos, methylation profile in common peak and common peak summary.
#' @param commonPeaks Required. commonPeaks() output, a matrix of CommonPeaksMM class objects.
#' @param return_common_peak_sites Either TRUE of FALSE (default). If TRUE, a list of data.frames containing common peak sets.
#' @param save_MethMotif_logo Either TRUE of FALSE (default). If TRUE, MethMotif logos for the common peak sets will be saved.
#' @param return_methylation_profile Either TRUE of FALSE (default). If TRUE, methylation profile in common peak sets will be returned.
#' @param return_summary Either TRUE of FALSE (default). If TRUE, a numeric matrix containing the percentage of peaks as common will be returned.
#' @param logo_type Logo type for the MethMotif logo to be saved, either "entropy" (default) or "frequency",
#' @param meth_level Methylation level to be plotted for the (Meth)Motif logo, and it should be one of the values, "all" (default), "methylated", and "unmethylated".
#' @return  a list of data.frames, a numeric matrix or (Meth)Motif logo PDF files depending on the options.
#' @keywords commonPeakResult
#' @export
#' @examples
#' target_id <- c("MM1_HSA_K562_CEBPB")
#' compared_id <- c("MM1_HSA_HepG2_CEBPB")
#' commonPeaks_output <- commonPeaks(target_peak_id=target_id,
#'                                   motif_only_for_target_peak=T,
#'                                   compared_peak_id=compared_id,
#'                                   motif_only_for_compared_peak=T,
#'                                   methylation_profile_in_narrow_region=T)
#' commonPeaks_result <- commonPeakResult(commonPeaks=commonPeaks_output,
#'                                        return_common_peak_sites=T,
#'                                        save_MethMotif_logo=T,
#'                                        return_methylation_profile=T,
#'                                        return_summary=T)

commonPeakResult <- function(commonPeaks, return_common_peak_sites = FALSE,
                             save_MethMotif_logo = FALSE, return_methylation_profile = FALSE,
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
  if (!is.logical(return_methylation_profile))
  {
    stop("'return_methylation_profile' should be either TRUE (T) or FALSE (F, default)!")
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
    if (return_summary == TRUE && return_methylation_profile == FALSE)
    {
      message("... ... You chose NOT to return methylation profile;")
      message("... ... You chose to return common peak summary;")
      message("... ... ... Both common peak sets and peak summary will be stored in a list, and named with 'common_peak_list' and 'peak_summary' in the list. Use 'names()' in the output for its detials.")
    }
    else if (return_summary == FALSE && return_methylation_profile == TRUE)
    {
      message("... ... You chose to return methylation profile;")
      message("... ... You chose NOT to return common peak summary;")
      message("... ... ... Both common peak sets and methylation profiles will be stored in a list, and named with 'common_peak_list' and 'methylation_profile' in the list. Use 'names()' in the output for its detials.")
    }
    else if (return_summary == FALSE && return_methylation_profile == FALSE)
    {
      message("... ... You chose NOT to return methylation profile;")
      message("... ... You chose NOT to return common peak summary;")
      message("... ... ... Only common peak sets will be returned in a list and named with 'common_peak_list'.")
    }
    else
    {
      message("... ... You chose to return methylation profile;")
      message("... ... You chose to return common peak summary;")
      message("... ... ... ALL of common peak sets, methylation profiles and peak summary will be stored in a list, and named with 'common_peak_list', 'methylation_profile' and 'peak_summary' in the list. Use 'names()' in the output for its detials.")
    }
  }
  else
  {
    message("... ... You chose NOT to return common peak sites;")
    if (return_summary == TRUE && return_methylation_profile == FALSE)
    {
      message("... ... You chose NOT to return methylation profile;")
      message("... ... You chose to return common peak summary;")
      message("... ... ... Only peak summary will be returned in a list and named with 'peak_summary'.")
    }
    else if (return_summary == FALSE && return_methylation_profile == TRUE)
    {
      message("... ... You chose to return methylation profile;")
      message("... ... You chose NOT to return common peak summary;")
      message("... ... ... Only methylation_profiles will be returned in a list, and named with 'methylation_profile'.")
    }
    else if (return_summary == FALSE && return_methylation_profile == FALSE)
    {
      message("... ... You chose NOT to return methylation profile;")
      message("... ... You chose NOT to return common peak summary;")
      message("... ... ... None will be returned.")
    }
    else
    {
      message("... ... You chose to return methylation profile;")
      message("... ... You chose to return common peak summary;")
      message("... ... ... Both methylation profiles and peak summary will be stored in a list, and named with 'methylation_profile' and 'peak_summary' in the list. Use 'names()' in the output for its detials.")
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

  if (return_common_peak_sites == FALSE && save_MethMotif_logo == FALSE && return_summary == FALSE && return_methylation_profile==FALSE)
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
    # store all the results in the following
    return_all <- list()
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
      return_all[["common_peak_list"]] <- common_peak_list
    }
    # if return methylation profile
    if (return_methylation_profile)
    {
      methylation_profile <- list()
      for (i in 1:nrow(commonPeaks))
      {
        meth_profile_i <- as.matrix(commonPeaks[i,1][[1]]@methylation_profile)
        peak_id <- commonPeaks[i,1][[1]]@id
        methylation_profile[[peak_id]] <- meth_profile_i
      }
      return_all[["methylation_profile"]] <- methylation_profile
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
      return_all[["peak_summary"]] <- summary_matrix
    }
    # return values
    if (!(return_common_peak_sites == F && return_methylation_profile == F && return_summary == F))
    {
      return(return_all)
    }
  }
}
