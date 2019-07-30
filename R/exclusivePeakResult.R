#' exclusivePeaks result
#'
#' This function allows you to get the reuslts from the exclusivePeaks() output, including a list of exclusive peak sets, (Meth)Motif logo, methylation profile in exclsive peaks and exclsive peak summary.
#' @param exclusivePeaks Required. exclusivePeaks() output, a matrix of ExclusivePeaksMM class objects.
#' @param return_exclusive_peak_sites Either TRUE of FALSE (default). If TRUE, a list of data.frames containing exclusive peak sets will be returned.
#' @param save_MethMotif_logo Either TRUE of FALSE (default). If TRUE, (Meth)Motif logos for the exclusive peak sets will be saved.
#' @param return_methylation_profile Either TRUE of FALSE (default). If TRUE, the methylation profiles in 200bp window around exclusive peak summits will be returned.
#' @param return_summary Either TRUE of FALSE (default). If TRUE, a numeric matrix containing the percentage of peaks as exclusive will be returned.
#' @param logo_type Logo type for the (Meth)Motif logo to be saved, either "entropy" (default) or "frequency".
#' @param meth_level Methylation level to be plotted for the (Meth)Motif logo, and it should be one of the values, "all" (default), "methylated", and "unmethylated".
#' @return  a list of data.frames, a numeric matrix or (Meth)Motif logo PDF files depending on the options.
#' @keywords exclusivePeakResult
#' @export
#' @examples
#' target_id <- "MM1_HSA_K562_CEBPB"
#' excluded_id <- c("MM1_HSA_HepG2_CEBPB", "MM1_HSA_HCT116_CEBPB")
#' excluPeak_output <- exclusivePeaks(target_peak_id=target_id,
#'                                    motif_only_for_target_peak=TRUE,
#'                                    excluded_peak_id=excluded_id,
#'                                    motif_only_for_excluded_peak=TRUE,
#'                                    methylation_profile_in_narrow_region=TRUE)
#' exclusivePeaks_result <- exclusivePeakResult(exclusivePeaks=excluPeak_output,
#'                                              return_exclusive_peak_sites=TRUE,
#'                                              save_MethMotif_logo=TRUE,
#'                                              return_summary=TRUE)

exclusivePeakResult <- function(exclusivePeaks,
                                return_exclusive_peak_sites = FALSE,
                                save_MethMotif_logo = FALSE,
                                return_methylation_profile = FALSE,
                                return_summary = FALSE,
                                logo_type="entropy",
                                meth_level="all")
{
  # check input arguments
  if (missing(exclusivePeaks))
  {
    stop("please provide results from 'exclusivePeaks()' using 'exclusivePeaks ='!")
  }
  if (!is.logical(return_exclusive_peak_sites))
  {
    stop("'return_exclusive_peak_sites' should be either TRUE (T) or FALSE (F, default)!")
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
  message("Start getting the results of exclusivePeaks ...")
  if (return_exclusive_peak_sites == TRUE)
  {
    message("... ... You chose to return exclusive peak sites;")
    if (return_summary == TRUE && return_methylation_profile == FALSE)
    {
      message("... ... You chose NOT to return methylation profile;")
      message("... ... You chose to return exclusive peak summary;")
      message("... ... ... Both exclusive peak sets and peak summary will be stored in a list, and named with 'exclusive_peak_list' and 'peak_summary' in the list. Use 'names()' in the output for its details.")
    }
    else if (return_summary == FALSE && return_methylation_profile == TRUE)
    {
      message("... ... You chose to return methylation profile;")
      message("... ... You chose NOT to return exclusive peak summary;")
      message("... ... ... Both exclusive peak sets and methylation profiles will be stored in a list, and named with 'exclusive_peak_list' and 'methylation_profile' in the list. Use 'names()' in the output for its details.")
    }
    else if (return_summary == FALSE && return_methylation_profile == FALSE)
    {
      message("... ... You chose NOT to return methylation profile;")
      message("... ... You chose NOT to return exclusive peak summary;")
      message("... ... ... Only exclusive peak sets will be returned in a list and named with 'exclusive_peak_list'.")
    }
    else
    {
      message("... ... You chose to return methylation profile;")
      message("... ... You chose to return exclusive peak summary;")
      message("... ... ... ALL of exclusive peak sets, methylation profiles and peak summary will be stored in a list, and named with 'exclusive_peak_list', 'methylation_profile' and 'peak_summary' in the list. Use 'names()' in the output for its details.")
    }
  }
  else
  {
    message("... ... You chose NOT to return exclusive peak sites;")
    if (return_summary == TRUE && return_methylation_profile == FALSE)
    {
      message("... ... You chose NOT to return methylation profile;")
      message("... ... You chose to return exclusive peak summary;")
      message("... ... ... Only peak summary will be returned in a list and named with 'peak_summary'.")
    }
    else if (return_summary == FALSE && return_methylation_profile == TRUE)
    {
      message("... ... You chose to return methylation profile;")
      message("... ... You chose NOT to return exclusive peak summary;")
      message("... ... ... Only methylation_profiles will be returned in a list, and named with 'methylation_profile'.")
    }
    else if (return_summary == FALSE && return_methylation_profile == FALSE)
    {
      message("... ... You chose NOT to return methylation profile;")
      message("... ... You chose NOT to return exclusive peak summary;")
      message("... ... ... None will be returned.")
    }
    else
    {
      message("... ... You chose to return methylation profile;")
      message("... ... You chose to return exclusive peak summary;")
      message("... ... ... Both methylation profiles and peak summary will be stored in a list, and named with 'methylation_profile' and 'peak_summary' in the list. Use 'names()' in the output for its details.")
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

  if (return_exclusive_peak_sites == FALSE && save_MethMotif_logo == FALSE && return_summary == FALSE && return_methylation_profile==FALSE)
  {
    message("... ... You chose no action. EXIT!!")
    return(NULL)
  }
  else
  {
    # if save methmotif logo
    if (save_MethMotif_logo)
    {
      for (i in seq(1, nrow(exclusivePeaks), 1))
      {
        logo_id <- exclusivePeaks[i,1][[1]]@id
        is_TFregulome <- exclusivePeaks[i,1][[1]]@isTFregulomeID
        nsites <- exclusivePeaks[i,1][[1]]@MethMotif@MMmotif@nsites
        if (is_TFregulome == TRUE && nsites > 0)
        {
          plotLogo(MM_object = exclusivePeaks[i,1][[1]]@MethMotif, logo_type = logo_type, meth_level = meth_level)
        }
        else
        {
          message(paste0("... ... ... The input peak set for the results '",logo_id,"' was not originated from TFregulomeR or the number of direct binding sites in the exclusive peaks is 0, so no motif logo available."))
        }
      }
    }
    # store all the results in the following
    return_all <- list()
    # if return exclusive peaks
    if (return_exclusive_peak_sites)
    {
      exclusive_peak_list <- list()
      for (i in seq(1,nrow(exclusivePeaks), 1))
      {
        exclusive_peak_i <- as.data.frame(exclusivePeaks[i,1][[1]]@exclusive_peak)
        peak_id <- exclusivePeaks[i,1][[1]]@id
        exclusive_peak_list[[peak_id]] <- exclusive_peak_i
      }
      return_all[["exclusive_peak_list"]] <- exclusive_peak_list
    }
    # if return methylation profile
    if (return_methylation_profile)
    {
      methylation_profile <- list()
      for (i in seq(1,nrow(exclusivePeaks), 1))
      {
        meth_profile_i <- as.matrix(exclusivePeaks[i,1][[1]]@methylation_profile)
        peak_id <- exclusivePeaks[i,1][[1]]@id
        methylation_profile[[peak_id]] <- meth_profile_i
      }
      return_all[["methylation_profile"]] <- methylation_profile
    }
    # if return summary
    if (return_summary)
    {
      summary_matrix <- matrix(nrow = nrow(exclusivePeaks), ncol = 1)
      id_list <- c()
      for (i in seq(1, nrow(exclusivePeaks), 1))
      {
        summary_id <- exclusivePeaks[i,1][[1]]@id
        id_list <- c(id_list, summary_id)
        summary_matrix[i,1] <- exclusivePeaks[i,1][[1]]@exclusive_percentage
      }
      rownames(summary_matrix) <- id_list
      colnames(summary_matrix) <- c("percentage_in_original_inputs(%)")
      return_all[["peak_summary"]] <- summary_matrix
    }
    # return values
    if (!(return_exclusive_peak_sites == FALSE && return_methylation_profile == FALSE && return_summary == FALSE))
    {
      return(return_all)
    }
  }
}



