#' intersectPeakMatrix result
#'
#' This function allows you to get the reuslts from the intersectPeakMatrix() output, including a matrix of the pair-wise intersecting percentages between two lists of peak sets, DNA methylation profiles in the intersected regions (for peaks from TFregulomeR) and (Meth)Motif logos for each pair of intersections (for peaks from TFregulomeR).
#' @param intersectPeakMatrix Required. intersectPeakMatrix() output, a marix of IntersectPeakMatrix class objects.
#' @param return_intersection_matrix Either TRUE of FALSE (default). If TRUE, a matrix of the pair-wise intersecting percentages between two lists of peak sets will be returned.
#' @param angle_of_matrix Either "x" (default) or "y". If "x", a matrix denoting the percentages of peak sets in "peak_list_x" intersected with "peak_list_y" will be returned; if "y", a matrix denoting the percentages of peak sets in "peak_list_y" intersected with "peak_list_x" will be returned.
#' @param return_tag_density Either TRUE of FALSE (default). If TRUE, a matrix of tag density values in intersected peaks between "peak_list_x" and "peak_list_y" will be returned.
#' @param angle_of_tag_density Either "x" (default) or "y". If "x", a matrix denoting tag density values in "peak_list_x" intersected with "peak_list_y" will be returned; if "y", a matrix denoting tag density values in "peak_list_y" intersected with "peak_list_x" will be returned.
#' @param tag_density_value The value of tag density in the intersected peaks. It should be one of the following values: "median" (default),"mean","SD","quartile_25","quartile_75".
#' @param return_external_source Either TRUE of FALSE (default). If TRUE, a matrix of external source values in intersected peaks between "peak_list_x" and "peak_list_y" will be returned.
#' @param angle_of_external_source Either "x" (default) or "y". If "x", a matrix denoting external source values in "peak_list_x" intersected with "peak_list_y" will be returned; if "y", a matrix denoting tag density values in "peak_list_y" intersected with "peak_list_x" will be returned.
#' @param external_source_value The value of external source signal in the intersected peaks. It should be one of the following values: "median" (default),"mean","SD","quartile_25","quartile_75".
#' @param return_methylation_profile Either TRUE of FALSE (default). If TRUE, a matrix of DNA methylation state data in the intersected regions will be returned.
#' @param angle_of_methylation_profile Either "x" (default) or "y". If "x", a matrix denoting DNA methylation state data in "peak_list_x" intersected with "peak_list_y" will be returned; if "y", a matrix denoting DNA methylation state data in "peak_list_y" intersected with "peak_list_x" will be returned.
#' @param save_MethMotif_logo Either TRUE of FALSE (default). If TRUE, (Meth)Motif logos for the intersected peaks will be saved.
#' @param angle_of_logo Either "x" (default) or "y". If "x", (Meth)Motif logos for the peak sets in "peak_list_x" intersected with "peak_list_y" will be saved; if "y", (Meth)Motif logos for the peak sets in "peak_list_y" intersected with "peak_list_x" will be saved.
#' @param logo_type Logo type for the (Meth)Motif logo to be saved, either "entropy" (default) or "frequency".
#' @param meth_level Methylation level to be plot for the (Meth)Motif logo, and it should be one of the values, "all" (default), "methylated", and "unmethylated".
#' @param saving_MethMotif_logo_x_id Either "all" (default) or a subset of "peak_id_x". If a subset of "peak_id_x" is provided, only the (Meth)Motif logos for them will be saved.
#' @param saving_MethMotif_logo_y_id Either "all" (default) or a subset of "peak_id_y". If a subset of "peak_id_y" is provided, only the (Meth)Motif logos for them will be saved.
#' @return  a matrix of pair-wise intersecting percantages between two lists of peak sets, a matrix of DNA methylation data in the intersected regions or (Meth)Motif PDF files depending on the options
#' @keywords intersectPeakMatrixResult
#' @export
#' @examples
#' peak_id_x <- c("MM1_HSA_K562_CEBPB", "MM1_HSA_HCT116_CEBPB")
#' peak_id_y <- c("MM1_HSA_HepG2_CEBPB", "MM1_HSA_HCT116_CEBPB")
#' intersect_output <- intersectPeakMatrix(peak_id_x=peak_id_x,
#'                                                   motif_only_for_id_x=TRUE,
#'                                                   peak_id_y=peak_id_y,
#'                                                   motif_only_for_id_y=TRUE)
#' intersect_matrix <- intersectPeakMatrixResult(intersectPeakMatrix=intersect_output,
#'                                               return_intersection_matrix=TRUE,
#'                                               save_MethMotif_logo=TRUE,
#'                           saving_MethMotif_logo_x_id=c("MM1_HSA_K562_CEBPB"))

intersectPeakMatrixResult <- function(intersectPeakMatrix,
                                      return_intersection_matrix = FALSE,
                                      angle_of_matrix = "x",
                                      return_tag_density = FALSE,
                                      angle_of_tag_density = "x",
                                      tag_density_value = "median",
                                      return_external_source = FALSE,
                                      angle_of_external_source = "x",
                                      external_source_value = "median",
                                      return_methylation_profile = FALSE,
                                      angle_of_methylation_profile = "x",
                                      save_MethMotif_logo = FALSE,
                                      angle_of_logo="x",
                                      logo_type="entropy",
                                      meth_level="all",
                                      saving_MethMotif_logo_x_id = "all",
                                      saving_MethMotif_logo_y_id = "all")
{
  # check input arguments
  if (missing(intersectPeakMatrix))
  {
    stop("Please provide results from 'intersectPeakMatrix()' using 'intersectPeakMatrix ='!")
  }
  if (!is.logical(return_intersection_matrix))
  {
    stop("'return_intersection_matrix' should be either TRUE (T) or FALSE (F, default)!")
  }
  if (!is.logical(return_methylation_profile))
  {
    stop("'return_methylation_profile' should be either TRUE (T) or FALSE (F, default)!")
  }
  if (!is.logical(return_tag_density))
  {
    stop("'return_tag_density' should be either TRUE (T) or FALSE (F, default)!")
  }
  if (!is.logical(return_external_source))
  {
    stop("'return_external_source' should be either TRUE (T) or FALSE (F, default)!")
  }
  if (angle_of_matrix != "x" && angle_of_matrix != "y")
  {
    stop("'angle_of_matrix' should be either 'x' (default) or 'y'!")
  }
  if (angle_of_methylation_profile != "x" && angle_of_methylation_profile != "y")
  {
    stop("'angle_of_methylation_profile' should be either 'x' (default) or 'y'!")
  }
  if (angle_of_tag_density != "x" && angle_of_tag_density != "y")
  {
    stop("'angle_of_tag_density' should be either 'x' (default) or 'y'!")
  }
  if (angle_of_external_source != "x" && angle_of_external_source != "y")
  {
    stop("'angle_of_external_source' should be either 'x' (default) or 'y'!")
  }
  if (!(tag_density_value %in% c("median","mean","SD","quartile_25","quartile_75")))
  {
    stop("'tag_density_value' should be one of the following values: 'median','mean','SD','quartile_25','quartile_75'")
  }
  if (!(external_source_value %in% c("median","mean","SD","quartile_25","quartile_75")))
  {
    stop("'external_source_value' should be one of the following values: 'median','mean','SD','quartile_25','quartile_75'")
  }
  if (!is.logical(save_MethMotif_logo))
  {
    stop("'save_MethMotif_logo' should be either TRUE (T) or FALSE (F, default)!")
  }
  if (angle_of_logo != "x" && angle_of_logo != "y")
  {
    stop("'angle_of_logo' should be either 'x' (default) or 'y'!")
  }
  if (logo_type != "entropy" && logo_type != "frequency")
  {
    stop("'logo_type' should be either 'entropy' (default) or 'frequency'!")
  }
  if (meth_level != "all" && meth_level != "methylated" && meth_level != "unmethylated")
  {
    stop("'meth_level' should be one of 'all' (default), 'methylated' and 'unmethylated'!")
  }
  if (length(saving_MethMotif_logo_x_id) > nrow(intersectPeakMatrix))
  {
    stop("number of x ids input in 'saving_MethMotif_logo_x_id' is larger than number of rows in 'intersectPeakMatrix'!!")
  }
  if (length(saving_MethMotif_logo_y_id) > ncol(intersectPeakMatrix))
  {
    stop("number of y ids input in 'saving_MethMotif_logo_y_id' is larger than number of columns in 'intersectPeakMatrix'!!")
  }

  # output the arguments
  message("Start getting the results of intersectPeakMatrix ...")
  if (return_intersection_matrix == TRUE)
  {
    message("... ... You chose to return intersection matrix;")
    if (angle_of_matrix == "x")
    {
      message("... ... ... You chose x-wise intersection matrix;")
    }
    else
    {
      message("... ... ... You chose y-wise intersection matrix;")
    }
  }
  else
  {
    message("... ... You chose NOT to return intersection matrix;")
  }

  if (return_tag_density == TRUE)
  {
    message("... ... You chose to return tag density;")
    message(paste0("... ... ... the tag density value you chose to return is ",
                   tag_density_value))
    if (angle_of_tag_density == "x")
    {
      message("... ... ... You chose to return tag density for peak list x;")
    }
    else
    {
      message("... ... ... You chose to return tag density for peak list y;")
    }
  }
  else
  {
    message("... ... You chose NOT to return tag density;")
  }

  if (return_external_source == TRUE)
  {
    message("... ... You chose to return external source signal;")
    message(paste0("... ... ... the external source signal value you chose to return is ",
                   external_source_value))
    if (angle_of_external_source == "x")
    {
      message("... ... ... You chose to return external source signal for peak list x;")
    }
    else
    {
      message("... ... ... You chose to return external source signal for peak list y;")
    }
  }
  else
  {
    message("... ... You chose NOT to return external source signal;")
  }

  if (return_methylation_profile == TRUE)
  {
    message("... ... You chose to return methylation profile;")
    if (angle_of_methylation_profile == "x")
    {
      message("... ... ... You chose to return methylation profile for peak list x;")
    }
    else
    {
      message("... ... ... You chose to return methylation profile for peak list y;")
    }
  }
  else
  {
    message("... ... You chose NOT to return methylation profile;")
  }

  if (save_MethMotif_logo == TRUE)
  {
    message("... ... You chose to save MethMotif logo in PDF if any;")
    if (angle_of_logo == "x")
    {
      message("... ... ... You chose x-wise MethMotif logo;")
    }
    else
    {
      message("... ... ... You chose y-wise MethMotif logo;")
    }
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


  if (return_intersection_matrix == FALSE && save_MethMotif_logo==FALSE &&
      return_methylation_profile ==FALSE && return_tag_density == FALSE &&
      return_external_source == FALSE)
  {
    message("... ... You chose no action. EXIT!!")
    return(NULL)
  }
  else
  {
    # if save methmotif logo
    if (save_MethMotif_logo)
    {
      for (i in seq(1, nrow(intersectPeakMatrix), 1))
      {
        if (saving_MethMotif_logo_x_id != "all" && !(rownames(intersectPeakMatrix)[i] %in% saving_MethMotif_logo_x_id))
        {
          next
        }
        for (j in seq(1, ncol(intersectPeakMatrix), 1))
        {
          if (saving_MethMotif_logo_y_id != "all" && !(colnames(intersectPeakMatrix)[j] %in% saving_MethMotif_logo_y_id))
          {
            next
          }
          if (angle_of_logo == "x")
          {
            logo_id <- intersectPeakMatrix[i,j][[1]]@id_x
            is_TFregulome <- intersectPeakMatrix[i,j][[1]]@isxTFregulomeID
            nsites <- intersectPeakMatrix[i,j][[1]]@MethMotif_x@MMmotif@nsites
            if (is_TFregulome==TRUE && nsites>0)
            {
              plotLogo(MM_object = intersectPeakMatrix[i,j][[1]]@MethMotif_x, logo_type = logo_type, meth_level = meth_level)
            }
            else
            {
              message(paste0("No (Meth)Motif logo for ", logo_id, " will be generated, because its ID does match any record of existing TFregulomeR IDs or number of the interected motif is zero."))
            }
          }
          else
          {
            logo_id <- intersectPeakMatrix[i,j][[1]]@id_y
            is_TFregulome <- intersectPeakMatrix[i,j][[1]]@isyTFregulomeID
            nsites <- intersectPeakMatrix[i,j][[1]]@MethMotif_y@MMmotif@nsites
            if (is_TFregulome==TRUE && nsites>0)
            {
              plotLogo(MM_object = intersectPeakMatrix[i,j][[1]]@MethMotif_y, logo_type = logo_type, meth_level = meth_level)
            }
            else
            {
              message(paste0("No (Meth)Motif logo for ", logo_id, " will be generated, because its ID does match any record of existing TFregulomeR IDs or number of the interected motif is zero."))
            }
          }
        }
      }
    }
    return_all <- list()
    # if return intersection matrix
    if (return_intersection_matrix)
    {
      intersection_matrix <- data.frame(matrix(nrow = nrow(intersectPeakMatrix), ncol = ncol(intersectPeakMatrix)))
      rownames(intersection_matrix) <- rownames(intersectPeakMatrix)
      colnames(intersection_matrix) <- colnames(intersectPeakMatrix)
      for (i in seq(1, nrow(intersectPeakMatrix), 1))
      {
        for (j in seq(1, ncol(intersectPeakMatrix), 1))
        {
          if (angle_of_matrix == "x")
          {
            intersection_matrix[i,j] <- intersectPeakMatrix[i,j][[1]]@overlap_percentage_x
          }
          else
          {
            intersection_matrix[i,j] <- intersectPeakMatrix[i,j][[1]]@overlap_percentage_y
          }
        }
      }
      return_all[["intersection_matrix"]] <- intersection_matrix
    }
    # if return tag density
    if (return_tag_density)
    {
      tag_density_matrix <- data.frame(matrix(nrow = nrow(intersectPeakMatrix),
                                              ncol = ncol(intersectPeakMatrix)))
      rownames(tag_density_matrix) <- rownames(intersectPeakMatrix)
      colnames(tag_density_matrix) <- colnames(intersectPeakMatrix)
      for (i in seq(1, nrow(intersectPeakMatrix), 1))
      {
        for (j in seq(1, ncol(intersectPeakMatrix), 1))
        {
          if (angle_of_tag_density == "x")
          {
            tag_density_matrix[i,j] <- intersectPeakMatrix[i,j][[1]]@tag_density_x[tag_density_value]
          }
          else
          {
            tag_density_matrix[i,j] <- intersectPeakMatrix[i,j][[1]]@tag_density_y[tag_density_value]
          }
        }
      }
      return_all[["tag_density_matrix"]] <- tag_density_matrix
    }
    # if return external source signal
    if (return_external_source)
    {
      external_source_matrix <- data.frame(matrix(nrow = nrow(intersectPeakMatrix),
                                                  ncol = ncol(intersectPeakMatrix)))
      rownames(external_source_matrix) <- rownames(intersectPeakMatrix)
      colnames(external_source_matrix) <- colnames(intersectPeakMatrix)
      for (i in seq(1, nrow(intersectPeakMatrix), 1))
      {
        for (j in seq(1, ncol(intersectPeakMatrix), 1))
        {
          if (angle_of_tag_density == "x")
          {
            external_source_matrix[i,j] <- intersectPeakMatrix[i,j][[1]]@external_signal_x[external_source_value]
          }
          else
          {
            external_source_matrix[i,j] <- intersectPeakMatrix[i,j][[1]]@external_signal_y[external_source_value]
          }
        }
      }
      return_all[["external_source_matrix"]] <- external_source_matrix
    }
    # if return methylation profile
    if (return_methylation_profile)
    {
      methylation_profile_list <- list()
      count = 1
      for (i in seq(1, nrow(intersectPeakMatrix), 1))
      {
        for (j in seq(1, ncol(intersectPeakMatrix), 1))
        {

          if (angle_of_methylation_profile == "x")
          {
            id <- intersectPeakMatrix[i,j][[1]]@id_x
            methylation_profile_list[[count]] <- intersectPeakMatrix[i,j][[1]]@methylation_profile_x
            count <- count +1
          }
          else
          {
            id <- intersectPeakMatrix[i,j][[1]]@id_y
            methylation_profile_list[[count]] <- intersectPeakMatrix[i,j][[1]]@methylation_profile_y
            count <- count +1
          }
        }
      }
      methylation_profile_matrix <- matrix(methylation_profile_list, nrow = nrow(intersectPeakMatrix),
                                           ncol = ncol(intersectPeakMatrix), byrow = TRUE)
      rownames(methylation_profile_matrix) <- rownames(intersectPeakMatrix)
      colnames(methylation_profile_matrix) <- colnames(intersectPeakMatrix)
      return_all[["methylation_profile_matrix"]] <- methylation_profile_matrix
    }
    #return value
    if (!(return_intersection_matrix==FALSE && return_methylation_profile==FALSE &&
          return_tag_density==FALSE && return_external_source==FALSE))
    {
      return(return_all)
    }
  }
}
