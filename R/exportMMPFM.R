#' export motif position frequency matrix and beta score matrix
#'
#' This function allows you to export motif position frequency matrix and beta score matrix from the output of "searchMotif", "commonPeaks", "exclusivePeaks" or "intersectPeakMatrix".
#' @param fun_output Required. Output from "searchMotif", "commonPeaks", "exclusivePeaks" or "intersectPeakMatrix".
#' @param fun Required. The function that was used to get the output and should be one of the options, 'searchMotif', 'commonPeaks', 'exclusivePeaks' and 'intersectPeakMatrix'.
#' @param save_motif_PFM Either TRUE or FALSE (default). If "TRUE", the motif position frequency matrix will be saved.
#' @param save_betaScore_matrix Either TRUE or FALSE (default). If "TRUE", the beta score matrix will be saved.
#' @param angle_of_matrix_for_intersectPeakMatrix Only applicable when "fun='intersectPeakMatrix'". Either "x" (default) or "y". If "x", motif PFM for the peak sets in "peak_list_x" intersected with "peak_list_y" will be saved; if "y", motif PFM for the peak sets in "peak_list_y" intersected with "peak_list_x" will be saved.
#' @param saving_id_x_for_intersectPeakMatrix Only applicable when "fun='intersectPeakMatrix'". Either "all" (default) or a subset of "peak_id_x". If a subset of "peak_id_x" is provided, only the MethMotif logos for them will be saved.
#' @param saving_id_y_for_intersectPeakMatrix Only applicable when "fun='intersectPeakMatrix'". Either "all" (default) or a subset of "peak_id_y". If a subset of "peak_id_y" is provided, only the MethMotif logos for them will be saved.
#' @return  motif position frequency matrix file and beta score matrix file
#' @keywords exportMMPFM
#' @export
#' @examples
#' methmotif_cebpb <- searchMotif(id = "MM1_HSA_K562_CEBPB")
#' exportMMPFM(fun_output = methmotif_cebpb, fun = "searchMotif",
#'             save_motif_PFM = TRUE, save_betaScore_matrix = TRUE)

exportMMPFM <- function(fun_output, fun, save_motif_PFM = FALSE,
                        save_betaScore_matrix = FALSE,
                        angle_of_matrix_for_intersectPeakMatrix = "x",
                        saving_id_x_for_intersectPeakMatrix = "all",
                        saving_id_y_for_intersectPeakMatrix = "all")
{
  # check input arguments
  if (missing(fun_output))
  {
    stop("Please provide a TFregulomeR function output using 'fun_output ='!")
  }
  if (missing(fun))
  {
    stop("Please tell us which TFregulomeR function you used to get the output using 'fun ='! It should be one of them, 'searchMotif', 'commonPeaks', 'exclusivePeaks'and 'intersectPeakMatrix'!")
  }
  if (!(fun %in% c('searchMotif', 'commonPeaks', 'exclusivePeaks', 'intersectPeakMatrix')))
  {
    stop("Sorry, we cannot recognise the function you input using 'fun ='! exportMMPFM only works for the output from 'searchMotif', 'commonPeaks', 'exclusivePeaks' or 'intersectPeakMatrix', and please choose one of them as input!")
  }
  if (!is.logical(save_motif_PFM))
  {
    stop("'save_motif_PFM' should be either TRUE (T) or FALSE (F, default)!")
  }
  if (!is.logical(save_betaScore_matrix))
  {
    stop("'save_betaScore_matrix' should be either TRUE (T) or FALSE (F, default)!")
  }
  if (angle_of_matrix_for_intersectPeakMatrix != "x" && angle_of_matrix_for_intersectPeakMatrix != "y")
  {
    stop("'angle_of_matrix_for_intersectPeakMatrix' should be either 'x' or 'y'!")
  }
  if (fun == "intersectPeakMatrix" && length(saving_id_x_for_intersectPeakMatrix) > nrow(fun_output))
  {
    stop("You chose to export the result from intersectPeakMatrix(), but the number of x ids you input in 'saving_id_x_for_intersectPeakMatrix' is larger than the row number in your 'fun_output'")
  }
  if (fun == "intersectPeakMatrix" && length(saving_id_y_for_intersectPeakMatrix) > ncol(fun_output))
  {
    stop("You chose to export the result from intersectPeakMatrix(), but the number of x ids you input in 'saving_id_y_for_intersectPeakMatrix' is larger than the column number in your 'fun_output'")
  }

  message("Start exporting ... ...")
  if (save_motif_PFM == TRUE && save_betaScore_matrix == TRUE)
  {
    message("... ... You chose to save motif PFM and beta score matrix.")
  }
  else if (save_motif_PFM == TRUE && save_betaScore_matrix == FALSE)
  {
    message("... ... You chose to save motif PFM only.")
  }
  else if (save_motif_PFM == FALSE && save_betaScore_matrix == TRUE)
  {
    message("... ... You chose to save motif beta score matrix only.")
  }
  else
  {
    message("... ... You chose NOT to save motif PFM OR beta score matrix.")
    message("... ... No action. EXIT!!")
    return(NULL)
  }

  # start exporting
  ## for searchMotif
  if (fun == "searchMotif")
  {
    message("... ... export searchMotif")
    if (class(fun_output)[1] != "MethMotif")
    {
      message(paste0("... ... Error: searchMotif output should be MethMotif class objects."))
      message("... ... Please make sure your input is searchMotif output and you have provided the correct information in 'fun ='!")
      message("EXIT!!")
      return(NULL)
    }
    # if save betaScore matrix
    if (save_betaScore_matrix)
    {
      beta_score_matrix <- fun_output@MMBetaScore
      if (is.na(beta_score_matrix[1,1]))
      {
        message("... ... ... ... No beta score matrix is available. Skip!")
      }
      else
      {
        saved_BetaScore <- save_BetaScore(MethMotif_object = fun_output)
        message(paste0("... ... ... ... Beta score matrix has been saved as '", saved_BetaScore,"'."))
      }
    }
    # if save motif PFM
    if (save_motif_PFM)
    {
      saved_motifPFM <- save_motifPFM(MethMotif_object = fun_output)
      message(paste0("... ... ... ... Motif PFM has been saved as '", saved_motifPFM,"'."))
    }
  }
  ## for commonPeaks
  if (fun == "commonPeaks")
  {
    message("... ... export commonPeaks...")
    for (i in 1:nrow(fun_output))
    {
      fun_output_i <- fun_output[i, 1][[1]]
      if (class(fun_output_i)[1] != "CommonPeaksMM")
      {
        message(paste0("Error: commonPeaks output is a matrix of CommonPeaksMM class objects. We found the element row No.", i," is not CommonPeaksMM class in your input!"))
        message("... ... Please make sure your input is commonPeaks output and you have provided the correct information in 'fun ='!")
        message("EXIT!!")
        return(NULL)
      }
      id_i <- fun_output_i@id
      message(paste0("... ... ... export id = ", id_i))
      isTFregulomeID <- fun_output_i@isTFregulomeID
      nsites <- fun_output_i@MethMotif@MMmotif@nsites
      if (isTFregulomeID == TRUE && nsites > 0)
      {
        # if save betaScore matrix
        if (save_betaScore_matrix)
        {
          beta_score_matrix <- fun_output_i@MethMotif@MMBetaScore
          if (is.na(beta_score_matrix[1,1]))
          {
            message("... ... ... ... No beta score matrix is available. Skip!")
          }
          else
          {
            saved_BetaScore <- save_BetaScore(MethMotif_object = fun_output_i@MethMotif)
            message(paste0("... ... ... ... Beta score matrix has been saved as '", saved_BetaScore,"'."))
          }
        }
        # if save motif PFM
        if (save_motif_PFM)
        {
          saved_motifPFM <- save_motifPFM(MethMotif_object = fun_output_i@MethMotif)
          message(paste0("... ... ... ... Motif PFM has been saved as '", saved_motifPFM,"'."))
        }
      }
      else
      {
        message("... ... ... the original peaks of ",id_i, " is not loaded from TFregulome database, or in the common peak the number of TFBS is zero. Hence no further action for this id!")
      }
    }
  }
  ## for exclusivePeaks
  if (fun == "exclusivePeaks")
  {
    message("... ... export exclusivePeaks...")
    for (i in 1:nrow(fun_output))
    {
      fun_output_i <- fun_output[i, 1][[1]]
      if (class(fun_output_i)[1] != "ExclusivePeaksMM")
      {
        message(paste0("Error: exclusivePeaks output is a matrix of ExclusivePeaksMM class objects. We found the element row No.", i," is not ExclusivePeaksMM class in your input!"))
        message("... ... Please make sure your input is exclusivePeaks output and you have provided the correct information in 'fun ='!")
        message("EXIT!!")
        return(NULL)
      }
      id_i <- fun_output_i@id
      message(paste0("... ... ... export id = ", id_i))
      isTFregulomeID <- fun_output_i@isTFregulomeID
      nsites <- fun_output_i@MethMotif@MMmotif@nsites
      if (isTFregulomeID == TRUE && nsites > 0)
      {
        # if save betaScore matrix
        if (save_betaScore_matrix)
        {
          beta_score_matrix <- fun_output_i@MethMotif@MMBetaScore
          if (is.na(beta_score_matrix[1,1]))
          {
            message("... ... ... ... No beta score matrix is available. Skip!")
          }
          else
          {
            saved_BetaScore <- save_BetaScore(MethMotif_object = fun_output_i@MethMotif)
            message(paste0("... ... ... ... Beta score matrix has been saved as '", saved_BetaScore,"'."))
          }
        }
        # if save motif PFM
        if (save_motif_PFM)
        {
          saved_motifPFM <- save_motifPFM(MethMotif_object = fun_output_i@MethMotif)
          message(paste0("... ... ... ... Motif PFM has been saved as '", saved_motifPFM,"'."))
        }
      }
      else
      {
        message("... ... ... the original peaks of ",id_i, " is not loaded from MethMotif database, or in the exclusive peaks the number of TFBS recorded in MethMotif database is zero. Hence no further action for this id!")
      }
    }
  }
  ## for intersectPeakMatrix
  if (fun == "intersectPeakMatrix")
  {
    message("... ... export intersectPeakMatrix...")
    message(paste0("... ... we will export in the ", angle_of_matrix_for_intersectPeakMatrix," wide of intersectPeakMatrix output since the input angle_of_matrix_for_intersectPeakMatrix = '", angle_of_matrix_for_intersectPeakMatrix,"'"))
    for (i in 1:nrow(fun_output))
    {
      if (saving_id_x_for_intersectPeakMatrix != "all" && !(rownames(fun_output)[i] %in% saving_id_x_for_intersectPeakMatrix))
      {
        next
      }
      for (j in 1:ncol(fun_output))
      {
        if (saving_id_y_for_intersectPeakMatrix != "all" && !(colnames(fun_output)[j] %in% saving_id_y_for_intersectPeakMatrix))
        {
          next
        }
        fun_output_i <- fun_output[i, j][[1]]
        if (class(fun_output_i)[1] != "IntersectPeakMatrix")
        {
          message(paste0("Error: intersectPeakMatrix output is a matrix of IntersectPeakMatrix class objects. We found the element row No.", i,", column No.", j," is not IntersectPeakMatrix class in your input!"))
          message("... ... Please make sure your input is intersectPeakMatrix output and you have provided the correct information in 'fun ='!")
          message("EXIT!!")
          return(NULL)
        }
       if (angle_of_matrix_for_intersectPeakMatrix == "x")
       {
         id_i <- fun_output_i@id_x
         message(paste0("... ... ... export id = ", id_i))
         isTFregulomeID <- fun_output_i@isxTFregulomeID
         nsites <- fun_output_i@MethMotif_x@MMmotif@nsites
         if (isTFregulomeID == TRUE && nsites > 0)
         {
           # if save betaScore matrix
           if (save_betaScore_matrix)
           {
             beta_score_matrix <- fun_output_i@MethMotif_x@MMBetaScore
             if (is.na(beta_score_matrix[1,1]))
             {
               message("... ... ... ... No beta score matrix is available. Skip!")
             }
             else
             {
               saved_BetaScore <- save_BetaScore(MethMotif_object = fun_output_i@MethMotif_x)
               message(paste0("... ... ... ... Beta score matrix has been saved as '", saved_BetaScore,"'."))
             }
           }
           # if save motif PFM
           if (save_motif_PFM)
           {
             saved_motifPFM <- save_motifPFM(MethMotif_object = fun_output_i@MethMotif_x)
             message(paste0("... ... ... ... Motif PFM has been saved as '", saved_motifPFM,"'."))
           }
         }
         else
         {
           message("... ... ... the original peaks of ",id_i, " is not loaded from MethMotif database, or in the intersected peaks the number of TFBS recorded in MethMotif database is zero. Hence no further action for this id!")
         }
       }
       else
       {
         id_i <- fun_output_i@id_y
         message(paste0("... ... ... export id = ", id_i))
         isTFregulomeID <- fun_output_i@isyTFregulomeID
         nsites <- fun_output_i@MethMotif_y@MMmotif@nsites
         if (isTFregulomeID == TRUE && nsites > 0)
         {
           # if save betaScore matrix
           if (save_betaScore_matrix)
           {
             beta_score_matrix <- fun_output_i@MethMotif_y@MMBetaScore
             if (is.na(beta_score_matrix[1,1]))
             {
               message("... ... ... ... No beta score matrix is available. Skip!")
             }
             else
             {
               saved_BetaScore <- save_BetaScore(MethMotif_object = fun_output_i@MethMotif_y)
               message(paste0("... ... ... ... Beta score matrix has been saved as '", saved_BetaScore,"'."))
             }
           }
           # if save motif PFM
           if (save_motif_PFM)
           {
             saved_motifPFM <- save_motifPFM(MethMotif_object = fun_output_i@MethMotif_y)
             message(paste0("... ... ... ... Motif PFM has been saved as '", saved_motifPFM,"'."))
           }
         }
         else
         {
           message("... ... ... the original peaks of ",id_i, " is not loaded from MethMotif database, or in the intersected peaks the number of TFBS recorded in MethMotif database is zero. Hence no further action for this id!")
         }
        }
      }
    }
  }
}

save_BetaScore <- function(MethMotif_object)
{
  id <- MethMotif_object@MMmotif@id
  beta_score <- as.data.frame(MethMotif_object@MMBetaScore)
  colnames_original <- sapply(seq(1, ncol(beta_score), 1), as.character)
  colnames(beta_score) <- colnames_original
  beta_score$position <- c("beta score < 10%", "beta score in between", "beta score > 90%")
  beta_score <- rbind(beta_score, colnames(beta_score))
  beta_score <- beta_score[c(4,1:3), ]
  beta_score <- beta_score[, c("position", colnames_original)]
  write.table(beta_score, paste0(id, "-methScore.txt"), sep = "\t", col.names = F, row.names = F, quote = F)
  return(paste0(id, "-methScore.txt"))
}

save_motifPFM <- function(MethMotif_object)
{
  motif_format <- MethMotif_object@MMmotif@motif_format
  if (motif_format == "MEME")
  {
    line1 <- paste0("MEME version ", MethMotif_object@MMmotif@version)
    line2 <- ""
    line3 <- paste0("ALPHABET= ", MethMotif_object@MMmotif@alphabet)
    line4 <- ""
    line5 <- paste0("strands: ", MethMotif_object@MMmotif@strand)
    line6 <- ""
    line7 <- "Background letter frequencies"
    background <- MethMotif_object@MMmotif@background
    line8 <- paste0("A ", background[1], " C ", background[2], " G ", background[3], " T ", background[4], " ")
    line9 <- ""
    line10 <- paste0("MOTIF ", MethMotif_object@MMmotif@id, " ", MethMotif_object@MMmotif@alternate_name)
    line11 <- ""
    line12 <- paste0("letter-probability matrix: alength= 4 w= ",MethMotif_object@MMmotif@width, " nsites= ",
                     MethMotif_object@MMmotif@nsites, " E= ", MethMotif_object@MMmotif@evalue)
    motif_matrix <- MethMotif_object@MMmotif@motif_matrix
    line_matrix <- ""
    for (i in 1:nrow(motif_matrix))
    {
      if (i < nrow(motif_matrix)){
        line_matrix <- paste0(line_matrix, "  ", specify_decimal(as.numeric(motif_matrix[i,1]), 6), "\t  ", specify_decimal(as.numeric(motif_matrix[i,2]), 6),
                              "\t  ", specify_decimal(as.numeric(motif_matrix[i,3]), 6), "\t  ", specify_decimal(as.numeric(motif_matrix[i,4]), 6), "\t\n")
      }
      else
      {
        line_matrix <- paste0(line_matrix, "  ", specify_decimal(as.numeric(motif_matrix[i,1]), 6), "\t  ", specify_decimal(as.numeric(motif_matrix[i,2]), 6),
                              "\t  ", specify_decimal(as.numeric(motif_matrix[i,3]), 6), "\t  ", specify_decimal(as.numeric(motif_matrix[i,4]), 6), "\t")
      }
    }
    file_name <- paste0(MethMotif_object@MMmotif@id, "-motif-MEME.txt")
    content_all <- c(line1, line2, line3, line4, line5, line6, line7, line8, line9, line10, line11, line12, line_matrix)
    fileConn <- file(file_name)
    writeLines(content_all, fileConn)
    close(fileConn)
    return(file_name)
  }
  else
  {
    id <- MethMotif_object@MMmotif@id
    line1 <- paste0("AC ", id)
    line2 <- "XX"
    line3 <- paste0("ID ", id)
    line4 <- "XX"
    line5 <- paste0("DE ", id, " ", MethMotif_object@MMmotif@alternate_name, " ; from TFregulomeR")
    line6 <- paste0("PO\tA\tC\tG\tT")
    motif_matrix <- MethMotif_object@MMmotif@motif_matrix
    line_matrix <- ""
    for (i in 1:nrow(motif_matrix))
    {
      if(i < nrow(motif_matrix))
      {
        line_matrix <- paste0(line_matrix, i, "\t", motif_matrix[i,1], "\t", motif_matrix[i,2], "\t", motif_matrix[i,3], "\t", motif_matrix[i,4], "\n")
      }
      else
      {
        line_matrix <- paste0(line_matrix, i, "\t", motif_matrix[i,1], "\t", motif_matrix[i,2], "\t", motif_matrix[i,3], "\t", motif_matrix[i,4])
      }
    }
    line_aftermatrix_1 <- "XX"
    line_aftermatrix_2 <- "CC program: TFregulomeR"
    line_aftermatrix_3 <- "XX"
    line_aftermatrix_4 <- "//"
    file_name <- paste0(id, "-motif-TRANSFAC.txt")
    content_all <- c(line1, line2, line3, line4, line5, line6, line_matrix, line_aftermatrix_1, line_aftermatrix_2, line_aftermatrix_3, line_aftermatrix_4)
    fileConn <- file(file_name)
    writeLines(content_all, fileConn)
    close(fileConn)
    return(file_name)
  }
}

specify_decimal <- function(x, k) trimws(format(round(x, k), nsmall=k))


