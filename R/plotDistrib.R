#' plot TFBS distribution from motifDistrib() output
#'
#' This function allows you to plot TFBS distributions in a given list of peak sets from the output of motifDistrib().
#' @param motifDistrib Required. motifDistrib() output.
#' @return  TFBS distribution PDF file.
#' @keywords plotDistrib
#' @export
#' @examples
#' CEBPB_peaks <- loadPeaks(id = "MM1_HSA_K562_CEBPB")
#' motifDistrib_output <- motifDistrib(id = "MM1_HSA_K562_CEBPB", peak_list = list(CEBPB_peaks),
#'                                     peak_id = c("MM1_HSA_K562_CEBPB"))
#' plotDistrib(motifDistrib = motifDistrib_output)

plotDistrib <- function(motifDistrib)
{
  if (missing(motifDistrib))
  {
    stop("Please provided output from motifDistrib() function using 'motifDistrib ='!")
  }
  for (i in 1:length(motifDistrib))
  {
    target_id <- motifDistrib[[i]]$target_id
    peak_id <- names(motifDistrib)[i]
    input_peak_num <- motifDistrib[[i]]$input_peak_number
    peak_with_motif <- motifDistrib[[i]]$peak_with_motif
    motif_occurrence <- motifDistrib[[i]]$motif_occurrence
    pdf_name <- paste0("motif_", target_id, "_in_peak_set_", peak_id, ".pdf")
    plot_title <- paste0("motif ", target_id, " in peak set ", peak_id,"\n input peak number: ", input_peak_num, "\n peak with motif: ", peak_with_motif)
    window_size <- length(motif_occurrence)-1
    pdf(pdf_name)
    plot(x=seq((-1*window_size)/2, window_size/2, 1), y = motif_occurrence, type="l",
         xlab="Position relative to center (bp)", ylab="Motif percentage (%)",
         main=plot_title, cex.main=0.8)
    dev.off()
    message(paste0("Distribution of motif ", target_id," in peak set ", peak_id, " has been saved!"))
  }
}
