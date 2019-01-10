#' load peaks from MethMotif database
#'
#' This function allows you to obtain the peaks using TFregulome ID.
#' @param id Required. TFregulome ID
#' @param includeMotifOnly Either TRUE or FALSE (default). If TRUE, only peaks with motif will be returned
#' @param TFregulome_url TFregulome server is implemented in MethMotif server. If the MethMoitf url is NO more "http://bioinfo-csi.nus.edu.sg/methmotif/", please use a new url.
#' @return  a data.frame containing peak coordinates
#' @keywords loadPeaks
#' @export
#' @examples
#' CEBPB_peaks <- loadPeaks(id = "MM1_HSA_K562_CEBPB")

loadPeaks <- function(id, includeMotifOnly = F, TFregulome_url)
{
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

  # check input argument id
  if (missing(id))
  {
    stop("Please input a TFregulome id using 'id = '.")
  }
  else
  {
    query_url <- paste0("listTFBS.php?AllTable=F&id=",id)
    #parse JSON from API endpoint
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
    if (!is.null(request_content_json))
    {
      request_content_df <- as.data.frame(request_content_json$TFBS_records)
      if (nrow(request_content_df)==0)
      {
        message(request_content_json$message)
        return(NULL)
      }
      else
      {
        if (includeMotifOnly)
        {
          peak_file <- request_content_df[1,c("peak_with_motif_file")]
        }
        else
        {
          peak_file <- request_content_df[1,c("all_peak_file")]
        }
        # read peak file
        peak_df <- tryCatch({
          read.delim(peak_file, sep = "\t", header = F)
        },
        warning = function(w)
        {
          message(paste0("Warning: No peak file for id=", id))
          stop()
        },
        error = function(err)
        {
          message(paste0("Erro: No peak file for id=", id))
        })
        peak_df = peak_df[,c(1,2,3)]
        # use new id for each peak region
        if (includeMotifOnly)
        {
          peak_df$V4 = paste0(id,"_peaks_with_motif_", as.vector(rownames(peak_df)))
        }
        else
        {
          peak_df$V4 = paste0(id,"_all_peaks_", as.vector(rownames(peak_df)))
        }
        colnames(peak_df) = c("chr","start","end","id")
        message("Success: peak file has been returned in a data frame!")
        return(peak_df)
      }
    }
  }
}
