#' load peaks from TFregulomeR
#'
#' This function allows you to obtain the peaks from TFregulomeR using TFregulomeR ID.
#' @param id Required. TFregulomeR ID
#' @param includeMotifOnly Either TRUE or FALSE (default). If TRUE, only peaks with motif will be returned
#' @param server server localtion to be linked, either 'sg' or 'ca'.
#' @param TFregulome_url TFregulomeR server is implemented in MethMotif server. If the MethMotif url is NO more "https://bioinfo-csi.nus.edu.sg/methmotif/" or "https://methmotif.org", please use a new url.
#' @return  a data.frame containing peak coordinates
#' @keywords loadPeaks
#' @export
#' @examples
#' CEBPB_peaks <- loadPeaks(id = "MM1_HSA_K562_CEBPB")

loadPeaks <- function(id, includeMotifOnly = FALSE,
                      server = "ca", TFregulome_url)
{
  # check server location
  if (server != "sg" && server != "ca")
  {
    stop("server should be either 'sg' (default) or 'ca'!")
  }

  # make an appropriate API url
  if (missing(TFregulome_url)){
    if(server == 'sg')
    {
      TFregulome_url <- "https://bioinfo-csi.nus.edu.sg/methmotif/api/table_query/"
    }
    else
    {
      TFregulome_url <- "https://methmotif.org/api/table_query/"
    }
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
    stop("Please input a TFregulomeR id using 'id = '.")
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
      message("There is a warning to connect TFregulomeR API!")
      message("Advice:")
      message("1) Check internet access;")
      message("2) Check dependent package 'jsonlite';")
      message("3) Current TFregulomeR server is implemented in MethMotif database, whose homepage is 'https://bioinfo-csi.nus.edu.sg/methmotif/' or 'https://methmotif.org'. If MethMotif homepage url is no more valid, please Google 'MethMotif', and input the valid MethMotif homepage url using 'TFregulome_url = '.")
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
          read.delim(peak_file, sep = "\t", header = FALSE)
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
        colnames(peak_df) <- c("chr","start","end","id","tag_fold_change")
        # use new id for each peak region
        if (includeMotifOnly)
        {
          peak_df$id = paste0(id,"_peaks_with_motif_", as.vector(rownames(peak_df)))
        }
        else
        {
          peak_df$id = paste0(id,"_all_peaks_", as.vector(rownames(peak_df)))
        }
        message("Success: peak file has been returned in a data frame!")
        return(peak_df)
      }
    }
  }
}
