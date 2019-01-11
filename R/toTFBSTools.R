#' convert motif PFM in TFregulome database into PFMatrix class object in TFBSTools package
#'
#' This function allows you to retrieve and convert motif PFM in TFregulome database into PFMatrix class object, which can be used in TFBSTools package.
#' @param id Required. TFregulome ID.
#' @param TFregulome_url TFregulome server is implemented in MethMotif server. If the MethMoitf url is NO more "http://bioinfo-csi.nus.edu.sg/methmotif/", please use a new url.
#' @return  An object of class PFMatrix
#' @keywords toTFBSTools
#' @export
#' @examples
#' CEBPB_pfm <- toTFBSTools(id = "MM1_HSA_K562_CEBPB")

toTFBSTools <- function(id, TFregulome_url)
{
  if (missing(id))
  {
    stop("Please provide a TFregulome ID using 'id ='")
  }
  # make an appropriate API url
  if (missing(TFregulome_url)){
    TFregulome_url <- "http://localhost:8888/api/table_query/"
    # store TFregulome_url as TFregulome_url_bk for searchMotif() later
    TFregulome_url_bk <- "http://localhost:8888"
  } else if (endsWith(TFregulome_url, suffix = "/index.php")==TRUE){
    # store TFregulome_url as TFregulome_url_bk for searchMotif() later
    TFregulome_url_bk <- TFregulome_url
    TFregulome_url <- gsub("index.php", "", TFregulome_url)
    TFregulome_url <- paste0(TFregulome_url, "api/table_query/")
  } else if (endsWith(TFregulome_url, suffix = "/")==TRUE){
    # store TFregulome_url as TFregulome_url_bk for searchMotif() later
    TFregulome_url_bk <- TFregulome_url
    TFregulome_url <- paste0(TFregulome_url, "api/table_query/")
  } else {
    # store TFregulome_url as TFregulome_url_bk for searchMotif() later
    TFregulome_url_bk <- TFregulome_url
    TFregulome_url <- paste0(TFregulome_url, "/api/table_query/")
  }

  methmotif_output <- suppressMessages(searchMotif(id = id, motif_format = "TRANSFAC",
                                                   TFregulome_url=TFregulome_url_bk))
  if (is.null(methmotif_output))
  {
    message(paste0("No record for id ", id, " in TFregulome database!"))
    return(NULL)
  }
  else
  {
    methmotif_output_transfac <- methmotif_output
    pfm <- PFMatrix(ID = methmotif_output_transfac@MMmotif@id, name = methmotif_output_transfac@MMmotif@alternate_name, strand = "*",
                   bg = methmotif_output_transfac@MMmotif@background, profileMatrix = t(methmotif_output_transfac@MMmotif@motif_matrix))
    return(pfm)
  }
}
