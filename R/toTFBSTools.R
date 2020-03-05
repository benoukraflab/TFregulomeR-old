#' convert motif PFM in TFregulomeR into PFMatrix class object in TFBSTools package
#'
#' This function allows you to retrieve and convert motif PFM in TFregulomeR database into PFMatrix class object, which can be used in TFBSTools package.
#' @param id Required. TFregulomeR ID.
#' @param server server localtion to be linked, either 'sg' or 'ca'.
#' @param TFregulome_url TFregulomeR server is implemented in MethMotif server. If the MethMotif url is NO more "http://bioinfo-csi.nus.edu.sg/methmotif/" or "http://methmotif.org", please use a new url.
#' @return  An object of class PFMatrix
#' @keywords toTFBSTools
#' @export
#' @examples
#' require(TFBSTools)
#' CEBPB_pfm <- toTFBSTools(id = "MM1_HSA_K562_CEBPB")

toTFBSTools <- function(id, server = "sg",TFregulome_url)
{
  if (missing(id))
  {
    stop("Please provide a TFregulomeR ID using 'id ='")
  }

  # check server location
  if (server != "sg" && server != "ca")
  {
    stop("server should be either 'sg' (default) or 'ca'!")
  }
  # make an appropriate API url
  if (missing(TFregulome_url)){
    if(server == 'sg')
    {
      TFregulome_url <- "http://bioinfo-csi.nus.edu.sg/methmotif/api/table_query/"
      # store TFregulome_url as TFregulome_url_bk for searchMotif() later
      TFregulome_url_bk <- "http://bioinfo-csi.nus.edu.sg/methmotif"
    }
    else
    {
      TFregulome_url <- "http://methmotif.org/api/table_query/"
      # store TFregulome_url as TFregulome_url_bk for searchMotif() later
      TFregulome_url_bk <- "http://methmotif.org"
    }
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
    message(paste0("No record for id ", id, " in TFregulomeR!"))
    return(NULL)
  }
  else
  {
    methmotif_output_transfac <- methmotif_output
    pfm <- TFBSTools::PFMatrix(ID = methmotif_output_transfac@MMmotif@id, name = methmotif_output_transfac@MMmotif@alternate_name, strand = "*",
                   bg = methmotif_output_transfac@MMmotif@background, profileMatrix = t(methmotif_output_transfac@MMmotif@motif_matrix))
    return(pfm)
  }
}
