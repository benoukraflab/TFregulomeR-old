#' convert motif PFM in MethMotif database into PFMatrix class object in TFBSTools package
#'
#' This function allows you to retrieve and convert motif PFM in MethMotif database into PFMatrix class object, which can be used in TFBSTools package.
#' @param id Required. MethMotif ID.
#' @param methmotif_url Optional. If the MethMoitf url is NO more "http://bioinfo-csi.nus.edu.sg/methmotif/", please use a new url.
#' @return  An object of class PFMatrix
#' @keywords toTFBSTools
#' @export
#' @examples
#' CEBPB_pfm <- toTFBSTools(id = "MM1_HSA_K562_CEBPB")

toTFBSTools <- function(id, methmotif_url)
{
  if (missing(id))
  {
    stop("Please provide a MethMotif id using 'id ='")
  }
  # make an appropriate API url
  if (missing(methmotif_url)){
    methmotif_url <- "http://bioinfo-csi.nus.edu.sg/methmotif/api/table_query/"
    # store methmotif_url as methmotif_url_bk for searchMethMotif() later
    methmotif_url_bk <- "http://bioinfo-csi.nus.edu.sg/methmotif"
  } else if (endsWith(methmotif_url, suffix = "/index.php")==TRUE){
    # store methmotif_url as methmotif_url_bk for readMMmotif() later
    methmotif_url_bk <- methmotif_url
    methmotif_url <- gsub("index.php", "", methmotif_url)
    methmotif_url <- paste0(methmotif_url, "api/table_query/")
  } else if (endsWith(methmotif_url, suffix = "/")==TRUE){
    # store methmotif_url as methmotif_url_bk for readMMmotif() later
    methmotif_url_bk <- methmotif_url
    methmotif_url <- paste0(methmotif_url, "api/table_query/")
  } else {
    # store methmotif_url as methmotif_url_bk for readMMmotif() later
    methmotif_url_bk <- methmotif_url
    methmotif_url <- paste0(methmotif_url, "/api/table_query/")
  }

  methmotif_output <- suppressMessages(searchMethMotif(id = id, motif_format = "transfac"))
  if (is.null(methmotif_output))
  {
    message(paste0("No record for id ", id, " in MethMotif database!"))
    return(NULL)
  }
  else
  {
    methmotif_output_transfac <- methmotif_output[[1]]
    pfm <- PFMatrix(ID = methmotif_output_transfac@MMmotif@id, name = methmotif_output_transfac@MMmotif@alternate_name, strand = "*",
                   bg = methmotif_output_transfac@MMmotif@background, profileMatrix = t(methmotif_output_transfac@MMmotif@motif_matrix))
    return(pfm)
  }
}
