#' Search motif PFM and beta score matrix (if source is MethMotif) for a given TFregulomeR ID in TFregulomeR
#'
#' This function allows you to obtain motif PFM matrix and beta score matrix (if source is MethMotif) for a given TFregulomeR ID in TFregulomeR
#' @param id Required. TFregulomeR ID.
#' @param motif_format Motif PFM format, either in MEME by default or TRANSFAC.
#' @param server server localtion to be linked, either 'sg' or 'ca'.
#' @param TFregulome_url TFregulomeR server is implemented in MethMotif server. If the MethMotif url is NO more "http://bioinfo-csi.nus.edu.sg/methmotif/" or "http://methmotif.org", please use a new url.
#' @return MethMotif class object
#' @keywords MethMotif
#' @export
#' @examples
#' K562_CEBPB <- searchMotif(id = "MM1_HSA_K562_CEBPB")
#' K562_CEBPB_transfac <- searchMotif(id = "MM1_HSA_K562_CEBPB",
#'                                    motif_format = "TRANSFAC")

searchMotif <- function(id, motif_format = "MEME",
                        server = "sg", TFregulome_url)
{
  # check motif_format MEME and TRANSFAC.
  motif_format = toupper(motif_format)
  if (motif_format != "MEME" & motif_format != "TRANSFAC")
  {
    stop("Please check motif_format! Currently we only support MEME (default) and TRANSFAC formats!")
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
    }
    else
    {
      TFregulome_url <- "http://methmotif.org/api/table_query/"
    }
  } else if (endsWith(TFregulome_url, suffix = "/index.php")==TRUE){
    TFregulome_url <- gsub("index.php", "", TFregulome_url)
    TFregulome_url <- paste0(TFregulome_url, "api/table_query/")
  } else if (endsWith(TFregulome_url, suffix = "/")==TRUE){
    TFregulome_url <- paste0(TFregulome_url, "api/table_query/")
  } else {
    TFregulome_url <- paste0(TFregulome_url, "/api/table_query/")
  }

  if(missing(id))
  {
    stop("Please input regulome ID using 'id = '!")
  }
  else
  {
    #query url
    query_url <- paste0("listTFBS.php?AllTable=F&id=", id)
  }

  # start to query
  request_content_json <- tryCatch({
    fromJSON(paste0(TFregulome_url,query_url))
  },
  error = function(cond)
  {
    message("There is a warning to connect TFregulomeR API!")
    message("Advice:")
    message("1) Check internet access;")
    message("2) Check dependent package 'jsonlite';")
    message("3) Current TFregulomeR server is implemented in MethMotif database, whose homepage is 'http://bioinfo-csi.nus.edu.sg/methmotif/' or 'http://methmotif.org'. If MethMotif homepage url is no more valid, please Google 'MethMotif', and input the valid MethMotif homepage url using 'TFregulome_url = '.")
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
      message("There are a matched record exported in a MethMotif object.")
      # using the previous backup TFregulome_url to parse into readMMmotif()
      if (motif_format == "MEME")
      {
        motif_file_path <- request_content_df[1, "motif_MEME"]
      }
      else
      {
        motif_file_path <- request_content_df[1, "motif_TRANSFAC"]
      }
      betaScore_file_path <- request_content_df[1, "beta_score_matrix"]
      num_peak <- request_content_df[1,"peak_with_motif_num"]
      nsites <- request_content_df[1,"TFBS_num"]

      # MEME file path for background extraction when motif formati is TRANSFAC
      motif_file_path_MEME <- request_content_df[1, "motif_MEME"]
      methmotif_item_motif <- readMMmotif(motif_file_path = motif_file_path,
                                          motif_format = motif_format,
                                          id = id,
                                          num_peak = as.integer(num_peak),
                                          nsites = as.integer(nsites),
                                          motif_file_path_MEME = motif_file_path_MEME)
      methmotif_item_betaScore <- readBetaScoreMatrix(betaScore_file_path = betaScore_file_path)
      methmotif_item <- new("MethMotif")
      methmotif_item <- updateMethMotif(methmotif_item,
                                        MMBetaScore = methmotif_item_betaScore,
                                        MMmotif = methmotif_item_motif)

      return(methmotif_item)
    }
  }
  else
  {
    message("Empty output for the request!")
    return(NULL)
  }
}
