#' browse existing TFBSs available in TFregulome server
#'
#' This function allows you to get existing TFBSs in TFregulome database
#' @param species The species of interset
#' @param organ The organ of interset
#' @param sample_type The sample type of interset
#' @param cell_tissue_name The name of tissue or cell of interset
#' @param tf The TF of interset
#' @param disease_state The disease state of interset
#' @param source The source of interset
#' @param TFregulome_url TFregulome server is implemented in MethMotif server. If the MethMoitf url is NO more "http://bioinfo-csi.nus.edu.sg/methmotif/", please use a new url.
#' @return  data.frame giving the existing TFBSs in MethMotif database for a given cell type
#' @keywords TFBS
#' @export
#' @examples
#' TFBS_brain <- TFBSBrowser(organ = "brain")

TFBSBrowser <- function(species, organ, sample_type, cell_tissue_name,
                        tf, disease_state, source, TFregulome_url)
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

  query_index <- rep(0,7)
  query_value <- rep("",7)
  names(query_index) <- c("species","organ","sample_type","cell_or_tissue_name",
                          "tf","disease_state","source")
  names(query_value) <- c("species","organ","sample_type","cell_or_tissue_name",
                          "tf","disease_state","source")
  if (!missing(species))
  {
    query_index["species"] <- 1
    query_value["species"] <- paste0("species=",species)
  }
  if (!missing(organ))
  {
    query_index["organ"] <- 1
    query_value["organ"] <- paste0("organ=",organ)
  }
  if (!missing(sample_type))
  {
    query_index["sample_type"] <- 1
    query_value["sample_type"] <- paste0("sample_type=",sample_type)
  }
  if (!missing(cell_tissue_name))
  {
    query_index["cell_tissue_name"] <- 1
    query_value["cell_tissue_name"] <- paste0("cell_tissue_name=",
                                               cell_tissue_name)
  }
  if (!missing(tf))
  {
    query_index["tf"] <- 1
    query_value["tf"] <- paste0("TF=",tf)
  }
  if (!missing(disease_state))
  {
    query_index["disease_state"] <- 1
    query_value["disease_state"] <- paste0("disease_state=",disease_state)
  }
  if (!missing(source))
  {
    query_index["source"] <- 1
    query_value["source"] <- paste0("source=",source)
  }

  if (sum(query_index)==0)
  {
    query_url <- "listTFBS.php?AllTable=T"
  }
  else
  {
    query_url <- paste0("listTFBS.php?AllTable=F&",
                        paste0(query_value[query_index==1], collapse = "&"))
  }
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
      request_content_df_output <- request_content_df[,c("ID", "species", "organ",
                                                        "sample_type",
                                                        "cell_tissue_name",
                                                        "description",
                                                        "disease_state",
                                                        "TF","source","source_ID",
                                                        "peak_num","peak_with_motif_num")]
      tfbs_num <- nrow(request_content_df_output)
      species_result <- unique(request_content_df_output$species)
      species_num <- length(species_result)
      organ_result <- unique(request_content_df_output$organ)
      organ_num <- length(organ_result)
      sample_type_result <- unique(request_content_df_output$sample_type)
      sample_type_num <- length(sample_type_result)
      cell_tissue_name_num <- length(unique(request_content_df_output$cell_tissue_name))
      disease_state_result <- unique(request_content_df_output$disease_state)
      disease_state_num <- length(disease_state_result)
      tf_num <- length(unique(request_content_df_output$TF))
      source_result <- unique(request_content_df_output$source)
      message(paste0(tfbs_num," TFBS(s) founded: ..."))
      message(paste0("... covering ", tf_num, " TF(s)"))
      message(paste0("... from ", species_num," species:"))
      message(paste0("... ...", paste0(species_result, collapse = ", ")))
      message(paste0("... from ", organ_num," organ(s):"))
      message(paste0("... ... ", paste0(organ_result, collapse = ", ")))
      message(paste0("... in ", sample_type_num, " sample type(s):"))
      message(paste0("... ... ", paste0(sample_type_result, collapse = ", ")))
      message(paste("... in ", cell_tissue_name_num, " different cell(s) or tissue(s)"))
      message(paste0("... in ", disease_state_num," type(s) of disease state(s):"))
      message(paste0("... ... ", paste0(disease_state_result, collapse = ", ")))
      message(paste0("... from the source(s): ", paste0(source_result, collapse = ", ")))
      return(request_content_df_output)
    }
  }
  else
  {
    message("Empty output for the request!")
    return(NULL)
  }
}
