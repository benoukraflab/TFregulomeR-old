APIconnection <- function(){
  test_id <- "http://bioinfo-csi.nus.edu.sg/methmotif/api/table_query/listTFBS.php?AllTable=F&id=MM1_HSA_K562_CEBPB"
  api_request <- fromJSON(test_id)
  api_request_res <- as.data.frame(api_request$TFBS_records)
  RUnit::checkTrue(nrow(api_request_res)>0)
}
