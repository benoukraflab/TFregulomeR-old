#' greatAnnotate
#'
#' This function allows you to analyse the gene ontologies of targeting genes by cis-regulatory regions.
#' @param peaks Required. A bed-format genomic regions in data frame.
#' @param assembly Optional. The assembly the input regions are using, currently supporting 'hg19', 'hg38' (default), 'mm9' and 'mm10'.
#' @param return_annotation Optional. Either TRUE of FALSE (default). If TRUE, a data.frame containing annotation results will be returned.
#' @param return_html_report Optional. Either TRUE of FALSE (default). If TRUE, a dynamic HTML report will be saved.
#' @param pvalue Optional. The adjusted p-value which is applied to filter the results, by default 0.01.
#' @param test Optional. The statistical test used in GREAT analysis, either 'binomial' (default) or 'hypergeometric'.
#' @param great_rule Optional. Equivalent to the rGREAT input 'rule', 'basalPlusExt' (default, basal plus extension), 'twoClosest' (two nearest genes), or 'oneClosest' (single nearest gene).
#' @param great_adv_upstream Optional. Equivalent to the rGREAT input 'adv_upstream' (upstream extension, kb). Only applicable when 'great_rule' is  'basalPlusExt', by default 5.0 (kb).
#' @param great_adv_downstream Optional. Equivalent to the rGREAT input 'adv_downstream' (downstream extension, kb). Only applicable when 'great_rule' is  'basalPlusExt', by default 1.0 (kb).
#' @param great_adv_span Optional. Equivalent to the rGREAT input 'adv_span' (maximal distal region). Only applicable when 'great_rule' is  'basalPlusExt', by default 1000.0 (kb).
#' @param great_adv_twoDistance Optional. Equivalent to the rGREAT input 'adv_twoDistance' (region range to be considered). Only applicable when 'great_rule' is  'twoClosest', by default 1000.0 (kb).
#' @param great_adv_oneDistance Optional. Equivalent to the rGREAT input 'adv_oneDistance' (region range to be considered). Only applicable when 'great_rule' is  'oneClosest', by default 1000.0 (kb).
#' @param request_interval Optional. The minimal gap time between two requests using greatAnnotate, by default 60 (s).
#' @param great_version Optional. Equivalent to the rGREAT input 'version', by default 3.0.
#' @return  a data.frame, or an HTML report depending on the options.
#' @keywords greatAnnotate
#' @export
#' @examples
#' K562_CEBPB_regions <- loadPeaks(id = "MM1_HSA_K562_CEBPB")
#' K562_CEBPB_regions_annotation <- greatAnnotate(peaks = K562_CEBPB_regions,
#' return_annotation = T, return_html_report = T)

greatAnnotate <- function(peaks, assembly = "hg38", return_annotation = FALSE,
                          return_html_report = FALSE, pvalue = 0.01, test = "binomial",
                          great_rule = "basalPlusExt", great_adv_upstream = 5.0,
                          great_adv_downstream = 1.0, great_adv_span = 1000.0,
                          great_adv_twoDistance = 1000.0, great_adv_oneDistance = 1000.0,
                          request_interval = 60, great_version = 3.0)
{
  # check input arguments
  if (missing(peaks))
  {
    stop("please provide peak regions using 'peaks ='!")
  }
  if (class(peaks) != "data.frame")
  {
    stop("The 'peaks' should be a BED-format data.frame!")
  }
  if (!(assembly %in% c("hg38","hg19","mm10","mm9")))
  {
    stop("Currently greatAnnotate only supports hg19, hg38, mm9 and mm10.")
  }
  if (!is.logical(return_annotation))
  {
    stop("'return_annotation' should be either TRUE (T) or FALSE (F, default)")
  }
  if (!is.logical(return_html_report))
  {
    stop("'return_html_report' should be either TRUE (T) or FALSE (F, default)")
  }
  if (!is.numeric(pvalue))
  {
    stop("'pvalue' should be a numeric!")
  }
  if (!(test %in% c("binomial", "hypergeometric")))
  {
    stop("'test' should be either 'binomial' (default) or 'hypergeometric'!")
  }
  if (!(great_rule %in% c("basalPlusExt", "twoClosest", "oneClosest")))
  {
    stop("According to 'rGREAT', 'great_rule' shoule be 'basalPlusExt' (default), 'twoClosest', or 'oneClosest'!")
  }
  # check loaded package
  if (!("rGREAT" %in% (.packages())))
  {
    stop("GREAT R package 'rGREAT' (>=1.14.0) is NOT loaded yet!")
  }
  if (assembly == "hg38" && !("liftOver" %in% (.packages())))
  {
    stop("Your input is hg38. Currently GREAT doesn't support hg38. We need to convert hg38 to hg19 using liftOver (>=1.4.0). But liftOver package is NOT loaded yet!")
  }
  if (return_html_report==TRUE && !("rbokeh" %in% (.packages())))
  {
    stop("A dynamic html report requires 'rbokeh' pakcage (>=0.5.0). 'rbokeh' is NOT loaded yet!")
  }

  #message
  message("Start greatAnnotate ...")
  if (return_annotation)
  {
    message("... ... You chose to return annotated results in a data.frame.")
  }
  else
  {
    message("... ... You chose NOT to return annotated results in a data.frame.")
  }
  if (return_html_report)
  {
    message("... ... You chose to return an HTML report.")
  }
  else
  {
    message("... ... You chose NOT to return an HTML report.")
  }
  if (!return_annotation && ! return_html_report)
  {
    message("... ... You chose no action. EXIT!!")
    return(NULL)
  }

  peaks <- peaks[,1:3]
  colnames(peaks) <- c("chr","start", "end")
  peaks$id <- paste0("greatAnnotate_peak_", as.vector(rownames(peaks)))
  # hg38 to hg19
  if (assembly == "hg38")
  {
    message("... ... assembly is hg38. Now converting to hg19 using liftOver...")
    message(paste0("... ... number of the original input regions is ", nrow(peaks)))
    peaks_grange <- with(peaks, GRanges(chr, IRanges(start, end), id=id))
    peaks <- hg38Tohg19(peaks_grange)
    assembly <- "hg19"
    message(paste0("... ... number of the regions successfully converted to hg19 is ", nrow(peaks)))
  }
  # great analysis
  message("... ... start GREAT analysis")
  if (great_rule == "basalPlusExt")
  {
    jobs <- submitGreatJob(gr = peaks, species = assembly, rule = great_rule,
                          adv_upstream = great_adv_upstream, adv_downstream = great_adv_downstream,
                          adv_span = great_adv_span, request_interval = request_interval, version = great_version)
  }
  else if (great_rule == "twoClosest")
  {
    jobs <- submitGreatJob(gr = peaks, species = assembly, rule = great_rule,
                          adv_twoDistance = great_adv_twoDistance, request_interval = request_interval,
                          version = great_version)
  }
  else
  {
    jobs <- submitGreatJob(gr = peaks, species = assembly, rule = great_rule,
                          adv_oneDistance = great_adv_oneDistance, request_interval = request_interval,
                          version = great_version)
  }
  tb <- getEnrichmentTables(jobs)
  go_MF <- tb$`GO Molecular Function`
  go_BP <- tb$`GO Biological Process`
  go_CC <- tb$`GO Cellular Component`
  if (test == "binomial")
  {
    if (nrow(go_MF) > 0)
    {
      go_MF_sorted <- go_MF[order(go_MF$Binom_Adjp_BH),]
      go_MF_sorted_pass <- go_MF_sorted[which(go_MF_sorted$Binom_Adjp_BH<pvalue),]
    }
    else
    {
      go_MF_sorted_pass <- go_MF
    }
    if (nrow(go_BP) > 0)
    {
      go_BP_sorted <- go_BP[order(go_BP$Binom_Adjp_BH),]
      go_BP_sorted_pass <- go_BP_sorted[which(go_BP_sorted$Binom_Adjp_BH<pvalue),]
    }
    else
    {
      go_BP_sorted_pass <- go_BP
    }
    if (nrow(go_CC) > 0)
    {
      go_CC_sorted <- go_CC[order(go_CC$Binom_Adjp_BH),]
      go_CC_sorted_pass <- go_CC_sorted[which(go_CC_sorted$Binom_Adjp_BH<pvalue),]
    }
    else
    {
      go_CC_sorted_pass <- go_CC
    }
  }
  else
  {
    if (nrow(go_MF) > 0)
    {
      go_MF_sorted <- go_MF[order(go_MF$Hyper_Adjp_BH),]
      go_MF_sorted_pass <- go_MF_sorted[which(go_MF_sorted$Hyper_Adjp_BH<pvalue),]
    }
    else
    {
      go_MF_sorted_pass <- go_MF
    }
    if (nrow(go_BP) > 0)
    {
      go_BP_sorted <- go_BP[order(go_BP$Hyper_Adjp_BH),]
      go_BP_sorted_pass <- go_BP_sorted[which(go_BP_sorted$Hyper_Adjp_BH<pvalue),]
    }
    else
    {
      go_BP_sorted_pass <- go_BP
    }
    if (nrow(go_CC) > 0)
    {
      go_CC_sorted <- go_CC[order(go_CC$Hyper_Adjp_BH),]
      go_CC_sorted_pass <- go_CC_sorted[which(go_CC_sorted$Hyper_Adjp_BH<pvalue),]
    }
    else
    {
      go_CC_sorted_pass <- go_CC
    }
  }
  # no result for MF?
  if (nrow(go_MF_sorted_pass) > 0)
  {
    go_MF_sorted_pass$go_id <- "MF"
  }
  else
  {
    go_MF_sorted_pass$go_id <- character(0)
  }
  # no result for BP?
  if (nrow(go_BP_sorted_pass) > 0)
  {
    go_BP_sorted_pass$go_id <- "BP"
  }
  else
  {
    go_BP_sorted_pass$go_id <- character(0)
  }
  # no result for CC?
  if (nrow(go_CC_sorted_pass) > 0)
  {
    go_CC_sorted_pass$go_id <- "CC"
  }
  else
  {
    go_CC_sorted_pass$go_id <- character(0)
  }
  all <- rbind(go_MF_sorted_pass, go_BP_sorted_pass, go_CC_sorted_pass)
  if (nrow(all) == 0){
    message("... ... No result for the current settings!")
    return(NULL)
  }
  else{
    if (return_html_report)
    {
      if (test == "binomial")
      {
        all$log10p <- -log10(all$Binom_Adjp_BH)
        all$Term <- all$name
        all$gene_number <- all$Hyper_Observed_Gene_Hits
        all$adjusted_pvalue <- all$Binom_Adjp_BH
      }
      else
      {
        all$log10p <- -log10(all$Hyper_Adjp_BH)
        all$Term <- all$name
        all$gene_number <- all$Hyper_Observed_Gene_Hits
        all$adjusted_pvalue <- all$Hyper_Adjp_BH
      }
      html_output <- formHTMLoutput(all, test)
      write(html_output, "greatAnnotate_result.html")
      message("... ... An html report has been generated as 'greatAnnotate_result.html'!")
    }
    if (return_annotation)
    {
      if (test == "binomial")
      {
        all_output <- all[,c("go_id","ID","name","Hyper_Observed_Gene_Hits", "Binom_Adjp_BH")]
      }
      else
      {
        all_output <- all[,c("go_id", "ID","name","Hyper_Observed_Gene_Hits", "Hyper_Adjp_BH")]
      }
      colnames(all_output) <- c("category","ID", "name", "number_of_targeting_genes", "adjusted_pvalue")
      message("... ... The annotation results have been returned in a data.frame!")
      return(all_output)
    }
  }
}

hg38Tohg19 <- function(peaks)
{
  seqlevelsStyle(peaks) <- "UCSC"
  path <- system.file(package="liftOver", "extdata", "hg38ToHg19.over.chain")
  ch <- import.chain(path)
  cur19 <- suppressWarnings(liftOver(peaks, ch))
  cur19 <- unlist(cur19)
  genome(cur19) <- "hg19"
  cur19_df <- data.frame(cur19)
  return(cur19_df[,c("seqnames","start","end","id")])
}

formHTMLoutput <- function(all, test)
{
  html_contents_CC <- rbokehPlot(inputdata = all[which(all$go_id=="CC"),], inputtype = "CC")
  modelid_CC <- html_contents_CC[[1]]
  docid_CC <- html_contents_CC[[2]]
  docs_json_CC <- html_contents_CC[[3]]
  html_contents_BP <- rbokehPlot(inputdata = all[which(all$go_id=="BP"),], inputtype = "BP")
  modelid_BP <- html_contents_BP[[1]]
  docid_BP <- html_contents_BP[[2]]
  docs_json_BP <- html_contents_BP[[3]]
  html_contents_MF <- rbokehPlot(inputdata = all[which(all$go_id=="MF"),], inputtype = "MF")
  modelid_MF <- html_contents_MF[[1]]
  docid_MF <- html_contents_MF[[2]]
  docs_json_MF <- html_contents_MF[[3]]

  go_BP_sorted_pass <- all[which(all$go_id=="BP"),]
  go_MF_sorted_pass <- all[which(all$go_id=="MF"),]
  go_CC_sorted_pass <- all[which(all$go_id=="CC"),]
  table_html <- ""
  table_html <- paste0(table_html, "
        <tbody id='table_BP' style='display:;'>")
  if (nrow(go_BP_sorted_pass)>0){
    for (elm in 1:nrow(go_BP_sorted_pass)){
      go_id = go_BP_sorted_pass[elm, "ID"]
      go_name = go_BP_sorted_pass[elm, "name"]
      if (test == "binomial")
      {
        adj_pvalue = go_BP_sorted_pass[elm, "Binom_Adjp_BH"]
      }
      else
      {
        adj_pvalue = go_BP_sorted_pass[elm, "Hyper_Adjp_BH"]
      }
      number_of_gene = go_BP_sorted_pass[elm, "Hyper_Observed_Gene_Hits"]
      table_html <- paste0(table_html, "
          <tr>
            <td>", go_id,"</td>
            <td>", go_name,"</td>
            <td>", adj_pvalue,"</td>
            <td>", number_of_gene,"</td>
          </tr>")
    }
    } else{
      table_html <- paste0(table_html, "
          <tr>
            <td>No result</td>
            <td>No result</td>
            <td>No result</td>
            <td>No result</td>
          </tr>")
    }
  table_html <- paste0(table_html, "
        </tbody>
        <tbody id='table_MF' style='display: none;'>")
  if (nrow(go_MF_sorted_pass)>0){
    for (elm in 1:nrow(go_MF_sorted_pass)){
      go_id = go_MF_sorted_pass[elm, "ID"]
      go_name = go_MF_sorted_pass[elm, "name"]
      if (test == "binomial")
      {
        adj_pvalue = go_MF_sorted_pass[elm, "Binom_Adjp_BH"]
      }
      else
      {
        adj_pvalue = go_MF_sorted_pass[elm, "Hyper_Adjp_BH"]
      }
      number_of_gene = go_MF_sorted_pass[elm, "Hyper_Observed_Gene_Hits"]
      table_html <- paste0(table_html, "
          <tr>
            <td>", go_id,"</td>
            <td>", go_name,"</td>
            <td>", adj_pvalue,"</td>
            <td>", number_of_gene,"</td>
          </tr>")
    }
    } else{
      table_html <- paste0(table_html, "
          <tr>
            <td>No result</td>
            <td>No result</td>
            <td>No result</td>
            <td>No result</td>
          </tr>")
    }
  table_html <- paste0(table_html, "
        </tbody>
        <tbody id='table_CC' style='display: none;'>")
  if (nrow(go_CC_sorted_pass)>0){
    for (elm in 1:nrow(go_CC_sorted_pass)){
      go_id = go_CC_sorted_pass[elm, "ID"]
      go_name = go_CC_sorted_pass[elm, "name"]
      if (test == "binomial")
      {
        adj_pvalue = go_CC_sorted_pass[elm, "Binom_Adjp_BH"]
      }
      else
      {
        adj_pvalue = go_CC_sorted_pass[elm, "Hyper_Adjp_BH"]
      }
      number_of_gene = go_CC_sorted_pass[elm, "Hyper_Observed_Gene_Hits"]
      table_html <- paste0(table_html, "
          <tr>
            <td>", go_id,"</td>
            <td>", go_name,"</td>
            <td>", adj_pvalue,"</td>
            <td>", number_of_gene,"</td>
          </tr>")
    }
    } else{
      table_html <- paste0(table_html, "
          <tr>
            <td>No result</td>
            <td>No result</td>
            <td>No result</td>
            <td>No result</td>
          </tr>")
    }
  table_html <- paste0(table_html, "
        </tbody>")

  output_html <- paste0("<!DOCTYPE html>
<html>
<head>
  <script src='https://code.jquery.com/jquery-3.2.1.min.js' ></script>
  <script src='https://cdn.pydata.org/bokeh/release/bokeh-0.12.2.min.js'></script>
  <link href='https://cdn.pydata.org/bokeh/release/bokeh-0.12.2.min.css' rel='stylesheet'>
  <link rel='stylesheet' href='https://maxcdn.bootstrapcdn.com/bootstrap/3.3.7/css/bootstrap.min.css'>
  <style>
    /* width */
    ::-webkit-scrollbar {
      width: 10px;
    }

    /* Track */
    ::-webkit-scrollbar-track {
      background: #f1f1f1;
    }

    /* Handle */
    ::-webkit-scrollbar-thumb {
      background: #888;
    }

    /* Handle on hover */
    ::-webkit-scrollbar-thumb:hover {
      background: #555;
    }
  </style>
</head>
<body>
<div align='center'>
  <h4>TFregulomeR - GreatAnnotate Result</h4>
</div>
<hr align='center' style='width: 80%'>
<br>
<div align='center'>
  Figure - Gene ontology analyses of targeting genes:
  <select id='change_go_term' onchange='changeTerm(this.value)'>
    <option value='BP' selected>Biological Process</option>
    <option value='MF'>Molecular Fcuntion</option>
    <option value='CC'>Cellular Component</option>
  </select>
</div>
<div id='great_result' align='center'>
  <div class='bk-root' class='plotdiv' style=' width: 500px; height: 500px;' align='left'>
    <div id='great_analysis' class='plotdiv' ></div>
  </div>
</div>
<div class='container' align='center'>
  <br>
  Table - Gene ontology analyses of targeting genes: <select id='change_go_term' onchange='changeTermTable(this.value)'>
    <option value='BP' selected>Biological Process</option>
    <option value='MF'>Molecular Fcuntion</option>
    <option value='CC'>Cellular Component</option>
  </select>
  <button type='button' class='btn btn-primary btn-sm' onclick=\"exportTableToCSV('greatAnnotate_results.csv')\">Export All To CSV File</button>
  <div class='row' style='width:90%;height:400px;overflow:auto;'>
    <div class='table-responsive'>
      <table class='table table-bordred table-striped'>
        <thead>
          <th>GO ID</th>
          <th>GO name</th>
          <th>adjusted p-value</th>
          <th>number of targeting genes</th>
        </thead>", table_html,"
      </table>
    </div>
  </div>
</div>
<br>
<hr align='center' style='width: 80%'>
<script type='text/javascript'>
Bokeh.$(function() {
  ", modelid_BP,"
  var elementid = 'great_analysis';
  ", docid_BP, "
  ", docs_json_BP, "
  var refkey = Object.keys(docs_json)[0]
  var refs = docs_json[refkey].roots.references
  function traverseObject(obj) {
    for (var key in obj) {
      if (obj[key].constructor === Object) {
        traverseObject(obj[key]);
      } else if (obj[key].constructor === Array) {
        for (var i = 0; i < obj[key].length; i++) {
          if (obj[key][i] === null)
          obj[key][i] = NaN;
        };
      }
    };
  }
  for (var i = 0; i < refs.length; i++) {
    if (refs[i].type === 'ColumnDataSource')
      traverseObject(refs[i].attributes.data);
  };
  var render_items = [{
    'docid': docid,
    'elementid': elementid,
    'modelid': modelid
  }];
  Bokeh.set_log_level('info');
  Bokeh.embed.embed_items(docs_json, render_items);
});
function changeTerm(opt){
  $('#great_result').html(\"<div class='bk-root' class='plotdiv' style=' width: 500px; height: 500px;' align='left'><div id='great_analysis' class='plotdiv'></div>\");
  var elementid = 'great_analysis';
  if (opt == 'CC'){
    ", modelid_CC, "
    ", docid_CC, "
    ", docs_json_CC, "
  }
  if (opt == 'MF'){
    ", modelid_MF, "
    ", docid_MF, "
    ", docs_json_MF,"
  }
  if (opt=='BP'){
    ", modelid_BP, "
    ", docid_BP, "
    ", docs_json_BP, "
  }
  var refkey = Object.keys(docs_json)[0]
  var refs = docs_json[refkey].roots.references
  function traverseObject(obj) {
    for (var key in obj) {
      if (obj[key].constructor === Object) {
        traverseObject(obj[key]);
      } else if (obj[key].constructor === Array) {
        for (var i = 0; i < obj[key].length; i++) {
          if (obj[key][i] === null)
            obj[key][i] = NaN;
        };
      }
    };
  }
  for (var i = 0; i < refs.length; i++) {
    if (refs[i].type === 'ColumnDataSource')
      traverseObject(refs[i].attributes.data);
  };
  var render_items = [{
    'docid': docid,
    'elementid': elementid,
    'modelid': modelid
  }];
  Bokeh.set_log_level('info');
  Bokeh.embed.embed_items(docs_json, render_items);
}

function changeTermTable(opt){
  if (opt == 'BP'){
    $('#table_BP').css('display', '');
    $('#table_MF').css('display', 'none');
    $('#table_CC').css('display', 'none');
  }
  if (opt == 'MF'){
    $('#table_BP').css('display', 'none');
    $('#table_MF').css('display', '');
    $('#table_CC').css('display', 'none');
  }
  if (opt == 'CC'){
    $('#table_BP').css('display', 'none');
    $('#table_MF').css('display', 'none');
    $('#table_CC').css('display', '');
  }
}
function downloadCSV(csv, filename) {
    var csvFile;
    var downloadLink;
    // CSV file
    csvFile = new Blob([csv], {type: 'text/csv'});
    // Download link
    downloadLink = document.createElement('a');
    // File name
    downloadLink.download = filename;
    // Create a link to the file
    downloadLink.href = window.URL.createObjectURL(csvFile);
    // Hide download link
    downloadLink.style.display = 'none';
    // Add the link to DOM
    document.body.appendChild(downloadLink);
    // Click download link
    downloadLink.click();
}
function exportTableToCSV(filename) {
    var csv = [];
    var rows = document.querySelectorAll('table tr');
    for (var i = 0; i < rows.length; i++) {
      var row = [], cols = rows[i].querySelectorAll('td, th');
      for (var j = 0; j < cols.length; j++)
      row.push(cols[j].innerText);
      csv.push(row.join(','));
    }
    // Download CSV file
    downloadCSV(csv.join('\\n'), filename);
}
</script>
</body>
</html>")
  return(output_html)
}

rbokehPlot <- function(inputdata, inputtype)
{
  if (inputtype == "BP")
  {
    col = "blue"
  }
  else if (inputtype == "MF")
  {
    col = "red"
  }
  else
  {
    col = "green"
  }
  p <- suppressWarnings(figure(xgrid = F, ygrid = F, xaxes = "below",yaxes = "left",
                               width = 500, height = 500, title = "GREAT analysis",
                               ylab = "Number of genes",xlab = "-log10(adjusted p-value)",
                               logo = NULL, tools = c( "pan", "wheel_zoom", "box_zoom",
                                                       "reset"))
                        %>% ly_points(data=inputdata, x = log10p,
                                      y = Hyper_Observed_Gene_Hits, color = col, alpha = 0.5,
                                    size = 10,hover = list(Term,adjusted_pvalue,gene_number)))
  suppressMessages(suppressWarnings(rbokeh2html(p, file = "rbokeh_temp.html")))
  html_contents <- get_html_content("rbokeh_temp.html")
  file.remove("rbokeh_temp.html")
  return(html_contents)
}

get_html_content <- function(filename){
  index_html = file(filename, "r")
  modelid = ""
  docid = ""
  docs_json = ""
  while (TRUE){
    line = readLines(index_html, n=1, warn = FALSE)
    if ( length(line) == 0 ) {
      break
    }
    line = trimws(line)
    if(startsWith(line, "var modelid")){
      modelid = gsub("[\r\n]", "", line)
    }
    if(startsWith(line, "var docid")){
      docid = gsub("[\r\n]", "", line)
    }
    if(startsWith(line, "var docs_json")){
      docs_json = gsub("[\r\n]", "", line)
    }
  }
  close(index_html)
  return(list(modelid, docid, docs_json))
}


