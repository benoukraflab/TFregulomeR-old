#' genomeAnnotate
#'
#' This function allows you to annotate genomic locations of cis-regulatory regions.
#' @param peaks Required. A bed-format genomic regions in data frame.
#' @param assembly Optional. The assembly the input regions are using, currently supporting 'hg19', 'hg38' (default), 'mm9' and 'mm10'.
#' @param return_annotation Optional. Either TRUE of FALSE (default). If TRUE, a data.frame containing annotation results will be returned.
#' @param return_html_report Optional. Either TRUE of FALSE (default). If TRUE, a dynamic HTML report will be saved.
#' @param promoter_range Optional. A numeric vector to define promoter range. By default, c(-1000, 100) defines promoters as 1000bp upstream and 100bp downstream of TSS.
#' @param TTS_range Optional. A numeric vector to define TTS range. By default, c(-100, 1000) defines promoters as 100bp upstream and 1000bp downstream of real TTS.
#' @param TFregulome_url Optional. If the MethMoitf url is NO more "http://bioinfo-csi.nus.edu.sg/methmotif/", please use a new url.
#' @return  a data.frame, or an HTML report depending on the options.
#' @keywords genomeAnnotate
#' @export
#' @examples
#' K562_CEBPB_regions <- loadPeaks(id = "MM1_HSA_K562_CEBPB")
#' K562_CEBPB_regions_annotation <- genomeAnnotate(peaks = K562_CEBPB_regions,
#' return_annotation = T, return_html_report = T)

genomeAnnotate <- function(peaks, assembly = "hg38", return_annotation = FALSE,
                          return_html_report = FALSE, promoter_range = c(-1000,100),
                          TTS_range = c(-100, 1000), TFregulome_url)
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
  if (!is.vector(promoter_range))
  {
    stop("'promoter_range' should be a numric vector, by default c(-1000, 100)")
  }
  if (!is.vector(TTS_range))
  {
    stop("'TTS_range' should be a numric vector, by default c(-100, 1000)")
  }
  # check loaded package
  if (assembly == "hg38" && !("TxDb.Hsapiens.UCSC.hg38.knownGene" %in% (.packages())))
  {
    stop("R package 'TxDb.Hsapiens.UCSC.hg38.knownGene' (>=3.4.0) is NOT loaded yet!")
  }
  if (assembly == "hg19" && !("TxDb.Hsapiens.UCSC.hg19.knownGene" %in% (.packages())))
  {
    stop("R package 'TxDb.Hsapiens.UCSC.hg19.knownGene' (>=3.2.2) is NOT loaded yet!")
  }
  if (assembly == "mm10" && !("TxDb.Mmusculus.UCSC.mm10.knownGene" %in% (.packages())))
  {
    stop("R package 'TxDb.Mmusculus.UCSC.mm10.knownGene' (>=3.4.4) is NOT loaded yet!")
  }
  if (assembly == "mm9" && !("TxDb.Mmusculus.UCSC.mm9.knownGene" %in% (.packages())))
  {
    stop("R package 'TxDb.Mmusculus.UCSC.mm9.knownGene' (>=3.2.2) is NOT loaded yet!")
  }

  #message
  message("Start genomeAnnotate ...")
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

  # make an appropriate API url
  if (missing(TFregulome_url)){
    TFregulome_url <- "http://bioinfo-csi.nus.edu.sg/methmotif/api/methmotifR/genomeAnnotate/"
  } else if (endsWith(TFregulome_url, suffix = "/index.php")==TRUE){
    TFregulome_url <- gsub("index.php", "", TFregulome_url)
    TFregulome_url <- paste0(TFregulome_url, "api/methmotifR/genomeAnnotate/")
  } else if (endsWith(TFregulome_url, suffix = "/")==TRUE){
    TFregulome_url <- paste0(TFregulome_url, "api/methmotifR/genomeAnnotate/")
  } else {
    TFregulome_url <- paste0(TFregulome_url, "/api/methmotifR/genomeAnnotate/")
  }

  #check existence of geneName conversion file in methmotif server
  name_conversion_file <- paste0(TFregulome_url, "hg38_UCSC_to_GeneName.txt")
  name_conversion <- tryCatch(read.table(name_conversion_file, sep = "\t"),
                              warning=function(w) data.frame())
  if (nrow(name_conversion) ==0)
  {
    message("There is a warning to connect MethMotif API!")
    message("Advice:")
    message("1) Check internet access;")
    message("2) Current default MethMotif homepage is 'http://bioinfo-csi.nus.edu.sg/methmotif/'. If MethMotif homepage url is no more valid, please Google 'MethMotif', and input the valid MethMotif homepage url using 'TFregulome_url = '.")
    message(paste0("warning: ",cond))
    return(NULL)
  }

  peaks <- peaks[,1:3]
  colnames(peaks) <- c("chr","start", "end")
  peaks$id <- paste0("genomeAnnotate_peak_", as.vector(rownames(peaks)))
  peaks_gr <- with(peaks, GRanges(chr, IRanges(start+1, end), id=id))
  #load knowgene
  if (assembly=="hg38")
  {
    txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
    seqlevels(txdb) <- paste0(rep("chr",times=24), c(seq(1,22,1),"X","Y"))
  } else if (assembly=="hg19")
  {
    txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
    seqlevels(txdb) <- paste0(rep("chr",times=24), c(seq(1,22,1),"X","Y"))
  } else if (assembly=="mm10")
  {
    txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
    seqlevels(txdb) <- paste0(rep("chr",times=21), c(seq(1,19,1),"X","Y"))
    name_conversion_file <- paste0(TFregulome_url, "mm10_UCSC_to_GeneName.txt")
    name_conversion <- tryCatch(read.table(name_conversion_file, sep = "\t"),
                                warning=function(w) data.frame())
  } else
  {
    txdb <- TxDb.Mmusculus.UCSC.mm9.knownGene
    seqlevels(txdb) <- paste0(rep("chr",times=21), c(seq(1,19,1),"X","Y"))
    name_conversion_file <- paste0(TFregulome_url, "mm10_UCSC_to_GeneName.txt")
    name_conversion <- tryCatch(read.table(name_conversion_file, sep = "\t"),
                                warning=function(w) data.frame())
  }
  colnames(name_conversion) <- c("UCSC","gene_name")
  UCSC_name <- as.character(name_conversion$UCSC)
  gene_name <- as.character(name_conversion$gene_name)

  all_TSS <- promoters(txdb, upstream = 0, downstream = 1)
  # promoter
  message(paste0("... ... annotating promoters defined as upstream ",
                 abs(promoter_range[1]), "bp and downstream ",
                 abs(promoter_range[2]), "bp"))
  all_promoters <- promoters(txdb, upstream = abs(promoter_range[1]),
                             downstream = abs(promoter_range[2]))
  suppressWarnings(promoter_overlapped <- subsetByOverlaps(peaks_gr, all_promoters))
  # combine overlapped motif seq information with peak
  hits <- findOverlaps(peaks_gr, all_promoters)
  feature_id <- CharacterList(split(names(all_promoters)[subjectHits(hits)],
                                    queryHits(hits)))
  mcols(promoter_overlapped) <- DataFrame(mcols(promoter_overlapped), feature_id)
  promoter_overlapped_df <- as.data.frame(promoter_overlapped)
  new_promoter <- addDistanceAndClean(promoter_overlapped_df, all_TSS, name_conversion,"promoter-TSS")

  # TTS
  message(paste0("... ... annotating TTS defined as upstream ",
                 abs(TTS_range[1]), "bp and downstream ",
                 abs(TTS_range[2]), "bp"))
  peaks <- peaks[!(peaks$id %in% new_promoter$id), ]
  if (nrow(peaks)>0)
  {
    peaks_gr <- with(peaks, GRanges(chr, IRanges(start+1, end), id=id))
    # get TTS
    all_transcript <- transcripts(txdb, use.names=T)
    all_transcript_df <- as.data.frame(all_transcript)
    all_transcript_df_pos <- all_transcript_df[all_transcript_df$strand=="+",]
    pos_newStart <- as.integer(all_transcript_df_pos$end) - abs(TTS_range[1])
    pos_newEnd <- as.integer(all_transcript_df_pos$end) + abs(TTS_range[2])
    all_transcript_df_pos$start <- pos_newStart
    all_transcript_df_pos$end <- pos_newEnd
    all_transcript_df_neg <- all_transcript_df[all_transcript_df$strand=="-",]
    neg_newStart <- as.integer(all_transcript_df_neg$start) - abs(TTS_range[2])
    neg_newEnd <- as.integer(all_transcript_df_neg$start) + abs(TTS_range[1])
    all_transcript_df_neg$start <- neg_newStart
    all_transcript_df_neg$end <- neg_newEnd
    all_transcript_df_new <- rbind(all_transcript_df_pos, all_transcript_df_neg)
    all_tts <- with(all_transcript_df_new, GRanges(seqnames, IRanges(start, end),
                                                   strand = strand, tx_id = tx_id,
                                                   tx_name=tx_name))
    names(all_tts) <- all_transcript_df_new$tx_name
    # annotating TTS
    suppressWarnings(tts_overlapped <- subsetByOverlaps(peaks_gr, all_tts))
    # combine overlapped motif seq information with peak
    hits <- findOverlaps(peaks_gr, all_tts)
    feature_id <- CharacterList(split(names(all_tts)[subjectHits(hits)],
                                      queryHits(hits)))
    mcols(tts_overlapped) <- DataFrame(mcols(tts_overlapped), feature_id)
    tts_overlapped_df <- as.data.frame(tts_overlapped)
    new_tts <- addDistanceAndClean(tts_overlapped_df, all_TSS, name_conversion,"TTS")
  } else
  {
    new_tts <- as.data.frame(matrix(nrow = 0,ncol = 8))
    colnames(new_tts) <- c("chr","start","end","id","annotation","geneName",
                           "transcript", "distanceToTSS")
  }

  # cds exons
  message("... ... annotating exons")
  peaks <- peaks[!(peaks$id %in% new_tts$id), ]
  if (nrow(peaks)>0)
  {
    peaks_gr <- with(peaks, GRanges(chr, IRanges(start+1, end), id=id))
    all_exons <- exonsBy(txdb, by = "tx",use.names=T)
    suppressWarnings(exons_overlapped <- subsetByOverlaps(peaks_gr, all_exons))
    # combine overlapped motif seq information with peak
    hits <- findOverlaps(peaks_gr, all_exons)
    feature_id <- CharacterList(split(names(all_exons)[subjectHits(hits)],
                                      queryHits(hits)))
    mcols(exons_overlapped) <- DataFrame(mcols(exons_overlapped), feature_id)
    exons_overlapped_df <- as.data.frame(exons_overlapped)
    new_exon <- addDistanceAndClean(exons_overlapped_df, all_TSS, name_conversion, "exon")
  } else
  {
    new_exon <- as.data.frame(matrix(nrow = 0,ncol = 8))
    colnames(new_exon) <- c("chr","start","end","id","annotation","geneName",
                           "transcript", "distanceToTSS")
  }


  # 5UTR exon
  message("... ... annotating 5' UTR")
  peaks <- peaks[!(peaks$id %in% new_exon$id), ]
  if (nrow(peaks)>0)
  {
    peaks_gr <- with(peaks, GRanges(chr, IRanges(start+1, end), id=id))
    all_5utr <- fiveUTRsByTranscript(txdb, use.names=T)
    suppressWarnings(fiveUtr_overlapped <- subsetByOverlaps(peaks_gr, all_5utr))
    # combine overlapped motif seq information with peak
    hits <- findOverlaps(peaks_gr, all_5utr)
    feature_id <- CharacterList(split(names(all_5utr)[subjectHits(hits)],
                                      queryHits(hits)))
    mcols(fiveUtr_overlapped) <- DataFrame(mcols(fiveUtr_overlapped), feature_id)
    fiveUtr_overlapped_df <- as.data.frame(fiveUtr_overlapped)
    new_5utr <- addDistanceAndClean(fiveUtr_overlapped_df, all_TSS, name_conversion, "5UTR")
  } else
  {
    new_5utr <- as.data.frame(matrix(nrow = 0,ncol = 8))
    colnames(new_5utr) <- c("chr","start","end","id","annotation","geneName",
                            "transcript", "distanceToTSS")
  }


  # 3UTR exon
  message("... ... annotating 3' UTR")
  peaks <- peaks[!(peaks$id %in% new_5utr$id), ]
  if (nrow(peaks)>0)
  {
    peaks_gr <- with(peaks, GRanges(chr, IRanges(start+1, end), id=id))
    all_3utr <- threeUTRsByTranscript(txdb, use.names=T)
    suppressWarnings(threeUtr_overlapped <- subsetByOverlaps(peaks_gr, all_3utr))
    # combine overlapped motif seq information with peak
    hits <- findOverlaps(peaks_gr, all_3utr)
    feature_id <- CharacterList(split(names(all_3utr)[subjectHits(hits)],
                                      queryHits(hits)))
    mcols(threeUtr_overlapped) <- DataFrame(mcols(threeUtr_overlapped), feature_id)
    threeUtr_overlapped_df <- as.data.frame(threeUtr_overlapped)
    new_3utr <- addDistanceAndClean(threeUtr_overlapped_df, all_TSS, name_conversion, "3UTR")
  } else
  {
    new_3utr <- as.data.frame(matrix(nrow = 0,ncol = 8))
    colnames(new_3utr) <- c("chr","start","end","id","annotation","geneName",
                            "transcript", "distanceToTSS")
  }

  #introns
  message("... ... annotating introns")
  peaks <- peaks[!(peaks$id %in% new_3utr$id), ]
  if (nrow(peaks)>0)
  {
    peaks_gr <- with(peaks, GRanges(chr, IRanges(start+1, end), id=id))
    all_introns <- intronsByTranscript(txdb, use.names=T)
    suppressWarnings(intron_overlapped <- subsetByOverlaps(peaks_gr, all_introns))
    # combine overlapped motif seq information with peak
    hits <- findOverlaps(peaks_gr, all_introns)
    feature_id <- CharacterList(split(names(all_introns)[subjectHits(hits)],
                                      queryHits(hits)))
    mcols(intron_overlapped) <- DataFrame(mcols(intron_overlapped), feature_id)
    intron_overlapped_df <- as.data.frame(intron_overlapped)
    new_intron <- addDistanceAndClean(intron_overlapped_df, all_TSS, name_conversion, "intron")
  } else
  {
    new_intron <- as.data.frame(matrix(nrow = 0,ncol = 8))
    colnames(new_intron) <- c("chr","start","end","id","annotation","geneName",
                            "transcript", "distanceToTSS")
  }


  #intergenic
  message("... ... annotating intergenic regions")
  peaks <- peaks[!(peaks$id %in% new_intron$id), ]
  new_intergenic <- as.data.frame(matrix(nrow = nrow(peaks),
                                         ncol = 8))
  if (nrow(peaks)>0)
  {
    peaks_gr <- with(peaks, GRanges(chr, IRanges(start+1, end), id=id))
    near_start <- follow(peaks_gr, unstrand(all_TSS))
    near_end <- precede(peaks_gr, unstrand(all_TSS))
    start_indx <- seq(1,length(near_start),1)
    distance_start <- rep(Inf, length(near_start))
    distance_start[!(is.na(near_start))] <- unlist(lapply(start_indx[!(is.na(near_start))],
                                                          function(x) distance(peaks_gr[x], all_TSS[near_start[x]])))
    transcript_start <- rep(NA, length(near_start))
    transcript_start[!(is.na(near_start))] <- names(all_TSS[near_start[!(is.na(near_start))]])
    genename_start <- rep(NA, length(near_start))
    genename_start[!(is.na(near_start))] <- unlist(lapply(transcript_start[!(is.na(near_start))],
                                                          function(x) paste0(gene_name[which(UCSC_name %in% x)], collapse = ",")))

    end_indx <- seq(1,length(near_end),1)
    distance_end <- rep(Inf, length(near_end))
    distance_end[!(is.na(near_end))] <- unlist(lapply(end_indx[!(is.na(near_end))],
                                                      function(x) distance(peaks_gr[x], all_TSS[near_end[x]])))
    transcript_end <- rep(NA, length(near_end))
    transcript_end[!(is.na(near_end))] <- names(all_TSS[near_end[!(is.na(near_end))]])
    genename_end <- rep(NA, length(near_end))
    genename_end[!(is.na(near_end))] <- unlist(lapply(transcript_end[!(is.na(near_end))],
                                                      function(x) paste0(gene_name[which(UCSC_name %in% x)], collapse = ";")))

    distance_both <- data.frame(start = distance_start, end=distance_end)
    transcript_both <- data.frame(start = transcript_start, end=transcript_end)
    genename_both <- data.frame(start = genename_start, end=genename_end)
    min_index <- apply(distance_both, 1, function(x) which.min(x))

    distance_final <- distance_start
    distance_final[min_index==2] <- distance_end[min_index==2]

    transcript_final <- transcript_start
    transcript_final[min_index==2] <- transcript_end[min_index==2]

    genename_final <- genename_start
    genename_final[min_index==2] <- genename_end[min_index==2]

    new_intergenic[,1:4] <- peaks[,1:4]
    new_intergenic[,5] <- "intergenic"
    new_intergenic[,6] <- genename_final
    new_intergenic[,7] <- transcript_final
    new_intergenic[,8] <- distance_final
  }
  colnames(new_intergenic) <- c("chr","start","end","id","annotation",
                                "geneName","transcript", "distanceToTSS")
  annotation_result <- rbind(new_promoter, new_tts, new_exon, new_5utr,
                             new_3utr, new_intron, new_intergenic)

  if (return_html_report)
  {
    output_html <- html_report(annotation_result)
    write(output_html, "genomeAnnotate_result.html")
    message("... ... An html report has been generated as 'genomeAnnotate_result.html'!")
  }
  if (return_annotation)
  {
    message("... ... The annotation results have been returned in a data.frame!")
    return(annotation_result)
  }
}

addDistanceAndClean <- function(feature_df, all_TSS, name_conversion, annotation)
{
  UCSC_name <- as.character(name_conversion$UCSC)
  gene_name <- as.character(name_conversion$gene_name)
  new_feature_df <- as.data.frame(matrix(nrow = nrow(feature_df), ncol = 8))
  if (nrow(feature_df)>0)
  {
    feature_df_peak <- feature_df[,c("seqnames","start","end","id")]
    feature_df_peak_gr <- with(feature_df_peak, GRanges(seqnames, IRanges(start, end), id=id))
    feature_df$distance_tss <- unlist(lapply(seq(1,nrow(feature_df),1),
                                             function(x) paste0(distance(feature_df_peak_gr[x],all_TSS[feature_df[x,"feature_id"][[1]]]), collapse = ";")))
    feature_df$transcript <- unlist(apply(as.data.frame(feature_df$feature_id),
                                          1,function(x) paste0(as.character(x[[1]]), collapse = ";")))
    feature_df$gene_name <- unlist(lapply(feature_df$feature_id,
                                          function(x) paste0(unique(gene_name[which(UCSC_name %in% x[[1]])]), collapse = ";")))
    new_feature_df[,1:4] <- feature_df[,c("seqnames","start","end","id")]
    new_feature_df[,5] <-  annotation
    new_feature_df[,6] <- feature_df[,"gene_name"]
    new_feature_df[,7] <- feature_df[,"transcript"]
    new_feature_df[,8] <- feature_df[,"distance_tss"]
  }
  colnames(new_feature_df) <- c("chr","start","end","id","annotation",
                                "geneName","transcript", "distanceToTSS")
  return(new_feature_df)
}

html_report <- function(annotation_result){
  promoter_num <- nrow(annotation_result[annotation_result$annotation=="promoter-TSS", ])
  tts_num <- nrow(annotation_result[annotation_result$annotation=="TTS", ])
  exon_num <- nrow(annotation_result[annotation_result$annotation=="exon", ])
  fiveUtr_num <- nrow(annotation_result[annotation_result$annotation=="5UTR", ])
  threeUtr_num <- nrow(annotation_result[annotation_result$annotation=="3UTR", ])
  intron_num <- nrow(annotation_result[annotation_result$annotation=="intron", ])
  intergenic_num <- nrow(annotation_result[annotation_result$annotation=="intergenic", ])
  total_num <- promoter_num+tts_num+exon_num+fiveUtr_num+threeUtr_num+intron_num+intergenic_num

  table_html <- formGenomeAnnotateTableBody(annotation_result)
  output_html <- paste0("<!DOCTYPE html>
<meta charset='utf-8'>
<head>
  <style type='text/css'>
    body {
        font-size: 100%;
        font-family: 'Helvetica Neue', Helvetica, Arial, sans-serif;
    }
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
    table td, table th
    {
       max-width: 150px;
       word-wrap: break-word;
    }
  </style>
  <link rel='stylesheet' href='https://maxcdn.bootstrapcdn.com/bootstrap/3.3.7/css/bootstrap.min.css'>
  <script src='https://code.jquery.com/jquery-3.2.1.min.js' ></script>
  <script src='https://d3js.org/d3.v3.min.js' charset='utf-8'></script>
</head>
<body>
<div align='center'>
  <h4>TFregulomeR - genomeAnnotate Result</h4>
</div>
<hr align='center' style='width: 80%'>
<div align='center' style='width: 80%'>
  <div id='donut-charts' align='center'></div>
</div>
<br>
<div class='container' align='center'>
  Table - genomic annotations of the peaks: <select id='change_go_term' onchange='changeTermTable(this.value)'>
    <option value='promoter' selected>promoter</option>
    <option value='TTS'>TTS</option>
    <option value='exon'>exon</option>
    <option value='5utr'>5 UTR</option>
    <option value='3utr'>3 UTR</option>
    <option value='intron'>intron</option>
    <option value='intergenic'>intergenic</option>
  </select>
  <button type='button' class='btn btn-primary btn-sm' onclick=\"exportTableToCSV('genomeAnnotate_results.csv')\">Export All To CSV File</button>
  <div class='row' style='width:90%;height:400px;overflow:auto;'>
    <div class='table-responsive'>
      <table class='table table-bordred table-striped'>
        <thead>
          <th>chr</th>
          <th>start</th>
          <th>end</th>
          <th>annotation</th>
          <th>geneName</th>
          <th>transcript</th>
          <th>distanceToTSS</th>
        </thead>",table_html,"
      </table>
    </div>
  </div>
</div>
<hr align='center' style='width: 80%'>
<script>
$(function() {
  var donutData = genData();
  var donuts = new DonutCharts();
  donuts.create(donutData);
});
function DonutCharts() {
  var charts = d3.select('#donut-charts');
  var legendWidth = 100;
  var width = 200;
  var height = 320;
  var legendRectSize = 20;
  var legendSpacing = 4;
  var chart_m,
      chart_r,
      color = d3.scale.category10();
  var getCatNames = function(dataset) {
    var catNames = new Array();
    for (var i = 0; i < dataset[0].data.length; i++) {
      catNames.push(dataset[0].data[i].cat);
    }
    return catNames;
  }
  var createLegend = function(catNames) {
    var legends = charts.select('.legend')
                        .selectAll('g')
                        .data(catNames)
                        .enter().append('g')
                        .attr('transform', function(d, i) {
                          var height = legendRectSize + legendSpacing;
                          var offset =  height * color.domain().length / 4;
                          var vert = i * height - offset;
                          return 'translate(' + width / 2 + ',' + vert + ')';
                        });
    legends.append('rect')
           .attr('width', legendRectSize)
           .attr('height', legendRectSize)
           .style('fill', function(d, i) {
                return color(i);
           });
    legends.append('text')
           .attr('x', legendRectSize + legendSpacing)
           .attr('y', legendRectSize - legendSpacing)
           .text(function(d) {
                return d;
           });
  }
  var createCenter = function(pie) {
    var eventObj = {
        'click': function(d, i) {
          var paths = charts.selectAll('.clicked');
          pathAnim(paths, 0);
          paths.classed('clicked', false);
          resetAllCenterText();
        }
    }
    var donuts = d3.selectAll('.donut');
    donuts.append(\"svg:circle\")
          .attr(\"r\", chart_r * 1)
          .style(\"fill\", \"#E7E7E7\")
          .on(eventObj);
    donuts.append('text')
          .attr('class', 'center-txt type')
          .attr('y', chart_r * -0.16)
          .attr('text-anchor', 'middle')
          .style('font-weight', 'bold')
          .text(function(d, i) {
             return d.type;
          });
    donuts.append('text')
          .attr('class', 'center-txt value')
          .attr('text-anchor', 'middle');
    donuts.append('text')
          .attr('class', 'center-txt percentage')
          .attr('y', chart_r * 0.16)
          .attr('text-anchor', 'middle')
          .style('fill', '#A2A2A2');
  }
  var setCenterText = function(thisDonut) {
    var sum = d3.sum(thisDonut.selectAll('.clicked').data(), function(d) {
      return d.data.val;
    });
    thisDonut.select('.value')
             .text(function(d) {
               return (sum)? sum.toFixed(0) + d.unit
                           : d.total.toFixed(0) + d.unit;
             });
    thisDonut.select('.percentage')
             .text(function(d) {
               return (sum)? (sum/d.total*100).toFixed(2) + '%'
                           : '';
             });
  }
  var resetAllCenterText = function() {
    charts.selectAll('.value')
          .text(function(d) {
             return d.total.toFixed(0) + d.unit;
           });
    charts.selectAll('.percentage')
          .text('');
  }
  var pathAnim = function(path, dir) {
    switch(dir) {
        case 0:
            path.transition()
                .duration(1000)
                .ease('bounce')
                .attr('d', d3.svg.arc()
                             .innerRadius(chart_r * 0.3)
                             .outerRadius(chart_r)
                );
            break;
        case 1:
            path.transition()
                .attr('d', d3.svg.arc()
                             .innerRadius(chart_r * 0.5)
                             .outerRadius(chart_r *1.15)
                 );
            break;
    }
  }
  var updateDonut = function() {
    var eventObj = {
        'mouseover': function(d, i, j) {
            pathAnim(d3.select(this), 1);
            var thisDonut = charts.select('.type' + j);
            thisDonut.select('.value').text(function(donut_d) {
                return d.data.val.toFixed(0) + donut_d.unit;
            });
            thisDonut.select('.percentage').text(function(donut_d) {
                return (d.data.val/donut_d.total*100).toFixed(2) + '%';
            });
        },
        'mouseout': function(d, i, j) {
            var thisPath = d3.select(this);
            if (!thisPath.classed('clicked')) {
                pathAnim(thisPath, 0);
            }
            var thisDonut = charts.select('.type' + j);
            setCenterText(thisDonut);
        },
        'click': function(d, i, j) {
            var thisDonut = charts.select('.type' + j);
            if (0 === thisDonut.selectAll('.clicked')[0].length) {
                thisDonut.select('circle').on('click')();
            }
            var thisPath = d3.select(this);
            var clicked = thisPath.classed('clicked');
            pathAnim(thisPath, ~~(!clicked));
            thisPath.classed('clicked', !clicked);
            setCenterText(thisDonut);
        }
    };
    var pie = d3.layout.pie()
                       .sort(null)
                       .value(function(d) {
                            return d.val;
                        });
    var arc = d3.svg.arc()
                    .innerRadius(chart_r * 0.3)
                    .outerRadius(function() {
                        return (d3.select(this).classed('clicked'))? chart_r * 1.08
                                                                   : chart_r;
                    });
    var paths = charts.selectAll('.donut')
                      .selectAll('path')
                      .data(function(d, i) {
                        return pie(d.data);
                      });
    paths
        .transition()
        .duration(1000)
        .attr('d', arc);
    paths.enter()
         .append('svg:path')
         .attr('d', arc)
         .style('fill', function(d, i) {
            return color(i);
          })
         .style('stroke', '#FFFFFF')
         .on(eventObj)
    paths.exit().remove();
    resetAllCenterText();
  }
  this.create = function(dataset) {
    var $charts = $('#donut-charts');
    chart_m = $charts.innerWidth() / dataset.length / 4 * 0.13;
    chart_r = $charts.innerWidth() / dataset.length / 4 * 0.85;
    charts.append('svg')
          .attr('class', 'legend')
          .attr('width', width+legendWidth)
          .attr('height', height)
          .attr('align','center')
          .attr('transform', 'translate(' + (width / 2) + ',' + (height / 2) + ')');
    var donut = charts.selectAll('.donut')
                      .data(dataset)
                      .enter().append('svg:svg')
                      .attr('width', (chart_r + chart_m) * 2)
                      .attr('height', (chart_r + chart_m) * 2)
                      .append('svg:g')
                      .attr('class', function(d, i) {
                        return 'donut type' + i;
                      })
                      .attr('transform', 'translate(' + (chart_r+chart_m) + ',' + (chart_r+chart_m) + ')');
    createLegend(getCatNames(dataset));
    createCenter();
    updateDonut();
  }
  this.update = function(dataset) {
    var donut = charts.selectAll(\".donut\")
                      .data(dataset);
    updateDonut();
  }
}
function genData() {
    var type = ['Peaks'];
    var unit = [''];
    var cat = ['promoter-TSS', 'TTS', 'exon', '5 UTR', '3 UTR', 'intron', 'intergenic'];
    var dataset = new Array();
    for (var i = 0; i < type.length; i++) {
        var data = new Array();
        var total = ",total_num,"
    data = [
        {cat: 'promoter-TSS', val: ",promoter_num,"},
        {cat: 'TTS', val: ",tts_num,"},
        {cat: 'exon', val: ",exon_num,"},
        {cat: '5 UTR', val: ",fiveUtr_num,"},
        {cat: '3 UTR', val: ",threeUtr_num,"},
        {cat: 'intron', val: ",intron_num,"},
        {cat: 'intergenic', val: ",intergenic_num,"}
        ];
    dataset.push({
        'type': type[i],
        'unit': unit[i],
        'data': data,
        'total': total
    });
    }
    return dataset;
}
function changeTermTable(opt){
  if (opt == 'promoter'){
    $('#promoter').css('display', '');
    $('#TTS').css('display', 'none');
    $('#exon').css('display', 'none');
    $('#5UTR').css('display', 'none');
    $('#3UTR').css('display', 'none');
    $('#intron').css('display', 'none');
    $('#intergenic').css('display', 'none');
  }
  if (opt == 'TTS'){
    $('#promoter').css('display', 'none');
    $('#TTS').css('display', '');
    $('#exon').css('display', 'none');
    $('#5UTR').css('display', 'none');
    $('#3UTR').css('display', 'none');
    $('#intron').css('display', 'none');
    $('#intergenic').css('display', 'none');
  }
  if (opt == 'exon'){
    $('#promoter').css('display', 'none');
    $('#TTS').css('display', 'none');
    $('#exon').css('display', '');
    $('#5UTR').css('display', 'none');
    $('#3UTR').css('display', 'none');
    $('#intron').css('display', 'none');
    $('#intergenic').css('display', 'none');
  }
  if (opt == '5utr'){
    $('#promoter').css('display', 'none');
    $('#TTS').css('display', 'none');
    $('#exon').css('display', 'none');
    $('#5UTR').css('display', '');
    $('#3UTR').css('display', 'none');
    $('#intron').css('display', 'none');
    $('#intergenic').css('display', 'none');
  }
  if (opt == '3utr'){
    $('#promoter').css('display', 'none');
    $('#TTS').css('display', 'none');
    $('#exon').css('display', 'none');
    $('#5UTR').css('display', 'none');
    $('#3UTR').css('display', '');
    $('#intron').css('display', 'none');
    $('#intergenic').css('display', 'none');
  }
  if (opt == 'intron'){
    $('#promoter').css('display', 'none');
    $('#TTS').css('display', 'none');
    $('#exon').css('display', 'none');
    $('#5UTR').css('display', 'none');
    $('#3UTR').css('display', 'none');
    $('#intron').css('display', '');
    $('#intergenic').css('display', 'none');
  }
  if (opt == 'intergenic'){
    $('#promoter').css('display', 'none');
    $('#TTS').css('display', 'none');
    $('#exon').css('display', 'none');
    $('#5UTR').css('display', 'none');
    $('#3UTR').css('display', 'none');
    $('#intron').css('display', 'none');
    $('#intergenic').css('display', '');
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
</body>")

  return(output_html)
}


formGenomeAnnotateTableBody <- function(annotation_result)
{
  annotation_result$tag = apply(annotation_result,1,
                                function(x) paste0("          <tr>\n",
                                                   "            <td>",x["chr"],
                                                   "</td>\n",
                                                   "            <td>",x["start"],
                                                   "</td>\n",
                                                   "            <td>",x["end"],
                                                   "</td>\n",
                                                   "            <td>",x["annotation"],
                                                   "</td>\n",
                                                   "            <td>",x["geneName"],
                                                   "</td>\n",
                                                   "            <td>",x["transcript"],
                                                   "</td>\n",
                                                   "            <td>",x["distanceToTSS"],
                                                   "</td>\n          </tr>\n"))
  all_annotation <- c("promoter-TSS", "TTS", "exon", "5UTR", "3UTR", "intron", "intergenic")
  table_html=""
  for (annotation_i in all_annotation){
    annotation_result_i <- annotation_result[annotation_result$annotation==annotation_i,]
    annotation_result_i_tag <- paste0(annotation_result_i$tag, collapse = "\n")
    if (annotation_i == "promoter-TSS"){
      table_html <- paste0(table_html, "
        <tbody id='promoter' style='display:;'>\n",annotation_result_i_tag,"
        </tbody>")
    } else{
      table_html <- paste0(table_html,"
        <tbody id='",annotation_i,"' style='display:none;'>\n",annotation_result_i_tag,"
        </tbody>")
    }
    }
  return(table_html)
}
