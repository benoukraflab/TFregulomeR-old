#' interactome3D
#'
#' This function allows you to get a html report of a 3D dynamic TF interactome with CpG methylation and external source signal.
#' @param intersectPeakMatrix Output of function 'intersectPeakMatrix()'.
#' @param return_interactome_with_mCpG Either TRUE of FALSE (default). If TRUE, html report of TF interactome with mCpG portion will be saved.
#' @param mCpG_threshold A mininum beta score to determine CpG methylation. Should be 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8 (default), or 0.9
#' @param return_interactome_with_external_source Either TRUE of FALSE (default). If TRUE, html report of TF interactome with external source signal will be saved.
#' @param external_source_value The value of external source signal in the intersected peaks. It should be one of the following values: "median" (default),"mean","SD","quartile_25","quartile_75".
#' @param angle_of_matrix Either "x" (default) or "y". If "x", will focus on the peak sets in "peak_list_x" intersected with "peak_list_y"; if "y", will focus on peak sets in "peak_list_y" intersected with "peak_list_x".
#' @return An html file
#' @keywords interactome3D
#' @export
#'

interactome3D <- function(intersectPeakMatrix,
                          return_interactome_with_mCpG = FALSE,
                          mCpG_threshold = 0.8,
                          return_interactome_with_external_source = FALSE,
                          external_source_value = "median",
                          angle_of_matrix = "x")
{
  # check input arguments
  if (missing(intersectPeakMatrix))
  {
    stop("Please provide output of 'intersectPeakMatrix()' using 'intersectPeakMatrix ='!")
  }
  # check the validity of input intersectPeakMatrix
  if (class(intersectPeakMatrix[1,1][[1]])[1] != "IntersectPeakMatrix")
  {
    stop("The input 'intersectPeakMatrix' is not valid. Please use the output of function 'intersectPeakMatrix()'")
  }
  if (!is.logical(return_interactome_with_mCpG))
  {
    stop("'return_interactome_with_mCpG' should be either TRUE (T) or FALSE (F, default)!")
  }
  if (!(mCpG_threshold %in% seq(0.1,0.9,0.1)))
  {
    stop("'mCpG_threshold' should be one of the following values: 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8 (default), and 0.9")
  }
  if (!is.logical(return_interactome_with_external_source))
  {
    stop("'return_interactome_with_external_source' should be either TRUE (T) or FALSE (F, default)!")
  }
  if (angle_of_matrix != "x" && angle_of_matrix != "y")
  {
    stop("'angle_of_matrix' should be either 'x' (default) or 'y'!")
  }
  if (!(external_source_value %in% c("median","mean","SD","quartile_25","quartile_75")))
  {
    stop("'external_source_value' should be one of the following values: 'median','mean','SD','quartile_25','quartile_75'")
  }

  # message
  message("Start interactome3D ...")
  if (return_interactome_with_mCpG)
  {
    message("... You chose to report TF interactome coupled with DNA methylation ...")
    message(paste0("... ... The mCpG threshold you chose is ", mCpG_threshold))
  }
  else
  {
    message("... You chose NOT to report TF interactome coupled with DNA methylation ...")
  }
  if (return_interactome_with_external_source)
  {
    message("... You chose to report TF interactome coupled with external source signal ...")
    message(paste0("... ... The external source signal value you chose is ", external_source_value))
  }
  else
  {
    message("... You chose NOT to report TF interactome coupled with external source signal ...")
  }
  if (return_interactome_with_mCpG == FALSE && return_interactome_with_external_source == FALSE)
  {
    message("... You chose no action. EXIT!!")
    return(NULL)
  }

  # start reporting
  # get intersectPeakMatrix result
  suppressMessages(intersectPeakMatrix_res <- intersectPeakMatrixResult(intersectPeakMatrix=intersectPeakMatrix,
                                                       return_intersection_matrix = TRUE,
                                                       angle_of_matrix = angle_of_matrix,
                                                       return_methylation_profile = TRUE,
                                                       angle_of_methylation_profile = angle_of_matrix,
                                                       return_external_source = TRUE,
                                                       external_source_value = external_source_value))
  TF_intersection_matrix <- intersectPeakMatrix_res$intersection_matrix
  TF_methylation_matrix <- intersectPeakMatrix_res$methylation_profile_matrix
  TF_external_source_matrix <- intersectPeakMatrix_res$external_source_matrix
  # get mCpG portion matrix
  TF_mCpG_matrix = as.data.frame(matrix(nrow = nrow(TF_intersection_matrix),
                                        ncol = ncol(TF_intersection_matrix)))
  rownames(TF_mCpG_matrix) <- rownames(TF_intersection_matrix)
  colnames(TF_mCpG_matrix) <- colnames(TF_intersection_matrix)
  for (i in 1:nrow(TF_intersection_matrix))
  {
    for (j in 1:ncol(TF_intersection_matrix))
    {
      meth_matrix_i <- TF_methylation_matrix[i,j][[1]]
      if (!(is.na(meth_matrix_i[1,1])) && sum(meth_matrix_i) >0 )
      {
        TF_mCpG_matrix[i,j] <- 100*sum(meth_matrix_i[(mCpG_threshold*10+1):10])/sum(meth_matrix_i)
      }
    }
  }

  # hierarchical clustering of TF intersection matrix if more than one column and one row
  if (nrow(intersectPeakMatrix)>1 && ncol(intersectPeakMatrix)>1)
  {
    TF_intersection_matrix_heatmap <- heamap.invisible(as.matrix(TF_intersection_matrix))
    TF_intersection_matrix_hc <- TF_intersection_matrix[rev(TF_intersection_matrix_heatmap$rowInd),
                                                        TF_intersection_matrix_heatmap$colInd]
    TF_mCpG_matrix_hc <- TF_mCpG_matrix[rev(TF_intersection_matrix_heatmap$rowInd),
                                        TF_intersection_matrix_heatmap$colInd]
    TF_external_source_matrix_hc <- TF_external_source_matrix[rev(TF_intersection_matrix_heatmap$rowInd),
                                                              TF_intersection_matrix_heatmap$colInd]
  }
  else
  {
    TF_intersection_matrix_hc <- TF_intersection_matrix
    TF_mCpG_matrix_hc <- TF_mCpG_matrix
    TF_external_source_matrix_hc <- TF_external_source_matrix
  }

  if (return_interactome_with_mCpG)
  {
    # TF interactome with mCpG
    TF_intersection_with_mCpG_json <- get.jsondata(TF_intersection_matrix_hc, TF_mCpG_matrix_hc)
    x_value <- get.names(rownames(TF_intersection_matrix_hc))
    y_value <- get.names(colnames(TF_intersection_matrix_hc))
    html_report_mCpG <- html_3d.report(TF_intersection_with_mCpG_json, x_value,
                                       y_value, "mCpG")
    write(html_report_mCpG, "TF_interactome_with_mCpG.html")
    message("... report of TF interactome with mCpG portion has been saved as 'TF_interactome_with_mCpG.html'")
  }

  if (return_interactome_with_external_source)
  {
    # TF interactome with external source signal jsondata
    TF_intersection_with_external_source_json <- get.jsondata(TF_intersection_matrix_hc,
                                                              TF_external_source_matrix_hc)
    x_value <- get.names(rownames(TF_intersection_matrix_hc))
    y_value <- get.names(colnames(TF_intersection_matrix_hc))
    html_report_external_source <- html_3d.report(TF_intersection_with_external_source_json,
                                                  x_value, y_value, "external_source")
    write(html_report_external_source, "TF_interactome_with_external_source.html")
    message("... report of TF interactome with external source signal has been saved as 'TF_interactome_with_external_source.html'")
  }
}

heamap.invisible <- function(data_matrix)
{
  ff <- tempfile()
  png(filename=ff)
  res <- gplots::heatmap.2(data_matrix)
  dev.off()
  unlink(ff)
  return(res)
}

get.jsondata <- function(intersection_table, z_table)
{
  jsondata <- "[\n"
  for (i in seq(0, nrow(intersection_table)-1,1))
  {
    for (j in seq(0, ncol(intersection_table)-1,1))
    {
      intersection_value <- intersection_table[i+1,j+1]
      z_value <- z_table[i+1,j+1]
      if (is.na(z_value))
      {
        z_value <- -1
      }
      jsondata <- paste0(jsondata,"{x:",i,",y:",j,",z:",z_value,",style:",intersection_value,"},\n")
    }
  }
  jsondata <- paste0(jsondata, "];\n")
}

get.names <- function(namelist)
{
  namelist_new <- unlist(lapply(namelist,
                         function(x) tail(unlist(strsplit(x, split = "_")),1)))
  namelist_new_str <- paste0("'", paste(namelist_new, collapse = "','"), "'")
}

html_3d.report <- function(jsondata, x_value, y_value, data_type)
{
  if (data_type == "mCpG")
  {
    z_value_lable <- "mCpGs(%)"
  }
  else
  {
    z_value_lable <- "external source signal"
  }
  html_3d_res <- paste0("<!DOCTYPE>
<html>
<head>
  <style>
    html, body {
      font: 10pt arial;
      padding: 0;
      margin: 0;
      width: 100%;
      height: 100%;
    }

  </style>
  <script type=\"text/javascript\" src=\"https://bioinfo-csi.nus.edu.sg/methmotif/interactome/js/vis.min.js\"></script>

  <script type=\"text/javascript\">
    var data = null;
    var graph = null;
    function drawVisualization() {
      var data = ", jsondata,"
      // specify options
      var options = {
        width:  '90%',
        height: '90%',
        style: 'bar-color',
        showPerspective: true,
        showGrid: false,
        showShadow: true,
        keepAspectRatio: true,
        verticalRatio: 0.4,
        showLegend:true,
        xLabel: \"TF-x\",
        yLabel: \"TF-y\",
        zLabel: '", z_value_lable,"',
        legendLabel: \"co-binding(%)\",
        xStep: 5,
        yStep: 5,
        zMin: 0,
        xValueLabel: function(value) {
           var x_name = [",x_value,"];
           return x_name[value];
        },
        yValueLabel: function(value) {
           var y_name = [",y_value,"];
           return y_name[value];
        },
      tooltip: function (point) {
          var x_name = [", x_value,"];
          var y_name = [", y_value,"];
          var output = 'TF-x: <b>'+x_name[point.x]+'</b><br>TF-y: <b>'+y_name[point.y]+'</b><br>",z_value_lable,": '+point.z;
          return output;
        },
      };
      var container = document.getElementById('mygraph');
      graph = new vis.Graph3d(container, data, options);

    }
  </script>

</head>
<body onresize=\"graph.redraw();\" onload=\"drawVisualization()\">
<div id=\"mygraph\" align=\"center\"></div>
</body>
</html>

")
  return(html_3d_res)
}


