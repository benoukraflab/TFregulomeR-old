#' cofactorReport
#'
#' This function allows you to get a PDF report of top cofactors along with DNA methylation for a TF.
#' @param intersectPeakMatrix Output of function 'intersectPeakMatrix()'.
#' @param top_num Number of highest co-binding factors to report for each TF (up to 10). By default the number is 10.
#' @param cobinding_threshold Only the co-factors with co-binding percentages more than this threshold value will be reported. By default the threshold is 0.05.
#' @return A PDF file
#' @keywords cofactorReport
#' @export
#' @examples
#' peak_id_x <- c("MM1_HSA_K562_CEBPB")
#' peak_id_y <- c("MM1_HSA_K562_CEBPD", "MM1_HSA_K562_ATF4")
#' intersect_output <- intersectPeakMatrix(peak_id_x=peak_id_x,
#'                                         motif_only_for_id_x=TRUE,
#'                                         peak_id_y=peak_id_y)
#' cofactorReport(intersectPeakMatrix = intersect_output)

cofactorReport <- function(intersectPeakMatrix,
                           top_num = 10,
                           cobinding_threshold = 0.05)
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
    # check other input arguments
    if (top_num <= 0)
    {
        stop("Invalid input for 'top_num'. It should be > 0 but <=10")
    }
    if (top_num > 10)
    {
        stop("This function only reports no more than 10 cofactors for each TF.")
    }
    if (cobinding_threshold <= 0 || cobinding_threshold > 1)
    {
        stop("'cobinding_threshold' should be >0 but <=1.")
    }

    # start reporting
    message("Start cofactorReport ...")
    message(paste0("... The maximum number of cofactors to be reported is ", top_num))
    message(paste0("... The minimum percent of co-binding peaks for a cofactor is ", cobinding_threshold*100,"%"))
    message("... Each peak set derived from TFregulomeR compendium in 'peak list x' will be reported in an individual PDF file")

    for (i in seq(1,nrow(intersectPeakMatrix),1))
    {
        intersectPeakMatrix_i <- intersectPeakMatrix[i, ,drop=FALSE]
        id_i <- rownames(intersectPeakMatrix_i)
        # judge whether peak set id_i is derived from TFregulomeR compendium
        is_from_TFregulomeR <- intersectPeakMatrix_i[1,1][[1]]@isxTFregulomeID
        if (!is_from_TFregulomeR)
        {
            message(paste0("... ... peak id '",id_i,"' in 'peak list x' is NOT a TFregulomeR ID. SKIP!!"))
            next
        }
        message(paste0("... ... Start reporting peak id '",id_i,"' ..."))
        suppressMessages(intersectPeakMatrix_res_i <- intersectPeakMatrixResult(intersectPeakMatrix_i,
                                                                                return_intersection_matrix = TRUE,
                                                                                angle_of_matrix = "x"))
        intersect_matrix_i <- intersectPeakMatrix_res_i$intersection_matrix
        intersect_matrix_t_i <- as.data.frame(t(intersect_matrix_i))
        intersect_matrix_filter_i <- intersect_matrix_t_i[(intersect_matrix_t_i[,1] >= cobinding_threshold*100),,drop=FALSE]

        if (nrow(intersect_matrix_filter_i) == 0)
        {
            message(paste0("... ... ... The number of cofactors passing 'cobinding_threshold' for peak id '", id_i,"' is 0. SKIP!!"))
            next
        }
        intersect_matrix_order_i <- intersect_matrix_filter_i[order(intersect_matrix_filter_i[,1])
                                                              ,,drop = FALSE]
        if (nrow(intersect_matrix_filter_i) > top_num)
        {
            message(paste0("... ... ... The number of cofactors passing 'cobinding_threshold' for peak id '",
                           id_i,"' is ", nrow(intersect_matrix_filter_i), ". Only top ",
                           top_num," will be selected."))
            intersect_matrix_order_i <- intersect_matrix_order_i[seq(nrow(intersect_matrix_filter_i)-top_num+1,
                                                                     nrow(intersect_matrix_filter_i),1)
                                                                 ,,drop=FALSE]
        }
        # cobinding heatmap
        intersect_matrix_heatmap_i <- as.data.frame(matrix(nrow=nrow(intersect_matrix_order_i),
                                                           ncol = 5))
        colnames(intersect_matrix_heatmap_i) <- c("x","new_x","y","new_y","value")
        intersect_matrix_heatmap_i$y <- rownames(intersect_matrix_order_i)
        intersect_matrix_heatmap_i$new_y <- unlist(lapply(rownames(intersect_matrix_order_i),
                                                          function(x) tail(unlist(strsplit(x,split = "_")),1)))
        intersect_matrix_heatmap_i$x <- colnames(intersect_matrix_order_i)
        intersect_matrix_heatmap_i$new_x <- unlist(lapply(colnames(intersect_matrix_order_i),
                                                          function(x) tail(unlist(strsplit(x,split = "_")),1)))
        intersect_matrix_heatmap_i$value <- intersect_matrix_order_i[,1]
        intersect_matrix_heatmap_i$new_y <- factor(intersect_matrix_heatmap_i$new_y,
                                                   levels = as.character(intersect_matrix_heatmap_i$new_y))
        cobinding_ylabel <- as.character(intersect_matrix_heatmap_i$new_y[rev(seq(1,nrow(intersect_matrix_heatmap_i),1))])
        cobinding_ylabel_new <- paste0(as.character(intersect_matrix_heatmap_i$new_x),
                                       "\n+\n",cobinding_ylabel)
        colors_cobinding <- colorRampPalette(c("white","#D46A6A", "#801515", "#550000"))(11)
        p1 <- ggplot(intersect_matrix_heatmap_i,aes(x=intersect_matrix_heatmap_i$new_x,
                                                    y=intersect_matrix_heatmap_i$new_y,
                                                    fill=intersect_matrix_heatmap_i$value))+
            geom_tile()+scale_fill_gradientn(colours=colors_cobinding,
                                             breaks=c(seq(0, 100, length=11)),
                                             limits=c(0,100))+
            theme(panel.background = element_blank(),
                  plot.background = element_blank(),
                  axis.title=element_blank(),
                  axis.ticks = element_blank(),
                  plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"),
                  axis.text.x = element_text(size=6),axis.text.y = element_text(size=6),
                  legend.position = "none")+
            scale_y_discrete(breaks=cobinding_ylabel,
                             labels=cobinding_ylabel_new)
        # cobinding heatmap legend
        legend1_matrix <- as.data.frame(matrix(nrow = 11, ncol = 2))
        colnames(legend1_matrix) <- c("cobinding","value")
        legend1_matrix[,1] <- "x"
        legend1_matrix[,2] <- seq(0,100,10)
        p_legend1 <- ggplot(legend1_matrix,aes(x=legend1_matrix$value,
                                               y=legend1_matrix$cobinding,
                                               fill=legend1_matrix$value))+
            geom_tile()+scale_fill_gradientn(colours=colors_cobinding,
                                             breaks=c(seq(0, 100, length=11)),
                                             limits=c(0,100))+
            theme(panel.background = element_blank(),
                  plot.background = element_blank(),axis.title=element_blank(),
                  axis.ticks = element_blank(),
                  plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"),
                  axis.text = element_text(size=10),legend.position = "none")+
            scale_y_discrete(breaks=c("x"),labels=c("co-binding(%)"))
        # motifs
        motif_plot_list_p <- list()
        count <- 1
        for (j in order(-seq(1,nrow(intersect_matrix_heatmap_i),1)))
        {
            x_j <- intersect_matrix_heatmap_i$x[j]
            y_j <- intersect_matrix_heatmap_i$y[j]
            motif_j <- t(intersectPeakMatrix_i[x_j,y_j][[1]]@MethMotif_x@MMmotif@motif_matrix)
            top <- 10
            bottom <- 0
            if (j==1)
            {
                bottom <- 10
                top <- 0
            }
            p_j <- ggplot() + geom_logo(data = motif_j, method = "bits") +
                theme(axis.title.y=element_blank(), axis.text.y=element_blank(),
                      axis.ticks.y=element_blank(),
                      axis.ticks.x=element_blank(),axis.text.x=element_blank(),
                      plot.margin = margin(t = top, r = 0, b = bottom, l = 0, unit = "pt"),
                      panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(),
                      panel.background = element_blank())
            motif_plot_list_p[[count]] <- p_j
            count <- count+1
        }
        if (nrow(intersect_matrix_heatmap_i)==1)
        {
            p2 <- arrangeGrob(grobs=motif_plot_list_p,
                              layout_matrix = as.matrix(c(NA,1,2)))
        }
        else if(nrow(intersect_matrix_heatmap_i)==2)
        {
            p2 <- arrangeGrob(grobs=motif_plot_list_p,
                              layout_matrix = as.matrix(c(NA,1,NA,NA,2,3)))
        }
        else if (nrow(intersect_matrix_heatmap_i)==3)
        {
            p2 <- arrangeGrob(grobs=motif_plot_list_p,
                              layout_matrix = as.matrix(c(NA,1,1,1,NA,
                                                          NA,2,2,2,NA,
                                                          NA,3,3,3,4)))
        }
        else
        {
            p2 <- arrangeGrob(grobs=motif_plot_list_p,
                              nrow=nrow(intersect_matrix_heatmap_i))
        }

        # methylation within motif heatmap
        if (!is.na(intersectPeakMatrix_i[1,1][[1]]@MethMotif_x@MMBetaScore[1,1]) )
        {
            methylation_matrix_inside <- as.data.frame(matrix(nrow = nrow(intersect_matrix_order_i)*3,
                                                              ncol = 5))
            colnames(methylation_matrix_inside) <- c("x","y","new_y","meth_interval","value")
            methylation_matrix_inside$x <- colnames(intersect_matrix_order_i)
            methylation_matrix_inside$y <- unlist(lapply(rownames(intersect_matrix_order_i),
                                                         function(x) rep(x,3)))
            methylation_matrix_inside$new_y <- unlist(lapply(methylation_matrix_inside$y,
                                                             function(x) tail(unlist(strsplit(x,split = "_")),1)))

            for(j in seq(1,nrow(intersect_matrix_heatmap_i),1))
            {
                x_j <- intersect_matrix_heatmap_i$x[j]
                y_j <- intersect_matrix_heatmap_i$y[j]
                meth_matrix <- intersectPeakMatrix_i[x_j,y_j][[1]]@MethMotif_x@MMBetaScore
                if (sum(meth_matrix)>0)
                {
                    unmeth_j <- 100*sum(meth_matrix[1,])/sum(meth_matrix)
                    between_j <- 100*sum(meth_matrix[2,])/sum(meth_matrix)
                    meth_j <- 100*sum(meth_matrix[3,])/sum(meth_matrix)
                    methylation_matrix_inside[(methylation_matrix_inside$x==x_j & methylation_matrix_inside$y==y_j),
                                              "meth_interval"] <- c("unmeth","in-between","meth")
                    methylation_matrix_inside[(methylation_matrix_inside$x==x_j & methylation_matrix_inside$y==y_j),
                                              "value"] <- c(unmeth_j,between_j,meth_j)
                }
                else
                {
                    methylation_matrix_inside[(methylation_matrix_inside$x==x_j & methylation_matrix_inside$y==y_j),
                                              "meth_interval"] <- c("unmeth","in-between","meth")
                    methylation_matrix_inside[(methylation_matrix_inside$x==x_j & methylation_matrix_inside$y==y_j),
                                              "value"] <- c(0,0,0)
                }
            }
            methylation_matrix_inside$new_y <- factor(methylation_matrix_inside$new_y,
                                                      levels = as.character(intersect_matrix_heatmap_i$new_y))

            methylation_matrix_inside$meth_interval <- factor(methylation_matrix_inside$meth_interval,
                                                              levels = as.character(c("unmeth","in-between","meth")))
            colors_meth <- colorRampPalette(c("white","#7887AB", "#4F628E", "#2E4172","#061539"))(11)
            p3 <- ggplot(methylation_matrix_inside,aes(x=methylation_matrix_inside$meth_interval,
                                                       y=methylation_matrix_inside$new_y,
                                                       fill=methylation_matrix_inside$value))+
                geom_tile()+scale_fill_gradientn(colours=colors_meth,
                                                 breaks=c(seq(0, 100, length=11)),
                                                 limits=c(0,100))+
                theme(legend.position = "none" ,panel.background = element_blank(),
                      plot.background = element_blank(),axis.title=element_blank(),
                      axis.ticks = element_blank(),
                      plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"),
                      axis.text.y = element_blank(),
                      axis.text.x = element_text(size=6))
        }
        else
        {
            message(paste("... ... ... No methylation information within motif for id '",
                          id_i,"', since it is not from MethMotif."))
            p3 <- ggplot() + theme(panel.grid.major = element_blank(),
                                   panel.grid.minor = element_blank(),
                                   panel.background = element_blank())
        }

        # methylation surrounding heatmap
        if (!is.na(intersectPeakMatrix_i[1,1][[1]]@methylation_profile_x[1,1]))
        {
            methylation_matrix <- as.data.frame(matrix(nrow = nrow(intersect_matrix_order_i)*3,
                                                       ncol = 5))
            colnames(methylation_matrix) <- c("x","y","new_y","meth_interval","value")
            methylation_matrix$x <- colnames(intersect_matrix_order_i)
            methylation_matrix$y <- unlist(lapply(rownames(intersect_matrix_order_i),
                                                  function(x) rep(x,3)))
            methylation_matrix$new_y <- unlist(lapply(methylation_matrix$y,
                                                      function(x) tail(unlist(strsplit(x,split = "_")),1)))

            for (j in seq(1,nrow(intersect_matrix_order_i),1))
            {
                x_j <- intersect_matrix_heatmap_i$x[j]
                y_j <- intersect_matrix_heatmap_i$y[j]
                meth_matrix <- intersectPeakMatrix_i[x_j,y_j][[1]]@methylation_profile_x
                if (sum(meth_matrix)>0)
                {
                    unmeth_j <- 100*meth_matrix[1]/sum(meth_matrix)
                    between_j <- 100*sum(meth_matrix[2:9])/sum(meth_matrix)
                    meth_j <- 100*meth_matrix[10]/sum(meth_matrix)
                    methylation_matrix[(methylation_matrix$x==x_j & methylation_matrix$y==y_j),
                                       "meth_interval"] <- c("unmeth","in-between","meth")
                    methylation_matrix[(methylation_matrix$x==x_j & methylation_matrix$y==y_j),
                                       "value"] <- c(unmeth_j,between_j,meth_j)
                }
                else
                {
                    methylation_matrix[(methylation_matrix$x==x_j & methylation_matrix$y==y_j),
                                       "meth_interval"] <- c("unmeth","in-between","meth")
                    methylation_matrix[(methylation_matrix$x==x_j & methylation_matrix$y==y_j),
                                       "value"] <- c(0,0,0)
                }
            }
            methylation_matrix$new_y <- factor(methylation_matrix$new_y,
                                               levels = as.character(intersect_matrix_heatmap_i$new_y))

            methylation_matrix$meth_interval <- factor(methylation_matrix$meth_interval,
                                                       levels = as.character(c("unmeth","in-between","meth")))
            p4 <- ggplot(methylation_matrix,aes(x=methylation_matrix$meth_interval,
                                                y=methylation_matrix$new_y,
                                                fill=methylation_matrix$value))+
                geom_tile()+scale_fill_gradientn(colours=colors_meth,
                                                 breaks=c(seq(0, 100, length=11)),
                                                 limits=c(0,100))+
                theme(panel.background = element_blank(),
                      plot.background = element_blank(),
                      axis.title=element_blank(),
                      axis.ticks = element_blank(),
                      plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"),
                      axis.text.y = element_blank(),
                      axis.text.x = element_text(size=6), legend.position = "none")
        }
        else
        {
            message(paste("... ... ... No methylation information in 200bp peak regions for id '",
                          id_i,"', since it is not from MethMotif or 'methylation_profile_in_narrow_region=FALSE' during intersectPeakMatrix()."))
            p4 <- ggplot() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                   panel.background = element_blank())
        }

        # methylation legend
        p_legend2 <- ggplot(legend1_matrix,aes(x=legend1_matrix$value,
                                               y=legend1_matrix$cobinding,
                                               fill=legend1_matrix$value))+
            geom_tile()+scale_fill_gradientn(colours=colors_meth,
                                             breaks=c(seq(0, 100, length=11)),
                                             limits=c(0,100))+
            theme(panel.background = element_blank(),
                  plot.background = element_blank(),
                  axis.title=element_blank(),
                  axis.ticks = element_blank(),
                  plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"),
                  axis.text = element_text(size=10),legend.position = "none")+
            scale_y_discrete(breaks=c("x"),labels=c("CG percent(%)"))

        # read enrichment score
        suppressMessages(tag_median_res_i <- intersectPeakMatrixResult(intersectPeakMatrix_i,
                                                                       return_tag_density = TRUE,
                                                                       tag_density_value = "median"))
        tag_median_i <- tag_median_res_i$tag_density_matrix
        tag_median_i_t <- as.data.frame(t(tag_median_i))
        tag_median_i_t_filter <- tag_median_i_t[intersect_matrix_heatmap_i$y,
                                                unique(intersect_matrix_heatmap_i$x),drop=FALSE]
        colnames(tag_median_i_t_filter) = "y"
        tag_median_i_t_filter$x = seq(1,nrow(tag_median_i_t_filter),1)

        suppressMessages(tag_q1_res_i <- intersectPeakMatrixResult(intersectPeakMatrix_i,
                                                                   return_tag_density = TRUE,
                                                                   tag_density_value = "quartile_25"))
        tag_q1_i <- tag_q1_res_i$tag_density_matrix
        tag_q1_i_t <- as.data.frame(t(tag_q1_i))
        tag_q1_i_t_filter <- tag_q1_i_t[intersect_matrix_heatmap_i$y,
                                        unique(intersect_matrix_heatmap_i$x),drop=FALSE]
        colnames(tag_q1_i_t_filter) = "y"
        tag_q1_i_t_filter$x = seq(1,nrow(tag_q1_i_t_filter),1)

        suppressMessages(tag_q3_res_i <- intersectPeakMatrixResult(intersectPeakMatrix_i,
                                                                   return_tag_density = TRUE,
                                                                   tag_density_value = "quartile_75"))
        tag_q3_i <- tag_q3_res_i$tag_density_matrix
        tag_q3_i_t <- as.data.frame(t(tag_q3_i))
        tag_q3_i_t_filter <- tag_q3_i_t[intersect_matrix_heatmap_i$y,
                                        unique(intersect_matrix_heatmap_i$x),drop=FALSE]
        colnames(tag_q3_i_t_filter) = "y"
        tag_q3_i_t_filter$x = seq(1,nrow(tag_q3_i_t_filter),1)

        tag_max_value <- max(c(tag_median_i_t_filter[,"y"],
                               tag_q1_i_t_filter[,"y"],
                               tag_q3_i_t_filter[,"y"]))
        tag_max_value_int <- (as.integer(tag_max_value/10)+1)*10
        tag_x_max_list <- c(1,0.7,0.5,0.5,0.4,0.4,0.3,0.3,0,0)
        tag_x_min <- 1-tag_x_max_list[nrow(tag_median_i_t_filter)]
        tag_x_max <- nrow(tag_median_i_t_filter) +
            tag_x_max_list[nrow(tag_median_i_t_filter)]

        tag_quartile_total1 <- rbind(tag_q1_i_t_filter,tag_q3_i_t_filter)
        colnames(tag_q1_i_t_filter) <- c("y1","x1")
        colnames(tag_q3_i_t_filter) <- c("y3","x3")
        tag_quartile_total2 <- cbind(tag_q1_i_t_filter,tag_q3_i_t_filter)
        p_tag <- ggplot(tag_median_i_t_filter, aes(x=tag_median_i_t_filter$x,
                                                   y = tag_median_i_t_filter$y))+
            geom_point(shape=19, size=3)+
            ylim(0,tag_max_value_int)+
            xlim(tag_x_min,tag_x_max)+
            theme(panel.background = element_blank(),
                  plot.background = element_blank(),
                  plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"),
                  axis.title.y = element_blank(),
                  axis.text.y=element_blank(),
                  axis.ticks.y=element_blank(),
                  axis.title.x = element_text(size=6),
                  axis.text.x = element_text(size=6))+
            geom_point(data = tag_quartile_total1, shape = "|", size = 3,
                       mapping = aes(x=tag_quartile_total1$x,
                                     y=tag_quartile_total1$y))+
            geom_segment(data = tag_quartile_total2,
                         mapping=aes(x=tag_quartile_total2$x1,
                                     y=tag_quartile_total2$y1,
                                     xend=tag_quartile_total2$x3,
                                     yend=tag_quartile_total2$y3))+
            labs(y="Read enrichment score")+
            coord_flip()



        text1 <- textGrob("cofactors")
        text2 <- textGrob("motif")
        text3 <- textGrob("methylation in motif")
        text4 <- textGrob("methylation in\n 200bp regions")
        text5 = textGrob("read enrichment")

        pdf_i_name <- paste0(id_i,"_cofactor_report.pdf")
        pdf(pdf_i_name)
        grid.arrange(text1,text2,text3,text4,text5,p1, p2,p3,p4,p_tag,p_legend1,p_legend2,
                     layout_matrix = rbind(c(1,2,2,3,3,4,4,5,5),
                                           c(6,7,7,8,8,9,9,10,10),
                                           c(6,7,7,8,8,9,9,10,10),
                                           c(6,7,7,8,8,9,9,10,10),
                                           c(6,7,7,8,8,9,9,10,10),
                                           c(6,7,7,8,8,9,9,10,10),
                                           c(6,7,7,8,8,9,9,10,10),
                                           c(6,7,7,8,8,9,9,10,10),
                                           c(6,7,7,8,8,9,9,10,10),
                                           c(6,7,7,8,8,9,9,10,10),
                                           c(6,7,7,8,8,9,9,10,10),
                                           c(6,7,7,8,8,9,9,10,10),
                                           c(6,7,7,8,8,9,9,10,10),
                                           c(NA,11,11,11,NA,12,12,12,NA)))
        dev.off()
        message(paste0("... ... ... Cofactor report for id '", id_i,"' has been saved as ", pdf_i_name))
    }
}

