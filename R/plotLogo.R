#' plot (Meth)Motif logo
#'
#' This function allows you to plot (Meth)Motif logo.
#' @param MM_object Required. MethMotif class object
#' @param logo_type Logo type for the (Meth)Motif logo to be saved, either "entropy" (default) or "frequency".
#' @param meth_level Methylation level to be plot for the (Meth)Motif logo, and it should be one of the values, "all" (default), "methylated", and "unmethylated".
#' @return  (Meth)Motif logo pdf file. If the TFregulomeR peak source is from MethMotif, MethMotif logo will be saved; if the source is GTRD, only motif logo will be saved.
#' @keywords plotLogo
#' @export
#' @examples
#' K562_CEBPB <- searchMotif(id = "MM1_HSA_K562_CEBPB")
#' plotLogo(MM_object = K562_CEBPB)

plotLogo <- function(MM_object, logo_type="entropy", meth_level="all")
{
  # check logo_type
  if (logo_type != "entropy" && logo_type != "frequency")
  {
    stop("Please check your input argument 'logo_type'! Please choose either 'entropy' (default) or 'frequency'!")
  }
  else if (logo_type == "entropy")
  {
    logo_method <- "bits"
    y_max <- 2
  }
  else if (logo_type == "frequency")
  {
    logo_method <- "prob"
    y_max <- 1
  }

  # check logo_type
  if (meth_level != "all" && meth_level != "methylated" && meth_level != "unmethylated")
  {
    stop("Please check your input argument 'meth_level'! Please choose one of 'all' (default), 'methylated' or 'unmethylated'!")
  }

  # check input argumen
  if (missing(MM_object))
  {
    stop("Please input an MethMotif_object using 'MM_object = '!")
  }
  else if (class(MM_object)[1] != "MethMotif")
  {
    stop("Your input is not a MethMotif obejct. Try to use searchMotif() function to return a MethMotif object for the TFBS of interest!")
  }
  else
  {
    # get beta score matrix and motif matrix
    MMBetaScore <- MM_object@MMBetaScore
    MMmotif <- MM_object@MMmotif
    ID <- MMmotif@id
    motif_matrix <- t(MMmotif@motif_matrix)

    motif_length <- ncol(motif_matrix)
    nPeaks <- MM_object@MMmotif@nPeaks
    title_for_nPeaks <- paste0("Number of peaks with motif = ", nPeaks)

    if (!(is.na(MMBetaScore[1,1])))
    {
      # generate a dataframe for beta score plotting
      plot_beta_score <- matrix(rep(0,length(MMBetaScore)*3), ncol = 3)
      colnames(plot_beta_score) <- c("number","pos","meth")
      plot_beta_score[seq(1,motif_length,1),1] <- as.vector(MMBetaScore[3,])
      plot_beta_score[(motif_length+1):(2*motif_length),1] <- as.vector(MMBetaScore[2,])
      plot_beta_score[(2*motif_length+1):(3*motif_length),1] <- as.vector(MMBetaScore[1,])
      plot_beta_score[seq(1,motif_length,1),2] <- seq(1, motif_length, 1)
      plot_beta_score[(motif_length+1):(2*motif_length),2] <- seq(1, motif_length, 1)
      plot_beta_score[(2*motif_length+1):(3*motif_length),2] <- seq(1, motif_length,1)
      plot_beta_score[seq(1,motif_length,1), 3] <- "beta score>90%"
      plot_beta_score[(motif_length+1):(2*motif_length),3] <- "beta score 10-90%"
      plot_beta_score[(2*motif_length+1):(3*motif_length),3] <- "beta score<10%"
      plot_beta_score <- as.data.frame(plot_beta_score)

      # plot all, methylated or unmethylated
      if (meth_level == "all")
      {
        # make levels in beta score plotting matrix
        plot_beta_score$meth <- factor(plot_beta_score$meth,levels = c("beta score>90%",  "beta score 10-90%","beta score<10%"))
        plot_beta_score$pos <- factor(plot_beta_score$pos, levels = seq(1,motif_length,1))
        ylim <- round(max(as.vector(apply(MMBetaScore,2,sum)))/1000+1)*1000+500
        barplot_color <- c("darkorange1","darkgreen", "dodgerblue1")
        pdf_name <- paste0(ID, "-logo-", logo_type, ".pdf")
      }
      else if (meth_level == "methylated")
      {
        plot_beta_score <- plot_beta_score[which(plot_beta_score$meth == "beta score>90%"), ]
        plot_beta_score$pos <- factor(plot_beta_score$pos, levels = seq(1,motif_length,1))
        ylim <- round(max(as.vector(MMBetaScore[3, ]))/1000+1)*1000+500
        barplot_color <- c("darkorange1")
        pdf_name <- paste0(ID, "-logo-", logo_type, "-methylated-only.pdf")
      }
      else
      {
        plot_beta_score <- plot_beta_score[which(plot_beta_score$meth == "beta score<10%"), ]
        plot_beta_score$pos <- factor(plot_beta_score$pos, levels = seq(1,motif_length,1))
        ylim <- round(max(as.vector(MMBetaScore[1, ]))/1000+1)*1000+500
        barplot_color <- c("dodgerblue1")
        pdf_name <- paste0(ID, "-logo-", logo_type, "-unmethylated-only.pdf")
      }

      sum_of_pos <- aggregate(as.numeric(as.character(plot_beta_score$number)),
                              by=list(pos=plot_beta_score$pos),
                              FUN=sum)
      colnames(sum_of_pos) <- c("pos", "sum")
      #plot beta score
      p1 <- ggplot(data = plot_beta_score[order(plot_beta_score$meth, decreasing = FALSE),],
                   aes(x=pos,y=as.numeric(as.character(plot_beta_score$number)),
                       fill=plot_beta_score$meth)) +
        geom_bar(colour="black", stat="identity") +
        scale_fill_manual(values = barplot_color) + ylim(0, ylim) +
        theme(axis.title.y=element_blank(), axis.title.x=element_blank(), axis.text.y=element_blank(),
              axis.ticks.y=element_blank(), axis.ticks.x=element_blank(), axis.text.x=element_blank(),
              legend.title = element_blank(), legend.background = element_blank(),
              legend.box.background = element_rect(colour = "black"), legend.key.size = unit(0.8,"line"),
              legend.position=c(0.9,0.9), plot.margin = margin(t = 10, r = 20, b = 0, l = 35, unit = "pt"),
              panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              panel.background = element_blank(),
              plot.title = element_text(hjust = 0.5, size = 10)) +
        stat_summary(fun.y = sum, aes(label = stat(sum_of_pos$sum), group = pos), geom = "text",vjust = -0.5)+
        ggtitle(title_for_nPeaks)
    }
    else
    {
      p1 <- ggplot() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                panel.background = element_blank(),
                             plot.title = element_text(hjust = 0.5, size = 10))+
        ggtitle(title_for_nPeaks)
      pdf_name = paste0(ID, "-logo-", logo_type, ".pdf")
    }

    # motif logo position size
    #size xlab
    if (motif_length>40)
    {
      xlab_size <- 4
    }
    else
    {
      xlab_size <- -0.5*motif_length+24
    }
    #plot motif logo
    p2 <- ggplot() + geom_logo(data = motif_matrix, method = logo_method) +
        theme(axis.title.y=element_blank(),
              axis.ticks.x=element_blank(),axis.text.x=element_text(size=xlab_size),
              plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit =  "pt"),
              panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              panel.background = element_blank())+ylim(0,y_max)

    #combine p1 and p2 together
    p3 <- grid.arrange(p1,p2, nrow=2)
    pdf(pdf_name)
    plot(p3)
    dev.off()
    message(paste0("Success: a PDF named '", pdf_name,"' has been saved!"))
  }
}
