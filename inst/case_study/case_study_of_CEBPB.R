library(TFregulomeR)
library(gplots)

# CEBPB motifs in TFregulomeR
CEBPB_record <- TFBSBrowser(tf = "CEBPB")

################### K562 shared CEBPB targets (Figure 2A, B, C and D left)###############################
K562_commonPeak <- commonPeaks(target_peak_id = "MM1_HSA_K562_CEBPB",
                                 motif_only_for_target_peak = T,  
                                 compared_peak_id = CEBPB_record$ID, 
                                 motif_only_for_compared_peak = T)
K562_commonPeak_result <- commonPeakResult(commonPeaks = K562_commonPeak,
                                           return_common_peak_sites = T, 
                                           save_MethMotif_logo = T)
K562_commonPeak_result_peak <- K562_commonPeak_result$common_peak_list$MM1_HSA_K562_CEBPB_common_peaks

# cofactor profiling across K562 shared CEBPB targets
K562_TFBS <- TFBSBrowser(cell_tissue_name = "K562")
K562_commonPeak_cofactor <- intersectPeakMatrix(user_peak_list_x = list(K562_commonPeak_result_peak),
                                                user_peak_x_id = "MM1_HSA_K562_CEBPB",
                                                peak_id_y = K562_TFBS$ID, 
                                                motif_only_for_id_y = T, 
                                                methylation_profile_in_narrow_region = T)
K562_commonPeak_cofactor_res <- intersectPeakMatrixResult(intersectPeakMatrix = K562_commonPeak_cofactor, 
                                                          return_intersection_matrix = T, 
                                                          return_methylation_profile = T)
K562_commonPeak_cofactor_matrix <- K562_commonPeak_cofactor_res$intersection_matrix
K562_commonPeak_cofactor_matrix_t <- as.data.frame(t(K562_commonPeak_cofactor_matrix))
attach(K562_commonPeak_cofactor_matrix_t)
K562_commonPeak_cofactor_matrix_order <- as.data.frame(K562_commonPeak_cofactor_matrix_t[order(-MM1_HSA_K562_CEBPB),,drop = FALSE])
detach(K562_commonPeak_cofactor_matrix_t)

color <- colorRampPalette(c("white","#D46A6A", "#801515", "#550000"))
bk <- c(seq(0, 100, length=100))
pdf("K562_shared_CEBPB_targets_co-factor_profile.pdf")
heatmap.2(as.matrix(cbind(K562_commonPeak_cofactor_matrix_order,K562_commonPeak_cofactor_matrix_order)),
          col=color(99), trace = "none", Colv = NULL, Rowv = NULL,
          dendrogram = "none", density.info = "none",
          key.xlab = "percentage (%)", keysize = 1.2, cexRow = .3, 
          key.title = "",labCol = NA, breaks = bk)
dev.off()


K562_commonPeak_cofactor_meth <- K562_commonPeak_cofactor_res$methylation_profile_matrix[,c("MM1_HSA_K562_CEBPB",
                                                                                            "MM1_HSA_K562_CEBPD",
                                                                                            "MM1_HSA_K562_CTCF",
                                                                                            "MM1_HSA_K562_E4F1",
                                                                                            "MM1_HSA_K562_JUND")]
K562_commonPeak_cofactor_meth_matrix <- matrix(nrow = 5, ncol = 10)
for (i in 1:length(K562_commonPeak_cofactor_meth)){
    meth_i <- K562_commonPeak_cofactor_meth[i][[1]]
    K562_commonPeak_cofactor_meth_matrix[i, ] <- 100*meth_i/sum(meth_i)
}
rownames(K562_commonPeak_cofactor_meth_matrix) <- c("CEBPB",
                                                    "CEBPD",
                                                    "CTCF",
                                                    "E4F1",
                                                    "JUND")
colnames(K562_commonPeak_cofactor_meth_matrix) <- rownames(K562_commonPeak_cofactor_meth[1][[1]])
color <- colorRampPalette(c("white", "#00B2D3", "#012694"))
bk <- c(seq(0, 100, length=100))
pdf("methylation_profile_of_K562_shared_CEBPB_targets_co-factors.pdf")
heatmap.2(as.matrix(K562_commonPeak_cofactor_meth_matrix),
          col=color(99), trace = "none", Colv = NULL, Rowv = NULL,
          dendrogram = "none", density.info = "none",
          key.xlab = "CpG percentage (%)", keysize = 1.2, cexRow = 1,  
          key.title = "",srtCol=45, breaks = bk)
dev.off()

# MethMotif logo of CEBPB/CEBPD in shared CEBPB targets
intersectPeakMatrixResult(intersectPeakMatrix = K562_commonPeak_cofactor, 
                          save_MethMotif_logo = T, 
                          saving_MethMotif_logo_y_id = "MM1_HSA_K562_CEBPD")
# K562 shared CEBPB targets without CEBPD co-binding regions
K562_common_no_cebpd <- exclusivePeaks(user_target_peak_list = list(K562_commonPeak_result_peak), 
                                       user_target_peak_id = "MM1_HSA_K562_CEBPB",
                                       excluded_peak_id = "MM1_HSA_K562_CEBPD", 
                                       motif_only_for_excluded_peak = T, 
                                       methylation_profile_in_narrow_region = T)
K562_common_no_cebpd_res <- exclusivePeakResult(exclusivePeaks = K562_common_no_cebpd, 
                                                save_MethMotif_logo = T)
################### K562 shared CEBPB targets (Figure 2A, B, C and D left) #######################

################### K562 exclusive CEBPB targets (Figure 2A, B, C and D right) #####################
#remove MM1_HSA_K562_CEBPB ID from all CEBPB TFregulomeR IDs
CEBPB_record_ID_noK562 <- CEBPB_record$ID[!(CEBPB_record$ID %in% "MM1_HSA_K562_CEBPB")]

K562_exclusivePeak_output <- exclusivePeaks(target_peak_id = "MM1_HSA_K562_CEBPB",
                                            motif_only_for_target_peak = T,
                                            excluded_peak_id = CEBPB_record_ID_noK562,
                                            motif_only_for_excluded_peak = T)
K562_exclusivePeak_result <- exclusivePeakResult(exclusivePeaks = K562_exclusivePeak_output, 
                                                 return_exclusive_peak_sites = T, 
                                                 save_MethMotif_logo = T)
K562_exclusivePeak_peak <- K562_exclusivePeak_result$exclusive_peak_list$MM1_HSA_K562_CEBPB_exclusive_peaks

# cofactor profiling across K562 exclusive CEBPB targets
K562_exclusivePeak_cofactor <- intersectPeakMatrix(user_peak_list_x = list(K562_exclusivePeak_peak),
                                                   user_peak_x_id = "MM1_HSA_K562_CEBPB",
                                                   peak_id_y = K562_TFBS$ID, 
                                                   motif_only_for_id_y = T,
                                                   methylation_profile_in_narrow_region = T)
K562_exclusivePeak_cofactor_res <- intersectPeakMatrixResult(intersectPeakMatrix = K562_exclusivePeak_cofactor,
                                                             return_intersection_matrix = T, 
                                                             return_methylation_profile = T)
K562_exclusivePeak_cofactor_matrix <- K562_exclusivePeak_cofactor_res$intersection_matrix

K562_exclusivePeak_cofactor_matrix_t <- as.data.frame(t(K562_exclusivePeak_cofactor_matrix))
attach(K562_exclusivePeak_cofactor_matrix_t)
K562_exclusivePeak_cofactor_order <- as.data.frame(K562_exclusivePeak_cofactor_matrix_t[order(-MM1_HSA_K562_CEBPB),,drop = FALSE])
detach(K562_exclusivePeak_cofactor_matrix_t)

color <- colorRampPalette(c("white","#D46A6A", "#801515", "#550000"))
bk <- c(seq(0, 100, length=100))
pdf("K562_exclusive_CEBPB_targets_co-factor_profile.pdf")
heatmap.2(as.matrix(cbind(K562_exclusivePeak_cofactor_order,K562_exclusivePeak_cofactor_order)),
          col=color(99), trace = "none", Colv = NULL, Rowv = NULL,
          dendrogram = "none", density.info = "none",
          key.xlab = "percentage (%)", keysize = 1.2, cexRow = .3,  
          key.title = "",labCol = NA, breaks = bk)
dev.off()


K562_exclusivePeak_cofactor_meth <- K562_exclusivePeak_cofactor_res$methylation_profile_matrix[,c("MM1_HSA_K562_CEBPB",
                                                                                                  "MM1_HSA_K562_ATF4",
                                                                                                  "MM1_HSA_K562_TCF12",
                                                                                                  "MM1_HSA_K562_CBFA2T3",
                                                                                                  "MM1_HSA_K562_TAL1",
                                                                                                  "MM1_HSA_K562_GATA1")]
K562_exclusivePeak_cofactor_meth_matrix <- matrix(nrow = 6, ncol = 10)
for (i in 1:length(K562_exclusivePeak_cofactor_meth)){
    meth_i <- K562_exclusivePeak_cofactor_meth[i][[1]]
    K562_exclusivePeak_cofactor_meth_matrix[i, ] <- 100*meth_i/sum(meth_i)
}

rownames(K562_exclusivePeak_cofactor_meth_matrix) <- c("CEBPB",
                                                       "ATF4",
                                                       "TCF12",
                                                       "CBFA2T3",
                                                       "TAL1",
                                                       "GATA1")

colnames(K562_exclusivePeak_cofactor_meth_matrix) <- rownames(K562_exclusivePeak_cofactor_meth[1][[1]])
color <- colorRampPalette(c("white", "#00B2D3", "#012694"))
bk <- c(seq(0, 100, length=100))
pdf("methylation_profile_of_K562_exclusive_CEBPB_targets_co-factors.pdf")
heatmap.2(as.matrix(K562_exclusivePeak_cofactor_meth_matrix),
          col=color(99), trace = "none", Colv = NULL, Rowv = NULL,
          dendrogram = "none", density.info = "none",
          key.xlab = "CpG percentage (%)", keysize = 1.2, cexRow = 1,  
          key.title = "",srtCol=45, breaks = bk)
dev.off()

# MethMotif logo of CEBPB/ATF4 in exclusive CEBPB targets
intersectPeakMatrixResult(intersectPeakMatrix = K562_exclusivePeak_cofactor, 
                          save_MethMotif_logo = T, saving_MethMotif_logo_y_id = "MM1_HSA_K562_ATF4")
# K562 exclusive CEBPB targets without ATF4 co-binding regions
K562_exclusive_no_atf4 <- exclusivePeaks(user_target_peak_list = list(K562_exclusivePeak_peak),
                                         user_target_peak_id = "MM1_HSA_K562_CEBPB",
                                         excluded_peak_id = c("MM1_HSA_K562_ATF4"),
                                         motif_only_for_excluded_peak = T,
                                         methylation_profile_in_narrow_region = T)
K562_exclusive_no_atf4_res <- exclusivePeakResult(exclusivePeaks = K562_exclusive_no_atf4, 
                                                  save_MethMotif_logo = T)
################### K562 exclusive CEBPB targets (Figure 2A, B, C and D right) #####################


# CEBPB Meth(Motif) in all cell types (Supplementary Figure1)
for (i in CEBPB_record$ID){
    motif_i <- searchMotif(id = i)
    plotLogo(MM_object = motif_i)
}

# motif in shared and exclusive CEBPB targets in all other 15 cell types (Supplementary Figure2)
for (i in CEBPB_record$ID){
    common_i <- commonPeaks(target_peak_id = i,
                           motif_only_for_target_peak = T, 
                           compared_peak_id = CEBPB_record$ID,
                           motif_only_for_compared_peak = T)
    common_i_res <- commonPeakResult(commonPeaks = common_i,save_MethMotif_logo = T)
    
    # exclude ID i from all CEBPB IDs
    cebpd_id_no_i <- CEBPB_record$ID[!(CEBPB_record$ID %in% i)]
    exclusive_i <- exclusivePeaks(target_peak_id = i,
                                 motif_only_for_target_peak = T,
                                 excluded_peak_id = cebpd_id_no_i, 
                                 motif_only_for_excluded_peak = T)
    exclusive_i_res <- exclusivePeakResult(exclusivePeaks = exclusive_i, save_MethMotif_logo = T)
}

# functions of CEBPB/CEBPD and CEBPB/ATF4 targets in K562 (Supplementary Figure3)
# load required package for GREAT annotation
library(rGREAT)
# load required package for genomic conversion from hg38 to hg19
library(liftOver)
# all CEBPB/CEBPD co-binding regions in K562
K562_CEBPB_CEBPD <- commonPeaks(target_peak_id = "MM1_HSA_K562_CEBPB",
                                motif_only_for_target_peak = T,
                                compared_peak_id = "MM1_HSA_K562_CEBPD", 
                                motif_only_for_compared_peak = T)
K562_CEBPB_CEBPD_res <- commonPeakResult(commonPeaks = K562_CEBPB_CEBPD,
                                         return_common_peak_sites = T, 
                                         save_MethMotif_logo = T, 
                                         return_summary = T)
K562_CEBPB_CEBPD_res$peak_summary
#>                                 percentage_in_original_inputs(%)
#> MM1_HSA_K562_CEBPB_common_peaks                         6.532727
K562_CEBPB_CEBPD_peak <- K562_CEBPB_CEBPD_res$common_peak_list$MM1_HSA_K562_CEBPB_common_peaks
K562_CEBPB_CEBPD_great <- greatAnnotate(peaks = K562_CEBPB_CEBPD_peak, return_annotation = T)
K562_CEBPB_CEBPD_great_bp <- K562_CEBPB_CEBPD_great[which(K562_CEBPB_CEBPD_great$category=="BP"),]

# all CEBPB/ATF4 co-binding regions in K562
K562_CEBPB_ATF4 <- commonPeaks(target_peak_id = "MM1_HSA_K562_CEBPB",
                               motif_only_for_target_peak = T,
                               compared_peak_id = "MM1_HSA_K562_ATF4",
                               motif_only_for_compared_peak = T)
K562_CEBPB_ATF4_res <- commonPeakResult(commonPeaks = K562_CEBPB_ATF4, 
                                       return_common_peak_sites = T, 
                                       save_MethMotif_logo = T, 
                                       return_summary = T)
K562_CEBPB_ATF4_res$peak_summary
#>                                 percentage_in_original_inputs(%)
#> MM1_HSA_K562_CEBPB_common_peaks                         32.34774
K562_CEBPB_ATF4_peak <- K562_CEBPB_ATF4_res$common_peak_list$MM1_HSA_K562_CEBPB_common_peaks
K562_CEBPB_ATF4_great <- greatAnnotate(peaks = K562_CEBPB_ATF4_peak, return_annotation = T)
K562_CEBPB_ATF4_great_bp <- K562_CEBPB_ATF4_great[which(K562_CEBPB_ATF4_great$category=="BP"),]
