library(TFregulomeR)
library(gplots)

MAFF_TFBS <- TFBSBrowser(tf = "MAFF")
# plot (Meth)Motif logos for MAFF in all cell types (Figure 3A)
for (i in MAFF_TFBS$ID){
    motif_i <- searchMotif(id = i)
    plotLogo(MM_object = motif_i)
}

# shared MAFF targets across all three cell types (Figure 3B)
MAFF_common_output <- commonPeaks(target_peak_id = MAFF_TFBS$ID,
                                  motif_only_for_target_peak = T,
                                  compared_peak_id = MAFF_TFBS$ID,
                                  motif_only_for_compared_peak = T)
MAFF_common_result <- commonPeakResult(commonPeaks = MAFF_common_output,
                                       return_common_peak_sites = T)
K562_MAFF_commonPeaks <- MAFF_common_result$common_peak_list$MM1_HSA_K562_MAFF_common_peaks
HeLa_MAFF_commonPeaks <- MAFF_common_result$common_peak_list$`MM1_HSA_HeLa-S3_MAFF_common_peaks`
HepG2_MAFF_commonPeaks <- MAFF_common_result$common_peak_list$MM1_HSA_HepG2_MAFF_common_peaks

############### K562 exclusive MAFF targets (Figure 3B, C and D) #################################
K562_MAFF_exclu <- exclusivePeaks(target_peak_id = "MM1_HSA_K562_MAFF",
                                  motif_only_for_target_peak = T,
                                  excluded_peak_id = c("MM1_HSA_HeLa-S3_MAFF",
                                                       "MM1_HSA_HepG2_MAFF"),
                                  motif_only_for_excluded_peak = T)
K562_MAFF_exclu_res <- exclusivePeakResult(exclusivePeaks = K562_MAFF_exclu,
                                           return_exclusive_peak_sites = T,
                                           save_MethMotif_logo = T)
K562_MAFF_exclu_peaks <- K562_MAFF_exclu_res$exclusive_peak_list$MM1_HSA_K562_MAFF_exclusive_peaks

# cofactors in K562 exclusive MAFF targets
K562_TFBS <- TFBSBrowser(cell_tissue_name = "K562")
K562_MAFF_exclu_peaks_intersect <- intersectPeakMatrix(user_peak_list_x = list(K562_MAFF_exclu_peaks),
                                                       user_peak_x_id = "MM1_HSA_K562_MAFF",
                                                       peak_id_y = K562_TFBS$ID,
                                                       motif_only_for_id_y = T,
                                                       methylation_profile_in_narrow_region = T)

K562_MAFF_exclu_peaks_intersect_res <- intersectPeakMatrixResult(intersectPeakMatrix = K562_MAFF_exclu_peaks_intersect,
                                                                 return_intersection_matrix = T,
                                                                 return_methylation_profile = T)

K562_MAFF_exclu_cofactor_matrix <- K562_MAFF_exclu_peaks_intersect_res$intersection_matrix
K562_MAFF_exclu_cofactor_matrix_t <- as.data.frame(t(K562_MAFF_exclu_cofactor_matrix))
attach(K562_MAFF_exclu_cofactor_matrix_t)
K562_MAFF_exclu_cofactor_matrix_order <- as.data.frame(K562_MAFF_exclu_cofactor_matrix_t[order(-MM1_HSA_K562_MAFF),,drop = FALSE])
detach(K562_MAFF_exclu_cofactor_matrix_t)

color <- colorRampPalette(c("white","#D46A6A", "#801515", "#550000"))
bk <- c(seq(0, 100, length=100))
pdf("cofactor_profile_in_K562_excluisve_MAFF_targets.pdf")
heatmap.2(as.matrix(cbind(K562_MAFF_exclu_cofactor_matrix_order,K562_MAFF_exclu_cofactor_matrix_order)),
          col=color(99), trace = "none", Colv = NULL, Rowv = NULL,
          dendrogram = "none", density.info = "none",
          key.xlab = "percentage (%)", keysize = 1.2, cexRow = .3,  
          key.title = "",labCol = NA, breaks = bk)
dev.off()

K562_MAFF_exclu_cofactor_meth <- K562_MAFF_exclu_peaks_intersect_res$methylation_profile_matrix[,c("MM1_HSA_K562_MAFF",
                                                                                                   "MM1_HSA_K562_C11orf30",
                                                                                                   "MM1_HSA_K562_NFE2",
                                                                                                   "MM1_HSA_K562_MAFK",
                                                                                                   "MM1_HSA_K562_ZNF316",
                                                                                                   "MM1_HSA_K562_NFE2L2",
                                                                                                   "MM1_HSA_K562_JUND",
                                                                                                   "MM1_HSA_K562_DPF2",
                                                                                                   "MM1_HSA_K562_BACH1",
                                                                                                   "MM1_HSA_K562_JUN",
                                                                                                   "MM1_HSA_K562_FOSL1",
                                                                                                   "MM1_HSA_K562_ATF3")]
K562_MAFF_exclu_cofactor_meth_matrix <- matrix(nrow = 12, ncol = 10)
for (i in 1:length(K562_MAFF_exclu_cofactor_meth)){
    meth_i <- K562_MAFF_exclu_cofactor_meth[i][[1]]
    K562_MAFF_exclu_cofactor_meth_matrix[i, ] <- 100*meth_i/sum(meth_i)
}
rownames(K562_MAFF_exclu_cofactor_meth_matrix) <- c("MAFF",
                                                    "C11orf30",
                                                    "NFE2",
                                                    "MAFK",
                                                    "ZNF316",
                                                    "NFE2L2",
                                                    "JUND",
                                                    "DPF2",
                                                    "BACH1",
                                                    "JUN",
                                                    "FOSL1",
                                                    "ATF3")
colnames(K562_MAFF_exclu_cofactor_meth_matrix) <- rownames(K562_MAFF_exclu_cofactor_meth[1][[1]])

color <- colorRampPalette(c("white", "#00B2D3", "#012694"))
bk <- c(seq(0, 100, length=100))
pdf("methylation_profile_of_cofactors_in_K562_exclusive_MAFF_targets.pdf")
heatmap.2(as.matrix(K562_MAFF_exclu_cofactor_meth_matrix),
          col=color(99), trace = "none", Colv = NULL, Rowv = NULL,
          dendrogram = "none", density.info = "none",
          key.xlab = "CpG percentage (%)", keysize = 1.2, cexRow = 1,  
          key.title = "",srtCol=45, breaks = bk)
dev.off()

# K562 exclusive MAFF targets without NFE2 and NFE2L2 co-binding
K562_MAFF_exclu_no_NFE <- exclusivePeaks(user_target_peak_list = list(K562_MAFF_exclu_peaks), 
                                         user_target_peak_id = "MM1_HSA_K562_MAFF",
                                         excluded_peak_id = c("MM1_HSA_K562_NFE2",
                                                              "MM1_HSA_K562_NFE2L2"), 
                                         motif_only_for_excluded_peak = T)
K562_MAFF_exclu_no_NFE_res <- exclusivePeakResult(exclusivePeaks = K562_MAFF_exclu_no_NFE, 
                                                  return_exclusive_peak_sites = T,
                                                  save_MethMotif_logo = T)
K562_MAFF_exclu_no_NFE_peaks <- K562_MAFF_exclu_no_NFE_res$exclusive_peak_list$MM1_HSA_K562_MAFF_exclusive_peaks

# cofactors in K562 exclusive MAFF targets without NFE2 and NFE2L2 co-binding
K562_MAFF_exclu_no_NFE_cofactor <- intersectPeakMatrix(user_peak_list_x = list(K562_MAFF_exclu_no_NFE_peaks),
                                                        user_peak_x_id = "MM1_HSA_K562_MAFF",
                                                        peak_id_y = K562_TFBS$ID,
                                                        motif_only_for_id_y = T, 
                                                        methylation_profile_in_narrow_region = T)

K562_MAFF_exclu_no_NFE_cofactor_res <- intersectPeakMatrixResult(intersectPeakMatrix = K562_MAFF_exclu_no_NFE_cofactor, 
                                                                  return_intersection_matrix = T, 
                                                                  return_methylation_profile = T)

K562_MAFF_exclu_no_NFE_cofactor_matrix <-K562_MAFF_exclu_no_NFE_cofactor_res$intersection_matrix

K562_MAFF_exclu_no_NFE_cofactor_matrix_t <- as.data.frame(t(K562_MAFF_exclu_no_NFE_cofactor_matrix))
attach(K562_MAFF_exclu_no_NFE_cofactor_matrix_t)
K562_MAFF_exclu_no_NFE_cofactor_matrix_order <- as.data.frame(K562_MAFF_exclu_no_NFE_cofactor_matrix_t[order(-MM1_HSA_K562_MAFF),,drop = FALSE])
detach(K562_MAFF_exclu_no_NFE_cofactor_matrix_t)

color <- colorRampPalette(c("white","#D46A6A", "#801515", "#550000"))
bk <- c(seq(0, 100, length=100))
pdf("cofactors_in_K562_exclusive_MAFF_targets_without_NFE2_and_NFE2L2.pdf")
heatmap.2(as.matrix(cbind(K562_MAFF_exclu_no_NFE_cofactor_matrix_order,K562_MAFF_exclu_no_NFE_cofactor_matrix_order)),
          col=color(99), trace = "none", Colv = NULL, Rowv = NULL,
          dendrogram = "none", density.info = "none",
          key.xlab = "percentage (%)", keysize = 1.2, cexRow = .3,  
          key.title = "",labCol = NA, breaks = bk)
dev.off()

K562_MAFF_exclu_no_NFE_cofactor_meth <- K562_MAFF_exclu_no_NFE_cofactor_res$methylation_profile_matrix[,c("MM1_HSA_K562_MAFF",
                                                                                                         "MM1_HSA_K562_ZNF316",
                                                                                                         "MM1_HSA_K562_C11orf30",
                                                                                                         "MM1_HSA_K562_MAFK")]
K562_MAFF_exclu_no_NFE_cofactor_meth_matrix <- matrix(nrow = 4, ncol = 10)
for (i in 1:length(K562_MAFF_exclu_no_NFE_cofactor_meth)){
    meth_i <- K562_MAFF_exclu_no_NFE_cofactor_meth[i][[1]]
    K562_MAFF_exclu_no_NFE_cofactor_meth_matrix[i, ] <- 100*meth_i/sum(meth_i)
}
rownames(K562_MAFF_exclu_no_NFE_cofactor_meth_matrix) <- c("MAFF","ZNF316","C11orf30","MAFK")
colnames(K562_MAFF_exclu_no_NFE_cofactor_meth_matrix) <- rownames(K562_MAFF_exclu_no_NFE_cofactor_meth[1][[1]])

color <- colorRampPalette(c("white", "#00B2D3", "#012694"))
bk <- c(seq(0, 100, length=100))
pdf("methylation_profile_of_cofactors_in_K562_exclusive_MAFF_targets_without_NFE2_and_NFE2L2_co-binding.pdf")
heatmap.2(as.matrix(K562_MAFF_exclu_no_NFE_cofactor_meth_matrix),
          col=color(99), trace = "none", Colv = NULL, Rowv = NULL,
          dendrogram = "none", density.info = "none",
          key.xlab = "CpG percentage (%)", keysize = 1.2, cexRow = 1,  
          key.title = "",srtCol=45, breaks = bk)
dev.off()
############### K562 exclusive MAFF targets (Figure 3B, C and D) #################################

############### HeLa-S3 exclusive MAFF targets (Figure 3B and C) #################################
HeLa_MAFF_exclu <- exclusivePeaks(target_peak_id = "MM1_HSA_HeLa-S3_MAFF",
                                  motif_only_for_target_peak = T,
                                  excluded_peak_id = c("MM1_HSA_K562_MAFF","MM1_HSA_HepG2_MAFF"),
                                  motif_only_for_excluded_peak = T)
HeLa_MAFF_exclu_res <- exclusivePeakResult(exclusivePeaks = HeLa_MAFF_exclu,
                                           return_exclusive_peak_sites = T,
                                           save_MethMotif_logo = T)
HeLa_MAFF_exclu_res_peaks <- HeLa_MAFF_exclu_res$exclusive_peak_list$`MM1_HSA_HeLa-S3_MAFF_exclusive_peaks`

# cofactors in HeLa-S3 exclusive MAFF targets
HeLa_TFBS <- TFBSBrowser(cell_tissue_name = "HeLa-S3")
HeLa_MAFF_exclu_peaks_cofactor <- intersectPeakMatrix(user_peak_list_x = list(HeLa_MAFF_exclu_res_peaks),
                                                       user_peak_x_id = "MM1_HSA_HeLa-S3_MAFF",
                                                       peak_id_y = HeLa_TFBS$ID,
                                                       motif_only_for_id_y = T,
                                                       methylation_profile_in_narrow_region = T)
HeLa_MAFF_exclu_peaks_cofactor_res <- intersectPeakMatrixResult(intersectPeakMatrix = HeLa_MAFF_exclu_peaks_cofactor,
                                                                 return_intersection_matrix = T,
                                                                 return_methylation_profile = T)

HeLa_MAFF_exclu_peaks_cofactor_matrix <- HeLa_MAFF_exclu_peaks_cofactor_res$intersection_matrix
HeLa_MAFF_exclu_peaks_cofactor_matrix_t <- as.data.frame(t(HeLa_MAFF_exclu_peaks_cofactor_matrix))
attach(HeLa_MAFF_exclu_peaks_cofactor_matrix_t)
HeLa_MAFF_exclu_peaks_cofactor_matrix_order <- as.data.frame(HeLa_MAFF_exclu_peaks_cofactor_matrix_t[order(-`MM1_HSA_HeLa-S3_MAFF`),,drop = FALSE])
detach(HeLa_MAFF_exclu_peaks_cofactor_matrix_t)

color <- colorRampPalette(c("white","#D46A6A", "#801515", "#550000"))
bk <- c(seq(0, 100, length=100))
pdf("cofactors_in_HeLa-S3_exclusive_MAFF_targets.pdf")
heatmap.2(as.matrix(cbind(HeLa_MAFF_exclu_peaks_cofactor_matrix_order,HeLa_MAFF_exclu_peaks_cofactor_matrix_order)),
          col=color(99), trace = "none", Colv = NULL, Rowv = NULL,
          dendrogram = "none", density.info = "none",
          key.xlab = "percentage (%)", keysize = 1.2, cexRow = .3,  
          key.title = "",labCol = NA, breaks = bk)
dev.off()

HeLa_MAFF_exclu_cofactor_meth <- HeLa_MAFF_exclu_peaks_cofactor_res$methylation_profile_matrix[,c("MM1_HSA_HeLa-S3_MAFF",
                                                                                                  "MM1_HSA_HeLa-S3_MAFK",
                                                                                                 "MM1_HSA_HeLa-S3_JUND",
                                                                                                 "MM1_HSA_HeLa-S3_FOS",
                                                                                                 "MM1_HSA_HeLa-S3_NFE2L2",
                                                                                                 "MM1_HSA_HeLa-S3_CEBPB")]
HeLa_MAFF_exclu_cofactor_meth_matrix <- matrix(nrow = 6, ncol = 10)
for (i in 1:length(HeLa_MAFF_exclu_cofactor_meth)){
    meth_i <- HeLa_MAFF_exclu_cofactor_meth[i][[1]]
    HeLa_MAFF_exclu_cofactor_meth_matrix[i, ] <- 100*meth_i/sum(meth_i)
}
rownames(HeLa_MAFF_exclu_cofactor_meth_matrix) <- c("MAFF",
                                                    "MAFK",
                                                  "JUND",
                                                  "FOS",
                                                  "NFE2L2",
                                                  "CEBPB")
colnames(HeLa_MAFF_exclu_cofactor_meth_matrix) <- rownames(HeLa_MAFF_exclu_cofactor_meth[1][[1]])

color <- colorRampPalette(c("white", "#00B2D3", "#012694"))
bk <- c(seq(0, 100, length=100))
pdf("methylation_profile_of_cofactors_in_HeLa-S3_exclusive_MAFF_targets.pdf")
heatmap.2(as.matrix(HeLa_MAFF_exclu_cofactor_meth_matrix),
          col=color(99), trace = "none", Colv = NULL, Rowv = NULL,
          dendrogram = "none", density.info = "none",
          key.xlab = "CpG percentage (%)", keysize = 1.2, cexRow = 1,  
          key.title = "",srtCol=45, breaks = bk)
dev.off()

############### HeLa-S3 exclusive MAFF targets (Figure 3B and C) #################################


############### HepG2 exclusive MAFF targets (Figure 3B and C) #################################

HepG2_MAFF_exclu <- exclusivePeaks(target_peak_id = "MM1_HSA_HepG2_MAFF",
                                   motif_only_for_target_peak = T,
                                   excluded_peak_id = c("MM1_HSA_K562_MAFF",
                                                        "MM1_HSA_HeLa-S3_MAFF"),
                                   motif_only_for_excluded_peak = T)
HepG2_MAFF_exclu_res <- exclusivePeakResult(exclusivePeaks = HepG2_MAFF_exclu,
                                            return_exclusive_peak_sites = T,
                                            save_MethMotif_logo = T)
HepG2_MAFF_exclu_peaks <- HepG2_MAFF_exclu_res$exclusive_peak_list$MM1_HSA_HepG2_MAFF_exclusive_peaks

# cofactors in HepG2 exclusive MAFF targets
HepG2_TFBS <- TFBSBrowser(cell_tissue_name = "HepG2")
HepG2_MAFF_exclu_peaks_cofactor <- intersectPeakMatrix(user_peak_list_x = list(HepG2_MAFF_exclu_peaks),
                                                        user_peak_x_id = "MM1_HSA_HepG2_MAFF",
                                                        peak_id_y = HepG2_TFBS$ID,
                                                        motif_only_for_id_y = T,
                                                        methylation_profile_in_narrow_region = T)
HepG2_MAFF_exclu_peaks_cofactor_res <- intersectPeakMatrixResult(intersectPeakMatrix = HepG2_MAFF_exclu_peaks_cofactor,
                                                                  return_intersection_matrix = T,
                                                                  return_methylation_profile = T)

HepG2_MAFF_exclu_peaks_cofactor_matrix <- HepG2_MAFF_exclu_peaks_cofactor_res$intersection_matrix
HepG2_MAFF_exclu_peaks_cofactor_matrix_t <- as.data.frame(t(HepG2_MAFF_exclu_peaks_cofactor_matrix))
attach(HepG2_MAFF_exclu_peaks_cofactor_matrix_t)
HepG2_MAFF_exclu_peaks_cofactor_matrix_order <- as.data.frame(HepG2_MAFF_exclu_peaks_cofactor_matrix_t[order(-MM1_HSA_HepG2_MAFF),,drop = FALSE])
detach(HepG2_MAFF_exclu_peaks_cofactor_matrix_t)

color <- colorRampPalette(c("white","#D46A6A", "#801515", "#550000"))
bk <- c(seq(0, 100, length=100))
pdf("cofactors_in_HepG2_exclusive_MAFF_targets.pdf")
heatmap.2(as.matrix(cbind(HepG2_MAFF_exclu_peaks_cofactor_matrix_order,HepG2_MAFF_exclu_peaks_cofactor_matrix_order)),
          col=color(99), trace = "none", Colv = NULL, Rowv = NULL,
          dendrogram = "none", density.info = "none",
          key.xlab = "percentage (%)", keysize = 1.2, cexRow = .3,  
          key.title = "",labCol = NA, breaks = bk)
dev.off()

HepG2_MAFF_exclu_cofactor_meth <- HepG2_MAFF_exclu_peaks_cofactor_res$methylation_profile_matrix[,c("MM1_HSA_HepG2_MAFF",
                                                                                                    "MM1_HSA_HepG2_MAFK")]
HepG2_MAFF_exclu_cofactor_meth_matrix <- matrix(nrow = 2, ncol = 10)
for (i in 1:length(HepG2_MAFF_exclu_cofactor_meth)){
    meth_i <- HepG2_MAFF_exclu_cofactor_meth[i][[1]]
    HepG2_MAFF_exclu_cofactor_meth_matrix[i, ] <- 100*meth_i/sum(meth_i)
}
rownames(HepG2_MAFF_exclu_cofactor_meth_matrix) <- c("MAFF","MAFK")
colnames(HepG2_MAFF_exclu_cofactor_meth_matrix) <- rownames(HepG2_MAFF_exclu_cofactor_meth[1][[1]])

color = colorRampPalette(c("white", "#00B2D3", "#012694"))
bk <- c(seq(0, 100, length=100))
pdf("methylation_profile_of_cofactors_in_HepG2_exclusive_MAFF_targets.pdf")
heatmap.2(as.matrix(HepG2_MAFF_exclu_cofactor_meth_matrix),
          col=color(99), trace = "none", Colv = NULL, Rowv = NULL,
          dendrogram = "none", density.info = "none",
          key.xlab = "CpG percentage (%)", keysize = 1.2, cexRow = 1,  
          key.title = "",srtCol=45, breaks = bk)
dev.off()

############### HepG2 exclusive MAFF targets (Figure 3B and C) #################################

