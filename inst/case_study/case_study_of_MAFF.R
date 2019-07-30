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
                                  motif_only_for_target_peak = TRUE,
                                  compared_peak_id = MAFF_TFBS$ID,
                                  motif_only_for_compared_peak = TRUE)
MAFF_common_result <- commonPeakResult(commonPeaks = MAFF_common_output,
                                       return_common_peak_sites = TRUE)
K562_MAFF_commonPeaks <- MAFF_common_result$common_peak_list$MM1_HSA_K562_MAFF_common_peaks
HeLa_MAFF_commonPeaks <- MAFF_common_result$common_peak_list$`MM1_HSA_HeLa-S3_MAFF_common_peaks`
HepG2_MAFF_commonPeaks <- MAFF_common_result$common_peak_list$MM1_HSA_HepG2_MAFF_common_peaks

############### K562 exclusive MAFF targets (Figure 3B, C and D) #################################
K562_MAFF_exclu <- exclusivePeaks(target_peak_id = "MM1_HSA_K562_MAFF",
                                  motif_only_for_target_peak = TRUE,
                                  excluded_peak_id = c("MM1_HSA_HeLa-S3_MAFF",
                                                       "MM1_HSA_HepG2_MAFF"),
                                  motif_only_for_excluded_peak = TRUE)
K562_MAFF_exclu_res <- exclusivePeakResult(exclusivePeaks = K562_MAFF_exclu,
                                           return_exclusive_peak_sites = TRUE,
                                           save_MethMotif_logo = TRUE)
K562_MAFF_exclu_peaks <- K562_MAFF_exclu_res$exclusive_peak_list$MM1_HSA_K562_MAFF_exclusive_peaks

# cofactors in K562 exclusive MAFF targets
K562_TFBS <- TFBSBrowser(cell_tissue_name = "K562")
K562_MAFF_exclu_peaks_intersect <- intersectPeakMatrix(user_peak_list_x = list(K562_MAFF_exclu_peaks),
                                                       user_peak_x_id = "MM1_HSA_K562_MAFF",
                                                       peak_id_y = K562_TFBS$ID,
                                                       motif_only_for_id_y = TRUE,
                                                       methylation_profile_in_narrow_region = TRUE)

# simply using cofactorReport() can report the top cofactors along with motif logos, DNA methylation and read enrichments
cofactorReport(intersectPeakMatrix = K562_MAFF_exclu_peaks_intersect, cobinding_threshold = 0.1)

# K562 exclusive MAFF targets without NFE2 and NFE2L2 co-binding
K562_MAFF_exclu_no_NFE <- exclusivePeaks(user_target_peak_list = list(K562_MAFF_exclu_peaks), 
                                         user_target_peak_id = "MM1_HSA_K562_MAFF",
                                         excluded_peak_id = c("MM1_HSA_K562_NFE2",
                                                              "MM1_HSA_K562_NFE2L2"), 
                                         motif_only_for_excluded_peak = TRUE)
K562_MAFF_exclu_no_NFE_res <- exclusivePeakResult(exclusivePeaks = K562_MAFF_exclu_no_NFE, 
                                                  return_exclusive_peak_sites = TRUE,
                                                  save_MethMotif_logo = TRUE)
K562_MAFF_exclu_no_NFE_peaks <- K562_MAFF_exclu_no_NFE_res$exclusive_peak_list$MM1_HSA_K562_MAFF_exclusive_peaks

# cofactors in K562 exclusive MAFF targets without NFE2 and NFE2L2 co-binding
K562_MAFF_exclu_no_NFE_cofactor <- intersectPeakMatrix(user_peak_list_x = list(K562_MAFF_exclu_no_NFE_peaks),
                                                        user_peak_x_id = "MM1_HSA_K562_MAFF",
                                                        peak_id_y = K562_TFBS$ID,
                                                        motif_only_for_id_y = TRUE, 
                                                        methylation_profile_in_narrow_region = TRUE)
cofactorReport(K562_MAFF_exclu_no_NFE_cofactor, cobinding_threshold = 0.1)

############### K562 exclusive MAFF targets (Figure 3B, C and D) #################################

############### HeLa-S3 exclusive MAFF targets (Figure 3B and C) #################################
HeLa_MAFF_exclu <- exclusivePeaks(target_peak_id = "MM1_HSA_HeLa-S3_MAFF",
                                  motif_only_for_target_peak = TRUE,
                                  excluded_peak_id = c("MM1_HSA_K562_MAFF","MM1_HSA_HepG2_MAFF"),
                                  motif_only_for_excluded_peak = TRUE)
HeLa_MAFF_exclu_res <- exclusivePeakResult(exclusivePeaks = HeLa_MAFF_exclu,
                                           return_exclusive_peak_sites = TRUE,
                                           save_MethMotif_logo = TRUE)
HeLa_MAFF_exclu_res_peaks <- HeLa_MAFF_exclu_res$exclusive_peak_list$`MM1_HSA_HeLa-S3_MAFF_exclusive_peaks`

# cofactors in HeLa-S3 exclusive MAFF targets
HeLa_TFBS <- TFBSBrowser(cell_tissue_name = "HeLa-S3")
HeLa_MAFF_exclu_peaks_cofactor <- intersectPeakMatrix(user_peak_list_x = list(HeLa_MAFF_exclu_res_peaks),
                                                       user_peak_x_id = "MM1_HSA_HeLa-S3_MAFF",
                                                       peak_id_y = HeLa_TFBS$ID,
                                                       motif_only_for_id_y = TRUE,
                                                       methylation_profile_in_narrow_region = TRUE)

cofactorReport(intersectPeakMatrix = HeLa_MAFF_exclu_peaks_cofactor, cobinding_threshold = 0.1)

############### HeLa-S3 exclusive MAFF targets (Figure 3B and C) #################################


############### HepG2 exclusive MAFF targets (Figure 3B and C) #################################

HepG2_MAFF_exclu <- exclusivePeaks(target_peak_id = "MM1_HSA_HepG2_MAFF",
                                   motif_only_for_target_peak = TRUE,
                                   excluded_peak_id = c("MM1_HSA_K562_MAFF",
                                                        "MM1_HSA_HeLa-S3_MAFF"),
                                   motif_only_for_excluded_peak = TRUE)
HepG2_MAFF_exclu_res <- exclusivePeakResult(exclusivePeaks = HepG2_MAFF_exclu,
                                            return_exclusive_peak_sites = TRUE,
                                            save_MethMotif_logo = TRUE)
HepG2_MAFF_exclu_peaks <- HepG2_MAFF_exclu_res$exclusive_peak_list$MM1_HSA_HepG2_MAFF_exclusive_peaks

# cofactors in HepG2 exclusive MAFF targets
HepG2_TFBS <- TFBSBrowser(cell_tissue_name = "HepG2")
HepG2_MAFF_exclu_peaks_cofactor <- intersectPeakMatrix(user_peak_list_x = list(HepG2_MAFF_exclu_peaks),
                                                        user_peak_x_id = "MM1_HSA_HepG2_MAFF",
                                                        peak_id_y = HepG2_TFBS$ID,
                                                        motif_only_for_id_y = TRUE,
                                                        methylation_profile_in_narrow_region = TRUE)

cofactorReport(intersectPeakMatrix = HepG2_MAFF_exclu_peaks_cofactor, cobinding_threshold = 0.1)
HepG2_MAFF_exclu_peaks_cofactor_res <- intersectPeakMatrixResult(intersectPeakMatrix = HepG2_MAFF_exclu_peaks_cofactor,
                                                                  return_intersection_matrix = TRUE,
                                                                  return_methylation_profile = TRUE)
############### HepG2 exclusive MAFF targets (Figure 3B and C) #################################

