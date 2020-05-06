library(TFregulomeR)
library(gplots)

# CEBPB motifs in TFregulomeR
CEBPB_record <- dataBrowser(tf = "CEBPB", species = "human")

# (Meth)Motif logo for all CEBPB in TFregulomeR compendium for Supplementary Figure 1
for (id in CEBPB_record$ID){
    motif_matrix <- searchMotif(id = id)
    plotLogo(MM_object = motif_matrix)
}

# segregate K562 CEBPB binding sites into 16 sub-ensembles according to the number
# of cell types where they are enriched
# 1) load K562 CEBPB peaks
K562_CEBPB_all_peaks <- loadPeaks(id = "MM1_HSA_K562_CEBPB",
                                  includeMotifOnly = TRUE)
# 2) segregate them into 16 sub-ensembles
for (i in seq(1,16,1)){
    cell_i <- CEBPB_record$cell_tissue_name[i]
    id_i <- CEBPB_record$ID[i]
    common_peak_i <- commonPeaks(target_peak_id = "MM1_HSA_K562_CEBPB",
                                 motif_only_for_target_peak = TRUE,
                                 compared_peak_id = id_i,
                                 motif_only_for_compared_peak = TRUE)
    common_peak_res_i <- commonPeakResult(commonPeaks = common_peak_i,
                                          return_common_peak_sites = TRUE)
    common_peak_regions_i <- common_peak_res_i$common_peak_list$MM1_HSA_K562_CEBPB_common_peaks
    K562_CEBPB_all_peaks[,cell_i] <- 0
    K562_CEBPB_all_peaks[(K562_CEBPB_all_peaks$id %in% common_peak_regions_i$id),cell_i] <- 1
}
K562_CEBPB_all_peaks$sum <- apply(K562_CEBPB_all_peaks[,CEBPB_record$cell_tissue_name],1,sum)

# Peak numbers, MethMotif logs and read enrichments for 16 K562 CEBPB sub-ensembles for Supplementary Figure 2 and Figure 1A
sub_ensemble_peak_num <- c()
sub_ensemble_read_score <- as.data.frame(matrix(nrow = nrow(K562_CEBPB_all_peaks),
                                                ncol = 2))
colnames(sub_ensemble_read_score) <- c("cell_type_num", "read_fold_change")
sub_ensemble_read_score$cell_type_num <- K562_CEBPB_all_peaks$sum
for (i in seq(1,16,1)){
    peak_subset_i <- K562_CEBPB_all_peaks[which(K562_CEBPB_all_peaks$sum==i),]
    sub_ensemble_peak_num <- c(sub_ensemble_peak_num, nrow(peak_subset_i))
    sub_ensemble_read_score[which(sub_ensemble_read_score$cell_type_num==i),'read_fold_change'] <-
        peak_subset_i$tag_fold_change

    # MethMotif logo
    common_peak_i <- commonPeaks(user_target_peak_list = list(peak_subset_i[,1:5]),
                                 user_target_peak_id = "MM1_HSA_K562_CEBPB",
                                 compared_peak_id = "MM1_HSA_K562_CEBPB",
                                 motif_only_for_compared_peak = TRUE)
    common_peak_res_i <- commonPeakResult(commonPeaks = common_peak_i,
                                          save_MethMotif_logo = TRUE)
    file.rename("MM1_HSA_K562_CEBPB_common_peaks-logo-entropy.pdf",
                paste0("K562_CEBPB_MethMotif_in_peaks_shared_by_",i,"_cell_types.pdf"))
}
sub_ensemble_read_score$cell_type_num <- factor(sub_ensemble_read_score$cell_type_num,
                                                levels = seq(1,16,1))
pdf("read_enrichment_scores_across_16_K562_CEBPB_sub-ensembles.pdf")
boxplot(read_fold_change~cell_type_num,sub_ensemble_read_score,
        xlab = "number of cell types sharing the peaks",
        ylab = "Read enrichment score", outline=FALSE)
dev.off()

pdf("peak_numbers_across_16_K562_CEBPB_sub-ensembles.pdf")
plot(x=seq(1,16,1), y = sub_ensemble_peak_num, type = "l",
     xlab = "number of cell types sharing the peaks",
     ylab = "number of peaks")
dev.off()

# cofactor analysis of 16 K562 CEBPB sub-ensembles for Figure 2B
K562_TFBS <- dataBrowser(cell_tissue_name = "K562")
K562_CEBPB_16_subsets_list <- list()
for (i in seq(1,16,1)){
    K562_CEBPB_16_subsets_i <- K562_CEBPB_all_peaks[which(K562_CEBPB_all_peaks$sum==i),]
    K562_CEBPB_16_subsets_list[[i]] <- K562_CEBPB_16_subsets_i
}
cofactor_16_subsets <- intersectPeakMatrix(user_peak_list_x = K562_CEBPB_16_subsets_list,
                                           peak_id_y = K562_TFBS$ID,
                                           motif_only_for_id_y = TRUE)
cofactor_16_subsets_res <- intersectPeakMatrixResult(intersectPeakMatrix = cofactor_16_subsets,
                                                     return_intersection_matrix = TRUE,
                                                     angle_of_matrix = "x")
cofactor_16_subsets_matrix <- cofactor_16_subsets_res$intersection_matrix
cofactor_16_subsets_matrix_t <- as.data.frame(t(cofactor_16_subsets_matrix))
# filter out cofactor whose binding percents are less than 5 in all sub-ensembles
cofactor_16_subsets_matrix_filtered <- cofactor_16_subsets_matrix_t[!(cofactor_16_subsets_matrix_t$user_peak_x1<=5 &
                                                                         cofactor_16_subsets_matrix_t$user_peak_x2<=5 &
                                                                         cofactor_16_subsets_matrix_t$user_peak_x3<=5 &
                                                                         cofactor_16_subsets_matrix_t$user_peak_x4<=5 &
                                                                         cofactor_16_subsets_matrix_t$user_peak_x5<=5 &
                                                                         cofactor_16_subsets_matrix_t$user_peak_x6<=5 &
                                                                         cofactor_16_subsets_matrix_t$user_peak_x7<=5 &
                                                                         cofactor_16_subsets_matrix_t$user_peak_x8<=5 &
                                                                         cofactor_16_subsets_matrix_t$user_peak_x9<=5 &
                                                                         cofactor_16_subsets_matrix_t$user_peak_x10<=5 &
                                                                         cofactor_16_subsets_matrix_t$user_peak_x11<=5 &
                                                                         cofactor_16_subsets_matrix_t$user_peak_x12<=5 &
                                                                         cofactor_16_subsets_matrix_t$user_peak_x13<=5 &
                                                                         cofactor_16_subsets_matrix_t$user_peak_x14<=5 &
                                                                         cofactor_16_subsets_matrix_t$user_peak_x15<=5 &
                                                                         cofactor_16_subsets_matrix_t$user_peak_x16<=5),]
color <- colorRampPalette(c("white","#D46A6A", "#801515", "#550000"))
bk <- c(seq(0, 100, length=100))
pdf("cofactors_in_16_K562_CEBPB_sub-ensembles.pdf")
heatmap.2(as.matrix(cofactor_16_subsets_matrix_filtered), dendrogram = "row",
          col=color(99), trace = "none", Colv = NULL, density.info = "none",
          key.xlab = "percentage (%)", keysize = 1.2, cexRow = .3,
          key.title = "",labCol = NA, breaks = bk)
dev.off()

# mCG percentage and read enrichment of CEBPB-CEBPD, CEBPB-ATF4 combinations across 16 sub-ensembles for Figure 2C
CEBPB_CEBPD_meth_value <- c()
CEBPB_ATF4_meth_value <- c()
CEBPB_CEBPD_tag_median <- c()
CEBPB_ATF4_tag_median <- c()

CEBPB_CEBPD_tag_q1 <- c()
CEBPB_ATF4_tag_q1 <- c()

CEBPB_CEBPD_tag_q3 <- c()
CEBPB_ATF4_tag_q3 <- c()
for (i in seq(1,16,1)){
    K562_CEBPB_16_subsets_i <- K562_CEBPB_16_subsets_list[[i]]
    cofactor_in_subsets_i <- intersectPeakMatrix(user_peak_list_x = list(K562_CEBPB_16_subsets_i),
                                                 user_peak_x_id = "MM1_HSA_K562_CEBPB",
                                                 peak_id_y = c("MM1_HSA_K562_CEBPD",
                                                               "MM1_HSA_K562_ATF4"),
                                                 motif_only_for_id_y = TRUE,
                                                 methylation_profile_in_narrow_region = TRUE)
    cofactor_in_subsets_res_i <- intersectPeakMatrixResult(intersectPeakMatrix = cofactor_in_subsets_i,
                                                           return_tag_density = TRUE,
                                                           angle_of_tag_density = "x",
                                                           tag_density_value = "median",
                                                           return_methylation_profile = TRUE,
                                                           angle_of_methylation_profile = "x")

    meth_matrix_i <- cofactor_in_subsets_res_i$methylation_profile_matrix
    CEBPB_CEBPD_meth_i <- meth_matrix_i["MM1_HSA_K562_CEBPB","MM1_HSA_K562_CEBPD"][[1]]
    CEBPB_CEBPD_meth_value_i <- sum(CEBPB_CEBPD_meth_i[9:10])*100/sum(CEBPB_CEBPD_meth_i)
    CEBPB_CEBPD_meth_value <- c(CEBPB_CEBPD_meth_value, CEBPB_CEBPD_meth_value_i)

    CEBPB_ATF4_meth_i <- meth_matrix_i["MM1_HSA_K562_CEBPB","MM1_HSA_K562_ATF4"][[1]]
    CEBPB_ATF4_meth_value_i <- sum(CEBPB_ATF4_meth_i[9:10])*100/sum(CEBPB_ATF4_meth_i)
    CEBPB_ATF4_meth_value <- c(CEBPB_ATF4_meth_value, CEBPB_ATF4_meth_value_i)

    tag_density_median_i <- cofactor_in_subsets_res_i$tag_density_matrix
    CEBPB_CEBPD_tag_median <- c(CEBPB_CEBPD_tag_median,
                                tag_density_median_i["MM1_HSA_K562_CEBPB",
                                                     "MM1_HSA_K562_CEBPD"])
    CEBPB_ATF4_tag_median <- c(CEBPB_ATF4_tag_median,
                               tag_density_median_i["MM1_HSA_K562_CEBPB",
                                                    "MM1_HSA_K562_ATF4"])

    cofactor_tag_q1_i <- intersectPeakMatrixResult(intersectPeakMatrix = cofactor_in_subsets_i,
                                                   return_tag_density = TRUE,
                                                   angle_of_tag_density = "x",
                                                   tag_density_value = "quartile_25")
    tag_density_q1_i <- cofactor_tag_q1_i$tag_density_matrix
    CEBPB_CEBPD_tag_q1 <- c(CEBPB_CEBPD_tag_q1,
                            tag_density_q1_i["MM1_HSA_K562_CEBPB",
                                             "MM1_HSA_K562_CEBPD"])
    CEBPB_ATF4_tag_q1 <- c(CEBPB_ATF4_tag_q1,
                           tag_density_q1_i["MM1_HSA_K562_CEBPB",
                                            "MM1_HSA_K562_ATF4"])
    cofactor_tag_q3_i <- intersectPeakMatrixResult(intersectPeakMatrix = cofactor_in_subsets_i,
                                                   return_tag_density = TRUE,
                                                   tag_density_value = "quartile_75",
                                                   angle_of_tag_density = "x")
    tag_density_q3_i <- cofactor_tag_q3_i$tag_density_matrix
    CEBPB_CEBPD_tag_q3 <- c(CEBPB_CEBPD_tag_q3,
                            tag_density_q3_i["MM1_HSA_K562_CEBPB",
                                             "MM1_HSA_K562_CEBPD"])
    CEBPB_ATF4_tag_q3 <- c(CEBPB_ATF4_tag_q3,
                           tag_density_q3_i["MM1_HSA_K562_CEBPB",
                                            "MM1_HSA_K562_ATF4"])
}
pdf("mCG_percentage_of_CEBPB-CEPBD_and_CEBPB-ATF4_in_16_subsets.pdf")
plot(x = seq(1,16,1),
     y = CEBPB_CEBPD_meth_value,
     type="l", ylim = c(0,30), xlim = c(0,17), ylab = "5mC percentage (%)",
     xlab = "number of shared cell types", col="blue")
lines(x = seq(1,16,1),
      y = CEBPB_ATF4_meth_value, col="red")
dev.off()

pdf("read_enrichments_of_CEBPB-CEPBD_and_CEBPB-ATF4_in_16_subsets.pdf")
plot(x = seq(1,16,1),
     y = CEBPB_CEBPD_tag_median,
     type="l", ylim = c(0,100), xlim = c(0,17), ylab = "tag density median",
     xlab = "number of shared cell types", col="blue")
points(x = seq(1,16,1),
       y = CEBPB_CEBPD_tag_median, pch=20, col="blue")
points(x = seq(1,16,1),
       y = CEBPB_CEBPD_tag_q1, pch="-", col="blue")
points(x = seq(1,16,1),
       y = CEBPB_CEBPD_tag_q3, pch="-", col="blue")
segments(x0=seq(1,16,1), x1=seq(1,16,1),
         y0=CEBPB_CEBPD_tag_q1, y1=CEBPB_CEBPD_tag_q3, col="blue")
lines(x = seq(1,16,1),
      y = CEBPB_ATF4_tag_median, col="red")
points(x = seq(1,16,1),
       y = CEBPB_ATF4_tag_median, pch=20,col="red")
points(x = seq(1,16,1),
       y = CEBPB_ATF4_tag_q1, pch="-", col="red")
points(x = seq(1,16,1),
       y = CEBPB_ATF4_tag_q3, pch="-", col="red")
segments(x0=seq(1,16,1), x1=seq(1,16,1),
         y0=CEBPB_ATF4_tag_q1, y1=CEBPB_ATF4_tag_q3, col="red")
dev.off()

# MethMotif logos in K562 exclusive CEBPB peaks with and without ATF4 for Figure 2D
#remove MM1_HSA_K562_CEBPB ID from all CEBPB TFregulomeR IDs
CEBPB_record_ID_noK562 <- CEBPB_record$ID[!(CEBPB_record$ID %in% "MM1_HSA_K562_CEBPB")]

K562_exclusivePeak_output <- exclusivePeaks(target_peak_id = "MM1_HSA_K562_CEBPB",
                                            motif_only_for_target_peak = TRUE,
                                            excluded_peak_id = CEBPB_record_ID_noK562,
                                            motif_only_for_excluded_peak = TRUE)
K562_exclusivePeak_result <- exclusivePeakResult(exclusivePeaks = K562_exclusivePeak_output,
                                                 return_exclusive_peak_sites = TRUE)
K562_exclusivePeak_peak <- K562_exclusivePeak_result$exclusive_peak_list$MM1_HSA_K562_CEBPB_exclusive_peaks

K562_exclusivePeak_with_ATF4_output <- commonPeaks(user_target_peak_list = list(K562_exclusivePeak_peak),
                                                   user_target_peak_id = "MM1_HSA_K562_CEBPB",
                                                   compared_peak_id = "MM1_HSA_K562_ATF4",
                                                   motif_only_for_compared_peak = TRUE)
K562_exclusivePeak_with_ATF4_res <- commonPeakResult(commonPeaks = K562_exclusivePeak_with_ATF4_output,
                                                     return_common_peak_sites = TRUE,
                                                     save_MethMotif_logo = TRUE)
K562_exclusivePeak_with_ATF4_peaks <- K562_exclusivePeak_with_ATF4_res$common_peak_list$MM1_HSA_K562_CEBPB_common_peaks
K562_exclusivePeak_without_ATF4_peaks <- K562_exclusivePeak_peak[!(K562_exclusivePeak_peak$id %in% K562_exclusivePeak_with_ATF4_peaks$id),]

K562_exclusivePeak_without_ATF4_output <- commonPeaks(user_target_peak_list = list(K562_exclusivePeak_without_ATF4_peaks),
                                                      user_target_peak_id = "MM1_HSA_K562_CEBPB",
                                                      compared_peak_id = "MM1_HSA_K562_CEBPB",
                                                      motif_only_for_compared_peak = TRUE)
commonPeakResult(commonPeaks = K562_exclusivePeak_without_ATF4_output,
                 save_MethMotif_logo = TRUE)

# MethMotif logos in K562 shared CEBPB peaks with and without CEBPD for Figure 2D
K562_commonPeak_output <- commonPeaks(target_peak_id = "MM1_HSA_K562_CEBPB",
                                      motif_only_for_target_peak = TRUE,
                                      compared_peak_id = CEBPB_record$ID,
                                      motif_only_for_compared_peak = TRUE)
K562_commonPeak_result <- commonPeakResult(commonPeaks = K562_commonPeak_output,
                                           return_common_peak_sites = TRUE)
K562_commonPeak_peak <- K562_commonPeak_result$common_peak_list$MM1_HSA_K562_CEBPB_common_peaks

K562_commonPeak_with_CEBPB_output <- commonPeaks(user_target_peak_list = list(K562_commonPeak_peak),
                                                 user_target_peak_id = "MM1_HSA_K562_CEBPB",
                                                 compared_peak_id = "MM1_HSA_K562_CEBPD",
                                                 motif_only_for_compared_peak = TRUE)
K562_commonPeak_with_CEBPB_res <- commonPeakResult(commonPeaks = K562_commonPeak_with_CEBPB_output,
                                                   return_common_peak_sites = TRUE,
                                                   save_MethMotif_logo = TRUE)
K562_commonPeak_with_CEBPB_peaks <- K562_commonPeak_with_CEBPB_res$common_peak_list$MM1_HSA_K562_CEBPB_common_peaks
K562_commonPeak_without_CEBPB_peaks <- K562_commonPeak_peak[!(K562_commonPeak_peak$id %in% K562_commonPeak_with_CEBPB_peaks$id),]

K562_commonPeak_without_ATF4_output <- commonPeaks(user_target_peak_list = list(K562_commonPeak_without_CEBPB_peaks),
                                                   user_target_peak_id = "MM1_HSA_K562_CEBPB",
                                                   compared_peak_id = "MM1_HSA_K562_CEBPB",
                                                   motif_only_for_compared_peak = TRUE)
commonPeakResult(commonPeaks = K562_commonPeak_without_ATF4_output,
                 save_MethMotif_logo = TRUE)


# motif in shared and exclusive CEBPB targets in all other 15 cell types (Supplementary Figure 3)
for (i in CEBPB_record$ID){
    common_i <- commonPeaks(target_peak_id = i,
                            motif_only_for_target_peak = TRUE,
                            compared_peak_id = CEBPB_record$ID,
                            motif_only_for_compared_peak = TRUE)
    common_i_res <- commonPeakResult(commonPeaks = common_i,
                                     save_MethMotif_logo = TRUE)

    # exclude ID i from all CEBPB IDs
    cebpd_id_no_i <- CEBPB_record$ID[!(CEBPB_record$ID %in% i)]
    exclusive_i <- exclusivePeaks(target_peak_id = i,
                                  motif_only_for_target_peak = TRUE,
                                  excluded_peak_id = cebpd_id_no_i,
                                  motif_only_for_excluded_peak = TRUE)
    exclusive_i_res <- exclusivePeakResult(exclusivePeaks = exclusive_i,
                                           save_MethMotif_logo = TRUE)
}


# functions of CEBPB/CEBPD and CEBPB/ATF4 targets in K562 (Supplementary Figure 4)
# load required package for GREAT annotation
library(rGREAT)
# load required package for genomic conversion from hg38 to hg19
library(liftOver)
# all CEBPB-CEBPD co-binding regions in K562
K562_CEBPB_CEBPD <- commonPeaks(target_peak_id = "MM1_HSA_K562_CEBPB",
                                motif_only_for_target_peak = TRUE,
                                compared_peak_id = "MM1_HSA_K562_CEBPD",
                                motif_only_for_compared_peak = TRUE)
K562_CEBPB_CEBPD_res <- commonPeakResult(commonPeaks = K562_CEBPB_CEBPD,
                                         return_common_peak_sites = TRUE,
                                         save_MethMotif_logo = TRUE,
                                         return_summary = TRUE)
K562_CEBPB_CEBPD_res$peak_summary
#>                                 percentage_in_original_inputs(%)
#> MM1_HSA_K562_CEBPB_common_peaks                         6.532727
K562_CEBPB_CEBPD_peak <- K562_CEBPB_CEBPD_res$common_peak_list$MM1_HSA_K562_CEBPB_common_peaks
K562_CEBPB_CEBPD_great <- greatAnnotate(peaks = K562_CEBPB_CEBPD_peak,
                                        return_annotation = TRUE)
K562_CEBPB_CEBPD_great_bp <- K562_CEBPB_CEBPD_great[which(K562_CEBPB_CEBPD_great$category=="BP"),]

# all CEBPB-ATF4 co-binding regions in K562
K562_CEBPB_ATF4 <- commonPeaks(target_peak_id = "MM1_HSA_K562_CEBPB",
                               motif_only_for_target_peak = TRUE,
                               compared_peak_id = "MM1_HSA_K562_ATF4",
                               motif_only_for_compared_peak = TRUE)
K562_CEBPB_ATF4_res <- commonPeakResult(commonPeaks = K562_CEBPB_ATF4,
                                        return_common_peak_sites = TRUE,
                                        save_MethMotif_logo = TRUE,
                                        return_summary = TRUE)
K562_CEBPB_ATF4_res$peak_summary
#>                                 percentage_in_original_inputs(%)
#> MM1_HSA_K562_CEBPB_common_peaks                         32.34774
K562_CEBPB_ATF4_peak <- K562_CEBPB_ATF4_res$common_peak_list$MM1_HSA_K562_CEBPB_common_peaks
K562_CEBPB_ATF4_great <- greatAnnotate(peaks = K562_CEBPB_ATF4_peak, return_annotation = TRUE)
K562_CEBPB_ATF4_great_bp <- K562_CEBPB_ATF4_great[which(K562_CEBPB_ATF4_great$category=="BP"),]



# identify CEBPB motif in K562 ATF4 peaks(Supplementary Figure 5A)
# 1) plot ATF4 overal motif logo
ATF4_motif <- searchMotif(id="MM1_HSA_K562_ATF4")
plotLogo(MM_object = ATF4_motif)
# 2) get ATF4 peaks co-bound by CEBPB
ATF4_with_CEBPB <- commonPeaks(target_peak_id = "MM1_HSA_K562_ATF4",
                              motif_only_for_target_peak = TRUE,
                              compared_peak_id = "MM1_HSA_K562_CEBPB",
                              motif_only_for_compared_peak = TRUE)
ATF4_with_CEBPB_res <- commonPeakResult(commonPeaks = ATF4_with_CEBPB,
                                       return_common_peak_sites = TRUE,
                                       save_MethMotif_logo = TRUE)
ATF4_peaks_with_CEBPB <- ATF4_with_CEBPB_res$common_peak_list$MM1_HSA_K562_ATF4_common_peaks

K562_TFBS <- dataBrowser(cell_tissue_name = "K562")
# 3) get all K562 PWM IDs except CEBPB and ATF4
K562_TFBS_no_CEBPB_ATF4 <- K562_TFBS$ID[(K562_TFBS$ID != "MM1_HSA_K562_ATF4" &
                                          K562_TFBS$ID != "MM1_HSA_K562_CEBPB")]
# 4) filter out the peaks also co-bound by other TFs in the ATF4-CEBPB co-binding peaks obtained at step 2.
ATF4_peaks_with_CEBPB_no_other <- exclusivePeaks(user_target_peak_list = list(ATF4_peaks_with_CEBPB),
                                               user_target_peak_id = "MM1_HSA_K562_ATF4",
                                               excluded_peak_id = K562_TFBS_no_CEBPB_ATF4,
                                               motif_only_for_excluded_peak = TRUE)

ATF4_peaks_with_CEBPB_no_other_res <- exclusivePeakResult(exclusivePeaks = ATF4_peaks_with_CEBPB_no_other,
                                                        save_MethMotif_logo = TRUE)


