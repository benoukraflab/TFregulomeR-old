library(TFregulomeR)
library(gplots)

# ATF3 (Meth)Motif in five cell types (Figure 4A)

ATF3_TFBS <- dataBrowser(tf = "ATF3")
for (i in ATF3_TFBS$ID){
    motif_i <- searchMotif(id = i)
    plotLogo(MM_object = motif_i)
}

######################## cofactors in ATF3 binding loci and genomic annotations (Figure 4B) ########################

# HCT116
hct116_tfbs <- dataBrowser(cell_tissue_name = "HCT116")
hct116_tfbs_no_ATF3 <- hct116_tfbs$ID[!(hct116_tfbs$ID %in% "MM1_HSA_HCT116_ATF3")]
atf3_intersect <- intersectPeakMatrix(peak_id_x = "MM1_HSA_HCT116_ATF3",
                                     motif_only_for_id_x =FALSE,
                                     peak_id_y = hct116_tfbs_no_ATF3,
                                     motif_only_for_id_y = TRUE)
cofactorReport(atf3_intersect, cobinding_threshold = 0.1)

# HCT116 ATF3 peak genomic annotation
HCT116_ATF3_peaks <- loadPeaks(id = "MM1_HSA_HCT116_ATF3", includeMotifOnly = TRUE)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
HCT116_ATF3_peaks_genome <- genomeAnnotate(peaks = HCT116_ATF3_peaks, return_annotation = TRUE)
out_3utr <- nrow(HCT116_ATF3_peaks_genome[which(HCT116_ATF3_peaks_genome$annotation=="3UTR"),])
out_5utr <- nrow(HCT116_ATF3_peaks_genome[which(HCT116_ATF3_peaks_genome$annotation=="5UTR"),])
out_exon <- nrow(HCT116_ATF3_peaks_genome[which(HCT116_ATF3_peaks_genome$annotation=="exon"),])
out_intergenic <- nrow(HCT116_ATF3_peaks_genome[which(HCT116_ATF3_peaks_genome$annotation=="intergenic"),])
out_intron <- nrow(HCT116_ATF3_peaks_genome[which(HCT116_ATF3_peaks_genome$annotation=="intron"),])
out_promoter <- nrow(HCT116_ATF3_peaks_genome[which(HCT116_ATF3_peaks_genome$annotation=="promoter-TSS"),])
out_TTS <- nrow(HCT116_ATF3_peaks_genome[which(HCT116_ATF3_peaks_genome$annotation=="TTS"),])

pie <- c(out_3utr, out_5utr, out_exon, out_intergenic, out_intron, out_promoter, out_TTS)

pdf("HCT116_ATF3_genomic_locations.pdf")
pie(pie, labels = c("3' UTR exon","5' UTR exon","exon","intergenic","intron","promoter","TTS"),
    col = c("#009662","#ff44ec","#648c00","#3b7cff","#c4a778","#00313f","#f2006d"))
dev.off()

### get the distribution plot for TSS distance
HCT116_distances_to_TSS <- apply(as.data.frame(HCT116_ATF3_peaks_genome$distanceToTSS),1,
                         function(x) min(as.numeric(unlist(strsplit(x, ";")))))
pdf("HCT116_ATF3_peak_distance_to_TSS.pdf")
plot(density(log10(HCT116_distances_to_TSS+1)), xlab="log10(distance to TSS)", main="HCT116 ATF3 peaks")
abline(v=log10(1000), col="red")
dev.off()

# K562

k562_tfbs <- dataBrowser(cell_tissue_name = "K562")
k562_tfbs_no_atf3 <- k562_tfbs$ID[!(k562_tfbs$ID %in% "MM1_HSA_K562_ATF3") ]
k562_atf3_intersect <- intersectPeakMatrix(peak_id_x = "MM1_HSA_K562_ATF3",
                                          motif_only_for_id_x =FALSE,
                                          peak_id_y = k562_tfbs_no_atf3,
                                          motif_only_for_id_y = TRUE)
cofactorReport(intersectPeakMatrix = k562_atf3_intersect, cobinding_threshold = 0.1)

K562_ATF3_peaks <- loadPeaks(id = "MM1_HSA_K562_ATF3", includeMotifOnly = TRUE)
K562_ATF3_peaks_genome <- genomeAnnotate(peaks = K562_ATF3_peaks, return_annotation = TRUE)
out_3utr <- nrow(K562_ATF3_peaks_genome[which(K562_ATF3_peaks_genome$annotation=="3UTR"),])
out_5utr <- nrow(K562_ATF3_peaks_genome[which(K562_ATF3_peaks_genome$annotation=="5UTR"),])
out_exon <- nrow(K562_ATF3_peaks_genome[which(K562_ATF3_peaks_genome$annotation=="exon"),])
out_intergenic <- nrow(K562_ATF3_peaks_genome[which(K562_ATF3_peaks_genome$annotation=="intergenic"),])
out_intron <- nrow(K562_ATF3_peaks_genome[which(K562_ATF3_peaks_genome$annotation=="intron"),])
out_promoter <- nrow(K562_ATF3_peaks_genome[which(K562_ATF3_peaks_genome$annotation=="promoter-TSS"),])
out_TTS <- nrow(K562_ATF3_peaks_genome[which(K562_ATF3_peaks_genome$annotation=="TTS"),])

pie <- c(out_3utr, out_5utr, out_exon, out_intergenic, out_intron, out_promoter, out_TTS)

pdf("K562_ATF3_genomic_locations.pdf")
pie(pie, labels = c("3' UTR exon","5' UTR exon","exon","intergenic","intron","promoter","TTS"), col = c("#009662","#ff44ec","#648c00","#3b7cff","#c4a778","#00313f","#f2006d"))
dev.off()

### get the distribution plot for TSS distance
K562_distances_to_TSS <- apply(as.data.frame(K562_ATF3_peaks_genome$distanceToTSS),1,
                                function(x) min(as.numeric(unlist(strsplit(x, ";")))))
pdf("K562_ATF3_peak_distance_to_TSS.pdf")
plot(density(log10(K562_distances_to_TSS+1)), xlab="log10(distance to TSS)", main="K562 ATF3 peaks")
abline(v=log10(1000), col="red")
dev.off()

# GM12878
gm_tfbs <- dataBrowser(cell_tissue_name = "GM12878")
gm_tfbs_no_atf3 <- gm_tfbs$ID[!(gm_tfbs$ID %in% "MM1_HSA_GM12878_ATF3")]
gm_atf3_intersect <- intersectPeakMatrix(peak_id_x = "MM1_HSA_GM12878_ATF3",
                                        motif_only_for_id_x = FALSE,
                                        peak_id_y = gm_tfbs_no_atf3,
                                        motif_only_for_id_y = TRUE)
cofactorReport(intersectPeakMatrix = gm_atf3_intersect, cobinding_threshold = 0.1)

GM_ATF3_peaks <- loadPeaks(id = "MM1_HSA_GM12878_ATF3", includeMotifOnly = TRUE)
GM_ATF3_peaks_genome <- genomeAnnotate(peaks = GM_ATF3_peaks, return_annotation = TRUE)

out_3utr <- nrow(GM_ATF3_peaks_genome[which(GM_ATF3_peaks_genome$annotation=="3UTR"),])
out_5utr <- nrow(GM_ATF3_peaks_genome[which(GM_ATF3_peaks_genome$annotation=="5UTR"),])
out_exon <- nrow(GM_ATF3_peaks_genome[which(GM_ATF3_peaks_genome$annotation=="exon"),])
out_intergenic <- nrow(GM_ATF3_peaks_genome[which(GM_ATF3_peaks_genome$annotation=="intergenic"),])
out_intron <- nrow(GM_ATF3_peaks_genome[which(GM_ATF3_peaks_genome$annotation=="intron"),])
out_promoter <- nrow(GM_ATF3_peaks_genome[which(GM_ATF3_peaks_genome$annotation=="promoter-TSS"),])
out_TTS <- nrow(GM_ATF3_peaks_genome[which(GM_ATF3_peaks_genome$annotation=="TTS"),])

pie <- c(out_3utr, out_5utr, out_exon, out_intergenic, out_intron, out_promoter, out_TTS)

pdf("GM12878_ATF3_genomic_locations.pdf")
pie(pie, labels = c("3' UTR exon","5' UTR exon","exon","intergenic","intron","promoter","TTS"), col = c("#009662","#ff44ec","#648c00","#3b7cff","#c4a778","#00313f","#f2006d"))
dev.off()

### get the distribution plot for TSS distance
GM_distances_to_TSS <- apply(as.data.frame(GM_ATF3_peaks_genome$distanceToTSS),1,
                               function(x) min(as.numeric(unlist(strsplit(x, ";")))))
pdf("GM12878_ATF3_peak_distance_to_TSS.pdf")
plot(density(log10(GM_distances_to_TSS+1)), xlab="log10(distance to TSS)", main="GM12878 ATF3 peaks")
abline(v=log10(1000), col="red")
dev.off()

# H1-hESC
h1_tfbs <- dataBrowser(cell_tissue_name = "H1-hESC")
h1_tfbs_no_atf3 <- h1_tfbs$ID[!(h1_tfbs$ID %in% "MM1_HSA_H1-hESC_ATF3")]
h1_tfbs_intersect <- intersectPeakMatrix(peak_id_x = "MM1_HSA_H1-hESC_ATF3",
                                        motif_only_for_id_x = FALSE,
                                        peak_id_y = h1_tfbs_no_atf3,
                                        motif_only_for_id_y = TRUE)
cofactorReport(intersectPeakMatrix = h1_tfbs_intersect, cobinding_threshold = 0.1)

H1_ATF3_peaks <- loadPeaks(id = "MM1_HSA_H1-hESC_ATF3", includeMotifOnly = TRUE)
H1_ATF3_peaks_genome <- genomeAnnotate(peaks = H1_ATF3_peaks, return_annotation = TRUE)

out_3utr <- nrow(H1_ATF3_peaks_genome[which(H1_ATF3_peaks_genome$annotation=="3UTR"),])
out_5utr <- nrow(H1_ATF3_peaks_genome[which(H1_ATF3_peaks_genome$annotation=="5UTR"),])
out_exon <- nrow(H1_ATF3_peaks_genome[which(H1_ATF3_peaks_genome$annotation=="exon"),])
out_intergenic <- nrow(H1_ATF3_peaks_genome[which(H1_ATF3_peaks_genome$annotation=="intergenic"),])
out_intron <- nrow(H1_ATF3_peaks_genome[which(H1_ATF3_peaks_genome$annotation=="intron"),])
out_promoter <- nrow(H1_ATF3_peaks_genome[which(H1_ATF3_peaks_genome$annotation=="promoter-TSS"),])
out_TTS <- nrow(H1_ATF3_peaks_genome[which(H1_ATF3_peaks_genome$annotation=="TTS"),])

pie <- c(out_3utr, out_5utr, out_exon, out_intergenic, out_intron, out_promoter, out_TTS)

pdf("H1-hESC_ATF3_genomic_locations.pdf")
pie(pie, labels = c("3' UTR exon","5' UTR exon","exon","intergenic","intron","promoter","TTS"), col = c("#009662","#ff44ec","#648c00","#3b7cff","#c4a778","#00313f","#f2006d"))
dev.off()

### get the distribution plot for TSS distance
H1_distances_to_TSS <- apply(as.data.frame(H1_ATF3_peaks_genome$distanceToTSS),1,
                             function(x) min(as.numeric(unlist(strsplit(x, ";")))))
pdf("H1-hESC_ATF3_peak_distance_to_TSS.pdf")
plot(density(log10(H1_distances_to_TSS+1)), xlab="log10(distance to TSS)", main="H1-hESC ATF3 peaks")
abline(v=log10(1000), col="red")
dev.off()


# HepG2
hegp2_tfbs <- dataBrowser(cell_tissue_name = "HepG2")
hegp2_tfbs_no_atf3 <- hegp2_tfbs$ID[!(hegp2_tfbs$ID %in% "MM1_HSA_HepG2_ATF3")]
hegp2_atf3_intersect <- intersectPeakMatrix(peak_id_x = "MM1_HSA_HepG2_ATF3",
                                           motif_only_for_id_x = FALSE,
                                           peak_id_y = hegp2_tfbs_no_atf3,
                                           motif_only_for_id_y = TRUE)
cofactorReport(intersectPeakMatrix = hegp2_atf3_intersect, cobinding_threshold = 0.1)

Hepg2_ATF3_peaks <- loadPeaks(id = "MM1_HSA_HepG2_ATF3", includeMotifOnly = TRUE)
Hepg2_ATF3_peaks_genome <- genomeAnnotate(peaks = Hepg2_ATF3_peaks, return_annotation = TRUE)

out_3utr <- nrow(Hepg2_ATF3_peaks_genome[which(Hepg2_ATF3_peaks_genome$annotation=="3UTR"),])
out_5utr <- nrow(Hepg2_ATF3_peaks_genome[which(Hepg2_ATF3_peaks_genome$annotation=="5UTR"),])
out_exon <- nrow(Hepg2_ATF3_peaks_genome[which(Hepg2_ATF3_peaks_genome$annotation=="exon"),])
out_intergenic <- nrow(Hepg2_ATF3_peaks_genome[which(Hepg2_ATF3_peaks_genome$annotation=="intergenic"),])
out_intron <- nrow(Hepg2_ATF3_peaks_genome[which(Hepg2_ATF3_peaks_genome$annotation=="intron"),])
out_promoter <- nrow(Hepg2_ATF3_peaks_genome[which(Hepg2_ATF3_peaks_genome$annotation=="promoter-TSS"),])
out_TTS <- nrow(Hepg2_ATF3_peaks_genome[which(Hepg2_ATF3_peaks_genome$annotation=="TTS"),])

pie <- c(out_3utr, out_5utr, out_exon, out_intergenic, out_intron, out_promoter, out_TTS)

pdf("HepG2_ATF3_genomic_locations.pdf")
pie(pie, labels = c("3' UTR exon","5' UTR exon","exon","intergenic","intron","promoter","TTS"),
    col = c("#009662","#ff44ec","#648c00","#3b7cff","#c4a778","#00313f","#f2006d"))
dev.off()

### get the distribution plot for TSS distance
Hepg2_distances_to_TSS <- apply(as.data.frame(Hepg2_ATF3_peaks_genome$distanceToTSS),1,
                             function(x) min(as.numeric(unlist(strsplit(x, ";")))))
pdf("HepG2_ATF3_peak_distance_to_TSS.pdf")
plot(density(log10(Hepg2_distances_to_TSS+1)), xlab="log10(distance to TSS)", main="HepG2 ATF3 peaks")
abline(v=log10(1000), col="red")
dev.off()

######################## cofactors in ATF3 binding loci and genomic annotations (Figure 4B) ########################

######################## Peak overlapping and GREAT analysis (Figure 4C and D) ########################

#### HCT116 and K562 ####

# common peaks between HCT116 and K562
HCT116_K562_common <- commonPeaks(target_peak_id = "MM1_HSA_HCT116_ATF3",
                                 motif_only_for_target_peak = TRUE,
                                 compared_peak_id = "MM1_HSA_K562_ATF3",
                                 motif_only_for_compared_peak = TRUE)
HCT116_K562_common_res <- commonPeakResult(commonPeaks = HCT116_K562_common,
                                          return_common_peak_sites = TRUE)
HCT116_K562_common_peak <- HCT116_K562_common_res$common_peak_list$MM1_HSA_HCT116_ATF3_common_peaks
# great analysis of common peaks between HCT116 and K562
# load GREAT API package
library(rGREAT)
# load liftover to convert hg38 to hg19
library(liftOver)
HCT116_K562_common_GREAT <- greatAnnotate(peaks = HCT116_K562_common_peak,
                                         return_annotation = TRUE)
HCT116_K562_common_GREAT_bp <- HCT116_K562_common_GREAT[which(HCT116_K562_common_GREAT$category=="BP"),]

# HCT116 exclusive ATF3 peaks
HCT116_exclu <- exclusivePeaks(target_peak_id = "MM1_HSA_HCT116_ATF3",
                              motif_only_for_target_peak = TRUE,
                              excluded_peak_id = "MM1_HSA_K562_ATF3",
                              motif_only_for_excluded_peak = TRUE)
HCT116_exclu_res <- exclusivePeakResult(exclusivePeaks = HCT116_exclu,
                                       return_exclusive_peak_sites = TRUE)
HCT116_exclu_peak <- HCT116_exclu_res$exclusive_peak_list$MM1_HSA_HCT116_ATF3_exclusive_peaks
HCT116_exclu_GREAT <- greatAnnotate(peaks = HCT116_exclu_peak,
                                   return_annotation = TRUE)
HCT116_exclu_GREAT_bp <- HCT116_exclu_GREAT[which(HCT116_exclu_GREAT$category=="BP"),]

# K562 exclusive ATF3 peaks
K562_exclu <- exclusivePeaks(target_peak_id = "MM1_HSA_K562_ATF3",
                            motif_only_for_target_peak = TRUE,
                            excluded_peak_id = "MM1_HSA_HCT116_ATF3",
                            motif_only_for_excluded_peak = TRUE)
K562_exclu_res <- exclusivePeakResult(exclusivePeaks = K562_exclu,
                                     return_exclusive_peak_sites = TRUE)
K562_exclu_peak <- K562_exclu_res$exclusive_peak_list$MM1_HSA_K562_ATF3_exclusive_peaks
K562_exclu_GREAT <- greatAnnotate(peaks = K562_exclu_peak,
                                 return_annotation = TRUE)
K562_exclu_GREAT_bp <- K562_exclu_GREAT[which(K562_exclu_GREAT$category=="BP"),]

#### GM12878, H1-hESC and HepG2 ####
# common peaks amongst GM12878, H1-hESC and HepG2
gm_h1_hepg2_common <- commonPeaks(target_peak_id = c("MM1_HSA_GM12878_ATF3",
                                                    "MM1_HSA_H1-hESC_ATF3",
                                                    "MM1_HSA_HepG2_ATF3"),
                                 motif_only_for_target_peak = TRUE,
                                 compared_peak_id = c("MM1_HSA_GM12878_ATF3",
                                                      "MM1_HSA_H1-hESC_ATF3",
                                                      "MM1_HSA_HepG2_ATF3"),
                                 motif_only_for_compared_peak = TRUE)
gm_h1_hepg2_common_res <- commonPeakResult(commonPeaks = gm_h1_hepg2_common,
                                          return_common_peak_sites = TRUE)
gm_common_peaks <- gm_h1_hepg2_common_res$common_peak_list$MM1_HSA_GM12878_ATF3_common_peaks
h1_common_peaks <- gm_h1_hepg2_common_res$common_peak_list$`MM1_HSA_H1-hESC_ATF3_common_peaks`
hepg2_common_peaks <- gm_h1_hepg2_common_res$common_peak_list$MM1_HSA_HepG2_ATF3_common_peaks
gm_h1_hepg2_common_GREAT <- greatAnnotate(peaks = gm_common_peaks,
                                         return_annotation = TRUE)
gm_h1_hepg2_common_GREAT_bp <- gm_h1_hepg2_common_GREAT[which(gm_h1_hepg2_common_GREAT$category=="BP"),]

# GM12878 exclusive ATF3 targets
GM_exclu <- exclusivePeaks(target_peak_id = c("MM1_HSA_GM12878_ATF3"),
                          motif_only_for_target_peak = TRUE,
                          excluded_peak_id = c("MM1_HSA_H1-hESC_ATF3",
                                               "MM1_HSA_HepG2_ATF3"),
                          motif_only_for_excluded_peak = TRUE)
GM_exclu_res <- exclusivePeakResult(exclusivePeaks = GM_exclu,
                                   return_exclusive_peak_sites = TRUE)
GM_exclu_peak <- GM_exclu_res$exclusive_peak_list$MM1_HSA_GM12878_ATF3_exclusive_peaks
GM_exclu_peak_GREAT <- greatAnnotate(peaks = GM_exclu_peak,
                                    return_annotation = TRUE)
GM_exclu_peak_GREAT_bp <- GM_exclu_peak_GREAT[which(GM_exclu_peak_GREAT$category=="BP"),]

# H1-hESC exclusive ATF3 targets
h1_exclu <- exclusivePeaks(target_peak_id = c("MM1_HSA_H1-hESC_ATF3"),
                          motif_only_for_target_peak = TRUE,
                          excluded_peak_id = c("MM1_HSA_HepG2_ATF3",
                                               "MM1_HSA_GM12878_ATF3"),
                          motif_only_for_excluded_peak = TRUE)
h1_exclu_res <- exclusivePeakResult(exclusivePeaks = h1_exclu,
                                   return_exclusive_peak_sites = TRUE)
h1_exclu_peak <- h1_exclu_res$exclusive_peak_list$`MM1_HSA_H1-hESC_ATF3_exclusive_peaks`
h1_exclu_peak_GREAT <- greatAnnotate(peaks = h1_exclu_peak, return_annotation = TRUE)
h1_exclu_peak_GREAT_bp <- h1_exclu_peak_GREAT[which(h1_exclu_peak_GREAT$category=="BP"),]

# HepG2 exclusive ATF3 targets
hepg2_exclu <- exclusivePeaks(target_peak_id = c("MM1_HSA_HepG2_ATF3"),
                             motif_only_for_target_peak = TRUE,
                             excluded_peak_id = c("MM1_HSA_GM12878_ATF3",
                                                  "MM1_HSA_H1-hESC_ATF3"),
                             motif_only_for_excluded_peak = TRUE)
hepg2_exclu_res <- exclusivePeakResult(exclusivePeaks = hepg2_exclu,
                                      return_exclusive_peak_sites = TRUE)
hepg2_exclu_peak <- hepg2_exclu_res$exclusive_peak_list$MM1_HSA_HepG2_ATF3_exclusive_peaks
hepg2_exclu_GREAT <- greatAnnotate(peaks = hepg2_exclu_peak, return_annotation = TRUE)
hepg2_exclu_GREAT_bp <- hepg2_exclu_GREAT[which(hepg2_exclu_GREAT$category=="BP"),]

######################## Peak overlapping and GREAT analysis (Figure 4C and D) ########################
