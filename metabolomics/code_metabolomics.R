#Set working directory
setwd()

# Load in packages
library(tidyverse)
library(mixOmics)
library(ggpubr)
library(vegan)
library(caret)
library(limma) 
library(patchwork)
library(rstatix)
library(Spectra) 
library(MsBackendMgf) 
library(homologueDiscoverer)
library(ggrepel)

###########################################################################################
# Function for PEG removal
plotAnnotatedStatic <- function(annotated, legend_setting = "bottom"){
  annotated <- mutate(annotated,
                      homologue_id = if_else(is.na(homologue_id), 0L, homologue_id),
                      homologue_id = as.factor(homologue_id))
  colourCount = length(unique(annotated$homologue_id))
  getPalette = colorRampPalette(RColorBrewer::brewer.pal(9, "Set1"))
  ncolor = length(unique(annotated$homologue_id))
  annotated <- arrange(annotated, desc(homologue_id))
  g <- ggplot(annotated, aes(group = homologue_id)) +
    geom_point(aes(x = rt, y = mz, color = homologue_id, size = homologue_id,
                   shape = homologue_id, alpha = homologue_id)) +
    geom_line(data = filter(annotated, homologue_id != 0),
              aes(x = rt, y = mz, group = homologue_id), alpha = 0.4, size = 0.1) +
    ggtitle("Annotated Peak Table") +
    scale_colour_manual(values = c("black", getPalette(colourCount-1))) +
    scale_alpha_manual(values =  c(0.2, rep(0.85, times = ncolor))) +
    scale_size_manual(values = c(0.1, rep(1.5, times = ncolor))) +
    scale_shape_manual(values = c(1, rep(19, times = ncolor))) +
    xlab("Retention Time (s)") +
    ylab("Mass to Charge Ratio") +
    theme(legend.position=legend_setting, text = element_text(family="mono"))
  return(g)
}
###########################################################################################

# Read in data
feature_table <- read_csv("gnps_quant.csv") 
colnames(feature_table)[3] <- "RT"
metadata <- read.csv("metadata_VEOIBD.csv")
sample_order <- read.csv("sequence.csv") %>% dplyr::select(-1)
annotations <- read.delim("merged_results_with_gnps_all_libs.tsv")
annotations$X.Scan. <- as.character(annotations$X.Scan.)

info_feature <- feature_table %>% dplyr::select(1:3,7)
colnames(info_feature) <- c("Feature", "mz", "RT", "Corr_ID")
info_feature$Feature <- as.character(info_feature$Feature)

# Filter out features eluting in the dead volume 
info_feature_filter <- info_feature %>% dplyr::filter(RT >= 0.70) 

info_feature_complete <- info_feature_filter %>% 
  left_join(annotations, by = c("Feature" = "X.Scan.")) %>% 
  dplyr::select(1:4,18)

# Data table
data <- feature_table %>%
  dplyr::filter(RT >= 0.70) %>% 
  column_to_rownames("row ID") %>% 
  dplyr::select(contains("Peak")) %>% t() %>% as.data.frame() %>%
  rownames_to_column("SampleID") %>% arrange(SampleID) %>% 
  distinct(SampleID, .keep_all = TRUE) %>%
  dplyr::filter(!(SampleID %in% c("PT_CHANG_00.mzML Peak area", "PT_RODRI_001.mzML Peak area"))) # excluded due to clinical phenotype 

data$SampleID <- gsub(".mzML Peak area", "", data$SampleID)

# Metadata
metadata_metabolomics <- metadata %>% 
  dplyr::left_join(data %>% dplyr::select(SampleID), by = c("sample_ID" = "SampleID")) %>%
  left_join(sample_order, by = c("sample_ID" = "File.Name"))

# Characteristics of the study population
metadata %>% count(Cohort, Sex)
metadata %>%
  group_by(Cohort) %>%
  summarise(mean_age = mean(Age, na.rm = TRUE))
metadata %>%
  group_by(Cohort) %>%
  summarise(
    min_age = min(Age, na.rm = TRUE),
    max_age = max(Age, na.rm = TRUE))
wilcox.test(Age ~ Cohort, data = metadata)
chisq.test(table(metadata$Cohort, metadata$Sex))
time_from_diagnosis <- metadata %>%
  filter(Cohort == "VEOIBD") %>%
  mutate(Time_From_Diagnosis = Age - Age_IBD_diagnosis) %>%
  summarise(
    mean_time = mean(Time_From_Diagnosis, na.rm = TRUE),
    sd_time = sd(Time_From_Diagnosis, na.rm = TRUE))

# Discard first blanks and pools used for conditioning
data_final <- data %>% 
  dplyr::filter(!(SampleID %in% c("Blank_1", "Blank_2", "sixmix_1", "pool_qc_veoIBD_1", "pool_qc_veoIBD_2", "pool_qc_veoIBD_3")))

# Investigate total peak area
data_TIC <- data.frame(TIC = rowSums(data_final %>% column_to_rownames("SampleID"))) %>% 
  rownames_to_column("SampleID") %>% left_join(sample_order, by = c("SampleID" = "File.Name"))

data_TIC %>% ggscatter("Run_Order", "TIC", add = "reg.line", label = "SampleID") +
  stat_cor() + ylim(0, 2.5e10) # pool and some samples have high TIC

# Samples only
data_TIC %>% dplyr::filter(!(str_detect(pattern = "Blank|six|pool", SampleID))) %>%
  ggscatter("Run_Order", "TIC", add = "reg.line") + ylim(0, 2.5e10) +
  stat_cor() # possibly 4 outliers with very high TIC

data_TIC %>% dplyr::filter(!(str_detect(pattern = "Blank|six|pool", SampleID))) %>%
  dplyr::filter(TIC < 4e9) %>%
  ggscatter("Run_Order", "TIC", add = "reg.line") + ylim(0, 4e9) +
  stat_cor() # possibly 4 outliers with very high TIC

# Check Blank, QCpool and QCmix
data_TIC %>% dplyr::filter(str_detect(pattern = "pool_qc", SampleID)) %>% 
  ggscatter("Run_Order", "TIC", add = "reg.line") +
  stat_cor() + ylim(0, 2e10)

data_TIC %>% dplyr::filter(str_detect(pattern = "sixmix", SampleID)) %>% 
  ggscatter("Run_Order", "TIC", add = "reg.line") +
  stat_cor() + ylim(0, 6e9)

data_TIC %>% dplyr::filter(str_detect(pattern = "Blank", SampleID)) %>% 
  ggscatter("Run_Order", "TIC", add = "reg.line", label = "SampleID") +
  stat_cor() + ylim(0, 1e9)


# Check sample type
sample_tic <- data_TIC %>% dplyr::filter(!(str_detect(pattern = "Blank|pool|mix", SampleID))) %>% summarise(mean(TIC))
pool_tic <- data_TIC %>% dplyr::filter(str_detect(pattern = "pool", SampleID)) %>% summarise(mean(TIC))
blank_tic <- data_TIC %>% dplyr::filter(str_detect(pattern = "Blank", SampleID)) %>% summarise(mean(TIC))

# Ratios
ratio_tic_pb <- pool_tic/sample_tic #pools were concentrated or samples with higher TIC caused high concentration
ratio_tic_sb <- sample_tic/blank_tic


# Check internal standard and remove sample where there was a shift
fbmn_IS <- annotations %>% dplyr::filter(str_detect(Compound_Name, regex("sulf", ignore_case = TRUE))) %>% 
  distinct(X.Scan., .keep_all = TRUE) %>% dplyr::filter(Organism != "BILELIB19")

# Extract IS feature from the table
table_IS <- data_final %>% column_to_rownames("SampleID") %>% t() %>% as.data.frame() %>% rownames_to_column("ID") %>% 
  dplyr::filter(ID %in% fbmn_IS$X.Scan.) %>% column_to_rownames("ID") %>% t() %>% as.data.frame() %>% 
  rownames_to_column("filename") %>% dplyr::filter(!(str_detect(filename, pattern = "Blank|pool|mix"))) %>%
  left_join(metadata_metabolomics, by = c("filename" = "sample_ID")) %>%
  dplyr::select(filename, `15938`, Run_Order)

colnames(table_IS)[2] <- "sulfadimethoxine"

table_IS %>% ggscatter(x = "Run_Order", y = "sulfadimethoxine", add = c("reg.line")) +
  stat_cor() + ylim(0, 7e6)

table_IS %>% ggbarplot(x = "Run_Order", y = "sulfadimethoxine", xlab = "Run Order", 
                       ylab = "Peak Area sulfadimethoxine", title = "Internal Standard Acquisition") +
  geom_hline(yintercept = mean(table_IS$sulfadimethoxine, na.rm = TRUE), linetype = "dashed", color = "blue")

cv_is <- sd(table_IS$sulfadimethoxine)/mean(table_IS$sulfadimethoxine)
# CV is 20% is okish


# Check features per sample type
data_blank <- data_final %>% dplyr::filter(str_detect(pattern = "Blank", SampleID))
data_pool <- data_final %>% dplyr::filter(str_detect(pattern = "pool", SampleID))
data_sixmix <- data_final %>% dplyr::filter(str_detect(pattern = "sixmix", SampleID))

# Blanks
blanks_feature_info <- data.frame(Feature = colnames(data_blank)[-1],
                                  Mean_blank = data_blank %>% column_to_rownames("SampleID") %>% colMeans(), 
                                  SD_blank =  data_blank %>% column_to_rownames("SampleID") %>% apply(2, sd)) %>%
  dplyr::mutate(CV_blank = SD_blank/Mean_blank) %>% left_join(annotations, by = c("Feature" = "X.Scan.")) %>% 
  left_join(info_feature) %>% dplyr::select(Feature, mz, RT, Compound_Name, SpectrumID, 
                                            Precursor_MZ, Mean_blank, SD_blank, CV_blank) %>% 
  dplyr::filter(Mean_blank > 0) %>% arrange(desc(Mean_blank))

# QCpool
pools_feature_info <- data.frame(Feature = colnames(data_pool)[-1],
                                 Mean_pool = data_pool %>% column_to_rownames("SampleID") %>% colMeans(), 
                                 SD_pool =  data_pool %>% column_to_rownames("SampleID") %>% apply(2, sd)) %>%
  dplyr::mutate(CV_pool = SD_pool/Mean_pool) %>% left_join(annotations, by = c("Feature" = "X.Scan.")) %>% 
  left_join(info_feature) %>% dplyr::select(Feature, mz, RT, Compound_Name, SpectrumID, 
                                            Precursor_MZ, Mean_pool, SD_pool, CV_pool) %>% 
  dplyr::filter(Mean_pool > 0) %>% arrange(desc(Mean_pool))

# QCmix
sixmix_feature_info <- data.frame(Feature = colnames(data_sixmix)[-1],
                                  Mean_sixmix = data_sixmix %>% column_to_rownames("SampleID") %>% colMeans(), 
                                  SD_sixmix = data_sixmix %>% column_to_rownames("SampleID") %>% apply(2, sd)) %>%
  dplyr::mutate(CV_sixmix = SD_sixmix/Mean_sixmix) %>% left_join(annotations, by = c("Feature" = "X.Scan.")) %>% 
  left_join(info_feature) %>% dplyr::select(Feature, mz, RT, Compound_Name, SpectrumID, 
                                            Precursor_MZ, Mean_sixmix, SD_sixmix, CV_sixmix) %>% 
  dplyr::filter(Mean_sixmix > 0) %>% arrange(desc(Mean_sixmix))
# CVs of the six standards are ranging between 8% to 17%


# Features to be removed Pools/Blank < 5
feature_to_remove <- blanks_feature_info %>% left_join(pools_feature_info) %>%
  dplyr::filter(Mean_blank > 0) %>% 
  dplyr::mutate(Pool_Blank = Mean_pool/Mean_blank) %>% 
  dplyr::filter(Pool_Blank < 5 | is.na(Pool_Blank))

# Data with blank removal
data_clean <- data_final %>% dplyr::select(-c(feature_to_remove$Feature)) %>%
  column_to_rownames("SampleID") %>%
  select_if(~sum(.) != 0) %>% rownames_to_column("SampleID")

# PCA raw data
PCA_raw <- mixOmics::pca(data_clean %>% select_at(vars(-one_of(nearZeroVar(., names = TRUE)))) %>%
                           column_to_rownames("SampleID"), ncomp = 2, center = TRUE, scale = TRUE)
PCA_raw_scores <- data.frame(PCA_raw$variates$X) %>% rownames_to_column("sample_ID") %>% left_join(metadata_metabolomics)

i <- "Cohort"

PCA_raw_plot <- PCA_raw_scores %>%
  ggscatter(x = "PC1", y = "PC2", color = i, alpha = 0.6, label = "sample_ID",
            title = paste("PCA - Fecal Metabolome", i, sep = " "),
            xlab = paste("PC1 (", round(PCA_raw$prop_expl_var$X[1]*100, digits = 1),"%)", sep = ""), 
            ylab = paste("PC2 (", round(PCA_raw$prop_expl_var$X[2]*100, digits = 1),"%)", sep = ""),
            ggtheme = theme_classic()) +
  geom_point(data = PCA_raw_scores %>% group_by((!!sym(i))) %>% 
               summarise_at(vars(matches("PC")), mean), aes(PC1, PC2, color = (!!sym(i))), size = 4, shape = 8) +
  theme(plot.title = element_text(size = 10), axis.title = element_text(size = 8),
        axis.text = element_text(size = 4)) + coord_fixed()

# there are 4 samples to potentially be removed due to outliers (PT_S83, PT_S97, PT_AC_095, PT_AC_175)


# Remove plastic contamination
feature_check <- info_feature %>% dplyr::select(1:3) %>% 
  dplyr::filter(Feature %in% colnames(data_clean)) %>% as_tibble()
colnames(feature_check) <- c("peak_id", "mz", "rt")
feature_check$mz <- as.double(feature_check$mz)
feature_check$rt <- as.double(feature_check$rt)

feature_pool <- data_clean %>% 
  dplyr::filter(SampleID == "pool_qc_veoIBD_5") %>% 
  column_to_rownames("SampleID") %>%
  t() %>% as.data.frame() %>% rownames_to_column("peak_id")

feature_check_pool <- feature_check %>% left_join(feature_pool)

feature_check_pool$peak_id <- as.integer(feature_check_pool$peak_id)
colnames(feature_check_pool)[4] <- "intensity"

# Remove PEGs
peak_table_PEG <- detectHomologues(feature_check_pool, mz_steps = c(44.0262, 88.0524, 132.0786),
                                   rt_min = 0.01, rt_max = 100, ppm_tolerance = 5, 
                                   min_series_length = 4, search_mode = "targeted", 
                                   step_mode = "increment", verbose = TRUE)

#plotAnnotatedStatic(peak_table_PEG)

sdb <- sdbCreate(peak_table_PEG, sample_origin = "data")
sdb1_info <- info_feature_complete %>% dplyr::filter(Feature %in% sdb$sample_peak_id)

# Remove untargeted 
peak_table_untargeted <- detectHomologues(feature_check_pool, mz_min = 10, mz_max = 50, 
                                          rt_min = 0.1, rt_max = 100, ppm_tolerance = 5, 
                                          min_series_length = 5, search_mode = "untargeted", 
                                          step_mode = "increment", verbose = TRUE)

#plotAnnotatedStatic(peak_table_untargeted)

sdb2 <- sdbCreate(peak_table_untargeted, sample_origin = "data")
sdb2_info <- info_feature_complete %>% dplyr::filter(Feature %in% sdb2$sample_peak_id)

# Combine tables
overlapping_ids <- sdbCheckContained(sdb, peak_table_untargeted, ppm_tolerance = 0.5, rt_tolerance = 0.5)
print(overlapping_ids)
sdb <- sdbPush(sdb, filter(peak_table_untargeted, !(homologue_id %in% overlapping_ids)))
summary_sdb <- sdbSummarize(sdb)
sdb$sample_peak_id <- as.character(sdb$sample_peak_id)
features_contaminat <- info_feature_complete %>% 
  dplyr::filter(Feature %in% sdb$sample_peak_id) %>%
  dplyr::filter(is.na(Compound_Name) | str_detect(pattern = "glycol", Compound_Name))

# Check one PEG
peg <- data_clean %>% dplyr::select(SampleID, `5185`) %>% 
  left_join(metadata_metabolomics, by = c("SampleID" = "sample_ID")) 
colnames(peg)[2] <- "PEG"

peg %>% dplyr::mutate(LogPEG = log2(PEG+1)) %>% 
  ggscatter(x = "SampleID", y = "PEG") + scale_color_viridis_c()

# Clean data from contaminants
data_clean2 <- data_clean %>% 
  dplyr::select(-c(features_contaminat$Feature))

# Features to be removed Pool/QCmix < 5
feature_to_remove_mix <- sixmix_feature_info %>% left_join(pools_feature_info) %>% 
  dplyr::filter(Mean_sixmix > 0) %>% 
  dplyr::mutate(Pool_Mix = Mean_pool/Mean_sixmix) %>% 
  dplyr::filter(Pool_Mix < 5 | is.na(Pool_Mix)) %>% 
  dplyr::filter(!(Feature %in% feature_to_remove$Feature)) %>% 
  dplyr::filter(!(Feature %in% features_contaminat$Feature))

# Data with blank2 removal
data_clean3 <- data_clean2 %>% dplyr::select(-c(feature_to_remove_mix$Feature))

# Remove feature before 0.70 min and after 9.5 min
feature_to_rt <- info_feature_complete %>% dplyr::filter(RT < 0.70 | RT > 9.5) %>%
  dplyr::filter(!(Feature %in% feature_to_remove$Feature))  %>%
  dplyr::filter(!(Feature %in% feature_to_remove_mix$Feature)) %>% 
  dplyr::filter(!(Feature %in% features_contaminat$Feature)) %>% 
  dplyr::filter(!(Feature %in% c(42452, 41947, 43385, 43414)))

# Final cleaned table
data_clean4 <- data_clean3 %>% dplyr::select(-c(feature_to_rt$Feature)) 

# Clean contaminant found via MN
contaminants <- read_csv("edge_table_deltamz.csv") 

contaminants_2 <- contaminants %>% 
  dplyr::filter((deltamz >= 44.02 & deltamz <= 44.03) | 
                  (deltamz >= -44.03 & deltamz <= -44.03) |
                  (deltamz >= -88.06 & deltamz <= -88.05) |
                  (deltamz >= 88.05 & deltamz <= 88.06))

contaminants_3 <- contaminants_2 %>%
  separate(`shared name`, into = c("ID1", "ID2"), sep = " \\(-\\) ", remove = FALSE, extra = "drop") %>%
  pivot_longer(cols = c("ID1", "ID2"), names_to = "ID_type", values_to = "sample_peak_id") %>%
  dplyr::select(-ID_type, -`shared name`) %>% 
  dplyr::select(deltamz, sample_peak_id) %>%
  distinct(sample_peak_id, .keep_all = TRUE) %>%
  dplyr::filter(sample_peak_id %in% colnames(data_clean4))

data_clean5 <- data_clean4 %>%
  dplyr::select(-c(contaminants_3$sample_peak_id)) 

# PCA raw data
PCA_raw <- mixOmics::pca(data_clean5 %>% select_at(vars(-one_of(nearZeroVar(., names = TRUE)))) %>%
                           column_to_rownames("SampleID"), ncomp = 2, center = TRUE, scale = TRUE)
PCA_raw_scores <- data.frame(PCA_raw$variates$X) %>% rownames_to_column("sample_ID") %>% left_join(metadata_metabolomics)

i <- "Cohort"

PCA_raw_plot <- PCA_raw_scores %>%
  ggscatter(x = "PC1", y = "PC2", color = i, alpha = 0.6,
            title = paste("PCA - Fecal Metabolome", i, sep = " "),
            xlab = paste("PC1 (", round(PCA_raw$prop_expl_var$X[1]*100, digits = 1),"%)", sep = ""), 
            ylab = paste("PC2 (", round(PCA_raw$prop_expl_var$X[2]*100, digits = 1),"%)", sep = ""),
            ggtheme = theme_classic()) +
  geom_point(data = PCA_raw_scores %>% group_by((!!sym(i))) %>% 
               summarise_at(vars(matches("PC")), mean), aes(PC1, PC2, color = (!!sym(i))), size = 4, shape = 8) +
  theme(plot.title = element_text(size = 10), axis.title = element_text(size = 8),
        axis.text = element_text(size = 4)) + coord_fixed()

# Remove 3 contaminated samples (PT_S83, PT_S97, PT_AC_095) and 3 samples with low TIC (HC_98, PT_PEDRZ_012, PT_BHATT_018)
data_filter <- data_clean5 %>% dplyr::filter(!(SampleID %in% c("PT_S83", "PT_S97", "PT_AC_095", 
                                                               "HC_98", "PT_PEDRZ_012", "PT_BHATT_018"))) # these 3 have low tic

# RCLR transformation
data_final_clr <- decostand(data_filter %>% column_to_rownames("SampleID"), method = "rclr")

# PCA
PCA_whole <- mixOmics::pca(data_final_clr %>% select_at(vars(-one_of(nearZeroVar(., names = TRUE)))),
                           ncomp = 2, center = TRUE, scale = TRUE)
PCA_whole_scores <- data.frame(PCA_whole$variates$X) %>% 
  rownames_to_column("sample_ID") %>% left_join(metadata_metabolomics)

i <- "Cohort"

PCA_plot <- PCA_whole_scores %>%
  ggscatter(x = "PC1", y = "PC2", color = i, alpha = 0.6,
            title = paste("PCA - Fecal"),
            xlab = paste("PC1 (", round(PCA_whole$prop_expl_var$X[1]*100, digits = 1),"%)", sep = ""), 
            ylab = paste("PC2 (", round(PCA_whole$prop_expl_var$X[2]*100, digits = 1),"%)", sep = ""),
            ggtheme = theme_classic()) +
  geom_point(data = PCA_whole_scores %>% group_by((!!sym(i))) %>% 
               summarise_at(vars(matches("PC")), mean), aes(PC1, PC2, color = (!!sym(i))), size = 4, shape = 8) +
  theme(plot.title = element_text(size = 10), axis.title = element_text(size = 8),
        axis.text = element_text(size = 4)) + coord_fixed()

# Keep sample only
data_sample <- data_filter %>% 
  dplyr::filter(!(str_detect(pattern = "six|Blank|pool", SampleID)))

# CLR transformation
data_final_clr <- decostand(data_sample %>% column_to_rownames("SampleID"), method = "rclr")


#write_csv(x = data_final_clr %>% rownames_to_column("SampleID"), file = "table_metabolomics_clr.csv")

# PCA
PCA_whole <- mixOmics::pca(data_final_clr %>% select_at(vars(-one_of(nearZeroVar(., names = TRUE)))),
                           ncomp = 2, center = TRUE, scale = TRUE)
PCA_whole_scores <- data.frame(PCA_whole$variates$X) %>% rownames_to_column("sample_ID") %>% left_join(metadata_metabolomics)

i <- "Cohort"

PCA_plot <- PCA_whole_scores %>%
  ggscatter(x = "PC1", y = "PC2", color = i, alpha = 0.6,
            title = paste("PCA - Fecal Metabolome"),
            xlab = paste("PC1 (", round(PCA_whole$prop_expl_var$X[1]*100, digits = 1),"%)", sep = ""), 
            ylab = paste("PC2 (", round(PCA_whole$prop_expl_var$X[2]*100, digits = 1),"%)", sep = ""),
            ggtheme = theme_classic()) +
  geom_point(data = PCA_whole_scores %>% group_by((!!sym(i))) %>% 
               summarise_at(vars(matches("PC")), mean), aes(PC1, PC2, color = (!!sym(i))), size = 4, shape = 8) +
  scale_color_manual(values = c("VEOIBD" = "#ff7f00", "Healthy" = "#6a3d9a")) + 
  theme(plot.title = element_text(size = 10), axis.title = element_text(size = 8),
        axis.text = element_text(size = 6)) + coord_fixed()

# PERMANOVA
dist_metabolites <- vegdist(data_final_clr, method = "euclidean")
disper_donor <- betadisper(dist_metabolites, PCA_whole_scores$Cohort)
anova(disper_donor)
permanova <- adonis2(dist_metabolites ~ Cohort + Age, PCA_whole_scores, na.action = na.omit, by = "terms")

#ggsave(filename = "VEOIBD_PCA_fecal.svg", plot = PCA_plot, device = "svg", width = 4, height = 3, dpi = "retina")


# PLSDA - Cohort
PLSDA_cohort <- mixOmics::plsda(data_final_clr %>% select_at(vars(-one_of(nearZeroVar(., names = TRUE)))),
                                     PCA_whole_scores$Cohort, ncomp = 2, scale = TRUE)
PLSDA_cohort_scores <- data.frame(PLSDA_cohort$variates$X) %>% 
  rownames_to_column("sample_ID") %>% left_join(metadata_metabolomics)

PLSDA_cohort_plot <- PLSDA_cohort_scores %>%
  ggscatter(x = "comp1", y = "comp2", color = "Cohort", alpha = 0.6, title = "PLSDA - Cohort",
            xlab = paste("Component 1 (", round(PLSDA_cohort$prop_expl_var$X[1]*100, digits = 1),"%)", sep = ""), 
            ylab = paste("Component 2 (", round(PLSDA_cohort$prop_expl_var$X[2]*100, digits = 1),"%)", sep = ""),
            legend.title = "Group", ggtheme = theme_classic()) +
  geom_point(data = PLSDA_cohort_scores %>% group_by(Cohort) %>% 
               summarise_at(vars(matches("comp")), mean), aes(comp1, comp2, color = Cohort), size = 3, shape = 8) +
  scale_color_manual(values = c("VEOIBD" = "#ff7f00", "Healthy" = "#6a3d9a")) +
  theme(plot.title = element_text(size = 10),axis.title = element_text(size = 8),
        axis.text = element_text(size = 5)) + coord_fixed()

#ggsave(filename = "VEOIBD_PLSDA_cohort_scoreplot.svg", plot = PLSDA_cohort_plot, device = "svg", width = 5, height = 3, dpi = "retina")

Loadings_cohort <- plotLoadings(PLSDA_cohort, plot = FALSE, contrib = "max") 
Loadings_cohort <- Loadings_cohort$X %>%
  rownames_to_column() %>% dplyr::select(rowname, GroupContrib)

perf_plsda_cohort <- perf(PLSDA_cohort, validation = "Mfold", folds = 4, nrepeat = 99, progressBar = TRUE) 
plot(perf_plsda_cohort, legend = TRUE)
perf_plsda_cohort$error.rate
#ggsave(filename = "VEOIBD_performance_plsda_cohort.svg", plot = perf_plsda_cohort, device = "svg", width = 5, height = 3, dpi = "retina")

#pdf("VEOIBD_PLS_DA_performance.pdf", width = 4.5, height = 3.5)
#plot(perf_plsda_cohort, legend = FALSE)
#dev.off()

# Extract VIPs
VIPs_cohort <- as.data.frame(mixOmics::vip(PLSDA_cohort))
VIPs_cohort_filter <- dplyr::filter(VIPs_cohort, VIPs_cohort$comp1 > 1)
VIPs_cohort_filter$ID <- rownames(VIPs_cohort_filter)
VIPs_cohort_select <- VIPs_cohort_filter %>% dplyr::select(ID, comp1)
VIPs_cohort_Load <- VIPs_cohort_select %>% 
  left_join(Loadings_cohort, by = c("ID" = "rowname")) %>% 
  left_join(info_feature_complete, by = c("ID" = "Feature")) %>% arrange(desc(comp1))

#write_csv(x = VIPs_cohort_Load, file = "VIPs_cohort_Load.csv")

# Selected features based on VIP
vip_healthy <- VIPs_cohort_Load %>% dplyr::filter(GroupContrib == "Healthy") %>% head(500)
vip_ibd <- VIPs_cohort_Load %>% dplyr::filter(GroupContrib == "VEOIBD") %>% head(500)

data_vip <- data_sample %>%
  dplyr::select("SampleID", vip_healthy$ID, vip_ibd$ID) %>%
  dplyr::mutate(Healthy = rowSums(select(., vip_healthy$ID))) %>%
  dplyr::mutate(IBD = rowSums(select(., vip_ibd$ID))) %>%
  dplyr::mutate(Ratio = log(Healthy/IBD)) %>%
  dplyr::select(SampleID, Ratio) %>%
  left_join(metadata_metabolomics, by = c("SampleID" = "sample_ID"))

# Plot ratio
plot_ratio <- data_vip %>%
  ggboxplot(y = "Ratio", x = "Cohort", add = "jitter", size = 0.35, ylab = "Log(Healthy/VEO-IBD)",
            add.params = list(color = "Cohort", size = 1.5), legend = "none",
            fill = "Cohort", 
            palette = c("#6a3d9a", "#ff7f00"),
            alpha = 0.5,
            title = "Differential features from PLS-DA models") +
  stat_compare_means() +
  theme(plot.title = element_text(size = 8), axis.title = element_text(size = 7),
        axis.text = element_text(size = 6))
  
#ggsave(filename = "ratio_features_PLSDA.svg", plot = plot_ratio, device = "svg", width = 1.7, height = 3.5, dpi = "retina")

##### N ACYL LIPIDS #####
nacyl_all <- info_feature_complete %>% 
  dplyr::filter(str_detect(pattern = "C.*:", Compound_Name)) %>%
  dplyr::filter(!(str_detect(pattern = "bile|Collision|IIN|CHEBI", Compound_Name))) %>% 
  dplyr::filter(!(Feature %in% c("40624", "38982",    # remove not N-acyl lipids
                                 "13298", "14784", "7808", "10766", "5104", "13439", "12847", "15138"))) # remove misannotations

acyl_interest <-  VIPs_cohort_Load %>% 
  dplyr::filter(str_detect(pattern = "C.*:", Compound_Name)) %>%
  dplyr::filter(!(str_detect(pattern = "bile|Collision|IIN|CHEBI", Compound_Name)))

acyl_short_chain <- acyl_interest %>% 
  dplyr::filter(str_detect(pattern = "C2|C3|C4|C5|C6", Compound_Name)) %>% 
  dplyr::filter(!ID %in% c("3318", "11337"))

meth_C3 <- nacyl_all %>% dplyr::filter(Feature == "7255") 
meth_C4 <- nacyl_all %>% dplyr::filter(Feature == "10409") 
meth_C5 <- nacyl_all %>% dplyr::filter(Feature == "13798") 
phe_C2 <- nacyl_all %>% dplyr::filter(Feature == "10149") 
phe_C3 <- nacyl_all %>% dplyr::filter(Feature == "12767") 
phe_C4 <- nacyl_all %>% dplyr::filter(Feature == "15419") 
phe_C5 <- nacyl_all %>% dplyr::filter(Feature == "19169") 
leu_C4 <- nacyl_all %>% dplyr::filter(Feature == "13758") 
leu_C5 <- nacyl_all %>% dplyr::filter(Feature == "17502") 
tyr_C3 <- nacyl_all %>% dplyr::filter(Feature == "7167") 
tyr_C4 <- nacyl_all %>% dplyr::filter(Feature == "9593") 

data_acyl_sum <- data_final_clr %>%
  rownames_to_column() %>%
  dplyr::select(rowname, '7255', '10409', '13798', '10149', '12767', '15419', '19169',
                '13758', '17502', '7167', '9593') %>%
  pivot_longer(
    cols = -rowname,
    names_to = "Feature",
    values_to = "Abundance"
  ) %>%
  left_join(metadata, by = c("rowname" = "sample_ID"))

p_values <- compare_means(Abundance ~ Cohort, group.by = "Feature", 
                          data = data_acyl_sum, method = "wilcox.test")

plot_data_acyl <- ggboxplot(data_acyl_sum, x = "Feature", y = "Abundance",
                            add = "jitter", ylab = "rclr Peak area", legend = "none",
                            add.params = list(color = "Cohort"),
                            fill = "Cohort", 
                            palette = c("#6a3d9a", "#ff7f00"),
                            alpha = 0.5,
                            title = "N acyl lipids") +
  theme(plot.title = element_text(size = 8),axis.title = element_text(size = 7),
        axis.text = element_text(size = 6))

#ggsave(filename = "plot_data_acyl.svg", plot = plot_data_acyl, device = "svg", width = 7.5, height = 3, dpi = "retina")

### Make volcano plot for the N ACYL LIPIDS
data_sum_all <- data_sample %>% 
  column_to_rownames("SampleID") %>% 
  dplyr::select(where(~ sum(.) != 0)) %>%
  dplyr::mutate(row_sum = rowSums(.)) %>%
  dplyr::mutate(across(everything(), ~ . / row_sum)) %>%
  dplyr::select(-row_sum) %>% 
  rownames_to_column("SampleID") %>%
  left_join(metadata_metabolomics %>% dplyr::select(sample_ID, Cohort), by = c("SampleID" = "sample_ID")) %>% 
  dplyr::select(-SampleID) %>% 
  group_by(Cohort) %>%
  summarise(across(everything(), mean)) 

data_nacyl_log2fc <- data_sum_all %>%
  pivot_longer(cols = -Cohort, names_to = "Feature", values_to = "Value") %>% 
  dplyr::filter(Feature %in% nacyl_all$Feature) %>%
  pivot_wider(names_from = Cohort, values_from = Value) %>%
  mutate(across(c("VEOIBD", "Healthy"), ~ifelse(. == 0, 1e-8, .))) %>%
  mutate(Fold_Change = VEOIBD / Healthy) %>%
  mutate(Log2FC = log2(Fold_Change)) %>%
  arrange(desc(Log2FC))

p_values_nacyl <- data_sample %>%
  column_to_rownames("SampleID") %>% 
  dplyr::select(where(~ sum(.) != 0)) %>%
  dplyr::mutate(row_sum = rowSums(.)) %>%
  dplyr::mutate(across(everything(), ~ . / row_sum)) %>%
  dplyr::select(-row_sum) %>% 
  rownames_to_column("SampleID") %>%
  left_join(metadata_metabolomics %>% dplyr::select(sample_ID, Cohort), by = c("SampleID" = "sample_ID")) %>%
  pivot_longer(cols = -c(SampleID, Cohort), names_to = "Feature", values_to = "Value") %>%
  dplyr::filter(Feature %in% nacyl_all$Feature) %>% 
  group_by(Feature) %>%
  summarise(p_value = wilcox.test(Value ~ Cohort)$p.value) %>%
  mutate(adj_p_value = p.adjust(p_value, method = "BH"))

volcano_data <- data_nacyl_log2fc %>%
  left_join(p_values_nacyl, by = "Feature") %>%
  mutate(
    neg_log10_p = -log10(adj_p_value),
    significance = ifelse(adj_p_value < 0.05 & abs(Log2FC) > 1, "Yes", "No")) %>% 
  left_join(info_feature_complete) %>% 
  mutate(C_number = as.numeric(str_extract(Compound_Name, "(?<=C)\\d+(?=:)")),
         Type = case_when(
           C_number >= 2 & C_number <= 6  ~ "short",
           C_number >= 7 & C_number <= 12 ~ "medium",
           C_number >= 13 & C_number <= 30 ~ "long",
           TRUE ~ NA_character_
         )) %>%
  select(-C_number)

top_features_with_id <- volcano_data %>%
  filter(significance == "Yes")

volcano_plot <- ggplot(volcano_data, aes(x = Log2FC, y = neg_log10_p, color = Type)) +
  geom_point(alpha = 0.75, size = 3) +  
  scale_color_manual(values = c("long" = "#C9CBA3", "medium" = "#A31621", "short" = "#1F7A8C")) +
  geom_text_repel(data = top_features_with_id, aes(label = Feature), size = 3,
                  box.padding = 0.3, point.padding = 0.3, max.overlaps = 1000) +
  theme_minimal() +
  labs(title = "Volcano Plot - N acyl lipids", x = "log2FC", y = "-log10(adjusted P-value)") +
  theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        panel.grid.major = element_blank(),        
        panel.grid.minor = element_blank(),        
        axis.line = element_line(color = "black")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
  scale_x_continuous(limits = c(-10, 5)) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black")
#ggsave(filename = "VEOIBD_volcano_nacyl.svg", plot = volcano_plot, device = "svg", width = 4.5, height = 3.5, dpi = "retina")



##### DIPEPTIDES #####
# Only select discriminant dipeptides from the PLS-DA model
dipeptides_interest <- VIPs_cohort_Load %>% 
  dplyr::filter(str_detect(pattern = "-Phe|-Leu|-Ile|-Met|-Trp|ASP-PHE", Compound_Name)) %>% 
  dplyr::filter(!ID %in% c("8739", "9022"))

data_sum_TIC_di <- data.frame(TIC = rowSums(data_sample %>% column_to_rownames("SampleID"))) %>% 
    rownames_to_column("SampleID")

data_vip_dipeptides <- data_sample %>%
  dplyr::select("SampleID", dipeptides_interest$ID) %>%
  dplyr::mutate(PeakArea = rowSums(select(., dipeptides_interest$ID))) %>%
  left_join(data_sum_TIC_di, by = "SampleID") %>% 
  dplyr::mutate(Relative_abundance = PeakArea / TIC) %>% 
  dplyr::select(SampleID, Relative_abundance) %>%
  mutate(log_Relative_abundance = log(Relative_abundance)) %>% 
  left_join(metadata_metabolomics, by = c("SampleID" = "sample_ID"))

# Plot ratio dipeptides
plot_ratio_dipeptides <- data_vip_dipeptides %>%
  ggboxplot(y = "log_Relative_abundance", x = "Cohort", add = "jitter", size = 0.35, ylab = "Log(relative Abundance)",
            add.params = list(color = "Cohort", size = 1.5), legend = "none",
            fill = "Cohort", 
            palette = c("#6a3d9a", "#ff7f00"),
            alpha = 0.5,
            title = "Dipeptides") +
  stat_compare_means() +
  theme(plot.title = element_text(size = 8),axis.title = element_text(size = 9),
        axis.text = element_text(size = 7))
#ggsave(filename = "plot_ratio_dipeptides.svg", plot = plot_ratio_dipeptides, device = "svg", width = 1.7, height = 3.5, dpi = "retina")



##### TRIPEPTIDES #####
tripeptides_interest <- VIPs_cohort_Load %>% 
  dplyr::filter(str_detect(pattern = "-Phe|-Leu|-Ile|-Met|-Trp|ASP-PHE", Compound_Name)) %>% 
  dplyr::filter(ID %in% c("9022", "8739")) %>% 
  dplyr::select(ID, comp1, GroupContrib, Compound_Name) %>% 
  mutate(Type = "Tri")

data_sum_TIC_tri <- data.frame(TIC = rowSums(data_sample %>% column_to_rownames("SampleID"))) %>% 
  rownames_to_column("SampleID")

data_vip_tripeptides <- data_sample %>%
  dplyr::select("SampleID", tripeptides_interest$ID) %>%
  dplyr::mutate(PeakArea = rowSums(select(., tripeptides_interest$ID))) %>%
  left_join(data_sum_TIC_tri, by = "SampleID") %>% 
  dplyr::mutate(Relative_abundance = PeakArea / TIC) %>% 
  dplyr::select(SampleID, Relative_abundance) %>%
  mutate(log_Relative_abundance = log(Relative_abundance)) %>% 
  left_join(metadata_metabolomics, by = c("SampleID" = "sample_ID"))

# Plot ratio dipeptides
plot_ratio_tripeptides <- data_vip_tripeptides %>%
  ggboxplot(y = "log_Relative_abundance", x = "Cohort", add = "jitter", size = 0.35, ylab = "Log(relative Abundance)",
            add.params = list(color = "Cohort", size = 1.5), legend = "none",
            fill = "Cohort", 
            palette = c("#6a3d9a", "#ff7f00"),
            alpha = 0.5,
            title = "Tripeptides") +
  stat_compare_means() +
  theme(plot.title = element_text(size = 8),axis.title = element_text(size = 9),
        axis.text = element_text(size = 7))
#ggsave(filename = "plot_ratio_tripeptides.svg", plot = plot_ratio_tripeptides, device = "svg", width = 1.7, height = 3.5, dpi = "retina")

# Include features predicted to be tripeptides by SIRIUS/Canopus
canopus <- read.delim("canopus_formula_summary.tsv")
canopus$mappingFeatureId <- as.character(canopus$mappingFeatureId)

class_predictions <- VIPs_cohort_Load %>% left_join(canopus, by = c("ID" = "mappingFeatureId")) %>% 
  dplyr::select(ID, comp1, GroupContrib, Compound_Name, NPC.pathway, NPC.pathway.Probability, NPC.superclass,
                NPC.superclass.Probability, NPC.class, NPC.class.Probability) %>% 
  dplyr::mutate(NPC.class = ifelse(NPC.class.Probability >= 0.7, NPC.class, NA))

tripeptides_canopus <- class_predictions %>% 
  dplyr::filter(str_detect(pattern = "Tripe", NPC.class)) %>% 
  dplyr::select(ID, comp1, GroupContrib, Compound_Name) %>% 
  mutate(Type = "Tri")

tripeptides_interest_filter <- tripeptides_interest %>% 
  dplyr::filter(ID == "9022")

tripeptides_all <- rbind(tripeptides_canopus, tripeptides_interest_filter)

data_sum_TIC_tri <- data.frame(TIC = rowSums(data_sample %>% column_to_rownames("SampleID"))) %>% 
  rownames_to_column("SampleID")

data_vip_tripeptides_all <- data_sample %>%
  dplyr::select("SampleID", tripeptides_all$ID) %>%
  dplyr::mutate(PeakArea = rowSums(select(., tripeptides_all$ID))) %>%
  left_join(data_sum_TIC_tri, by = "SampleID") %>% 
  dplyr::mutate(Relative_abundance = PeakArea / TIC) %>% 
  dplyr::select(SampleID, Relative_abundance) %>%
  mutate(log_Relative_abundance = log(Relative_abundance)) %>% 
  left_join(metadata_metabolomics, by = c("SampleID" = "sample_ID"))

# Plot ratio SIRIUS predicted tripeptides
plot_ratio_tripeptides_all <- data_vip_tripeptides_all %>%
  ggboxplot(y = "log_Relative_abundance", x = "Cohort", add = "jitter", size = 0.35, ylab = "Log(relative Abundance)",
            add.params = list(color = "Cohort", size = 1.5), legend = "none",
            fill = "Cohort", 
            palette = c("#6a3d9a", "#ff7f00"),
            alpha = 0.5,
            title = "Tripeptides all") +
  stat_compare_means() +
  theme(plot.title = element_text(size = 8),axis.title = element_text(size = 9),
        axis.text = element_text(size = 7))
#ggsave(filename = "plot_ratio_tripeptides_all_SIRIUS.svg", plot = plot_ratio_tripeptides_all, device = "svg", width = 1.7, height = 3.5, dpi = "retina")



##### BILE ACIDS #####

# Pof molecular networking validation of bile acid candidates
BA_query <- read.delim("BA_query_library_final.tsv")
names(BA_query)[names(BA_query) == "X.Scan."] <- "ID"
BA_query_filter <- BA_query %>% 
  dplyr::filter(str_detect(pattern = "Did not pass", query_validation))

ba_all <- info_feature_complete %>% 
  dplyr::filter(str_detect(pattern = "bile|cholic|-UDCA|-CA|-HDCA|-DCA|GLYCOCHENODEOXYCHOLIC|Tauro|tauro|TAUROCHENODEOXYCHOLIC",
                           Compound_Name)) %>% 
  dplyr::filter(!(Feature %in% c("8825", "28614", "23482")))

# Check bile acids ratio associated with one or the other group
# Only select discriminant BAs - PLS-DA model
ba_interst <- VIPs_cohort_Load %>% 
  dplyr::filter(str_detect(pattern = "bile|cholic|-UDCA|-CA|-HDCA|GLYCOCHENODEOXYCHOLIC", Compound_Name)) %>% 
  dplyr::filter(!(ID=="8825"))

# Identify BAs that did not pass the validation
matched_ids <- BA_query_filter[BA_query_filter$ID %in% ba_interst$ID, ] %>% 
  dplyr::filter(!ID %in% c(32275, 33374, 21568, 21571, 33375)) # these are true BAs

# Remove features annotated as BAs that did not pass validation
ba_interst_filter <- ba_interst %>% 
  dplyr::filter(!ID %in% matched_ids$ID)

# Selected features based on VIP
vip_healthy_ba <- VIPs_cohort_Load %>% dplyr::filter(GroupContrib == "Healthy") %>% 
  dplyr::filter(ID %in% ba_interst_filter$ID)
vip_ibd_ba <- VIPs_cohort_Load %>% dplyr::filter(GroupContrib == "VEOIBD") %>% 
  dplyr::filter(ID %in% ba_interst_filter$ID)

#vip_all <- rbind(vip_healthy_ba, vip_ibd_ba)
#write_csv(vip_all, "VEOIBD_vip_all_ba.csv")

data_vip_ba <- data_sample %>%
  dplyr::select("SampleID", vip_healthy_ba$ID, vip_ibd_ba$ID) %>%
  dplyr::mutate(Healthy = rowSums(select(., vip_healthy_ba$ID))) %>%
  dplyr::mutate(IBD = rowSums(select(., vip_ibd_ba$ID))) %>%
  dplyr::mutate(Ratio = log(IBD/Healthy)) %>%
  dplyr::select(SampleID, Ratio) %>%
  left_join(metadata_metabolomics, by = c("SampleID" = "sample_ID"))

# Plot ratio bile acids
plot_ratio_ba <- data_vip_ba %>%
  ggboxplot(y = "Ratio", x = "Cohort", add = "jitter", size = 0.35, ylab = "Log(VEO-IBD/Healthy)",
            add.params = list(color = "Cohort", size = 1.5), legend = "none",
            fill = "Cohort", 
            palette = c("#6a3d9a", "#ff7f00"),
            alpha = 0.5,
            title = "Bile Acid Ratio") +
  stat_compare_means() +
  theme(plot.title = element_text(size = 8),axis.title = element_text(size = 9),
        axis.text = element_text(size = 7))

#ggsave(filename = "plot_ratio_ba_updated.svg", plot = plot_ratio_ba, device = "svg", width = 1.7, height = 3.5, dpi = "retina")

plot_cor_ba_age <- data_vip_ba %>%
  ggscatter(x = "Age", y = "Ratio", add = "reg.line", color = "Cohort", legend = "none",
            alpha=0.7, ylab = "Log(VEO-IBD/Healthy)") +
  stat_cor(aes(color = Cohort), method = "pearson") + 
  scale_color_manual(
    values = c("Healthy" = "#6a3d9a", "VEOIBD" = "#ff7f00")) +
  scale_x_continuous(breaks = c(3, 6, 9, 12), limits = c(1, 14),
                     labels = c("3", "6", "9", "12"))
#ggsave(filename = "VEOIBD_plot_cor_ba_age.svg", plot = plot_cor_ba_age, device = "svg", width = 3, height = 4, dpi = "retina")

# LM to predict ratio based on age
model_lm <- data_vip_ba %>% 
  dplyr::mutate_at("Cohort", as.factor) %>%
  dplyr::mutate_at("Sex", as.factor) %>%
  lm(formula = Ratio ~ Age + Cohort)

summary(model_lm)

# GLM with age as confounding factor
model_glm <- data_vip_ba %>% 
  dplyr::mutate_at("Cohort", as.factor) %>%
  dplyr::mutate_at("Sex", as.factor) %>%
  glm(formula = Cohort ~ Ratio + Age, family = "binomial")

summary(model_glm)

# Look at the primary BA (CA and CDCA) and Gly and Tau conjugates
cholic_acid <- ba_all %>% 
  dplyr::filter(Feature == "20570")

DCA <- ba_all %>% 
  dplyr::filter(Feature == "35038")

GCA <- ba_all %>% 
  dplyr::filter(Feature == "22666")

GCDCA <- ba_all %>% 
  dplyr::filter(Feature == "29237")

TCA <- ba_all %>% 
  dplyr::filter(Feature == "19768")

TCDCA <- ba_all %>% 
  dplyr::filter(Feature == "23929")

data_primary_ba_sum <- data_sample %>% 
  dplyr::select("SampleID", cholic_acid$Feature, DCA$Feature, GCA$Feature, GCDCA$Feature, TCA$Feature, TCDCA$Feature) %>%
  dplyr::mutate(cholic_acid_sum = rowSums(dplyr::select(., cholic_acid$Feature))) %>% 
  dplyr::mutate(DCA_sum = rowSums(dplyr::select(., DCA$Feature))) %>% 
  dplyr::mutate(GCA_sum = rowSums(dplyr::select(., GCA$Feature))) %>% 
  dplyr::mutate(GCDCA_sum = rowSums(dplyr::select(., GCDCA$Feature))) %>%
  dplyr::mutate(TCA_sum = rowSums(dplyr::select(., TCA$Feature))) %>%
  dplyr::mutate(TCDCA_sum = rowSums(dplyr::select(., TCDCA$Feature))) %>%
  dplyr::select(SampleID, cholic_acid_sum, DCA_sum, GCA_sum, GCDCA_sum, TCA_sum, TCDCA_sum) %>%
  dplyr::mutate(Log_cholic_acid = log(cholic_acid_sum+1)) %>%
  dplyr::mutate(Log_DCA = log(DCA_sum+1)) %>%
  dplyr::mutate(Log_GCA = log(GCA_sum+1)) %>%
  dplyr::mutate(Log_GCDCA = log(GCDCA_sum+1)) %>%
  dplyr::mutate(Log_TCA = log(TCA_sum+1)) %>%
  dplyr::mutate(Log_TCDCA = log(TCDCA_sum+1)) %>%
  left_join(metadata_metabolomics, by = c("SampleID" = "sample_ID"))

data_primary_ba_sum_long <- data_primary_ba_sum %>%
  dplyr::select(SampleID, Log_cholic_acid, Log_DCA, Log_GCA, Log_GCDCA, Log_TCA, Log_TCDCA) %>%
  pivot_longer(
    cols = -SampleID,
    names_to = "Feature",
    values_to = "Abundance"
  ) %>%
  left_join(metadata, by = c("SampleID" = "sample_ID"))

p_values <- compare_means(Abundance ~ Cohort, group.by = "Feature", 
                          data = data_primary_ba_sum_long, method = "wilcox.test")

plot_primary_ba_conjugated <- ggboxplot(data_primary_ba_sum_long, x = "Feature", y = "Abundance", 
                                add = "jitter", ylab = "Log(sum Peak areas)", legend = "none",
                                add.params = list(color = "Cohort"),
                                fill = "Cohort", 
                                palette = c("#6a3d9a", "#ff7f00"),
                                alpha = 0.5,
                                title = "Conjugated BA") +
  theme(plot.title = element_text(size = 8),axis.title = element_text(size = 7),
        axis.text = element_text(size = 6))

#ggsave(filename = "plot_primary_ba_conjugated.svg", plot = plot_primary_ba_conjugated, device = "svg", width = 6, height = 3, dpi = "retina")



##### OXO BAs ######

keto_ba <- ba_all %>% 
  dplyr::filter(str_detect(Compound_Name, regex("oxo", ignore_case = TRUE)))

keto_3712oxo <- ba_all %>% dplyr::filter(Feature == "21568")
keto_37oxo12oxo <- ba_all %>% dplyr::filter(Feature == "21571")

data_keto_ba_sum <- data_final_clr %>%
  rownames_to_column() %>%
  dplyr::select(rowname, '21568', '21571') %>%
  pivot_longer(
    cols = -rowname,
    names_to = "Feature",
    values_to = "rclr_Peak_area"
  ) %>%
  left_join(metadata, by = c("rowname" = "sample_ID"))

p_values <- compare_means(rclr_Peak_area ~ Cohort, group.by = "Feature", 
                          data = data_keto_ba_sum, method = "wilcox.test")

plot_data_keto_ba <- ggboxplot(data_keto_ba_sum, x = "Feature", y = "rclr_Peak_area", size = 0.35,
                               add = "jitter", ylab = "rclr Peak area", legend = "none",
                               add.params = list(color = "Cohort", size = 1.5),
                               fill = "Cohort", 
                               palette = c("#6a3d9a", "#ff7f00"),
                               alpha = 0.5,
                               title = "Oxo BA") +
  theme(plot.title = element_text(size = 8),axis.title = element_text(size = 7),
        axis.text = element_text(size = 6))

#ggsave(filename = "plot_data_keto_ba1.svg", plot = plot_data_keto_ba, device = "svg", width = 2.8, height = 3.5, dpi = "retina")

keto_ba_1 <- ba_all %>% dplyr::filter(Feature == "19831")
keto_ba_2 <- ba_all %>% dplyr::filter(Feature == "24594")
keto_ba_3 <- ba_all %>% dplyr::filter(Feature == "32275")
keto_ba_4 <- ba_all %>% dplyr::filter(Feature == "33374")

data_keto_ba_sum_SI <- data_final_clr %>%
  rownames_to_column() %>%
  dplyr::select(rowname, '33374', '19831', '24594', '32275') %>%
  pivot_longer(
    cols = -rowname,
    names_to = "Feature",
    values_to = "rclr_Peak_area"
  ) %>%
  left_join(metadata, by = c("rowname" = "sample_ID"))

p_values <- compare_means(rclr_Peak_area ~ Cohort, group.by = "Feature", 
                          data = data_keto_ba_sum_SI, method = "wilcox.test")

plot_data_keto_ba_SI <- ggboxplot(data_keto_ba_sum_SI, x = "Feature", y = "rclr_Peak_area", size = 0.35,
                               add = "jitter", ylab = "rclr Peak area", legend = "none",
                               add.params = list(color = "Cohort", size = 1.5),
                               fill = "Cohort", 
                               palette = c("#6a3d9a", "#ff7f00"),
                               alpha = 0.5,
                               title = "Oxo BA") +
  theme(plot.title = element_text(size = 8),axis.title = element_text(size = 7),
        axis.text = element_text(size = 6))

#ggsave(filename = "plot_data_keto_ba_SI.svg", plot = plot_data_keto_ba_SI, device = "svg", width = 4, height = 3.5, dpi = "retina")

glu_CA <- data_final_clr %>% rownames_to_column("sample_ID") %>% 
  dplyr::select(sample_ID, '21984') %>% 
  left_join(metadata_metabolomics, by = "sample_ID")
colnames(glu_CA)[2] <- "feat"

plot_glu_CA <- glu_CA %>%
  ggboxplot(y = "feat", x = "Cohort", add = "jitter", size = 0.35, ylab = "rclr Peak area",
            add.params = list(color = "Cohort", size = 1.5), legend = "none",
            fill = "Cohort", 
            palette = c("#6a3d9a", "#ff7f00"),
            alpha = 0.5,
            title = "Glu-CA") +
  stat_compare_means() +
  theme(plot.title = element_text(size = 8),axis.title = element_text(size = 9),
        axis.text = element_text(size = 7))

#ggsave(filename = "plot_glu_CA.svg", plot = plot_glu_CA, device = "svg", width = 1.7, height = 3.5, dpi = "retina")



##### BAR PLOT - OTHER FEATURES DRIVING SEPARATION #####
data_other_feat <- VIPs_cohort_Load %>% 
  dplyr::filter(ID %in% c("5185", "5158", "4622", "7191", "37125", "4744",
                          "21984", "9857", "6895", "8673", "5518", "16870", "10060", "15321"))


data_other_feat_sum <- data_sample %>% 
  column_to_rownames("SampleID") %>% 
  dplyr::select(where(~ sum(.) != 0)) %>%
  dplyr::mutate(row_sum = rowSums(.)) %>%
  dplyr::mutate(across(everything(), ~ . / row_sum)) %>%
  dplyr::select(-row_sum) %>% 
  rownames_to_column("SampleID") %>%
  left_join(metadata_metabolomics %>% dplyr::select(sample_ID, Cohort), by = c("SampleID" = "sample_ID")) %>% 
  dplyr::select(-SampleID) %>% group_by(Cohort) %>%
  summarise(across(everything(), mean))

data_fold <- data_other_feat_sum %>%
  pivot_longer(cols = -Cohort, names_to = "Feature", values_to = "Value") %>%
  pivot_wider(names_from = Cohort, values_from = Value) %>%
  mutate(across(c("VEOIBD", "Healthy"), ~ifelse(. == 0, 1e-9, .))) %>%
  mutate(Fold_Change = VEOIBD/Healthy) %>%
  mutate(Log2FC = log2(Fold_Change)) %>%
  arrange(desc(Log2FC)) %>% 
  dplyr::filter(Feature %in% data_other_feat$ID) %>%
  left_join(info_feature_complete) %>%
  dplyr::filter(abs(Log2FC) > 1)

data_ra <- data_sample %>% 
  column_to_rownames("SampleID") %>% 
  dplyr::select(where(~ sum(.) != 0)) %>%
  dplyr::mutate(row_sum = rowSums(.)) %>%
  dplyr::mutate(across(everything(), ~ . / row_sum)) %>%
  dplyr::select(-row_sum) %>%
  dplyr::select(data_other_feat$ID) %>%
  #dplyr::mutate(across(everything(), ~ log(. + 1))) %>%
  rownames_to_column("SampleID") %>%
  left_join(metadata_metabolomics %>% dplyr::select("sample_ID", "Cohort"), by = c("SampleID" = "sample_ID")) %>%
  dplyr::select(-SampleID) %>%
  rename_with(~ paste0("Feature_", .), -Cohort) 

data_ra_long <- data_ra %>%
  pivot_longer(cols = -Cohort, names_to = "Feature", values_to = "Value") # Convert to long format

# Wilcoxon test for each feature
wilcox_results <- data_ra_long %>%
  group_by(Feature) %>%
  summarise( p_value = wilcox.test(Value[Cohort == "Healthy"], Value[Cohort == "VEOIBD"], exact = FALSE)$p.value,
             .groups = "drop") %>%
  dplyr::mutate(adj_p_value = p.adjust(p_value, method = "BH")) %>%
  dplyr::mutate(Feature = gsub("Feature_", "", Feature)) %>%
  left_join(data_fold %>% dplyr::select(Feature, Log2FC, Compound_Name))

VEOIBD_barplot <- ggplot(wilcox_results, aes(x = reorder(Feature, Log2FC), y = Log2FC)) +
  geom_bar(stat = "identity", fill = "#D55E00", width = 0.7) + 
  theme_minimal() +
  coord_flip() +
  labs(
    title = "Log2 Fold Change for Features",
    x = "Feature",
    y = "Log2 Fold Change"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"), axis.title = element_text(size = 14),
    axis.text = element_text(size = 12), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"), axis.ticks = element_line(color = "black"),
    axis.ticks.length = unit(0.15, "cm"))

#ggsave(filename = "VEOIBD_barplot.svg", plot = VEOIBD_barplot, device = "svg", width = 3, height = 3.5, dpi = "retina")


##### Boxplot sulfasalazine #####


sulfa <- data_sample %>%
  dplyr::select(SampleID, `21282`)
colnames(sulfa)[2] <- "feat"
sulfa$SampleID <- as.character(sulfa$SampleID)
sulfa$sulfa_use_known <- ifelse(sulfa$SampleID %in% c("PT_S8", "PT_AC_020", "PT_AC_001", "PT_AC_175"), "yes", "no")
sulfa <- dplyr::left_join(sulfa, metadata_metabolomics, by = c("SampleID" = "sample_ID")) 

# N-Acetyl Mesalazine 21282
sulfa_feat <- data_final_clr %>% rownames_to_column("sample_ID") %>% 
  dplyr::select(sample_ID, '21282') %>% 
  left_join(metadata_metabolomics, by = "sample_ID")
colnames(sulfa_feat)[2] <- "feat"
plot_sulfa <- sulfa_feat %>%
  ggboxplot(y = "feat", x = "Cohort", add = "jitter", ylab = "rclr",
            fill = "Cohort",
            add.params = list(color = "Cohort", alpha = 0.7), legend = "none",
            palette = c("#6a3d9a", "#ff7f00"),
            alpha = 0.5,
            title = "Sulfasalazine") +
  stat_compare_means() + 
  theme(plot.title = element_text(size = 10), axis.title.y = element_text(size = 8), 
        axis.text.y = element_text(size = 6), axis.title.x = element_blank())

#ggsave(filename = "VEOIBD_plot_sulfasalazine.svg", plot = plot_sulfa, device = "svg", width = 1.5, height = 3, dpi = "retina")



##############################################################################################
# Export mgf for SIRIUS predictions
vip_200 <- VIPs_cohort_Load %>% head(250)

# Import full mgf and filter it
# mgf for microbeMASST 
# mgf for Sirius
dda <- Spectra("metabolomics/mzmine/gnps.mgf", source = MsBackendMgf())
dda_ids <- data.frame(ID = dda@backend@spectraData@listData$FEATURE_ID) %>%
  dplyr::mutate(Interest = ID %in% vip_200$ID)
dda_filtered <- dda[dda_ids$Interest]
export(dda_filtered, MsBackendMgf(), file = "ibd_interest_200.mgf", exportTitle = FALSE)
###############################################################################################

##### SIRIUS PREDICTIONS #####
canopus <- read.delim("canopus_formula_summary.tsv")
canopus$mappingFeatureId <- as.character(canopus$mappingFeatureId)

class_predictions <- VIPs_cohort_Load %>% left_join(canopus, by = c("ID" = "mappingFeatureId")) %>% 
  dplyr::select(ID, comp1, GroupContrib, Compound_Name, NPC.pathway, NPC.pathway.Probability, NPC.superclass,
                NPC.superclass.Probability, NPC.class, NPC.class.Probability) 
  
pathway_predictions_filtered <- class_predictions %>% 
  dplyr::mutate(NPC.pathway = ifelse(NPC.pathway.Probability >= 0.7, NPC.pathway, NA)) %>%
  dplyr::select(ID, comp1, GroupContrib, Compound_Name, NPC.pathway) %>% 
  group_by(NPC.pathway) %>% summarise(count = n())
#write_csv(pathway_predictions_filtered, "pathway_predictions_filtered_SIRIUS.csv")

superclass_predictions_filtered <- class_predictions %>% 
  dplyr::mutate(NPC.superclass = ifelse(NPC.superclass.Probability >= 0.7, NPC.superclass, NA)) %>% 
  dplyr::select(ID, comp1, GroupContrib, Compound_Name, NPC.superclass) %>% 
  group_by(NPC.superclass) %>% summarise(count = n())
#write_csv(superclass_predictions_filtered, "superclass_predictions_filtered_SIRIUS.csv")
class_predictions_filtered <- class_predictions %>% 
  dplyr::mutate(NPC.class = ifelse(NPC.class.Probability >= 0.7, NPC.class, NA)) %>%
  dplyr::select(ID, comp1, GroupContrib, Compound_Name, NPC.class) %>% 
  group_by(NPC.class) %>% summarise(count = n())

### Make volcano plot based on SIRIUS predictions - class

vip_unknown <- VIPs_cohort_Load %>% 
  filter(is.na(Compound_Name)) %>% 
  filter(comp1 > 2.5) %>% 
  left_join(canopus, by = c("ID" = "mappingFeatureId")) %>% 
  mutate(NPC.class = if_else(NPC.class.Probability < 0.70, NA_character_, NPC.class))

colnames(vip_unknown)[colnames(vip_unknown) == "ID"] <- "Feature"

data_sum_all_vip_unknown <- data_sample %>% 
  column_to_rownames("SampleID") %>% 
  dplyr::select(where(~ sum(.) != 0)) %>%
  dplyr::mutate(row_sum = rowSums(.)) %>%
  dplyr::mutate(across(everything(), ~ . / row_sum)) %>%
  dplyr::select(-row_sum) %>% 
  rownames_to_column("SampleID") %>%
  left_join(metadata_metabolomics %>% dplyr::select(sample_ID, Cohort), by = c("SampleID" = "sample_ID")) %>% 
  dplyr::select(-SampleID) %>% 
  group_by(Cohort) %>%
  summarise(across(everything(), mean)) 

data_sirius_log2fc <- data_sum_all_vip_unknown %>%
  pivot_longer(cols = -Cohort, names_to = "Feature", values_to = "Value") %>% 
  dplyr::filter(Feature %in% vip_unknown$Feature) %>%
  pivot_wider(names_from = Cohort, values_from = Value) %>%
  mutate(across(c("VEOIBD", "Healthy"), ~ifelse(. == 0, 1e-8, .))) %>%
  mutate(Fold_Change = VEOIBD / Healthy) %>%
  mutate(Log2FC = log2(Fold_Change)) %>%
  arrange(desc(Log2FC))

p_values_sirius <- data_sample %>%
  column_to_rownames("SampleID") %>% 
  dplyr::select(where(~ sum(.) != 0)) %>%
  dplyr::mutate(row_sum = rowSums(.)) %>%
  dplyr::mutate(across(everything(), ~ . / row_sum)) %>%
  dplyr::select(-row_sum) %>% 
  rownames_to_column("SampleID") %>%
  left_join(metadata_metabolomics %>% dplyr::select(sample_ID, Cohort), by = c("SampleID" = "sample_ID")) %>%
  pivot_longer(cols = -c(SampleID, Cohort), names_to = "Feature", values_to = "Value") %>%
  dplyr::filter(Feature %in% vip_unknown$Feature) %>% 
  group_by(Feature) %>%
  summarise(p_value = wilcox.test(Value ~ Cohort)$p.value) %>%
  mutate(adj_p_value = p.adjust(p_value, method = "BH"))

vip_pathway <- vip_unknown %>% select(Feature, NPC.class)

volcano_data <- data_sirius_log2fc %>%
  left_join(p_values_sirius, by = "Feature") %>%
  mutate(
    neg_log10_p = -log10(adj_p_value),
    significance = ifelse(adj_p_value < 0.05 & abs(Log2FC) > 1, "Yes", "No")) %>% 
  left_join(vip_pathway)

top_features_with_id <- volcano_data %>%
  filter(significance == "Yes")

volcano_plot <- ggplot(volcano_data, aes(x = Log2FC, y = neg_log10_p, color = NPC.class)) +
  geom_point(alpha = 0.75, size = 3) +  
  #scale_color_manual(values = c("long" = "#C9CBA3", "medium" = "#A31621", "short" = "#1F7A8C")) +
  theme_minimal() +
  labs(title = "Volcano Plot - SIRIUS pathway predictions", x = "log2FC", y = "-log10(adjusted P-value)") +
  theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        panel.grid.major = element_blank(),        
        panel.grid.minor = element_blank(),        
        axis.line = element_line(color = "black")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
  scale_x_continuous(limits = c(-10, 5)) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black")
#ggsave(filename = "VEOIBD_volcano_SIRIUS_pathway.svg", plot = volcano_plot, device = "svg", width = 4.5, height = 3.5, dpi = "retina")



############################################################################################################
##### TISSUE MASST SEARCH OUTPUT #####

# Unknown feature
table_feat_7024 <- matrix(c(15, 28, 75, 3024),
                          nrow = 2, byrow = TRUE,
                          dimnames = list(c("Presence", "Absence"), c("IBD", "Healthy")))

print(table_feat_7024)
chisq.test(table_feat_7024)

table_feat_7024 <- matrix(c(5, 28, 46, 3024),
                          nrow = 2, byrow = TRUE,
                          dimnames = list(c("Presence", "Absence"), c("IBD", "Healthy")))

print(table_feat_7024)
chisq.test(table_feat_7024)

### Feature 5158 - Glu-Phe ###
table_feat_5158 <- matrix(c(119, 167, 304, 2885),
                          nrow = 2, byrow = TRUE,
                          dimnames = list(c("Presence", "Absence"), c("IBD", "Healthy")))

print(table_feat_5158)
chisq.test(table_feat_5158)

data <- data.frame(
  Status = c("Presence", "Absence"),
  IBD = c(119, 304),
  Healthy = c(167, 2885)
)

# Convert to long format
data_long <- pivot_longer(data, cols = c("IBD", "Healthy"),
                          names_to = "Group", values_to = "Count")

# Plot as percentages
plot_Glu_Phe <- ggplot(data_long, aes(x = Group, y = Count, fill = Status)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_fill_manual(values = c("Presence" = "#21409A", "Absence" = "#FFDE17")) +
  labs(title = "Glu-Phe",
       x = "Group", y = "Percentage", fill = "Status") +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),  
    panel.grid.minor = element_blank(),  
    axis.line = element_line(color = "black", linewidth = 0.6), 
    axis.ticks = element_line(color = "black"))
#ggsave(filename = "Glu-Phe_tissuemasst.svg", plot = plot_Glu_Phe, device = "svg", width = 3, height = 3.5, dpi = "retina")


### Feature 4622 - Asn-Phe ###
table_feat_4622 <- matrix(c(52, 6, 423, 3046),
                          nrow = 2, byrow = TRUE,
                          dimnames = list(c("Presence", "Absence"), c("IBD", "Healthy")))

print(table_feat_4622)
chisq.test(table_feat_4622)

data <- data.frame(
  Status = c("Presence", "Absence"),
  IBD = c(52, 423),
  Healthy = c(6, 3046)
)

# Convert to long format
data_long <- pivot_longer(data, cols = c("IBD", "Healthy"),
                          names_to = "Group", values_to = "Count")

# Plot as percentages
plot_Asn_Phe <- ggplot(data_long, aes(x = Group, y = Count, fill = Status)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_fill_manual(values = c("Presence" = "#21409A", "Absence" = "#FFDE17")) +
  labs(title = "Asn-Phe ",
       x = "Group", y = "Percentage", fill = "Status") +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),   
    axis.line = element_line(color = "black", linewidth = 0.6),  
    axis.ticks = element_line(color = "black"))
#ggsave(filename = "Asn-Phe_tissuemasst.svg", plot = plot_Asn_Phe, device = "svg", width = 3, height = 3.5, dpi = "retina")


### Feature 8739 - Thr-Val-Leu ###
table_feat_8739 <- matrix(c(142, 161, 281, 2891),
                          nrow = 2, byrow = TRUE,
                          dimnames = list(c("Presence", "Absence"), c("IBD", "Healthy")))

print(table_feat_8739)
chisq.test(table_feat_8739)

data <- data.frame(
  Status = c("Presence", "Absence"),
  IBD = c(142, 281),
  Healthy = c(161, 2891)
)

# Convert to long format
data_long <- pivot_longer(data, cols = c("IBD", "Healthy"),
                          names_to = "Group", values_to = "Count")

# Plot as percentages
plot_Thr_Val_Leu <- ggplot(data_long, aes(x = Group, y = Count, fill = Status)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_fill_manual(values = c("Presence" = "#21409A", "Absence" = "#FFDE17")) +
  labs(title = "Thr-Val-Leu ",
       x = "Group", y = "Percentage", fill = "Status") +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),   
    panel.grid.minor = element_blank(),   
    axis.line = element_line(color = "black", linewidth = 0.6), 
    axis.ticks = element_line(color = "black"))
#ggsave(filename = "Thr-Val-Leu_tissuemasst.svg", plot = plot_Thr_Val_Leu, device = "svg", width = 3, height = 3.5, dpi = "retina")


### Feature 9022 - Ile-Pro-Ile ###
table_feat_9022 <- matrix(c(35, 42, 388, 3010),
                          nrow = 2, byrow = TRUE,
                          dimnames = list(c("Presence", "Absence"), c("IBD", "Healthy")))

print(table_feat_9022)
chisq.test(table_feat_9022)

data <- data.frame(
  Status = c("Presence", "Absence"),
  IBD = c(35, 388),
  Healthy = c(42, 3010)
)

# Convert to long format
data_long <- pivot_longer(data, cols = c("IBD", "Healthy"),
                          names_to = "Group", values_to = "Count")

# Plot as percentages
plot_Ile_Pro_Ile <- ggplot(data_long, aes(x = Group, y = Count, fill = Status)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_fill_manual(values = c("Presence" = "#21409A", "Absence" = "#FFDE17")) +
  labs(title = "Ile-Pro_Ile",
       x = "Group", y = "Percentage", fill = "Status") +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),   
    panel.grid.minor = element_blank(),   
    axis.line = element_line(color = "black", linewidth = 0.6),  
    axis.ticks = element_line(color = "black"))
#ggsave(filename = "Ile-Pro-Ile_tissuemasst.svg", plot = plot_Ile_Pro_Ile, device = "svg", width = 3, height = 3.5, dpi = "retina")


### Feature 9471 - m/z 360.2129 ###
table_feat_9471 <- matrix(c(46, 104, 180, 2948),
                          nrow = 2, byrow = TRUE,
                          dimnames = list(c("Presence", "Absence"), c("IBD", "Healthy")))

print(table_feat_9471)
chisq.test(table_feat_9471)

data <- data.frame(
  Status = c("Presence", "Absence"),
  IBD = c(46, 180),
  Healthy = c(104, 2948)
)

# Convert to long format
data_long <- pivot_longer(data, cols = c("IBD", "Healthy"),
                          names_to = "Group", values_to = "Count")

# Plot as percentages
plot_9471 <- ggplot(data_long, aes(x = Group, y = Count, fill = Status)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_fill_manual(values = c("Presence" = "#21409A", "Absence" = "#FFDE17")) +
  labs(title = "Feature Presence in IBD vs Healthy (Proportion)",
       x = "Group", y = "Percentage", fill = "Status") +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),   
    panel.grid.minor = element_blank(),  
    axis.line = element_line(color = "black", linewidth = 0.6),  
    axis.ticks = element_line(color = "black")  
  )
#ggsave(filename = "Feat_9471_tissuemasst.svg", plot = plot_9471, device = "svg", width = 3, height = 3.5, dpi = "retina")


