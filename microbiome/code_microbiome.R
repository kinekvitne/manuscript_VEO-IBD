#Set working directory
setwd()

# Load in packages
library(tidyverse)
library(phyloseq)
library(vegan)
library(caret)
library(ggpubr)
library(mixOmics)
library(ALDEx2)
library(ANCOMBC)

# Read data
otu_table <- read_tsv("table_microbiome.tsv")
taxonomy_table <- read_tsv("taxonomy.tsv")
metadata <- read_csv("metadata_VEOIBD.csv")

# Fix taxonomy
taxonomy_table_split <- taxonomy_table %>% 
  separate(Taxon, into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep = ";") %>%
  column_to_rownames("Feature ID")

taxonomy_table_split$Kingdom <- gsub("k__", "", taxonomy_table_split$Kingdom)
taxonomy_table_split$Phylum <- gsub("p__", "", taxonomy_table_split$Phylum)
taxonomy_table_split$Class <- gsub("c__", "", taxonomy_table_split$Class)
taxonomy_table_split$Order <- gsub("o__", "", taxonomy_table_split$Order)
taxonomy_table_split$Family <- gsub("f__", "", taxonomy_table_split$Family)
taxonomy_table_split$Genus <- gsub("g__", "", taxonomy_table_split$Genus)
taxonomy_table_split$Species <- gsub("s__", "", taxonomy_table_split$Species)

taxonomy_table_split$Kingdom <- gsub("__", "", taxonomy_table_split$Kingdom)
taxonomy_table_split$Phylum <- gsub("__", "", taxonomy_table_split$Phylum)
taxonomy_table_split$Class <- gsub("__", "", taxonomy_table_split$Class)
taxonomy_table_split$Order <- gsub("__", "", taxonomy_table_split$Order)
taxonomy_table_split$Family <- gsub("__", "", taxonomy_table_split$Family)
taxonomy_table_split$Genus <- gsub("__", "", taxonomy_table_split$Genus)
taxonomy_table_split$Species <- gsub("__", "", taxonomy_table_split$Species)

taxonomy_table_split[taxonomy_table_split == "" | taxonomy_table_split == " "] <- NA

# Fill unknown taxa
taxonomy_table_fill <- taxonomy_table_split %>% dplyr::filter(!(is.na(Phylum)))

taxonomy_table_fill$Class[is.na(taxonomy_table_fill$Class)] <- as.character(taxonomy_table_fill$Phylum[is.na(taxonomy_table_fill$Class)])
taxonomy_table_fill$Order[is.na(taxonomy_table_fill$Order)] <- as.character(taxonomy_table_fill$Class[is.na(taxonomy_table_fill$Order)])
taxonomy_table_fill$Family[is.na(taxonomy_table_fill$Family)] <- as.character(taxonomy_table_fill$Order[is.na(taxonomy_table_fill$Family)])
taxonomy_table_fill$Genus[is.na(taxonomy_table_fill$Genus)] <- as.character(taxonomy_table_fill$Family[is.na(taxonomy_table_fill$Genus)])
taxonomy_table_fill$Species[is.na(taxonomy_table_fill$Species)] <- as.character(taxonomy_table_fill$Genus[is.na(taxonomy_table_fill$Species)])

# Fix OTU table
otu_table_t <- t(otu_table %>% column_to_rownames("OTU_ID")) %>% as.data.frame() %>%
  rownames_to_column("sample_ID") %>%
  dplyr::mutate(sample_ID = gsub("15748.M.", "", sample_ID)) %>%
  dplyr::mutate(sample_ID = gsub("\\.", "_", sample_ID)) %>% 
  column_to_rownames("sample_ID") %>% as.matrix()

metadata_phylo <- metadata

# Filter taxonomy 
aa <- colnames(otu_table_t) %>% as.data.frame() 
colnames(aa) <- "ID"
bb <- taxonomy_table_fill %>% rownames_to_column("ID")
cc <- aa %>% left_join(bb) %>% dplyr::filter(is.na(Kingdom))
dd <- otu_table %>% 
  dplyr::filter(OTU_ID %in% cc$ID) %>% 
  dplyr::mutate(RowSum = rowSums(across(where(is.numeric))))

otu_table_final <- otu_table_t %>% as.data.frame() %>% 
  dplyr::select(-(cc$ID)) %>% t() %>% as.matrix()
taxonomy_table_final <- taxonomy_table_fill %>% as.matrix()

metadata_ps <- data_frame(sample_ID = colnames(otu_table_final)) %>%
  left_join(metadata_phylo) %>% column_to_rownames("sample_ID")


# Generate PS object
ps <- phyloseq(otu_table(otu_table_final, taxa_are_rows = TRUE), sample_data(metadata_ps),
               tax_table(taxonomy_table_final))
ps_genus <- ps %>% tax_glom(taxrank = "Genus")


# Extract count table
table_genus <- as.data.frame(ps_genus@otu_table) %>% as.data.frame()
table_taxa <- as.data.frame(ps_genus@tax_table) %>% as.data.frame() %>%
  rownames_to_column("ID") %>% dplyr::select(ID, Family, Genus) %>%
  dplyr::mutate(taxa = paste(Family, Genus, sep = "_")) %>% 
  dplyr::mutate(taxa = gsub(" ", "", taxa)) %>%
  dplyr::select(ID, taxa)

# Check counts
table_genus_sum <- data_frame(sample_ID = colnames(table_genus), Abundance = colSums(table_genus)) %>%
  dplyr::filter(Abundance > 5000) %>%
  arrange(desc(Abundance))

table_genus_final <- table_genus %>% dplyr::select(table_genus_sum$sample_ID) %>%
  rownames_to_column("ID") %>% left_join(table_taxa) %>% dplyr::select(-ID) %>%
  column_to_rownames("taxa") %>% t() %>% as.data.frame() %>% rownames_to_column("sample_ID")

#write_csv(x = table_genus_final, file = "table_microbiome.csv")

# Robust CLR Transformation
OTU_clr <- table_genus_final %>% column_to_rownames("sample_ID") %>%
  decostand(method = "rclr")

#write_csv(x = OTU_clr %>% rownames_to_column("SampleID"), file = "table_microbiome_clr.csv")

# Remove zero variance
OTU_nvz <- OTU_clr %>%
  select_at(vars(-one_of(caret::nearZeroVar(., names = TRUE))))

# PCA
PCA_whole <- mixOmics::pca(OTU_nvz, ncomp = 2, center = TRUE, scale = TRUE)
PCA_whole_scores <- data.frame(PCA_whole$variates$X) %>% 
  rownames_to_column("sample_ID") %>% left_join(metadata)

i <- "Cohort"

PCA_plot <- PCA_whole_scores %>%
  ggscatter(x = "PC1", y = "PC2", color = i, alpha = 0.6, 
            title = paste("PCA - Fecal microbiome"),
            xlab = paste("PC1 (", round(PCA_whole$prop_expl_var$X[1]*100, digits = 1),"%)", sep = ""), 
            ylab = paste("PC2 (", round(PCA_whole$prop_expl_var$X[2]*100, digits = 1),"%)", sep = ""),
            ggtheme = theme_classic()) +
  geom_point(data = PCA_whole_scores %>% group_by((!!sym(i))) %>% 
               summarise_at(vars(matches("PC")), mean), aes(PC1, PC2, color = (!!sym(i))), size = 4, shape = 8) +
  scale_color_manual(values = c("VEOIBD" = "#ff7f00", "Healthy" = "#6a3d9a")) + 
  theme(plot.title = element_text(size = 10), axis.title = element_text(size = 8),
        axis.text = element_text(size = 4)) + coord_fixed()


# PERMANOVA
dist_clr <- vegdist(OTU_nvz, method = "euclidean")
disper <- betadisper(dist_clr, PCA_whole_scores$Cohort)
anova(disper)
permanova_cohort <- adonis2(dist_clr ~ Cohort * Age * Sex, PCA_whole_scores, na.action = na.omit, by = "terms")

#ggsave(filename = "VEOIBD_PCA_fecal_microbiome.svg", plot = PCA_plot, device = "svg", width = 4, height = 3, dpi = "retina")


# PLSDA
PLSDA_cohort <- mixOmics::plsda(OTU_nvz, PCA_whole_scores$Cohort, ncomp = 3, scale = TRUE)
PLSDA_cohort_scores <- data.frame(PLSDA_cohort$variates$X) %>% 
  rownames_to_column("sample_ID") %>% left_join(metadata)

PLSDA_cohort_plot <- PLSDA_cohort_scores %>%
  ggscatter(x = "comp1", y = "comp2", color = "Cohort", alpha = 0.6, title = "PLSDA - Microbiome",
            xlab = paste("Component 1 (", round(PLSDA_cohort$prop_expl_var$X[1]*100, digits = 1),"%)", sep = ""), 
            ylab = paste("Component 2 (", round(PLSDA_cohort$prop_expl_var$X[2]*100, digits = 1),"%)", sep = ""),
            legend.title = "Group", ggtheme = theme_classic()) +
  geom_point(data = PLSDA_cohort_scores %>% group_by(Cohort) %>% 
               summarise_at(vars(matches("comp")), mean), aes(comp1, comp2, color = Cohort), size = 3, shape = 8) +
  scale_color_manual(values = c("VEOIBD" = "#ff7f00", "Healthy" = "#6a3d9a")) +
  theme(plot.title = element_text(size = 10),axis.title = element_text(size = 8),
        axis.text = element_text(size = 5)) + coord_fixed()

#ggsave(filename = "VEOIBD_PLSDA_cohort_scoreplot_microbiome.svg", plot = PLSDA_cohort_plot, device = "svg", width = 5, height = 5.5, dpi = "retina")


Loadings_cohort <- plotLoadings(PLSDA_cohort, plot = FALSE, contrib = "max")$X %>% 
  as.data.frame() %>% rownames_to_column() %>% dplyr::select(rowname, GroupContrib)

perf_PLSDA_cohort <- perf(PLSDA_cohort,  validation = "Mfold", folds = 4, nrepeat = 99, progressBar = TRUE) 
plot(perf_PLSDA_cohort, legend = FALSE)

VIPs_cohort <- as.data.frame(mixOmics::vip(PLSDA_cohort))
VIPs_cohort_filter <- dplyr::filter(VIPs_cohort, VIPs_cohort$comp1 > 1)
VIPs_cohort_filter$ID <- rownames(VIPs_cohort_filter)
VIPs_cohort_select <- VIPs_cohort_filter %>% dplyr::select(ID, comp1)
VIPs_cohort_Load <- VIPs_cohort_select %>% 
  left_join(Loadings_cohort, by = c("ID" = "rowname")) %>% arrange(desc(comp1))

#write_csv(x = VIPs_cohort_Load, file = "VIPs_cohort_Load_microbiome.csv")


# Selected features based on VIP
vip_healthy <- VIPs_cohort_Load %>% dplyr::filter(GroupContrib == "Healthy") 
vip_ibd <- VIPs_cohort_Load %>% dplyr::filter(GroupContrib == "VEOIBD")

data_vip <- table_genus_final %>%
  dplyr::select("sample_ID", vip_healthy$ID, vip_ibd$ID) %>%
  dplyr::mutate(Healthy = rowSums(dplyr::select(., vip_healthy$ID))) %>%
  dplyr::mutate(IBD = rowSums(dplyr::select(., vip_ibd$ID))) %>%
  dplyr::mutate(Ratio = log(Healthy/IBD)) %>%
  dplyr::select(sample_ID, Ratio) %>%
  left_join(metadata)

data_vip$Cohort <- factor(data_vip$Cohort, levels = c("Healthy", "VEOIBD"))

# Plot ratio
plot_ratio <- data_vip %>%
  ggboxplot(y = "Ratio", x = "Cohort", add = "jitter", ylab = "Log(Healthy/IBD)",
            add.params = list(color = "Cohort"), legend = "none",
            fill = "Cohort", 
            palette = c("Healthy" = "#6a3d9a", "VEOIBD" = "#ff7f00"),
            alpha = 0.5,
            title = "Differential features from PLS-DA models") +
  stat_compare_means() +
  theme(plot.title = element_text(size = 8), axis.title = element_text(size = 7),
        axis.text = element_text(size = 6))

#ggsave(filename = "VEOIBD_ratio_features_PLSDA_microbiome.svg", plot = plot_ratio, device = "svg", width = 1.7, height = 3.5, dpi = "retina")


######################################
# Differential abundance with ALDEx2 #
######################################

all(table_genus_final$sample_ID == metadata$sample_ID)

metadata_aldex <- data_frame(sample_ID = table_genus_final$sample_ID) %>% 
  left_join(metadata)

all(table_genus_final$sample_ID == metadata_aldex$sample_ID)

# ALDEx2 Analysis
cond_cohort <- as.vector(metadata_aldex$Cohort)
ALDE_cohort <- aldex.clr(t(table_genus_final %>% column_to_rownames("sample_ID")), cond_cohort, denom = "all") 
ALDE_cohort_test <- aldex.ttest(ALDE_cohort)
ALDE_cohort_eff <- aldex.effect(ALDE_cohort, useMC = TRUE)

# Store results
aldex_output <- data.frame(rownames(ALDE_cohort_eff), ALDE_cohort_eff, ALDE_cohort_test)
aldex_output_filter <- aldex_output[which(aldex_output$wi.eBH < 0.1),]
colnames(aldex_output_filter)[1] <- "Name"

# Baloon plots
# Add a column to classify effect size as positive or negative
aldex_output_filter$Group <- ifelse(aldex_output_filter$effect > 0, "VEOIBD", "Healthy")

otu_aldex_plot <- ggdotchart(aldex_output_filter, x = "Name", y = "effect",
                             sorting = "ascending",  dot.size = 2, 
                             add = "segments", ylab = "Effect Size",                      
                             add.params = list(color = "lightgray", size = 0.1),
                             title = "OTUs - Healthy vs VEO-IBD",
  color = "Group", palette = c("Healthy" = "#6a3d9a", "VEOIBD" = "#ff7f00")) +
  geom_hline(yintercept = 0, linetype = 2, color = "lightgray") +
  coord_flip() + theme_bw() + theme(plot.title = element_text(size = 10),
                                    axis.title = element_text(size = 8),
                                    axis.text = element_text(size = 6),
                                    legend.title = element_text(size = 8),
                                    legend.text = element_text(size = 6),
                                    panel.grid.major = element_blank(), 
                                    panel.grid.minor = element_blank()) +
  theme(legend.position = "none")

#ggsave(filename = "VEOIBD_baloon_plot_microbiome.svg", plot = otu_aldex_plot, device = "svg", width = 3.5, height = 3, dpi = "retina")


########################################
# Differential abundance with ANCOM-BC #
########################################
set.seed(123)
ancom_output <- ancombc2(data = ps_genus, fix_formula = "Cohort + Age + Sex",
                         p_adj_method = "holm", pseudo_sens = TRUE, prv_cut = 0.10,
                         group = "Cohort", struc_zero = TRUE, neg_lb = TRUE,
                         alpha = 0.05, n_cl = 4, verbose = TRUE, global = TRUE,
                         pairwise = TRUE, dunnet = TRUE, trend = TRUE,
                         iter_control = list(tol = 1e-2, max_iter = 20, verbose = TRUE),
                         em_control = list(tol = 1e-5, max_iter = 100),
                         mdfdr_control = list(fwer_ctrl_method = "holm", B = 100),
                         trend_control = list(contrast = list(matrix(c(1, 0, -1, 1), nrow = 2, byrow = TRUE),
                                                              matrix(c(-1, 0, 1, -1), nrow = 2, byrow = TRUE),
                                                              matrix(c(1, 0, 1, -1), nrow = 2, byrow = TRUE)),
                                              node = list(2, 2, 1), solver = "ECOS", B = 100))

ancom_output_filter <- ancom_output$res %>%
  left_join(taxonomy_table_final %>% as.data.frame() %>% rownames_to_column("taxon")) %>%
  dplyr::filter(p_CohortVEOIBD < 0.05) %>%
  dplyr::filter(abs(lfc_CohortVEOIBD) > 0.5) %>%
  dplyr::select(34,35, 3:5, 15:17, 19:21) %>% 
  arrange(desc(lfc_CohortVEOIBD)) %>%
  dplyr::mutate(Name = paste(Family, Genus, sep = "_")) %>%
  dplyr::mutate(Name = gsub(" ", "", Name))

#write_csv(x = ancom_output_filter, file = "ancom_bc_output.csv")


# Intersect ALDEx2 and ANCOM-BC2
common_otu <- ancom_output_filter %>% full_join(aldex_output_filter, by = "Name")




###################################################################################
# Plot genera for Figure 2
genera_pls_da <- OTU_clr %>% rownames_to_column("sample_ID") %>%
  dplyr::select(sample_ID, Bifidobacteriaceae_Bifidobacterium, Lachnospiraceae_Blautia, 
                Lachnospiraceae_Lachnospira, Alcaligenaceae_Sutterella) %>%
  pivot_longer(
    cols = -sample_ID,
    names_to = "Feature",
    values_to = "rclr_Peak_area"
  ) %>%
  left_join(metadata) 

p_values <- compare_means(rclr_Peak_area ~ Cohort, group.by = "Feature", 
                          data = genera_pls_da, method = "wilcox.test")

plot_genera <- ggboxplot(genera_pls_da, x = "Feature", y = "rclr_Peak_area", size = 0.35,
                               add = "jitter", ylab = "rclr Peak area", legend = "none",
                               add.params = list(color = "Cohort", size = 1.5),
                               fill = "Cohort", 
                               palette = c("#6a3d9a", "#ff7f00"),
                               alpha = 0.5,
                               title = "Genera - PLS-DA") +
  theme(plot.title = element_text(size = 8),axis.title = element_text(size = 7),
        axis.text = element_text(size = 6))

#ggsave(filename = "plot_genera.svg", plot = plot_genera, device = "svg", width = 4.3, height = 3.5, dpi = "retina")

# Plot genera for SI Figure 2
genera_pls_da <- OTU_clr %>% rownames_to_column("sample_ID") %>%
  dplyr::select(sample_ID, Lachnospiraceae_Coprococcus, Coriobacteriaceae_Collinsella, Lachnospiraceae_Dorea, 
                Lachnospiraceae_Clostridium, Bacteroidaceae_Bacteroides, Rikenellaceae_Alistipes, Enterobacteriaceae_Enterobacteriaceae) %>%
  pivot_longer(
    cols = -sample_ID,
    names_to = "Feature",
    values_to = "rclr_Peak_area"
  ) %>%
  left_join(metadata) 

p_values <- compare_means(rclr_Peak_area ~ Cohort, group.by = "Feature", 
                          data = genera_pls_da, method = "wilcox.test")

plot_genera <- ggboxplot(genera_pls_da, x = "Feature", y = "rclr_Peak_area", size = 0.35,
                         add = "jitter", ylab = "rclr Counts", legend = "none",
                         add.params = list(color = "Cohort", size = 1.5),
                         fill = "Cohort", 
                         palette = c("#6a3d9a", "#ff7f00"),
                         alpha = 0.5,
                         title = "Genera - PLS-DA") +
  theme(plot.title = element_text(size = 8),axis.title = element_text(size = 7),
        axis.text = element_text(size = 6))

#ggsave(filename = "plot_genera_SI.svg", plot = plot_genera, device = "svg", width = 6, height = 3.5, dpi = "retina")



# Alpha diversity
ps_rar <- rarefy_even_depth(ps,
                            sample.size = 5000,
                            rngseed = 1234,
                            replace = FALSE,
                            trimOTUs = TRUE)

alpha_diversity <- ps_rar %>%
  estimate_richness(split = TRUE, measures = c("Observed", "Shannon", "Chao1")) %>%
  dplyr::mutate(sample_ID = sample_names(ps_rar))

alpha_diversity_info <- alpha_diversity %>% left_join(metadata)

alpha_diversity_info$Cohort <- factor(alpha_diversity_info$Cohort, levels = c("Healthy", "VEOIBD"))

alpha_diversity_shannon_plot <- alpha_diversity_info %>%
  ggboxplot(y = "Shannon", x = "Cohort", add = "jitter", ylab = "Shannon Diversity Index", 
            add.params = list(color = "Cohort"), xlab = FALSE, legend = "none",
            fill = "Cohort",
            palette = c("Healthy" = "#6a3d9a", "VEOIBD" = "#ff7f00"),
            alpha = 0.5,
            title = "Alpha Diversity") +
  stat_compare_means(method = "wilcox.test") +
  theme(plot.title = element_text(size = 10), axis.title = element_text(size = 8),
        axis.text = element_text(size = 6))

#ggsave(filename = "VEOIBD_alpha_diversity.svg", plot = alpha_diversity_shannon_plot, device = "svg", width = 1.7, height = 3.5, dpi = "retina")

# Alpha diversity over time
plot_alpha_time <- alpha_diversity_info %>%
  ggscatter(x = "Age", y = "Shannon", add = "reg.line", color = "Cohort", legend = "none",
            alpha=0.7, ylab = "Shannon Diversity Index", xlab = "Age (years)") +
  stat_cor(aes(color = Cohort), method = "pearson") + 
  scale_color_manual(
    values = c("Healthy" = "#6a3d9a", "VEOIBD" = "#ff7f00"))  +
  scale_x_continuous(breaks = c(3, 6, 9, 12), limits = c(1, 14),
                     labels = c("3", "6", "9", "12"))
#ggsave(filename = "plot_alpha_time.svg", plot = plot_alpha_time, device = "svg", width = 3.6, height = 4, dpi = "retina")


