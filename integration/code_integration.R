# Set working directory
setwd()

# Read packages
library(tidyverse)
library(ggpubr)
library(mixOmics)
library(caret)

# Read data 
metabolomics <- read_csv("table_metabolomics_clr.csv")
microbiome <- read_csv("table_microbiome_clr.csv")
metadata <- read_csv("metadata_VEOIBD.csv")

# Keep samples for which we have both metabolomics and microbiome data
id_sample <- metabolomics %>% dplyr::select(SampleID) %>%
  inner_join(microbiome %>% dplyr::select(SampleID)) # HC_97 and HC_99 seems to cluster with IBD

metabolomics_filter <- metabolomics %>% dplyr::filter(SampleID %in% id_sample$SampleID) %>%
  arrange(SampleID) %>% select_at(vars(-one_of(nearZeroVar(., names = TRUE))))
microbiome_filter <- microbiome %>% dplyr::filter(SampleID %in% id_sample$SampleID) %>%
  arrange(SampleID) %>% select_at(vars(-one_of(nearZeroVar(., names = TRUE))))
metadata_filter <- metadata %>% dplyr::filter(sample_ID %in% id_sample$SampleID) %>%
  dplyr::rename(SampleID = sample_ID) %>%
  arrange(SampleID)

##########################
# Integration via DIABLO #
##########################

# PLS between the two datasets. Extract correlation structure between components for DIABLO design matrix 
pls_model <- mixOmics::pls(metabolomics_filter %>% column_to_rownames("SampleID"), 
                           microbiome_filter %>% column_to_rownames("SampleID"), 
                           ncomp = 2, scale = TRUE)
plotIndiv(pls_model, comp = 1:2, rep.space= 'XY-variate', group = metadata_filter$Cohort, ind.names = FALSE,
          legend = TRUE, title = "PLS comp 1 - 2, XY-space", pch = 20, style = "lattice", centroid = TRUE, ellipse = TRUE)
cor(pls_model$variates$X, pls_model$variates$Y) %>% diag() # correlation close to 0.9

# Prepare datasets
data <- list(Metabolomics = as.matrix(metabolomics_filter %>% column_to_rownames("SampleID")), 
             Microbiome = as.matrix(microbiome_filter %>% column_to_rownames("SampleID")))
Y <- metadata_filter$Cohort
design <- matrix(0.8, ncol = length(data), nrow = length(data), 
                 dimnames = list(names(data), names(data)))
diag(design) = 0

# Preliminary DIABLO 
DIABLO <- block.splsda(X = data, Y = Y, ncomp = 2, design = design)
perf_DIABLO <- perf(DIABLO, validation = 'loo')
plot(perf_DIABLO)  # 3 components and centroid distance

# Tune number of variables to retain per component
test_keepX <- list(Metabolomics = c(seq(10, 20, 2)),
                   Microbiome = c(seq(10, 20, 2)))

tune_DIABLO <- tune.block.splsda(X = data, Y = Y, ncomp = 2, 
                                 test.keepX = test_keepX, design = design,
                                 validation = 'loo', dist = "centroids.dist")

# Final DIABLO model
DIABLO <- block.splsda(X = data, Y = Y, ncomp = 2, 
                       keepX = tune_DIABLO$choice.keepX, design = design)
perf_DIABLO <- perf(DIABLO, validation = 'loo')
plotDiablo(DIABLO)

# Plot 
plotIndiv(DIABLO, ind.names = FALSE, legend = TRUE, title = 'DIABLO', col.per.group = c("#3A383F", "#85BEDC"),
          style = "lattice", ellipse = TRUE, centroid = TRUE, pch = 20, X.label = "Component 1", Y.label = "Component 2")

# Circos Plot - Evaluate correlation between variables in the two omics datasets
#library(circlize)
circosPlot(DIABLO, cutoff = 0.5, line = FALSE, size.labels = 0.5, comp = 1,
           color.blocks = c("#4F6D7A", "#3A383F"), showIntraLinks = FALSE, legend = TRUE,
           size.variables = 0.5, size.legend = 0.5)

svg("circosPlot_DIABLO.svg", width = 10, height = 10)  # Square dimensions
circosPlot(DIABLO, cutoff = 0.5, line = FALSE, size.labels = 0.5, comp = 1,
           color.blocks = c("#4F6D7A", "#56535C"), showIntraLinks = FALSE, legend = TRUE,
           size.variables = 0.5, size.legend = 0.5)
dev.off()



