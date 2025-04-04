# Banksy MULTI-SAMPLE CLUSTERING - qMSI data (RPCA)
# https://github.com/prabhakarlab/Banksy/blob/bioc/vignettes/multi-sample.Rmd
# Benedetta Manzato - 26-09-2024

install.packages("BiocManager")
BiocManager::install(version = "3.20")

BiocManager::install("SpatialExperiment")#,version = "3.18")
BiocManager::install("ggspavis")
install.packages('remotes')
remotes::install_github("prabhakarlab/Banksy", force = TRUE)
remotes::install_github("satijalab/seurat", "seurat5", force = TRUE)
install.packages('Seurat')
remotes::install_github("mojaveazure/seurat-disk")
remotes::install_version("Matrix", version = "1.6-4")
install.packages("GenomicRanges",version='2.37.1')

library(Banksy)
library(SeuratDisk)
library(Seurat)
library(SpatialExperiment)

wd <- "/exports/humgen/bmanzato/nieromics_dir"

############################# prep data

# read RDS

rpca_object <- readRDS(paste(wd,"/qMSI_data/seurat_rpca.rds",sep=""))
rpca_data_layer <- GetAssayData(rpca_object, slot = "data", assay = "integrated")
rpca_df <- as.data.frame(t(as.matrix(rpca_data_layer)))

write.csv(rpca_df, file = paste(wd,'/qMSI_data/rpca_df.csv',sep=""), row.names = TRUE)

ct_ann <- read.csv(paste(wd,"/rosalie/qmsi_analysis/ct_annotations.csv",sep=""),row.names = 1)
ct_ann <- subset(ct_ann, x != "gaps")




# SUBSET the features of the object to keep only the subset of lipids decided on 06-11-2024
subset_lipids <- read.csv(paste(wd,"/qMSI_data/subset_lipids_06112024.csv",sep=""))

target_columns <- as.character(subset_lipids$x)
rpca_df <- rpca_df[, colnames(rpca_df) %in% target_columns]
write.csv(rpca_df, file = paste(wd,'/qMSI_data/rpca_df_130lipids.csv',sep=""), row.names = TRUE)


# biri
coord_biri1 <- read.csv(paste(wd,"/qMSI_data/coord/20240720 mouse IRI kidney bIRI_1_xycoord.csv", sep = ""),sep=";",comment.char = "#",row.names = 1)
coord_biri2 <- read.csv(paste(wd,"/qMSI_data/coord/20240720 mouse IRI kidney bIRI_2_xycoord.csv", sep = ""),sep=";",comment.char = "#",row.names = 1)
coord_biri3 <- read.csv(paste(wd,"/qMSI_data/coord/20240720 mouse IRI kidney bIRI_3_xycoord.csv", sep = ""),sep=";",comment.char = "#",row.names = 1)

# sham
coord_sham1 <- read.csv(paste(wd,"/qMSI_data/coord/20240720 mouse IRI kidney sham_1_xycoord.csv", sep = ""),sep=";",comment.char = "#",row.names = 1)
coord_sham2 <- read.csv(paste(wd,"/qMSI_data/coord/20240720 mouse IRI kidney sham_2_xycoord.csv", sep = ""),sep=";",comment.char = "#",row.names = 1)
coord_sham3 <- read.csv(paste(wd,"/qMSI_data/coord/20240720 mouse IRI kidney sham_3_xycoord.csv", sep = ""),sep=";",comment.char = "#",row.names = 1)

# Add a new column to the coord data frames to indicate the source
coord_biri1$kidney <- 'biri1'
coord_biri2$kidney <- 'biri2'
coord_biri3$kidney <- 'biri3'
coord_sham1$kidney <- 'sham1'
coord_sham2$kidney <- 'sham2'
coord_sham3$kidney <- 'sham3'

# Concatenate 
coord_biri <- rbind(coord_biri1, coord_biri2, coord_biri3)
coord_sham <- rbind(coord_sham1, coord_sham2, coord_sham3)


# separate sham and biri
rpca_biri <- rpca_df[1:nrow(coord_biri), ]

start_row <- nrow(rpca_df) - nrow(coord_sham) + 1
rpca_sham <- rpca_df[start_row:nrow(rpca_df), ]

# rename rows
rownames(coord_biri) <- rownames(rpca_biri)
rownames(coord_sham) <- rownames(rpca_sham)

# filter rpca_biri and rpca_sham and coord to keep only spots != gaps
rpca_biri <- rpca_biri[rownames(rpca_biri) %in% rownames(ct_ann), ]
coord_biri <- coord_biri[rownames(coord_biri) %in% rownames(ct_ann), ]

rpca_sham <- rpca_sham[rownames(rpca_sham) %in% rownames(ct_ann), ]
coord_sham <- coord_sham[rownames(coord_sham) %in% rownames(ct_ann), ]

######################### BANKSY on bIRI tissues

# Create SpatialExperiment object; it requires two assays: counts and logcounts
spe <- SpatialExperiment(
  assays = list(counts = t(rpca_biri), 
                logcounts = t(rpca_biri)), 
  colData = coord_biri,
  spatialCoordsNames = c("x", "y"))

gc()


# Subset the object
# change the elements in sample_list with the names of your samples (same elements in metadata$sample_id)
sample_names <- c('biri1','biri2','biri3')
spe_list <- lapply(sample_names, function(x) spe[, spe$kidney == x])
spe_list
table(spe$kidney)


# Preprocess the data
# no need to preprocess the data (already normalized and n features is just 300 so no need to find hvg)


# Running BANKSY

# To run BANKSY across multiple samples, we first compute the BANKSY neighborhood
# feature matrices for each sample separately. 
# k-geoms = 24 means two order of neighbors with stereo-seq
compute_agf <- TRUE
k_geom <- 8
aname <- 'logcounts'
spe_list <- lapply(spe_list, computeBanksy, assay_name = aname, compute_agf = compute_agf, k_geom = k_geom)
gc()

# We then merge the samples to perform joint dimensional reduction and clustering:

spe_joint <- do.call(cbind, spe_list)


# When running multi-sample BANKSY PCA, the `group` argument may be provided. 
# This specifies the grouping variable for the cells or spots across the samples.
# Features belonging to cells or spots corresponding to each level of the 
# grouping variable will be z-scaled separately. 
# Here the column that indicates this is sample_id in colData(spe)

lambda <- 1
use_agf <- TRUE
spe_joint <- runBanksyPCA(spe_joint, use_agf = use_agf, lambda = lambda, group = "kidney", scale = FALSE, seed = 1000)
gc()


# Run UMAP on the BANKSY embedding:
spe_joint <- runBanksyUMAP(spe_joint, use_agf = use_agf, lambda = lambda, seed = 1000)
invisible(gc())

# Finally, obtain cluster labels for spots across all 3 samples.
res <- 1
spe_joint <- clusterBanksy(spe_joint, use_agf = use_agf, lambda = lambda, resolution = res, seed = 1000)


clusters <- colData(spe_joint)


wd <- "/exports/humgen/bmanzato/banksy_data/"

write.csv(clusters,paste(wd,'output_data/bIRI_banksy_QMSI_RPCA_lambda1_kgeom8_res1_agftrue_scalefalse_130lipids_NOgaps.csv',sep=''))





