library(dplyr)
library(ggplot2)
wd_data <- '/exports/nieromics-hpc/spatial_multiomics_data_MS/'

steseq_meta <- read.csv(paste0(wd_data, 'transcriptomics/metadata_complete.csv')) %>%
  filter(condition == "IRI") %>%
  rename("ss_x" = x) %>%
  rename("ss_y" = y)

qmsi_correct_banksy <- read.csv(paste0(wd_data, 'metabolomics/banksy/bIRI_banksy_QMSI_RPCA_lambda1_kgeom8_res1_agftrue_scalefalse_130lipids_NOgaps.csv')) %>%
  rename("banksy_correct" = clust_M1_lam1_k50_res1) %>%
  select(X, banksy_correct)

qmsi_transformed <- read.csv(paste0(wd_data, 'overlay/qmsi_transformed_biri1.csv')) %>%
  rename("qm_x" = x) %>%
  rename("qm_y" = y) %>%
  left_join(qmsi_correct_banksy, by = "X")

biri_ss_sub <- steseq_meta %>% filter(sample == "IRI1")
qmsi_rescaled <- qmsi_transformed

d <- 50  # threshold distance

# Query the kd-tree to get the nearest neighbor distances for each point in qmsi_rescaled
nn_result <- nn2(
  data = biri_ss_sub[, c("ss_x", "ss_y")],
  query = qmsi_rescaled[, c("qm_x", "qm_y")],
  k = 1
)

# Extract distances and closest point indices
distances <- nn_result$nn.dists[, 1]
closest_indices <- nn_result$nn.idx[, 1]

# Filter msi to keep only points within the distance threshold
qmsi_filtered <- qmsi_rescaled[distances <= d, ]
valid_indices <- closest_indices[distances <= d]  # Keep only indices within the threshold
qmsi_filtered$ss_index <- biri_ss_sub$X[valid_indices]

# Filter for SDs of interest using correct banksy annotation
q234 <- qmsi_filtered %>% filter(banksy_correct %in% c("2", "3", "4")) %>%
  mutate(q_overlay = "SD234")
q7 <- qmsi_filtered %>% filter(banksy_correct %in% c("7")) %>%
  mutate(q_overlay = "SD7")
q9 <- qmsi_filtered %>% filter(banksy_correct %in% c("9")) %>%
  mutate(q_overlay = "SD9")
qoverlay <- rbind(q234, q7, q9)

# Check and filter for pixels that are assigned to both
inconsistent_rows <- qoverlay %>%
  group_by(ss_index) %>%
  summarise(unique_q_overlay = n_distinct(q_overlay)) %>%
  filter(unique_q_overlay > 1) %>%
  pull(ss_index)  # Extract the list of problematic indices
qoverlay_clean <- qoverlay %>%
  filter(!ss_index %in% inconsistent_rows)

#write.csv(qoverlay_clean, paste0(wd_data, 'overlay/overlay_biri2.csv'))

# Prepare the combined df with stereoseq
qoverlay <- read.csv(paste0(wd_data, 'overlay/overlay_biri1.csv'))
steseq <- steseq_meta %>% filter(sample == "IRI1") %>% #Put iri3 here 
  rename('ss_index' = "X")
overlay <- steseq %>%
  left_join(qoverlay, by = 'ss_index') %>%
  select(ss_index, ss_x, ss_y, sample, figure_idents, X, q_overlay) %>%
  mutate(q_overlay = ifelse(is.na(q_overlay), "other", q_overlay))

proj_palette <- c("#EAEAEA", "#739BD7", "#E8BE49","#9F1E22")
ggplot(overlay, aes(x = ss_x, y = ss_y, color = q_overlay)) +
  geom_point(size = 0.45, shape = 15) +
  scale_color_manual(values = proj_palette) +
  labs(color = "SD") +
  theme_void() +
  theme(legend.key.size = unit(1, 'cm')) +  
  guides(color = guide_legend(override.aes = list(size = 4))) +
  coord_fixed()  # Fix the aspect ratio

towrite <- na.omit(overlay) %>%
  filter(figure_idents %in% c("PT-S1/S2", "PT-S3")) %>%
  mutate(projection = paste(figure_idents, q_overlay, sep = "_")) %>%
  group_by(ss_index) %>%
  summarise(projection = first(projection)) %>%
  ungroup() %>%
  select(ss_index, projection)

write.csv(towrite, paste0(wd_data, 'overlay/overlay_pts_biri1.csv') )