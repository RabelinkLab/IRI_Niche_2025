# Based on https://katrionagoldmann.github.io/volcano3D/index.html
library(volcano3D)
library(dplyr)
library(tidyr)
library(ggplot2)


# Load metabolite data
imputed_iri <- read.csv("//vf-nieromics.researchlumc.nl/NIERomics$/qMSI_data/imputed_qMSI_data/qMSI_countmatrix_imputed.csv")
iri_annotated <- read.csv( "W:/data/rrietjens-Rosalie/StereoSeq_qMSI/annotated_iri.csv")
iri <- merge(imputed_iri, iri_annotated, by = "spot")

data <- iri %>% 
  mutate(niche = case_when(
    ban_idents %in% c(2, 3, 4) ~ "Healthy niche",
    ban_idents == 7 ~ "Injured 1 niche", 
    ban_idents == 9 ~ "Injured 2 niche"
  )) %>%
  na.omit()

# Reformat data according to volcano3D requirements
data$grouping <- paste(data$niche, data$fig_idents)

data <- data %>%
  select(spot, kidney, ban_idents, fig_idents, niche, grouping, everything(), -sample_id, -x, -y, -X) %>%
  filter(grouping %in% c("Healthy niche PT-S1/S2", "Injured 1 niche PT-S1/S2", "Injured 2 niche PT-S1/S2"))

rld <- data %>%
  select(spot, everything(), -kidney, -ban_idents, -fig_idents, -niche, -grouping,) 
rownames(rld) <- rld$spot
rld[,1] <- NULL

# Cleaning up of lipids not taken along for DE analysis
lipids_toremove <- c(
  "CPA.18.0", "CPA.18.1", "LPA.16.0", "LPA.18.0", "LPA.18.1", 
  "LPA.18.2", "LPA.O.18.1", "PA.34.1", "PA.34.2", "PA.36.3", 
  "PA.36.4", "PA.38.6", "PE.34.1", "PE.34.2", "PE.36.1", 
  "PE.36.2", "PE.36.3", "PE.36.4.PE.O.36.5.O", "PE.38.3", 
  "PE.38.4.PE.O.38.5.O", "PE.38.5.PE.O.38.6.O", "PI.34.1", "PI.36.2"
)
rld <- rld %>%
  select(-all_of(lipids_toremove))
rld <- t(rld) # To match the input of the example

metadata <- data[,1:6]
rownames(metadata) <- metadata$spot
metadata$grouping <- factor(metadata$grouping, 
                            levels = c("Healthy niche PT-S1/S2", "Injured 1 niche PT-S1/S2", "Injured 2 niche PT-S1/S2"))

# Change the default (anova & t-test) since the data is not normally distributed
# Tested normal distribution visually with qqnorm/line and ks.test
polar <- polar_coords(outcome = metadata$grouping,
                      data = t(rld), 
                      group_test = "kruskal.test", 
                      pairwise_test = "wilcoxon"
)

# Subset and view individually
pvals_subset <- significance_subset(polar, 
                                    significance = c("I2nP+"), 
                                    output="pvals")
dat <- head(pvals_subset) %>%
  kable()

# Outcome in numbers:
# S12: H = 16, I12 = 3, I2 = 0, HI2 = 9, HI1 = 10
# S3: H = 5, I12 = 3, I2 = 2, HI2 = 15, HI1 = 14

# Visualization 3D interactive:
p <- volcano3D(polar,
               label_rows = c("Succinic.acid",  "FFA.18.2."),
               label_size = 10,
               height = 500)
p

add_animation(p)

# Static ggplot
# For PT-S1/S2 labels: c(4,5,8,10,11,18,19,33,35)
# For PT-S3 labels: c(3,4,6,12,13,19,22,26,31,33)

radial_ggplot(polar = polar,
              label_rows = c(3,4,6,12,13,19,22,26,31,33),
              marker_size = 6,
              marker_outline_width = 0,
              legend_size = 10) +
  theme(legend.position = "right")


