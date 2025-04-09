library(dplyr)
library(readxl)
library(ggplot2)


data <- iri %>% 
  filter(kidney == 'biri1') %>%
  mutate(niche = case_when(
    ban_idents %in% c(2, 3, 4) ~ "Healthy niche",
    ban_idents == 7 ~ "Injured niche 1",     
    ban_idents == 9 ~ "Injured niche 2", 
    TRUE ~ "Other")
  )

#Visualization
create_spatial_visualization <- function(data, metab, nic) {
  # Normalize the data for the selected metabolite
  subdata <- data
  subdata$normalized_percentage <- (subdata[[metab]] / max(subdata[[metab]], na.rm = TRUE)) * 100 * 1.65
  
  # Hotspot removal
  quantile_99 <- quantile(subdata$normalized_percentage, 0.99, na.rm = TRUE)
  subdata$normalized_percentage <- ifelse(subdata$normalized_percentage > quantile_99, quantile_99, subdata$normalized_percentage)
  
  # Compute convex hull for the outline
  hull <- subdata %>%
    filter(normalized_percentage > 0) %>%  # Only include non-black data points
    slice(chull(x, y))  # Compute convex hull
  
  # Create the ggplot visualization
  q <- ggplot(subdata, aes(x = x, y = y, fill = normalized_percentage)) +
    geom_tile() +  # Main tiles
    geom_tile(data = subset(subdata, niche != nic),
              fill = "black", alpha = 0.80) +  # Black overlay
    geom_polygon(data = hull, aes(x = x, y = y), fill = NA, color = "white", linewidth = 0.3) +  # Convex hull outline
    scale_fill_viridis_c(option = "viridis", name = "Intensity (%)", limits = c(0, quantile_99)) +  
    labs(fill = "Intensity (%)") +
    ggtitle(paste("Spatial Visualization of", metab, "in", nic)) +
    theme_minimal() +
    theme(
      panel.grid = element_blank(), 
      panel.background = element_rect(fill = "black", color = NA), 
      plot.background = element_rect(fill = "black", color = NA), 
      text = element_text(color = "white"), 
      axis.text = element_text(color = "white"), 
      axis.title = element_text(color = "white"), 
      plot.title = element_text(color = "white")
    ) +
    coord_fixed()
  
  return(q)
}

result_plot <- create_spatial_visualization(data, "FFA(18:2)", "Healthy niche")
print(result_plot)