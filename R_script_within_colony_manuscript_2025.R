---
title: "Exploring Within-Colony Thermal Tolerance Variability in Reef Building Corals"
author: "Marlen Schlotheuber"
format: html
editor: visual
---

# R code for Manuscript

### Install needed packages: CBASS package

```{r}
if(!require(devtools)){
   install.packages("devtools")
}

devtools::install_github("reefgenomics/CBASSED50@0.1.5.1", force=TRUE)
```

### Other packages

```{r}
install.packages("ggmap")
install.packages("osmdata")
install.packages("ggspatial")
library(sf)
library(rnaturalearth)
library(ggmap)
library(dplyr)
library(tidyr)
library(ggplot2)
library(readxl)
library(rstudioapi)
library(RColorBrewer)
library(CBASSED50)
library(tidyverse)
library(patchwork)
library (scales)
library(car)
library(mgcv)
library(stringr)
```

## **Fig.1 Study regions and the sampling locations within massive Porites colonies**

```{r}
# Study 
# Get country boundaries
world <- ne_countries(scale = "medium", returnclass = "sf")

# Subset for Indonesia and Malaysia
countries <- world[world$name %in% c("Indonesia", "Malaysia"), ]

# plot
map_plot <-  ggplot() +
  geom_sf(data = world, fill = "white", color = "black") +
  geom_sf(data = subset(countries, name == "Indonesia"), fill = "white", color = "black") +
  geom_sf(data = subset(countries, name == "Malaysia"), fill = "white", color = "black") +
  ggtitle("Sampling sites: Indonesia and Malaysia") +
  coord_sf(xlim = c(90, 155), ylim = c(-35, 35), expand = FALSE) +
  scale_x_continuous(name = "Longitude (°E)", breaks = seq(90, 155, 10)) +
  scale_y_continuous(name = "Latitude (°N)", breaks = seq(-35, 35, 10)) +
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "grey90", color = NA),  # ocean
    plot.background  = element_rect(fill = "white", color = NA),   # outer plot bg
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

map_plot

# Save map as PDF
getwd()

ggsave("/Users/MarlenSchlotheuber/Nextcloud/10. projects/11. PhD/10. projects/Intraspecific bleaching project/Project planning/MYS_2024/Fieldwork/Data_26042024/Manuscript/Intra-colonial/Final_Figures_01092025/Fig.1_Indonesia_Malaysia_Map.png",
       plot = map_plot, width = 7, height = 10, units = "in")
```

# Thermal tolerance data from Indonesia and Malaysia

To examine the photosynthetic efficiency data from Malaysia and Indonesia, the R-package from Iakovleva, Colin & Voolstra, 2025 was used ([https://doi.org/10.5281/zenodo.8370644](https://doi.org/10.5281/zenodo.8370644.)).

The photosynthetic efficiency data from the four temperature profiles of each sampling location from Indonesia and Malysia was used to fit log-logistic dose-response curves, using the ‘drc’ R package. These curves were used to determine the effective dose 50 (ED50) for each sampling location (Evensen et al. 2021; Voolstra et al. 2021). The ED50 value defines the temperature at which photosynthetic efficiency has declined by 50% relative to the baseline measurement, and is therefore used as an empirically derived proxy for heat stress tolerance (Evensen et al. 2022). 

### Load Data ED50 data

```{r}
# Indonesia
# Upload Indonesia data
input_data_path <- selectFile(
  caption = "Select XLSX or CSV Input File")

# ID
cbass_dataset_ID <- read_data(input_data_path)

# To specify the prefix for output files
output_prefix <- tools::file_path_sans_ext(input_data_path)

rlog::log_info(paste("Your current directory is", getwd()))
rlog::log_info(paste("Your input filename is", basename(input_data_path)))
rlog::log_info(paste("The output files will be written into", output_prefix))

# Malaysia data
# upload Malaysia data
input_data_path <- selectFile(
  caption = "Select XLSX or CSV Input File")

# MYS
cbass_dataset_MYS <- read_data(input_data_path)

# To specify the prefix for output files
output_prefix <- tools::file_path_sans_ext(input_data_path)

rlog::log_info(paste("Your current directory is", getwd()))
rlog::log_info(paste("Your input filename is", basename(input_data_path)))
rlog::log_info(paste("The output files will be written into", output_prefix))

```

### Analysis and Preprocessing

```{r}
# ID
cbass_dataset_ID <- preprocess_dataset(cbass_dataset_ID)
validate_cbass_dataset(cbass_dataset_ID)

# MYS
cbass_dataset_MYS <- preprocess_dataset(cbass_dataset_MYS)
validate_cbass_dataset(cbass_dataset_MYS)
```

### Explore ED50

```{r}
# ID
# grouping properties Indonesia
grouping_properties <- c("Country", "Condition")

# PAM data Indonesia
drm_formula <- "PAM ~ Temperature"
models <- fit_drms(cbass_dataset_ID, grouping_properties, drm_formula, is_curveid = TRUE)
str(cbass_dataset_ID)

# PAM data 
ed50_ID <- get_ed50_by_grouping_property(models)
cbass_dataset_ID <- define_grouping_property(cbass_dataset_ID, grouping_properties) %>%
  mutate(GroupingProperty = paste(GroupingProperty, Genotype, sep = "_"))

ed50_df_ID <- 
  left_join(ed50_ID, cbass_dataset_ID, by = "GroupingProperty") %>%
  dplyr::select(names(ed50_ID), Depth, Genotype, all_of(grouping_properties)) %>%
  distinct()

# MYS
# grouping properties
grouping_properties <- c("Country", "Condition")
# important: keep gouping property: Country, Condition to get the right EDS50 plot !!!!

# PAM data Malaysia
drm_formula <- "PAM ~ Temperature"
models <- fit_drms(cbass_dataset_MYS, grouping_properties, drm_formula, is_curveid = TRUE)

# PAM data
ed50_MYS <- get_ed50_by_grouping_property(models)
cbass_dataset_MYS <- define_grouping_property(cbass_dataset_MYS, grouping_properties) %>%
  mutate(GroupingProperty = paste(GroupingProperty, Genotype, sep = "_"))

ed50_df_MYS <- 
  left_join(ed50_MYS, cbass_dataset_MYS, by = "GroupingProperty") %>%
  dplyr::select(names(ed50_MYS), Depth, Genotype, all_of(grouping_properties)) %>%
  distinct()

```

### Adjust dataset for plotting

```{r}
# ID
# Create the Location column and name it as the sampling locazions based on the Condition column
cbass_dataset_ID <- cbass_dataset_ID %>%
  mutate(Location = case_when(
    Condition %in% c("Top") ~ "top",
    Condition %in% c("Bottom") ~ "bottom"))

ed50_df_ID <- ed50_df_ID %>%
  mutate(Location = case_when(
    Condition %in% c("Top") ~ "top",
    Condition %in% c("Bottom") ~ "bottom"))

# ID
str(cbass_dataset_ID)
cbass_dataset_ID$Genotype <- as.numeric(cbass_dataset_ID$Genotype)
cbass_dataset_ID$Location <- as.factor(cbass_dataset_ID$Location)
cbass_dataset_ID$Temperature <- as.factor(cbass_dataset_ID$Temperature)
cbass_dataset_ID$PAM <- as.factor(cbass_dataset_ID$PAM)



# MYS
# add top, middle, bottom to the cbass MYS dataset as sampling locations based on the Condition column 
cbass_dataset_MYS <- cbass_dataset_MYS %>%
  mutate(Location = case_when(
    Condition %in% c("E1", "W1") ~ "top",    
    Condition %in% c("E2", "W2") ~"upper_middle",
    Condition %in% c("E3", "W3") ~ "lower_middle",
    Condition %in% c("E4", "W4") ~ "bottom", 
    TRUE ~ NA_character_  
  ))

# also add to the ed50_df dataset
ed50_df_MYS <- ed50_df_MYS %>%
  mutate(Location = case_when(
    Condition %in% c("E1", "W1") ~ "top",    
    Condition %in% c("E2", "W2") ~"upper_middle",
    Condition %in% c("E3", "W3") ~ "lower_middle",
    Condition %in% c("E4", "W4") ~ "bottom", 
    TRUE ~ NA_character_  
  ))

# MYS
# change genotype numbers to numeric so its the same as for the MYS dataset
str(cbass_dataset_MYS)
cbass_dataset_MYS$Genotype <- as.numeric(cbass_dataset_MYS$Genotype)
cbass_dataset_MYS$Location <- as.factor(cbass_dataset_MYS$Location)
```

### Temperature response curve

```{r}
# ID
# set order of the Locations
cbass_dataset_ID$Location <- factor(cbass_dataset_ID$Location, 
                                   levels = c("bottom", "lower_middle", "upper_middle", "top"))

# visualize temperature response curve 
exploratory_curve <-
  ggplot(data = cbass_dataset_ID,
       aes(
         x = Temperature,
         y = PAM,
         # You can play around with the group value (e.g., Species, Site, Condition)
         group = GroupingProperty,
         color = Location)) +
  geom_smooth(
    method = drc::drm,
    method.args = list(
      fct = drc::LL.3()),
    se = FALSE,
    size = 0.7
  ) +
  geom_point(size = 1.5) +
  facet_wrap( ~ Genotype, nrow = 2) +
  scale_color_brewer(palette = "Paired") # Colorblind-friendly palette

#plot
exploratory_curve


# MYS
# set order of the Locations
cbass_dataset_MYS$Location <- factor(cbass_dataset_MYS$Location, 
                                   levels = c("bottom", "lower_middle", "upper_middle", "top"))

# visualize temperature response curve
exploratory_curve <-
  ggplot(data = cbass_dataset_MYS,
       aes(
         x = Temperature,
         y = PAM,
         # You can play around with the group value (e.g., Species, Site, Condition)
         group = GroupingProperty,
         color = Location)) +
  geom_smooth(
    method = drc::drm,
    method.args = list(
      fct = drc::LL.3()),
    se = FALSE,
    size = 0.7
  ) +
  geom_point(size = 1.5) +
  facet_wrap( ~ Genotype, nrow = 2) +
  scale_color_brewer(palette = "Paired") # Colorblind-friendly palette

# plot
exploratory_curve
```

## **Fig 2. Thermal tolerance across the surface of massive Porites sp. colonies from Indonesia and Malaysia.**

```{r}
# Indonesia

# only select needed columns from the indonesia dataset
ed50_df_ID_filtered <- ed50_df_ID %>%
  dplyr::select(c("ED50", "Depth", "Location", "Genotype", "Country"))

# calculate mean ED50 per colony/genotype and calculate differences
ED50_ID_means <- ed50_df_ID_filtered %>%
  group_by(Genotype) %>%
  summarise(mean_ED50 = mean(ED50))

Final_ED50_ID <- ed50_df_ID_filtered %>%
  left_join(ED50_ID_means, by = c("Genotype")) %>%
  mutate(ED50_diff = ED50 - mean_ED50) %>%
  select(-contains(".x")) %>% # remove duplicates
  distinct()  # Re

# change genotype numbre from 01 to 1
Final_ED50_ID$Genotype <- as.numeric(Final_ED50_ID$Genotype)

# prepare for plotting
# set the order of the locations
Final_ED50_ID$Location <- factor(Final_ED50_ID$Location, 
                                   levels = c("bottom", "top"))

# set the order of the coral colonies
Final_ED50_ID <- Final_ED50_ID %>%
  mutate(
    Genotype = factor(Genotype, levels = c("1", "2", "3", "4", "5", "6", "7"))  # Order genotypes
  )

# we don't have ITS2 and 16S data for colony 5 so remove colony 5 also from the ED50 data
# Remove Genotype 5 from the dataset
Final_ED50_ID_filtered <- Final_ED50_ID[Final_ED50_ID$Genotype != "5", ]

 # set the colors               
genotype_colors <- c("1" = "pink4", 
                     "2" = "orchid4", 
                     "3" = "pink2", 
                     "4" = "mistyrose2", 
                     "6" = "lightcoral",
                     "7" = "hotpink4")

# ED50 Difference Bar Plot
ED50_diff_barplot_ID <- ggplot(Final_ED50_ID_filtered, aes(x = ED50_diff, y = Location, fill = Genotype)) +
  geom_vline(aes(xintercept = 0), color = "black", linetype = "dashed", size = 1) +
  geom_vline(xintercept = -0.5, color = "grey40", linetype = "dashed") +
  geom_vline(xintercept = 0.5, color = "grey40", linetype = "dashed") +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~Genotype, nrow = 6) +
  xlim(-1, 1) +
  scale_y_discrete(labels = c("Bottom", "Top")) +
  scale_fill_manual(values = genotype_colors) +
  labs(
    x = "Deviation from Mean ED50 [°C]", 
    y = "Location within Colony",
    fill = "Genotype"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    text = element_text(size = 13),
    strip.text = element_blank(),
    axis.title = element_text(size = 13, face = "bold"),
    axis.text = element_text(size = 13),
    legend.position = "none",
    panel.spacing = unit(0.3, "lines"),
    strip.background = element_blank()
  )


# ED50 Violin Plot
ED50_violin_plot_ID <- ggplot(Final_ED50_ID_filtered, aes(x = 1, y = ED50, fill = Genotype)) +
  geom_violin() +
  facet_wrap(~Genotype, nrow = 7) +
  stat_summary(fun = "mean", geom = "point", color = "black", size = 3, shape = 18, show.legend = FALSE) +
  stat_summary(fun = "mean", geom = "text", color = "black", size = 4,
               aes(label = round(..y.., 2)), vjust = -1) +
  scale_y_continuous(breaks = c(37, 38), limits = c(36.5, 38.5)) +
  scale_fill_manual(values = genotype_colors) +
  labs(
    x = NULL,
    y = "ED50 [°C]",
    fill = "Genotype"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    text = element_text(size = 13),
    strip.text = element_blank(),
    axis.title = element_text(size = 13, face = "bold"),
    axis.text = element_text(size = 13),
    legend.position = "none",
    panel.spacing = unit(0.3, "lines"),
    strip.background = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )


# Combine plots
combined_plot_ID <- ED50_violin_plot_ID + ED50_diff_barplot_ID + 
  plot_layout(ncol = 2, widths = c(0.4, 0.8))

# Display combined plot
combined_plot_ID

# save as image
ggsave("/Users/MarlenSchlotheuber/Nextcloud/10. projects/11. PhD/10. projects/Intraspecific bleaching project/Project planning/MYS_2024/Fieldwork/Data_26042024/Manuscript/Final_Figures/Fig3_ED50_ID_within_colonies.png", plot = combined_plot_ID, width = 12, height = 10, dpi = 300, bg = "white")

# save as pdf
ggsave("/Users/MarlenSchlotheuber/Nextcloud/10. projects/11. PhD/10. projects/Intraspecific bleaching project/Project planning/MYS_2024/Fieldwork/Data_26042024/Manuscript/Final_Figures/Fig3_ED50_ID_within_colonies.pdf", 
       plot = combined_plot_ID, 
       width = 12, 
       height = 10, 
       dpi = 300, 
       bg = "white")



# Malayisa
# merge datasets from MY to combine ed50 data with general sampling data
Final_ED50_MYS  <- merge(cbass_dataset_MYS, ed50_df_MYS)

# create the column Orientation
Final_ED50_MYS <- Final_ED50_MYS %>%
  mutate(Orientation = case_when(
    Site %in% c("E") ~ "East",
    Site %in% c("W") ~ "West"))

# set order of the Locations
Final_ED50_MYS$Location <- factor(Final_ED50_MYS$Location, 
                                   levels = c("bottom", "lower_middle", "upper_middle", "top"))
 
# set order of the colonies with increasing depth
Final_ED50_MYS <- Final_ED50_MYS %>%
  mutate(
    Genotype = factor(Genotype, levels = c("1", "6", "2", "3", "5", "4")))

# adjust Orientation column
Final_ED50_MYS <- Final_ED50_MYS %>%
  mutate(
    Orientation = factor(Orientation, levels = c("West", "East")))

# Compute mean ED50 per Genotype
ED50_MYS_means <- Final_ED50_MYS %>%
  group_by(Genotype) %>%
  summarise(mean_ED50 = mean(ED50))

# Merge back to main dataset (on Genotype *and* Orientation)
ED50_MYS_means$Genotype <- as.factor(ED50_MYS_means$Genotype)
str(ED50_MYS_means)

# merge datasets
Final_ED50_MYS <- Final_ED50_MYS %>%
  left_join(ED50_MYS_means, by = c("Genotype")) 


Final_ED50_MYS <- Final_ED50_MYS %>%
  mutate(ED50_diff = ED50 - mean_ED50) %>%  # Compute difference
  select(-contains(".x")) %>%  # Remove duplicate columns (like Orientation.x)
  distinct()  # Re
           
# set the colors    
 genotype_colors <- c("1" = "skyblue3", 
                     "2" = "lightblue", 
                     "3" = "dodgerblue4", 
                     "4" = "cornflowerblue", 
                     "5" = "lightskyblue", 
                     "6" = "dodgerblue3")


# Bar plot
ED50_diff_barplot_MYS <- ggplot(Final_ED50_MYS, aes(x = ED50_diff, y = Location, fill = Genotype)) +
  geom_vline(aes(xintercept = 0), color = "black", linetype = "dashed", size = 1) +
  geom_vline(xintercept = -0.5, color = "grey40", linetype = "dashed") +
  geom_vline(xintercept = 0.5, color = "grey40", linetype = "dashed") +
  geom_bar(stat = "identity", position = "dodge") +
  facet_grid(Genotype ~ Orientation, scales = "free_y") +
  xlim(-1, 1) +
  scale_y_discrete(labels = c("Bottom", "Lower Middle", "Upper Middle", "Top")) +
  scale_fill_manual(values = genotype_colors) +
  labs(
    x = "Deviation from Mean ED50 [°C]", 
    y = "Location within Colony",
    fill = "Genotype"
  ) +
  theme_minimal(base_size = 13) +  # Global text size
  theme(
    text = element_text(size = 13),
    strip.text = element_blank(),
    axis.title = element_text(size = 13, face = "bold"),  
    axis.text = element_text(size = 13),  
    legend.position = "none",
    panel.spacing = unit(0.4, "lines"),
    strip.background = element_blank()
  )


# Violin plot
ED50_violin_plot_MYS <- ggplot(Final_ED50_MYS, aes(x = 1, y = ED50, fill = Genotype)) +
  geom_violin() +
  facet_wrap(~Genotype, nrow = 6) +
  stat_summary(fun = "mean", geom = "point", color = "black", size = 3, shape = 18, show.legend = FALSE) + 
  stat_summary(fun = "mean", geom = "text", color = "black", size = 4,  # Keep this a little larger for visibility
               aes(label = round(..y.., 2)), vjust = -1) +
  scale_fill_manual(values = genotype_colors) +
  guides(fill = guide_legend(title = "Coral Colony")) +
  labs(
    x = NULL,
    y = "ED50 [°C]",
    fill = "Genotype"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    text = element_text(size = 13),
    strip.text = element_blank(),
    axis.title = element_text(size = 13, face = "bold"),
    axis.text = element_text(size = 13),
    legend.position = "none",
    panel.spacing = unit(0.4, "lines"),
    strip.background = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )


# Combine
combined_plot_MYS <- ED50_violin_plot_MYS + ED50_diff_barplot_MYS +
  plot_layout(ncol = 2, widths = c(0.4, 0.8))

combined_plot_MYS

# save as image
ggsave("/Users/MarlenSchlotheuber/Nextcloud/10. projects/11. PhD/10. projects/Intraspecific bleaching project/Project planning/MYS_2024/Fieldwork/Data_26042024/Manuscript/Final_Figures/Fig2_ED50_MYS_within_colonies.png", plot = combined_plot_MYS, width = 12, height = 10, dpi = 300, bg = "white")

# save as pdf
ggsave("/Users/MarlenSchlotheuber/Nextcloud/10. projects/11. PhD/10. projects/Intraspecific bleaching project/Project planning/MYS_2024/Fieldwork/Data_26042024/Manuscript/Final_Figures/Fig2_ED50_MYS_within_colonies.pdf", 
       plot = combined_plot_MYS, 
       width = 12, 
       height = 10, 
       dpi = 300, 
       bg = "white")


# plot together
# Define each row (violin left, barplot right)
row1 <- ED50_violin_plot_ID     + ED50_diff_barplot_ID     + plot_layout(ncol = 2, widths = c(0.4, 0.8))
row2 <- ED50_violin_plot_MYS + ED50_diff_barplot_MYS + plot_layout(ncol = 2, widths = c(0.4, 0.8))

# Combine plots into one row
final_combined_plot <-  ED50_diff_barplot_ID + ED50_diff_barplot_MYS +
  plot_layout(ncol = 4, widths = c(0.4, 0.8, 0.4, 0.8)) +
  plot_annotation(
    theme = theme(
      plot.margin = margin(10, 10, 10, 10)
    )
  )

# Display
final_combined_plot

# figure was safed as PDF and edited in Affinity designer
# Added the mean ED50 from the violin plot manually to the image

# Save as PNG: A4 landscape width, 1/3 height
ggsave(
  filename = "ITS2_MYS_ID.png",
  plot = final_combined_plot,
  width = 10,      # A4 width (landscape)
  height = 8,     # 1/3 of A4 height
  dpi = 300,
  units = "in",
  bg = "white"
)

# Save as PDF with same dimensions
ggsave(
  filename = "ITS2_MYS_ID.pdf",
  plot = final_combined_plot,
  width = 11.7,
  height = 8,
  units = "in"
)
getwd()
```

### Testing for significant differences in thermal tolerance

```{r}
# Indonesia
# check for significant differences
# data distribution
hist(Final_ED50_ID$ED50) # looks like its not normally distributed, but also not many data points
shapiro.test(Final_ED50_ID$ED50) # p-value: 0.8591, not signfiicant, so data is normaly distributed
qqnorm(Final_ED50_ID$ED50) # normal distributed

# data distribution
hist(Final_ED50_ID$ED50_diff) # looks like its not normally distributed
shapiro.test(Final_ED50_ID$ED50_diff) # p-value: 0.5761, normaly distributed
qqnorm(Final_ED50_ID$ED50_diff) # normal distributed

# across colonies
# can't test for significant differences across genotypes because for each genotpye I only have one datapoint
ggplot(Final_ED50_ID, aes(x = Genotype, y = ED50, color = Location )) +
  geom_point(size = 3)+
  facet_wrap(~Location)

# within colonies ED50
# check for differences in variance between the top and bottom sampling location
leveneTest(ED50 ~ Location, data = Final_ED50_ID)
# p-value = 0.9143, not significant difference in the variance of groups, can use ANOVA

# two groups only top and bottom
t.test(ED50 ~ Location, data = Final_ED50_ID)
# p-value = 0.797 not significant 


# Malaysia

# filter data
Final_ED50_MYS <- Final_ED50_MYS %>% 
  filter(Temperature == 30)

# check data distribution: normality and homogeneity
hist(Final_ED50_MYS$ED50) # could be normally distirbuted
# normallity test
shapiro.test(Final_ED50_MYS$ED50) # p-value: 0.1506, so we do not reject the 0-hypthesis --> data normally distributed 
# plot normallity
qqnorm(Final_ED50_MYS$ED50) # looks like its normally distributed

# Ed50 across depth
# check for corrleation of ED50 and Depth, data is normally distributed but check what relationship deoth and ED50 have (linear, non-linear)

# what model fits best: linear or gam?
lm_model <- lm(ED50 ~ Depth, data = Final_ED50_MYS)
summary(lm_model) # p-value = <2e-16 ***
gam_model <- gam(ED50 ~ s(Depth), data = Final_ED50_MYS)
summary(gam_model) # <2e-16 ***

# check for the best fitting model
AIC(lm_model, gam_model)

# gam fits better has lower AIC, so data is rather not linear -> use spearman to check for correlation
cor.test(Final_ED50_MYS$ED50, Final_ED50_MYS$Depth, method = "spearman")
# p-value < 7.013e-06
# rho: -0.598 --> -0.60

# within colonies
# check that variance is normally distributed before using a normaly distributed test (e.g. ANOVA)
leveneTest(ED50 ~ Location, data = Final_ED50_MYS)
# p-value = 0.2624 --> not significan means the variance in groups is equal so variance in groups is equal we can use ANOVA

# test for significance using ANOVA (4 locations within colonies)
anova_result <- aov(ED50 ~ Location, data = Final_ED50_MYS)
summary(anova_result)
# p-value =0.677

# check for East
Final_ED50_East <- subset(Final_ED50_MYS, Orientation == "East")

# Check for West
Final_ED50_West <- subset(Final_ED50_MYS, Orientation == "West")

# ANOVA for East
# test homogeneity
leveneTest(ED50 ~ Location, data = Final_ED50_East)
# p-value: 0.725, not significant use ANOVA
anova_east <- aov(ED50 ~ Location, data = Final_ED50_East)
summary(anova_east)
# p-value= 0.84, not signf

# ANOVA for West
# test homogeneity
leveneTest(ED50 ~ Location, data = Final_ED50_West)
# p-value = 0.1565, not signficant use ANOVA
anova_west <- aov(ED50 ~ Location, data = Final_ED50_West)
summary(anova_west)
# p-value: 0.884


###########
# check for signficant differences in the ED50_Diff
# data distirbution
hist(Final_ED50_MYS$ED50_diff) #  normally distirbuted, or slightly left skewed
# normality
shapiro.test(Final_ED50_MYS$ED50_diff) # P-value above 0.8281, data normally distirbuted
# plot normallity
qqnorm(Final_ED50_MYS$ED50_diff, main = "Q-Q Plot of ED50") # looks like its normally distribute
# test homogeneity
leveneTest(ED50_diff ~ Location, data = Final_ED50_MYS)
# p-value = 0.418, not signficant use ANOVA

# normally distributed data
anova_result <- aov(ED50_diff ~ Location, data = Final_ED50_MYS)
summary(anova_result)
# p-value = 0.251 not significant

# ANOVA for East
anova_east <- aov(ED50_diff ~ Location, data = Final_ED50_East)
summary(anova_east)
# p-value = 0.595

# ANOVA for West
anova_west <- aov(ED50_diff ~ Location, data = Final_ED50_West)
summary(anova_west)
# p-value = 0.880

```

# ITS2 type profiles - Symbiodiniaceae community composition

ITS2 sequencing data was analysed using SymPortal (Hume et al., 2019). ITS2 type profiles being representative of putative Symbiodiniaceae genotypes were characterized by specific sets of defining intragenomic ITS2 sequence variants (DIVs) (Hume et al. 2019) and visualized using R studio.

### Malaysia: Load Data and Preprocessing

```{r}
# upload ITS2 type profile meta data
input_data_path <- selectFile(
  caption = "Select XLSX or CSV Input File")

its2_metadata <- read_data(input_data_path)

#upload ITS2 profile absolute abundance
its2_profile <- read.table("your_pathname_/560_20250219T125105_DBV_20250221T224139.profiles.relative.abund_and_meta.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

colnames(its2_profile)[1] <- "ITS2_profile"
colnames(its2_profile)[2] <- "sample_name"

# merge dataframes
its2_data<- merge(its2_metadata, its2_profile, by = "sample_name", all.x = TRUE)

# prepare for plotting
# Convert to long format for ggplot
its2_data_long <- its2_data %>%
  gather(key = "ITS2_profile", value = "value", C15.C15ev.C15dt.C15kl.C15pm:C15.C15by.C15ai)

# just select data with, top middle bottom and middle 1
# Colony 2 was so large that we were able to sample at 6 sampling locations within the colonies (top, upper middle 1, upper middle 2, lower middle 1, lower middle 2, bottom) but since no differences within the colony were detected we removed the two additional middle sampling locations upper and lower middle 2).
its2_data_long_filtered <- its2_data_long %>%
  filter(collection_depthrange %in% c("top", "upper_middle", "upper_middle_1", "lower_middle","lower_middle_1", "bottom"))

# rename upper and lower middle 1 to just lower and upper middle
its2_data_long_filtered <- its2_data_long_filtered %>%
  mutate(collection_depthrange = recode(collection_depthrange,
                                      "lower_middle_1" = "lower_middle",
                                      "upper_middle_1" = "upper_middle"))

# First, ensure your collection_location is a factor with appropriate labels.
its2_data_long_filtered$collection_location <- factor(its2_data_long_filtered$collection_location,
                                              levels = c(4, 3, 2, 1),
                                              labels = c("bottom", "lower_middle", "upper_middle", "top"))

# change the naming in collection_location
# rename upper and lower middle 1 to just lower and upper midlee for colony 2
its2_data_long_filtered <- its2_data_long_filtered %>%
  mutate(collection_location = recode(collection_location,
                                      "lower_middle" = "Lower Middle",
                                      "upper_middle" = "Upper Middle",
                                       "top" = "Top",
                                      "bottom" = "Bottom"))

# Ensure 'collection_depthrange' is ordered and 'collection_orientation' is a factor with "West" and "East"
its2_data_long_filtered <- its2_data_long_filtered %>%
  mutate(
    collection_location = factor(collection_location, levels = c("Bottom", "Lower Middle","Upper Middle", "Top")),
    collection_orientation = factor(collection_orientation, levels = c("West", "East"))
  )

```

### Plot ITS2 Type Profile

```{r}
# organize data
# set order of the colonies with decreasing depth
its2_data_long_filtered <- its2_data_long_filtered %>%
  mutate(
    collection_colony = factor(collection_colony, levels = c("1", "6", "2", "3", "5", "4")))

#plot
ITS2_plot_MYS <- ggplot(its2_data_long_filtered, aes(x = value, y = collection_location, fill = ITS2_profile)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_x_continuous(labels = scales::percent_format(), expand = c(0, 0)) +
  facet_grid(collection_colony ~ collection_orientation, scales = "free_y") +
  theme_minimal(base_size = 13) +
  guides(fill = guide_legend(ncol = 1)) +
  labs(
    x = "Relative Abundance",
    y = "Location within Colony",
    fill = "ITS2 Type Profiles"
  ) +
  scale_fill_manual(
    values = c(
      "C15.C15ev.C15dt.C15kl.C15pm" = "#66CC66",
         "C15.C15ad.C15ai.C116.C116aa" = "#6699CC",
      "C15.C15ev.C15by" = "#66CC99",
      "C15.C15ev.C15dt" = "#009999",
      "C15.C15ev" = "#99CCCC",
        "C15.C15by.C15ai" = "#006666"
    ),
    labels = c(
      "C15.C15ev.C15dt.C15kl.C15pm" = "C15-C15ev-C15dt-C15kl-C15pm",
      "C15.C15ad.C15ai.C116.C116aa" = "C15-C15ad-C15ai-C116-C116aa",
      "C15.C15ev.C15by" = "C15-C15ev-C15by",
      "C15.C15ev.C15dt" = "C15-C15ev-C15dt",
      "C15.C15ev" = "C15-C15ev",
      "C15.C15by.C15ai" = "C15-C15by-C15ai"
    )
  ) +
  theme(
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.text = element_text(size = 13),
    legend.title = element_text(size = 13),
    axis.text = element_text(size = 13),
    axis.title = element_text(size = 13),
    strip.text = element_text(size = 13),
    panel.grid = element_blank(),
    panel.border = element_blank()
  )

ITS2_plot_MYS

# save plot as image
ggsave("/Users/MarlenSchlotheuber/Nextcloud/10. projects/11. PhD/10. projects/Intraspecific bleaching project/Project planning/MYS_2024/Fieldwork/Data_26042024/Manuscript/Final_Figures/FIG4_ITS2_type_profiles_MY.png", plot = ITS2_plot_MYS, width = 8, height = 8, dpi = 300, bg = "white")

# save as pdf
ggsave("/Users/MarlenSchlotheuber/Nextcloud/10. projects/11. PhD/10. projects/Intraspecific bleaching project/Project planning/MYS_2024/Fieldwork/Data_26042024/Manuscript/Final_Figures/FIG4_ITS2_type_profiles_MY.pdf", 
       plot = ITS2_plot_MYS, 
       width = 8, 
       height = 8, 
       dpi = 300, 
       bg = "white")


```

### Indonesia: Load Data and Preprocessing

```{r}
# ITS2 type profile metadata
input_data_path <- selectFile(
  caption = "Select XLSX or CSV Input File")

# ITS2 profile data
its2_metadata <- read_data(input_data_path)


#upload ITS2 profile absolute abundance
its2_profile_ID <- read.table('your_pathname/582_20250509T135802_DBV_20250510T135051.profiles.relative.abund_and_meta.txt', header = TRUE, sep = "\t", stringsAsFactors = FALSE)

colnames(its2_profile)[1] <- "ITS2_profile"
colnames(its2_profile_ID)[2] <- "sample_name"

# merge dataframes
its2_data_ID <- merge(its2_metadata, its2_profile, by = "sample_name", all.x = TRUE)

# prepare for plotting
# Convert to long format for ggplot
its2_profile_ID_long <- its2_profile_ID %>%
  gather(key = "ITS2_profile", value = "value", C15.C15aae:C15.C15ad.C15ai.C116.C116aa)

# split colunm into seperate columns to extract the inforamtion
# Ensure sample_name is a character vector
its2_profile_ID_long$sample_name <- as.character(its2_profile_ID_long$sample_name)

# Split the sample_naem column and extract the components
split_info <- do.call(rbind, strsplit(its2_profile_ID_long$sample_name, "_"))

# Assign the components to new columns
its2_profile_ID_long$Coral_species <- split_info[, 1]
its2_profile_ID_long$Genotype <- split_info[, 2]
its2_profile_ID_long$Location <- split_info[, 3]

its2_data_ID_final <- its2_profile_ID_long

# Ensure order
its2_profile_ID_long <- its2_profile_ID_long %>%
  mutate(
    collection_location = factor(Location, levels = c("Bottom", "Top")))
```

### Plot ITS2 Type Profiles

```{r}

ITS2_plot_ID <- ggplot(its2_profile_ID_long, aes(x = value, y = Location, fill = ITS2_profile)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_x_continuous(labels = scales::percent_format(), expand = c(0, 0)) +
  facet_grid(Genotype ~ ., scales = "free_y") +
  theme_minimal(base_size = 13) +
    guides(fill = guide_legend(ncol = 1)) +
  labs(
    x = "Relative Abundance",
    y = "Location within Colony",
    fill = "ITS2 Type Profiles"
  ) +
  scale_fill_manual(
    values = c(
      "C15.C15aae" = '#006633',
      "C15.C15bn.C15by" = "#99CC33",
      "C15.C116aa.C15ad" = "#339900",
      "C57d.C57a" = "#000099",
      "C15.C15kl.C15he" = "#009966",
      "C15.C15ad.C15ai.C116.C116aa" = "#6699CC"
    ),
    labels = c(
      "C15.C15aae" = "C15-C15aae",
      "C15.C15bn.C15by" = "C15-C15bn-C15by",
      "C15.C116aa.C15ad" = "C15/C116aa-C15ad",
      "C57d.C57a" = "C57d/C57a",
      "C15.C15kl.C15he" = "C15-C15kl-C15he",
      "C15.C15ad.C15ai.C116.C116aa" = "C15-C15ad-C15ai-C116-C116aa"
    )
  ) +
  theme(
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.text = element_text(size = 13),
    legend.title = element_text(size = 13),
    axis.text = element_text(size = 13),
    axis.title = element_text(size = 13),
    strip.text = element_text(size = 13),
    panel.grid = element_blank(),
    panel.border = element_blank()
  )

ITS2_plot_ID

# save plot as image
ggsave("/Users/MarlenSchlotheuber/Nextcloud/10. projects/11. PhD/10. projects/Intraspecific bleaching project/Project planning/MYS_2024/Fieldwork/Data_26042024/Manuscript/Final_Figures/Fig5_ID_ITS2.png", plot = ITS2_plot_ID, width = 8, height = 8, dpi = 300, bg = "white")

# save as pdf
ggsave("/Users/MarlenSchlotheuber/Nextcloud/10. projects/11. PhD/10. projects/Intraspecific bleaching project/Project planning/MYS_2024/Fieldwork/Data_26042024/Manuscript/Final_Figures/Fig5_ID_ITS2.pdf", 
       plot = ITS2_plot_ID, 
       width = 12, 
       height = 10, 
       dpi = 300, 
       bg = "white")

```

## **Fig. 3 ITS2 type profiles across the surface of massive Porites sp. colonies from Indonesia and Malaysia.**

```{r}
# Combine plots side-by-side with equal width
final_combined_plot <- ITS2_plot_ID + ITS2_plot_MYS +
  plot_layout(ncol = 2, widths = c(1, 1)) +  # Equal widths for both plots
  plot_annotation(
    theme = theme(
      plot.margin = margin(10, 10, 10, 10)
    )
  )
final_combined_plot

# plot as saved and edited with Affinity designer.

# Save to A4 size (portrait orientation) as PDF
ggsave(
  filename = "ITS2_MYS_ID_new.pdf",
  plot = final_combined_plot,
  width = 11.7,      # A4 width in inches
  height = 10,     # A4 height in inches
  units = "in"
)
```

### Supplementary Fig. 1 ITS2 sequences - absolute abundance

```{r}
# plot absolute abundance of ITS2 sequences
# Indonesia
# upload data
its2_sequences <- read.table('/Users/MarlenSchlotheuber/Nextcloud/10. projects/11. PhD/10. projects/Intraspecific bleaching project/Project planning/MYS_2024/Fieldwork/Data_26042024/SymPortal/Sequencing_data_ITS2_SymPortal/ITS2_SymPortal_data/202408_ID_ITS2_data/post_med_seqs/582_20250509T135802_DBV_20250510T135051.seqs.absolute.abund_and_meta.txt', header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# remove columns with meta data only keep sample name and ITS2 sequences
its2_sequences_filtered <- its2_sequences[, -c(1, 3:40)]

# remove rows with NA
its2_sequences_filtered <- its2_sequences_filtered %>%
  filter(sample_name != "" & !is.na(sample_name))

# Convert to long format for ggplot
its2_sequences_filtered_long <- its2_sequences_filtered %>%
  pivot_longer(
    cols = -sample_name,          # all columns except SampleID
    names_to = "Sequence",
    values_to = "AbsoluteAbundance"
  )

# make sure the order of sequences is correct (highest abundnace at the bottom of the barplot)
seq_order <- its2_sequences_filtered_long %>%
  group_by(Sequence) %>%
  summarize(total_abundance = sum(AbsoluteAbundance, na.rm = TRUE)) %>%
  arrange(desc(total_abundance)) %>%      # highest abundance first
  pull(Sequence)

# Set factor levels so that highest abundance is at the bottom
its2_sequences_filtered_long$Sequence <- factor(
  its2_sequences_filtered_long$Sequence,
  levels = seq_order
)

# color palette
palette_88 <- colorRampPalette(c("magenta1", "darkorange1", "green", "blue"))(88)

# plot
ID_ITS2seq <- ggplot(its2_sequences_filtered_long, aes(x = sample_name, y = AbsoluteAbundance, fill = Sequence)) +
  geom_bar(stat = "identity", position = position_stack(reverse = TRUE)) +
  scale_fill_manual(values = palette_88) +
  labs(x = "Sample", y = "Relative Abundance") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none"   # hide legend
  )
ID_ITS2seq 

# Malaysia
# upload data
its2_sequences <- read.table('/Users/MarlenSchlotheuber/Nextcloud/10. projects/11. PhD/10. projects/Intraspecific bleaching project/Project planning/MYS_2024/Fieldwork/Data_26042024/SymPortal/Sequencing_data_ITS2_SymPortal/ITS2_SymPortal_data/202404_MYS_ITS2_data/post_med_seqs/560_20250219T125105_DBV_20250221T224139.seqs.absolute.abund_and_meta.txt', header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# remove columns with meta data only keep sample name and ITS2 sequences
its2_sequences_filtered <- its2_sequences[, -c(1, 3:39)]

# remove rows with NA
its2_sequences_filtered <- its2_sequences_filtered %>%
  filter(sample_name != "" & !is.na(sample_name))

# Convert to long format for ggplot
its2_sequences_filtered_long <- its2_sequences_filtered %>%
  pivot_longer(
    cols = -sample_name,          # all columns except SampleID
    names_to = "Sequence",
    values_to = "AbsoluteAbundance"
  )

# make sure the order of sequences is correct (highest abundnace at the bottom of the barplot)
seq_order <- its2_sequences_filtered_long %>%
  group_by(Sequence) %>%
  summarize(total_abundance = sum(AbsoluteAbundance, na.rm = TRUE)) %>%
  arrange(desc(total_abundance)) %>%      # highest abundance first
  pull(Sequence)

# Set factor levels so that highest abundance is at the bottom
its2_sequences_filtered_long$Sequence <- factor(
  its2_sequences_filtered_long$Sequence,
  levels = seq_order
)

# color palette
palette_88 <- colorRampPalette(c("magenta1", "darkorange1", "green", "blue"))(88)

# plot
ggplot(its2_sequences_filtered_long, aes(x = sample_name, y = AbsoluteAbundance, fill = Sequence)) +
  geom_bar(stat = "identity", position = position_stack(reverse = TRUE)) +
  scale_fill_manual(values = palette_88) +
  labs(x = "Sample", y = "Relative Abundance") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none"   # hide legend
  )
```

# 16S - bacterial community composition

16S sequencing data was processed in R studio using the DADA2 pipeline. DADA2 (Callahan et al., 2016) filtered and trimmed reads, calculated error rates, inferred sequence variants, merged paired reads, and removed chimeras. Taxonomy was assigned by using the Silva v138 dataset from the Silva reference database (Benjamin, 2024). The top 20 bacterial families across the colony surface were visualised using R studio. <https://benjjneb.github.io/dada2/tutorial.html>

### Malaysia: DADA2

```{r}
install.packages("seqinr")
library(seqinr)

# Save as FASTA file
write.fasta(sequences = as.list(raw_fasta),
            names = names(raw_fasta),
            file.out = "top_ASVs.fasta")



# load packages
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("dada2")

BiocManager::install("dada2", force = TRUE)
install.packages(c("cubature", "ape", "reshape2", "edgeR", "plyr"))

install.packages("dada2")
install.packages("Rcpp")
library(dada2)
library(Rcpp)

getwd()
# set working directory
# read file
cat("Reading files") # prints a message to indicate the script is starting to read files
path <- "." # current working directory
list.files(path) # lists all files in the current directory 

# Forward and reverse fastq filenames have format: 
fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))
# searches for files matching specific patterns _R1_001...and sorts files, get a sorted list of file paths for forward (fnFs) and reverse reads (fnRs)

# Extract sample names, assuming filenames have format
sample.names <- sapply(strsplit(basename(fnFs), "__|_0"), function(x) x[2])

# basename(fnFs): Removes the directory path and keeps only the file names. Example: "/path/to/SAMPLE1_R1_001.fastq" → "SAMPLE1_R1_001.fastq"
# strsplit(..., "_"): Splits each filename into parts using the underscore (_) as the delimiter. Example: "SAMPLE1_R1_001.fastq" → c("SAMPLE1", "R1", "001.fastq")
# sapply(..., [, 1): Extracts the first element from each split result (the part before the first underscore). Example: c("SAMPLE1", "R1", "001.fastq") → "SAMPLE1"

cat("Processing",length(sample.names),"samples:", sample.names)
# outputs a summary of the number of samples detected and their names, legth(): counts the number, cat():prints the message

#Inspect read quality profiles
plotQualityProfile(fnFs[1:2])
# generates a quality profile plot for the first two forward read files in the fnFs list, DADA2 package in R for processing high-throughput sequence data, to visualize the quality score of sequencing reads across all base positions in a FASTQ file
# selected the first two 1:2 forward read files in fnFs 

# gray-scale: heat map of the frequency of each quality score at each base position
# green line: mean quality score at each position
# orange line: quartiles of the quality score distribution
# red line: scaled proportion of reads that extend to at least that position

# forwards reads
# --> Truncate the forward reads at position 280

plotQualityProfile(fnRs[1:2])

# reversed reads
# worse quality, especially at the end which is common in Illumina sequencing
# DADA2 incorporates quality information into its error model which makes the algorithm robust to lower quality sequence

# reversed reads
# --> Truncate the forward reads at position 200


```

### Filtering and Trimming

This pipeline ensures that only high-quality reads are retained for downstream analysis, which is crucial for the accuracy of biological insights derived from sequencing data.

```{r}
# filter and trimm 
cat("Filtering and trimming")

# prints a message to indicate that the filtering and triming steps are about to begin in the sequencing data processing pipeline

# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

### set parameters  maxN=0 (DADA2 requires no Ns), truncQ=2, rm.phix=TRUE and maxEE=2. The maxEE parameter sets the maximum number of “expected errors” allowed in a read, which is a better filter than simply averaging quality scores.
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(280,200),maxN=0, maxEE=c(2,2), truncQ=5, rm.phix=TRUE, compress=TRUE, multithread=TRUE) 
head(out)

```

### Learn error rates & sample interference

```{r}
# learn error rates
# The error model is learned using the filtered and trimmed sequences, which helps in determining the probability of errors at each base position across the reads.Forward (errF) and reverse (errR) errors are learned separately because the forward and reverse reads often have different error profiles.

cat("Learning error rates")
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)
# This function in DADA2 is used to learn the error rates for sequencing reads.
# multithread=TRUE: Enables multi-threading to speed up the error rate learning process if the system supports it.

# visualize the esstimated error rates as a sanity check
plotErrors(errF, nominalQ = T)
# everything looks resonable (black line follows the points, the model fits) and we proceed with confidence

```

### Sample Inference

```{r}
# This step is essential for transforming the raw sequencing data into biologically meaningful information by correcting errors and identifying the true sequences present in your samples
# apply the core sample inference algorithm to the filtered and trimmed sequence data.
dadaFs <- dada(filtFs, err=errF, multithread=TRUE, pool = T)
dadaRs <- dada(filtRs, err=errR, multithread=TRUE, pool = T)
# It uses the error rates learned from the reverse reads (errR) to apply the error-correction model and infer the true biological sequences.

#Inspecting the returned dada-class object:
dadaFs[[1]]
dadaRs[[1]]
```

### Merging paired reads and checking on read length

```{r}
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
# merges forward (F) and reveres (R) sequences and filtered data in one column
# Inspect the merger data.frame from the first sample
head(mergers[[1]])

# Construct sequencing table
seqtab <- makeSequenceTable(mergers)
# creates a sequence table from the merged paired-end reads. A sequence table is a matrix where: rows represent individual samples, columns represent unique DNA sequences, entries contain the counts of each sequence in each sample
dim(seqtab)
# returns the dimensions of the sequence table (1st value number of rows (sample), 2nd value is the number of columns (unique sequences))

# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))
# looks good
# check visually how the districution of read length is
# Get the read length distribution
length_table <- table(nchar(getSequences(seqtab)))

# Convert to data frame for plotting
length_df <- as.data.frame(length_table)
colnames(length_df) <- c("Length", "Count")

# Make sure Length is numeric
length_df$Length <- as.numeric(as.character(length_df$Length))

ggplot(length_df, aes(x = Length, y = Count)) +
  geom_col(fill = "steelblue") +  # vertical bars
  labs(title = "Distribution of Read Lengths",
       x = "Read Length (bp)",
       y = "Number of Reads") +
  xlim(280, 400)+
  theme_minimal()

```

### Remove Chimeras

```{r}
# Remove chimeras  
# Chimeric sequences are artifacts of the PCR amplification process. These sequences do not represent real biological sequences but rather result from errors in the molecular processes used to prepare and sequence DNA.
# Chimeric sequences can inflate diversity estimates. Removing them ensures that the analysis focuses only on true biological sequences.
# This code removes chimeric sequences from the sequence table, calculates the dimensions of the updated sequence table and determines the proportion of non-chimeric sequences retained.

seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
# removeBim..: identifies and removes chimeric sequences from the sequence table
# method: chimera removal method that identifies chimeras based on agreement across samples
# multithread=T: enables multi-threading for faster computation
# verbose=T: provides detailed output about the chimera removal process including the number of chimeric sequences identified.
dim(seqtab.nochim)

## calculate the  frequency of chimeric sequences 
sum(seqtab.nochim)/sum(seqtab)
# 0.885 = 88.5% of your reads were retained, only 11.5% of reads were discarded as chimeras


# Track reads through the pipeline 
# This block of code tracks the number of reads at each stage of the DADA2 pipeline for each sample. This tracking helps monitor data loss, identify problematic samples, and ensure the pipeline is functioning as expected.
getN <- function(x) sum(getUniques(x))

# custom function that calculates the total number of unique reads (amplicon sequences)
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# combines data from multiple steps into a single table
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)

# Purpose
# Quality Control: This table allows you to monitor read retention at each step of the pipeline, ensuring that no unexpected data loss occurs.
# Troubleshooting: If a sample shows a significant drop in reads at a specific step, it may i indicate a problem (e.g., poor quality reads or insufficient overlap for merging).
# Summary: Provides an overview of pipeline performance, enabling comparison across samples.
```

### Assign Taxonomy

```{r}
# Assign taxonomy
# This block of code assigns taxonomic classifications to the non-chimeric sequences (ASVs) in your dataset, enabling the identification of the taxa represented by the sequences

taxa <- assignTaxonomy(seqtab.nochim,"/Users/MarlenSchlotheuber/Nextcloud/10. projects/11. PhD/10. projects/Intraspecific bleaching project/Project planning/MYS_2024/Fieldwork/Data_26042024/16S_Bacteria_data/16S_files/silva_nr_v138_train_set.fa.gz" , multithread=TRUE)
# uses a reference database to classify each unique sequence (ASV) in the seqtab.nochim table to a taxonomic level
# rows: correspond to unique sequences (ASV)

##optional: make species level assignments based on exact matching
taxa <- addSpecies(taxa, "/Users/MarlenSchlotheuber/Nextcloud/10. projects/11. PhD/10. projects/Intraspecific bleaching project/Project planning/MYS_2024/Fieldwork/Data_26042024/16S_Bacteria_data/16S_files/silva_species_assignment_v138.fa.gz")
                   
#inspect the taxonomic assignments:
taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)
```

### Prepare Data for plotting

plots the 20 most abundant bacterial families

```{r}
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

# Install Biostrings
BiocManager::install("Biostrings")
# Install phyloseq
BiocManager::install("phyloseq")

library(phyloseq)
library(Biostrings)
library(reshape2)
library(ggplot2)
library(scales)
library(dplyr)
library(gridExtra)

# add needed info to data frame
samples.out_MYS <- rownames(seqtab.nochim)
info_split_MYS <- strsplit(samples.out_MYS, "_")
Coral_Species <- sapply(info_split_MYS, function(x) sub("[0-9]+", "", x[1]))
Genotype <- sapply(info_split_MYS, function(x) gsub("[^0-9]", "", x[1]))
Orientation <- sapply(info_split_MYS, function(x) sub("[0-9]+", "", x[2]))
Location <- sapply(info_split_MYS, function(x) gsub("[^0-9]", "", x[2]))
samdf_MYS <- data.frame(Coral_Species, Genotype, Orientation, Location, stringsAsFactors = FALSE)
rownames(samdf_MYS) <- samples.out_MYS

# also add the column Location
# Convert Location to numeric for mapping
samdf_MYS$Location <- as.numeric(samdf_MYS$Location)

# Map numeric sampling location to descriptive Location
location_map_MYS <- c("top", "upper middle", "lower middle", "bottom")

samdf_MYS$Location <- location_map_MYS[samdf_MYS$Location]

rownames(samdf_MYS) <- samples.out_MYS

# colony 2 has two more sampling regions but for comparison I will remove it for now
# Remove rows with NA in any column
samdf_MYS <- na.omit(samdf_MYS)

# construct a phyloseq object directly from the dada2 outputs
ps_MYS <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               sample_data(samdf_MYS), 
               tax_table(taxa))

# its easier to work with the short ASV names than with the full DNA sequence
# store DNA sequences of our ASC in the refseq slot of the phyloseq object
dna_MYS <- Biostrings::DNAStringSet(taxa_names(ps_MYS))
names(dna_MYS) <- taxa_names(ps_MYS)
ps_MYS <- merge_phyloseq(ps_MYS, dna_MYS)
taxa_names(ps_MYS) <- paste0("ASV", seq(ntaxa(ps_MYS)))
ps_MYS
```

### Indonesia: DADA2

Data processing with [DADA2](https://github.com/benjjneb/dada2) to infer ASVs

```{r}
install.packages("seqinr")
library(seqinr)

# Save as FASTA file
write.fasta(sequences = as.list(raw_fasta),
            names = names(raw_fasta),
            file.out = "top_ASVs.fasta")



# load packages
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("dada2")

BiocManager::install("dada2", force = TRUE)
install.packages(c("cubature", "ape", "reshape2", "edgeR", "plyr"))

install.packages("dada2")
install.packages("Rcpp")
library(dada2)
library(Rcpp)

getwd()
# set working directory
# read file
cat("Reading files") # prints a message to indicate the script is starting to read files
path <- "." # current working directory
list.files(path) # lists all files in the current directory 

# Forward and reverse fastq filenames have format: 
fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))
# searches for files matching specific patterns _R1_001...and sorts files, get a sorted list of file paths for forward (fnFs) and reverse reads (fnRs)

# Extract sample names, assuming filenames have format
sample.names <- sapply(strsplit(basename(fnFs), "__|_0"), function(x) x[2])

# basename(fnFs): Removes the directory path and keeps only the file names. Example: "/path/to/SAMPLE1_R1_001.fastq" → "SAMPLE1_R1_001.fastq"
# strsplit(..., "_"): Splits each filename into parts using the underscore (_) as the delimiter. Example: "SAMPLE1_R1_001.fastq" → c("SAMPLE1", "R1", "001.fastq")
# sapply(..., [, 1): Extracts the first element from each split result (the part before the first underscore). Example: c("SAMPLE1", "R1", "001.fastq") → "SAMPLE1"

cat("Processing",length(sample.names),"samples:", sample.names)
# outputs a summary of the number of samples detected and their names, legth(): counts the number, cat():prints the message

#Inspect read quality profiles
plotQualityProfile(fnFs[1:2])
# generates a quality profile plot for the first two forward read files in the fnFs list, DADA2 package in R for processing high-throughput sequence data, to visualize the quality score of sequencing reads across all base positions in a FASTQ file
# selected the first two 1:2 forward read files in fnFs 

# gray-scale: heat map of the frequency of each quality score at each base position
# green line: mean quality score at each position
# orange line: quartiles of the quality score distribution
# red line: scaled proportion of reads that extend to at least that position

# forwards reads
# --> Truncate the forward reads at position 280

plotQualityProfile(fnRs[1:2])

# reversed reads
# worse quality, especially at the end which is common in Illumina sequencing
# DADA2 incorporates quality information into its error model which makes the algorithm robust to lower quality sequence

# reversed reads
# --> Truncate the forward reads at position 200


```

### Filtering and Trimming

This pipeline ensures that only high-quality reads are retained for downstream analysis, which is crucial for the accuracy of biological insights derived from sequencing data.

```{r}
# filter and trimm 
cat("Filtering and trimming")

# prints a message to indicate that the filtering and triming steps are about to begin in the sequencing data processing pipeline

# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

### set parameters  maxN=0 (DADA2 requires no Ns), truncQ=2, rm.phix=TRUE and maxEE=2. The maxEE parameter sets the maximum number of “expected errors” allowed in a read, which is a better filter than simply averaging quality scores.
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(280,200),maxN=0, maxEE=c(2,2), truncQ=5, rm.phix=TRUE, compress=TRUE, multithread=TRUE) 

head(out)
```

### Learn error rates & sample interference

```{r}
# learn error rates
# The error model is learned using the filtered and trimmed sequences, which helps in determining the probability of errors at each base position across the reads.Forward (errF) and reverse (errR) errors are learned separately because the forward and reverse reads often have different error profiles.

cat("Learning error rates")
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)
# This function in DADA2 is used to learn the error rates for sequencing reads.
# multithread=TRUE: Enables multi-threading to speed up the error rate learning process if the system supports it.

# visualize the esstimated error rates as a sanity check
plotErrors(errF, nominalQ = T)
# everything looks resonable (black line follows the points, the model fits) and we proceed with confidence
```

### Sample Inference

```{r}
# This step is essential for transforming the raw sequencing data into biologically meaningful information by correcting errors and identifying the true sequences present in your samples
# apply the core sample inference algorithm to the filtered and trimmed sequence data.
dadaFs <- dada(filtFs, err=errF, multithread=TRUE, pool = T)
dadaRs <- dada(filtRs, err=errR, multithread=TRUE, pool = T)
# It uses the error rates learned from the reverse reads (errR) to apply the error-correction model and infer the true biological sequences.

#Inspecting the returned dada-class object:
dadaFs[[1]]
dadaRs[[1]]
```

### Merging paired reads and checking on read length

```{r}
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
# merges forward (F) and reveres (R) sequences and filtered data in one column
# Inspect the merger data.frame from the first sample
head(mergers[[1]])

# Construct sequencing table
seqtab <- makeSequenceTable(mergers)
# creates a sequence table from the merged paired-end reads. A sequence table is a matrix where: rows represent individual samples, columns represent unique DNA sequences, entries contain the counts of each sequence in each sample
dim(seqtab)
# returns the dimensions of the sequence table (1st value number of rows (sample), 2nd value is the number of columns (unique sequences))

# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))
# looks good
# check visually how the districution of read length is
# Get the read length distribution
length_table <- table(nchar(getSequences(seqtab)))

# Convert to data frame for plotting
length_df <- as.data.frame(length_table)
colnames(length_df) <- c("Length", "Count")

# Make sure Length is numeric
length_df$Length <- as.numeric(as.character(length_df$Length))

ggplot(length_df, aes(x = Length, y = Count)) +
  geom_col(fill = "steelblue") +  # vertical bars
  labs(title = "Distribution of Read Lengths",
       x = "Read Length (bp)",
       y = "Number of Reads") +
  xlim(280, 400)+
  theme_minimal()
```

### Remove Chimeras

```{r}
# Remove chimeras  
# Chimeric sequences are artifacts of the PCR amplification process. These sequences do not represent real biological sequences but rather result from errors in the molecular processes used to prepare and sequence DNA.
# Chimeric sequences can inflate diversity estimates. Removing them ensures that the analysis focuses only on true biological sequences.
# This code removes chimeric sequences from the sequence table, calculates the dimensions of the updated sequence table and determines the proportion of non-chimeric sequences retained.

seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
# removeBim..: identifies and removes chimeric sequences from the sequence table
# method: chimera removal method that identifies chimeras based on agreement across samples
# multithread=T: enables multi-threading for faster computation
# verbose=T: provides detailed output about the chimera removal process including the number of chimeric sequences identified.
dim(seqtab.nochim)

## calculate the  frequency of chimeric sequences 
sum(seqtab.nochim)/sum(seqtab)
# 0.654 = 65.4% of your reads were retained, only the rest% of reads were discarded as chimeras


# Track reads through the pipeline 
# This block of code tracks the number of reads at each stage of the DADA2 pipeline for each sample. This tracking helps monitor data loss, identify problematic samples, and ensure the pipeline is functioning as expected.
getN <- function(x) sum(getUniques(x))

# custom function that calculates the total number of unique reads (amplicon sequences)
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# combines data from multiple steps into a single table
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)

# Purpose
# Quality Control: This table allows you to monitor read retention at each step of the pipeline, ensuring that no unexpected data loss occurs.
# Troubleshooting: If a sample shows a significant drop in reads at a specific step, it may i indicate a problem (e.g., poor quality reads or insufficient overlap for merging).
# Summary: Provides an overview of pipeline performance, enabling comparison across samples.
```

### Assign Taxonomy

```{r}
{r}
# Assign taxonomy
# This block of code assigns taxonomic classifications to the non-chimeric sequences (ASVs) in your dataset, enabling the identification of the taxa represented by the sequences

taxa <- assignTaxonomy(seqtab.nochim,"/Users/MarlenSchlotheuber/Nextcloud/10. projects/11. PhD/10. projects/Intraspecific bleaching project/Project planning/MYS_2024/Fieldwork/Data_26042024/16S_Bacteria_data/16S_files/silva_nr_v138_train_set.fa.gz" , multithread=TRUE)
# uses a reference database to classify each unique sequence (ASV) in the seqtab.nochim table to a taxonomic level
# rows: correspond to unique sequences (ASV)

##optional: make species level assignments based on exact matching
# taxa <- addSpecies(taxa, "/Users/MarlenSchlotheuber/Nextcloud/10. projects/11. PhD/10. projects/Intraspecific bleaching project/Project planning/MYS_2024/Fieldwork/Data_26042024/16S_Bacteria_data/16S_files/silva_species_assignment_v138.fa.gz")
                   
#inspect the taxonomic assignments:
taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)
```

### Prepare Data for plotting

```{r}
{r}
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

# Install Biostrings
BiocManager::install("Biostrings")
# Install phyloseq
BiocManager::install("phyloseq")

library(phyloseq)
library(Biostrings)
library(reshape2)
library(ggplot2)
library(scales)
library(dplyr)
library(gridExtra)

# add needed info to data frame
samples.out <- rownames(seqtab.nochim)
info_split <- strsplit(samples.out, "_")
Species <- sapply(info_split, function(x) sub("[0-9]+", "", x[2]))
Genotype <- sapply(info_split, function(x) gsub("[^0-9]", "", x[2]))
SamplingLocation <- sapply(info_split, function(x) x[3])
samdf <- data.frame(Species, Genotype, SamplingLocation, stringsAsFactors = FALSE)
rownames(samdf) <- samples.out

# construct a phyloseq object directly from the dada2 outputs
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               sample_data(samdf), 
               tax_table(taxa))

# its easier to work with the short ASV names than with the full DNA sequence
# store DNA sequences of our ASC in the refseq slot of the phyloseq object
dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps_ID <- merge_phyloseq(ps, dna)
taxa_names(ps_ID) <- paste0("ASV", seq(ntaxa(ps_ID)))
ps_ID
```

### Barplot: Bacterial Family

```{r}
# Malaysia
# plot it again but with all other bacteria families in grey
# Groups ASVs by Family level, Counts are summed up for each Family
ps_MYS_family <- tax_glom(ps_MYS, taxrank = "Family")

# check what happens if I work with the df
df_ps_MYS_family <- psmelt(ps_MYS_family)

# clean sample column name
df_ps_MYS_family$Sample <- sub("_S.*", "", df_ps_MYS_family$Sample)

# 3. Normalize abundance to relative abundance per sample
df_ps_MYS_family <- df_ps_MYS_family %>%
  group_by(Sample) %>%
  mutate(RelAbundance = Abundance / sum(Abundance)) %>%
  ungroup()

# total number of reads
total_reads_MYS <- df_ps_MYS_family %>%
  summarise(TotalReads = sum(Abundance, na.rm = TRUE))

print(total_reads_MYS)
# total reads: 198.725.1	

# how many ASV do we have?
ntaxa(ps_MYS)
# 4417 MYS
length(unique(df_ps_MYS_family$OTU))
# 274

# how many families do we have?
length(unique(tax_table(ps)[, "Family"]))
#  43

# top 20
top20_families <- df_ps_MYS_family %>%
  group_by(Family) %>%
  summarise(TotalRelAbundance = sum(RelAbundance)) %>%
  arrange(desc(TotalRelAbundance)) %>%
  slice_head(n = 20) %>%
  pull(Family)

# Group all non-top 20 families as "Other"
df_ps_MYS_top20_family <- df_ps_MYS_family %>%
  mutate(Family = ifelse(Family %in% top20_families, as.character(Family), "Other"))

# check number of top 20 families
df_ps_MYS_top20_family %>% distinct(Family) %>% count()
# 21 (including others)
unique(df_ps_MYS_top20_family$Family)

# check on % of each family of the top 20
# Step 1: Sum abundance per family
family_sums <- aggregate(RelAbundance ~ Family, data = df_ps_MYS_top20_family, sum)

# Step 2: Compute total abundance
total_abundance <- sum(family_sums$RelAbundance)

# Step 3: Calculate percent abundance
family_sums$Percent_Abundance <- (family_sums$RelAbundance / total_abundance) * 100

# Set factor levels so "Other" appears last in stack
df_ps_MYS_top20_family <- df_ps_MYS_top20_family %>%
  mutate(Family = factor(Family, levels = c("Other", top20_families)))

# PLOT
# organize Genotype
df_ps_MYS_top20_family <- df_ps_MYS_top20_family %>%
  mutate(
    Genotype = factor(Genotype, levels = c("1", "6", "2", "3", "5", "4")))

# Orientation
df_ps_MYS_top20_family <- df_ps_MYS_top20_family %>%
  mutate(
    Orientation = factor(Orientation, levels = c("W", "E")))

# write the y-axis labesl in capital letters
df_ps_MYS_top20_family <- df_ps_MYS_top20_family %>%
  mutate(Location = str_to_title(Location))

# Location
df_ps_MYS_top20_family <- df_ps_MYS_top20_family %>%
  mutate(
    Location = factor(Location, levels = c("Bottom", "Lower Middle", "Upper Middle", "Top")))

# wehn combining MY and ID make sure the colorcode aligns than use this code if not continue as below
# for plotting MY and ID together asign each fmaily to a colors for plotting later
# 9 overlapping bacterial family names in total
# 31 unique family names, so in total we need 31 unique colors

# first list bacterial family names
# 1. Define all unique bacterial families from both MY and ID (31 unique)
all_families <- c(
  # 9 overlapping families first (same order for consistency)
  "Halomonadaceae",
  "Marinobacteraceae",
  "Chlorobiaceae",
  "Alcanivoracaceae1",
  "Idiomarinaceae",
  "Alteromonadaceae",
  "Rhodobacteraceae",
  "Vibrionaceae",
  "Flavobacteriaceae",
  
  # Families unique to MY
  "Endozoicomonadaceae",
  "Moraxellaceae",
  "Spirochaetaceae",
  "Marinilabiliaceae",
  "Hungateiclostridiaceae",
  "Cyclobacteriaceae",
  "Propionibacteriaceae",
  "Desulfosarcinaceae",
  "Desulfobacteraceae",
  "Clade_I",
  "Pseudomonadaceae",
  
  # Families unique to ID
  "Pseudoalteromonadaceae",
  "Colwelliaceae",
  "Marinomonadaceae",
  "Peptostreptococcales-Tissierellales_fa",
  "Arcobacteraceae",
  "Nitrincolaceae",
  "Kordiimonadaceae",
  "Stappiaceae",
  "Cryomorphaceae",
  "Marinifilaceae",
  "Saprospiraceae"
)

# assign colors to each family
kelly_colors <- c(
  "#00538A", # Strong Blue
  "#990000", # Vivid Red
  "#FF6800", # Vivid Orange
  "#007D35", # Vivid Green
  "#FFB300", # Vivid Yellow
  "#803E75", # Strong Purple
  "#A6BDD7", # Very Light Blue
  "#CC6699", # Medium Gray
  "#F6768E", # Strong Purplish Pink
  "#53377A", # Strong Violet
  "#FF0000", # Strong Purplish Red
  "#FFEA00", # Vivid Greenish Yellow
  "#99CCFF", # Light Blue
  "#93AA00", # Yellowish Green
  "#593318", # Deep Yellowish Brown
  "#3399CC", # Vivid Reddish Orange
  "#003300", # Dark Olive Green
  "#33CC99", # Aquamarine
  "#CC99FF", # Lavender
  "#CC3333", # Soft Red
  "#66CC66", # Soft Green
  "#FFCC99", # Warm Yellow
  "#6699CC", # Muted Blue
  "#996633", # Brown
  "#FF6699", # Light Pink
  "#1B9E77", # Teal
  "#D95F02", # Orange
  "#7570B3", # Indigo
  "#E7298A", # Magenta
  "#66A61E", # Olive Green
  "#1F78B4"  # Medium Blue
)

# 3. Create a named vector mapping families to colors
family_colors <- setNames(kelly_colors, all_families)

# 4. Add 'Other' color for families not listed explicitly
family_colors <- c(family_colors, Other = "grey70")

# colors by family
df_ps_MYS_top20_family <- df_ps_MYS_top20_family %>%
  mutate(Family = factor(Family, levels = c("Other", all_families)))

# plot
bacteria_family_plot_MYS <- ggplot(df_ps_MYS_top20_family, aes(x = RelAbundance * 100, y = Location, fill = Family)) +
  geom_bar(stat = "identity", position = "stack") +
  facet_grid(Genotype ~ Orientation, scales = "free_y") +
  scale_x_continuous(labels = scales::percent_format(scale = 1)) +
  scale_fill_manual(values = family_colors) +
  theme_minimal(base_size = 13) +
  theme(
    strip.text = element_text(face = "bold", size = 13),
    legend.position = "bottom",
    legend.text = element_text(size = 13),
    legend.title = element_text(size = 13),
    axis.text = element_text(size = 13),
    axis.title = element_text(size = 13),
    panel.grid = element_blank(),
    panel.border = element_blank()
  ) +
  guides(fill = guide_legend(ncol = 3)) +
  labs(
    x = "Relative Abundance (%)",
    y = "Sampling Location",
    fill = "Bacterial Family"
  )

bacteria_family_plot_MYS

# save plots of MY
# save as image
ggsave("/Users/MarlenSchlotheuber/Nextcloud/10. projects/11. PhD/10. projects/Intraspecific bleaching project/Project planning/MYS_2024/Fieldwork/Data_26042024/Manuscript/R_code/Bacterial_MY.png", plot = bacteria_family_plot_MYS, width = 8, height = 8, dpi = 300, bg = "white")

# save as pdf
ggsave("/Users/MarlenSchlotheuber/Nextcloud/10. projects/11. PhD/10. projects/Intraspecific bleaching project/Project planning/MYS_2024/Fieldwork/Data_26042024/Manuscript/R_code/Bacterial_MYS.pdf", 
       plot = bacteria_family_plot_MYS, 
       width = 12, 
       height = 10, 
       dpi = 300, 
       bg = "white")


# Indonesia
# check overall numbers 
# Groups ASVs by Family level, Counts are summed up for each Family
ps_ID_family <- tax_glom(ps_ID, taxrank = "Family")

# check what happens if I work with the df
df_ps_ID_family <- psmelt(ps_ID_family)

#  Normalize abundance to relative abundance per sample
df_ps_ID_family <- df_ps_ID_family %>%
  group_by(Sample) %>%
  mutate(RelAbundance = Abundance / sum(Abundance)) %>%
  ungroup()

# total number of reads
total_reads <- df_ps_ID_family %>%
  summarise(TotalReads = sum(Abundance, na.rm = TRUE))

print(total_reads)
# total reads:502.311	

# how many ASV do we have?
ntaxa(ps)
# 678

# how many families do we have?
length(unique(tax_table(ps)[, "Family"]))
# 44

# top 20
top20_families <- df_ps_ID_family %>%
  group_by(Family) %>%
  summarise(TotalRelAbundance = sum(RelAbundance)) %>%
  arrange(desc(TotalRelAbundance)) %>%
  slice_head(n = 20) %>%
  pull(Family)

# check top 20 families
top20_families

# Group all non-top 20 families as "Other"
df_ID_top20_family <- df_ps_ID_family %>%
  mutate(Family = ifelse(Family %in% top20_families, Family, "Other"))

# check number of bacteria families
df_ID_top20_family %>%
  distinct(Family) %>%
  count()
# 21 which is correct because I have 20 families and the rest is other

# check on % of each family of the top 20
# Step 1: Sum abundance per family
family_sums <- aggregate(RelAbundance ~ Family, data = df_ID_top20_family, sum)

# Step 2: Compute total abundance
total_abundance <- sum(family_sums$RelAbundance)

# Step 3: Calculate percent abundance
family_sums$Percent_Abundance <- (family_sums$RelAbundance / total_abundance) * 100

# Plot
# organize Genotype
df_ID_top20_family <- df_ID_top20_family %>%
  mutate(
    Genotype = factor(Genotype, levels = sort(unique(Genotype))))

# order of colonies
df_ID_top20_family <- df_ID_top20_family %>%
  mutate(
    Genotype = factor(Genotype, levels = c("01", "02", "03", "04", "05", "06", "07")))

# Location
df_ID_top20_family <- df_ID_top20_family %>%
  mutate(
    SamplingLocation = factor(SamplingLocation, levels = c("bottom", "top")))

df_ID_top20_family <- df_ID_top20_family %>%
  mutate(
    SamplingLocation = factor(SamplingLocation,
                              levels = c("bottom", "top"),
                              labels = c("Bottom", "Top"))
  )

# Set factor levels so "Other" appears last in stack
df_ID_top20_family <- df_ID_top20_family %>%
  mutate(Family = factor(Family, levels = c("Other", all_families)))

#plot
bacteria_family_plot_ID <- ggplot(df_ID_top20_family, aes(x = RelAbundance * 100, y = SamplingLocation, fill = Family)) +
  geom_bar(stat = "identity", position = "stack") +
  facet_wrap(~Genotype, scales = "free_y", ncol = 1) +
  scale_x_continuous(labels = scales::percent_format(scale = 1)) +
  scale_fill_manual(values = family_colors) +
  theme_minimal(base_size = 13) +
  theme(
    strip.text = element_blank(),  # Remove genotype labels
    legend.position = "bottom",    # Move legend underneath the plot
    legend.text = element_text(size = 13),
    legend.title = element_text(size = 13),
    axis.text = element_text(size = 13),
    axis.title = element_text(size = 13),
    panel.grid = element_blank(),
    legend.box = "horizontal"      # Layout legend horizontally
  ) +
  guides(fill = guide_legend(ncol = 3)) +  # Adjust number of rows as needed
  labs(
    x = "Relative Abundance[%]",
    y = "Sampling Location",
    fill = "Bacterial Family"
  )

bacteria_family_plot_ID

# save plots
# ID
# save as image
ggsave("/Users/MarlenSchlotheuber/Nextcloud/10. projects/11. PhD/10. projects/Intraspecific bleaching project/Project planning/MYS_2024/Fieldwork/Data_26042024/Manuscript/R_code/Fig7_ID_bacteria_family_16S_data.png", plot = bacteria_family_plot_ID, width = 8, height = 8, dpi = 300, bg = "white")

# save as pdf
ggsave("/Users/MarlenSchlotheuber/Nextcloud/10. projects/11. PhD/10. projects/Intraspecific bleaching project/Project planning/MYS_2024/Fieldwork/Data_26042024/Manuscript/R_code/Fig7_ID_bacteria_family_16S_data.pdf", 
       plot = bacteria_family_plot_ID, 
       width = 12, 
       height = 10, 
       dpi = 300, 
       bg = "white")


```

## **Fig. 4 Bacterial community composition across the surface of massive Porites sp. colonies from Indonesia and Malaysia**

```{r}
# combine ID and MYS
# Combine plots side-by-side with equal width
final_16S_combined_plot <-  bacteria_family_plot_ID + bacteria_family_plot_MYS +
  plot_layout(ncol = 2, widths = c(1, 1)) +
  plot_annotation(
    theme = theme(
      plot.margin = margin(10, 10, 10, 10)
    )
  )

final_16S_combined_plot

# save plot and edit in affinity designer
# Save as PDF
ggsave(
  filename = "16S_MYS_ID.pdf",
  plot = final_16S_combined_plot,
  width = 11.7,      # A4 width in inches
  height = 10,     # A4 height in inches
  units = "in"
)
```
