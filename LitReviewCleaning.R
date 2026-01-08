##Installand load  packages

install.packages("sf")
install.packages("ggplot2")
install.packages("dplyr")
install.packages("ggnewscale")
install.packages("tidyr")
library(dplyr)
library(sf)
library(ggplot2)
library(ggnewscale)
library(tidyr)

#Load data
# 1. Path and File List
folder_path <- "C:/Users/alexa/OneDrive - University of Guelph/Documents/WEB lab/North American Bat Polygons"
shp_files <- list.files(folder_path, pattern = "\\.shp$", full.names = TRUE)

# 2. Optimized Import: Simplify and Standardize DURING the loop
all_bats_list <- lapply(shp_files, function(f) {
  # Read file
  s <- st_read(f, quiet = TRUE)
  
  # A. Simplify immediately (dTolerance = 1000 ignores details smaller than 1km)
  # This reduces the number of "points" R has to draw by up to 90%
  s <- st_simplify(s, preserveTopology = TRUE, dTolerance = 1000)
  
  # B. Standardize CRS now (WGS84)
  # ggplot is MUCH faster if it doesn't have to convert coordinates on the fly
  s <- st_transform(s, 4326)
  
  # C. Keep only what we need (The name and the shape)
  s$species_name <- basename(f)
  return(s[, c("species_name", "geometry")])
})

# 3. Combine the simplified layers
all_bats <- do.call(rbind, all_bats_list)

# 4. Fast Plotting
# We'll use a single color to speed up rendering
ggplot(data = all_bats) +
  geom_sf(fill = "blue", alpha = 0.2, color = NA) + 
  theme_minimal() +
  labs(title = "Simplified Overlap Map")

head(all_bats)

##Adding study points
## 2. Load and Clean Study Points
diet_data <- read.csv("Diet_database.csv")

# Clean coordinates immediately to prevent crashes
# Make sure Location is included here!
study_points_cleaned <- diet_cleaned %>%
  select(Article_Title, Bat_spp, Location, Easting, Northing, Zone) %>%
  distinct() %>%
  filter(!is.na(Easting)) # Remove any rows with missing coordinates

## 3. Convert to Spatial Points (study_sf)
# We use all_bats directly now since you don't need to crop
study_sf <- lapply(unique(study_points_cleaned$Zone), function(z) {
  zone_subset <- filter(study_points_cleaned, Zone == z)
  
  st_as_sf(zone_subset, coords = c("Easting", "Northing"), crs = 32600 + z) %>%
    st_transform(4326) 
}) %>% bind_rows()

## 4. Plot to Verify
ggplot() +
  # Plot all bat ranges
  geom_sf(data = all_bats, aes(fill = species_name), alpha = 0.3) +
  new_scale_fill() + 
  # Plot study sites
  geom_sf(data = study_sf, aes(fill = Bat_spp), 
          color = "black", size = 2, shape = 21, stroke = 0.5) + 
  theme_minimal() +
  labs(title = "Study Sites and Bat Range Overlap", 
       subtitle = "Using full North American range maps")
##Spatial analysis 

# 1. Spatial Join
intersections <- st_intersects(study_sf, all_bats_cropped)

# 2. Add the count as a new column to your study data
study_sf$species_richness <- lengths(intersections)

# 3. Check the results
summary(study_sf$species_richness)

#plot
ggplot() +
  # 1. Bat Ranges
  geom_sf(data = all_bats_cropped, aes(fill = species_name), alpha = 0.2) +
  labs(fill = "Range Map Species") +
  
  new_scale_fill() +
  
  # 2. Study Sites
  # We now map 'size' to the count we just created
  geom_sf(data = study_sf, 
          aes(fill = Bat_spp, size = species_richness), 
          color = "black", 
          shape = 21, 
          stroke = 0.5,
          alpha = 0.8) +
  
  # 3. Adjust size scale so the dots are visible but distinct
  scale_size_continuous(range = c(1, 6), name = "Species Richness\n(Range Count)") +
  
  theme_minimal() +
  labs(
    title = "Study Sites vs. Local Species Richness",
    subtitle = "Larger dots indicate areas where more bat species ranges overlap",
    fill = "Species Studied"
  )


###########################Dietary Comparison data cleaning

# 1. Load your fixed data
diet_raw <- read.csv("Diet_database.csv", stringsAsFactors = FALSE, check.names = FALSE)

# 2. Filter out Redundant Family Rows
diet_cleaned <- diet_raw %>%
  # Remove ghost columns
  select(where(~ !all(is.na(.)))) %>%
  mutate(
    Order = tolower(trimws(Order)),
    Bat_spp = trimws(Bat_spp),
    Prop_Numeric = as.numeric(as.character(Proportion)),
    # Scale percentages to decimals
    Prop_Final = ifelse(Prop_Numeric > 1, Prop_Numeric / 100, Prop_Numeric)
  ) %>%
  # --- THE FIX: IGNORE FAMILY ROWS ---
  # Only keep rows where Family is empty to avoid double-counting the Order total
  filter(is.na(Family) | Family == "" | Family == " ") %>%
  filter(!is.na(Prop_Final), !is.na(Order), Order != "")

# 3. Aggregate by Species
diet_composition <- diet_cleaned %>%
  group_by(Bat_spp, Order) %>%
  summarise(Mean_Prop = mean(Prop_Final, na.rm = TRUE), .groups = 'drop') %>%
  # Final normalization check: force the mean to sum to 100% per species
  group_by(Bat_spp) %>%
  mutate(Mean_Prop = Mean_Prop / sum(Mean_Prop)) %>%
  ungroup()

# 4. Plot
ggplot(diet_composition, aes(x = reorder(Bat_spp, Mean_Prop, sum), y = Mean_Prop, fill = Order)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  theme_minimal() +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1.001)) +
  scale_fill_viridis_d(option = "turbo") +
  labs(title = "Average Dietary Composition",
       subtitle = "Aggregated using Order-level totals only",
       x = "Bat Species", y = "Proportion of Diet", fill = "Prey Order")


#####Merging species richness map data with diet data
library(dplyr)
library(sf)

# 1. Prepare the richness data so there is exactly ONE value per unique site/species
richness_collapsed <- study_sf %>% 
  st_drop_geometry() %>% 
  select(Article_Title, Bat_spp, Location, species_richness) %>%
  mutate(across(c(Article_Title, Bat_spp, Location), trimws)) %>%
  # This ensures we have one richness value per Location/Species
  group_by(Article_Title, Bat_spp, Location) %>%
  summarize(map_species_richness = first(species_richness), .groups = "drop")

# 2. Join to diet_cleaned
# Every matching row in diet_cleaned will now get that richness value
diet_cleaned <- diet_cleaned %>%
  mutate(across(c(Article_Title, Bat_spp, Location), trimws)) %>%
  left_join(richness_collapsed, by = c("Article_Title", "Bat_spp", "Location"))

####Adding weight dataset

library(dplyr)

# 1. Load weights
bat_weights <- read.csv("bat_weights.csv", stringsAsFactors = FALSE)

# 2. Add weight to your existing dataset
# This maps 'Bat_spp' from your data to 'species' in the CSV
diet_cleaned <- diet_family_richness %>%
  left_join(
    bat_weights %>% select(species, weight), 
    by = c("Bat_spp" = "species")
  )

# 3. Quick Check
if(any(is.na(diet_family_richness$weight))) {
  missing <- unique(diet_family_richness$Bat_spp[is.na(diet_family_richness$weight)])
  warning("Weights missing for: ", paste(missing, collapse = ", "))
}
###Dietary analysis

library(dplyr)
library(tidyr)
library(vegan)
library(tibble)

# 1. Prepare the Data (Strict Aggregation)
diet_matrix <- diet_cleaned %>%
  # Create the ID
  mutate(unique_id = paste(Article_Title, Location, Bat_spp, sep = "_")) %>%
  
  # Select ONLY the columns needed for the matrix
  select(unique_id, Order, Prop_Final) %>%
  
  # Group strictly by ID and Order to sum any duplicates
  group_by(unique_id, Order) %>%
  summarize(total_prop = sum(Prop_Final, na.rm = TRUE), .groups = "drop") %>%
  
  # Now Pivot - this is guaranteed to work because duplicates are gone
  pivot_wider(names_from = Order, values_from = total_prop, values_fill = 0) %>%
  
  # Move ID to row names
  column_to_rownames("unique_id")


metadata <- diet_cleaned %>%
  mutate(unique_id = paste(Article_Title, Location, Bat_spp, sep = "_")) %>%
  
  # Keep your predictors
  select(unique_id, Article_Title, Bat_spp, species_richness, Northing, weight) %>%
  
  # Remove duplicate rows so we have 1 row per unique_id
  distinct(unique_id, .keep_all = TRUE) %>%
  
  # Ensure the metadata rows are in the EXACT same order as the matrix rows
  filter(unique_id %in% rownames(diet_matrix)) %>%
  arrange(match(unique_id, rownames(diet_matrix))) %>%
  column_to_rownames("unique_id")

# FINAL CHECK: This must be TRUE before you proceed!
all(rownames(diet_matrix) == rownames(metadata))

# 1. Calculate Dissimilarity (Bray-Curtis)
diet_dist <- vegdist(diet_matrix, method = "bray")

# 2. Run PERMANOVA
# strata = Article_Title ensures we control for the study effect
permanova_result <- adonis2(diet_dist ~ Bat_spp + species_richness + Northing + weight, 
                            data = metadata, 
                            strata = metadata$Article_Title,
                            permutations = 999)

print(permanova_result)

# This will split the "Model" row into 3 separate rows
permanova_detailed <- adonis2(diet_dist ~ Bat_spp + species_richness + Northing + weight, 
                              data = metadata, 
                              strata = metadata$Article_Title,
                              permutations = 999,
                              by = "margin") # <--- This is the magic switch

print(permanova_detailed)


library(ggplot2)
library(vegan)

# 1. Run the NMDS (Non-metric Multidimensional Scaling)
# This squashes the complex diet data into 2D space
nmds_result <- metaMDS(diet_matrix, distance = "bray", k = 2, trymax = 100)

# 2. Extract coordinates for plotting
nmds_scores <- as.data.frame(scores(nmds_result, display = "sites"))
nmds_scores$Bat_spp <- metadata$Bat_spp # Add species names back

# 3. Calculate "Centroids" (The average diet center for each species)
species_centers <- aggregate(cbind(NMDS1, NMDS2) ~ Bat_spp, data = nmds_scores, FUN = mean)

# 4. Plot
ggplot(nmds_scores, aes(x = NMDS1, y = NMDS2, color = Bat_spp)) +
  # Draw the points (individual populations)
  geom_point(alpha = 0.6) +
  # Draw the Species Labels at the center of their cluster
  geom_label(data = species_centers, aes(label = Bat_spp), 
             color = "black", fill = "white", alpha = 0.8, size = 3) +
  theme_minimal() +
  labs(title = "Bat Diet Composition by Species (NMDS)",
       subtitle = "Closer points = More similar diets") +
  theme(legend.position = "none") # Hide legend since labels are on the plot


###Cleaned figure
library(ggplot2)
library(dplyr)
library(vegan)

# 1. Clean up the data labels
# Create the Genus column again if needed
nmds_scores$Genus <- sub("_.*", "", nmds_scores$Bat_spp)

# 2. Calculate Centers for the Genus labels
genus_centers <- aggregate(cbind(NMDS1, NMDS2) ~ Genus, data = nmds_scores, FUN = mean)

library(ggplot2)

ggplot(nmds_scores, aes(x = NMDS1, y = NMDS2, color = Genus, fill = Genus)) +
  
  # 1. The Blobs (Confidence Ellipses)
  # Draws a shape around the 'territory' of each Genus
  stat_ellipse(geom = "polygon", alpha = 0.2, level = 0.95) +
  
  # 2. The Points
  geom_point(size = 2, alpha = 0.6) +
  
  # 3. Styling
  theme_minimal() +
  labs(title = "Bat Dietary Niches by Genus",
       x = "NMDS Axis 1", 
       y = "NMDS Axis 2") +
  theme(legend.position = "right") # Legend provides the key


#niche overlap

# Function to calculate Pianka's Niche Overlap
pianka_index <- function(matrix) {
  # matrix: rows = species, columns = insect orders
  n <- nrow(matrix)
  overlap_mat <- matrix(NA, n, n)
  rownames(overlap_mat) <- rownames(matrix)
  colnames(overlap_mat) <- rownames(matrix)
  
  for(i in 1:n) {
    for(j in 1:n) {
      p1 <- matrix[i, ]
      p2 <- matrix[j, ]
      # Pianka's Formula
      overlap_mat[i, j] <- sum(p1 * p2) / sqrt(sum(p1^2) * sum(p2^2))
    }
  }
  return(overlap_mat)
}

# 1. Prepare the average diet per species
# (Using the diet_composition object from your previous successful code)
species_diet_matrix <- diet_composition %>%
  pivot_wider(names_from = Order, values_from = Mean_Prop, values_fill = 0) %>%
  tibble::column_to_rownames("Bat_spp") %>%
  as.matrix()

# 2. Run the calculation
overlap_results <- pianka_index(species_diet_matrix)

# 3. View the results (Pairwise overlap)
print(round(overlap_results, 2))

# Observed average overlap across all species
obs_mean_overlap <- mean(overlap_results[lower.tri(overlap_results)])

# Permutation Test (Null Model)
set.seed(123)
null_overlaps <- replicate(999, {
  # Shuffle the diet proportions within each species
  null_matrix <- t(apply(species_diet_matrix, 1, sample))
  res <- pianka_index(null_matrix)
  mean(res[lower.tri(res)])
})

# Calculate P-value
p_value <- sum(null_overlaps <= obs_mean_overlap) / (999 + 1)

cat("Observed Mean Overlap:", obs_mean_overlap, "\n")
cat("P-value for Niche Partitioning:", p_value, "\n")

library(reshape2)
library(ggplot2)

# 1. Create a copy of your overlap results
overlap_tri <- overlap_results

# 2. Set the upper triangle and the diagonal (self-comparison) to NA
# This removes the x vs y / y vs x redundancy
overlap_tri[upper.tri(overlap_tri, diag = TRUE)] <- NA

# 3. Melt the data, telling R to ignore the NA values
overlap_melted <- melt(overlap_tri, na.rm = TRUE)

# 4. Plot the cleaned Heatmap
ggplot(overlap_melted, aes(Var1, Var2, fill = value)) +
  geom_tile(color = "white") + # Adds a thin white border to tiles
  scale_fill_gradient2(low = "white", high = "red", mid = "orange", 
                       midpoint = 0.5, limit = c(0,1), space = "Lab", 
                       name="Pianka\nOverlap") +
  theme_minimal() + 
  theme(
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
    panel.grid.major = element_blank() # Cleans up the background
  ) +
  labs(
    title = "Pairwise Niche Overlap (Pianka's Index)",
    x = "", y = ""
  ) +
  coord_fixed() # Makes the tiles perfectly square

##Analaysis and cleanning at the family level

library(dplyr)
library(tidyr)
library(vegan)
library(ggplot2)
library(reshape2)

# RE-RUN AGGREGATION (The Fix: Added 'Northing' to group_by)
diet_fam_matrix_df <- diet_family_richness %>%
  # Group by ALL the variables we want to keep
  group_by(Article_Title, Bat_spp, Location, Northing, species_richness, Family, weight.y) %>%
  summarise(Prop = sum(Prop_Final), .groups = "drop") %>%
  
  # Renormalize to ensure 100% sum per sample
  group_by(Article_Title, Bat_spp, Location, Northing) %>%
  mutate(Prop = Prop / sum(Prop)) %>%
  ungroup() %>%
  
  # Pivot to Matrix format
  pivot_wider(names_from = Family, values_from = Prop, values_fill = 0)

# Create the Matrix and Metadata
# 1. The Community Matrix (Only the insect families)
fam_matrix <- as.matrix(diet_fam_matrix_df %>% 
                          select(-Article_Title, -Bat_spp, -Location, -Northing, -species_richness))

# 2. The Metadata (Now includes Northing!)
fam_metadata <- diet_fam_matrix_df %>% 
  select(Article_Title, Bat_spp, Location, Northing, species_richness, weight.y)

# Check to make sure Northing is there now
head(fam_metadata$Northing)

##Permanova

library(vegan)

# 2. Run the Expanded PERMANOVA
# Formula: Diet depends on Species + Richness + Northing
set.seed(123)
permanova_expanded <- adonis2(fam_matrix ~ Bat_spp + species_richness + Northing + weight.y, 
                              data = fam_metadata, 
                              method = "bray",
                              by = "margin") # 'margin' tests the unique effect of each

# 3. Print the Breakdown
print("--- PERMANOVA RESULTS: Species, Richness, and Latitude ---")
print(permanova_expanded)

# --- DIAGNOSTIC: CHECK FOR COLLINEARITY ---
# If Richness and Northing are perfectly correlated (e.g., Richness drops as you go North),
# the model might struggle to tell them apart. Let's check:
correlation <- cor.test(fam_metadata$species_richness, fam_metadata$Northing)
print(paste("Correlation between Richness and Northing:", round(correlation$estimate, 2)))

# Test Weight alone to see its influence without Species 'blocking' it
set.seed(123)
weight_test <- adonis2(fam_matrix ~ weight.y + species_richness + Northing, 
                       data = fam_metadata, 
                       method = "bray", 
                       by = "margin")
print(weight_test)

##Making weight categorical to try and parse apart idenity and size
library(dplyr)
library(vegan)

# 1. Create the categories using case_when
fam_metadata <- fam_metadata %>%
  mutate(size_class = case_when(
    weight.y < 8 ~ "Small",
    weight.y >= 8.1 & weight.y <= 15 ~ "Medium",
    weight.y > 15.1 ~ "Large"
  )) %>%
  # Convert to factor and set the order so the legend is logical
  mutate(size_class = factor(size_class, levels = c("Small", "Medium", "Large")))

# 2. Run the Combined Model
# 'by = "margin"' lets us see if Species explains anything NEW after Size is accounted for
set.seed(123)
permanova_combined <- adonis2(fam_matrix ~ size_class + Bat_spp + Northing + size_class, 
                              data = fam_metadata, 
                              method = "bray", 
                              by = "margin")

print(permanova_combined)

##Plot (might be refunendant)

library(vegan)
library(ggplot2)

# --- 1. Run NMDS and Prepare Data ---
set.seed(123)
# Running NMDS on the family-level matrix
nmds <- metaMDS(fam_matrix, distance = "bray", k = 2, trymax = 100, autotransform = FALSE)

# Extract the coordinates for the samples (sites)
nmds_sites <- as.data.frame(scores(nmds, display = "sites"))

# Merge metadata for plotting
nmds_sites$Bat_spp <- fam_metadata$Bat_spp
nmds_sites$Size_Class <- fam_metadata$size_class

# Define specific colors for size clouds
size_cloud_colors <- c("Small" = "#1b9e77", "Medium" = "#d95f02", "Large" = "#7570b3")

# --- 2. Create the Plot ---
ggplot(nmds_sites, aes(x = NMDS1, y = NMDS2)) +
  
  # Layer 1: The "Clouds" (95% Confidence Ellipses by Size Class)
  stat_ellipse(geom = "polygon", aes(fill = Size_Class), 
               alpha = 0.2, level = 0.95, show.legend = FALSE) +
  
  # Layer 2: The Points (Colored by Species, Shaped by Size)
  geom_point(aes(color = Bat_spp, shape = Size_Class), 
             size = 3, alpha = 0.8) +
  
  # --- 3. Styling ---
  theme_bw() +
  
  # Use Turbo for distinct species colors in the legend
  scale_color_viridis_d(option = "turbo", name = "Species") +
  
  # Use the custom colors for the filled clouds
  scale_fill_manual(values = size_cloud_colors) +
  
  # FIXED: Assigning simple integers to shapes (16=circle, 17=triangle, 15=square)
  scale_shape_manual(values = c("Small" = 16, "Medium" = 17, "Large" = 15), 
                     name = "Size Class") +
  
  # Labels and Theme
  labs(
    title = "Family-Level Dietary Niche Space",
    subtitle = paste("Stress:", round(nmds$stress, 3), "| Clouds: 95% CI for Size Classes"),
    x = "NMDS Axis 1",
    y = "NMDS Axis 2"
  ) +
  theme(
    panel.grid = element_blank(),
    legend.position = "right",
    legend.box = "vertical",
    legend.text = element_text(face = "italic") # Italicize species names
  )

table(fam_metadata$size_class)

# 1. Extract the scores for the Insect Families (the "Species" in vegan terms)
family_scores <- as.data.frame(scores(nmds, display = "species"))
family_scores$Family <- rownames(family_scores)

# 2. Find the "Top Drivers" for Axis 1 and Axis 2
# We look for the largest absolute values (the families furthest from the center)
top_nmds1 <- family_scores %>%
  arrange(desc(abs(NMDS1))) %>%
  head(10)

top_nmds2 <- family_scores %>%
  arrange(desc(abs(NMDS2))) %>%
  head(10)

# 3. Print the results to see the "meaning" of your axes
print("--- Top Drivers for NMDS Axis 1 ---")
print(top_nmds1[, c("Family", "NMDS1")])

print("--- Top Drivers for NMDS Axis 2 ---")
print(top_nmds2[, c("Family", "NMDS2")])


# --- Pianka's Overlap Test (Family Level) ---

library(tidyr)
library(dplyr)

# --- 1. Re-Define the Calculation Function ---
pianka_calc <- function(mat) {
  n <- nrow(mat)
  res <- matrix(0, n, n)
  rownames(res) <- rownames(mat); colnames(res) <- rownames(mat)
  for(i in 1:n) for(j in 1:n) {
    p1 <- mat[i,]; p2 <- mat[j,]
    res[i,j] <- sum(p1*p2) / sqrt(sum(p1^2)*sum(p2^2))
  }
  return(res)
}

# --- 2. Re-Create the Species Profile Matrix ---
# This averages the diet for each bat species
species_fam_profile <- diet_fam_matrix_df %>%
  pivot_longer(cols = -c(Article_Title:species_richness), names_to = "Family", values_to = "Prop") %>%
  group_by(Bat_spp, Family) %>%
  summarise(Mean_Prop = mean(Prop), .groups = "drop") %>%
  # Renormalize to ensure sum = 1
  group_by(Bat_spp) %>%
  mutate(Mean_Prop = Mean_Prop / sum(Mean_Prop)) %>%
  pivot_wider(names_from = Family, values_from = Mean_Prop, values_fill = 0) %>%
  tibble::column_to_rownames("Bat_spp") %>%
  as.matrix()

# --- 3. Re-Calculate the Missing Object 'fam_overlap' ---
fam_overlap <- pianka_calc(species_fam_profile)

# --- 4. Run the Permutation Test ---
obs_mean_overlap <- mean(fam_overlap[lower.tri(fam_overlap)])

set.seed(123)
null_overlaps <- replicate(999, {
  # Shuffle the diet proportions within each species row
  null_matrix <- t(apply(species_fam_profile, 1, sample))
  res <- pianka_calc(null_matrix)
  mean(res[lower.tri(res)])
})

p_value <- sum(null_overlaps <= obs_mean_overlap) / (999 + 1)

# --- 5. Print Results ---
print(paste("Observed Mean Family Overlap:", round(obs_mean_overlap, 3)))
print(paste("P-value:", p_value))

if(p_value > 0.95) {
  print("CONCLUSION: Significant Niche Aggregation (Diets are highly similar)")
} else if(p_value < 0.05) {
  print("CONCLUSION: Significant Niche Partitioning (Diets are distinct)")
} else {
  print("CONCLUSION: Random Overlap (No strong pattern)")
}

##Heat map
library(reshape2)
library(ggplot2)

# 1. Clean the Matrix for Plotting
# We create a copy so we don't mess up the original math object
fam_overlap_plot <- fam_overlap

# Set upper triangle and diagonal to NA (removes redundancy & self-comparison)
fam_overlap_plot[upper.tri(fam_overlap_plot, diag = TRUE)] <- NA

# Melt into long format for ggplot
melted_fam_plot <- melt(fam_overlap_plot, na.rm = TRUE)

# 2. Plot
ggplot(melted_fam_plot, aes(Var1, Var2, fill = value)) +
  geom_tile(color = "white") + # White borders make tiles distinct
  
  # Color Scale: 0 = White, 1 = Red
  scale_fill_gradient2(
    low = "white", 
    high = "#b2182b", # Deep Red for high overlap
    mid = "#fddbc7", # Light Salmon for moderate overlap
    midpoint = 0.5, 
    limit = c(0, 1), 
    name = "Pianka's\nOverlap"
  ) +
  
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 10),
    axis.text.y = element_text(size = 10),
    panel.grid = element_blank(), # Removes distracting background grid
    legend.position = "right"
  ) +
  
  # Force square tiles
  coord_fixed() +
  
  labs(
    title = "Pairwise Niche Overlap (Family Level)",
    subtitle = "Redder tiles indicate species pairs with highly similar diets",
    x = "", 
    y = ""
  )

# SIMPER Analysis (Similarity Percentage)
# This asks: "Which insect families contribute most to the difference between bat species?"

library(vegan)

# 1. Run SIMPER on your family-level matrix
# (Using the 'fam_matrix' from the PERMANOVA step)
simper_results <- simper(fam_matrix, group = fam_metadata$Bat_spp)

# 2. View the summary
# This produces a LOT of text (every pair compared), so let's look at the top contributors
summary(simper_results)

# 3. Visualize the "Top Drivers" of difference
# We can extract the summary to find the #1 insect family driving differences
# (Note: This helps you write sentences like "Differences were driven primarily by consumption of X")


###figure of top drivers of dietary differences
library(dplyr)
library(ggplot2)
library(tidyr)
library(viridis)

# 1. Prepare Data (STRICT VERSION)
top_diet_contributors <- diet_fam_matrix_df %>%
  pivot_longer(
    cols = -c(Article_Title, Bat_spp, Location, Northing, species_richness), 
    names_to = "Family", 
    values_to = "Proportion"
  ) %>%
  group_by(Bat_spp, Family) %>%
  summarise(Mean_Prop = mean(Proportion), .groups = "drop") %>%
  
  # Remove pure zeros before ranking
  filter(Mean_Prop > 0) %>%
  
  group_by(Bat_spp) %>%
  # STRICT: Grab exactly 5, even if there are ties
  slice_max(order_by = Mean_Prop, n = 5, with_ties = FALSE) %>%
  ungroup() %>%
  
  # Final Factor cleanup
  mutate(Family = as.character(Family)) # Convert to text to kill old factor levels

# 2. Define Palette based on the FINAL list
families_in_plot <- sort(unique(top_diet_contributors$Family))
my_colors <- viridis::viridis_pal(option = "turbo")(length(families_in_plot))
names(my_colors) <- families_in_plot

# 3. Plot
ggplot(top_diet_contributors, aes(x = Bat_spp, y = Mean_Prop, fill = Family)) +
  geom_bar(stat = "identity", position = "stack", width = 0.8, color = "white") +
  
  scale_fill_manual(
    values = my_colors, 
    name = "Insect Family",
    breaks = families_in_plot # <--- Strictly limits legend to this list
  ) +
  
  theme_minimal() +
  labs(
    title = "Top 5 Dietary Contributors by Bat Species",
    subtitle = "Mean proportion of diet (Family Level)",
    x = "", 
    y = "Proportion of Diet"
  ) +
  coord_flip() + 
  theme(
    axis.text.y = element_text(face = "bold", size = 10, color = "black"),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "right",
    legend.text = element_text(size = 9) # Adjust size if legend is too long
  )


######General summary stats about study characteristics 

library(dplyr)
library(ggplot2)

# 1. Prepare the Data
study_resolution_stats <- diet_data %>%
  mutate(Level = tolower(trimws(`Level.of.classification`))) %>%
  mutate(Level = case_when(
    Level %in% c("genus/species") ~ "species",
    Level == "famly" ~ "family",
    TRUE ~ Level
  )) %>%
  filter(Level %in% c("order", "family", "genus", "species")) %>%
  group_by(Level) %>%
  summarise(Study_Count = n_distinct(Article_Title)) %>%
  ungroup() %>%
  # Force Biological Order
  mutate(Level = factor(Level, levels = c("order", "family", "genus", "species")))

# 2. Plot (Clean bars, no numbers)
ggplot(study_resolution_stats, aes(x = Level, y = Study_Count, fill = Level)) +
  geom_bar(stat = "identity", color = "black", width = 0.7) +
  
  # Styling
  theme_minimal() +
  scale_fill_viridis_d(option = "mako", begin = 0.3, end = 0.8, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  
  labs(
    title = "Taxonomic Resolution of Diet Studies",
    x = "Taxonomic Level",
    y = "Number of Studies"
  ) +
  
  theme(
    axis.text.x = element_text(size = 12, face = "bold", vjust = 0.5),
    axis.title = element_text(size = 11, face = "bold"),
    panel.grid.major.x = element_blank()
  )

####sptial temporal and avaialbility considerations 

library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)

# 1. Prepare the Data
study_metadata_counts <- diet_data %>%
  # Select and clean column names
  select(
    Article_Title, 
    Insect_Ref = matches("insect.*reference|reference.*collection", ignore.case = TRUE),
    Temporal = matches("temporal.*consideration", ignore.case = TRUE),
    Spatial = matches("spatial.*consideration", ignore.case = TRUE)
  ) %>%
  distinct() %>% # Count each study only once
  pivot_longer(
    cols = c(Insect_Ref, Temporal, Spatial),
    names_to = "Category",
    values_to = "Response"
  ) %>%
  mutate(
    Response = tolower(trimws(Response)),
    Response = case_when(
      str_starts(Response, "y") ~ "Yes",
      str_starts(Response, "n") ~ "No",
      TRUE ~ "Not Reported"
    ),
    Category = recode(Category,
                      "Insect_Ref" = "Insect Reference\nCollection Used?",
                      "Temporal" = "Temporal Factors\nConsidered?",
                      "Spatial" = "Spatial Factors\nConsidered?"
    )
  ) %>%
  filter(Response %in% c("Yes", "No")) %>%
  group_by(Category, Response) %>%
  summarise(Count = n(), .groups = "drop")

# 2. Create the Grouped Bar Chart
ggplot(study_metadata_counts, aes(x = Category, y = Count, fill = Response)) +
  # Use position = dodge to place bars side-by-side
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7, color = "black") +
  
  # Styling
  theme_minimal() +
  scale_fill_manual(values = c("No" = "#d73027", "Yes" = "#1a9850")) + # High contrast Red/Green
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) + 
  
  labs(
    title = "Methodological Considerations in Diet Studies",
    subtitle = "Analysis of unique studies in the database",
    x = "",
    y = "Number of Studies",
    fill = "Response"
  ) +
  
  theme(
    axis.text.x = element_text(size = 11, face = "bold", color = "black"),
    axis.title.y = element_text(face = "bold"),
    legend.position = "top",
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank()
  )
