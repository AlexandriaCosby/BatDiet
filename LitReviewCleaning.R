# ==============================================================================
# Project: North American Bat Diet Synthesis
# Description: Spatial analysis, dietary niche partitioning, and literature review stats.
# ==============================================================================

# --- 1. SETUP & LIBRARIES ---
library(dplyr)
library(sf)
library(ggplot2)
library(ggnewscale)
library(tidyr)
library(vegan)
library(reshape2)
library(stringr)
library(viridis)
library(ape)

# Set global theme
theme_set(theme_minimal())

# --- 2. CUSTOM FUNCTIONS ---

# Pianka's Niche Overlap Index
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

# --- 3. DATA IMPORT ---

# Paths
shp_path <- "C:/Users/alexa/OneDrive - University of Guelph/Documents/WEB lab/North American Bat Polygons"
diet_file <- "Diet_database.csv"
weight_file <- "bat_weights.csv"

# Load Raw Data
diet_raw <- read.csv(diet_file, stringsAsFactors = FALSE, check.names = FALSE)
bat_weights <- read.csv(weight_file, stringsAsFactors = FALSE)
shp_files <- list.files(shp_path, pattern = "\\.shp$", full.names = TRUE)

# Load data
diet_raw <- read.csv("Diet_database.csv", stringsAsFactors = FALSE, check.names = FALSE)

# --- FIX: Remove Ghost Columns ---
# This keeps only columns where the name is NOT empty and NOT NA
diet_raw <- diet_raw[, !is.na(names(diet_raw)) & names(diet_raw) != ""]

# Check if it worked
print("Cleaned Column Names:")
print(names(diet_raw))

# ==============================================================================
# PART A: SPATIAL ANALYSIS (Bat Ranges & Richness)
# ==============================================================================

# 1. Process Range Maps (Simplify & Standardize)
all_bats_list <- lapply(shp_files, function(f) {
  st_read(f, quiet = TRUE) %>%
    st_simplify(preserveTopology = TRUE, dTolerance = 1000) %>% # Speed optimization
    st_transform(4326) %>%
    mutate(species_name = basename(f)) %>%
    select(species_name, geometry)
})
all_bats <- do.call(rbind, all_bats_list)

# 2. Process Study Sites
# Extract coordinates from raw diet data
study_points <- diet_raw %>%
  select(Article_Title, Bat_spp, Location, Easting, Northing, Zone) %>%
  distinct() %>%
  filter(!is.na(Easting))

# Convert to sf (handling different UTM zones)
study_sf <- lapply(unique(study_points$Zone), function(z) {
  study_points %>%
    filter(Zone == z) %>%
    st_as_sf(coords = c("Easting", "Northing"), crs = 32600 + z) %>%
    st_transform(4326)
}) %>% bind_rows()

# 3. Calculate Species Richness (Intersection)
# Overlap study sites with range maps
intersections <- st_intersects(study_sf, all_bats)
study_sf$species_richness <- lengths(intersections)

# 4. Create Richness Lookup Table (for merging later)
richness_lookup <- study_sf %>%
  st_drop_geometry() %>%
  mutate(across(c(Article_Title, Bat_spp, Location), trimws)) %>%
  group_by(Article_Title, Bat_spp, Location) %>%
  summarize(species_richness = first(species_richness), .groups = "drop")

###Fix island issue 
#Manual Fix for Haida Gwaii
# The spatial join missed this island site, so we force the correct richness (4).
study_sf <- study_sf %>%
  mutate(species_richness = case_when(
    # Look for "Haida" in the Location name (ignore case helps avoid typos)
    grepl("Haida", Location, ignore.case = TRUE) ~ 4,
    
    # Keep all other values the same
    TRUE ~ species_richness
  ))

# 5. Plot: Study Sites vs Bat Ranges (Trimmed to North America)
ggplot() +
  # 1. Bat Ranges (Polygons)
  geom_sf(data = all_bats, aes(fill = species_name), alpha = 0.2, show.legend = FALSE) +
  new_scale_fill() +
  
  # 2. Study Sites (Dots)
  geom_sf(data = study_sf, aes(fill = Bat_spp, size = species_richness), 
          shape = 21, stroke = 0.5, alpha = 0.8) +
  
  # 3. Scales & Labels
  scale_size_continuous(range = c(1, 6), name = "Range Overlap Count") +
  labs(title = "Study Sites & Local Species Richness", 
       subtitle = "Dots sized by number of overlapping bat ranges") +
  coord_sf(ylim = c(0, NA), expand = FALSE) +
  theme_minimal()
# 5.1. Generate a Heatmap Grid (Hexagons)

# 5. Generate High-Res Grid (Pixels)
richness_grid <- st_make_grid(all_bats, n = c(200, 200), square = TRUE) %>% 
  st_as_sf() %>% 
  st_make_valid() %>%
  rename(geometry = x)

# 6. Count Overlaps
grid_overlaps <- st_intersects(richness_grid, all_bats)
richness_grid$richness <- lengths(grid_overlaps)

# 7. Convert to Pixel Data (Raster-style)
pixel_data <- richness_grid %>%
  st_centroid() %>%
  mutate(
    X = st_coordinates(.)[,1],
    Y = st_coordinates(.)[,2]
  ) %>%
  st_drop_geometry() %>%
  filter(richness > 0, Y > 0) # Remove ocean & below equator

# 8. Plot: Continuous Surface + Visible Legend
ggplot() +
  # --- Layer 1: Continuous Surface ---
  geom_tile(data = pixel_data, aes(x = X, y = Y, fill = richness)) +
  scale_fill_viridis_c(option = "inferno", name = "Regional\nRichness") +
  
  # --- RESET SCALE for dots ---
  new_scale_fill() +
  
  # --- Layer 2: Study Sites ---
  geom_sf(data = study_sf, 
          aes(fill = Bat_spp, size = species_richness), 
          color = "white", # White borders for map visibility
          shape = 21, 
          stroke = 0.3, 
          alpha = 0.9) +
  
  # --- Color & Size Scales ---
  scale_fill_viridis_d(option = "viridis", name = "Study Species") +
  scale_size_continuous(range = c(2, 6), name = "Local Count") +
  
  # --- THE FIX: LEGEND OVERRIDES ---
  guides(
    # 1. Make Species dots big and colorful
    fill = guide_legend(override.aes = list(size = 5)), 
    
    # 2. Make "Local Count" dots visible (Grey fill instead of transparent/white)
    size = guide_legend(override.aes = list(fill = "grey50", color = "black"))
  ) + 
  
  # --- Theme & Labels ---
  labs(title = "North American Bat Species Richness",
       subtitle = "Background: Range Overlap Density | Dots: Study Sites") +
  theme_minimal() +
  theme(
    axis.text = element_blank(), 
    axis.title = element_blank(),
    panel.grid = element_blank()
  )
# ==============================================================================
# PART B: DIET DATA CLEANING (The Master Merge)
# ==============================================================================

# 1. Clean Proportions & Text
diet_cleaned <- diet_raw %>%
  # Remove completely empty columns/rows
  select(where(~ !all(is.na(.)))) %>% 
  mutate(
    across(c(Article_Title, Bat_spp, Location, Order, Family), ~trimws(.)),
    Order = tolower(Order),
    
    # Clean up non-numeric characters (like "<1")
    Proportion = gsub("<", "", Proportion),
    Prop_Numeric = as.numeric(as.character(Proportion)),
    
    # --- LOGIC FIX ---
    # If >= 1, treat as Percentage (divide by 100).
    # If < 1, treat as Proportion (keep as is).
    # Note: This treats "1" as 1% (0.01). If you have 100% diets entered as "1", this will need adjustment.
    Prop_Final = ifelse(Prop_Numeric >= 1, Prop_Numeric / 100, Prop_Numeric)
  ) %>%
  # Remove rows that didn't parse or are empty
  filter(!is.na(Prop_Final), !is.na(Order), Order != "")

# 2. Join Richness and Weights (The Master Merge)
diet_master <- diet_cleaned %>%
  left_join(richness_lookup, by = c("Article_Title", "Bat_spp", "Location")) %>%
  left_join(select(bat_weights, species, weight), by = c("Bat_spp" = "species")) %>%
  mutate(
    # Create Size Classes
    size_class = case_when(
      weight < 8 ~ "Small",
      weight >= 8.1 & weight <= 15 ~ "Medium",
      weight > 15.1 ~ "Large"
    ),
    # Set Factor Levels
    size_class = factor(size_class, levels = c("Small", "Medium", "Large"))
  )

# Quick Check: Did the logic work?
summary(diet_master$Prop_Final)
# ==============================================================================
# PART C: ORDER-LEVEL ANALYSIS
# ==============================================================================

# 1. Aggregation (Exclude Family-level rows to prevent double counting)
diet_order_comp <- diet_master %>%
  filter(is.na(Family) | Family == "") %>% 
  group_by(Bat_spp, Order) %>%
  summarise(Mean_Prop = mean(Prop_Final, na.rm = TRUE), .groups = 'drop') %>%
  group_by(Bat_spp) %>%
  mutate(Mean_Prop = Mean_Prop / sum(Mean_Prop)) # Renormalize to 100%

# 2. Plot: Order Composition
ggplot(diet_order_comp, aes(x = reorder(Bat_spp, Mean_Prop, sum), y = Mean_Prop, fill = Order)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1.001)) +
  scale_fill_viridis_d(option = "turbo") +
  labs(title = "Average Dietary Composition (Order Level)", x = "Bat Species", y = "Proportion")

# ==============================================================================
# PART G (Corrected): ORDER-LEVEL ANALYSIS (With Northing and easting)
# ==============================================================================

# 1. Prepare Data: Roll up ALL data to Order level
diet_order_matrix_df <- diet_master %>%
  # Ensure Order is present
  filter(!is.na(Order), Order != "") %>%
  
  # Safety Filter: Drop rows with missing predictors (Added Easting here)
  filter(!is.na(size_class), !is.na(species_richness), 
         !is.na(weight), !is.na(Northing), !is.na(Easting)) %>%
  
  # --- FIX: Added 'Easting' to the group_by so it is preserved ---
  group_by(Article_Title, Bat_spp, Location, species_richness, 
           Northing, Easting, Order, weight, size_class) %>%
  summarise(Prop = sum(Prop_Final), .groups = "drop") %>%
  
  # Renormalize to ensure every site sums to 1 (100%)
  group_by(Article_Title, Bat_spp, Location) %>%
  mutate(Prop = Prop / sum(Prop)) %>%
  ungroup() %>%
  
  # Pivot to Matrix format
  pivot_wider(names_from = Order, values_from = Prop, values_fill = 0)

# 2. Split into Matrix and Metadata
# Math Matrix (Only Insect Orders)
# --- FIX: Exclude 'Easting' from the math matrix ---
order_matrix <- as.matrix(diet_order_matrix_df %>% 
                            select(-c(Article_Title, Bat_spp, Location, species_richness, 
                                      Northing, Easting, weight, size_class)))

# Metadata (Predictors)
# --- FIX: Include 'Easting' in the metadata ---
order_metadata <- diet_order_matrix_df %>% 
  select(Article_Title, Bat_spp, species_richness, Northing, Easting, weight, size_class)

# 3. Run PERMANOVA (Order Level)
set.seed(123)
permanova_order_final <- adonis2(order_matrix ~ size_class + Bat_spp + species_richness + Northing + Easting, 
                                 data = order_metadata, 
                                 method = "bray", 
                                 # strata = order_metadata$Article_Title, # Removed for cross-study comparison
                                 by = "terms") # Sequential testing

print(permanova_order_final)

# 4. Run PCoA (To Visualize)
dist_order <- vegdist(order_matrix, method = "bray")
pcoa_order <- pcoa(dist_order)

# Extract Scores
pcoa_scores_ord <- as.data.frame(pcoa_order$vectors[, 1:2])
colnames(pcoa_scores_ord) <- c("PCoA1", "PCoA2")
pcoa_scores_ord <- cbind(pcoa_scores_ord, order_metadata)

# Calculate Variance Explained
var_expl_ord <- round(pcoa_order$values$Relative_eig[1:2] * 100, 1)

# 5. Plot
ggplot(pcoa_scores_ord, aes(x = PCoA1, y = PCoA2)) +
  stat_ellipse(geom = "polygon", aes(fill = size_class), alpha = 0.2, level = 0.95, show.legend = FALSE) +
  geom_point(aes(color = Bat_spp, shape = size_class), size = 3, alpha = 0.8) +
  scale_color_viridis_d(option = "turbo", name = "Species") +
  scale_fill_manual(values = c("Small" = "#1b9e77", "Medium" = "#d95f02", "Large" = "#7570b3")) +
  scale_shape_manual(values = c(16, 17, 15), name = "Size Class") +
  
  labs(title = "Order-Level Dietary Differences (PCoA)", 
       subtitle = paste0("Including Northing | PCoA 1: ", var_expl_ord[1], "%"),
       x = paste0("PCoA 1 (", var_expl_ord[1], "%)"), 
       y = paste0("PCoA 2 (", var_expl_ord[2], "%)")) +
  
  theme_bw() + 
  theme(panel.grid = element_blank())



# ==============================================================================
# CHECK: NUMBER OF STUDIES PER SPECIES
# ==============================================================================

# 1. Calculate Counts
species_study_counts <- diet_raw %>%
  # Clean names to ensure "Myotis lucifugus" and "Myotis lucifugus " don't count separately
  mutate(Bat_spp = trimws(Bat_spp)) %>%
  group_by(Bat_spp) %>%
  # Count unique Article Titles per species
  summarise(Study_Count = n_distinct(Article_Title)) %>%
  arrange(desc(Study_Count))

# 2. Print the Table
print("--- Number of Studies per Species ---")
print(species_study_counts, n = 50) # Print up to 50 rows

# 3. Plot the Results
ggplot(species_study_counts, aes(x = reorder(Bat_spp, Study_Count), y = Study_Count)) +
  geom_bar(stat = "identity", fill = "steelblue", color = "black", width = 0.7) +
  
  # Add the number labels on top of the bars
  geom_text(aes(label = Study_Count), hjust = -0.2, size = 3.5) +
  
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) + # Add space for labels
  coord_flip() + # Flip to make names readable
  
  labs(title = "Research Effort by Bat Species",
       x = "",
       y = "Number of Unique Studies") +
  
  theme_minimal() +
  theme(
    axis.text.y = element_text(face = "italic", size = 10), # Italicize species names
    panel.grid.major.y = element_blank()
  )


# ==============================================================================
# PART D: FAMILY-LEVEL ANALYSIS (Multivariate)
# ==============================================================================

# 1. Prepare Data Matrix (Strict Aggregation)
diet_fam_prep <- diet_master %>%
  filter(!is.na(Family), Family != "") %>%
  
  # Safety Filter: Drop rows with missing predictors (Added Easting)
  filter(!is.na(size_class), !is.na(species_richness), 
         !is.na(weight), !is.na(Northing), !is.na(Easting)) %>%
  
  # --- FIX: Added 'Easting' to group_by so it stays with the data ---
  group_by(Article_Title, Bat_spp, Location, species_richness, 
           Northing, Easting, Family, weight, size_class) %>%
  summarise(Prop = sum(Prop_Final), .groups = "drop") %>%
  
  # Renormalize per sample
  group_by(Article_Title, Bat_spp, Location) %>%
  mutate(Prop = Prop / sum(Prop)) %>%
  ungroup() 

# Pivot to Matrix
fam_matrix_df <- diet_fam_prep %>%
  pivot_wider(names_from = Family, values_from = Prop, values_fill = 0)

# Split into Matrix (Species) and Metadata (Env)
# --- FIX: Exclude 'Easting' from the math matrix ---
fam_matrix <- as.matrix(select(fam_matrix_df, -c(Article_Title, Bat_spp, Location, species_richness, 
                                                 Northing, Easting, weight, size_class)))

# --- FIX: Include 'Easting' in the metadata ---
fam_metadata <- select(fam_matrix_df, Article_Title, Bat_spp, species_richness, 
                       Northing, Easting, weight, size_class)

# 2. PERMANOVA (Cross-Study Analysis)
set.seed(123)
permanova_fam_final <- adonis2(fam_matrix ~ size_class + Bat_spp + species_richness + Northing + Easting, 
                               data = fam_metadata, 
                               method = "bray", 
                               # strata = fam_metadata$Article_Title, # Removed for cross-study comparison
                               by = "terms") # Sequential testing

print(permanova_fam_final)

# ==============================================================================
# ALTERNATIVE: PCoA (Principal Coordinates Analysis)
# ==============================================================================

library(vegan)
library(ape) # Needed for PCoA (pcoa function)
library(ggplot2)

# 1. Calculate the Distance Matrix (Bray-Curtis)
# We use the 'clean' matrix from the previous step
dist_matrix <- vegdist(fam_matrix_clean, method = "bray")

# 2. Run PCoA
pcoa_res <- pcoa(dist_matrix)

# 3. Extract Axis Scores (PC1 and PC2)
pcoa_scores <- as.data.frame(pcoa_res$vectors[, 1:2])
colnames(pcoa_scores) <- c("PCoA1", "PCoA2")

# 4. Bind Metadata
pcoa_scores <- cbind(pcoa_scores, fam_metadata_clean)

# 5. Calculate "Variance Explained" (To add to axes)
# This tells you how much of the diet difference is shown in the plot
var_explained <- round(pcoa_res$values$Relative_eig[1:2] * 100, 1)

# 6. Plot PCoA
ggplot(pcoa_scores, aes(x = PCoA1, y = PCoA2)) +
  # Ellipses
  stat_ellipse(geom = "polygon", aes(fill = size_class), alpha = 0.2, level = 0.95, show.legend = FALSE) +
  
  # Points
  geom_point(aes(color = Bat_spp, shape = size_class), size = 3, alpha = 0.8) +
  
  # Styling
  scale_color_viridis_d(option = "turbo", name = "Species") +
  scale_fill_manual(values = c("Small" = "#1b9e77", "Medium" = "#d95f02", "Large" = "#7570b3")) +
  scale_shape_manual(values = c(16, 17, 15), name = "Size Class") +
  
  # Labels with Variance Explained
  labs(title = "Family-Level Dietary Differences (PCoA)", 
       subtitle = "Using Bray-Curtis Dissimilarity",
       x = paste0("PCoA 1 (", var_explained[1], "%)"), 
       y = paste0("PCoA 2 (", var_explained[2], "%)")) +
  
  theme_bw() + 
  theme(panel.grid = element_blank())


# ==============================================================================
# PART E: NICHE OVERLAP (Pianka's Index)
# ==============================================================================

# 1. Create Average Diet Matrix
species_fam_profile <- diet_fam_prep %>%
  group_by(Bat_spp, Family) %>%
  summarise(Mean_Prop = mean(Prop), .groups = "drop") %>%
  group_by(Bat_spp) %>%
  mutate(Mean_Prop = Mean_Prop / sum(Mean_Prop)) %>%
  pivot_wider(names_from = Family, values_from = Mean_Prop, values_fill = 0) %>%
  tibble::column_to_rownames("Bat_spp") %>%
  as.matrix()

# 2. Calculate & Permute
fam_overlap <- pianka_calc(species_fam_profile)
obs_mean <- mean(fam_overlap[lower.tri(fam_overlap)])

# Null model
set.seed(123)
null_overlaps <- replicate(999, {
  mean(pianka_calc(t(apply(species_fam_profile, 1, sample)))[lower.tri(fam_overlap)])
})
p_val <- sum(null_overlaps <= obs_mean) / 1000

cat("Observed Mean Overlap:", obs_mean, "| P-value:", p_val, "\n")

# 3. Heatmap
overlap_melt <- melt(fam_overlap)
overlap_melt$value[upper.tri(fam_overlap, diag=TRUE)] <- NA # Clean upper triangle
overlap_melt <- na.omit(overlap_melt)

ggplot(overlap_melt, aes(Var1, Var2, fill = value)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "white", high = "#b2182b", mid = "#fddbc7", midpoint = 0.5, limit = c(0,1)) +
  coord_fixed() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Pairwise Niche Overlap (Pianka)", x = "", y = "", fill = "Overlap")

# ==============================================================================
# PART F: LITERATURE REVIEW STATS
# ==============================================================================

# 1. Taxonomic Resolution Stats
res_stats <- diet_raw %>%
  mutate(Level = tolower(trimws(`Classification`))) %>%
  mutate(Level = case_when(Level %in% c("genus/species") ~ "species", Level == "famly" ~ "family", TRUE ~ Level)) %>%
  filter(Level %in% c("order", "family", "genus", "species")) %>%
  group_by(Level) %>% summarise(Count = n_distinct(Article_Title)) %>%
  mutate(Level = factor(Level, levels = c("order", "family", "genus", "species")))

ggplot(res_stats, aes(x = Level, y = Count, fill = Level)) +
  geom_bar(stat = "identity", color = "black", width = 0.7) +
  scale_fill_viridis_d(option = "mako", begin = 0.3, end = 0.8) +
  labs(title = "Taxonomic Resolution of Diet Studies", y = "Number of Studies")

# 2. Methodological Considerations
meta_stats <- diet_raw %>%
  select(Article_Title, Insect_Ref = matches("insect.*reference"), 
         Temporal = matches("temporal"), Spatial = matches("spatial")) %>%
  distinct() %>%
  pivot_longer(cols = -Article_Title, names_to = "Category", values_to = "Response") %>%
  mutate(Response = case_when(str_starts(tolower(Response), "y") ~ "Yes", 
                              str_starts(tolower(Response), "n") ~ "No", 
                              TRUE ~ "Not Reported")) %>%
  filter(Response %in% c("Yes", "No"))

ggplot(meta_stats, aes(x = Category, fill = Response)) +
  geom_bar(position = "dodge", color = "black") +
  scale_fill_manual(values = c("No" = "#d73027", "Yes" = "#1a9850")) +
  labs(title = "Methodological Considerations", x = "")

# ==============================================================================
# PRINT EXACT NUMBERS FOR RESULTS
# ==============================================================================

# --- 1. TAXONOMIC RESOLUTION NUMBERS ---
res_stats_final <- res_stats %>%
  ungroup() %>%
  mutate(
    Total_Studies = sum(Count),
    Percentage = round(Count / Total_Studies * 100, 1)
  )

print("--- TABLE 1: TAXONOMIC RESOLUTION ---")
print(res_stats_final)


# --- 2. METHODOLOGICAL CONSIDERATIONS NUMBERS ---
# (The plot counted these automatically, so we must calculate the counts manually here)

method_stats_final <- meta_stats %>%
  group_by(Category, Response) %>%
  summarise(Count = n(), .groups = "drop") %>%
  group_by(Category) %>%
  mutate(
    Total = sum(Count),
    Percentage = round(Count / Total * 100, 1)
  )

print("--- TABLE 2: METHODOLOGICAL CONSIDERATIONS ---")
print(method_stats_final)
