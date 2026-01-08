# ==============================================================================
# Project: North American Bat Diet Synthesis - Multi-Level Analysis
# Description: Integrated Spatial, Order-level, and Family-level analysis
# Author: Alexandria Cosby
# ==============================================================================

# --- 1. SETUP & LIBRARIES ---
library(dplyr)
library(sf)
library(ggplot2)
library(ggnewscale)
library(tidyr)
library(vegan)
library(viridis)

# --- 2. SPATIAL ANALYSIS: RANGE OVERLAP & RICHNESS ---

# 2.1 Load and Simplify Bat Ranges
folder_path <- "C:/Users/alexa/OneDrive - University of Guelph/Documents/WEB lab/North American Bat Polygons"
shp_files   <- list.files(folder_path, pattern = "\\.shp$", full.names = TRUE)

all_bats <- do.call(rbind, lapply(shp_files, function(f) {
  st_read(f, quiet = TRUE) %>%
    st_simplify(preserveTopology = TRUE, dTolerance = 1000) %>% 
    st_transform(4326) %>%
    mutate(species_name = basename(f)) %>%
    select(species_name, geometry)
}))

# 2.2 Load Diet and Clean "Ghost Columns"
diet_raw <- read.csv("Diet_database.csv", stringsAsFactors = FALSE, check.names = FALSE)
diet_raw <- diet_raw[, names(diet_raw) != "" & !is.na(names(diet_raw))]

# 2.3 Convert Study Sites to Spatial
# 'remove = FALSE' is critical to keep Northing/Easting columns for later use
study_points <- diet_raw %>%
  select(Article_Title, Bat_spp, Location, Easting, Northing, Zone) %>%
  mutate(across(c(Article_Title, Bat_spp, Location), trimws)) %>%
  distinct() %>% 
  filter(!is.na(Easting), !is.na(Northing), !is.na(Zone))

study_sf <- bind_rows(lapply(unique(study_points$Zone), function(z) {
  st_as_sf(filter(study_points, Zone == z), coords = c("Easting", "Northing"), 
           crs = 32600 + z, remove = FALSE) %>%
    st_transform(4326)
}))

# 2.4 Calculate Richness and Create Metadata Lookup
study_sf$species_richness <- lengths(st_intersects(study_sf, all_bats))

richness_lookup <- study_sf %>% 
  st_drop_geometry() %>%
  group_by(Article_Title, Bat_spp, Location) %>%
  summarize(species_richness = first(species_richness), 
            Northing = first(Northing), 
            .groups = "drop")

# --- 3. DATA HARMONIZATION: WEIGHTS & PROPORTIONS ---

bat_weights <- read.csv("bat_weights.csv", stringsAsFactors = FALSE)

diet_base <- diet_raw %>%
  mutate(across(c(Article_Title, Bat_spp, Location), ~trimws(as.character(.)))) %>%
  # Handle special characters and convert to 0-1 scale
  mutate(Proportion = gsub("<", "", Proportion),
         Prop_Numeric = as.numeric(as.character(Proportion)),
         Prop_Numeric = ifelse(is.na(Prop_Numeric), 0, Prop_Numeric),
         Prop_Final = ifelse(Prop_Numeric > 1, Prop_Numeric / 100, Prop_Numeric)) %>%
  filter(!is.na(Prop_Final)) %>%
  # Join Metadata (Spatial + Weight)
  left_join(richness_lookup, by = c("Article_Title", "Bat_spp", "Location")) %>%
  left_join(bat_weights %>% rename(Bat_spp = species) %>% mutate(Bat_spp = trimws(Bat_spp)), 
            by = "Bat_spp") %>%
  # Functional Size Classes
  mutate(size_class = factor(case_when(weight < 10 ~ "Small", weight <= 20 ~ "Medium", weight > 20 ~ "Large"), 
                             levels = c("Small", "Medium", "Large")))

# --- 4. PERMANOVA: ORDER LEVEL ANALYSIS ---

# 1. Prepare Matrix
diet_order <- diet_base %>%
  filter(is.na(Family) | Family == "" | Family == " ") %>%
  group_by(Article_Title, Location, Bat_spp, Order) %>%
  summarise(Prop = sum(Prop_Final), .groups = "drop") %>%
  pivot_wider(names_from = Order, values_from = Prop, values_fill = 0)

# --- RE-ESTABLISH RICHNESS_CLEAN ---
# This pulls the spatial data back into a joinable format
richness_clean <- study_sf %>% 
  st_drop_geometry() %>% # Convert from map data to a standard table
  select(Article_Title, Bat_spp, Location, species_richness, Northing) %>%
  group_by(Article_Title, Bat_spp, Location) %>%
  summarize(
    species_richness = first(species_richness), 
    Northing = first(Northing), 
    .groups = "drop"
  )

# Now your metadata build will work:
metadata_order <- diet_order %>% 
  select(Article_Title, Location, Bat_spp) %>%
  left_join(richness_clean, by = c("Article_Title", "Location", "Bat_spp")) %>%
  left_join(bat_weights %>% rename(Bat_spp = species) %>% mutate(Bat_spp = trimws(Bat_spp)), 
            by = "Bat_spp") %>%
  mutate(size_class = factor(case_when(weight < 10 ~ "Small", weight <= 20 ~ "Medium", weight > 20 ~ "Large"), 
                             levels = c("Small", "Medium", "Large"))) %>%
  filter(!is.na(weight), !is.na(Northing))

# 2. Build Metadata strictly including Northing
metadata_order <- diet_order %>% 
  select(Article_Title, Location, Bat_spp) %>%
  # Force join Northing and Richness directly from our clean spatial source
  left_join(richness_clean, by = c("Article_Title", "Location", "Bat_spp")) %>%
  # Force join Weights
  left_join(bat_weights %>% rename(Bat_spp = species) %>% mutate(Bat_spp = trimws(Bat_spp)), 
            by = "Bat_spp") %>%
  # Re-calculate size class to be safe
  mutate(size_class = factor(case_when(weight < 10 ~ "Small", weight <= 20 ~ "Medium", weight > 20 ~ "Large"), 
                             levels = c("Small", "Medium", "Large"))) %>%
  # Now filter for complete cases
  filter(!is.na(weight), !is.na(Northing))

# 3. Align Matrix
order_mat <- diet_order %>% 
  filter(paste(Article_Title, Location, Bat_spp) %in% 
           paste(metadata_order$Article_Title, metadata_order$Location, metadata_order$Bat_spp)) %>%
  select(-Article_Title, -Location, -Bat_spp)

print("--- ORDER LEVEL PERMANOVA ---")
print(adonis2(vegdist(order_mat) ~ size_class + Bat_spp + Northing, 
              data = metadata_order, by = "margin", strata = metadata_order$Article_Title))

# --- 5. PERMANOVA: FAMILY LEVEL ANALYSIS ---

# 1. Prepare Matrix
diet_fam <- diet_base %>%
  filter(!is.na(Family), Family != "", Family != " ") %>%
  group_by(Article_Title, Location, Bat_spp, Family) %>%
  summarise(Prop = sum(Prop_Final), .groups = "drop") %>%
  pivot_wider(names_from = Family, values_from = Prop, values_fill = 0)

# 2. Build Metadata strictly including Northing
metadata_fam <- diet_fam %>% 
  select(Article_Title, Location, Bat_spp) %>%
  left_join(richness_clean, by = c("Article_Title", "Location", "Bat_spp")) %>%
  left_join(bat_weights %>% rename(Bat_spp = species) %>% mutate(Bat_spp = trimws(Bat_spp)), 
            by = "Bat_spp") %>%
  mutate(size_class = factor(case_when(weight < 10 ~ "Small", weight <= 20 ~ "Medium", weight > 20 ~ "Large"), 
                             levels = c("Small", "Medium", "Large"))) %>%
  filter(!is.na(weight), !is.na(Northing))

# 3. Align Matrix
fam_mat <- diet_fam %>% 
  filter(paste(Article_Title, Location, Bat_spp) %in% 
           paste(metadata_fam$Article_Title, metadata_fam$Location, metadata_fam$Bat_spp)) %>%
  select(-Article_Title, -Location, -Bat_spp)

print("--- FAMILY LEVEL PERMANOVA ---")
print(adonis2(vegdist(fam_mat) ~ size_class + Bat_spp + Northing, 
              data = metadata_fam, by = "margin", strata = metadata_fam$Article_Title))


