################## 0. PACKAGES AND LOAD DATA############################

install.packages("NMF")
install.packages("caret")
install.packages("inflection")
install.packages("nnls")
install.packages("tidygeocoder")
install.packages("nasapower")
install.packages("urltools", dependencies = TRUE)
install.packages("elevatr")
install.packages("sf")
install.packages("patchwork")
install.packages("cowplot")

library(gridExtra)
library(grid)
library(cowplot)
library(patchwork)
library(reshape2)
library(ggplot2)
library(sf)
library(nnls)
library(NMF)        # for nmf()
library(caret)      # for train/test split + RMSE
library(inflection) # for uik()
library(dplyr)
library(tidygeocoder)
library(nasapower)
library(elevatr)

# Load main OSACC dataset
new_df <- read.delim("./Data/OSACC Data Registered Varieties 1996-2021 with variety names.txt", header = T)

# Preview data
colnames(new_df)

####################### 1. Data Cleaning ##########################
# Store phenotypic traits
trait_cols <- c("Yield", "DTM", "Lodging", "Height", "SeedWt", "Protein", "Oil")

# Drop any rows with missing values
complete_cases <- complete.cases(new_df[,trait_cols])
nmf_matrix_clean <- new_df[complete_cases, ]

# Drop locations with < 3 valid rows
loc_counts <- table(nmf_matrix_clean$Location)
keep_loc <- names(loc_counts[loc_counts >= 3])
nmf_matrix <- nmf_matrix_clean[nmf_matrix_clean$Location %in% keep_loc, ]
nmf_matrix$Variety.Name <- factor(nmf_matrix$Variety.Name)

cat("Remaining rows afer cleaning:", nrow(nmf_matrix), "\n")


####################### 2. ENVIRONMENTAL DATA ##########################

# include Ontario as the province in the location column 
nmf_matrix$Location <- paste(nmf_matrix$Location, "Ontario", sep = ", ")

# Get location year pairs
loc_year_df <- nmf_matrix |> 
  select(Location, Year) |> 
  distinct()

# Get the unique location names from the dataset
unique_locations <- data.frame(Location = unique(nmf_matrix$Location))

# use tinygeocoder package to get the latitude and longitude of the locations
geocoded_locations <- unique_locations |> 
  geocode(Location, method = "osm", lat = Latitude, long = Longitude)

head(geocoded_locations)

# write the lat/long data on a file
write.csv(geocoded_locations, "location_coordinates.csv", row.names = FALSE)

# Read in the lat/long data file
loc_coords <- read.csv("location_coordinates.csv")

# Get the weather data from NASA POWER
weather_list <- list()
 
for (i in 1:nrow(loc_year_df)) {
  loc <- loc_year_df$Location[i]
  year <- loc_year_df$Year[i]
  
  coords <- loc_coords[loc_coords$Location == loc, ]
  lat <- coords$Latitude
  lon <- coords$Longitude
  
  # The dates of growing season (?)
  start_date <- paste0(year, "-05-01")
  end_date <- paste0(year, "-09-30")
  
  
  success <- FALSE
  attempts <- 0
  max_attempts <- 5
  
  while (!success && attempts < max_attempts) {
    attempts <- attempts + 1
    cat("Attempt", attempts, "for:", loc, year, "\n")
    
    weather_df <- tryCatch({
      get_power(
        community = "AG",
        lonlat = c(lon, lat),
        pars = c("T2M", "PRECTOTCORR", "ALLSKY_SFC_SW_DWN", "RH2M"),
        dates = c(start_date, end_date),
        temporal_api = "DAILY"
      )
    }, error = function(e) {
      message(paste("Error:", e$message))
      NULL
    })
    
    if (!is.null(weather_df)) {
      # Summarise
      summarized <- weather_df |> 
        summarise(
          Temp_c_avg = mean(T2M, na.rm = TRUE),
          Precip_mm_Total = sum(PRECTOTCORR, na.rm = TRUE),
          solarRad_avg = mean(ALLSKY_SFC_SW_DWN, na.rm = TRUE),
          RH2M_avg = mean(RH2M, na.rm = TRUE)
        ) |> 
        mutate(Location = loc, Year = year, Latitude = lat, Longitude = lon)
      
      weather_list[[i]] <- summarized
      success <- TRUE
      cat("Success for:", loc, year, "\n")
    } else {
      cat("Failed attempt", attempts, "for:", loc, year, "- Retrying...\n")
      Sys.sleep(5) # wait 5 seconds before retrying
    }
  }
  
  if (!success) {
    weather_list[[i]] <- NULL
    cat("Giving up on:", loc, year, "\n")
  }
  
  Sys.sleep(2) # short delay between locations
}


list_null <- which(sapply(weather_list, is.null))
list_null


saveRDS(weather_list, file = "weather_list_backup.Rds")

# Combine all weather data summaries
env_data <- bind_rows(weather_list)


# Inspect environmental data
head(env_data)

# Add elevation data
coord_df <- env_data |> 
  distinct(Latitude, Longitude)

# convert coord_df to simple feature type so elevatr package can use it
coord_sf <- st_as_sf(coord_df, coords = c("Longitude", "Latitude"), crs = 4326)

elevations <- get_elev_point(
  locations = coord_sf
)


coord_df$Elevation_m <- elevations$elevation

# Add elevation data to master list of environment data
env_data <- env_data |> 
  left_join(coord_df, by = c("Latitude", "Longitude"))

# save env data to file 
write.csv(env_data, "Environmental_Summary.csv", row.names = FALSE)

env_data <- read.csv("Environmental_Summary.csv")

# Make sure both data frames have matching types
str(nmf_matrix[, c("Location", "Year")])
str(env_data[, c("Location", "Year")])

# Merge original cleaned dataset with the environmental data obtained from NASA POWER
merged_df <- left_join(nmf_matrix, env_data, by = c("Location", "Year"))

# Check if all rows were merged
nrow(nmf_matrix)
nrow(merged_df)

# Check if there are any NA values introduced after merge
sum(is.na(merged_df))



####################### 3. SPLITTING DATASET ##########################

merged_df <- merged_df |> 
  filter(HUZ != 2200)

unique(merged_df$HUZ)

# Create a new column that only has 2 test group names
merged_df <- merged_df |> 
  mutate(TestGroup = ifelse(Test == "RR", "RR", "Conventional"))

# split the dataset by HUZ and further by Test type
test_split_list <- merged_df |> 
  group_by(TestGroup) |> 
  group_split()

# Name each of the subsets
names(test_split_list) <- sapply(test_split_list, function(df) {
  paste0(df$TestGroup[1])
})

# Check the names
names(test_split_list)


# Lower case and trim white spaces of soil types 
for (i in seq_along(test_split_list)) {
  test_split_list[[i]]$Soil.Type <- tolower(trimws(test_split_list[[i]]$Soil.Type))
}

####################### 4. NMF ##########################

env_cols <- c(
  "Temp_c_avg",
  "Precip_mm_Total",
  "solarRad_avg",
  "RH2M_avg",
  "Elevation_m",
  "Latitude",
  "Longitude",
  "Height",
  "Lodging",
  "SeedWt"
)

# NMF settings
n_runs <- 10
final_nmf_runs <- 50

# store results
nmf_results <- list()

# Loop through each subset
for (name in names(test_split_list)) {
  
  df_subset <- test_split_list[[name]]
  
  # Get the environment columns of each subset
  X <- as.matrix(df_subset[, env_cols])
  
  # Min-max normalize (safely handle constant or NA columns)
  X_scaled <- apply(X, 2, function(col) {
    col_min <- min(col, na.rm = TRUE)
    col_max <- max(col, na.rm = TRUE)
    
    # If column is constant, return zeros 
    if ((col_max - col_min) == 0) {
      return(rep(0, length(col)))
    } else {
      return((col - col_min) / (col_max - col_min))
    }
  })
  
  # Drop constant or all-NA columns
  drop_cols <- which(apply(X_scaled, 2, function(col) all(is.na(col)) || sd(col, na.rm = TRUE) == 0))
  
  if (length(drop_cols) > 0) {
    cat("Dropping problematic columns:", colnames(X_scaled)[drop_cols], "\n")
    X_scaled <- X_scaled[, -drop_cols, drop = FALSE]
  }
  
  # Check dimensions and set k_list accordingly
  max_k_allowed <- min(dim(X_scaled)) - 1
  k_list <- 2:max_k_allowed
  
  # find optimal K: run nmf [n_runs] number of times on each [k_list] value
  rss_list <- numeric(length(k_list))
  
  for (i in seq_along(k_list)) {
    k <- k_list[i]
    
    set.seed(2025)
    fit <- nmf(
      X_scaled,
      rank = k,
      nrun = n_runs,
      .options = "N"
    )
    # Store the residuals of each nmf run
    rss_list[i] <- residuals(fit)
  }
  
  # Store the best k using UIK method
  best_k <- uik(k_list, rss_list)[1]
  
  # Create folder to save UIK Plots
  dir.create("UIK_Plots", showWarnings = FALSE)
  
  # File name for each plot
  filename <- paste0("UIK_Plots/UIK_", name, ".png")
  
  png(filename, width = 800, height = 600)
  
  
  # Plot UIK curve
  plot(k_list, rss_list, type = "b", 
       main = paste0("UIK - ", name),
       xlab = "Number of Factors (k)", 
       ylab = "Residual Sum of Squares")
  abline(v = best_k, col = "red", lty = 2)
  
  dev.off()
  
  cat("\nSubset:", name, " → Best k =", best_k, "\n")
  
  # Final NMF fit using best k
  set.seed(2025)
  nmf_final <- nmf(
    X_scaled,
    rank = best_k,
    nrun = final_nmf_runs,
    .options = "N"
  )
  
  # Store the W and H matrices
  W <- basis(nmf_final)
  colnames(W) <- paste0("Factor", seq_len(ncol(W)))
  
  H <- coef(nmf_final)
  rownames(H) <- paste0("Factor", seq_len(nrow(H)))
  colnames(H) <- env_cols
  
  
  # Store the results
  nmf_results[[name]] <- list(
    best_k = best_k,
    W = W,
    H = H,
    nmf_model = nmf_final
  )
  
}

####################### 4.5 INTERPRET NMF FACTORS ##########################

# Create folder to save Heatmaps
dir.create("../Project_Trial_4 - split_by_test_type_env_and_pheno/H_Heatmaps", showWarnings = FALSE)

heatmap_list <- list()

for (group_name in names(nmf_results)) {
  cat("\nPlotting H matrix for group:", group_name, "\n")
  
  H <- nmf_results[[group_name]]$H
  
  # Convert H to a tidy dataframe for ggplot
  H_df <- as.data.frame(H)
  H_df$Factor <- rownames(H_df)
  
  H_long <- melt(H_df, id.vars = "Factor",
                 variable.name = "Variable",
                 value.name = "Loading")
  
  # Create ggplot heatmap
  p <- ggplot(H_long, aes(x = Variable, y = Factor, fill = Loading)) +
    geom_tile(color = "white") +
    scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0,
                         limits = c(min(H_long$Loading, na.rm = TRUE),
                                    max(H_long$Loading, na.rm = TRUE))) +
    theme_minimal() +
    labs(
      title = group_name,
      x = "Environmental Variable",
      y = "Latent Factor"
    ) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8))
  
  heatmap_list[[group_name]] <- p
  
  cat("Stored heatmap for group:", group_name, "\n")
  
}

# Combine all heatmaps into one figure
combined_heatmap <- wrap_plots(heatmap_list, ncol = 2) +
  plot_annotation(title = "NMF H Matrices Across Subsets")

# Save combined heatmap
ggsave("../Project_Trial_4 - split_by_test_type_env_and_pheno/H_Heatmaps/Combined_H_Heatmaps.png",
       combined_heatmap,
       width = 15, height = 10, dpi = 300)

####################### 5. LM ##########################

dir.create("../Project_Trial_4 - split_by_test_type_env_and_pheno/LM_Plots", showWarnings = FALSE)

plot_list <- list()
lm_results <- list()

for (subset in names(nmf_results)) {
  cat("\nFitting LM for subset:", subset, "\n")
  
  W <- nmf_results[[subset]]$W
  
  df_subset <- test_split_list[[subset]]
  
  W_df <- as.data.frame(W)
  W_df$RowID <- seq_len(nrow(W))
  
  df_subset <- df_subset |> 
    mutate(RowID = seq_len(nrow(df_subset)))
  
  df_lm <- df_subset |> 
    left_join(W_df, by = "RowID")
  
  if (!is.factor(df_lm$Soil.Type)) {
    df_lm$Soil.Type <- as.factor(df_lm$Soil.Type)
  }
  
  if (nrow(df_lm) < 10) {
    cat("Subset too small for train/test split: ", subset, "\n")
    next
  }
  
  # Create train and test splits (70/30)
  set.seed(2025)
  train_index <- createDataPartition(df_lm$Yield, p = 0.7, list = FALSE)
  df_train <- df_lm[train_index, ]
  df_test <- df_lm[-train_index, ]
  
  # LM equation
  factor_cols <- colnames(W)[!colnames(W) %in% "RowID"]
  formula_str <- paste("Yield ~ Soil.Type +", paste(factor_cols, collapse = " + "))
  lm_formula <- as.formula(formula_str)
  
  # Fit the model on training set
  lm_fit <- lm(lm_formula, data = df_train)
  
  # Predict on test set
  preds <- predict(lm_fit, newdata = df_test)
  
  # Calculate RMSE and R-squared
  rmse_val <- RMSE(preds, df_test$Yield)
  r2_val <- cor(preds, df_test$Yield)^2
  
  cat("Subset:", subset, "- RMSE = ", rmse_val, "- R-squared = ", r2_val, "\n")
  
  # plot predict vs actual for test set
  plot_df <- data.frame(
    actual = df_test$Yield,
    predicted = preds
  )
  
  p <- ggplot(plot_df, aes(x = actual, y = predicted)) +
    geom_point(color = "steelblue", alpha = 0.7) +
    geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
    theme_minimal() +
    labs(
      title = subset,
      x = NULL,
      y = NULL
    ) +
    theme(
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.text.x = element_blank(),
      axis.text.y = element_blank()
    )
  
  plot_list[[subset]] <- p
  
  # save results
  lm_results[[subset]] <- list(
    lm_model = lm_fit,
    rmse = rmse_val,
    r_squared = r2_val,
    model_summary = summary(lm_fit),
    plot_file = filename
  )
}

# Arrange plots into a grid
plots_grob <- arrangeGrob(
  grobs = plot_list,
  ncol = 2
)

# Now add global axis labels
final_plot <- arrangeGrob(
  plots_grob,
  left = textGrob("Predicted Yield", rot = 90, gp = gpar(fontsize = 16)),
  bottom = textGrob("Actual Yield", gp = gpar(fontsize = 16)),
  top = textGrob("Predicted vs Actual Yield Across Subsets (NMF+LM)", gp = gpar(fontsize = 18, fontface = "bold"))
)

# Save to file
ggsave("../Project_Trial_4 - split_by_test_type_env_and_pheno/LM_Plots/Combined_Pred_vs_Actual_with_labels.png",
       plot = final_plot,
       width = 15, height = 10, dpi = 300)

####################### 6. LM WITHOUT NMF ##########################

dir.create("../Project_Trial_4 - split_by_test_type_env_and_pheno/LM_Plots_NoNMF", showWarnings = FALSE)

# list to store results
plot_list_noNMF <- list()
lm_results_noNMF <- list()


for (subset in names(test_split_list)) {
  cat("\nRunning LM without NMF for subset", subset, "\n")
  
  df_subset <- test_split_list[[subset]]
  
  
  set.seed(2025)
  train_index <- createDataPartition(df_subset$Yield, p = 0.7, list = FALSE)
  df_train <- df_subset[train_index, ]
  df_test <- df_subset[-train_index, ]
  
  # LM equation
  formula_str <- paste("Yield ~ Soil.Type +", paste(env_cols, collapse = " + "))
  lm_formula <- as.formula(formula_str)
  
  # Fit the model on training set
  lm_fit <- lm(lm_formula, data = df_train)
  
  # Predict on test set
  preds <- predict(lm_fit, newdata = df_test)
  
  # Calculate RMSE and R-squared
  rmse_val <- RMSE(preds, df_test$Yield)
  r2_val <- cor(preds, df_test$Yield)^2
  
  cat("Subset:", subset, "- RMSE = ", rmse_val, "- R-squared = ", r2_val, "\n")
  
  # plot predict vs actual for test set
  plot_df <- data.frame(
    actual = df_test$Yield,
    predicted = preds
  )
  
  p <- ggplot(plot_df, aes(x = actual, y = predicted)) +
    geom_point(color = "steelblue", alpha = 0.7) +
    geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
    theme_minimal() +
    labs(
      title = subset,
      x = NULL,
      y = NULL
    ) +
    theme(
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.text.x = element_blank(),
      axis.text.y = element_blank()
    )
  
  plot_list_noNMF[[subset]] <- p
  
  cat("Saved no-NMF Plot for subset:", subset, "\n")
  
  # save results
  lm_results_noNMF[[subset]] <- list(
    lm_model = lm_fit,
    rmse = rmse_val,
    r_squared = r2_val,
    model_summary = summary(lm_fit),
    plot_file = filename
  )
}

# Arrange plots into a grid
plots_grob <- arrangeGrob(
  grobs = plot_list_noNMF,
  ncol = 2
)

# Now add global axis labels
final_plot <- arrangeGrob(
  plots_grob,
  left = textGrob("Predicted Yield", rot = 90, gp = gpar(fontsize = 16)),
  bottom = textGrob("Actual Yield", gp = gpar(fontsize = 16)),
  top = textGrob("Predicted vs Actual Yield Across Subsets (No NMF)", gp = gpar(fontsize = 18, fontface = "bold"))
)

# Save to file
ggsave("../Project_Trial_4 - split_by_test_type_env_and_pheno/LM_Plots_NoNMF/Combined_Pred_vs_Actual_NoNMF_with_labels.png",
       plot = final_plot,
       width = 15, height = 10, dpi = 300)
