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

####################### 1. Data Cleaning ##########################
# Phenotypic traits
trait_cols <- c("Yield", "DTM", "Lodging", "Height", "SeedWt", "Protein", "Oil")

# Drop missing values
nmf_matrix_clean <- new_df[complete.cases(new_df[, trait_cols]), ]

# Drop locations with fewer than 3 observations
loc_counts <- table(nmf_matrix_clean$Location)
keep_loc <- names(loc_counts[loc_counts >= 3])
nmf_matrix <- nmf_matrix_clean[nmf_matrix_clean$Location %in% keep_loc, ]
nmf_matrix$Variety.Name <- factor(nmf_matrix$Variety.Name)

cat("Rows after cleaning:", nrow(nmf_matrix), "\n")

# Convert Soil.Type to lowercase
nmf_matrix$Soil.Type <- tolower(nmf_matrix$Soil.Type)


####################### 2. ENVIRONMENTAL DATA ##########################

# Read in the lat/long data file
env_data <- read.csv("Environmental_Summary.csv")

# Append province for geocoding
nmf_matrix$Location <- paste(nmf_matrix$Location, "Ontario", sep = ", ")

# Merge original cleaned dataset with the environmental data obtained from NASA POWER
merged_df <- left_join(nmf_matrix, env_data, by = c("Location", "Year"))

# Check if all rows were merged
nrow(nmf_matrix)
nrow(merged_df)

# Check if there are any NA values introduced after merge
sum(is.na(merged_df))


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

# Environmental matrix
X <- as.matrix(merged_df[, env_cols])

# Min-max normalization
X_scaled <- apply(X, 2, function(col) (col - min(col))/(max(col)-min(col)))

# Optimal k (UIK method)
max_k_allowed <- min(dim(X_scaled)) - 1
k_list <- 2:max_k_allowed
rss_list <- numeric(length(k_list))

for (i in seq_along(k_list)) {
  set.seed(2025)
  fit <- nmf(X_scaled, rank = k_list[i], nrun = 10, .options = "N")
  rss_list[i] <- residuals(fit)
}

best_k <- uik(k_list, rss_list)[1]
cat("Optimal k:", best_k, "\n")

# Final NMF
set.seed(2025)
nmf_final <- nmf(X_scaled, rank = best_k, nrun = 50, .options = "N")

W <- basis(nmf_final)
colnames(W) <- paste0("Factor", seq_len(ncol(W)))
H <- coef(nmf_final)
rownames(H) <- paste0("Factor", 1:best_k)
colnames(H) <- env_cols

# Save UIK plot
png("../Project_Trial_6 - full_data_env_and_pheno/UIK_FullData.png", width=800, height=600)
plot(k_list, rss_list, type="b", main="UIK - Full Data", xlab="k", ylab="RSS")
abline(v=best_k, col="red", lty=2)
dev.off()

# Save Heatmap
png("../Project_Trial_6 - full_data_env_and_pheno/Heatmap_FullData.png", width=800, height=600)
heatmap(H, scale="column", margins=c(15,8), main="NMF Factors - Full Data")
dev.off()
  

####################### 5. LM ##########################

merged_df <- cbind(merged_df, W)

set.seed(2025)
train_index <- createDataPartition(merged_df$Yield, p=0.7, list=FALSE)
df_train <- merged_df[train_index, ]
df_test <- merged_df[-train_index, ]

colnames(W)

lm_formula <- as.formula(paste("Yield ~ Soil.Type +", paste(colnames(W), collapse=" + ")))

lm_fit <- lm(lm_formula, data=df_train)
preds <- predict(lm_fit, df_test)

rmse_nmf <- RMSE(preds, df_test$Yield)
r2_nmf <- cor(preds, df_test$Yield)^2

cat("NMF Pipeline - RMSE:", rmse_nmf, "- R²:", r2_nmf, "\n")

png("../Project_Trial_6 - full_data_env_and_pheno/LM_NMF_FullData.png", width=800, height=600)
plot(df_test$Yield, preds, xlab="Actual", ylab="Predicted", main="Predicted vs Actual (NMF+LM)", col="blue")
abline(0,1,col="red",lty=2)
dev.off()

####################### 6. LM WITHOUT NMF ##########################

predictors <- c("Soil.Type", env_cols)

df_noNMF <- merged_df[, c("Yield", predictors)]
df_noNMF <- na.omit(df_noNMF)

set.seed(2025)
train_index <- createDataPartition(df_noNMF$Yield, p=0.7, list=FALSE)
df_train <- df_noNMF[train_index, ]
df_test <- df_noNMF[-train_index, ]

lm_formula_noNMF <- as.formula(paste("Yield ~", paste(predictors, collapse=" + ")))

lm_fit_noNMF <- lm(lm_formula_noNMF, data=df_train)
preds_noNMF <- predict(lm_fit_noNMF, df_test)

rmse_noNMF <- RMSE(preds_noNMF, df_test$Yield)
r2_noNMF <- cor(preds_noNMF, df_test$Yield)^2

cat("No-NMF Pipeline - RMSE:", rmse_noNMF, "- R²:", r2_noNMF, "\n")

png("../Project_Trial_6 - full_data_env_and_pheno/LM_NoNMF_FullData.png", width=800, height=600)
plot(df_test$Yield, preds_noNMF, xlab="Actual", ylab="Predicted", main="Predicted vs Actual (No NMF)", col="blue")
abline(0,1,col="red",lty=2)
dev.off()

