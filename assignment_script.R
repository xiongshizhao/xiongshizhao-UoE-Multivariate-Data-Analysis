# Load required libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(factoextra)
library(cluster)
library(MASS)
library(RColorBrewer)

# Load Data
load("~/Desktop/Mathematics Project/xiongshizhao-UoE-Multivariate-Data-Analysis/Assignment3_Data.RData")


# --- Data Preparation ---

print("Original data dimensions (Genes x Samples):")
print(dim(data3T3))

# Transpose to Samples x Genes format
data <- as.data.frame(t(data3T3))
print("Transposed data dimensions (Samples x Genes):")
print(dim(data))


# --- Metadata Setup and Validation ---

# Retrieve sample ids from original data columns
if (is.null(colnames(data3T3))) {
  warning("data3T3 has no column names. Generating default SampleIDs.")
  sample_ids <- paste0("Sample_", 1:ncol(data3T3))
} else {
  sample_ids <- colnames(data3T3)
}

# Validate metadata lengths match sample count
expected_samples <- ncol(data3T3)
if (any(lengths(list(time = time.names, treatment = treatment.names, replicate = replicate.names)) != expected_samples)) {
  stop("One or more metadata vectors (time, treatment, replicate) do not match the number of samples.")
}

# Create metadata frame
data_meta <- data.frame(
  SampleID = sample_ids,
  Time = factor(time.names),
  Treatment = factor(treatment.names),
  Replicate = factor(replicate.names),
  stringsAsFactors = FALSE
)


# --- Combine Metadata and Transposed Data ---

# Verify rownames of transposed data match SampleIDs in metadata
if (is.null(rownames(data))) { 
  stop("Transposed expression data ('data') is missing rownames. Cannot verify alignment for cbind.")
}

# Use `identical()` to check values and order
if (!identical(data_meta$SampleID, rownames(data))) { # Compare meta SampleID to data rownames
  if (setequal(data_meta$SampleID, rownames(data))) { # Check if IDs match but order differs
    stop("Metadata SampleIDs and expression data rownames match but are in a different order. Use merge() for safe combination.")
  } else {
    stop("Mismatch between metadata SampleIDs and expression data rownames. Cannot safely use cbind(). Check input data consistency.")
  }
}

# Combine Data using `cbind`
message("SampleID alignment verified. Proceeding with cbind.")
data_full <- cbind(data_meta, data)

# Use default numeric row names for the final combined data frame
rownames(data_full) <- NULL


# --- Final Setup & Checks ---

# Display dimensions as a final check
print("Dimensions of final combined data frame (Samples x (Metadata + Genes)):")
print(dim(data_full))

# Define metadata and gene columns
meta_cols <- c("SampleID", "Time", "Treatment", "Replicate")
# Ensure meta_cols actually exist in data_full before using `setdiff`
if (!all(meta_cols %in% names(data_full))) {
  stop("One or more specified meta_cols are not present in the combined data_full dataframe.")
}
gene_columns <- setdiff(names(data_full), meta_cols)

print(paste("Identified", length(meta_cols), "metadata and", length(gene_columns), "gene columns.")) # Confirmation



# --- Task 1: Data Exploration ---
## Histogram of Log Gene Expression for a Specific Sample

# Example: First sample
sample_index <- 1 
if (sample_index > nrow(data_full) || sample_index < 1) {
  stop("Selected sample_index is out of bounds.")
}
sample_name <- data_full$SampleID[sample_index]

# Get data for the selected sample
sample_numeric <- as.numeric(data_full[sample_index, gene_columns])

# Calculate mean and median
mean_expression <- mean(sample_numeric, na.rm = TRUE)
median_expression <- median(sample_numeric, na.rm = TRUE)

# Generate the histogram
plot1 <- ggplot(data.frame(Expression = sample_numeric), aes(x = Expression)) +
  geom_histogram(aes(y = after_stat(density)), binwidth = 0.25,  
                 fill = "#69b3a2", color = "#e9ecef", alpha = 0.6) +  
  geom_density(color = "#ff7f0e", linewidth = 1.2) +
  geom_vline(xintercept = mean_expression,  
             color = "#1f77b4", linetype = "dashed", linewidth = 1) +
  geom_vline(xintercept = median_expression,  
             color = "darkred", linetype = "dotted", linewidth = 1) +
  
  # Add annotation for the MEAN value          
  annotate("text",  
           x = mean_expression,  
           y = 0,
           label = paste("Mean =", round(mean_expression, 2)),  
           hjust = -0.05,  
           vjust = -0.4, 
           color = "#1f77b4",
           size = 4.5,
           fontface = "bold") +
  
  # Add annotation for the MEDIAN value
  annotate("text",  
           x = median_expression,  
           y = 0,
           label = paste("Median =", round(median_expression, 2)),  
           hjust = -0.17,
           vjust = -2, 
           color = "darkred",
           size = 4.5,
           fontface = "bold") +
  
  ggtitle(sprintf("Distribution of Log Gene Expression for Sample '%s'", sample_name)) +
  xlab("Log Gene Expression") +
  ylab("Density") +  
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold", margin = margin(b = 30)),
    axis.title.x = element_text(size = 12, face = "bold", margin = margin(t = 15, r = 0, b = 0, l = 0)),
    axis.title.y = element_text(size = 12, face = "bold", margin = margin(t = 0, r = 15, b = 0, l = 0)),
    axis.text = element_text(size = 10, face = "bold"),  
    axis.line = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.7),
    plot.margin = unit(c(2, 2, 2, 2), "cm")
  ) 

# Print and save the boxplot
print(plot1)
ggsave("~/Desktop/plot1.png", plot = plot1, width = 12, height = 8, dpi = 900)


## Boxplot of Log Gene Expression for a Specific Gene

# Example: First Gene
gene_index <- 1 
if (gene_index > length(gene_columns) || gene_index < 1) {
  stop("Selected gene_index is out of bounds.")
}
gene_name <- gene_columns[gene_index] 

# Prepare data frame for the specific gene
gene_df_boxplot <- data.frame(
  Expression = data_full[[gene_name]],
  Time = data_full$Time,
  Treatment = data_full$Treatment,
  Condition = sprintf("%s (%sh)", data_full$Treatment, data_full$Time)
)

# Define the desired order for the x-axis conditions
condition_order <- c("Control (0h)", "Control (2h)", "Control (12h)",
                     "Insulin (2h)", "Insulin (12h)",
                     "Metformin (2h)", "Metformin (12h)")
gene_df_boxplot$Condition <- factor(gene_df_boxplot$Condition, levels = condition_order)

# Generate the boxplot
plot2 <- ggplot(gene_df_boxplot, aes(x = Condition, y = Expression, fill = Treatment)) +
  geom_boxplot(alpha = 0.7) +
  ggtitle(sprintf("Expression of Gene '%s' Across Conditions (Time & Treatment)", gene_name)) +
  xlab("Condition (Time.Treatment)") +
  ylab("Log Gene Expression") +  
  theme_classic() + 
  theme(
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold", margin = margin(b = 30)),
    axis.title.x = element_text(size = 12, face = "bold", margin = margin(t = 15, r = 0, b = 0, l = 0)),
    axis.title.y = element_text(size = 12, face = "bold", margin = margin(t = 0, r = 15, b = 0, l = 0)),
    axis.text = element_text(size = 10, face = "bold"),
    axis.text.x = element_text(angle = 20, hjust = 1, vjust = 1),
    axis.line = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.7),
    plot.margin = unit(c(2, 2, 2, 2), "cm"),
    legend.position = "top",
    legend.title = element_text(size=10, face="bold")
  ) +
  scale_fill_manual(values = c("Control" = "#69b3a2", "Insulin" = "#ff7f0e", "Metformin" = "#1f77b4"))

# Print and save the plot
print(plot2)
ggsave("~/Desktop/plot2.png", plot = plot2, width = 12, height = 8, dpi = 900)


## Time Course Plot for a Specific Gene

# Redefine data frame for the time course plot, requiring numeric time axis
gene_df_timecourse <- data.frame(
  Expression = data_full[[gene_name]],
  TimeNumeric = as.numeric(as.character(data_full$Time)),
  Treatment = data_full$Treatment
)

# Generate the time course plot
plot3 <- ggplot(gene_df_timecourse, aes(x = TimeNumeric, y = Expression, color = Treatment)) +
  geom_line(linewidth = 0.8) + 
  geom_point(size = 1.5) +
  
  # Add faceting, separate panels for each Treatment, stacked vertically (ncol=1)
  facet_wrap(~ Treatment, ncol = 1) +
  scale_x_continuous(breaks = c(0, 2, 12), labels = c("0h", "2h", "12h")) + 
  
  ggtitle(sprintf("Expression Trend of Gene '%s' Over Time (by Treatment)", gene_name)) +
  xlab("Time Point") +
  ylab("Log Gene Expression") +
  theme_classic() + 
  theme(
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold", margin = margin(b = 30)), 
    axis.title.x = element_text(size = 12, face = "bold", margin = margin(t = 15, r = 0, b = 0, l = 0)),
    axis.title.y = element_text(size = 12, face = "bold", margin = margin(t = 0, r = 15, b = 0, l = 0)),
    axis.text = element_text(size = 10, face = "bold"),
    axis.line = element_line(), 
    plot.margin = unit(c(2, 2, 2, 2), "cm"), 
    strip.background = element_rect(fill="grey85", color=NA),
    strip.text = element_text(face="bold", size=11),
    legend.position = "none" 
  ) +
  scale_color_manual(values = c("Control" = "#69b3a2", "Insulin" = "#ff7f0e", "Metformin" = "#1f77b4")) 

# Print and save the plot
print(plot3)
ggsave("~/Desktop/plot3.png", plot = plot3, width = 12, height = 8, dpi = 900)


## Outlier Sample Check

# Reshape data to long format for sample-wise boxplots
if (!exists("data_long")) {
  message("Reshaping data to long format.")
  data_long <- data_full %>%
    dplyr::select(all_of(meta_cols), all_of(gene_columns)) %>% 
    pivot_longer(
      cols = all_of(gene_columns),  
      names_to = "Gene",  
      values_to = "Expression"
    )
  message("Data reshaping complete.")
} else {
  message("Using existing 'data_long'.")
}

# Generate the boxplot
plot4 <- ggplot(data_long, aes(x = SampleID, y = Expression, fill = Treatment)) +
  geom_boxplot(alpha = 0.6) +
  ggtitle("Boxplot of Log Gene Expression Across Different Samples and Treatment Conditions") +
  xlab("Sample") +
  ylab("Log Gene Expression") +
  theme_classic() + 
  theme(
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold", margin = margin(b = 30)), 
    axis.title.x = element_text(size = 12, face = "bold", margin = margin(t = 15, r = 0, b = 0, l = 0)),
    axis.title.y = element_text(size = 12, face = "bold", margin = margin(t = 0, r = 15, b = 0, l = 0)),
    axis.text = element_text(size = 9, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5),
    axis.line = element_blank(), 
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.7),
    plot.margin = unit(c(1.5, 1.5, 1.5, 1.5), "cm"),
    legend.position = "top",
    legend.title = element_text(size=10, face="bold")
  ) +
  scale_fill_manual(values = c("Control" = "#69b3a2", "Insulin" = "#ff7f0e", "Metformin" = "#1f77b4"))

# Print and save the plot
print(plot4)
ggsave("~/Desktop/plot4.png", plot = plot4, width = 14, height = 8, dpi = 900) 



# --- Task 2: Principal Component Analysis ---

# Extract numeric gene data
pca_data <- data[, gene_columns] 
# Perform scaled PCA
pca_results <- prcomp(pca_data, scale. = TRUE, center = TRUE)

# Extract variance explained data for ALL PCs
num_pcs_total <- length(pca_results$sdev)
variance_data <- data.frame(
  PC = factor(paste0("PC", 1:num_pcs_total), levels = paste0("PC", 1:num_pcs_total)),
  VariancePercent = (pca_results$sdev^2 / sum(pca_results$sdev^2)) * 100
)

# Visualise the percentage of variance explained by each principal component
plot5 <- ggplot(variance_data, aes(x = PC, y = VariancePercent)) +
  geom_bar(stat = "identity", fill = "#69b3a2", color = "#e9ecef") +
  geom_line(aes(group = 1), color = "#ff7f0e", linewidth = 1) +
  geom_point(color = "#ff7f0e", size = 2) +
  geom_text(aes(label = paste0(round(VariancePercent, 1), "%")),
            vjust = -0.5,
            color = "black",
            size = 3) + 
  ggtitle("Scree Plot: Variance Explained by All Principal Components") +
  xlab("Principal Component") +
  ylab("% of Explained Variance") +
  ylim(0, max(variance_data$VariancePercent) * 1.1) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold", margin = margin(b = 30)),
    axis.title.x = element_text(size = 12, face = "bold", margin = margin(t = 15, r = 0, b = 0, l = 0)),
    axis.title.y = element_text(size = 12, face = "bold", margin = margin(t = 0, r = 15, b = 0, l = 0)),
    axis.text.x = element_text(size = 10, angle = 20, hjust = 0.5, vjust = 0.5),
    axis.text.y = element_text(size = 10),
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.7),
    plot.margin = unit(c(2, 2, 2, 2), "cm")
  )

# Print and save the plot
print(plot5)
ggsave("~/Desktop/plot5.png", plot = plot5, width = 12, height = 7, dpi = 900)


## Investigate how samples cluster on the first two principal components

# Create a data frame combining PC scores and sample metadata for plotting
pca_plot_df <- data.frame(
  PC1 = pca_results$x[, 1],
  PC2 = pca_results$x[, 2],
  SampleID = data_meta$SampleID,
  Time = data_meta$Time,
  Treatment = data_meta$Treatment,
  Replicate = data_meta$Replicate,
  Condition = factor(sprintf("%s (%sh)", data_meta$Treatment, data_meta$Time), levels = condition_order) # Use previously defined order
)

# Generate PCA plot: PC1 vs PC2 (Color points by Treatment, shape points by Time)
plot6 <- ggplot(pca_plot_df, aes(x = PC1, y = PC2, color = Treatment, shape = Time)) +
  geom_point(size = 3.5, alpha = 0.7) +
  scale_color_manual(values = c("Control" = "#69b3a2", "Insulin" = "#ff7f0e", "Metformin" = "#1f77b4")) +
  scale_shape_manual(values = c("0" = 16, "2" = 17, "12" = 15)) +
  ggtitle("PCA of Samples: PC1 vs PC2 (Colored by Treatment, Shaped by Time)") +
  xlab(paste0("PC1 (", round(summary(pca_results)$importance[2, 1] * 100, 1), "%)")) +
  ylab(paste0("PC2 (", round(summary(pca_results)$importance[2, 2] * 100, 1), "%)")) +
  theme_classic() + # Start with theme_classic as a base
  theme(
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold", margin = margin(b = 30)),
    axis.title.x = element_text(size = 12, face = "bold", margin = margin(t = 15, r = 0, b = 0, l = 0)),
    axis.title.y = element_text(size = 12, face = "bold", margin = margin(t = 0, r = 15, b = 0, l = 0)),
    axis.text = element_text(size = 10, face = "bold"),
    plot.margin = unit(c(2, 2, 2, 2), "cm"),
    axis.line = element_blank(), 
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.7),
    legend.title = element_text(size=12, face="bold")
  ) +
  labs(color = "Treatment", shape = "Time (h)")

# print and save the plot
print(plot6)
ggsave("~/Desktop/plot6.png", plot = plot6, width = 10, height = 8, dpi = 900)


## Select a small number of PCs that describe variability in the data

# Keep 2 PCs based on the previous scree plot
num_pcs_to_keep <- 2
# Print message indicating the choice
message(paste("Selecting the first", num_pcs_to_keep, "PCs for downstream analysis."))

# Extract the scores for the selected number of PCs
pc_scores_selected <- pca_results$x[, 1:num_pcs_to_keep]

# Calculate and report the cumulative variance explained by the selected PCs
variance_explained <- summary(pca_results)$importance[3, num_pcs_to_keep] * 100
message(paste("The selected", num_pcs_to_keep, "PCs explain", round(variance_explained, 2), "% of the total variance."))



# --- Task 3: Cluster Analysis ---

# Use the selected PC scores from Task 2
data_for_clustering <- pc_scores_selected

## Determine the Optimal Number of Clusters

# Elbow Method (Within-cluster Sum of Squares)
plot7 <- fviz_nbclust(data_for_clustering, kmeans, method = "wss", k.max = 10) +
  ggtitle("Estimating Optimal k for k-means using the Elbow Method") +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold", margin = margin(b = 30)),
    axis.title.x = element_text(size = 12, face = "bold", margin = margin(t = 15, r = 0, b = 0, l = 0)),
    axis.title.y = element_text(size = 12, face = "bold", margin = margin(t = 0, r = 15, b = 0, l = 0)),
    axis.text = element_text(size = 10, face = "bold"),
    axis.line = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.7),
    plot.margin = unit(c(2, 2, 2, 2), "cm")
  )

# Print and save the plot
print(plot7)
ggsave("~/Desktop/plot7.png", plot = plot7, width = 8, height = 6, dpi = 900)

# Silhouette Method (Average Silhouette Width)
plot8 <- fviz_nbclust(data_for_clustering, kmeans, method = "silhouette", k.max = 10) +
  ggtitle("Estimating Optimal k for k-means using the Silhouette Method") +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold", margin = margin(b = 30)),
    axis.title.x = element_text(size = 12, face = "bold", margin = margin(t = 15, r = 0, b = 0, l = 0)),
    axis.title.y = element_text(size = 12, face = "bold", margin = margin(t = 0, r = 15, b = 0, l = 0)),
    axis.text = element_text(size = 10, face = "bold"),
    axis.line = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.7),
    plot.margin = unit(c(2, 2, 2, 2), "cm")
  ) +
  
  annotate(
    geom = "text", x = 2.2, y = 0.48,
    label = "Optimal k = 2\n(Max Silhouette)",
    hjust = 0, vjust = 1,
    color = "darkred", size = 3.5, fontface = "bold"
  )

print(plot8)
ggsave("~/Desktop/plot8.png", plot = plot8, width = 10, height = 8, dpi = 900)

# Based on the plots, choose an optimal k
optimal_k <- 3
message(paste("Selected optimal number of clusters: k =", optimal_k))

## Perform and Compare Two Clustering Methods

set.seed(123) # for reproducibility

# Method 1: k-means Clustering
kmeans_clusters <- kmeans(data_for_clustering, centers = optimal_k, nstart = 25)
# Extract cluster assignments for each sample
clusters_kmeans <- kmeans_clusters$cluster

# Method 2: Hierarchical Clustering
# Calculate distance matrix (Euclidean distance on PC scores is appropriate)
dist_matrix <- dist(data_for_clustering, method = "euclidean")
hclust_results <- hclust(dist_matrix, method = "ward.D2")

# Cut the dendrogram to get the desired number of clusters
clusters_hclust <- cutree(hclust_results, k = optimal_k)

# Compare k-means assignments vs hierarchical assignments
cluster_comparison_table <- table(KMeans = clusters_kmeans, Hierarchical = clusters_hclust)
message("Contingency Table: k-means vs Hierarchical Clustering Assignments")
print(cluster_comparison_table)

# Visualise Hierarchical Clustering
plot9 <- fviz_dend(
  hclust_results,
  k = optimal_k,
  cex = 0.7,
  k_colors = c("#69b3a2", "#ff7f0e", "#1f77b4")[1:optimal_k],
  main = paste("Hierarchical Clustering Dendrogram (k=", optimal_k, ")"),
  xlab = "Height (Ward.D2 Distance)",
  horiz = TRUE,
  rect = TRUE,
  rect_fill = TRUE,
  ggtheme = theme_classic()
) +
  scale_x_continuous(breaks = c(0, 100, 200, 300, 400),
                     expand = expansion(mult=c(0, 0.05))) + 
  theme(
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold", margin = margin(b = 30)),
    axis.title.x = element_text(size = 12, face = "bold", margin = margin(t = 15, r = 0, b = 0, l = 0)),
    axis.title.y = element_blank(),
    axis.text.x = element_text(size = 10, face = "bold"),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.x = element_line(color = "black"),
    axis.line.y = element_blank(),
    panel.border = element_blank(), 
    panel.background = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.x = element_line(color = "grey70", linetype = "dashed", linewidth = 0.5),
    panel.grid.minor.x = element_blank(),
    plot.margin = unit(c(2, 2, 2, 2), "cm")
  )

# Print and save the plot
print(plot9)
ggsave("~/Desktop/plot9.png", plot = plot9, width = 10, height = 8, dpi = 900)


## Compare Clustering Results to Experimental Covariates

# Add cluster assignments to the PCA plot data frame
pca_plot_df$Cluster_KMeans <- factor(clusters_kmeans)
pca_plot_df$Cluster_Hclust <- factor(clusters_hclust)

# Plot PCA colored by k-means clusters
plot10 <- ggplot(pca_plot_df, aes(x = PC1, y = PC2, color = Cluster_KMeans, shape = Condition)) +
  geom_point(size = 3.5, alpha = 0.7) +
  scale_shape_manual(values = c("Control (0h)" = 15, "Control (2h)" = 16, "Control (12h)" = 17,
                                "Insulin (2h)" = 8, "Insulin (12h)" = 9,
                                "Metformin (2h)" = 10, "Metformin (12h)" = 18 )) +
  ggtitle(paste("PCA Plot Colored by k-means Clusters (k=", optimal_k, ")")) +
  xlab(paste0("PC1 (", round(summary(pca_results)$importance[2, 1] * 100, 1), "%)")) +
  ylab(paste0("PC2 (", round(summary(pca_results)$importance[2, 2] * 100, 1), "%)")) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold", margin = margin(b = 30)),
    axis.title.x = element_text(size = 12, face = "bold", margin = margin(t = 15, r = 0, b = 0, l = 0)),
    axis.title.y = element_text(size = 12, face = "bold", margin = margin(t = 0, r = 15, b = 0, l = 0)),
    axis.text = element_text(size = 10, face = "bold"),
    plot.margin = unit(c(2, 2, 2, 2), "cm"),
    axis.line = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.7),
    legend.title = element_text(size=12, face="bold"),
    legend.position = "right"
  ) +
  labs(color = paste("k-means (k=", optimal_k, ")", sep=""), shape = "Condition")

# Print and save the plot
print(plot10)
ggsave("~/Desktop/plot10.png", plot = plot10, width = 10, height = 8, dpi = 900)

# Plot PCA colored by Hierarchical clusters
plot11 <- ggplot(pca_plot_df, aes(x = PC1, y = PC2, color = Cluster_Hclust, shape = Condition)) +
  geom_point(size = 3.5, alpha = 0.7) +
  scale_shape_manual(values = c("Control (0h)" = 15, "Control (2h)" = 16, "Control (12h)" = 17,
                                "Insulin (2h)" = 8, "Insulin (12h)" = 9,
                                "Metformin (2h)" = 10, "Metformin (12h)" = 18 )) + 
  ggtitle(paste("PCA Plot Colored by Hierarchical Clusters (k=", optimal_k, ")")) +
  xlab(paste0("PC1 (", round(summary(pca_results)$importance[2, 1] * 100, 1), "%)")) +
  ylab(paste0("PC2 (", round(summary(pca_results)$importance[2, 2] * 100, 1), "%)")) +
  # Control legend order
  guides(
    color = guide_legend(order = 1),
    shape = guide_legend(order = 2)
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold", margin = margin(b = 30)),
    axis.title.x = element_text(size = 12, face = "bold", margin = margin(t = 15, r = 0, b = 0, l = 0)),
    axis.title.y = element_text(size = 12, face = "bold", margin = margin(t = 0, r = 15, b = 0, l = 0)),
    axis.text = element_text(size = 10, face = "bold"),
    plot.margin = unit(c(2, 2, 2, 2), "cm"),
    axis.line = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.7),
    legend.title = element_text(size=12, face="bold"),
    legend.position = "right"
  ) +
  # Define legend labels
  labs(color = paste("H-Clust (k=", optimal_k, ")", sep=""), shape = "Condition")

# Print and save the plot
print(plot11)
ggsave("~/Desktop/plot11.png", plot = plot11, width = 10, height = 8, dpi = 900)

# Create contingency tables to formally compare clustering and conditions
message("Contingency Table: k-means vs Condition")
table_kmeans_condition <- table(KMeans_Cluster = pca_plot_df$Cluster_KMeans, Condition = pca_plot_df$Condition)
print(table_kmeans_condition)

message("Contingency Table: Hierarchical vs Condition")
table_hclust_condition <- table(Hclust_Cluster = pca_plot_df$Cluster_Hclust, Condition = pca_plot_df$Condition)
print(table_hclust_condition)



# --- Task 4: Discriminant Analysis ---

# Prepare data: selected PC scores and the factor defining the groups (Condition)
lda_data <- data.frame(pc_scores_selected, Condition = pca_plot_df$Condition) # Use the combined Condition factor

## Check Parameter Estimation Feasibility

# Define number of groups (g), number of predictor variables (p = #PCs), sample size (n)
g <- length(levels(lda_data$Condition))  # Number of conditions = 7
p <- num_pcs_to_keep                     # Number of PCs used = 2
n <- nrow(lda_data)                      # Number of samples = 28

# Calculate approximate number of parameters LDA needs to estimate
# Formula assumes common covariance matrix across groups
num_params_lda <- g * p + p * (p + 1) / 2 + (g - 1)
# Calculate the ratio of samples to parameters
ratio_n_params <- n / num_params_lda

# Print parameters and ratio for the report
message(paste("LDA Setup: Number of samples (n):", n))
message(paste("LDA Setup: Number of PCs used (p):", p))
message(paste("LDA Setup: Number of groups (g):", g))
message(paste("LDA Setup: Approx. number of parameters to estimate:", round(num_params_lda)))
message(paste("LDA Setup: Ratio of samples to parameters (n/params):", round(ratio_n_params, 2)))

# Check if the ratio is sufficient (should be > 1, ideally > 3-5 for stability)
if (ratio_n_params <= 1) {
  warning("Ratio of samples to parameters is low (<= 1). LDA results may be unstable. Consider using fewer PCs if possible or interpreting results with caution.")
}


## Build and Evaluate the LDA Model

# Build the LDA model to predict 'Condition' using the selected PC scores
# CV = FALSE builds a single model on all data
lda_model <- lda(Condition ~ ., data = lda_data, CV = FALSE)
# Print the LDA model object summary (includes group means, coefficients)
print("LDA Model Summary:")
print(lda_model)

# Evaluate classification accuracy using Leave-One-Out Cross-Validation (LOOCV)
# This trains the model n times, each time leaving one sample out for prediction
lda_cv <- lda(Condition ~ ., data = lda_data, CV = TRUE)

# Create the Confusion Matrix from LOOCV results
confusion_matrix_lda <- table(Predicted = lda_cv$class, Actual = lda_data$Condition)
message("Confusion Matrix (Predicted vs Actual) from LOOCV:")
print(confusion_matrix_lda)

# Calculate Overall Accuracy from the confusion matrix
accuracy_lda <- sum(diag(confusion_matrix_lda)) / sum(confusion_matrix_lda)
message(paste("Overall Accuracy (LOOCV):", round(accuracy_lda * 100, 2), "%"))

# Calculate Overall Misclassification Rate
error_rate_lda <- 1 - accuracy_lda
message(paste("Overall Misclassification Rate (LOOCV):", round(error_rate_lda * 100, 2), "%"))


## Visualise LDA Results

# Get the predicted sample scores on the Linear Discriminants (LDs) from the model
# These scores represent the projection of samples onto the axes that best separate groups
lda_predictions <- predict(lda_model)

# Create a data frame for plotting samples on the first two LDs
lda_plot_df <- data.frame(
  LD1 = lda_predictions$x[, 1],
  LD2 = lda_predictions$x[, 2], # Only if p >= 2 and g >= 3
  Condition = lda_data$Condition
)

# Generate the LDA plot (LD1 vs LD2)
plot12 <- ggplot(lda_plot_df, aes(x = LD1, y = LD2, color = Condition, shape = Condition)) +
  geom_point(size = 3.5, alpha = 0.7) +
  scale_color_manual(values = brewer.pal(7, "Set2")) +
  scale_shape_manual(values = c("Control (0h)" = 15, "Control (2h)" = 16, "Control (12h)" = 17,
                                "Insulin (2h)" = 8, "Insulin (12h)" = 9,
                                "Metformin (2h)" = 10, "Metformin (12h)" = 18 )) + 
  ggtitle("LDA Plot: Samples Projected onto LD1 vs LD2") +
  xlab(paste0("LD1 (Linear Discriminant 1 using PC", 1:p,")")) +
  ylab(paste0("LD2 (Linear Discriminant 2 using PC", 1:p,")")) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold", margin = margin(b = 30)),
    axis.title.x = element_text(size = 12, face = "bold", margin = margin(t = 15, r = 0, b = 0, l = 0)),
    axis.title.y = element_text(size = 12, face = "bold", margin = margin(t = 0, r = 15, b = 0, l = 0)),
    axis.text = element_text(size = 10, face = "bold"),
    plot.margin = unit(c(2, 2, 2, 2), "cm"),
    axis.line = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.7),
    legend.title = element_text(size=12, face="bold"),
    legend.position = "right"
  ) +
  labs(color = "Condition", shape = "Condition")

# Print and save the plot
print(plot12)
ggsave("~/Desktop/plot12.png", plot = plot12, width = 10, height = 8, dpi = 900)