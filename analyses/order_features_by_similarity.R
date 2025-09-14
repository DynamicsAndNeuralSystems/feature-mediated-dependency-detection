# Generating features for random noise process and compute feature correlation
feature_values <- data.frame(matrix(ncol=24))
names(feature_values) <- unique(catch22_all(rnorm(10),catch24=T)$names)
for (seed in 1:50){
  set.seed(seed)
  ts <- rnorm(1000)
  new_ts_feature_values <- catch22_all(ts,catch24 = T)
  new_ts_feature_values <- pivot_wider(new_ts_feature_values,names_from = names, values_from = values)
  feature_values <- rbind.data.frame(feature_values,new_ts_feature_values)
}
feature_values <- feature_values[-1,]

feature_correlation_matrix <- cor(feature_values,method="spearman")

# Perform hierarchical clustering on the feature correlation matrix
feature_cluster <- hclust(as.dist(1 - abs(feature_correlation_matrix)))

dendrogram <- as.dendrogram(feature_cluster)
plot(dendrogram, main = "Feature Clustering Dendrogram", xlab = "Features")

# Get the IDs for each cluster
cluster_ids <- cutree(feature_cluster, k = 5)  

# Order the features based on cluster IDs
ordered_feature_names <- names(feature_values)[order(cluster_ids)]

file_path <- paste0(here("ordered_feature_names.txt"))

# Write the feature names to the text file
write.table(ordered_feature_names, file = file_path, sep = "\t", col.names = FALSE, row.names = FALSE)

