library(monocle3)
# Provide the path to the Cell Ranger output and generate a cell_data_set from 10X output
cds <- load_cellranger_data("./")

# Read in the cluster information generated from cellranger
Cluster <- read.csv("./outs/analysis/clustering/graphclust/clusters.csv", header = TRUE, stringsAsFactors = FALSE)
# Add the cluster information to the 3rd column of cds
colData(cds)[, 3] <- as.factor(Cluster[, 2])
# Read in the LibraryID of each cell
LibraryID <- read.csv("./outs/LibraryID.csv", header = TRUE, stringsAsFactors = FALSE)
# Add the LibraryID to the 4th column of cds
colData(cds)[, 4] <- as.factor(LibraryID[, 2])
# Assign column names of the 2 new added columns
colnames(colData(cds))[3:4] <- c("Cluster", "LibraryID")

# ------subsetting cds------
head(colData(cds))
head(rownames(colData(cds)))
# Remove cells from cluster 5 and 6
valid_cells <- rownames(colData(cds))[!(colData(cds)[, 3] %in% c(5, 6))]
# new cds after removing cells from cluster 5 and 6
cds <- cds[, valid_cells]
# ------done------

# ------violin plot------
library(ggplot2)
library(ggpubr)
# generate a new cds just including interesting genes
cds_subset1 <- cds[row.names(subset(rowData(cds), gene_short_name %in% c("ITGB4", "EPCAM"))),]
cds_subset2 <- cds[row.names(subset(rowData(cds), gene_short_name %in% c("ERBB2", "ERBB3", "EGFR"))),]

plot_genes_violin(cds_subset1, group_cells_by = "Cluster", ncol = 2) + stat_compare_means(method = "anova", label = "p.format", label.x.npc = "center")
plot_genes_violin(cds_subset2, group_cells_by = "Cluster", ncol = 3) + stat_compare_means(method = "anova", label = "p.format", label.x.npc = "center")

# pairwise comparisons with multiple groups
my_comparisons <- list(c("1", "2"), c("1", "3"), c("1", "4"), c("1", "7"), c("1", "8"), c("2", "3"), c("2", "4"), c("2", "7"), c("2", "8"), c("3", "4"), c("3", "7"), c("3", "8"), c("4", "7"), c("4", "8"), c("7", "8"))
plot_genes_violin(cds_subset1, group_cells_by = "Cluster", ncol = 2) + stat_compare_means(method = "wilcox.test", label = "p.format", comparisons = my_comparisons)
plot_genes_violin(cds_subset2, group_cells_by = "Cluster", ncol = 3) + stat_compare_means(method = "wilcox.test", label = "p.format", comparisons = my_comparisons)
# ------done------

# Using PCA to normalize the data, you should specify the number of principal components you want Monocle to compute (here 100).
cds <- preprocess_cds(cds, num_dim = 100)
# To check that you're using enough PCs to capture most of the variation in gene expression across all the cells in the data set. You can look at the fraction of variation explained by each PC using below function
plot_pc_variance_explained(cds)

# To reduce the dimensionality of the data down into the X, Y plane so we can plot it easily
cds <- reduce_dimension(cds)
# plot the data
plot_cells(cds, label_groups_by_cluster = FALSE, show_trajectory_graph = FALSE, color_cells_by = "Cluster", alpha = 0.8, graph_label_size = 2, cell_size = 1) 

# Grouping cells into clusters
cds <- cluster_cells(cds)
# visualization
plot_cells(cds)
# each cell is assigned not only to a cluster but also to a partition. When you are learning trajectories, each partition will eventually become a separate trajectory.
plot_cells(cds, color_cells_by = "partition")

# fit a principal graph within each parition 
cds <- learn_graph(cds)
plot_cells(cds, 
           color_cells_by = "LibraryID", 
           label_cell_groups = FALSE, 
           label_leaves = TRUE,
           label_branch_points = TRUE,
           graph_label_size = 1.5)

library(ggplot2)
plot_cells(cds, 
           color_cells_by = "LibraryID", 
           label_cell_groups = FALSE, 
           label_leaves = FALSE,
           label_branch_points = FALSE,
           graph_label_size = 1.5) + facet_wrap(~LibraryID, nrow = 1)

# Once we've learned a graph, we are ready to order the cells according to their progress through the developmental program. Monocle measures this progress in pseudotime.                      
cds <- order_cells(cds)
plot_cells(cds, 
           color_cells_by = "pseudotime", 
           label_cell_groups = FALSE, 
           label_leaves = FALSE,
           label_branch_points = FALSE,
           graph_label_size = 1.5)

# Retrieve the pseudotime information            
pseudotime <- pseudotime(cds, reduction_method = "UMAP")
pseudotime_df <- data.frame(Barcode = names(pseudotime), pseudotime = as.numeric(pseudotime))
head(pseudotime_df)

write.table(pseudotime_df, "./outs/pseudotime.csv", quote = FALSE, sep = ",", row.names = FALSE, col.names = TRUE)

