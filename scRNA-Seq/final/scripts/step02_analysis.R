# =================================================================
# second step of the analysis
# analysis of Pten PIP-Seq data
# =================================================================

library (Seurat)
library (ggplot2)
library (dplyr)
library (data.table)
library (clusterProfiler)
library (org.Mm.eg.db)
library (SeuratWrappers)
library (openxlsx)
library (RColorBrewer)
library (cowplot)
library (speckle)
library (pheatmap)
library (slingshot)
library (pheatmap)
library (purrr)
library (igraph)
library (ggraph)
library (ggbeeswarm)

# =================
# helper functions
# =================

# ggplot2 style color choice
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

# function to plot the fraction of features in a pool
prep_frac_plot <- function(meta_df = meta_df, sample = "orig.ident", genotype = "genotype", cell_type = "cell_type") {
  
  out_df <- list()
  
  frac_sample <- meta_df %>%
    group_by(.data[[sample]], .data[[genotype]], .data[[cell_type]]) %>%
    summarise(n = n(), .groups = "drop") %>%
    group_by(.data[[sample]]) %>%
    mutate(fraction = n / sum(n))
  
  mean_genotype <- frac_sample %>%
    group_by(.data[[genotype]], .data[[cell_type]]) %>%
    summarise(mean = mean(fraction), .groups = "drop")
  
  frac_plot <- ggplot() +
    geom_bar(
      data = mean_genotype,
      aes(x = .data[[genotype]], y = mean, fill = .data[[cell_type]]),
      stat = "identity"
    ) +
    geom_point(
      data = frac_sample,
      aes(x = .data[[genotype]], y = fraction)
    ) +
    facet_wrap(vars(.data[[cell_type]]), scales = "free_y") +
    theme_classic()
  
  # report plot and base data
  out_df[["fractions"]] <- frac_sample
  out_df[["means"]] <- mean_genotype
  out_df[["plot"]] <- frac_plot
  
  return (out_df)
}

# function to run slingshot on a data set
do_slingshot <- function (x, start_cluster = NA, end_cluster = NA, cluster_label = "seurat_clusters", reduction = "umap", color_label = NA ) {
  
  # initiate container for plots
  plot_list <- list()
  
  # NOTE: clusters with a single cell can't be processed by slingshot!
  
  # slingshot analysis starts here
  sce <- as.SingleCellExperiment(x)
  reducedDim(sce, "UMAP") <- Embeddings(x, reduction)
  
  if (is.na(start_cluster)) {
    sds <- slingshot(sce,
                     clusterLabels = cluster_label,   # column in colData
                     reducedDim    = "UMAP" )
  } else if (any(is.na(end_cluster))) {
    sds <- slingshot(sce,
                     clusterLabels = cluster_label,   # column in colData
                     reducedDim    = "UMAP",
                     start.clus    = start_cluster)
  } else {
    sds <- slingshot(sce,
                     clusterLabels = cluster_label,   # column in colData
                     reducedDim    = "UMAP",
                     start.clus    = start_cluster,
                     end.clus      = end_cluster)     # your progenitor cluster
  }
  
  # Each lineage is an ordered vector of clusters, e.g. c("0", "2", "5")
  lineages <- slingLineages(sds)
  
  umap_df <- as.data.frame( Embeddings(x, reduction = reduction) )
  colnames(umap_df) <- paste("umap", c(1:ncol(umap_df)), sep="_")
  cluster_idx <- which(colnames(x@meta.data) == cluster_label)
  umap_df$cluster <- as.character( x@meta.data[,cluster_idx] )
  
  # define a color for UMAP
  if (!is.na(color_label)) {
    umap_df$color <- unlist( x@meta.data[ color_label ] )
  } else {
    umap_df$color <- unlist( x@meta.data[ cluster_label ] )
  }
  
  centroids <- umap_df |>
    group_by(cluster) |>
    summarise(UMAP_1 = median(umap_1),
              UMAP_2 = median(umap_2),
              .groups = "drop")
  
  edges <- imap_dfr(lineages, function(clusters, lin_name) {
    tibble(
      from    = clusters[-length(clusters)],   # all but last
      to      = clusters[-1],                  # all but first
      Lineage = lin_name
    )
  }) |>
    distinct(from, to, .keep_all = TRUE) |>    # deduplicate shared edges
    left_join(centroids, by = c("from" = "cluster")) |>
    left_join(centroids, by = c("to" = "cluster"))
  
  colnames( edges ) <- c("from", "to", "Lineage", "x", "y", "xend", "yend")
  
  # this part uses ggplot to illustrate the identified trajectories
  
  plot_list[["cluster"]] <- ggplot() +
    geom_point(data = umap_df,
               aes(umap_1, umap_2, color = color),
               size = 0.6 ) +
    
    # Edges between cluster centroids
    geom_segment(data = edges,
                 aes(x = x, y = y, xend = xend, yend = yend),
                 linewidth = 1, color = "black",
                 arrow = arrow(length = unit(0.25, "cm"), type = "closed")) +
    
    # Centroid points on top
    geom_point(data = centroids,
               aes(UMAP_1, UMAP_2),
               size = 3, color = "black", shape = 21, fill = "white", stroke = 1) +
    
    # Cluster labels at centroids
    ggrepel::geom_label_repel(data = centroids,
                              aes(UMAP_1, UMAP_2, label = cluster),
                              size = 3, label.size = 0.2) +
    theme_classic() +
    guides(color = guide_legend(override.aes = list(size = 3)))
  
  curves <- slingCurves(sds, as.df = TRUE)
  
  plot_list[["smooth"]] <- ggplot() +
    # Cells, colored by cluster
    geom_point(data = umap_df,
               aes(x = umap_1, y = umap_2, color = cluster),
               size = 0.6) +
    
    # Trajectory curves
    geom_path(data = curves |> arrange(Order),
              aes(x = umap_1, y = umap_2, group = Lineage),
              linewidth = 1.2, color = "black") +
    
    theme_classic() +
    labs(x = "UMAP 1", y = "UMAP 2")
  
  # plot lineages as tree
  
  # slingPseudotime returns a matrix (one column per lineage)
  pt <- slingPseudotime(sds)
  
  # Add each lineage's pseudotime to Seurat metadata
  for (i in seq_len(ncol(pt))) {
    x <- AddMetaData(x, pt[, i], col.name = paste0("slingPT_lineage", i))
  }
  
  cluster_stats <- umap_df |>
    mutate(cluster = as.character(Idents(x))) |>
    group_by(cluster) |>
    summarise(n_cells = n(), .groups = "drop")
  
  # Mean pseudotime per cluster (across all lineages a cell belongs to)
  mean_pt <- rowMeans(pt, na.rm = TRUE)
  umap_df$mean_pt <- mean_pt
  
  pt_per_cluster <- umap_df |>
    mutate(cluster = as.character(Idents(x))) |>
    group_by(cluster) |>
    summarise(mean_pt = mean(mean_pt, na.rm = TRUE), .groups = "drop")
  
  # Node table — mark root and terminal clusters
  all_nodes <- unique(c(edges$from, edges$to))
  root      <- slingParams(sds)$start.clus
  terminals <- sapply(lineages, tail, 1) |> unique()
  
  nodes <- tibble(name = all_nodes) |>
    mutate(type = case_when(
      name == root        ~ "root",
      name %in% terminals ~ "terminal",
      TRUE                ~ "intermediate"
    ))
  
  nodes <- nodes |>
    left_join(cluster_stats, by = c("name" = "cluster")) |>
    left_join(pt_per_cluster, by = c("name" = "cluster"))
  
  g <- graph_from_data_frame(edges, vertices = nodes, directed = TRUE)
  
  # tree type output
  plot_list[["tree_1"]] <- ggraph(g, layout = "tree") +
    
    geom_edge_link(
      arrow     = arrow(length = unit(0.3, "cm"), type = "closed"),
      end_cap   = circle(8, "mm"),
      color     = "grey40",
      linewidth = 0.8
    ) +
    
    geom_node_point(aes(color = mean_pt, size = n_cells)) +
    
    geom_node_text(aes(label = paste0("C", name, "\nn=", n_cells)),
                   fontface = "bold", color = "white", size = 3) +
    
    scale_color_viridis_c(name = "Mean\npseudotime") +
    scale_size_continuous(name = "# cells", range = c(8, 18)) +
    
    theme_graph(base_family = "sans") +
    labs(title = "Slingshot lineage topology",
         subtitle = "Node size = cluster size · color = mean pseudotime")
 
  return(plot_list)
}

# =================================
# define base folders
# =================================

base <- "./scRNA-Seq/final"
no_save_base <- "./Miranda_et_al_noSave"
count_base <- paste( base, "PIPSeeker_output/", sep="" ) # folder where PIPSeeker output is located

set.seed(2401)

# ===============================
# Define a unified color scheme
# ===============================
color_vec <- c( gg_color_hue(5), "grey50" )
names(color_vec) <- c("aIP","RGP", "Astro", "OBNB", "Oligo", "Neuron")

# =======================================================================
# convert sample idx to genotype
# KO: sparse KO, i.e. green cells in a Pten-MADM
# DKO: double KO, sparse Pten ko in a Egfr cKO background
# =======================================================================
conv_vec <- c("KO1", "KO2", "DKO1", "WT1", "KO3", "DKO2", "WT2", "KO4")

# =================================
# read raw data and do basic QC
# and filtering
# =================================

seurat_obj_list <- list()
# for nice plotting order later
for (sample_idx in c(1:8) ) {
  
  message( paste("sample", sample_idx, sep=""))
  
  # we use sensitivity3 analysis option here
  file_path <- paste( count_base, "sample", sample_idx, "/filtered_matrix/sensitivity_3", sep="" )
  message( file_path )  
  
  counts <- Read10X( data.dir = file_path,  )
  seurat_obj_list[[ sample_idx ]] <- CreateSeuratObject(counts = counts, project = conv_vec[ sample_idx ], min.cells = 3, min.features = 200)
}

# merging objects
seurat.obj <- merge( x = seurat_obj_list[[ 1 ]], y = seurat_obj_list[ -1 ] )

# this par reproduces QC images from initial analysis
message("preparing QC plots")

# first 2 QC plots are just for basic information, not for the paper
# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
seurat.obj[["percent.mt"]] <- PercentageFeatureSet(seurat.obj, pattern = "^mt-")

# Visualize QC metrics as a violin plot
vln_plot <- VlnPlot(seurat.obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
ggsave( plot=vln_plot, file=paste (base, "/plots/QC/QC_basicFeatures.Vln.raw.pdf", sep=""), width = 10, height = 7)

message("filtering")
seurat.obj <- subset(seurat.obj, subset = nFeature_RNA > 500 & nFeature_RNA < 6000 & nCount_RNA < 40000 & percent.mt < 15)
vln_plot <- VlnPlot(seurat.obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

ggsave( plot=vln_plot, file=paste (base, "/plots/QC/QC_basicFeatures.Vln.filtered.pdf", sep=""), width = 12, height = 5)

plot1 <- FeatureScatter(seurat.obj, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(seurat.obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
QC_plot <- plot1 + plot2

ggsave( plot=QC_plot, file=paste (base, "/plots/QC/QC_combined.QC.pdf", sep=""), width = 10, height = 7)

# =================================
# standard Seurat pipeline
# =================================

seurat.obj <- NormalizeData(seurat.obj, normalization.method = "LogNormalize", scale.factor = 10000)
seurat.obj <- FindVariableFeatures(seurat.obj, selection.method = "vst", nfeatures = 2000)

seurat.obj <- ScaleData(seurat.obj)
seurat.obj <- RunPCA(seurat.obj, features = VariableFeatures(object = seurat.obj))

elbowPlot  <- ElbowPlot( object = seurat.obj, ndims = 30)
ggsave( plot=elbowPlot, file=paste (base, "/plots/QC/QC_ElbowPlot.QC.pdf", sep=""), width = 10, height = 7)

seurat.obj <- FindNeighbors(seurat.obj, dims = 1:25, reduction = "pca")
seurat.obj <- FindClusters(seurat.obj, resolution = 0.3)

seurat.obj <- RunUMAP(seurat.obj, dims = 1:25)

# for internal QC: minor batch effects visible
# replicates initial analysis
umap_replicates <- DimPlot( seurat.obj, group.by = c("orig.ident"), label=F, reduction = "umap" )
ggsave( filename = paste(base, "/plots/QC/QC_replicate_umap.pdf", sep=""), plot = umap_replicates, width=5.5, height=5 )

# ======================================
# removal of a small unspecific cluster
# ======================================

# for documentation purposes - one cluster will be removed, visualized here 
# replicates initial analysis
umap_clusters <- DimPlot( seurat.obj, group.by = c("seurat_clusters"), label=T) + NoLegend()
ggsave( filename = paste(base, "/plots/QC/QC_cluster_umap.pdf", sep=""), plot = umap_clusters, width=5, height=5 )

# prepare a plot to show that Seurat clusters are not majorly biased by batch effects
# exceptions are cluster 1, where WT2 shows an increased fraction and cluster 7 where DKO2 dominates
result <- seurat.obj@meta.data %>%
  group_by(seurat_clusters, orig.ident) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(seurat_clusters) %>%
  mutate(fraction = n / sum(n))

ggplot( result, aes(x=seurat_clusters, y=fraction, fill = orig.ident)) + geom_bar( stat="identity" ) + theme_classic()
ggsave( filename = paste(base, "/plots/QC/QC_cluster_replicate_fraction.pdf", sep=""), width=5, height=5 )

# write out seurat cluster proportions
write.csv( x = result, file = paste(base, "/plots/QC/QC_cluster_replicate_fraction.csv", sep=""))

# identify cell types
joined.obj <- JoinLayers( seurat.obj )
# create a variable that is ordered in a convenient way
joined.obj@meta.data$plot_id <- joined.obj@meta.data$orig.ident
joined.obj@meta.data$plot_id <- factor( joined.obj@meta.data$plot_id, levels = c("WT1", "WT2", "KO1", "KO2", "KO3", "KO4", "DKO1", "DKO2") )

# identify markers
# internal QC to rationalize removal of Cluster7
cluster_markers <- FindAllMarkers( joined.obj, only.pos = T)

cluster_markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10

heat_plot <- DoHeatmap( joined.obj, features = top10$gene ) + NoLegend()
ggsave( filename = paste(base, "/plots/QC/QC_marker_heatmap.pdf", sep=""), plot = heat_plot, width=8, height=9 )

# cluster 7 is not well localized and shows very few specific markers
# remove from analysis
seurat.obj <- subset(seurat.obj, seurat_clusters %in% c(7), invert=T)

# =======================
# Harmony integration
# =======================

seurat.obj <- IntegrateLayers(
  object = seurat.obj, method = HarmonyIntegration,
  orig.reduction = "pca", new.reduction = "harmony",
  verbose = FALSE
)

seurat.obj <- FindNeighbors(seurat.obj, reduction = "harmony", dims = 1:30)
seurat.obj <- FindClusters(seurat.obj, resolution = 0.3, cluster.name = "harmony_clusters")

seurat.obj <- RunUMAP(seurat.obj, reduction = "harmony", dims = 1:25, reduction.name = "umap.harmony")

DimPlot( seurat.obj, group.by = c("harmony_clusters"), reduction = "umap.harmony", label=T) + NoLegend()
ggsave( filename = paste(base, "/plots/Supplement/Sup14Edown_harmony_integrated_umap_clusters.pdf", sep=""), width=5, height=5 )

# ======================================================================================
# prepare a plot to show that Seurat clusters are not majorly affected by batch effects
# ======================================================================================

result <- seurat.obj@meta.data %>%
  group_by(harmony_clusters, orig.ident) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(harmony_clusters) %>%
  mutate(fraction = n / sum(n))

# write out seurat cluster proportions
write.csv( x = result, file = paste(base, "/plots/QC/QC_cluster_fraction_integrated.csv", sep=""))

ggplot( result, aes(x=harmony_clusters, y=fraction, fill = orig.ident)) + geom_bar( stat="identity" ) + theme_classic()
ggsave( filename = paste(base, "/plots/QC/QC_cluster_fraction_integrated.pdf", sep=""), width=5, height=5 )

# add a factor for nice sample order
seurat.obj@meta.data$sampleID <- factor( seurat.obj@meta.data$orig.ident, levels = c("WT1", "WT2", "KO1", "KO2", "KO3", "KO4", "DKO1", "DKO2") )

# ========================
# some QC for supplement
# ========================

# create a genotype variable
seurat.obj@meta.data$genotype <- seurat.obj@meta.data$orig.ident
seurat.obj@meta.data$genotype <- sapply( seurat.obj@meta.data$genotype, function (x) gsub(pattern = "1|2|3|4", replacement = "", x = x))

# create a color code for genotypes
#color_palette <- c("WT" = "Greys", "DKO" = "Blues", "KO" = "Oranges" )
color_palette <- c("WT" = "grey50", "DKO" = "cornflowerblue", "KO" = "white" )

colors <- sapply( unique(seurat.obj@meta.data$genotype), function (x) {
  n_col <- length( unique( seurat.obj@meta.data$orig.ident[ which(seurat.obj@meta.data$genotype == x) ] ) )
  color_list <- rep( color_palette[[x]], n_col ) 
  return(color_list)
})

vln_plot <- VlnPlot(seurat.obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, group.by = "sampleID")
vln_plot <- vln_plot & scale_fill_manual( values = unlist(colors) )

ggsave( plot=vln_plot, file=paste (base, "/plots/Supplement/Sup14B-D_basicFeatures.Vln.filtered.pdf", sep=""), width = 12, height = 5)

# ===============================================
# marker expression and cell type identification
# ===============================================

marker_plot <- FeaturePlot( seurat.obj, features = c("Dcx", "Neurod2", "Gad2", "Fabp7"), 
                            order = T, min.cutoff = "q15", reduction = "umap.harmony" ) & NoLegend()
ggsave( filename = paste(base, "/plots/Supplement/Sup14F_marker_expression_plots.pdf", sep=""), plot = marker_plot, width=9, height=9 )

non_proj_markers <- c("Mki67", "Pdgfra", "Olig2", "Cdk1", "Top2a", "Aldh1l1", "Sox10", "Olig1", "Gad2", "Dlx5", "Hes5", "Aldoc", "Fabp7", "Pcna", "Mcm3", "Sparc", "Ecrg4", "Gfap", "Slc1a3", "Aqp4")
proj_markers <- c("Bcl11b", "Foxp2", "Fezf2", "Satb2", "Eomes", "Neurod1","Cux1","Cux2" )

cluster_pseudobulk <- PseudobulkExpression( seurat.obj, group.by = "seurat_clusters")
pheatmap( cluster_pseudobulk[["RNA"]][c( proj_markers, non_proj_markers ) , ], scale="row",
          file = paste(base, "/plots/QC/QC_broad_marker_heatmap_clusters.pdf", sep=""))

# close drawing device if necessary, otherwise plotting goes funny
# wrap in tryCatch to allow for code to continue even if drawing device was closed before 
tryCatch({
  dev.off()
}, error = function(e) {
  message("An error occurred: ", e$message)
})

# give useful cell type names
cluster2celltype <- c("0" = "Neuron", "1" = "Neuron", "2" = "aIP", "3" = "RGP", "4" = "Oligo", "5" = "OBNB", "6" = "Neuron")

cluster2broad_celltype <- c("0" = "Projection Neurons", "1" = "Projection Neurons", 
                      "2" = "Non-Projection Neuron Cells", "3" = "Non-Projection Neuron Cells", 
                      "4" = "Non-Projection Neuron Cells", "5" = "Non-Projection Neuron Cells", 
                      "6" = "Projection Neurons")

seurat.obj@meta.data$broad_cell_type <- cluster2broad_celltype[ as.character(seurat.obj@meta.data$seurat_clusters) ]
seurat.obj@meta.data$cell_type <- cluster2celltype[ as.character(seurat.obj@meta.data$seurat_clusters) ]

DimPlot( seurat.obj, group.by = "broad_cell_type", reduction = "umap.harmony")
ggsave( filename = paste(base, "/plots/QC/QC_umap_integrated_broad_cell_types.pdf", sep=""), width=8, height=6 )

DimPlot( seurat.obj, group.by = "orig.ident", reduction = "umap.harmony")
ggsave( filename = paste(base, "/plots/Supplement/Sup14Eup_sample_umap_integrated.pdf", sep=""), width=7, height=6 )

DimPlot( seurat.obj, group.by = "cell_type", reduction = "umap.harmony") + scale_color_manual( values = color_vec)
ggsave( filename = paste(base, "/plots/Sup14G_umap_cell_types.pdf", sep=""), width=8, height=6 )

# write out the meta data
write.csv( x = seurat.obj@meta.data, file = paste(base, "/Supplement/scRNA_meta_data.csv", sep=""))

# ================================
# plot the cell number per sample
# ================================

total_cells <- table( seurat.obj@meta.data$orig.ident)

df2plot <- as.data.frame(total_cells)
df2plot$Var1 <- factor( df2plot$Var1, levels = c("WT1", "WT2", "KO1", "KO2", "KO3", "KO4", "DKO1", "DKO2") )
df2plot$group <- sapply( df2plot$Var1, function (x) gsub(pattern = "1|2|3|4", replacement = "", x = x))

df2plot %>%
  group_by( group) %>%
  summarize( mean = mean(Freq), sum = sum(Freq) ) -> av_group

df2plot$group <- factor( df2plot$group, levels = c("WT", "KO", "DKO") )

cell_number_plot <- ggplot( df2plot, aes(x=group, y=Freq, fill=group)) + 
  geom_bar( data = av_group, aes(x=group, y=mean), stat="identity", color="black") +
  geom_beeswarm( ) + theme_classic() + scale_fill_manual( values = c("WT" = "grey59", "DKO" = "cornflowerblue", "KO" = "white")) + 
  ylab( "# of cells" ) + NoLegend()
  
ggsave( filename = paste(base, "/plots/Supplement/Sup14A_cell_number_per_sample.pdf", sep=""), plot = cell_number_plot, width=3, height=5 )

# ==========================================================
# refine annotation
# ==========================================================

# in the next steps I extract and analyse projection neurons and non projection neurons separately
# goal is to refine annotation 

# ============================================================================================
# investigate neurons specifically to check for the layer specific neuron abundance phenotype
# ============================================================================================

# extract neuronal lineage and recluster
Idents( seurat.obj ) <- "broad_cell_type"
neuron.seurat.obj <- subset( seurat.obj, idents = c("Projection Neurons"))

neuron.seurat.obj <- FindNeighbors(neuron.seurat.obj, dims = 1:20, reduction = "harmony")
neuron.seurat.obj <- FindClusters(neuron.seurat.obj, resolution = 0.4)

neuron.seurat.obj <- RunUMAP(neuron.seurat.obj, dims = 1:20, reduction = "harmony")

neuron_umap_clusters <- DimPlot( neuron.seurat.obj, group.by = c("seurat_clusters"), label=T)
ggsave( filename = paste(base, "/plots/reviewer/Rev_Fig2A_neuron_umap_clusters.pdf", sep=""), plot = neuron_umap_clusters, width=6, height=5 )

# plot marker gene expression as heatmap to identify UL/DL/IP
cluster_pseudobulk <- PseudobulkExpression( neuron.seurat.obj, group.by = "seurat_clusters")
pheatmap( cluster_pseudobulk[["RNA"]][c("Bcl11b", "Foxp2", "Fezf2", "Satb2", "Eomes", "Neurod1","Cux1","Cux2" ) , ], 
          scale="row", file = paste(base, "/plots/reviewer/Rev_Fig2B_neuron_marker_heatmap_clusters.pdf", sep=""), silent = T )

# close drawing device if necessary, otherwise plotting goes funny
# wrap in tryCatch to allow for code to continue even if drawing device was closed before 
tryCatch({
  dev.off()
}, error = function(e) {
  message("An error occurred: ", e$message)
})

# assign cell types based on marker gene expression
cluster2cell_type <- c("0" = "UL", "1" = "UL", "2" = "UL", "3" = "IP", "4" = "DL", "5" = "UL" )
neuron.seurat.obj@meta.data$cell_type <- cluster2cell_type[ as.character(neuron.seurat.obj$seurat_clusters)]
DimPlot( neuron.seurat.obj, group.by = c("cell_type"), label=T)
ggsave( filename = paste(base, "/plots/reviewer/Rev_Fig2C_neuron_umap_cellTypes.pdf", sep=""), width=6, height=5 )

# sanity check of marker gene expression in cell types
neuron_pseudobulk <- PseudobulkExpression( neuron.seurat.obj, group.by = "cell_type")
pheatmap( neuron_pseudobulk[["RNA"]][c("Bcl11b", "Foxp2", "Fezf2", "Satb2", "Eomes", "Neurod1","Cux1","Cux2" ) , ], 
          scale="row", file = paste(base, "/plots/reviewer/Rev_QC_neuron_marker_heatmap_cellType.pdf", sep=""), silent = T )

# close drawing device if necessary, otherwise plotting goes funny
# wrap in tryCatch to allow for code to continue even if drawing device was closed before 
tryCatch({
  dev.off()
}, error = function(e) {
  message("An error occurred: ", e$message)
})

# do a Neuron specific analysis - excluding IPs
Idents(neuron.seurat.obj) <- "cell_type"
neurons_only <- subset( neuron.seurat.obj, idents = c("DL", "UL"))

# plot the cell type abundances
neuron_frac_df <- prep_frac_plot ( meta_df = neurons_only@meta.data )

prop_out  <- propeller(clusters = neurons_only$cell_type, 
                       sample = neurons_only$orig.ident, 
                       group = neurons_only$genotype, 
                       transform = "logit")

# Performing logit transformation of proportions
# group variable has > 2 levels, ANOVA will be performed
# Robust eBayes doesn't work with fewer than 3 cell types.
#                     Setting robust to FALSE
                    
prop_out$cell_type <- rownames( prop_out )
neuron_frac_df[["propeller"]] <- prop_out

# write out the raw data and anova for this plot
wb <- createWorkbook()

for (name in c("fractions", "means", "propeller")) {
  sheetName <- name
  tmp <- neuron_frac_df[[ name ]]
  addWorksheet(wb, sheetName)
  writeData(wb, sheetName, tmp)
  addFilter(wb, sheetName, row = 1, cols = 1:ncol(tmp))
  setColWidths(wb, sheetName, cols = 1:ncol(tmp), widths="auto")  
}

saveWorkbook(wb, file = paste(base, "/plots/reviewer/Rev_Fig2D_neuron_ULDLfrac_plot_data.xlsx", sep="/"), overwrite = T) 

ggsave( filename = paste(base, "/plots/reviewer/Rev_Fig2D_neuron_ULDLfrac_plot.pdf", sep=""), plot = neuron_frac_df[["plot"]], width=5, height=4 )

# ================================
# focus analysis on glia/OBNB
# ================================

# extract non projection neurons
seurat_glia_list <- lapply( seurat_obj_list, function (x) {

  # do a clean re-analysis of non projection neuron cells
  x <- subset(x, cells = rownames( seurat.obj@meta.data[which(seurat.obj@meta.data$broad_cell_type %in% c("Non-Projection Neuron Cells")),]))
  batch_vec <- x@meta.data$orig.ident
  
  x <- CreateSeuratObject( counts = GetAssay( object = x, slot = "counts", min.cells = 3, min.features = 200 ))
  x@meta.data$batch <- batch_vec
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x)
  
  return(x)
  
})

seurat.glia.obj <- merge( x = seurat_glia_list[[ 1 ]], y = seurat_glia_list[ -1 ] )

seurat.glia.obj <- AddMetaData( object = seurat.glia.obj, metadata = seurat.obj@meta.data[,c("cell_type", "genotype"), drop=F])

seurat.glia.obj <- ScaleData(seurat.glia.obj)

seurat.glia.obj <- RunPCA(seurat.glia.obj)
seurat.glia.obj <- RunUMAP(seurat.glia.obj, dims = 1:15, reduction.name = "non_integrated")

# integrate layers
seurat.glia.obj <- IntegrateLayers(
  object = seurat.glia.obj, method = HarmonyIntegration,
  orig.reduction = "pca", new.reduction = "harmony",
  verbose = FALSE
)

seurat.glia.obj <- FindNeighbors(seurat.glia.obj, dims = 1:10, reduction = "harmony")
seurat.glia.obj <- FindClusters(seurat.glia.obj, resolution = 1)
seurat.glia.obj <- RunUMAP(seurat.glia.obj, dims = 1:10, reduction.name = "umap", reduction = "harmony")

# ========================================================================================
# do a set of analyses to show that cluster 12 is not well defined and should be removed
# ========================================================================================

# determine cluster - cell type association
seurat.glia.obj@meta.data %>%
  group_by(seurat_clusters, cell_type) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(seurat_clusters) %>%
  mutate(fraction = n / sum(n)) -> cellType_cluster

# write out seurat cluster proportions
write.csv( x = cellType_cluster, file = paste(base, "/plots/QC/QC_glia_cell_type_fraction.csv", sep=""))

table(seurat.glia.obj@meta.data$seurat_clusters)

#   0   1   2   3   4   5   6   7   8   9  10  11  12  13 
# 288 279 264 240 223 199 138 130 127 114 109  74  47  22 

# all clusters have >80% association with a cell type
# except for cluster 12, where it is only 60%
cellType_cluster %>%
  group_by(seurat_clusters) %>%
  slice_max(fraction, n = 1) -> cluster_ct_frac

# write the data out
write.table( x = cluster_ct_frac, 
             file = paste(base, "/plots/QC_glia_ct_cluster_frac.tsv", sep=""), 
             row.names = F )

# plot the data
ggplot( cellType_cluster, aes(x=seurat_clusters, y=fraction, fill=cell_type)) + 
  geom_bar(stat="identity") + theme_classic()
ggsave( filename = paste(base, "/plots/QC/QC_glia_ct_cluster_frac.pdf", sep="") )

# =============================================================
# pseudobulk marker expression to refine cell type annotation
# =============================================================
glia_pseudobulk <- PseudobulkExpression( seurat.glia.obj, group.by = "seurat_clusters")
pheatmap( glia_pseudobulk[["RNA"]][c("Pdgfra", "Olig2", "Cdk1", "Top2a", "Aldh1l1", "Sox10", 
                                     "Olig1", "Gad2", "Dlx5", "Hes5", "Aldoc", "Fabp7", "Pcna", 
                                     "Mcm3", "Sparc", "Ecrg4", "Gfap", "Slc1a3", "Aqp4") , ], 
          scale="row", file = paste(base, "/plots/QC/QC_glia_marker_heatmap_clusters.pdf", sep="") )

# close drawing device if necessary, otherwise plotting goes funny
# wrap in tryCatch to allow for code to continue even if drawing device was closed before 
tryCatch({
  dev.off()
}, error = function(e) {
  message("An error occurred: ", e$message)
})

# assign new cell type labels
# note cluster 12 has no clear marker gene pattern
# note that 13 clearly expresses markers of oligos
glia_cluster_2_cellType <- c("0" = "RGP", "1" = "OBNB", "2" = "Astro", "3" = "Oligo", 
                             "4" = "aIP", "5" = "aIP", "6" = "aIP", 
                             "7" = "RGP", "8" = "Astro", "9" = "RGP", 
                             "10" = "Oligo", "11" = "RGP", "13" = "Oligo", "12" = "mixed")

seurat.glia.obj@meta.data$cell_type <- glia_cluster_2_cellType[ as.character(seurat.glia.obj@meta.data$seurat_clusters) ]

# ====================================================================
# detour: annotate refined neuron/glia cell types in complete dataset
# ====================================================================
seurat.obj@meta.data$cellID <- rownames(seurat.obj@meta.data)
seurat.glia.obj@meta.data$cellID <- rownames(seurat.glia.obj@meta.data)

seurat.obj@meta.data$cell_type[ match(seurat.glia.obj@meta.data$cellID, seurat.obj@meta.data$cellID) ] <- seurat.glia.obj$cell_type

# annotate neuron cells in complete dataset
neuron.seurat.obj@meta.data$cellID <- rownames(neuron.seurat.obj@meta.data)
seurat.obj@meta.data$cell_type[match(neuron.seurat.obj@meta.data$cellID, seurat.obj@meta.data$cellID)] <- neuron.seurat.obj@meta.data$cell_type

umap_plot <- DimPlot( seurat.obj, group.by = "cell_type", reduction = "umap.harmony")
ggsave( filename = paste(base, "/plots/QC/QC_refined_ct_umap.pdf", sep=""), plot = umap_plot, width=5, height=4 )

# save the Seurat object for further analysis
saveRDS( object = seurat.obj, file = paste(no_save_base, "/RDS_files/rev_seurat.obj.rds", sep="") )
# detour end

# ====================================================================
# analyse fraction of cells in glia data set
# ====================================================================

# remove mixed seurat cluster
seurat.glia.obj <- subset( seurat.glia.obj, idents = 12, invert = T)

# save the glia object for further analysis
saveRDS( object = seurat.glia.obj, file = paste(no_save_base, "/RDS_files/rev.seurat.glia.obj.rds", sep="") )

# draw marker expression maps
marker_plot <- FeaturePlot( seurat.glia.obj, features = c("Cdk1", "Mki67", "Aldh1l1", "Gfap", "Pdgfra", "Cspg4", "Dlx5", "Gad2"), 
                            order = T, min.cutoff = "q15", ncol = 4 ) & NoLegend()
ggsave( filename = paste(base, "/plots/Supplement/Sup14H_glia_marker_expression_plots.pdf", sep=""), 
        plot = marker_plot, width=15, height=9 )

# sanity check - pseuobulk expression of marker genes in cell types
glia_pseudobulk <- PseudobulkExpression( seurat.glia.obj, group.by = "cell_type")
pheatmap( glia_pseudobulk[["RNA"]][c("Pdgfra", "Olig2", "Cdk1", "Top2a", "Aldh1l1", "Sox10", "Olig1", "Gad2", "Dlx5", "Hes5", "Aldoc", "Fabp7", "Pcna", "Mcm3", "Sparc", "Ecrg4", "Gfap", "Slc1a3", "Aqp4") , ], 
          scale="row", file = paste(base, "/plots/main/Fig4B_glia_marker_heatmap_cellType.pdf", sep="") )

# close drawing device if necessary, otherwise plotting goes funny
# wrap in tryCatch to allow for code to continue even if drawing device was closed before 
tryCatch({
  dev.off()
}, error = function(e) {
  message("An error occurred: ", e$message)
})

DimPlot(seurat.glia.obj, reduction = "umap", group.by = "batch")
ggsave( filename = paste(base, "/plots/QC/QC_glia_umap_replicates.pdf", sep=""), width=5, height=4 )

DimPlot(seurat.glia.obj, reduction = "umap", group.by = "seurat_clusters", label = T)
ggsave( filename = paste(base, "/plots/QC/QC_glia_umap_clusters.pdf", sep=""), width=5, height=4 )

DimPlot(seurat.glia.obj, reduction = "umap", group.by = "cell_type", label = T) + scale_color_manual( values = color_vec )
ggsave( filename = paste(base, "/plots/main/Fig4B_glia_umap_cellType.pdf", sep=""), width=5, height=4 )

seurat.glia.obj@meta.data$cell_type <- factor( seurat.glia.obj@meta.data$cell_type, 
                                               levels = c("RGP", "aIP", "Astro", "Oligo", "OBNB"))

# ==============================================================
# add cell cycle information to the Seurat object and plot UMAP
# ==============================================================

# join kayers
seurat.glia.obj <- JoinLayers( seurat.glia.obj )

# read cell cyle genes
cell_cycle_phase_gene <- read.table(file = paste(base, "other_data/Mus_musculus.csv", sep=""), sep = ",", header = T)
gene_conversion <- bitr(geneID = cell_cycle_phase_gene$geneID, fromType = "ENSEMBL", toType = "SYMBOL", OrgDb = org.Mm.eg.db)
cell_cycle_phase_gene <- merge(cell_cycle_phase_gene, gene_conversion, by.x="geneID", by.y="ENSEMBL")
seurat.glia.obj <- CellCycleScoring (
   object = seurat.glia.obj,
   g2m.features = cell_cycle_phase_gene[which(cell_cycle_phase_gene$phase == "G2/M"), "SYMBOL"],
   s.features = cell_cycle_phase_gene[which(cell_cycle_phase_gene$phase == "S"), "SYMBOL"]
)

seurat.glia.obj@meta.data$genotype <- factor( seurat.glia.obj@meta.data$genotype, levels = c("WT", "KO", "DKO"))
DimPlot( seurat.glia.obj, group.by = "Phase", split.by = "genotype") + ggtitle("")
ggsave( filename = paste(base, "/plots/Supplement/Sup14I_glia_phase_plot_genotype.pdf", sep=""), width=9, height=4 )

# ====================================================
# prepare a slingshot analysis for broad trajectories
# ====================================================
Idents(seurat.glia.obj) <- "genotype"
seurat.WT <- subset(seurat.glia.obj, idents = "WT")
slingshot_out <- do_slingshot( seurat.WT, start_cluster = "RGP", end_cluster = NA, cluster_label = "cell_type" )
plot1 <- slingshot_out[[1]] + ggtitle("WT") + NoLegend() + scale_color_manual( values = color_vec )

seurat.KO <- subset(seurat.glia.obj, idents = "KO")
slingshot_out <- do_slingshot( seurat.KO, start_cluster = "RGP", end_cluster = NA, cluster_label = "cell_type" )
plot2 <- slingshot_out[[1]] + ggtitle("KO") + NoLegend() + scale_color_manual( values = color_vec )

seurat.DKO <- subset(seurat.glia.obj, idents = "DKO")
slingshot_out <- do_slingshot( seurat.DKO, start_cluster = "RGP", end_cluster = c("Oligo"), cluster_label = "cell_type" )
plot3 <- slingshot_out[[1]] + ggtitle("DKO") + NoLegend() + scale_color_manual( values = color_vec )

comb_traj_plot <- plot1 + plot2 + plot3
ggsave( filename = paste(base, "/plots/main/Fig4C_glia_broad_traj.pdf", sep=""), plot = comb_traj_plot, width = 15, height = 5)

# ====================================================
# prepare a slingshot analysis for fine trajectory
# ====================================================
slingshot_out <- do_slingshot( seurat.glia.obj, start_cluster = "7", end_cluster = "1", cluster_label = "seurat_clusters", color_label = "cell_type" )
slingshot_out[[1]] <- slingshot_out[[1]] + scale_color_manual( values = color_vec )
ggsave( filename = paste(base, "/plots/reviewer/Rev_Fig3A_glia_fine_traj.pdf", sep=""), plot = slingshot_out[[1]], width = 6, height = 5)

# =====================================================================
# add the glia seurat clusters to the Seurat object holding all cells
# this allows to calculate relative abundances relative to all cells 
# =====================================================================

# Neurons are treated as 1 entity
seurat.obj@meta.data$glia_clusters <- "Neuron"
seurat.obj@meta.data$glia_clusters[ match(seurat.glia.obj@meta.data$cellID, seurat.obj@meta.data$cellID) ] <- as.character( seurat.glia.obj$seurat_clusters )

# remove the mixed cell cluster
Idents( seurat.obj ) <- "cell_type"
seurat.obj <- subset(seurat.obj, idents = "mixed", invert = T)

# do the relative abundance plot and statistical testing
glia_cluster_frac_relTotal_df <- prep_frac_plot ( meta_df = seurat.obj@meta.data, cell_type = "glia_clusters", sample = "orig.ident" )
prop_out  <- propeller(clusters = seurat.obj$glia_clusters, 
                       sample = seurat.obj$orig.ident, 
                       group = seurat.obj$genotype, 
                       transform = "logit")

# Performing logit transformation of proportions
# group variable has > 2 levels, ANOVA will be performed

prop_out$cell_type <- rownames( prop_out )
glia_cluster_frac_relTotal_df[["propeller"]] <- prop_out

# write out the raw data and statistics for this plot
wb <- createWorkbook()

for (name in c("fractions", "means", "propeller")) {
  sheetName <- name
  tmp <- glia_cluster_frac_relTotal_df[[ name ]]
  addWorksheet(wb, sheetName)
  writeData(wb, sheetName, tmp)
  addFilter(wb, sheetName, row = 1, cols = 1:ncol(tmp))
  setColWidths(wb, sheetName, cols = 1:ncol(tmp), widths="auto")  
}

saveWorkbook(wb, file = paste(base, "/plots/reviewer/Rev_Fig3B_glia_cluster_plot_relTotal_data.xlsx", sep="/"), overwrite = T) 


glia_cluster_frac_relTotal_df[["plot"]] <- glia_cluster_frac_relTotal_df[["plot"]] + NoLegend()
ggsave( filename = paste(base, "/plots/reviewer/Rev_Fig3B_glia_cluster_plot_relTotal.pdf", sep=""), 
        glia_cluster_frac_relTotal_df[["plot"]], width=7, height=7 ) 

# relative abundance of glia cell types relative to all cells
seurat.obj@meta.data$glia_ct <- seurat.obj@meta.data$cell_type
seurat.obj@meta.data$glia_ct[ which( seurat.obj@meta.data$cell_type %in% c("UL", "DL", "IP")) ] <- "Neuron"

seurat.obj@meta.data$glia_ct <- factor( seurat.obj@meta.data$glia_ct, levels = c("RGP", "aIP", "Astro", "Oligo", "OBNB", "Neuron"))
glia_ct_frac_relTotal_df <- prep_frac_plot ( meta_df = seurat.obj@meta.data, cell_type = "glia_ct", sample = "orig.ident" )

prop_out  <- propeller(clusters = seurat.obj$glia_ct, 
                       sample = seurat.obj$orig.ident, 
                       group = seurat.obj$genotype, 
                       transform = "logit")

# Performing logit transformation of proportions
# group variable has > 2 levels, ANOVA will be performed

prop_out$cell_type <- rownames( prop_out )
glia_ct_frac_relTotal_df[["propeller"]] <- prop_out

# write out the raw data and anova for this plot
wb <- createWorkbook()

for (name in c("fractions", "means", "propeller")) {
  sheetName <- name
  tmp <- glia_ct_frac_relTotal_df[[ name ]]
  addWorksheet(wb, sheetName)
  writeData(wb, sheetName, tmp)
  addFilter(wb, sheetName, row = 1, cols = 1:ncol(tmp))
  setColWidths(wb, sheetName, cols = 1:ncol(tmp), widths="auto")  
}

saveWorkbook(wb, file = paste(base, "/plots/main/Fig4D_glia_cellType_plot_relTotal_data.xlsx", sep="/"), overwrite = T) 

glia_ct_frac_relTotal_df[["plot"]] <- glia_ct_frac_relTotal_df[["plot"]] + NoLegend() + scale_fill_manual( values = color_vec )
ggsave( filename = paste(base, "/plots/main/Fig4D_glia_cellType_plot_relTotal.pdf", sep=""), 
        glia_ct_frac_relTotal_df[["plot"]], width=7, height=7 ) 


# =================================================================================
# plot the relative abundances of seurat clusters based on trajectory information
# =================================================================================

trajectory_list <- list()
trajectory_list[["Oligo"]] <- c(7,0,11,10,3,13)
trajectory_list[["OBNB"]] <- c(7,1)
trajectory_list[["Astro"]] <- c(7,9,6,4,5,8,2)

traj_abundance_plot <- list()
for (cell_type in names(trajectory_list)) {
  
  df <- glia_cluster_frac_relTotal_df[[ "fractions" ]]
  colnames( df ) <- gsub( pattern = "glia_clusters", replacement = "seurat_clusters", x = colnames( df ) )  
  
  df <- df[which(df$seurat_clusters %in% trajectory_list[[ cell_type ]]),]
  
  df <- df %>%
    mutate(
      seurat_clusters = as.factor(seurat_clusters),
      genotype = as.factor(genotype)
    )
  
  wt_means <- df %>%
    filter(genotype == "WT") %>%
    group_by(seurat_clusters) %>%
    summarise(WT_mean_fraction = mean(fraction, na.rm = TRUE))
  
  df_rel <- df %>%
    left_join(wt_means, by = "seurat_clusters") %>%
    mutate(relative_fraction = log2(fraction / WT_mean_fraction) )
  
  summary_rel <- df_rel %>%
    filter(genotype != "WT") %>%         # we normalize to WT
    group_by(genotype, seurat_clusters) %>%
    summarise(
      mean_rel = mean(relative_fraction, na.rm = TRUE),
      sd_rel   = sd(relative_fraction, na.rm = TRUE),
      .groups = "drop"
    )
  
  # order the clusters based on the trajectory information
  summary_rel$seurat_clusters <- factor( summary_rel$seurat_clusters, levels = trajectory_list[[ cell_type ]])
  
  # do the plotting
  traj_abundance_plot[[ cell_type ]] <- ggplot(summary_rel, aes(x = seurat_clusters, y = mean_rel,
                                                                group = genotype, color = genotype, fill = genotype)) +
    geom_errorbar(aes(ymin = mean_rel - sd_rel, ymax = mean_rel + sd_rel)) +
    geom_line(size = 1.2) +
    geom_point(size = 2) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
    labs(
      title = paste( "log2 Relative cluster fractions", cell_type),
      x = "Seurat cluster",
      y = "Relative fraction (mean ± SD)"
    ) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
}

# save the plot
# NOTE: statistics for this plot come from Rev_Fig3B!
plot_grid( plotlist = traj_abundance_plot, ncol = 3)
ggsave( filename = paste(base, "/plots/reviewer/Rev_Fig3C_glia_cluster_traj_abundances.pdf", sep=""), width=15, height=4 )

# =============================================
# DEG analysis
# with MAST in preparation of iDEA analysis
# =============================================

de_res <- list()
cell_count <- list()

comp <- list()
comp[[ 1 ]] <- c("KO", "WT")
comp[[ 2 ]] <- c("DKO", "WT")
comp[[ 3 ]] <- c("DKO", "KO")

Idents( seurat.glia.obj ) <- 'genotype'

for ( idx in c(1:3) ) {
  
  # give a descriptor to the analysis
  del_gen  <- paste(comp [[ idx ]], collapse = "_")
  
  de_res[[del_gen]] <- list()
  cell_count[[del_gen]] <- list()
  
  test_diff <- subset(x = seurat.glia.obj, idents = c( comp[[ idx ]][1], comp[[ idx ]][2]), invert = FALSE)
  
  for ( i in unique( seurat.glia.obj$cell_type ) ) {
    
    message( del_gen, " ", i )
    
    Idents(test_diff) <- "cell_type"
    test_diff_cluster <- subset(x = test_diff, idents = i, invert = FALSE)
    Idents(test_diff_cluster) <- 'genotype'
    
    #testing for DE genes
    cell_count[[del_gen]][[i]] <- table(Idents(test_diff_cluster))
    
    #remove clusters with less than 9 cells
    if(cell_count[[del_gen]][[i]][1] > 9 & cell_count[[del_gen]][[i]][2] > 9) {
      
      de_res[[del_gen]][[i]] <- FindMarkers(test_diff_cluster, ident.1 = comp[[ idx ]][1], ident.2 = comp[[ idx ]][2], test.use = 'MAST', logfc.threshold = 0.01, assay = "RNA", slot = "data")
      
    } else {
      de_res[[del_gen]][[i]] <- NA
    }
  }
}

#write out DEGs - as Supplement
wb <- createWorkbook()

for (del_gen in names(de_res)) {
  for (comp in names(de_res[[del_gen]])) {
    if (!any(is.na(de_res[[del_gen]][[comp]]))) {
      
      tmp <- de_res[[del_gen]][[comp]]
      tmp$gene <- rownames(tmp)
      
      # modify cluster names to be compatible with EXEL tab naming conventions
      sheet_name_comp <- strsplit(x = comp, split = "\\ ")[[1]]
      sheet_name_comp <- sapply(sheet_name_comp, function (x) gsub(pattern = "\\(|\\)|\\[|\\]|-", replacement = "", x = x))
      sheet_name_comp <- paste(sheet_name_comp, collapse="_")
      
      sheetName <- paste(del_gen, sheet_name_comp, sep='_')
      addWorksheet(wb, sheetName)
      writeData(wb, sheetName, tmp)
      addFilter(wb, sheetName, row = 1, cols = 1:ncol(tmp))
      setColWidths(wb, sheetName, cols = 1:ncol(tmp), widths="auto")  
    }
  }
}

saveWorkbook(wb, file = paste(base, "/Supplement/DEG.analysis_MAST.xlsx", sep="/"), overwrite = T) 

# prepare a list suitable for iDEA

de_res <- list()

comparisons <- getSheetNames( paste(base, "/Supplement/DEG.analysis_MAST.xlsx", sep="") )
for (comp in comparisons) {
  comp_parts <- strsplit(x = comp, split = "_")[[1]]
  if (length(comp_parts) == 3) {
    comp_parts[1] <- paste(comp_parts[1], comp_parts[2], sep='_')
    comp_parts[2] <- comp_parts[3]
  } 
  
  if (length(de_res[[comp_parts[1]]]) == 0) {
    de_res[[comp_parts[1]]] <- list()
  }
  de_res[[comp_parts[1]]][[comp_parts[2]]] <- read.xlsx( xlsxFile = paste(base, "/Supplement/DEG.analysis_MAST.xlsx", sep=""), sheet = comp)
}

# this file is not included in the GitHub repo
saveRDS(object = de_res, file = paste(base, "/iDEA_analysis/DEG_iDEA.RDS", sep = ""))

sessionInfo()
# R version 4.3.2 (2023-10-31)
# Platform: x86_64-pc-linux-gnu (64-bit)
# Running under: Ubuntu 22.04.5 LTS
# 
# Matrix products: default
# BLAS:   /opt/R/4.3.2/lib/R/lib/libRblas.so 
# LAPACK: /usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3.10.0
# 
# locale:
#   [1] LC_CTYPE=en_GB.UTF-8       LC_NUMERIC=C               LC_TIME=de_AT.UTF-8        LC_COLLATE=en_GB.UTF-8    
# [5] LC_MONETARY=de_AT.UTF-8    LC_MESSAGES=en_GB.UTF-8    LC_PAPER=de_AT.UTF-8       LC_NAME=C                 
# [9] LC_ADDRESS=C               LC_TELEPHONE=C             LC_MEASUREMENT=de_AT.UTF-8 LC_IDENTIFICATION=C       
# 
# time zone: Europe/Vienna
# tzcode source: system (glibc)
# 
# attached base packages:
#   [1] stats4    stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] ggbeeswarm_0.7.2            ggraph_2.2.1                igraph_2.1.4                purrr_1.0.2                
# [5] slingshot_2.10.0            TrajectoryUtils_1.10.1      SingleCellExperiment_1.24.0 SummarizedExperiment_1.32.0
# [9] GenomicRanges_1.54.1        GenomeInfoDb_1.38.8         MatrixGenerics_1.14.0       matrixStats_1.2.0          
# [13] princurve_2.1.6             pheatmap_1.0.12             speckle_1.2.0               cowplot_1.1.3              
# [17] RColorBrewer_1.1-3          openxlsx_4.2.5.2            SeuratWrappers_0.3.5        org.Mm.eg.db_3.18.0        
# [21] AnnotationDbi_1.64.1        IRanges_2.36.0              S4Vectors_0.40.2            Biobase_2.62.0             
# [25] BiocGenerics_0.48.1         clusterProfiler_4.10.1      data.table_1.15.0           dplyr_1.1.4                
# [29] ggplot2_3.5.2               Seurat_5.0.1                SeuratObject_5.0.1          sp_2.1-3                   
# 
# loaded via a namespace (and not attached):
#   [1] fs_1.6.3                  spatstat.sparse_3.0-3     bitops_1.0-7              enrichplot_1.22.0        
# [5] HDO.db_0.99.1             httr_1.4.7                tools_4.3.2               sctransform_0.4.1        
# [9] utf8_1.2.4                R6_2.5.1                  lazyeval_0.2.2            uwot_0.1.16              
# [13] withr_3.0.0               prettyunits_1.2.0         gridExtra_2.3             progressr_0.14.0         
# [17] textshaping_0.3.7         cli_3.6.2                 spatstat.explore_3.2-6    fastDummies_1.7.3        
# [21] scatterpie_0.2.2          labeling_0.4.3            spatstat.data_3.0-4       ggridges_0.5.6           
# [25] pbapply_1.7-2             systemfonts_1.0.6         yulab.utils_0.1.4         gson_0.1.0               
# [29] DOSE_3.28.2               R.utils_2.12.3            harmony_1.2.0             parallelly_1.37.0        
# [33] limma_3.58.1              rstudioapi_0.16.0         RSQLite_2.3.6             generics_0.1.3           
# [37] gridGraphics_0.5-1        ica_1.0-3                 spatstat.random_3.2-2     zip_2.3.1                
# [41] GO.db_3.18.0              Matrix_1.6-5              fansi_1.0.6               abind_1.4-5              
# [45] R.methodsS3_1.8.2         lifecycle_1.0.4           edgeR_4.0.16              qvalue_2.34.0            
# [49] SparseArray_1.2.4         Rtsne_0.17                grid_4.3.2                blob_1.2.4               
# [53] promises_1.2.1            crayon_1.5.2              miniUI_0.1.1.1            lattice_0.21-9           
# [57] KEGGREST_1.42.0           pillar_1.9.0              fgsea_1.28.0              future.apply_1.11.1      
# [61] codetools_0.2-19          fastmatch_1.1-4           leiden_0.4.3.1            glue_1.7.0               
# [65] packrat_0.9.2             ggfun_0.1.4               remotes_2.5.0             vctrs_0.6.5              
# [69] png_0.1-8                 treeio_1.26.0             spam_2.10-0               gtable_0.3.6             
# [73] cachem_1.0.8              S4Arrays_1.2.1            mime_0.12                 tidygraph_1.3.1          
# [77] survival_3.5-7            statmod_1.5.0             ellipsis_0.3.2            fitdistrplus_1.1-11      
# [81] ROCR_1.0-11               nlme_3.1-163              ggtree_3.10.1             bit64_4.0.5              
# [85] progress_1.2.3            RcppAnnoy_0.0.22          irlba_2.3.5.1             vipor_0.4.7              
# [89] KernSmooth_2.23-22        colorspace_2.1-0          DBI_1.2.2                 ggrastr_1.0.2            
# [93] tidyselect_1.2.1          bit_4.0.5                 compiler_4.3.2            DelayedArray_0.28.0      
# [97] plotly_4.10.4             shadowtext_0.1.3          scales_1.4.0              lmtest_0.9-40            
# [101] stringr_1.5.1             digest_0.6.34             goftest_1.2-3             presto_1.0.0             
# [105] spatstat.utils_3.1-0      RhpcBLASctl_0.23-42       XVector_0.42.0            htmltools_0.5.7          
# [109] pkgconfig_2.0.3           sparseMatrixStats_1.14.0  fastmap_1.1.1             rlang_1.1.3              
# [113] htmlwidgets_1.6.4         DelayedMatrixStats_1.24.0 shiny_1.8.0               farver_2.1.1             
# [117] zoo_1.8-12                jsonlite_1.8.8            BiocParallel_1.36.0       GOSemSim_2.28.1          
# [121] R.oo_1.26.0               RCurl_1.98-1.14           magrittr_2.0.3            GenomeInfoDbData_1.2.11  
# [125] ggplotify_0.1.2           dotCall64_1.1-1           patchwork_1.2.0           Rcpp_1.0.12              
# [129] ape_5.8                   viridis_0.6.5             reticulate_1.44.1         stringi_1.8.3            
# [133] zlibbioc_1.48.2           MASS_7.3-60               MAST_1.28.0               plyr_1.8.9               
# [137] parallel_4.3.2            listenv_0.9.1             ggrepel_0.9.5             deldir_2.0-2             
# [141] Biostrings_2.70.3         graphlayouts_1.1.1        splines_4.3.2             tensor_1.5               
# [145] hms_1.1.3                 locfit_1.5-9.12           spatstat.geom_3.2-8       RcppHNSW_0.6.0           
# [149] reshape2_1.4.4            BiocManager_1.30.22       tweenr_2.0.3              httpuv_1.6.14            
# [153] RANN_2.6.1                tidyr_1.3.1               polyclip_1.10-6           future_1.33.1            
# [157] scattermore_1.2           ggforce_0.4.2             rsvd_1.0.5                xtable_1.8-4             
# [161] RSpectra_0.16-1           tidytree_0.4.6            later_1.3.2               ragg_1.3.0               
# [165] viridisLite_0.4.2         tibble_3.2.1              aplot_0.2.2               memoise_2.0.1            
# [169] beeswarm_0.4.0            cluster_2.1.4             globals_0.16.2    
