# ===============================
# second step of the analysis
# analysis of Pten PIP-Seq data
# ===============================

library (Seurat)
library (ggplot2)
library (dplyr)
library (data.table)
library (clusterProfiler)
library (org.Mm.eg.db)
library (monocle3)
library (SeuratWrappers)
library (openxlsx)
library (ggbeeswarm)
library (RColorBrewer)
library (cowplot)

base <- "./scRNA-Seq/biorxiv"
no_save_base <- "./Miranda_et_al_noSave" # folder for large files that are not necessary to backup
count_base <- paste( base, "PIPSeeker_output/", sep="" ) # folder where PIPSeeker output is located

set.seed(2401)

# =========================================================
# convert sample idx to genotype
# KO: sparse KO, i.e. green cells in a Pten-MADM
# DKO: double KO, sparse Pten ko in a Egfr cKO background
# =========================================================
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
  # PIPSeeker pipeline produces a number of outputs  -  only one used here
  file_path <- paste( count_base, "sample", sample_idx, "/filtered_matrix/sensitivity_3", sep="" )
  message( file_path )  
  
  counts <- Read10X( data.dir = file_path,  )
  seurat_obj_list[[ sample_idx ]] <- CreateSeuratObject(counts = counts, project = conv_vec[ sample_idx ], min.cells = 3, min.features = 200)
}

# merging objects
seurat.obj <- merge( x = seurat_obj_list[[ 1 ]], y = seurat_obj_list[ c(2:8) ] )

message("preparing QC plots")
# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
seurat.obj[["percent.mt"]] <- PercentageFeatureSet(seurat.obj, pattern = "^mt-")

# ==========
# QC plots
# ==========

vln_plot <- VlnPlot(seurat.obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
ggsave( plot=vln_plot, file=paste (base, "/plots/QC/basicFeatures.Vln.raw.pdf", sep=""), width = 10, height = 7)

message("filtering")
seurat.obj <- subset(seurat.obj, subset = nFeature_RNA > 500 & nFeature_RNA < 6000 & nCount_RNA < 40000 & percent.mt < 15)
vln_plot <- VlnPlot(seurat.obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

ggsave( plot=vln_plot, file=paste (base, "/plots/QC/basicFeatures.Vln.filtered.pdf", sep=""), width = 12, height = 5)

plot1 <- FeatureScatter(seurat.obj, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(seurat.obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
QC_plot <- plot1 + plot2

ggsave( plot=QC_plot, file=paste (base, "/plots/QC/combined.QC.pdf", sep=""), width = 10, height = 7)

seurat.obj <- NormalizeData(seurat.obj, normalization.method = "LogNormalize", scale.factor = 10000)
seurat.obj <- FindVariableFeatures(seurat.obj, selection.method = "vst", nfeatures = 2000)

seurat.obj <- ScaleData(seurat.obj)
seurat.obj <- RunPCA(seurat.obj, features = VariableFeatures(object = seurat.obj))
elbowPlot <- ElbowPlot( object = seurat.obj, ndims = 30)
seurat.obj <- FindNeighbors(seurat.obj, dims = 1:25)
seurat.obj <- FindClusters(seurat.obj, resolution = 0.3)

seurat.obj <- RunUMAP(seurat.obj, dims = 1:25)

# for internal QC: no batch effects visible
umap_replicates <- DimPlot( seurat.obj, group.by = c("orig.ident"), label=F )
ggsave( filename = paste(base, "/plots/QC/replicate_umap.pdf", sep=""), plot = umap_replicates, width=5.5, height=5 )

# for documentation purposes - one cluster will be removed, visualized here 
umap_clusters <- DimPlot( seurat.obj, group.by = c("seurat_clusters"), label=T) + NoLegend()
ggsave( filename = paste(base, "/plots/QC/cluster_umap.pdf", sep=""), plot = umap_clusters, width=5, height=5 )

# ====================
# cluster removal
# ====================

seurat.obj <- JoinLayers( seurat.obj )
# create a variable that is ordered in a convenient way
seurat.obj@meta.data$plot_id <- seurat.obj@meta.data$orig.ident
seurat.obj@meta.data$plot_id <- factor( seurat.obj@meta.data$plot_id, levels = c("WT1", "WT2", "KO1", "KO2", "KO3", "KO4", "DKO1", "DKO2") )

# identify markers
# internal QC to rationaize removal of Cluster7
cluster_markers <- FindAllMarkers( seurat.obj, only.pos = T)

cluster_markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10

heat_plot <- DoHeatmap( seurat.obj, features = top10$gene ) + NoLegend()
ggsave( filename = paste(base, "/plots/QC/marker_heatmap.pdf", sep=""), plot = heat_plot, width=8, height=9 )

# cluster 7 is not well localized and shows very few specific markers
# remove from analysis
seurat.obj <- subset(seurat.obj, seurat_clusters %in% c(7), invert=T)

# ====================
# QC for Supplement
# ====================

# add a factor for nice sample order
seurat.obj@meta.data$sampleID <- factor( seurat.obj@meta.data$orig.ident, levels = c("WT1", "WT2", "KO1", "KO2", "KO3", "KO4", "DKO1", "DKO2") )

# UMAP for supplement with replicates colored
umap_replicates <- DimPlot( seurat.obj, group.by = c("sampleID"), label=F)
ggsave( filename = paste(base, "/plots/Supplement/Sup8E_full_umap_filtered_replicates.pdf", sep=""), plot = umap_replicates, width=6, height=5 )

# create a genotype variable
seurat.obj@meta.data$genotype <- seurat.obj@meta.data$orig.ident
seurat.obj@meta.data$genotype <- sapply( seurat.obj@meta.data$genotype, function (x) gsub(pattern = "1|2|3|4", replacement = "", x = x))

# create a color code for genotypes
color_palette <- c("WT" = "grey50", "DKO" = "cornflowerblue", "KO" = "white" )

colors <- sapply( unique(seurat.obj@meta.data$genotype), function (x) {
  n_col <- length( unique( seurat.obj@meta.data$orig.ident[ which(seurat.obj@meta.data$genotype == x) ] ) )
  color_list <- rep( color_palette[[x]], n_col ) 
  
  return(color_list)
})

# basic QC
vln_plot <- VlnPlot(seurat.obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, group.by = "sampleID")
vln_plot <- vln_plot & scale_fill_manual( values = unlist(colors) )

ggsave( plot=vln_plot, file=paste (base, "/plots/Supplement/Sup8B-D_basicFeatures.Vln.filtered.pdf", sep=""), width = 12, height = 5)

# marker expression on UMAP
marker_plot <- FeaturePlot( seurat.obj, features = c("Dcx", "Neurod1", "Gad2", "Fabp7"), order = T, min.cutoff = "q10" ) & NoLegend()
ggsave( filename = paste(base, "/plots/Supplement/Sup8F_marker_expression_plots.pdf", sep=""), plot = marker_plot, width=9, height=9 )

# give useful cell type names
cluster2celltype <- c("iN", "Neurons", "aIP", "cProg", "oligo", "OBNB", "aIP")
names( cluster2celltype ) <- c("0", "1", "2", "3", "4", "5", "6")

seurat.obj@meta.data$cell_type <- cluster2celltype[ as.character(seurat.obj@meta.data$seurat_clusters) ]

# plot the cell number per sample
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
  
ggsave( filename = paste(base, "/plots/Supplement/Sup8A_cell_number_per_sample.pdf", sep=""), plot = cell_number_plot, width=3, height=5 )

# write out the raw data for this plot
wb <- createWorkbook()

tmp <- df2plot

sheetName <- "raw data point Fig5K"
addWorksheet(wb, sheetName)
writeData(wb, sheetName, tmp)
addFilter(wb, sheetName, row = 1, cols = 1:ncol(tmp))
setColWidths(wb, sheetName, cols = 1:ncol(tmp), widths="auto")  

tmp <- av_group

sheetName <- "averages Fig5K"
addWorksheet(wb, sheetName)
writeData(wb, sheetName, tmp)
addFilter(wb, sheetName, row = 1, cols = 1:ncol(tmp))
setColWidths(wb, sheetName, cols = 1:ncol(tmp), widths="auto")  

saveWorkbook(wb, file = paste(base, "Supplement/Sup8A_data.xlsx", sep="/"), overwrite = T) 

# ==========================================================
# cell type analysis
# ==========================================================

seurat.obj@meta.data$super_group <- "proj_neuron" 
seurat.obj@meta.data$super_group[ which( !seurat.obj@meta.data$cell_type %in% c("Neurons", "iN")) ] <- "non_proj_neuron"


plot2 <- DimPlot( seurat.obj, group.by = "super_group", label=T ) + NoLegend()

ggsave( filename = paste(base, "/plots/Supplement/Sup8G_broad_cluster_umap.pdf", sep=""), plot = plot2, width=6, height=5 )


# ================================
# focus analysis on glia/OBNB
# ================================

Idents(seurat.obj) <- "cell_type"
seurat_glia <- subset( seurat.obj, idents = c("aIP", "cProg", "OBNB", "oligo")) 

# remove all old analyses
seurat_glia@assays$RNA$data <- NULL
seurat_glia@assays$RNA$scale.data <- NULL
seurat_glia@reductions$pca <- NULL
seurat_glia@reductions$umap <- NULL

seurat_glia <- NormalizeData(seurat_glia, normalization.method = "LogNormalize", scale.factor = 10000)
seurat_glia <- FindVariableFeatures(seurat_glia, selection.method = "vst", nfeatures = 4000)

# partially regress out cell cycle effects
cell_cycle_phase_gene <- read.table(file = paste(base, "other_data/Mus_musculus.csv", sep=""), sep = ",", header = T)
gene_conversion <- bitr(geneID = cell_cycle_phase_gene$geneID, fromType = "ENSEMBL", toType = "SYMBOL", OrgDb = org.Mm.eg.db)
cell_cycle_phase_gene <- merge(cell_cycle_phase_gene, gene_conversion, by.x="geneID", by.y="ENSEMBL")

seurat_glia <- CellCycleScoring (
  object = seurat_glia,
  g2m.features = cell_cycle_phase_gene[which(cell_cycle_phase_gene$phase == "G2/M"), "SYMBOL"],
  s.features = cell_cycle_phase_gene[which(cell_cycle_phase_gene$phase == "S"), "SYMBOL"]
)

seurat_glia$CC.Difference <- seurat_glia$S.Score - seurat_glia$G2M.Score

all.genes <- rownames(seurat_glia)
seurat_glia <- ScaleData(seurat_glia, vars.to.regress = c("CC.Difference") )

seurat_glia <- RunPCA(seurat_glia, features = VariableFeatures(object = seurat_glia))

elbowPlot <- ElbowPlot( object = seurat_glia, ndims = 30)
seurat_glia <- FindNeighbors(seurat_glia, dims = 1:15)
seurat_glia <- FindClusters(seurat_glia, resolution = 0.5)

seurat_glia <- RunUMAP(seurat_glia, dims = 1:15)

seurat_glia@meta.data$gen2plot <- factor( seurat_glia@meta.data$genotype, levels = c("WT", "KO", "DKO") ) 

umap_df <- as.data.frame( Embeddings( object = seurat_glia, reduction = "umap") )
umap_df$genotype <- seurat_glia$genotype
umap_df$genotype <- factor( umap_df$genotype, levels = c("WT", "KO", "DKO") )
ggplot( ) + 
  geom_point(data = umap_df[ which(umap_df$genotype == "KO"), ], aes(x=umap_1, y=umap_2, fill = genotype), shape=21) +
  geom_point(data = umap_df[ which(umap_df$genotype == "DKO"), ], aes(x=umap_1, y=umap_2, fill = genotype), shape=21) +
  geom_point(data = umap_df[ which(umap_df$genotype == "WT"), ], aes(x=umap_1, y=umap_2, fill = genotype), shape=21) +
  scale_fill_manual(  values = c("WT" = "grey10", "DKO" = "cornflowerblue", "KO" = "white") ) +
  theme_classic()
ggsave( filename = paste(base, "/plots/QC/glia_umap_genotype.pdf", sep=""), width=6, height=5 )

FeaturePlot( seurat_glia, features = c("Cdk1", "Mki67", "Gad2", "Dlx5", "Pdgfra", "Cspg4", "Aldh1l1", "Aldoc", "Gfap"), order = T, min.cutoff = "q15" ) & NoLegend()
ggsave( filename = paste(base, "/plots/Supplement/Sup8H_glia_umap_marker_expression.pdf", sep=""), width=9, height=9 )

###
# refine cell type annotation
###

cluster2cellType <- c("cProg", "aIP", "iA", "iA", "iA", "iA", "oligo", "oligo", "OBNB")
names(cluster2cellType) <- as.character( c(0, 1, 2, 5, 7, 8, 4, 6, 3) )
seurat_glia@meta.data$cell_type <- cluster2cellType[ as.character(seurat_glia$seurat_clusters) ]

# sanity check for cell type assignment
VlnPlot( seurat_glia, features = c( "Mki67",  "Aldh1l1", "Gad2", "Pdgfra" ), group.by = "cell_type", ncol = 4)
ggsave( filename = paste(base, "/plots/QC/marker_expression_glia.png", sep=""), width=8, height=3 )

# figure plots
DimPlot( seurat_glia, group.by = c("cell_type", "Phase"), split.by = "gen2plot", label=F, label.size = 6)
ggsave( filename = paste(base, "/plots/Supplement/Sup8I_glia_umap_cellTypes.pdf", sep=""), width=9, height=8 )

DimPlot( seurat_glia, group.by = "cell_type", label=F)
ggsave( filename = paste(base, "/plots/main/Fig5B_umap.pdf", sep=""), width=6, height=5 )

# ======================
# trajectory analysis
# ======================

cds <- as.cell_data_set( seurat_glia )
cds <- cluster_cells(cds)
cds <- learn_graph(cds)
cds <- order_cells(cds)

cds@colData$cell_type <- cluster2cellType[ as.character(cds@colData$seurat_clusters) ]

traj_plot_ct <- list()
for (gen in c("WT", "KO", "DKO")) {

  traj_plot_ct[[ gen ]] <- plot_cells(cds[,rownames(cds@colData)[ which( cds@colData$genotype == gen ) ]], label_groups_by_cluster = FALSE, label_cell_groups = FALSE, cell_size = 1.5,
                             label_leaves = FALSE, label_branch_points = FALSE, label_roots = FALSE, 
                             color_cells_by = "cell_type", group_label_size = 5, trajectory_graph_segment_size = 3)
  traj_plot_ct[[ gen ]] <- traj_plot_ct[[ gen ]] + ggtitle( gen )
  
}

traj_plot <- plot_grid( plotlist = traj_plot_ct, ncol = 3 )

ggsave  ( traj_plot, filename = paste(base, "plots/QC/glia_traj_plot_genotypes.pdf", sep=""), width=15, height=5)

traj_plot_main <- plot_cells(cds, label_groups_by_cluster = FALSE, label_cell_groups = FALSE, cell_size = 1.5,
           label_leaves = FALSE, label_branch_points = FALSE, label_roots = FALSE, 
           color_cells_by = "cell_type", group_label_size = 5, trajectory_graph_segment_size = 3)
ggsave( plot = traj_plot_main, filename = paste(base, "/plots/main/Fig5C_traj_plot.pdf", sep=""), width=6, height=5 )

# analyse abundance changes of non projection neuron cell types
seurat_glia@meta.data %>%
  group_by( orig.ident, cell_type ) %>%
  summarize( n = n()) -> cluster_count

cluster_count$genotype <- sapply( cluster_count$orig.ident, function (x) gsub(pattern = "1|2|3|4", replacement = "", x = x))
cluster_count$total <- total_cells[ cluster_count$orig.ident ]
cluster_count$rel <- cluster_count$n / cluster_count$total

cluster_count %>%
  group_by( genotype, cell_type ) %>%
  summarize( mean = median(rel), sd = sd(rel)) -> cluster_count_comb

plot1 <- ggplot() + 
  geom_bar( data = cluster_count_comb, aes( x=genotype, y=mean ), stat="identity", position="dodge", fill="white", color="black") + 
  #geom_errorbar( data = cluster_count_comb, aes(x=genotype, ymin = mean-sd, ymax=mean+sd)) + 
  geom_beeswarm( data = cluster_count, aes( x=genotype, y=rel ) ) +
  facet_wrap( ~ cell_type, scales = "free_y" ) + theme_classic()

ggsave  ( filename = paste(base, "plots/QC/glia_cellTypes_abundance.pdf", sep=""), plot = plot1, width=6, height=4)

# arcsin transformation
cluster_count$asin_rel <- asin(sqrt( cluster_count$rel ))
anova <- aov( asin_rel ~ cell_type*genotype, data = cluster_count)
summary(anova)

#                    Df  Sum Sq Mean Sq F value   Pr(>F)    
# cell_type           4 0.07032 0.01758   5.685 0.002140 ** 
# genotype            2 0.07629 0.03814  12.334 0.000188 ***
# cell_type:genotype  8 0.04767 0.00596   1.927 0.100383    
# Residuals          25 0.07731 0.00309                     
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# do a chisquare analysis of pooled relative abundances
# due to the variation in the cell number direct comparisons might be less powerful
# much more sensitive than the anova with multiple testing
# do all pairwise comparisons

pvalues <- data.frame()

for (cell_type in unique(cluster_count$cell_type)) {
  tmp <- cluster_count[which(cluster_count$cell_type == cell_type),]
  
  tmp %>%
    group_by( genotype ) %>%
    summarize( sum_n = sum( n ), sum_total = sum( total ) ) -> tmp_comb
  
  ctrl_idx <- which(tmp_comb$genotype == "WT")
  
  for (genotype in c("KO", "DKO")) {
      idx <- which(tmp_comb$genotype == genotype)
      tmp_table <- as.table( rbind( c(as.numeric ( tmp_comb[ctrl_idx, "sum_n"] ) , as.numeric ( tmp_comb[ctrl_idx, "sum_total"] - tmp_comb[ctrl_idx, "sum_n"]) ),
                                    c(as.numeric ( tmp_comb[idx, "sum_n"] ), as.numeric ( tmp_comb[idx, "sum_total"] - tmp_comb[idx, "sum_n"] ) ) )
      )
      chisq_out <- chisq.test( tmp_table, correct = T, rescale.p = F, simulate.p.value = F )
      
      tmp_comb$rel <- tmp_comb$sum_n / tmp_comb$sum_total 
      tmp_comb <- as.data.frame( tmp_comb ) 
      rownames(tmp_comb) <- tmp_comb$genotype
      
      dir_change <- tmp_comb[ genotype, "rel"] / tmp_comb[ "WT", "rel"]
      
      if (dir_change > 1) {
        dir_change <- "pos" 
      } else {
        dir_change <- "neg"
      }
      pvalues <- rbind( pvalues, data.frame(cluster = cell_type, genotype = genotype, 
                                            abs_genotype = tmp_comb[ genotype, "sum_n"], rel_genotype = tmp_comb[ genotype, "rel"], 
                                            abs_wt = tmp_comb[ genotype, "sum_total"], rel_wt = tmp_comb[ "WT", "rel"], 
                                            pvalue = chisq_out$p.value, dir = dir_change ) )
    }
}

pvalues$p.adjust <- p.adjust( pvalues$pvalue, method = "bonferroni" )
pvalues$sig_stars <- sapply( pvalues$p.adjust, function (x) ifelse(x<0.05, "*", ""))

pvalues$score <- log10(pvalues$p.adjust) * -1
pvalues$score[ which(pvalues$score > 5 )] <- 5 
pvalues$score[ which( pvalues$dir == "neg") ] <- pvalues$score[ which( pvalues$dir == "neg") ] * -1

# main figure heatmap
ggplot( pvalues, aes(x=genotype, y=cluster, fill=dir, label=sig_stars)) + geom_tile( color = "black" ) + 
  geom_text( size = 15) + scale_fill_manual( values = c( "neg" = "cornflowerblue", "pos" = "yellow" ) ) + 
  theme_classic() + ylab( "cell type" )

ggsave  ( filename = paste(base, "/plots/main/Fig5D_abundance_change_p_plot.pdf", sep=""), width=4, height=4)

# save analyses to sup exel sheet
# save this analysis as exel
wb <- createWorkbook()

tmp <- pvalues

sheetName <- "raw data Figure 5D"
addWorksheet(wb, sheetName)
writeData(wb, sheetName, tmp)
addFilter(wb, sheetName, row = 1, cols = 1:ncol(tmp))
setColWidths(wb, sheetName, cols = 1:ncol(tmp), widths="auto")  

saveWorkbook(wb, file = paste(base, "Supplement/Fig5D_rel_abundance_analysis.xlsx", sep="/"), overwrite = T) 

# ===========================================
# DEG analysis
# with MAST in preparation of iDEA analysis
# ===========================================

de_res <- list()
cell_count <- list()

comp <- list()
comp[[ 1 ]] <- c("KO", "WT")
comp[[ 2 ]] <- c("DKO", "WT")
comp[[ 3 ]] <- c("DKO", "KO")

Idents( seurat_glia ) <- 'genotype'

for ( idx in c(1:3) ) {
  
  # give a descriptor to the analysis
  del_gen  <- paste(comp [[ idx ]], collapse = "_")
  
  de_res[[del_gen]] <- list()
  cell_count[[del_gen]] <- list()
  
  test_diff <- subset(x = seurat_glia, idents = c( comp[[ idx ]][1], comp[[ idx ]][2]), invert = FALSE)
  
  for ( i in unique( seurat_glia$cell_type ) ) {
    
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

saveWorkbook(wb, file = paste(base, "Supplement/DEG.analysis_MAST.xlsx", sep="/"), overwrite = T) 

# prepare a list suitable for iDEA

de_res <- list()

comparisons <- getSheetNames( paste(base, "Supplement/DEG.analysis_MAST.xlsx", sep="/") )
for (comp in comparisons) {
  comp_parts <- strsplit(x = comp, split = "_")[[1]]
  if (length(comp_parts) == 3) {
    comp_parts[1] <- paste(comp_parts[1], comp_parts[2], sep='_')
    comp_parts[2] <- comp_parts[3]
  } 
  
  if (length(de_res[[comp_parts[1]]]) == 0) {
    de_res[[comp_parts[1]]] <- list()
  }
  de_res[[comp_parts[1]]][[comp_parts[2]]] <- read.xlsx( xlsxFile = paste(base, "Supplement/DEG.analysis_MAST.xlsx", sep="/"), sheet = comp)
}

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
#   [1] LC_CTYPE=en_GB.UTF-8       LC_NUMERIC=C               LC_TIME=de_AT.UTF-8        LC_COLLATE=en_GB.UTF-8     LC_MONETARY=de_AT.UTF-8    LC_MESSAGES=en_GB.UTF-8   
# [7] LC_PAPER=de_AT.UTF-8       LC_NAME=C                  LC_ADDRESS=C               LC_TELEPHONE=C             LC_MEASUREMENT=de_AT.UTF-8 LC_IDENTIFICATION=C       
# 
# time zone: Europe/Vienna
# tzcode source: system (glibc)
# 
# attached base packages:
#   [1] stats4    stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] shiny_1.8.0                 openxlsx_4.2.5.2            SeuratWrappers_0.3.5        monocle3_1.3.7              SingleCellExperiment_1.24.0 SummarizedExperiment_1.32.0
# [7] GenomicRanges_1.54.1        GenomeInfoDb_1.38.8         MatrixGenerics_1.14.0       matrixStats_1.2.0           org.Mm.eg.db_3.18.0         AnnotationDbi_1.64.1       
# [13] IRanges_2.36.0              S4Vectors_0.40.2            Biobase_2.62.0              BiocGenerics_0.48.1         clusterProfiler_4.10.1      data.table_1.15.0          
# [19] dplyr_1.1.4                 ggplot2_3.5.0               Seurat_5.0.1                SeuratObject_5.0.1          sp_2.1-3                   
# 
# loaded via a namespace (and not attached):
#   [1] fs_1.6.3                spatstat.sparse_3.0-3   bitops_1.0-7            enrichplot_1.22.0       HDO.db_0.99.1           httr_1.4.7              RColorBrewer_1.1-3     
# [8] tools_4.3.2             sctransform_0.4.1       utf8_1.2.4              R6_2.5.1                lazyeval_0.2.2          uwot_0.1.16             withr_3.0.0            
# [15] prettyunits_1.2.0       gridExtra_2.3           progressr_0.14.0        textshaping_0.3.7       cli_3.6.2               spatstat.explore_3.2-6  fastDummies_1.7.3      
# [22] scatterpie_0.2.2        sass_0.4.8              labeling_0.4.3          spatstat.data_3.0-4     proxy_0.4-27            ggridges_0.5.6          pbapply_1.7-2          
# [29] systemfonts_1.0.6       yulab.utils_0.1.4       gson_0.1.0              DOSE_3.28.2             R.utils_2.12.3          parallelly_1.37.0       limma_3.58.1           
# [36] rstudioapi_0.16.0       RSQLite_2.3.6           generics_0.1.3          gridGraphics_0.5-1      ica_1.0-3               spatstat.random_3.2-2   zip_2.3.1              
# [43] GO.db_3.18.0            Matrix_1.6-5            ggbeeswarm_0.7.2        fansi_1.0.6             abind_1.4-5             R.methodsS3_1.8.2       lifecycle_1.0.4        
# [50] qvalue_2.34.0           SparseArray_1.2.4       Rtsne_0.17              grid_4.3.2              blob_1.2.4              promises_1.2.1          crayon_1.5.2           
# [57] miniUI_0.1.1.1          lattice_0.21-9          cowplot_1.1.3           KEGGREST_1.42.0         pillar_1.9.0            fgsea_1.28.0            boot_1.3-28.1          
# [64] future.apply_1.11.1     codetools_0.2-19        fastmatch_1.1-4         leiden_0.4.3.1          glue_1.7.0              packrat_0.9.2           leidenbase_0.1.27      
# [71] ggfun_0.1.4             remotes_2.5.0           vctrs_0.6.5             png_0.1-8               treeio_1.26.0           spam_2.10-0             gtable_0.3.4           
# [78] assertthat_0.2.1        cachem_1.0.8            S4Arrays_1.2.1          mime_0.12               tidygraph_1.3.1         survival_3.5-7          statmod_1.5.0          
# [85] ellipsis_0.3.2          fitdistrplus_1.1-11     ROCR_1.0-11             nlme_3.1-163            ggtree_3.10.1           bit64_4.0.5             progress_1.2.3         
# [92] RcppAnnoy_0.0.22        bslib_0.6.1             irlba_2.3.5.1           vipor_0.4.7             KernSmooth_2.23-22      colorspace_2.1-0        DBI_1.2.2              
# [99] ggrastr_1.0.2           tidyselect_1.2.1        bit_4.0.5               compiler_4.3.2          DelayedArray_0.28.0     plotly_4.10.4           shadowtext_0.1.3       
# [106] scales_1.3.0            lmtest_0.9-40           stringr_1.5.1           digest_0.6.34           goftest_1.2-3           presto_1.0.0            spatstat.utils_3.1-0   
# [113] minqa_1.2.6             XVector_0.42.0          htmltools_0.5.7         pkgconfig_2.0.3         lme4_1.1-35.3           fastmap_1.1.1           rlang_1.1.3            
# [120] htmlwidgets_1.6.4       jquerylib_0.1.4         farver_2.1.1            zoo_1.8-12              jsonlite_1.8.8          BiocParallel_1.36.0     GOSemSim_2.28.1        
# [127] R.oo_1.26.0             RCurl_1.98-1.14         magrittr_2.0.3          GenomeInfoDbData_1.2.11 ggplotify_0.1.2         dotCall64_1.1-1         patchwork_1.2.0        
# [134] munsell_0.5.0           Rcpp_1.0.12             ape_5.8                 viridis_0.6.5           reticulate_1.35.0       stringi_1.8.3           ggraph_2.2.1           
# [141] zlibbioc_1.48.2         MASS_7.3-60             MAST_1.28.0             plyr_1.8.9              parallel_4.3.2          listenv_0.9.1           ggrepel_0.9.5          
# [148] deldir_2.0-2            Biostrings_2.70.3       graphlayouts_1.1.1      splines_4.3.2           tensor_1.5              hms_1.1.3               igraph_2.0.2           
# [155] spatstat.geom_3.2-8     RcppHNSW_0.6.0          reshape2_1.4.4          BiocManager_1.30.22     nloptr_2.0.3            tweenr_2.0.3            httpuv_1.6.14          
# [162] RANN_2.6.1              tidyr_1.3.1             purrr_1.0.2             polyclip_1.10-6         future_1.33.1           scattermore_1.2         ggforce_0.4.2          
# [169] rsvd_1.0.5              xtable_1.8-4            RSpectra_0.16-1         tidytree_0.4.6          later_1.3.2             viridisLite_0.4.2       ragg_1.3.0             
# [176] tibble_3.2.1            aplot_0.2.2             memoise_2.0.1           beeswarm_0.4.0          cluster_2.1.4           globals_0.16.2
