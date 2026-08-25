# ===================================================================
# determine the abundance of RGPs in bulk RNA/Seq data
# Note: this script is the only one that was executed with R v4.5.3
# ===================================================================

library (Seurat)
library (SingleCellExperiment)
library (openxlsx)
library (dplyr)
library (tidyr)
library (ggplot2)
library (ggbeeswarm)

set.seed(2401)

# =================
# helper functions
# =================

# ggplot2 style color choice
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

# ===============================
# Define a unified color scheme
# ===============================

color_vec <- c( gg_color_hue(5), "grey50" )
names(color_vec) <- c("aIP","RGP", "Astro", "OBNB", "Oligo", "Neuron")

# ===============================
# read scRNA-Seq data
# ===============================
no_save_base <- "./Miranda_et_al_noSave"

seurat.obj <- readRDS( file = paste(no_save_base, "/RDS_files/rev.seurat.glia.obj.rds", sep="") )

# ===============================
# read bulk-Seq data
# ===============================
initial_base <- "./bulk_RNASeq/biorxiv" # I get some information from the initial analysis
base <- "./bulk_RNASeq/final"

#prepare environment hash table as lookup for ENSMUST - ENSMUSG - Symbol and ENTREZID
#this comes from STAR index folder and was prepared solely by using info from Gencode
#initialise variables
ensmusg_symbol <- new.env()
ensmusg_chr <- new.env()
ensmusg_entrezid <- new.env()

#read conversion table
ensmusg_symbol_chr <- read.table(paste(initial_base, "annotation/exchange_table_M27.tsv", sep="/"), fill = T, header = F, stringsAsFactors = F)
colnames(ensmusg_symbol_chr) <- c("ENSMUSG", "SYMBOL", "chr")

#fill conversion table - handled just like a list but apparently much faster
for (x in 1:nrow(ensmusg_symbol_chr)) {
  #deal here with multi mappers - one SYMBOL mapped to many ENSMUSG
  #a comma in the SYMBOL name indicates these multi-mappers
  if (is.null(ensmusg_symbol[[ensmusg_symbol_chr$ENSMUSG[x]]])) {
    ensmusg_symbol[[ensmusg_symbol_chr$ENSMUSG[x]]] <- ensmusg_symbol_chr$SYMBOL[x]
  } else {
    
    if (!ensmusg_symbol[[ensmusg_symbol_chr$ENSMUSG[x]]] == ensmusg_symbol_chr$SYMBOL[x]) {
      ensmusg_symbol[[ensmusg_symbol_chr$ENSMUSG[x]]] <- paste(ensmusg_symbol[[ensmusg_symbol_chr$ENSMUSG[x]]], ensmusg_symbol_chr$SYMBOL[x], sep=",")
    }
  }
}

#add an entry for the deletion I add manually
ensmusg_symbol[["Pten_ex5"]] <- "Pten_ex5"

#read sample list
sample_list_xls <- read.xlsx(xlsxFile =paste(initial_base, 'Supplement/final_sample_list.xlsx', sep="/"), sheet = 1, skipEmptyRows = F)
sample_list_xls <- sample_list_xls[ which(sample_list_xls$used == "Y"),]

#read count table
#initialise dataframe holding all read counts df
all <- data.frame(stringsAsFactors = F)

#read in all counts and merge in one dataframe - by gene_id
for (line in 1:nrow(sample_list_xls)) {
  this_count <- read.table(file = sample_list_xls$fn[line], header = F, sep="\t", stringsAsFactors = F)
  #remove first 4 lines - only basic stats of features
  this_count<-this_count[5:nrow(this_count),c(1,2)]
  
  colnames(this_count) <- c("gene_id", line)
  
  if (nrow(all) == 0) {
    all <- this_count
  } else {
    all<-merge(all, this_count, by="gene_id")
  }
}

# ============================================
# prepare count table for downstream analyses
# ============================================

colnames(all)<-c("gene_id", sample_list_xls$sample_name)
gn <- all$gene_id
symbol <- sapply(gn, function (x) ensmusg_symbol[[x]])
symbol <- as.character( symbol )
gn <- data.frame( ensembl = gn, symbol = symbol)
gn <- gn[!(duplicated(gn$symbol) | duplicated(gn$symbol, fromLast = TRUE)), ]
rownames(gn) <- gn$ensembl

rownames(all) <- all$gene_id

all<-all[gn$ensembl,]
rownames(all) <- gn$symbol

seurat.obj <- JoinLayers( seurat.obj )

# I use the cluster annotation - gave best results for me
Idents( seurat.obj ) <- "seurat_clusters"
RGP_cluster <- c("7", "0", "9", "11")
astro_cluster <- c("4", "6", "5", "2", "8")
seurat.obj <- subset( seurat.obj, idents = c(RGP_cluster, astro_cluster), invert = F)

# ===========
# BayesPrism
# ===========

library(BayesPrism)

sc_counts <- t(as.matrix(GetAssayData(seurat.obj, assay = "RNA", layer = "counts")))
# BayesPrism expects cells × genes

cell_type_labels    <- seurat.obj$cell_type
cell_state_labels   <- seurat.obj$seurat_clusters  # optional finer resolution

sc_dat_filtered <- cleanup.genes(
  input      = sc_counts,
  input.type = "count.matrix",
  species    = "hs",            # "hs" for human, "mm" for mouse
  gene.group = c("Rb", "Mrp", "other_Rb", "chrM", "MALAT1", "chrX", "chrY"),
  exp.cells  = 5                # exclude genes expressed in fewer than 5 cells
)

bulk_matrix <- t(as.matrix(all[,c(2:ncol(all))]))
# BayesPrism expects samples × genes

myPrism <- new.prism(
  reference         = sc_dat_filtered,
  mixture           = bulk_matrix,
  input.type        = "count.matrix",
  cell.type.labels  = cell_state_labels,
  cell.state.labels = cell_state_labels,
  key               = NULL,          
  outlier.cut       = 0.01,
  outlier.fraction  = 0.1
)

# number of cells in each cell state 
# cell.state.labels
# 11   9   8   7   6   5   4   2   0 
# 74 114 127 130 138 199 223 264 288 
# No tumor reference is speficied. Reference cell types are treated equally. 
# Number of outlier genes filtered from mixture = 8 
# Aligning reference and mixture... 
# Normalizing reference... 

# do the analysis - longest part of the script
bp_res <- run.prism(
  prism   = myPrism,
  n.cores = 6        # parallelise — this step is slow
)

# Cell type proportions (samples × cell types)
props_bp <- get.fraction(
  bp        = bp_res,
  which.theta = "final",    # "final" = after convergence, "first" = initialisation
  state.or.type = "type"    # "type" for coarse, "state" for fine-grained
)

props_bp <- as.data.frame(round(props_bp, 3))

rownames( sample_list_xls ) <- sample_list_xls$sample_name
props_bp$age <- sample_list_xls[ rownames(props_bp), "age" ]

# summarize the results and plot
df_plot <- props_bp %>%
  mutate(
    group_rgp_sum = rowSums(select(., all_of(RGP_cluster)), na.rm = TRUE),
    group_astro_sum = rowSums(select(., all_of(astro_cluster)), na.rm = TRUE)
  ) %>%
  select(age, group_rgp_sum, group_astro_sum)

df_long_BP <- df_plot %>%
  pivot_longer(cols = starts_with("group"), 
               names_to = "group", 
               values_to = "sum_value")

summary_long_BP <- df_long_BP %>%
  group_by(age, group) %>%
  summarize(
    mean_val = mean(sum_value),
    sd_val   = sd(sum_value),
    .groups = "drop"
  )

summary_long_BP$simple_group <- "RGP"
summary_long_BP$simple_group[ which( grepl(pattern = "astro", x = summary_long_BP$group) ) ] <- "Astro"

ggplot() + 
  geom_bar(data = summary_long_BP,
             aes(x = group, y = mean_val, fill = simple_group),
             stat="identity", position="dodge") +
  geom_beeswarm(data = df_long_BP,
              aes(x = group, y = sum_value),
              cex=3, size=2) +
  facet_grid( ~ age) + theme_classic() +
  scale_fill_manual( values = color_vec )  +
  theme(axis.text.x = element_text(angle = 90)) +
  ggtitle("BayesPrism deconvolution")
ggsave( file=paste (base, "/plots/Rev_Fig1A_bulk_deconv_BP.pdf", sep=""))

# ====================
# MuSiC deconvolution
# ====================
library (MuSiC)
seurat.obj$seurat_clusters <- factor(seurat.obj$seurat_clusters)
sce_ref <- as.SingleCellExperiment(seurat.obj, assay = "RNA")

# Create a separate ExpressionSet for bulk
bulk_matrix <- as.matrix(all[,c(2:ncol(all))])

deconv <- music_prop(
  bulk.mtx  = bulk_matrix,   # genes × samples raw count matrix
  sc.sce    = sce_ref,                       # SingleCellExperiment object
  clusters  = "seurat_clusters",
  samples   = "batch",
  verbose   = TRUE
)

props <- as.data.frame( deconv$Est.prop.weighted )
props$age <- sample_list_xls[ rownames(props), "age" ]

# summarize the results and plot
df_plot <- props %>%
  mutate(
    group_rgp_sum = rowSums(select(., all_of(RGP_cluster)), na.rm = TRUE),
    group_astro_sum = rowSums(select(., all_of(astro_cluster)), na.rm = TRUE)
  ) %>%
  select(age, group_rgp_sum, group_astro_sum)

df_long <- df_plot %>%
  pivot_longer(cols = starts_with("group"), 
               names_to = "group", 
               values_to = "sum_value")

summary_long <- df_long %>%
  group_by(age, group) %>%
  summarize(
    mean_val = mean(sum_value),
    sd_val   = sd(sum_value),
    .groups = "drop"
  )

summary_long$simple_group <- "RGP"
summary_long$simple_group[ which( grepl(pattern = "astro", x = summary_long$group) ) ] <- "Astro"

ggplot() + 
  geom_bar(data = summary_long,
           aes(x = group, y = mean_val, fill = simple_group),
           stat="identity", position="dodge") +
  geom_beeswarm(data = df_long,
                aes(x = group, y = sum_value),
                cex=3, size=2) +
  facet_grid( ~ age) + theme_classic() +
  theme(axis.text.x = element_text(angle = 90)) +
  scale_fill_manual( values = color_vec ) +
  ggtitle("MuSiC deconvolution")

ggsave( file=paste(base, "/plots/Rev_Fig1B_bulk_deconv_MuSiC.pdf", sep=""))

# ================================================
# write out the raw data and anova for this plot
# ================================================

raw_plot_data <- list()
raw_plot_data[["BP_average"]] <- summary_long_BP
raw_plot_data[["BP_datapoints"]] <- df_long_BP
raw_plot_data[["MuSiC_average"]] <- summary_long
raw_plot_data[["MuSiC_datapoints"]] <- df_long

wb <- createWorkbook()

for (name in names(raw_plot_data)) {
  sheetName <- name
  tmp <- raw_plot_data[[ name ]]
  addWorksheet(wb, sheetName)
  writeData(wb, sheetName, tmp)
  addFilter(wb, sheetName, row = 1, cols = 1:ncol(tmp))
  setColWidths(wb, sheetName, cols = 1:ncol(tmp), widths="auto")  
}

saveWorkbook(wb, file = paste(base, "/plots/deconv_plot_data.xlsx", sep="/"), overwrite = T) 

sessionInfo()
# R version 4.5.3 (2026-03-11)
# Platform: x86_64-pc-linux-gnu
# Running under: Ubuntu 22.04.5 LTS
# 
# Matrix products: default
# BLAS:   /opt/R/4.5.3/lib/R/lib/libRblas.so 
# LAPACK: /usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3.10.0  LAPACK version 3.10.0
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
# [1] MuSiC_1.0.0                 TOAST_1.24.0                quadprog_1.5-8              limma_3.66.0               
# [5] EpiDISH_2.26.0              nnls_1.6                    BayesPrism_2.2.3            NMF_0.28                   
# [9] cluster_2.1.8.2             rngtools_1.5.2              registry_0.5-1              snowfall_1.84-6.3          
# [13] snow_0.4-4                  ggbeeswarm_0.7.3            ggplot2_4.0.2               tidyr_1.3.2                
# [17] dplyr_1.2.0                 openxlsx_4.2.8.1            SingleCellExperiment_1.32.0 SummarizedExperiment_1.40.0
# [21] Biobase_2.70.0              GenomicRanges_1.62.1        Seqinfo_1.0.0               IRanges_2.44.0             
# [25] S4Vectors_0.48.1            BiocGenerics_0.56.0         generics_0.1.4              MatrixGenerics_1.22.0      
# [29] matrixStats_1.5.0           Seurat_5.4.0                SeuratObject_5.3.0          sp_2.2-1                   
# 
# loaded via a namespace (and not attached):
#  [1] RcppAnnoy_0.0.23       splines_4.5.3          later_1.4.8            bitops_1.0-9           tibble_3.3.1          
#  [6] polyclip_1.10-7        fastDummies_1.7.5      lifecycle_1.0.5        edgeR_4.8.2            doParallel_1.0.17     
# [11] globals_0.19.1         lattice_0.22-9         MASS_7.3-65            magrittr_2.0.4         plotly_4.12.0         
# [16] metapod_1.18.0         httpuv_1.6.17          otel_0.2.0             sctransform_0.4.3      zip_2.3.3             
# [21] spam_2.11-3            spatstat.sparse_3.1-0  reticulate_1.45.0      cowplot_1.2.0          pbapply_1.7-4         
# [26] RColorBrewer_1.1-3     abind_1.4-8            Rtsne_0.17             purrr_1.2.1            ggrepel_0.9.8         
# [31] irlba_2.3.7            listenv_0.10.1         spatstat.utils_3.2-2   goftest_1.2-3          MatrixModels_0.5-4    
# [36] RSpectra_0.16-2        dqrng_0.4.1            spatstat.random_3.4-5  fitdistrplus_1.2-6     parallelly_1.46.1     
# [41] codetools_0.2-20       DelayedArray_0.36.1    scuttle_1.20.0         tidyselect_1.2.1       locfdr_1.1-8          
# [46] farver_2.1.2           ScaledMatrix_1.18.0    viridis_0.6.5          spatstat.explore_3.8-0 jsonlite_2.0.0        
# [51] BiocNeighbors_2.4.0    e1071_1.7-17           progressr_0.18.0       ggridges_0.5.7         survival_3.8-6        
# [56] scater_1.38.1          iterators_1.0.14       systemfonts_1.3.2      foreach_1.5.2          tools_4.5.3           
# [61] ragg_1.5.2             ica_1.0-3              Rcpp_1.1.1             glue_1.8.0             gridExtra_2.3         
# [66] SparseArray_1.10.10    withr_3.0.2            BiocManager_1.30.27    fastmap_1.2.0          GGally_2.4.0          
# [71] bluster_1.20.0         SparseM_1.84-2         caTools_1.18.3         digest_0.6.39          rsvd_1.0.5            
# [76] R6_2.6.1               mime_0.13              textshaping_1.0.5      colorspace_2.1-2       scattermore_1.2       
# [81] gtools_3.9.5           tensor_1.5.1           spatstat.data_3.1-9    data.table_1.18.2.1    corpcor_1.6.10        
# [86] class_7.3-23           httr_1.4.8             htmlwidgets_1.6.4      S4Arrays_1.10.1        ggstats_0.13.0        
# [91] uwot_0.2.4             pkgconfig_2.0.3        gtable_0.3.6           lmtest_0.9-40          S7_0.2.1              
# [96] XVector_0.50.0         htmltools_0.5.9        dotCall64_1.2          MCMCpack_1.7-1         scales_1.4.0          
# [101] png_0.1-9              spatstat.univar_3.1-7  scran_1.38.1           rstudioapi_0.18.0      reshape2_1.4.5        
# [106] nlme_3.1-168           coda_0.19-4.1          proxy_0.4-29           zoo_1.8-15             stringr_1.6.0         
# [111] KernSmooth_2.23-26     parallel_4.5.3         miniUI_0.1.2           vipor_0.4.7            pillar_1.11.1         
# [116] grid_4.5.3             vctrs_0.7.2            gplots_3.3.0           RANN_2.6.2             promises_1.5.0        
# [121] BiocSingular_1.26.1    beachmat_2.26.0        xtable_1.8-8           beeswarm_0.4.0         locfit_1.5-9.12       
# [126] packrat_0.9.3          cli_3.6.5              compiler_4.5.3         rlang_1.1.7            future.apply_1.20.2   
# [131] labeling_0.4.3         plyr_1.8.9             stringi_1.8.7          gridBase_0.4-7         deldir_2.0-4          
# [136] viridisLite_0.4.3      BiocParallel_1.44.0    lazyeval_0.2.2         spatstat.geom_3.7-3    quantreg_6.1          
# [141] Matrix_1.7-5           RcppHNSW_0.6.0         patchwork_1.3.2        future_1.70.0          statmod_1.5.2         
# [146] mcmc_0.9-8             shiny_1.13.0           ROCR_1.0-12            igraph_2.2.2     
