# ==========================
# third step of the analysis 
# iDEA GO results
# ==========================

library (openxlsx)
library (org.Mm.eg.db)
library (pheatmap)
library (ggplot2)

#base folder to work with
base <- "./scRNA-Seq/biorxiv"
no_save_base <- "./Miranda_et_al_noSave"

set.seed(2401)

# iDEA analysis was performed on ISTA cluster
# detailes in iDEA_analysis folder

# ======================================
# read, filter and document the results
# ======================================

idea_GO_results <- readRDS(file = paste(base, "/iDEA_analysis/idea_GO.idea_results.RDS", sep="" ) )
  
# this list contains the GO terms ranked by p-value
topGO <- list()
for ( genotype in names(idea_GO_results) ) {
  topGO[[genotype]] <- list()
  for (cluster in names( idea_GO_results[[genotype]] )) {
    topGO[[genotype]][[cluster]] <- list()
    for (dir in names( idea_GO_results[[genotype]][[cluster]] )) {
      topGO[[genotype]][[cluster]][[dir]] <- idea_GO_results[[genotype]][[cluster]][[dir]]@gsea[order(idea_GO_results[[genotype]][[cluster]][[dir]]@gsea$pvalue_louis, decreasing = F), c("annot_id", "pvalue_louis", "pvalue")]
      topGO[[genotype]][[cluster]][[dir]]$pvalue_adj <- p.adjust ( topGO[[genotype]][[cluster]][[dir]]$pvalue )
    }
  }
}
  
# supplement to document all GO analyses results
  wb <- createWorkbook()
  
  for (comp in names(topGO)) {
    for (cluster in names(topGO[[comp]])) {
      for (dir in c("all")) {
        tmp <- topGO[[ comp ]][[ cluster ]][[ dir ]][,c("annot_id", "pvalue", "pvalue_adj")]
        tmp <- tmp[ order (tmp$pvalue, decreasing = F), ]
        sheetName <- paste( comp, cluster, dir, sep="_" )
        
        addWorksheet(wb, sheetName)
        writeData(wb, sheetName, tmp)
        addFilter(wb, sheetName, row = 1, cols = 1:ncol(tmp))
        setColWidths(wb, sheetName, cols = 1:ncol(tmp), widths="auto")  
      }
    }
  }
  
  saveWorkbook(wb, file = paste(base, "/Supplement/iDEA_GO_analysis.xlsx", sep="/"), overwrite = TRUE) 
  
# =========================
# read and plot MAST DEGs
# =========================

DEG_object <- readRDS(file = paste(base, "/iDEA_analysis/DEG_iDEA.RDS", sep="" ) )

# plot DEGs values
all_deg <- data.frame()
all_deg_symbol <- list()
gene_list <- list()
for (genotype in c("KO_WT", "DKO_KO")) {
  all_deg_symbol[[ genotype ]] <- c()
  
  for (cluster in names(DEG_object[[genotype]])) {
    tmp <- DEG_object[[genotype]][[cluster]][,c(5,6)]
    
    # change adjusted p to a score
    tmp$p_val_adj <- log10( tmp$p_val_adj )
    tmp$p_val_adj[ which( DEG_object[[genotype]][[cluster]]$avg_log2FC > 0 ) ] <- tmp$p_val_adj[ which( DEG_object[[genotype]][[cluster]]$avg_log2FC > 0 ) ] * -1
    
    # filter for significance
    tmp <- tmp[which(tmp$p_val_adj > 1 | tmp$p_val_adj < -1 ), ]
    
    #all_deg_symbol[[ genotype ]] <- c( all_deg_symbol[[ genotype ]], tmp$gene)
                                       
    all_deg_symbol[[ genotype ]] <- c( all_deg_symbol[[ genotype ]], 
                                        paste( tmp$gene, sapply( tmp[,1], function (x){
                                                                   ifelse (x < 0, "neg", "pos")
                                                                 }), sep="_" ))
    
    # give a meaningful name
    colnames(tmp) <- c(paste(genotype, cluster, sep="."), "gene")
    
    # combine to one df
    if (nrow(all_deg) == 0) {
      all_deg <- tmp
    } else {
      all_deg <- merge(tmp, all_deg, by="gene", all=T)
    }
  }
  gene_list[[genotype]] <- all_deg$gene
}

# =========================
# plot specificity of DEGs
# bar plot
# =========================

spec_plot_df2 <- as.data.frame( table(table(all_deg_symbol[["DKO_KO"]])) )
spec_plot_df2$perc <- spec_plot_df2$Freq *100 / sum(spec_plot_df2$Freq)
ggplot(spec_plot_df2, aes(x=Var1, y=perc)) + geom_bar(stat="identity", position="stack") + 
  theme_classic() + xlab("# cell type with DE") + ylab("% of all DEGs") + ggtitle("Pten-MADM-Egfr_cKO / Pten-MADM")
ggsave  ( filename = paste(base, "/plots/main/Fig5F_DEG_specificity_plot.pdf", sep=""), width=4, height=4)

wb <- createWorkbook()

tmp <- spec_plot_df2

sheetName <- "DEGspec Pten-Egfr vs Pten"
addWorksheet(wb, sheetName)
writeData(wb, sheetName, tmp)
addFilter(wb, sheetName, row = 1, cols = 1:ncol(tmp))
setColWidths(wb, sheetName, cols = 1:ncol(tmp), widths="auto")  

saveWorkbook(wb, file = paste(base, "Supplement/Fig5F_DEG_specificity.xlsx", sep="/"), overwrite = T) 

# ==========================================
# plot as heatmap
# ==========================================

# convert to a df that can be plotted by pheatmap
rownames ( all_deg ) <- all_deg$gene
all_deg <- all_deg[ , c(2:ncol(all_deg)) ]

# cut score for better plotting
all_deg[ all_deg > 5 ] <- 5
all_deg[ all_deg < -5 ] <- -5

# Pten mosaic / Pten mosaic Egfr cKO
sample_order <- c("DKO_KO.cProg", "DKO_KO.aIP", "DKO_KO.iA", "DKO_KO.oligo", "DKO_KO.OBNB")

  tmp_deg_mat <- all_deg[ gene_list[[ 2 ]] ]
  
  # for clustering no NAs are allowed
  tmp_deg_mat[ is.na(tmp_deg_mat) ] <- 0
  clust_out = hclust( dist(tmp_deg_mat) )
  # plotting needs NAs to give a white color
  tmp_deg_mat[ tmp_deg_mat == 0 ] <- NA
  
  out_file <- "/plots/main/Fig5E_DEG_heatmap.pdf"
  pheatmap( tmp_deg_mat[ clust_out$labels[clust_out$order], ], scale = "none", 
            show_rownames = F, cluster_cols = F, cluster_rows = F, na_col = "white",
            filename = paste(base, out_file, sep="") )

tryCatch({
  dev.off()
}, error = function(e) {
  #print("An error occurred:", e$message)
  print( e )
})

# ====================================
# plot overview of iDEA GO enrichment
# ====================================
library ("simplifyEnrichment")

# add a column with GOID - conversion table comes from preparation of GO terms
term2id <- readRDS( file = paste(base, "/iDEA_analysis/term2id.RDS", sep = "") )
term2id_mod <- names(term2id)
names(term2id_mod) <- gsub(pattern = " ", replacement = "_", x =  term2id)

cluster_df_list <- list()
sig_mat_list <- list()

for (comp in c( "KO_WT", "DKO_KO" )) {
  
  message(comp)
  
  # comparison specific parameters
  # this part links the clusters identified by simplifyEnrichment to broader groups for analysis
  
  if (comp == "KO_WT") {
    
    # GO annotations
    all_GO_clusters <- c( 4, 8, 1, 6, 3, 9, 2, 5 )
    all_GO_cluster_names <- c( "cell cycle", "gliogenesis", "cell death", "telomere", "proliferation", "growth", "adhesion", "signaling")
    
    # deregulated pathways
    sig_pw_cluster <- 5
    cluster_order <- c( 4, 3, 6, 2, 5 )
    cluster_names <- c("Tor", "Mapk", "Egfr", "Jak/Stat", "Notch" )
    
  } else if (comp == "DKO_WT") {
    
    # GO annotations
    all_GO_clusters <- c( 2, 7, 10, 6, 12, 11)
    all_GO_cluster_names <- c( "signaling", "cell cycle", "gliogenesis", "cell death / telomere", "glycolytic process", "NSC maintenance")
    
    # deregulated pathways
    sig_pw_cluster <- 2
    cluster_order <- c( 10, 3, 8, 6, 13 )
    cluster_names <- c("Tor", "Mapk", "Egfr", "Jak/Stat", "Notch" )
    
  } else if (comp == "DKO_KO") {
    
    all_GO_clusters <- c( 3, 5, 2, 6, 9, 11 )
    all_GO_cluster_names <- c( "signaling", "cell cycle", "gliogenesis", "cell death", "telomere", "stem cell maintenance")
    
    sig_pw_cluster <- 3
    cluster_order <- c(10, 3, 1, 13, 18)
    cluster_names <- c("Tor", "Mapk", "Egfr", "Jak/Stat", "Notch" )
    
  } else {
    
    stop ("unknown comparison!")
    
  }
  
  
  all_GO_clusters_label  <- all_GO_cluster_names
  names( all_GO_clusters_label ) <- all_GO_clusters
  sig_pw_cluster_label <- cluster_names
  names( sig_pw_cluster_label ) <- cluster_order
  
  # create list of GO terms for analysis
  list_of_gsea <- lapply ( names( topGO[[ comp ]] ), function (x) {
    tmp <- topGO[[ comp ]][[ x ]][[ "all" ]]
    tmp$GOID <- term2id_mod[ tmp$annot_id ]
    return(tmp)
  })
  
  names( list_of_gsea) <- names( topGO[[ comp ]] )
  
  # simplifyEnrichment analysis
  # saves the output in a df
  # save the heatmap with descriptions
  pdf( file = paste(base, "/iDEA_analysis/post_iDEA/", comp, "_all_terms_heatmap.pdf", sep="") )
  simpleGO <- simplifyGOFromMultipleLists( lt = list_of_gsea, go_id_column = 5, padj_column = 3, 
                                           padj_cutoff = 0.05, db = 'org.Mm.eg.db', ont = "BP", order_by_size = T, verbose = F )
  dev.off()
  
  # add additional information to the simplifyEnrichment analysis
  simpleGO$Description <- term2id[ simpleGO$id ]
  cluster_size <- table( simpleGO$cluster )
  simpleGO$cluster_size <- cluster_size[ as.character(simpleGO$cluster) ]
  # group names are added here
  simpleGO$cluster_label <- all_GO_clusters_label[ as.character(simpleGO$cluster) ]
  # save the simplifyEnrichment for later reporting
  cluster_df_list[[ paste(comp, "all", sep="_") ]] <- simpleGO[ order(simpleGO$cluster_size, simpleGO$cluster, decreasing = T), ]
  
  # filter the GO term lists from individual cell types based on simplifyEnrichment output
  list_of_gsea <- lapply ( names( topGO[[ comp ]] ), function (x) {
    tmp <- topGO[[ comp ]][[x]][["all"]]
    tmp$GOID <- term2id_mod[ tmp$annot_id ]
    tmp <- tmp[ which ( tmp$GOID %in% simpleGO[which( simpleGO$cluster %in% all_GO_clusters ), "id"] ), ]
    return(tmp)
  })
  
  names( list_of_gsea ) <- names( topGO[[ comp ]] )
  
  # prepare a matrix where rows are GO groups (simplifyEnrichment clusters) and columns are cell types
  # the most significant p-value is extracted and used for heatmap plotting
  mat <- sapply( list_of_gsea, function (x) {
    sapply ( all_GO_clusters, function (y) {
      idx <- which( x$GOID %in% simpleGO[which(simpleGO$cluster == y), "id"])
      if (length(idx) == 0) {
        return(1)
      } else {
        tmp <- x[ which( x$GOID %in% simpleGO[which(simpleGO$cluster == y), "id"]), ]
        return(min(tmp$pvalue_adj))
      }
    })
  })
  
  rownames( mat ) <- all_GO_cluster_names
  
  sig_mat_list[[ comp ]] <- mat
  
  mat[ which(mat == 0) ] <- min(mat[which(mat > 0)])
  mat <- log10(mat) * -1
  # cut the values for easier plotting
  mat[ which(mat > 6) ] <- 6
  pheatmap( mat, scale="none", cluster_cols = F, main = comp, filename = paste(base, "/iDEA_analysis/post_iDEA/", comp, "_allGO_heatmap.pdf", sep="") )
  
  ######
  ###
  # focus on signaling pathways
  ###
  ######
  
  # extract only GO IDs from signaling cluster
  list_of_gsea <- lapply ( names( topGO[[ comp ]] ), function (x) {
    tmp <- topGO[[ comp ]][[x]][["all"]]
    tmp$GOID <- term2id_mod[ tmp$annot_id ]
    tmp <- tmp[ which ( tmp$GOID %in% simpleGO[which( simpleGO$cluster == sig_pw_cluster ), "id"] ), ]
    return(tmp)
  })
  names( list_of_gsea ) <- names( topGO[[ comp ]] )
  
  # simplifyEnrichment analysis
  
  mat = GO_similarity( unique ( unlist ( sapply( list_of_gsea, function (x) x$GOID ) ) ), ont = "BP" )
  
  # saves the output in a df
  # save the heatmap with descriptions
  pdf( file = paste(base, "/iDEA_analysis/post_iDEA/", comp, "_pathway_overview_heatmap.pdf", sep="") )
  if ( comp %in% c("KO_WT", "DKO_WT" ) ) {
    df = simplifyGO(mat, control = list("cutoff" = 0.65), order_by_size = T, verbose = F) 
  } else {
    df = simplifyGO(mat, control = list("cutoff" = 0.6), order_by_size = T, verbose = F)
  }
  dev.off()
  
  # add additional information to the simplifyEnrichment analysis
  df$Description <- term2id[ df$id ]
  cluster_size <- table( df$cluster )
  df$cluster_size <- cluster_size[ as.character(df$cluster) ]
  
  pw_GO_clusters_label  <- cluster_names
  names( pw_GO_clusters_label ) <- as.character(cluster_order)
  
  df$cluster_label <- pw_GO_clusters_label[ as.character(df$cluster) ]
  
  cluster_df_list[[ paste(comp, "sig_pw", sep="_") ]] <- df[ order(df$cluster_size, df$cluster, decreasing = T), ]
  
  mat <- sapply( list_of_gsea, function (x) {
    sapply ( cluster_order, function (y) {
      idx <- which( x$GOID %in% df[which(df$cluster == y), "id"])
      if (length(idx) == 0) {
        return(1)
      } else {
        tmp <- x[ which( x$GOID %in% df[which(df$cluster == y), "id"]), ]
        return(min(tmp$pvalue_adj))
      }
    })
  })
  
  rownames( mat ) <- cluster_names
  
  # matrix is saved here
  sig_mat_list[[ comp ]] <- mat
  # deal with p-values == 0
  mat[ which(mat == 0) ] <- min(mat[which(mat > 0)])
  mat <- log10(mat) * -1
  mat[ which(mat > 4) ] <- 4
  pheatmap( mat[ c("Mapk", "Tor", "Egfr", "Jak/Stat", "Notch"), c("cProg", "aIP", "iA", "oligo", "OBNB") ], 
            scale="none", cluster_cols = F, main = comp, cluster_rows = F,
            filename = paste(base, "/iDEA_analysis/post_iDEA/", comp, "_pathway_heatmap.pdf", sep="") )
  
  # in some cases pheatmap makes troubles with the graphical output - fix here
  tryCatch({
    dev.off()
  }, error = function(e) {
    #print("An error occurred:", e$message)
    print( e )
  })
  
}


# write out cluster assignment - quality control
wb <- createWorkbook()

for (comp in names(cluster_df_list)) {
  tmp <- cluster_df_list[[ comp ]]
    
  sheetName <- comp
    
  addWorksheet(wb, sheetName)
  writeData(wb, sheetName, tmp)
  addFilter(wb, sheetName, row = 1, cols = 1:ncol(tmp))
  setColWidths(wb, sheetName, cols = 1:ncol(tmp), widths="auto")  
  
}

saveWorkbook(wb, file = paste(base, "/Supplement/GO_cluster_assignment.xlsx", sep="/"), overwrite = TRUE) 

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
#   [1] grid      stats4    stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] simplifyEnrichment_1.12.0 ggplot2_3.5.0             pheatmap_1.0.12           org.Mm.eg.db_3.18.0       AnnotationDbi_1.64.1      IRanges_2.36.0           
# [7] S4Vectors_0.40.2          Biobase_2.62.0            BiocGenerics_0.48.1       openxlsx_4.2.5.2         
# 
# loaded via a namespace (and not attached):
#   [1] tidyselect_1.2.1        dplyr_1.1.4             farver_2.1.1            blob_1.2.4              Biostrings_2.70.3       bitops_1.0-7            fastmap_1.1.1          
# [8] RCurl_1.98-1.14         digest_0.6.34           lifecycle_1.0.4         cluster_2.1.4           Cairo_1.6-2             NLP_0.3-0               KEGGREST_1.42.0        
# [15] RSQLite_2.3.6           magrittr_2.0.3          compiler_4.3.2          rlang_1.1.3             tools_4.3.2             utf8_1.2.4              labeling_0.4.3         
# [22] bit_4.0.5               xml2_1.3.6              RColorBrewer_1.1-3      withr_3.0.0             fansi_1.0.6             GOSemSim_2.28.1         tm_0.7-14              
# [29] colorspace_2.1-0        GO.db_3.18.0            scales_1.3.0            iterators_1.0.14        cli_3.6.2               crayon_1.5.2            ragg_1.3.0             
# [36] generics_0.1.3          rstudioapi_0.16.0       httr_1.4.7              rjson_0.2.23            commonmark_1.9.1        DBI_1.2.2               cachem_1.0.8           
# [43] stringr_1.5.1           zlibbioc_1.48.2         parallel_4.3.2          XVector_0.42.0          proxyC_0.4.1            matrixStats_1.2.0       vctrs_0.6.5            
# [50] yulab.utils_0.1.4       Matrix_1.6-5            slam_0.1-50             GetoptLong_1.0.5        bit64_4.0.5             clue_0.3-65             systemfonts_1.0.6      
# [57] packrat_0.9.2           foreach_1.5.2           glue_1.7.0              codetools_0.2-19        stringi_1.8.3           gtable_0.3.4            shape_1.4.6.1          
# [64] GenomeInfoDb_1.38.8     ComplexHeatmap_2.18.0   munsell_0.5.0           tibble_3.2.1            pillar_1.9.0            GenomeInfoDbData_1.2.11 circlize_0.4.16        
# [71] R6_2.5.1                textshaping_0.3.7       doParallel_1.0.17       lattice_0.21-9          markdown_1.13           png_0.1-8               gridtext_0.1.5         
# [78] memoise_2.0.1           Rcpp_1.0.12             zip_2.3.1               xfun_0.42               org.Hs.eg.db_3.18.0     fs_1.6.3                pkgconfig_2.0.3        
# [85] GlobalOptions_0.1.2   
