# ======================================
# bulk analysis of purified astrocytes
# ======================================

library (openxlsx)
library (DESeq2)
library (ggplot2)
library (clusterProfiler)
library (org.Mm.eg.db)
library (reshape2)
library (ggbeeswarm)
library (cowplot)
library (ggVennDiagram)

set.seed(2401)

base <- "./bulk_RNASeq/biorxiv"
figure_base <- paste(base, 'plots/main', sep='/')
sup_base <- paste(base, 'plots/Supplement', sep='/')
qc_base <- paste(base, 'plots/QC', sep='/')
gene_count_base <- paste(base, 'STAR_counts', sep='/') # count tables go here

# =====================================================================================
# prepare environment hash table as lookup for ENSMUST - ENSMUSG - Symbol and ENTREZID
# was prepared solely by using info from Gencode
# =====================================================================================

# initialise variables
  ensmusg_symbol <- new.env()
  ensmusg_chr <- new.env()
  ensmusg_entrezid <- new.env()
  
  #read conversion table
  ensmusg_symbol_chr <- read.table(paste(base, "annotation/exchange_table_M27.tsv", sep="/"), fill = T, header = F, stringsAsFactors = F)
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

# ==================================================================
# prepare sample list, combine count tables and add alignment info
# ==================================================================
 
#read sample list
sample_list_xls <- read.xlsx(xlsxFile =paste(base, 'other_data/R13415_Astrocytes.xlsx', sep="/"), sheet = 1, skipEmptyRows = F)

#match samples with read count files - done via sampleID

#recursive search for all files in the count folder - returns relative path
files <- list.files(path = gene_count_base, pattern="*ReadsPerGene.out.tab", recursive = T)
#make path absolute
full_fn_path_list <- lapply(files, function (x) {paste(gene_count_base, x, sep="/")})
#prepare dataframe from file list - matching ID with filename
files.DF <- data.frame(fn = unlist(full_fn_path_list), stringsAsFactors = F)
files.DF$sample_id <- apply(files.DF, 1, function (x) {
  #extract samle ID from filename
  parts <- strsplit(x[1], "\\.")
  return(parts[[1]][length(parts[[1]])-3])
})

sample_list_xls <- merge(sample_list_xls, files.DF, by="sample_id")

#read count tables
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

#give meaningful column names
colnames(all)<-c("gene_id", sample_list_xls$sample_id)

# log files are available as STAR_Logs.tar.gz
files.STAR <- c(list.files(path=gene_count_base, pattern=glob2rx("*.Log.final.out"), recursive = T))

STAR_df2plot <- data.frame()
for (file_name in files.STAR) {
  STAR.Log <- read.table(file = paste(gene_count_base, file_name, sep="/"), header = F, fill = T, sep="\t", stringsAsFactors = F)
  
  #isolate alignment folder name
  folders <- strsplit(x = file_name, split = "\\/")[[1]] 
  
  #if the folder name contains a dot this indicates that the sample ID needs to be extracted
  #otherwise the sample_id is in the filename
  sample_name <- strsplit(x = folders[[length(folders)]], split = "\\.")[[1]][2]
  
  STAR_df2plot <- rbind(STAR_df2plot, data.frame(sample_name=sample_name, STAR.unique.aligned=STAR.Log$V2[8], 
                                                 STAR.tot=STAR.Log$V2[5], STAR.perc = STAR.Log$V2[9], stringsAsFactors = F))
}

#extract some important infos from STAR alignment
#especially alignment rate is an important quality measure
sample_df_alignment_stats <- merge(STAR_df2plot, sample_list_xls, by.x="sample_name", by.y="sample_id")
sample_df_alignment_stats$STAR.tot <- as.numeric(sample_df_alignment_stats$STAR.tot)
sample_df_alignment_stats$STAR.unique.aligned <- as.numeric(sample_df_alignment_stats$STAR.unique.aligned)
sample_df_alignment_stats$STAR.perc <- sapply(sample_df_alignment_stats$STAR.perc, function (x) as.numeric(gsub(pattern = '%', replacement = '', x = x)))

#remove samples with too low alignment rate
#set an arbitrary cutoff
idx2rem <- sample_df_alignment_stats$sample_name[which(sample_df_alignment_stats$STAR.perc < 45)]

sample_list_xls[which(sample_list_xls$sample_id %in% idx2rem), c(1:6)]

#    sample_id genotype cell_type cell_number age sort_dat
# 6     192786       WT     Astro          22  P0 16.02.22
# 12    192803       WT     Astro         329  P4 22.03.22
# 20    192819       WT     Astro          93  P4 23.03.22

#create a group variable for correlation analysis of biological replicates
sample_list_xls$group <- paste( sample_list_xls$cell_type, sample_list_xls$genotype, sample_list_xls$age, sep="_")

# ==============================
# basic QC:
# PCA and correlation analysis
# ==============================

#prepare count table for later analysis
all <- all[, c("gene_id", as.character(sample_list_xls$sample_id))]

for (rem in c(1:3)) {
  
  #remove low quality samples
  sample_list_xls <- sample_list_xls[which(!sample_list_xls$sample_id %in% idx2rem),]
  
  #process all samples in comparisons
  all_filter <- all[, c("gene_id", sapply(sample_list_xls$sample_id, as.character))]
  all_filter <- all_filter[which(rowMeans(all[,2:ncol(all)]) > 10), ]
  
  #prepare countdata - change class to integer (read counts have to be integer)
  countData <- sapply(all_filter[,c(2:(ncol(all_filter)))], as.integer)
  #add rownames - name of genes
  rownames(countData) <- all_filter[,c(1)]
  
  #prepare DESeq2 comparison - group by "sample" date in the design formula
  vsd <- varianceStabilizingTransformation(countData, blind=T)
  
  rv <- rowVars(vsd)
  ntop <- 500
  select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
  pca <- prcomp(t(vsd[select, ]))
  
  d <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2])
  d$sample <- rownames(d)
  
  d$cell_type <- sapply(d$sample, function (x) {
    idx <- which(sample_list_xls$sample_id == x)
    out <- sample_list_xls$cell_type[idx]
    return(out)
  })
  
  d$age <- sapply(d$sample, function (x) {
    idx <- which(sample_list_xls$sample_id == x)
    out <- sample_list_xls$age[idx]
    return(out)
  })
  
  d$genotype <- sapply(d$sample, function (x) {
    idx <- which(sample_list_xls$sample_id == x)
    out <- sample_list_xls$genotype[idx]
    return(out)
  })
  
    
  d$sample_id <- rownames(d)
  plot1 <- ggplot(d, aes(x=PC1, y=PC2, color=age, shape=genotype)) + 
            geom_point(size=3) + theme_classic()
  plot2 <- ggplot(d, aes(x=PC1, y=PC2, color=cell_type, shape=genotype)) + 
    geom_point(size=3) + theme_classic() + geom_text(aes(label=sample_id), color="black")
  
  comb_plot <- plot_grid( plot1, plot2)
  
  ggsave(filename = paste(qc_base, "/PCA.", rem, ".pdf", sep=""), plot = comb_plot, width = 8, height = 4)

  #test for correlations between biological replicates
  df2plot <- data.frame()
  for (group in unique(sample_list_xls$group)) {
    
    ids <- sample_list_xls[which(sample_list_xls$group == group), "sample_id"]
    if (length(ids) == 1) {
      ids <- c(ids, ids)
    }
    
    pw_comp <- combn(ids, 2)
    
    for (i in 1: ncol(pw_comp)) {
      comp1 <- which(colnames(vsd) == pw_comp[1,i])
      comp2 <- which(colnames(vsd) == pw_comp[2,i])
      
      sample1 <- sample_list_xls[which(sample_list_xls$sample_id == pw_comp[1,i]), "sample_id"]
      sample2 <- sample_list_xls[which(sample_list_xls$sample_id == pw_comp[2,i]), "sample_id"]
      
      this_cor <- cor(vsd[,comp1], vsd[,comp2])
      
      df2plot <- rbind(df2plot, data.frame(group = group, cor=this_cor, sample1=sample1, sample2=sample2, stringsAsFactors = F))
    }
  }
  
  cor_plot <- ggplot(df2plot, aes(x=group, y=cor)) + 
            geom_boxplot() + 
            geom_point() + 
            theme_classic() + 
            ylim(0,1) + theme(axis.text.x = element_text(angle = 90, hjust = 1))
  
  ggsave(filename = paste(qc_base, "/corr_pw_comp_samples.", rem, ".pdf", sep=""), plot = cor_plot)
  
  if (rem == 1) {
    # low correlation with others
    idx2rem <- c(idx2rem, 192823, 192807 )
  }
  
  if (rem == 2) {
    # position on PCA plot
    idx2rem <- c(idx2rem, 192788, 192801, 192821 )
  }
}

# save the final PCA as Supplemental Figure
ggsave( filename = paste(sup_base, "/Sup6A_PCA_Astro.pdf", sep="/"), plot = plot1, width = 5, height = 5 )

# get the number of replicates
table(sample_list_xls$group)
# Astro_KO_P0 Astro_KO_P4 Astro_WT_P0 Astro_WT_P4 
#           5           3           3           3

#prepare a df for output of sample list
sample_df_alignment_stats$used <- 'N'
sample_df_alignment_stats$used[which(sample_df_alignment_stats$sample_name %in% sample_list_xls$sample_id)] <- 'Y'

#write out sample list and details
wb <- createWorkbook()

tmp <- sample_df_alignment_stats
sheetName <- 'sample list'
addWorksheet(wb, sheetName)
writeData(wb, sheetName, tmp)
addFilter(wb, sheetName, row = 1, cols = 1:ncol(tmp))
setColWidths(wb, sheetName, cols = 1:ncol(tmp), widths="auto")  

saveWorkbook(wb, file = paste(base, "Supplement/final_sample_list.xlsx", sep="/"), overwrite = T) 

###############################################################
#
# read information on deleted region of this analysis - Pten KO
# read all files with information on read coverage
#
###############################################################

# these count files are available as Pten_del_counts.tar.gz
files <- list.files(path = gene_count_base, pattern="Pten_del.counts.out", recursive = T)
#make path absolute
full_fn_path_list <- sapply(files, function (x) {paste(gene_count_base, x, sep="/")})

delcounts <- data.frame()
for (file in files) {
  #this number heavily relies on folder structure - only works in this set-up!
  seq_id <- strsplit(x = file, split = "/|\\.")[[1]][2]
  if (seq_id %in% sample_list_xls$sample_id) {
    tmp <- read.table(file = paste(gene_count_base, file, sep="/"), header = F, stringsAsFactors = F)
    tmp <- tmp[,c("V6", "V8")]
    colnames(tmp) <- c("gene_id", seq_id)
    if (nrow(delcounts) == 0){
      delcounts <- tmp
    } else {
      delcounts <- merge(delcounts, tmp, by="gene_id")
    }
  }
}

# ==========================
# plot Pten deletion counts
# ==========================

df2plot <- reshape2::melt(delcounts)
df2plot$group <- sapply(df2plot$variable, function (x) sample_list_xls[which(sample_list_xls$sample_id == x), 'group'])
df2plot <- merge(df2plot, STAR_df2plot, by.x='variable', by.y='sample_name')

#calculate reads per M uniquely aligned reads - nicer than the scaling done by DESeq2
df2plot$norm_reads <- df2plot$value / (as.numeric(df2plot$STAR.unique.aligned) / 1000000)
df2plot$cell_type <- sapply(df2plot$group, function (x) strsplit(x = x, split = '_')[[1]][1])
df2plot$genotype <- sapply(df2plot$group, function (x) strsplit(x = x, split = '_')[[1]][2])
df2plot$age <- sapply(df2plot$group, function (x) strsplit(x = x, split = '_')[[1]][3])

Pten_del_count_plot <- ggplot(df2plot, aes(x=genotype, y=norm_reads)) + 
  geom_boxplot() + 
  geom_beeswarm() + 
  theme_classic() + 
  ylab('normalised read count') +
  facet_wrap(~age, nrow=1) +
  ggtitle('Pten expression')

ggsave( filename = paste(sup_base, "/Sup6B_Pten_ex5_counts.pdf", sep="/"), plot = Pten_del_count_plot, width = 4 )

#write to xlsx
wb <- createWorkbook()
sheetName <- "Pten norm counts"
addWorksheet(wb, sheetName)
writeData(wb, sheetName, df2plot)
addFilter(wb, sheetName, row = 1, cols = 1:ncol(df2plot))
setColWidths(wb, sheetName, cols = 1:ncol(df2plot), widths="auto")
saveWorkbook(wb, file = paste(base, "Supplement/all_samples_Pten_normalised_counts.xlsx", sep=""), overwrite = T) 

#do DE analysis between genotypes of each age
res <- list()

for (age in c("P0", "P4"))  {
  
  #extract sub-group for analysis
  sl_txt_sub <- sample_list_xls[which(sample_list_xls$age == age),]

  #extract the columns from expression data fitting to the list of samples extracted above
  #also ectract gene_id!
  colnames_to_extract <- c("gene_id", sl_txt_sub$sample_id)
  all_filter <- all[ , colnames_to_extract ]
  
  #add delcounts data
  all_filter <- rbind(all_filter, delcounts[,colnames(all_filter)])
  all_filter <- all_filter[which(rowMeans(all_filter[,2:ncol(all_filter)]) > 10), ]
  
  cat(age, nrow(all_filter), "\n")   
  
  #prepare metadata for later DESeq2 analysis    
  meta<-data.frame(genotype = sl_txt_sub$genotype, sample_id = sl_txt_sub$sample_id)
  rownames(meta) <- colnames(all_filter[,2:(ncol(all_filter))])
  
  if (!any(rownames(meta) == meta$sample_id)) {
    warning("problem with: ", group)
  }
  
  #prepare countdata - change class to integer (read counts have to be integer)
  countData <- sapply(all_filter[,c(2:(ncol(all_filter)))], as.integer)
  #add rownames - name of genes
  rownames(countData) <- all_filter[,c(1)]

  #prepare DESeq2 comparison
  ddsMat<-DESeqDataSetFromMatrix(countData = countData,
                                 colData=meta,
                                 design= ~ genotype)
  
  #do DESeq2 analysis
  dds<-DESeq(ddsMat, fitType="local", quiet=T, betaPrior = F, parallel = F)
  
  #prepare a list of contrasts to be performed
  groups_contrasts <- list()
  groups_contrasts[[1]] <- c("genotype", "KO", "WT")

  #prepare a description to give meaningful names to the list entries
  descr <- c("KO_WT")
  gl <- c()
  
  #prepare lists of DE genes using contrasts
  for (x in 1:length(groups_contrasts)) {
    
    this_res <- results(dds, contrast = groups_contrasts[[x]])
    this_res <- this_res[order(this_res$padj),]
    
    #add gene symbols to DE results
    gn <- rownames(this_res)
    this_res$symbol <- sapply(gn, function (x) ensmusg_symbol[[x]])
    this_res$symbol <- as.character( this_res$symbol )
    
    #here the DE analysis is saved
    res[[paste(age, descr[x], sep="_")]] <- as.data.frame(this_res)
    gl <- c(gl, rownames(this_res))
  }

}    

wb <- createWorkbook()

for (cell_type in names(res)){
  tmp <- res[[cell_type]]
  tmp$symbol <- as.character(tmp$symbol)
  sheetName<-cell_type
  addWorksheet(wb, sheetName)
  writeData(wb, sheetName, tmp)
  addFilter(wb, sheetName, row = 1, cols = 1:ncol(tmp))
  setColWidths(wb, sheetName, cols = 1:ncol(tmp), widths="auto")  
}

saveWorkbook(wb, file = paste(base, "Supplement/DE_analysis.xlsx", sep="/"), overwrite = T) 

##################################
#
# plot and save number DE genes
#
##################################

glN <- data.frame()
deg_list <- list()
for (name in names(res)) {
  
  tmp <- res[[name]][which(res[[name]]$padj < 0.01),]
  
  deg_list[[ paste( gsub( pattern = "_KO_WT", replacement = "", x = name), "pos", sep="_" ) ]] <- tmp$symbol[ which(tmp$log2FoldChange > 0)]
  deg_list[[ paste( gsub( pattern = "_KO_WT", replacement = "", x = name), "neg", sep="_" ) ]] <- tmp$symbol[ which(tmp$log2FoldChange < 0)]
  
  nUP <- length(which(tmp$log2FoldChange > 0))
  nDOWN <- length(which(tmp$log2FoldChange < 0)) * -1
  
  glN <- rbind( glN, data.frame(n = c(nUP, nDOWN), dir = c('up', 'down'), group = name) )
  
}

glN$dir <- factor(glN$dir, levels=c('up', 'down'))

nDEG_plot <- ggplot(glN, aes(x=dir, y=n, fill=dir)) + 
  geom_bar(stat='identity', position='dodge') + 
  scale_fill_manual(values=c('up' = 'red', 'down' = 'blue')) +
  ylab('number of DEG') + ylim(-120, 300) + facet_wrap(~group) +
  ylim(-200, 300) + theme_classic()

ggsave(filename = paste( figure_base, "/Fig3B_DE_gene_number.pdf", sep="/"), plot = nDEG_plot, width = 3 )

list_of_plots <- list()
list_of_plots[[1]] <- ggVennDiagram( deg_list[ c("P0_pos", "P4_pos") ] ) + coord_flip(clip = "off") +
  scale_fill_gradient(low="white",high = "white") + theme(legend.position = "NONE")
list_of_plots[[2]] <- ggVennDiagram( deg_list[ c("P0_neg", "P4_neg") ] ) + coord_flip(clip = "off") +
  scale_fill_gradient(low="white",high = "white") + theme(legend.position = "NONE")

venn_plot <- plot_grid( plotlist = list_of_plots, ncol = 1, greedy = F )
ggsave( filename = paste(figure_base, "/Fig3C_DEG_venn_plot.pdf", sep="/"), plot = venn_plot, width = 4 )

# What genes are in the intersection?
intersect( deg_list[["P0_pos"]], deg_list[["P4_pos"]] )
# [1] "Luzp2"  "S100a1" "Scrg1"  "Cnp"    "Rgcc"   "Rnd2" 
intersect( deg_list[["P0_neg"]], deg_list[["P4_neg"]] )
#  [1] "Mef2c"   "Ssbp2"   "Dync1i1" "Gas7"    "Nwd2"    "Syt5"    "Lmo4"    "Caly"    "Trim67"  "Ly6h"    "Spint2"  "Camkv" 

# ====================
# GO term enrichment
# ====================

eGO <- list()

for (name in names(res)) {
 tmp <- res[[name]][which( res[[name]]$padj < 0.01 ),]
 tmp$symbol <- as.character( tmp$symbol )
 
 gl_all <- bitr(geneID = tmp$symbol, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)
 gl_all <-  gl_all[which(!is.na(gl_all$ENTREZID)),]
 
 tmp_eGO <- enrichGO( gene = gl_all$ENTREZID,
                   OrgDb         = org.Mm.eg.db,
                   ont           = "BP",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.05,
                   qvalueCutoff  = 0.1,
                   readable      = TRUE,
                   pool=F)
  
 eGO[[ name ]] <- tmp_eGO@result
 
}

wb <- createWorkbook()
for (comparison in names(eGO)){
  tmp <- eGO[[comparison]]
  sheetName<-comparison
  addWorksheet(wb, sheetName)
  writeData(wb, sheetName, tmp)
  addFilter(wb, sheetName, row = 1, cols = 1:ncol(tmp))
  setColWidths(wb, sheetName, cols = 1:ncol(tmp), widths="auto")  
}
saveWorkbook(wb, file = paste(base, "Supplement/GO_analysis.xlsx", sep="/"), overwrite = T) 

# plot the DEGs connected to gliogenesis - 3rd ranking in most significantly enriched GO term
eGO[[1]]$Description[3]
#[1] "gliogenesis"
gliogenesis_genes <- strsplit( x = eGO[[1]]$geneID[3], split = "/")[[1]]

lfc_df2plot <- res[[1]][ which( res[[1]]$symbol %in% gliogenesis_genes ), c("symbol", "log2FoldChange") ]
lfc_df2plot <- lfc_df2plot[ order(lfc_df2plot$log2FoldChange, decreasing = F), ]
lfc_df2plot$symbol <- factor( lfc_df2plot$symbol, levels = lfc_df2plot$symbol)
lfc_df2plot$color <- "red"
lfc_df2plot$color[which(lfc_df2plot$log2FoldChange < 0)] <- "blue"

ggplot( lfc_df2plot, aes(y=symbol, x=log2FoldChange, fill = color)) + geom_bar(stat="identity") + theme_classic() +
  scale_fill_manual( values = c("red" = "red", "blue" = "blue")) + ggtitle("Gliogenesis P0") + theme(legend.position = "none")

ggsave( filename = paste(figure_base, "/Fig3D_lfc_plot_gliogenesis.pdf", sep="/"), width = 5, height = 5 )

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
#   [1] ggbeeswarm_0.7.2            ggVennDiagram_1.5.2         cowplot_1.1.3               reshape2_1.4.4             
# [5] org.Mm.eg.db_3.18.0         AnnotationDbi_1.64.1        clusterProfiler_4.10.1      ggplot2_3.5.0              
# [9] DESeq2_1.42.1               SummarizedExperiment_1.32.0 Biobase_2.62.0              MatrixGenerics_1.14.0      
# [13] matrixStats_1.2.0           GenomicRanges_1.54.1        GenomeInfoDb_1.38.8         IRanges_2.36.0             
# [17] S4Vectors_0.40.2            BiocGenerics_0.48.1         openxlsx_4.2.5.2           
# 
# loaded via a namespace (and not attached):
#   [1] RColorBrewer_1.1-3      rstudioapi_0.16.0       jsonlite_1.8.8          magrittr_2.0.3          farver_2.1.1           
# [6] fs_1.6.3                zlibbioc_1.48.2         ragg_1.3.0              vctrs_0.6.5             memoise_2.0.1          
# [11] RCurl_1.98-1.14         ggtree_3.10.1           S4Arrays_1.2.1          SparseArray_1.2.4       gridGraphics_0.5-1     
# [16] plyr_1.8.9              cachem_1.0.8            igraph_2.1.4            lifecycle_1.0.4         pkgconfig_2.0.3        
# [21] Matrix_1.6-5            R6_2.5.1                fastmap_1.1.1           gson_0.1.0              GenomeInfoDbData_1.2.11
# [26] digest_0.6.34           aplot_0.2.2             enrichplot_1.22.0       colorspace_2.1-0        patchwork_1.2.0        
# [31] textshaping_0.3.7       RSQLite_2.3.6           labeling_0.4.3          fansi_1.0.6             httr_1.4.7             
# [36] polyclip_1.10-6         abind_1.4-5             compiler_4.3.2          bit64_4.0.5             withr_3.0.0            
# [41] BiocParallel_1.36.0     viridis_0.6.5           DBI_1.2.2               ggforce_0.4.2           MASS_7.3-60            
# [46] DelayedArray_0.28.0     HDO.db_0.99.1           tools_4.3.2             vipor_0.4.7             beeswarm_0.4.0         
# [51] ape_5.8                 scatterpie_0.2.2        zip_2.3.1               glue_1.7.0              nlme_3.1-163           
# [56] GOSemSim_2.28.1         grid_4.3.2              shadowtext_0.1.3        fgsea_1.28.0            generics_0.1.3         
# [61] gtable_0.3.4            tidyr_1.3.1             data.table_1.15.0       tidygraph_1.3.1         utf8_1.2.4             
# [66] XVector_0.42.0          ggrepel_0.9.5           pillar_1.9.0            stringr_1.5.1           yulab.utils_0.1.4      
# [71] splines_4.3.2           dplyr_1.1.4             tweenr_2.0.3            treeio_1.26.0           lattice_0.21-9         
# [76] bit_4.0.5               tidyselect_1.2.1        GO.db_3.18.0            locfit_1.5-9.12         Biostrings_2.70.3      
# [81] gridExtra_2.3           graphlayouts_1.1.1      stringi_1.8.3           lazyeval_0.2.2          ggfun_0.1.4            
# [86] codetools_0.2-19        ggraph_2.2.1            tibble_3.2.1            qvalue_2.34.0           ggplotify_0.1.2        
# [91] cli_3.6.2               systemfonts_1.0.6       munsell_0.5.0           Rcpp_1.0.12             png_0.1-8              
# [96] parallel_4.3.2          blob_1.2.4              DOSE_3.28.2             bitops_1.0-7            viridisLite_0.4.2      
# [101] tidytree_0.4.6          scales_1.3.0            packrat_0.9.2           purrr_1.0.2             crayon_1.5.2           
# [106] rlang_1.1.3             fastmatch_1.1-4         KEGGREST_1.42.0   
