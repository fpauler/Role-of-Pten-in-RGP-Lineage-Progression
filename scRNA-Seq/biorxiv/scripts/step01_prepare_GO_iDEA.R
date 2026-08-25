# ================================================================
# first step in the analysis: prepare GO terms for iDEA analysis
# ================================================================

library (GO.db)
library (data.table)
library (org.Mm.eg.db)
library (openxlsx)

set.seed(2401)

# define the base working folder
base <- "./scRNA-Seq/biorxiv"

# for iDEA analysis I have to prepare GO term <-> gene association:
# first define some search terms
GO_list <- read.csv(file = paste(base, "/other_data/Pip_Terms.csv", sep=""))

# cancer GO terms from here:
# https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-021-04105-8#availability-of-data-and-materials
# Additional file 9 has the consensus

consensus_terms <- read.xlsx( xlsxFile = paste(base, "other_data/Additional file 9.xlsx", sep="") )

# look up GO terms in my annotation
go_list <- Term(GOTERM)[ gsub( pattern = " ", replacement = "", x = consensus_terms$GO.terms ) ]
go_list <- go_list[ which(!is.na(go_list)) ]

# sanity check that search terms actually return GO terms
includeIDX <- sapply( GO_list$Include[which( nchar(GO_list$Include) > 0) ], function (y) {
  which(grepl(pattern =  y, x = Term(GOTERM), ignore.case = F, fixed = T ))
})

includeIDX <- c( unique( unlist(includeIDX) ), which( names(Term(GOTERM)) %in% names( go_list ) ) )
includeIDX <- unique( includeIDX )
                 
excludeIDX <- sapply( GO_list$Exclude[which( nchar(GO_list$Exclude) > 0) ], function (y) {
  which(grepl(pattern =  y, x = Term(GOTERM), ignore.case = F, fixed = T ))
})

term2id <- Term(GOTERM)[ setdiff( unique( unlist(includeIDX) ), unique( unlist(excludeIDX) ) ) ]
saveRDS( object = term2id, file = paste(base, "/iDEA_analysis/term2id.RDS", sep = "") )

lines_processed <- 0
#select GO terms of interest for iDEA analysis

# obtain all GO terms
go_id <- GOID( GOTERM[ setdiff( unique( unlist(includeIDX) ), unique( unlist(excludeIDX) ) ) ] )

# get genes
GO_GENE_DF <- data.frame()
for (i in 1:length(go_id)) {
  
  if (i %% 100 == 0){
    lines_processed <- lines_processed + 1
    message (lines_processed, " processed", dim(GO_GENE_DF))
  }
  allegs <- tryCatch( get(go_id[i], org.Mm.egGO2ALLEGS), 
                      error = function (error) {
                        # message(error, "\n")
                        return (NA)
                      }  )
  
  # get symbols
  if (!any(is.na(allegs))) {
    genes <- unique(unlist(mget(allegs,org.Mm.egSYMBOL)))
    tmp <- data.frame(gene = genes, term = 1)
    
    term_name <- Term(GOTERM)[[go_id[i]]]
    term_name <- gsub(pattern = " ", replacement = "_", x = term_name)
    
    colnames(tmp) <- c( "gene", term_name )
    
    if (nrow(tmp) > 10 & nrow(tmp < 1000)) {
      if ( nrow(GO_GENE_DF) == 0 ) {
        GO_GENE_DF <- as.data.table(tmp, key = "gene")
      } else {
        tmp <- as.data.table(tmp, key = "gene")
        GO_GENE_DF <- merge(GO_GENE_DF, tmp, all=T)
      }
      
    } else {
      message("removed ", term_name, " with ", nrow(tmp), " genes")
    }
  }
}

GO_GENE_DF <- as.data.frame(GO_GENE_DF)
rownames(GO_GENE_DF) <- GO_GENE_DF$gene
GO_GENE_DF <- GO_GENE_DF[,c(2:ncol(GO_GENE_DF))]
GO_GENE_DF[is.na(GO_GENE_DF)] <- 0

# save the data to do the iDEA analysis 
saveRDS(object = GO_GENE_DF, file = paste(base, "/iDEA_analysis/GO_anno_iDEA.RDS", sep = ""))

# sessionInfo()
# 
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
#   [1] openxlsx_4.2.5.2     org.Mm.eg.db_3.18.0  data.table_1.15.0    GO.db_3.18.0         AnnotationDbi_1.64.1 IRanges_2.36.0       S4Vectors_0.40.2     Biobase_2.62.0      
# [9] BiocGenerics_0.48.1 
# 
# loaded via a namespace (and not attached):
#   [1] zip_2.3.1               crayon_1.5.2            vctrs_0.6.5             httr_1.4.7              cli_3.6.2               rlang_1.1.3             packrat_0.9.2          
# [8] stringi_1.8.3           DBI_1.2.2               png_0.1-8               bit_4.0.5               RCurl_1.98-1.14         Biostrings_2.70.3       KEGGREST_1.42.0        
# [15] bitops_1.0-7            fastmap_1.1.1           GenomeInfoDb_1.38.8     memoise_2.0.1           compiler_4.3.2          RSQLite_2.3.6           blob_1.2.4             
# [22] Rcpp_1.0.12             pkgconfig_2.0.3         XVector_0.42.0          rstudioapi_0.16.0       R6_2.5.1                GenomeInfoDbData_1.2.11 tools_4.3.2            
# [29] bit64_4.0.5             zlibbioc_1.48.2         cachem_1.0.8 
