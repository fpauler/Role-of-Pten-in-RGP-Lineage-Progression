# analyses a list of DEGs with iDEA
# arguments from command line:

args = commandArgs(trailingOnly=TRUE)

if (length(args) < 3) {
	error ( "expected are 3 command line arguments ( in that exact order) : [annotation RDS] [DEG list RDS] [prefic char]" )
}

annoRDS <- args[1]
degRDS <- args[2]
prefix = args[3]

library (iDEA)

# a data.frame with gene names as rows and term names as columns
# cell have a 1 if a gene is associated with a trem or 0 if not
go_anno <- readRDS(file = annoRDS )

# this is a list of lists
# first level is gnotype, second level is t-type
de_res <- readRDS(file = degRDS )

idea_GO_results <- list()

for (genotype in names(de_res)) {

  idea_GO_results[[genotype]] <- list()

  for (cluster in  names(de_res[[genotype]])) {

    idea_GO_results[[genotype]][[cluster]] <- list()    
    
    for (dir in c("pos", "neg", "all")) {

       message(cluster, " ", genotype, " ", dir)

       ## Assume you have obtained the DE results from i.e. zingeR, edgeR or MAST with the data frame res_DE (column: pvalue and LogFC)
       if (dir == "pos") {
          deg_df <- de_res[[genotype]][[cluster]][ which(de_res[[genotype]][[cluster]]$avg_log2FC > 0), ]
       } else if (dir == "neg") {
          deg_df <- de_res[[genotype]][[cluster]][ which(de_res[[genotype]][[cluster]]$avg_log2FC < 0), ]
       } else {
          deg_df <- de_res[[genotype]][[cluster]]
       }

       pvalue <- deg_df$p_val #### the pvalue column
       zscore <- qnorm(pvalue/2.0, lower.tail=FALSE) #### convert the pvalue to z-score
       beta <- deg_df$avg_log2FC ## effect size
       se_beta <- abs(beta/zscore) ## to approximate the standard error of beta
       beta_var = se_beta^2  ### square 
       summary = data.frame(beta = beta,beta_var = beta_var)
       ## add the gene names as the rownames of summary
       rownames(summary) = deg_df$gene ### or the gene id column in the res_DE results
    
       idea <- CreateiDEAObject(summary, go_anno, max_var_beta = 100, min_precent_annot = 0.0025, num_core=20)
    
       idea <- iDEA.fit(idea,
                     fit_noGS=FALSE,
                     init_beta=NULL, 
                     init_tau=c(-2,0.5),
                     min_degene=5,
                     em_iter=15,
                     mcmc_iter=1000, 
                     fit.tol=1e-5,
                     modelVariant = F,
                     verbose=TRUE)

       # deal with potential empty list entries that cause errors for BMA calculations
       empty_list_ids <- which( sapply( idea@de, function (x) length(x)) == 0)
       filled_list_ids <- which( sapply( idea@de, function (x) length(x)) > 0)
       
       idea@de[ empty_list_ids ] = NULL
       idea@annot_id <- idea@annot_id[ filled_list_ids ]
       idea@annotation <- idea@annotation[ filled_list_ids ]
    
       idea <- iDEA.louis(idea)
       idea <- iDEA.BMA(idea)
    
       idea_GO_results[[genotype]][[cluster]][[dir]] <- idea
    }
  }
}

saveRDS( object = idea_GO_results, file = paste(prefix,"idea_results.RDS", sep="." ) )

sessionInfo()
