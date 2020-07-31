# Differential Expression Pipe
# Objective: Reduce dimensionality of RNAseq data by selecting subset of 
#            differentially expressed/methylated genes 

# If visualizations are desired, set see_viz = TRUE
see_viz <- FALSE

diff_exp_dimreduction <-function(analysis){

# Import libraries
library(limma)

# Load in data ----
if(.Platform[1] == "windows"){
  load("C:/Users/patterw/Box/data/DLBCL_multi_omics.rdata")
  setwd("C:/Users/patterw/Box/summer_research/SmCCNet-master")
} else if(.Platform[1] == "unix"){ 
  load("~/Box/data/DLBCL_multi_omics.rdata")
  setwd("~/Desktop/CCLCM/Research/dlbcl-multiomics/diff_exp_dimreduction")
} else {
  print("Error: This OS is not supported. You will need to
        specify the path to the data manually.")
}

# Perform differential expression analysis ----
## choose the corresponding data for DR and R
DRR <- which(wk.pheno$Status %in% c("Diagnostic","Relapsed"))
datDRR <- wk.pheno[DRR, ]
rnaDRR <- cpm.rna[ , DRR]

design <- model.matrix(~ Status + Patient, data = datDRR)
fitDRR <- lmFit(rnaDRR, design = design)
fitDRR <- eBayes(fitDRR)

# check the top genes using topTable()
top_genes_DRR <- topTable(fitDRR, coef = 2, number = Inf)
top_genes_DRR <- subset(top_genes_DRR, P.Value <= 0.05)

## For DR vs. DC
DRDC <- which(wk.pheno$Status %in% c("Diagnostic","Cured"))
datDRDC <- wk.pheno[DRDC, ]
rnaDRDC <- cpm.rna[ , DRDC]

design <- model.matrix(~ Status, data = datDRDC) # no patient here since these are independent samples
fitDRDC <- lmFit(rnaDRDC, design = design)
fitDRDC <- eBayes(fitDRDC)

top_genes_DRDC <- topTable(fitDRDC, coef = 2, number = Inf)
top_genes_DRDC <- subset(top_genes_DRDC, P.Value <= 0.05)


# Perform differential methylation analysis ----
## choose the corresponding data for DR and R
DRR <- which(wk.pheno$Status %in% c("Diagnostic","Relapsed"))
datDRR <- wk.pheno[DRR, ]
rnaDRR <- wk.methy[ , DRR]

design <- model.matrix(~ Status + Patient, data = datDRR)
fitDRR <- lmFit(rnaDRR, design = design)
fitDRR <- eBayes(fitDRR)

# check the top methylation sites using topTable()
top_sites_DRR <- topTable(fitDRR, coef = 2, number = Inf)
top_sites_DRR <- subset(top_sites_DRR, P.Value <= 0.05)


## For DR vs. DC
DRDC <- which(wk.pheno$Status %in% c("Diagnostic","Cured"))
datDRDC <- wk.pheno[DRDC, ]
rnaDRDC <- wk.methy[ , DRDC]

design <- model.matrix(~ Status, data = datDRDC) # no patient here since these are independent samples
fitDRDC <- lmFit(rnaDRDC, design = design)
fitDRDC <- eBayes(fitDRDC)

# check the top methylation sites using topTable()
top_sites_DRDC <- topTable(fitDRDC, coef = 2, p.value = 0.05, number = Inf)

# Return only the differential expresion data for the analysis
if(analysis == 1){
  return(list("TG_DRR" = top_genes_DRR, "TS_DRR" = top_sites_DRR))
} else if(analysis == 2){}
  return(list("TG_DRR" = top_genes_DRR, "TG_DRDC" = top_genes_DRDC, 
              "TS_DRR" = top_sites_DRR, "TS_DRDC" = top_sites_DRDC))
}

# Visualization ----
if (see_viz){
  a <- ggplot(top_genes_DRR, aes(x = P.Value)) + geom_histogram(color="black", fill="lightblue") +
    labs(title="Significance distribution (Expression, DRR)",x="Unadjusted p-value", y = "Count") + theme_minimal()
  b <- ggplot(top_genes_DRDC, aes(x = P.Value)) + geom_histogram(color="black", fill="lightblue") +
    labs(title="Significance distribution (Expression, DRDC)",x="Unadjusted p-value", y = "Count") + theme_minimal()
  c <- ggplot(top_sites_DRR, aes(x = P.Value)) + geom_histogram(color="black", fill="lightblue") +
    labs(title="Significance distribution (Methylation, DRR)",x="Unadjusted p-value", y = "Count") + theme_minimal()
  d <- ggplot(top_sites_DRDC, aes(x = P.Value)) + geom_histogram(color="black", fill="lightblue") +
    labs(title="Significance distribution (Methylation, DRDC)",x="Unadjusted p-value", y = "Count") + theme_minimal()
  
  plot_grid(a, b)
  plot_grid(c, d)
  
  
  test_alpha_vals <- function(gene_df){
    alpha_thresh_vals <- c(0.01, 0.05, 0.1)
    num_sig_genes <- c()
    for (alpha in alpha_thresh_vals){
      num_sig_genes <- append(num_sig_genes, nrow(subset(gene_df, P.Value <= alpha)))
    }
    return(data.frame("Counts" = num_sig_genes, "Alpha" = alpha_thresh_vals, 
                      "Comparison" = rep(deparse(substitute(gene_df)), length(alpha_thresh_vals))))
  }
  
  counts_geneDRR <- test_alpha_vals(top_genes_DRR)
  counts_geneDRDC <- test_alpha_vals(top_genes_DRDC) 
  counts_siteDRR <- test_alpha_vals(top_sites_DRR) 
  counts_siteDRDC <- test_alpha_vals(top_sites_DRDC) 
  alpha_counts_genes <- rbind(counts_geneDRR, counts_geneDRDC)
  alpha_counts_sites <- rbind(counts_siteDRR, counts_siteDRDC)
  
  ggplot(alpha_counts_genes, aes(x=as.character(Alpha), y=Counts, fill=Comparison)) +
    facet_grid(cols=vars(Comparison)) +
    geom_col() +
    geom_text(aes(label = Counts), vjust = -0.5) +
    labs(title = "Significance Counts for Genes by Comparison", x = "\u03b1",
         y = "Count less than \u03b1")
  
  ggplot(alpha_counts_sites, aes(x=as.character(Alpha), y=Counts, fill=Comparison)) +
    facet_grid(cols=vars(Comparison)) +
    geom_col() +
    geom_text(aes(label = Counts), vjust = -0.5) +
    labs(title = "Significance Counts for CpG Sites by Comparison", x = "\u03b1",
         y = "Count less than \u03b1")
}

