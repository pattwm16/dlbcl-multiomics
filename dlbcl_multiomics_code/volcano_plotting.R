# Generate volcano plots

library(EnhancedVolcano)
library(limma)

DRR <- which(wk.pheno$Status %in% c("Diagnostic","Relapsed"))
datDRR <- wk.pheno[DRR, ]
rnaDRR <- cpm.rna[ , DRR]

design <- model.matrix(~ Status + Patient, data = datDRR)
fitDRR <- lmFit(rnaDRR, design = design)
fitDRR <- eBayes(fitDRR)

# check the top genes using topTable()
# Default for topTable is Benjamini-Hochberg
top_genes_DRR <- topTable(fitDRR, coef = 2, number = Inf,
                          adjust.method = "BH")

q1 <- EnhancedVolcano(top_genes_DRR, 
                lab = rownames(top_genes_DRR), 
                x = "logFC", 
                y = "P.Value",
                title = NULL,
                subtitle = "Gene Expression")

DRDC <- which(wk.pheno$Status %in% c("Diagnostic","Cured"))
datDRDC <- wk.pheno[DRDC, ]
rnaDRDC <- cpm.rna[ , DRDC]

design <- model.matrix(~ Status, data = datDRDC) # no patient here since these are independent samples
fitDRDC <- lmFit(rnaDRDC, design = design)
fitDRDC <- eBayes(fitDRDC)

top_genes_DRDC <- topTable(fitDRDC, coef = 2, number = Inf,
                           adjust.method="BH")
p1 <- EnhancedVolcano(top_genes_DRDC, 
                lab = rownames(top_genes_DRDC), 
                x = "logFC", 
                y = "P.Value",
                title = NULL,
                subtitle = "Gene Expression")

### Methylation Sites -----
DRR <- which(wk.pheno$Status %in% c("Diagnostic","Relapsed"))
datDRR <- wk.pheno[DRR, ]
rnaDRR <- wk.methy[ , DRR]

design <- model.matrix(~ Status + Patient, data = datDRR)
fitDRR <- lmFit(rnaDRR, design = design)
fitDRR <- eBayes(fitDRR)

# check the top methylation sites using topTable()
top_sites_DRR <- topTable(fitDRR, coef = 2, number = Inf)
q2 <- EnhancedVolcano(top_sites_DRR, 
                lab = rownames(top_sites_DRR), 
                x = "logFC", 
                y = "P.Value",
                title = NULL,
                subtitle = "Methylation Sites")



DRDC <- which(wk.pheno$Status %in% c("Diagnostic","Cured"))
datDRDC <- wk.pheno[DRDC, ]
rnaDRDC <- wk.methy[ , DRDC]

design <- model.matrix(~ Status, data = datDRDC) # no patient here since these are independent samples
fitDRDC <- lmFit(rnaDRDC, design = design)
fitDRDC <- eBayes(fitDRDC)

# check the top methylation sites using topTable()
top_sites_DRDC <- topTable(fitDRDC, coef = 2, number = Inf)
p2 <- EnhancedVolcano(top_sites_DRDC, 
                lab = rownames(top_sites_DRDC), 
                x = "logFC", 
                y = "P.Value",
                title = NULL,
                subtitle = "Methylation Sites")
