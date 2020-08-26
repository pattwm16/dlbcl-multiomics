trim_out_unlabelled_data <- function(){
  # Gene trimming ----
  # collect the unlabelled genes
  unlabelled_genes <- wk.gene[which(wk.gene[,2] == ""),]
  
  # subset labelled genes
  labelled_genes <- subset(wk.gene, !wk.gene$id %in% unlabelled_genes$id)
  
  
  
  # Methylation Trimming ----
  # load in annotations from package
  library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
  
  # get the annotations and retrieve gene names
  FullAnnot = getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19) 
  FullAnnot = FullAnnot[,c("Name","UCSC_RefGene_Name")]
  
  # collect the unlabelled methyl sites
  unlabelled_methylsites <- FullAnnot[which(FullAnnot[,2] == ""),]
  
  # subset labelled genes
  labelled_methylsites <- subset(wk.methy, !rownames(wk.methy) %in% 
                                              unlabelled_methylsites$Name)
  
  
  return(list("Genes" = labelled_genes, "MethylSites" = labelled_methylsites))
}