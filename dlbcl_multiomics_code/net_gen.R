# Net Generator

anal <- 2

# Load in libraries
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
source("C:/Users/patterw/Documents/GitHub/dlbcl-multiomics/dlbcl_multiomics_code/R/SmCCNetSource.R")
source("C:/Users/patterw/Documents/GitHub/dlbcl-multiomics/dlbcl_multiomics_code/trim_out_unlabelled_data.R")

# get the annotations and retrieve gene names
load("C:/Users/patterw/Box/data/DLBCL_multi_omics.rdata")
labelled_data <- trim_out_unlabelled_data()
FullAnnot = getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19) 
FullAnnot = FullAnnot[,c("Name","UCSC_RefGene_Name")]

# Set working directory and load in analysis data
if(anal == 1){
  setwd("C:/Users/patterw/Documents/Summer Research/A1_data")
  load("SmCCNetWeightsA1.RData")
  load("CorrMatrixA1.RData")
  print("Analysis 1")
} else if(anal == 2){
  setwd("C:/Users/patterw/Documents/Summer Research/A2_data")
  load("SmCCNetWeightsA2.RData")
  load("CorrMatrixA2.RData")
  print("Analysis 2")
} else{
  print("Error: You have not set an analysis type")
}


# create AbarLabel object
AbarLabel <- c(colnames(Abar))

# subset methylsites in AbarLabel
AnnotMethylSites <- subset(FullAnnot, FullAnnot$Name %in% AbarLabel)

# subset genes in AbarLabel
gene_labels <- subset(wk.gene, wk.gene$id %in% AbarLabel)

# create new AbarLabel
methyl_labels <- paste0("MS: ", AnnotMethylSites$UCSC_RefGene_Name)
gene_labels   <- paste0("G: ", gene_labels$gene)

AbarLabel <- c(methyl_labels, gene_labels)

# Cut out edges with weight less than edgeCut
edgeCut <- 0.21

# Produce gene networks from adjacency network
for(idx in 1:length(Modules)){
  filename <- paste0("Net_", idx, ".pdf")
  plotMultiOmicsNetwork(Abar = Abar, CorrMatrix = bigCor, 
                        multiOmicsModule = Modules, ModuleIdx = idx, P1 = p1, 
                        EdgeCut = edgeCut, FeatureLabel = AbarLabel,
                        VertexLabelCex = 0.5, NetLayout = "lgl", 
                        VertexSize = 0.2, SaveFile = filename)
}

# get list of genes in network
for(idx in 1:length(Modules)){
  grp <- Modules[[idx]]
  grp.memb <- colnames(Abar)[grp]
  M.node <- grp.memb
  
  M <- as.matrix(Abar[M.node, M.node])
  M <- M * sign(bigCor[M.node, M.node])
  which(abs(M) >= 0.88, arr.ind = TRUE)
  net_verts <- unique(rownames(which(abs(M) >= edgeCut, arr.ind = TRUE)))
  
  genes_in_net <- subset(wk.gene$gene, wk.gene$id %in% net_verts)
  methylsites_in_net <- subset(FullAnnot$UCSC_RefGene_Name, 
                               FullAnnot$Name %in% net_verts)
  verts_list <- c(genes_in_net, methylsites_in_net)
  write.csv("Net_", idx, ".csv")
}
