# Net Generator

# Helper Functions ----
# Function: labelify
# Consumes a list of gene names (separated by ;), returns unique gene names
# with repeats spliced out (e.g. ["ZCCHC7", "ZCCHC7", "ZCCHC7"] returns 
# "ZCCHC7", but ["MIR219A2", "MIR1268A"] returns "MIR219A2;MIR1268A").
labelify <- function(names_list){
  lst <- c()
  for(elem in names_list){
    lst <- append(lst, paste(unique(unlist(strsplit(elem, ";"))), 
                             collapse = ";"))
  }
  return(lst)
}

# Analysis body ----

# Set which analysis is being performed
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
AnnotMethylSites$UCSC_RefGene_Name <- labelify(AnnotMethylSites$UCSC_RefGene_Name)

# subset genes in AbarLabel
gene_labels <- subset(wk.gene, wk.gene$id %in% AbarLabel)

# create new AbarLabel
methyl_labels <- paste0("MS: ", AnnotMethylSites$UCSC_RefGene_Name)
gene_labels   <- paste0("G: ", gene_labels$gene)

AbarLabel <- c(methyl_labels, gene_labels)

# Cut out edges with weight less than edgeCut
if(anal == 1){
  edgeCut <- 0.021
}else if(anal == 2){
  edgeCut <- 0.21
}


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
# TODO: Figure out why some genes shown in GeneList are not represented in
#       the networks produced by plotMultiomicsNetwork()
for(idx in 1:length(Modules)){
  # Copied next 5 lines from SmCCNetSource.R
  grp <- Modules[[idx]]
  grp.memb <- colnames(Abar)[grp]
  M.node <- grp.memb
  
  M <- as.matrix(Abar[M.node, M.node])
  M <- M * sign(bigCor[M.node, M.node])
  
  # Select the unique rows (since matrix is symmetric wrt names)
  net_verts <- unique(rownames(which(abs(M) >= edgeCut, arr.ind = TRUE)))
  
  # Select gene & methyl site names based on ids in net_verts
  genes_in_net <- subset(wk.gene$gene, wk.gene$id %in% net_verts)
  methylsites_in_net <- subset(AnnotMethylSites$UCSC_RefGene_Name, 
                               AnnotMethylSites$Name %in% net_verts)
  
  # Create a dataframe of these names and tidy
  verts_list <- as.data.frame(c(genes_in_net, methylsites_in_net))
  colnames(verts_list) <- c("GeneNames")
  
  # Check membership of these 
  TCGA <- read.csv("C:/Users/patterw/Documents/Summer Research/A1_data/publishedGenes.csv")
  verts_list$inTCGA <- verts_list$GeneNames %in% TCGA[,]
  
  # If not an empty graph, export list of genes as CSV
  if(nrow(verts_list >= 1)){
    write.csv(x = verts_list, 
              file = paste0("Net_", idx, "GeneList.csv")) 
  }
}
