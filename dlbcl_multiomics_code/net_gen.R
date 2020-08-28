# Net Generator

# Load in libraries
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)

source("R/SmCCNetSource_forcelayout.R")


# Load in data
load("SmCCNetWeights.RData")
load("CorrMatrix.RData")

# Cut out edges with weight less than edgeCut
edgeCut <- 0.4

# Add gene names to wk.methy
# TODO: Refactor to work with AbarLabel
row.names(wk.methy) <- subset(FullAnnot, 
                              FullAnnot$Name %in% row.names(wk.methy))$UCSC_RefGene_Name

# Before producing nets, use wk.gene to
# replace row/col names with Gene IDs
if(length(setdiff(colnames(Abar), rownames(Abar))) == 0){
  # There are ~556 NA values (and many other blank string values)
  # Replace NA with "" values
  wk.gene[which(is.na(wk.gene$gene)),]$gene <- ""
  # If gene ID is empty string, use the Ensemble ID
  wk.gene[which(wk.gene[,2] == ""),2] <- wk.gene[which(wk.gene[,2] == ""),1]
  revised_labels <- dplyr::recode(
    AbarLabel, 
    !!!setNames(as.character(wk.gene$gene), wk.gene$id))
}

# Produce gene networks from adjacency network
for(idx in 1:length(Modules)){
  filename <- paste0(CVDir, "Net_", idx, ".pdf")
  plotMultiOmicsNetwork(Abar = Abar, CorrMatrix = bigCor, 
                        multiOmicsModule = Modules, ModuleIdx = idx, P1 = p1, 
                        EdgeCut = edgeCut, FeatureLabel = AbarLabel, 
                        SaveFile = filename)
}


# TRASH ----
# Before producing nets, use wk.gene to
# replace row/col names with Gene IDs
if(length(setdiff(colnames(Abar), rownames(Abar))) == 0){
  # There are ~556 NA values (and many other blank string values)
  # TODO: Ask Bo about difference between NA vs. "" values in wk.gene
  
  # Replace NA with "" values
  wk.gene[which(is.na(wk.gene$gene)),]$gene <- ""
  # If gene ID is empty string, use the Ensemble ID
  wk.gene[which(wk.gene[,2] == ""),2] <- wk.gene[which(wk.gene[,2] == ""),1]
  revised_labels <- dplyr::recode(
    AbarLabel, 
    !!!setNames(as.character(wk.gene$gene), wk.gene$id))
}

# Produce gene networks from adjacency network
for(idx in 1:length(Modules)){
  filename <- paste0(CVDir, "Net_", idx, ".pdf")
  plotMultiOmicsNetwork(Abar = Abar, CorrMatrix = bigCor, 
                        multiOmicsModule = Modules, ModuleIdx = idx, P1 = p1, 
                        EdgeCut = edgeCut, FeatureLabel = revised_labels, SaveFile = filename)
}

plotMultiOmicsNetwork(Abar = Abar, CorrMatrix = bigCor,
                      multiOmicsModule = Modules, ModuleIdx = 1, P1 = p1,
                      EdgeCut = edgeCut, FeatureLabel = revised_labels,
                      VertexLabelCex = 0.2)
