# Net Generator

# Load in libraries
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
source("C:/Users/patterw/Documents/GitHub/dlbcl-multiomics/dlbcl_multiomics_code/R/SmCCNetSource_forcelayout.R")
source("C:/Users/patterw/Documents/GitHub/dlbcl-multiomics/dlbcl_multiomics_code/trim_out_unlabelled_data.R")

# get the annotations and retrieve gene names
load("C:/Users/patterw/Box/data/DLBCL_multi_omics.rdata")
labelled_data <- trim_out_unlabelled_data()
FullAnnot = getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19) 
FullAnnot = FullAnnot[,c("Name","UCSC_RefGene_Name")]

# Set working directory and load in analysis data
setwd("C:/Users/patterw/Documents/Summer Research/A1_data")
load("SmCCNetWeightsA1.RData")
load("CorrMatrixA1.RData")

# create AbarLabel object
AbarLabel <- c(colnames(Abar))

# subset methylsites in AbarLabel
AnnotMethylSites <- subset(FullAnnot, FullAnnot$Name %in% AbarLabel)

# subset genes in AbarLabel
gene_labels <- subset(wk.gene, wk.gene$id %in% AbarLabel)

# create new AbarLabel
methyl_labels <- paste0("MS: ", AnnotMethylSites$UCSC_RefGene_Name)
gene_labels   <- paste0("G: ", gene_labels$gene)

AbarLabel <- c(testy, testy2)

# Cut out edges with weight less than edgeCut
edgeCut <- 0.4

# Produce gene networks from adjacency network
for(idx in 1:length(Modules)){
  filename <- paste0(CVDir, "Net_", idx, ".pdf")
  plotMultiOmicsNetwork(Abar = Abar, CorrMatrix = bigCor, 
                        multiOmicsModule = Modules, ModuleIdx = idx, P1 = p1, 
                        EdgeCut = edgeCut, FeatureLabel = AbarLabel, 
                        SaveFile = filename)
}

plotMultiOmicsNetwork(Abar = Abar, CorrMatrix = bigCor,
                      multiOmicsModule = Modules, ModuleIdx = 1, P1 = p1,
                      EdgeCut = 0.021, FeatureLabel = AbarLabel1,
                      VertexLabelCex = 0.5, AddCorrSign = TRUE,
                      NetLayout = "lgl", VertexSize = 0.2)
