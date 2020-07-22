# DLBCL Multi-omics Script File
# Author: Will Patterson

# Setup ----
# First import dependencies used in SmCCNet paper
library(PMA)
library(pbapply)
library(Matrix)
library(igraph)

library(beepr) # Added beepr to alert when loop completes

# Set seed for reproducibility (TODO: delete in production code)
set.seed(12345)

# Specify path to code/functions from SmCCNet paper
# Stored in Box drive, but path is dependent on OS
if(.Platform[1] == "windows"){
  load("C:/Users/patterw/Box/data/DLBCL_multi_omics.rdata")
  setwd("C:/Users/patterw/Box/summer_research/SmCCNet-master")
} else if(.Platform[1] == "unix"){ 
  load("~/Box/data/DLBCL_multi_omics.rdata")
  setwd("~/Box/summer_research/SmCCNet-master")
} else {
  print("Error: This OS is not supported. You will need to
        specify the path to the data manually.")
}

# Import essential Rscripts from SmCCNet project
source("R/ModifiedPMA.R")
source("R/SmCCNetSource.R")

# Load in data ----

# First Analysis (Diagnostic samples only) ----

# Find the sample values and patient numbers that have
# diagnostic status (should be 13)
diagnostic_pts <- subset(wk.pheno, Status == "Diagnostic")

# Select only diagnostic pt data for wk.gene and wk.methy
## Subsetting cpm.rna matrix
cpm.rna_diagnostic <- subset(cpm.rna, select = as.matrix(diagnostic_pts[1]))

## Subsetting wk.methy matrix
sample_idx_wk.methy <- str_extract(colnames(wk.methy), "NG\\d+") %in%
  as.matrix(diagnostic_pts[1])
wk.methy_diagnostic <- subset(wk.methy, select = sample_idx_wk.methy)


# Load in data as symbolic form
X1 <- t(cpm.rna_diagnostic)
X2 <- t(wk.methy_diagnostic)

# Apply SmCCNet functions to multi-omics data ----

# Key Hyperparameters (Start)
# - # of folds for cross val (K)
# - number of subsamples (SubsamplingNum)

# Define num of features for each variate + num of subjects
p1 <- ncol(X1)
p2 <- ncol(X2)
n <- nrow(X1)
AbarLabel <- c(colnames(cbind(X1, X2)))

# Determining optimal sparsity penalties through cross val ----
K <- 3 # num folds in k-fold cross val
CCcoef <- NULL # unweighted version of SmCCNet
s1 <- 0.7; s2 <- 0.9 # feature sampling proportions 
SubsamplingNum <- 1000 # num of subsamples

# Create sparsity penalty options.
pen1 <- seq(.05, .3, by = .05) 
pen2 <- seq(.05, .3, by = .05) 
P1P2 <- expand.grid(pen1, pen2) # Map (l1, l2) to (c1, c2).
c1 <- sqrt(p1 * s1) * P1P2[ , 1]; c1[c1] <- 1
c2 <- sqrt(p2 * s2) * P1P2[ , 2]; c2[c2 < 1] <- 1

# Based on prior knowledge we may assume that there are at least as many
# genes as miRNAs in each network.
# NOTE: This may not hold for current project. Check with Bo.
P1P2 <- P1P2[which(c1>c2), ]

# Set a CV directory.
CVDir <- paste(as.character(K), "foldCV", format(Sys.time(), "%b.%d.%Y/"), 
               sep="")
dir.create(CVDir)

# Split up data into test / train sets
foldIdx <- split(1:n, sample(1:n, K))
for(i in 1:K){
  iIdx <- foldIdx[[i]]
  x1.train <- scale(X1[-iIdx, ]) 
  x2.train <- scale(X2[-iIdx, ]) 
  x1.test <- scale(X1[iIdx, ]) 
  x2.test <- scale(X2[iIdx, ]) 

  # Check if standardized data sets are valid.
  if(is.na(min(min(x1.train), min(x2.train), min(x1.test), min(x2.test)))){
    stop("Invalid scaled data. At least one of the data matrices include a 
         column with zero variance.")
  }
  subD <- paste0(CVDir, "CV_", i, "/")
  dir.create(subD)
  save(x1.train, x2.train, x1.test, x2.test,
       s1, s2, P1P2, p1, p2, SubsamplingNum, CCcoef, 
       file = paste0(subD, "Data.RData"))
}


# Running the K-fold CV (Parallelized) ----
library(parallel)

cl <- makeCluster(K, type = "FORK") # Create K parallel threads.
clusterExport(cl = cl, "CVDir") # Pass on variable CVDir to each thread.
parSapply(cl, 1:K, function(CVidx){
  # Reload source code files for each thread. 
  source("R/ModifiedPMA.R")
  source("R/SmCCNetSource.R")
  
  # Create a result directory for each thread.
  subD <- paste0(CVDir, "CV_", CVidx, "/")
  load(paste0(subD, "Data.RData"))
  dir.create(paste0(subD, "SmCCA/"))
  
  RhoTrain <- RhoTest <- DeltaCor <- rep(0, nrow(P1P2))
  for(idx in 1:nrow(P1P2)){
    # Consider one pair of sparsity penalties at a time.
    l1 <- P1P2[idx, 1]
    l2 <- P1P2[idx, 2]
    
    # Run SmCCA on the subsamples (Figure 1, Step II)
    Ws <- getRobustPseudoWeights(x1.train, x2.train, NULL, l1, l2, 
                                 s1, s2, NoTrait = TRUE,
                                 FilterByTrait = FALSE, 
                                 SubsamplingNum = SubsamplingNum, 
                                 CCcoef = CCcoef)
    
    # Aggregate pseudo-canonical weights from the subsamples.
    meanW <- rowMeans(Ws)
    v <- meanW[1:p1]
    u <- meanW[p1 + 1:p2]
    
    # Compute the prediction error for given CV fold and sparsity penalties.
    if(is.null(CCcoef)){CCcoef <- rep(1, 3)} # Unweighted SmCCA.
    rho.train <- cor(x1.train %*% v, x2.train %*% u) * CCcoef[1] + 
      cor(x1.train %*% v, yy.train) * CCcoef[2] + 
      cor(x2.train %*% u, yy.train) * CCcoef[3]
    rho.test <- cor(x1.test %*% v, x2.test %*% u) * CCcoef[1] +
      cor(x1.test %*% v, yy.test) * CCcoef[2] + 
      cor(x2.test %*% u, yy.test) * CCcoef[3]
    RhoTrain[idx] <- round(rho.train, digits = 5)
    RhoTest[idx] <- round(rho.test, digits = 5)
    DeltaCor[idx] <- abs(rho.train - rho.test)
    
    # Periodically save results in a temporary file.
    if(idx %% 10 == 0){
      save(P1P2, RhoTrain, RhoTest, DeltaCor, idx, 
           file = paste0(subD, "temp.RData"))
    }
  }
  
  # Record prediction errors for given CV fold and all sparsity penalty 
  # options.
  DeltaCor.all <- cbind(P1P2, RhoTrain, RhoTest, DeltaCor)
  colnames(DeltaCor.all) <- c("l1", "l2", "Training CC", "Test CC", 
                              "CC Pred. Error")
  write.csv(DeltaCor.all, 
            file = paste0(subD, "SmCCA/PredictionError.csv"))
  
  # Remove the temporary file.
  system(paste0("rm ", subD, "temp.RData"))
  return(CVidx)
})

# Notify me that it completed
beep(8)

# Close cluster
stopCluster(cl)

# Extract penalty pair with the smallest total prediction error ----
# Combine prediction errors from all K folds and compute the total prediction
# error for each sparsity penalty pair.
testCC <- predError <- NULL
for(j in 1:K){
  resultT <- paste0(CVDir, "CV_", j, "/SmCCA/PredictionError.csv")
  dCorT <- read.csv(resultT)[ , -1]
  testCC <- cbind(testCC, abs(dCorT[ , 4]))
  predError <- cbind(predError, dCorT[ , 5])
}

S1 <- rowMeans(testCC)
S2 <- rowMeans(predError)
T12 <- dCorT[ , -3]; T12[ , 3] <- S1; T12[ , 4] <- S2
write.csv(T12, file = paste0(CVDir, "TotalPredictionError.csv"))

# Visualization ----
library(plotly)
library(reshape2)

f1 <- list(
  family = "Arial, sans-serif",
  size = 20,
  color = "black"
)
f2 <- list(
  family = "Old Standard TT, serif",
  size = 20,
  color = "black"
)
a <- list(
  title = "l1",
  titlefont = f1,
  showticklabels = TRUE,
  tickfont = f2
)
b <- list(
  title = "l2",
  titlefont = f1,
  showticklabels = TRUE,
  tickfont = f2
)
hmelt <- melt(T12[ , -3], id.vars = c("l1", "l2"))
contourPlot <- plot_ly(hmelt, x = ~l1, y = ~l2, z = ~value, type = "contour") %>%
  layout(xaxis = a, yaxis = b, showlegend = TRUE, legend = f1)  
export(contourPlot, file = paste0(CVDir, "TotalPredictionError.pdf"))

pen <- which(S2 == min(S2))
l1 <- T12$l1[pen]; l2 <- T12$l2[pen]
print(paste0("Optimal penalty pair (l1, l2): (", l1, ",", l2, ")"))

# Integrate two omics data types and a quantitative phenotype ----
Ws <- getRobustPseudoWeights(X1, X2, Y, l1, l2, s1, s2, 
                             NoTrait = FALSE, FilterByTrait = FALSE, 
                             SubsamplingNum = SubsamplingNum, CCcoef = CCcoef)
Abar <- getAbar(Ws, FeatureLabel = AbarLabel) # had to add "FeatureLabel =" here, not sure if kosher...

# Get multi-omics modules and plot subnetworks
Modules <- getMultiOmicsModules(Abar, p1)
save(Ws, Abar, Modules, file = paste0(CVDir, "SmCCNetWeights.RData"))


# Cut out edges with weight less than edgeCut
bigCor <- cor(cbind(X1, X2))
edgeCut <- 0.1
for(idx in 1:length(Modules)){
  filename <- paste0(CVDir, "Net_", idx, ".pdf")
  plotMultiOmicsNetwork(Abar = Abar, CorrMatrix = bigCor, 
                        multiOmicsModule = Modules, ModuleIdx = idx, P1 = p1, 
                        EdgeCut = edgeCut, FeatureLabel = AbarLabel,
                        SaveFile = filename)
}


