#########################################################

# This script was developed to calculate mutual_exclusion 
# from MEGSA 

#########################################################

run_megsa <- function(df,detail=0){
  library(parallel)
  source("MEGSA.R")
  
  # mutationMat <- as.matrix(read.table(infile, header = TRUE, row.names = 1) != 0)
  ## Parallel processing step
  num_cores <-detectCores()-1
  cl <- makeCluster(num_cores)
  
  
  
  mutationMat <- as.matrix(df)
  # print(mutationMat)
  geneAll <- colnames(mutationMat)
  M <- ncol(mutationMat)
  tempCombn <- combn(M, 2)
  nPair <- ncol(tempCombn)
  pairS <- rep(NA, nPair); names(pairS) <- 1:nPair

  mx <- matrix(data=NA, nrow=nPair,ncol = 4)
  colnames(mx) <-list("gene1","gene2", "pvalue", "loglikelyhood")
  
  pb <- txtProgressBar(0, nPair, style = 3)
  for(i in 1:nPair){
    
    llk <-funEstimate(mutationMat[, tempCombn[, i]])$S
    pval <- 0.5 * pchisq(llk, 1, lower.tail = FALSE)

    g1 <-geneAll[tempCombn[, i]][1]
    g2 <-geneAll[tempCombn[, i]][2]
    
    mx[i,] <- c(g1,g2,pval,llk)
    if(detail)cat(i, "\t", geneAll[tempCombn[, i]],pval, "\n")
    setTxtProgressBar(pb, i)
  }
  
  return(data.frame(mx))
}
# 
# df <- read.table(file = 'mutationMat_LAML.txt', sep = '\t', header = TRUE, row.names = 1)
# 
# print(df)
# mx <- as.matrix(df)
# run_megsa(mx)