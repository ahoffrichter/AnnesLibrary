kk2foldchange <- function(kk){
  GO <- as.data.frame(kk)
  GeneRatio <- strsplit(GO$GeneRatio,"/")
  GRdf <- data.frame()
  for (i in seq_along(GeneRatio)){
    GRdf[i,1] <- GeneRatio[[i]][1]
    GRdf[i,2] <- GeneRatio[[i]][2]
  }
  GRdf$V1 <- as.numeric(GRdf$V1)
  GRdf$V2 <- as.numeric(GRdf$V2)
  GO$GeneRatio2 <- GRdf$V1/GRdf$V2
  BgRatio <- strsplit(GO$BgRatio,"/")
  BgRdf <- data.frame()
  for (i in seq_along(BgRatio)){
    BgRdf[i,1] <- BgRatio[[i]][1]
    BgRdf[i,2] <- BgRatio[[i]][2]
  }
  BgRdf$V1 <- as.numeric(BgRdf$V1)
  BgRdf$V2 <- as.numeric(BgRdf$V2)
  GO$BgRatio2 <- BgRdf$V1/BgRdf$V2
  GO$foldchange <- GO$GeneRatio2/GO$BgRatio2
  return(GO)
}
