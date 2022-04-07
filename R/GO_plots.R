GO_plots <- function(sigGenes, name="name"){
  kk <- enrichGO(gene = sigGenes, ont="BP", OrgDb = "org.Hs.eg.db")
  GO <- kk2foldchange(kk)
  l <- list()
  l[[paste("GO_BP_", name)]] <- GO
  if(length(GO$ID)>30){
    GOsub <- GO[1:30,]
  }else{
    GOsub <- GO[1:length(GO$ID),]
  }
  GO_ordered <- GOsub[order(GOsub$foldchange),]
  GO_ordered$foldchange <- as.numeric(GO_ordered$foldchange)
  GO_ordered$Description<- factor(GO_ordered$Description, levels = GO_ordered$Description[order(GO_ordered$foldchange, decreasing = FALSE)])
  ### GO plot - biological process
  if(length(GO$ID>0)){
    p1 <- ggplot(GO_ordered, aes(foldchange, Description)) +
      geom_point(aes(color=p.adjust), size=6)+
      scale_color_gradient(low="#00B0F0", high="lightgray")+
      theme_minimal(base_size = 16)+xlab("fold enrichment")+ylab("")+labs(color = "p.adjust")
  }else{
    p1 <- paste("No significant GO terms for Biological Process")
  }

  ### GO plot - biological process
  kk <- enrichGO(gene = sigGenes, ont="MF", OrgDb = "org.Hs.eg.db")
  GO <- kk2foldchange(kk)
  l[[paste("GO_MF_", name)]] <- GO
  if(length(GO$ID)>30){
    GOsub <- GO[1:30,]
  }else{
    GOsub <- GO[]
  }
  GO_ordered <- GOsub[order(GOsub$foldchange),]
  GO_ordered$foldchange <- as.numeric(GO_ordered$foldchange)
  GO_ordered$Description<- factor(GO_ordered$Description, levels = GO_ordered$Description[order(GO_ordered$foldchange, decreasing = FALSE)])

  if(length(GO$ID>0)){
    p2 <- ggplot(GO_ordered, aes(foldchange, Description)) +
      geom_point(aes(color=p.adjust), size=6)+
      scale_color_gradient(low="#00B0F0", high="lightgray")+
      theme_minimal(base_size = 16)+xlab("fold enrichment")+ylab("")+labs(color = "p.adjust")
  }else{
    p2 <- paste("No significant GO terms for Molecular Function")
  }

  ### GO plot - cellular component
  kk <- enrichGO(gene = sigGenes, ont="CC", OrgDb = "org.Hs.eg.db")
  GO <- kk2foldchange(kk)
  l[[paste("GO_CC_", name)]] <- GO
  if(length(GO$ID)>30){
    GOsub <- GO[1:30,]
  }else{
    GOsub <- GO[]
  }
  GO_ordered <- GOsub[order(GOsub$foldchange),]
  GO_ordered$foldchange <- as.numeric(GO_ordered$foldchange)
  GO_ordered$Description<- factor(GO_ordered$Description, levels = GO_ordered$Description[order(GO_ordered$foldchange, decreasing = FALSE)])
  if(length(GO$ID>0)){
    p3 <- ggplot(GO_ordered, aes(foldchange, Description)) +
      geom_point(aes(color=p.adjust), size=6)+
      scale_color_gradient(low="#00B0F0", high="lightgray")+
      theme_minimal(base_size = 16)+xlab("fold enrichment")+ylab("")+labs(color = "p.adjust")
  }else{
    p3 <- paste("No significant GO terms for Cellular Component")
  }
  plots <- list()
  # if("p1" %in%  ls()){
  #   plots[[1]] <- p1
  # }
  # if("p2" %in% ls()){
  #   plots[[2]] <- p2
  # }
  # if("p3" %in% ls()){
  #   plots[[3]] <- p3
  # }
  plots=list(BP=p1,MF=p2,CC=p3)
  my_outs <- list(l, plots)
  return(my_outs)
}


















