methyl_gene_plot_batch <- function(gene, chrom, start, end){
  minbase <- start - (0.25*(end-start))
  maxbase <- end + (0.25*(end-start))
  iTrack <- IdeogramTrack(genome = gen, chromosome = chrom, name="")
  gTrack <- GenomeAxisTrack(col="black", cex=1, name="", fontcolor="black")
  rTrack <- UcscTrack(genome=gen, chromosome=chrom, track="NCBI RefSeq",
                      from=minbase, to=maxbase, trackType="GeneRegionTrack",
                      rstarts="exonStarts", rends="exonEnds", gene="name",
                      symbol="name2", transcript="name", strand="strand",
                      fill="darkblue",stacking="squish", name="RefSeq",
                      showId=TRUE, geneSymbol=TRUE)
  annEPICOrd <- annEPICSub[order(annEPICSub$chr,annEPICSub$pos),]
  bValsOrd <- bVals[match(annEPICOrd$Name,rownames(bVals)),]
  cpgData <- GRanges(seqnames=Rle(annEPICOrd$chr),
                     ranges=IRanges(start=annEPICOrd$pos, end=annEPICOrd$pos),
                     strand=Rle(rep("*",nrow(annEPICOrd))),
                     betas=bValsOrd)
  methTrack <- DataTrack(range=cpgData, groups=targets$Sample_Group,genome = gen,
                         chromosome=chrom, ylim=c(-0.05,1.05), col=pal,
                         type=c("a","p"), name="DNA Meth. cond.\n(beta value)",
                         background.panel="white", legend=TRUE, cex.title=0.8,
                         cex.axis=0.8, cex.legend=0.8)
  methTrack_batch <- DataTrack(range=cpgData, groups=targets$batch,genome = gen,
                         chromosome=chrom, ylim=c(-0.05,1.05), col=pal,
                         type=c("a","p"), name="DNA Meth.batch\n(beta value)",
                         background.panel="white", legend=TRUE, cex.title=0.8,
                         cex.axis=0.8, cex.legend=0.8)
  dmrTrack <- AnnotationTrack(start=start, end=end, genome=gen, name=gene,
                              chromosome=chrom,fill="darkred")
  tracks <- list(iTrack, gTrack, methTrack,methTrack_batch, dmrTrack, rTrack)
  sizes <- c(1,1,7,7,1,1) # set up the relative sizes of the tracks
  plot <- plotTracks(tracks, from=minbase, to=maxbase, showTitle=TRUE, add53=TRUE,
             add35=TRUE, grid=TRUE, lty.grid=3, sizes = sizes, length(tracks))
  return(plot)
}
