library(RnaBindingProtein)
source("~/github/fimoService/batchMode/fimoBatchTools.R")
#----------------------------------------------------------------------------------------------------
toTrackTable <- function(tbl.fimo)
{
    browser()
    pval <- tbl.fimo$p.value
    score <- -log10(pval)
    strand <- tbl.fimo$strand
    motif <- tbl.fimo$motif_id

    loc.strings <- tbl.fimo$sequence_name
    starts <- tbl.fimo$start
    ends <- tbl.fimo$stop
    chroms <- unlist(lapply(strsplit(loc.strings, ":"), "[", 1))
    start.end.strings <- unlist(lapply(strsplit(loc.strings, ":"), "[", 2))
    region.starts <- as.numeric(unlist(lapply(strsplit(start.end.strings, "-"), "[", 1)))
    region.ends <- as.numeric(unlist(lapply(strsplit(start.end.strings, "-"), "[", 2)))
    bs.start <- region.starts + starts
    bs.end <- region.starts + ends

    data.frame(chrom=chroms, start=bs.start, end=bs.end, score=score, pval=pval,
               strand=strand, motif=motif, sequence=tbl.fimo$matched_sequence,
               stringsAsFactors=FALSE)


} # toTrackTable
#----------------------------------------------------------------------------------------------------
get.bach1.binding.sites <- function()
{
    eclip.file <- system.file(package="RnaBindingProtein", "extdata", "ENCFF565FNW.bigBed")
    rbp <- RnaBindingProtein$new("DDX3X", "BACH1", eclip.file, "K562")

    roi <- list(chrom="chr21", start=29274856, end=29392243)
    tbl.bach1.hits <- rbp$getBindingSites(roi) #, c("hg38_genes_3UTRs", "hg38_genes_5UTRs"))

    meme.file <- "ddx3x-ggc.meme"
    tbl.fimo <- fimoBatch(tbl.bach1.hits[, c("chrom", "start", "end")],
                          matchThreshold=1e-3, genomeName="hg38", pwmFile=meme.file,
                          expandResultsWithMotifDb=FALSE)
    tbl.track <- toTrackTable(tbl.fimo)
    track <- DataFrameQuantitativeTrack("enchriched motifs", tbl.track,  autoscale=TRUE, color="darkred")
    displayTrack(igv, track)

} # get.bach1.binding.sites
#----------------------------------------------------------------------------------------------------
viz <- function()
{
    igv <- start.igv("BACH1", "hg38")
    zoomOut(igv)
    track <- DataFrameQuantitativeTrack("DDX3X", tbl.bach1.hits[, c("chrom", "start", "end", "score")],
                                        autoscale=TRUE, color="darkgreen")
    displayTrack(igv, track)


} # viz
#----------------------------------------------------------------------------------------------------

