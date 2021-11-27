library(RUnit)
library(RnaBindingProtein)
#----------------------------------------------------------------------------------------------------
runTests <- function()
{
    test_ctor()
    test_getAnnotationTypes()
    test_findBindingSites()
    test_getUTRs()

} # runTests
#----------------------------------------------------------------------------------------------------
test_ctor <- function()
{
    message(sprintf("--- test_ctor"))

    eclip.file <- system.file(package="RnaBindingProtein", "extdata", "ENCFF565FNW.bigBed")
    checkTrue(file.exists(eclip.file))
    rbp <- RnaBindingProtein$new("DDX3X", "BACH1", eclip.file, "K562")
    checkTrue(all(c("R6", "RnaBindingProtein") %in% class(rbp)))

    tbl.rbp <- rbp$getBindingTable()
    checkEquals(dim(tbl.rbp), c(29606, 11))
    score.stats <- fivenum(tbl.rbp$score) # 1.000774   1.422952   2.526644   6.016320 400.000000
    checkEquals(score.stats[1], 1.0, tol=1)
    checkEquals(score.stats[5], 400.0, tol=1)


} # test_ctor
#----------------------------------------------------------------------------------------------------
test_getAnnotationTypes <- function()
{
    message(sprintf("--- test_getAnnotationTypes"))

    eclip.file <- system.file(package="RnaBindingProtein", "extdata", "ENCFF565FNW.bigBed")
    rbp <- RnaBindingProtein$new("DDX3X", "BACH1", eclip.file, "K562")

    expected <- c("hg38_basicgenes", "hg38_cpg_inter", "hg38_cpg_islands", "hg38_cpg_shelves",
                  "hg38_cpg_shores", "hg38_cpgs", "hg38_enhancers_fantom", "hg38_genes_1to5kb",
                  "hg38_genes_3UTRs", "hg38_genes_5UTRs", "hg38_genes_cds",
                  "hg38_genes_exonintronboundaries", "hg38_genes_exons", "hg38_genes_firstexons",
                  "hg38_genes_intergenic", "hg38_genes_intronexonboundaries", "hg38_genes_introns",
                  "hg38_genes_promoters", "hg38_lncrna_gencode")
    checkEquals(sort(rbp$getAnnotationTypes()), expected)

} # test_getAnnotationTypes
#----------------------------------------------------------------------------------------------------
test_start.igv <- function()
{
    if(!interactive()) return()

    message(sprintf("--- test_start.igv"))

    eclip.file <- system.file(package="RnaBindingProtein", "extdata", "ENCFF565FNW.bigBed")
    rbp <- RnaBindingProtein$new("DDX3X", "BACH1", eclip.file, "K562")
    igv <- rbp$start.igv()

} # test_start.igv
#----------------------------------------------------------------------------------------------------
test_getUTRs <- function()
{
    message(sprintf("--- test_getUTRs"))

    eclip.file <- system.file(package="RnaBindingProtein", "extdata", "ENCFF565FNW.bigBed")
    rbp <- RnaBindingProtein$new("DDX3X", "BACH1", eclip.file, "K562")
    tbl.gre <- rbp$getGenicRegions()
    checkEquals(dim(tbl.gre), c(26, 10))
    checkEquals(nrow(subset(tbl.gre, type=="hg38_genes_3UTRs")), 2)
    checkEquals(nrow(subset(tbl.gre, type=="hg38_genes_5UTRs")), 12)
    checkEquals(nrow(subset(tbl.gre, type=="hg38_genes_cds")), 12)

} # test_getUTRs
#----------------------------------------------------------------------------------------------------
test_findBindingSites <- function(viz=FALSE)
{
    message(sprintf("--- test_findBindingSites"))

    eclip.file <- system.file(package="RnaBindingProtein", "extdata", "ENCFF565FNW.bigBed")
    rbp <- RnaBindingProtein$new("DDX3X", "BACH1", eclip.file, "K562")

    roi <- list(chrom="chr21", start=29274856, end=29392243)
    tbl.bach1.hits <- rbp$getBindingSites(roi) #, c("hg38_genes_3UTRs", "hg38_genes_5UTRs"))
    checkEquals(dim(tbl.bach1.hits), c(7, 11))

} # test_findBindingSites
#----------------------------------------------------------------------------------------------------
# genic regions of interest: "hg38_genes_3UTRs", "hg38_genes_5UTRs"
test_getBindingSites.inGenicRegions <- function()
{
    message(sprintf("--- test_getBindingSites.inGenicRegions"))

    eclip.file <- system.file(package="RnaBindingProtein", "extdata", "ENCFF565FNW.bigBed")
    rbp <- RnaBindingProtein$new("DDX3X", "BACH1", eclip.file, "K562")

    x <- rbp$getBindingSites.inGenicRegions()
    checkEquals(lapply(x, dim), list(big=c(11, 12), small=c(7, 6)))
    tbl.small <- x$small   # unique DDX3X hits, sometimes in multiple regions
    tbl.big   <- x$big     # unique DDX3X hit/genic region pairs

    regionType.distribution <- as.list(table(tbl.small$regionType))
    checkEquals(regionType.distribution,
                list(hg38_genes_3UTRs=1,
                     hg38_genes_5UTRs=3,
                     hg38_genes_cds=2,
                    "hg38_genes_cds;hg38_genes_5UTRs"=1))

    viz <- function(){   # check these results visually
       igv <- start.igv("BACH1")
       zoomOut(igv)
         # for testing, get the genic regions associated with targetGene, independently, build roi
       tbl.gre <- rbp$getGenicRegions()
       shoulder <- 1000
       roi <- list(chrom=tbl.gre$chrom[1], start=min(tbl.gre$start)-shoulder, end=max(tbl.gre$end)+shoulder)
            #-------------------------
            # the DDX3X hits first
            #-------------------------
       tbl.hits <- rbp$getBindingSites(roi)
       track <- DataFrameQuantitativeTrack("DDX3X", tbl.hits[, c("chrom", "start", "end", "score")],
                                           autoscale=TRUE, color="brown")
       displayTrack(igv, track)

            #-------------------------
            # the 5'UTR regions
            #-------------------------
       track <- DataFrameAnnotationTrack("5'UTRs", subset(tbl.gre, type=="hg38_genes_5UTRs"), color="red")
       displayTrack(igv, track)

            #-------------------------
            # the 3'UTR regions
            #-------------------------
       track <- DataFrameAnnotationTrack("3'UTRs", subset(tbl.gre, type=="hg38_genes_3UTRs"), color="blue")
       displayTrack(igv, track)

            #-------------------------
            # CDS regions
            #-------------------------
       tbl.track <- subset(tbl.gre, type=="hg38_genes_cds")
       track <- DataFrameAnnotationTrack("CDS", tbl.track, color="darkGreen")
       displayTrack(igv, track)

            #--------------------------------------------------
            # DDX3X in genic regions, simplified in tbl.small
            #--------------------------------------------------
       region.types <- unique(sort(tbl.small$regionType))
       for(type in region.types){
         track <- DataFrameAnnotationTrack(type, subset(tbl.small, regionType==type))
         displayTrack(igv, track)
         }
       } # viz


} # test_getBindingSites.inUTRs
#----------------------------------------------------------------------------------------------------
if(!interactive())
    runTests()
