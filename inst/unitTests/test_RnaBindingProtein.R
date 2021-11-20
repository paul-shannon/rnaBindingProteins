library(RUnit)
library(RnaBindingProtein)
#----------------------------------------------------------------------------------------------------
runTests <- function()
{
    test_ctor()
    test_getAnnotationTypes()
    test_findBindingSites()

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
    x <- rbp$getUTRs()

} # test_getUTRs
#----------------------------------------------------------------------------------------------------
test_findBindingSites <- function(viz=FALSE)
{
    message(sprintf("--- test_findBindingSites"))

    eclip.file <- system.file(package="RnaBindingProtein", "extdata", "ENCFF565FNW.bigBed")
    rbp <- RnaBindingProtein$new("DDX3X", "BACH1", eclip.file, "K562")
    toi <- c("hg38_genes_3UTRs", "hg38_genes_5UTRs")
    checkEquals(sort(grep("UTR", rbp$getAnnotationTypes(), value=TRUE)), toi)
    roi <- list(chrom="chr21", start=29274856, end=29392243)
    tbl.bach1.hits <- rbp$getBindingSites(c("hg38_genes_3UTRs", "hg38_genes_5UTRs"), roi)
    checkEquals(dim(tbl.bach1.hits), c(7, 11))

    if(interactive() & viz){
       igv <- rbp$start.igv()
       #setTrackHeight(igv, "Refseq Genes", 100)
       rbp$showUTRs()
       }


} # test_findBindingSites
#----------------------------------------------------------------------------------------------------
if(!interactive())
    runTests()
