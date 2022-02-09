library(RUnit)
library(RnaBindingProtein)
#----------------------------------------------------------------------------------------------------
runTests <- function()
{
    test_ctor()
    test_getEncodeCatalog()
    test_getAnnotationTypes()
    test_getGeneRegions()
    test_getBindingSites()
    test_getGenicRegions()
    test_getAllGenicAnnotations()
    test_getBindingSites.inGenicRegions()
    test_writeFastaFile()
    test_runMeme()

} # runTests
#----------------------------------------------------------------------------------------------------
test_ctor <- function()
{
    message(sprintf("--- test_ctor"))

    rbp.tool <- RnaBindingProtein$new("DDX3X", "BACH1", "K562")
    checkTrue(all(c("R6", "RnaBindingProtein") %in% class(rbp.tool)))
    rbp.tool <- RnaBindingProtein$new("DDX3X", targetGene=NA, "K562")
    checkTrue(all(c("R6", "RnaBindingProtein") %in% class(rbp.tool)))

    tbl.rbp <- rbp.tool$getBindingTable()
    checkEquals(dim(tbl.rbp), c(29606, 11))
    score.stats <- fivenum(tbl.rbp$score) # 1.000774   1.422952   2.526644   6.016320 400.000000
    checkEquals(score.stats[1], 1.0, tol=1)
    checkEquals(score.stats[5], 400.0, tol=1)

    tbl.encodeMappings <- rbp.tool$getEncodeCatalog()

    tbl.rbp <- RnaBindingProtein$new("UPF1", NA, "HepG2")

    set.seed(17)
    i <- sample(seq_len(nrow(tbl.encodeMappings)), size=1)
    rbp.tool <- RnaBindingProtein$new(rbp=tbl.encodeMappings$target[i],
                                      targetGene=NA,
                                      cellType=tbl.encodeMappings$cellType[i])
    tbl.bindings <- rbp.tool$getBindingTable()
    checkTrue(nrow(tbl.bindings) > 1000)
    checkEquals(ncol(tbl.bindings), 11)

} # test_ctor
#----------------------------------------------------------------------------------------------------
test_getEncodeCatalog <- function()
{
    message(sprintf("--- test_getEncodeMappings"))

    rbp.tool <- RnaBindingProtein$new("DDX3X", targetGene=NA, "K562")
    tbl.cat <- rbp.tool$getEncodeCatalog()
    checkEquals(dim(tbl.cat), c(225, 3))

    colnames(tbl.cat)[3] <- "rbp"
    tbl.cat <- tbl.cat[order(tbl.cat$rbp),]
    rownames(tbl.cat) <- NULL
    tbl.cat[, c("rbp", "cellType", "id")]

} # test_getEncodeCatalog
#----------------------------------------------------------------------------------------------------
test_getAnnotationTypes <- function()
{
    message(sprintf("--- test_getAnnotationTypes"))

    rbp.tool <- RnaBindingProtein$new("DDX3X", "BACH1", "K562")

    expected <- c("hg38_basicgenes", "hg38_cpg_inter", "hg38_cpg_islands", "hg38_cpg_shelves",
                  "hg38_cpg_shores", "hg38_cpgs", "hg38_enhancers_fantom", "hg38_genes_1to5kb",
                  "hg38_genes_3UTRs", "hg38_genes_5UTRs", "hg38_genes_cds",
                  "hg38_genes_exonintronboundaries", "hg38_genes_exons", "hg38_genes_firstexons",
                  "hg38_genes_intergenic", "hg38_genes_intronexonboundaries", "hg38_genes_introns",
                  "hg38_genes_promoters", "hg38_lncrna_gencode")
    checkEquals(sort(rbp.tool$getAnnotationTypes()), expected)

} # test_getAnnotationTypes
#----------------------------------------------------------------------------------------------------
test_start.igv <- function()
{
    if(!interactive()) return()

    message(sprintf("--- test_start.igv"))

    eclip.file <- system.file(package="RnaBindingProtein", "extdata", "ENCFF565FNW.bigBed")
    rbp.tool <- RnaBindingProtein$new("DDX3X", "BACH1", eclip.file, "K562")
    igv <- rbp.tool$start.igv()

} # test_start.igv
#----------------------------------------------------------------------------------------------------
test_getGeneRegion <- function()
{
    goi <- c("CHD1", "CREBBP", "CTCF", "E2F4", "EGR1", "ELF1", "GATAD2A", "GTF2F1",
             "HDAC2", "HMGB3", "JUND", "KLF1", "KLF3", "LDB1", "NELFE", "NFE2",
             "NR3C1", "NRF1", "POLR2A", "POU2F1", "RAD21", "RCOR1", "RFX5", "SMC3",
             "TAL1", "TCF3", "WDR5", "ZBTB7A")

    for(gene in goi){
        rbp.tool <- RnaBindingProtein$new("DDX3X", gene, "K562")
        roi <- rbp.tool$getGeneRegion()
        checkEquals(names(roi), c("chrom", "start", "end"))
        }

} # testGetGeneRegion
#----------------------------------------------------------------------------------------------------
test_getGenicRegions <- function()
{
    message(sprintf("--- test_getGenicRegions"))

    rbp.tool <- RnaBindingProtein$new("DDX3X", "BACH1", "K562")
    tbl.gre <- rbp.tool$getGenicRegions()
    checkEquals(dim(tbl.gre), c(219, 10))
    checkEquals(nrow(subset(tbl.gre, type=="hg38_genes_3UTRs")), 3)
    checkEquals(nrow(subset(tbl.gre, type=="hg38_genes_5UTRs")), 22)
    checkEquals(nrow(subset(tbl.gre, type=="hg38_genes_cds")), 23)

} # test_getGenicRegions
#----------------------------------------------------------------------------------------------------
test_getAllGenicAnnotations <- function()
{
    message(sprintf("--- test_getAllGenicAnnotations"))

    #eclip.file <- system.file(package="RnaBindingProtein", "extdata", "ENCFF565FNW.bigBed")
    rbp.tool <- RnaBindingProtein$new("DDX3X", NA, "K562")
    tbl.anno <- rbp.tool$getAllGenicAnnotations()
    checkTrue(nrow(tbl.anno) > 500000)
    checkEquals(ncol(tbl.anno), 10)
    genes <- sort(unique(tbl.anno$symbol))
    checkTrue(length(genes) > 19000)

} # test_getAllGenicAnnotations
#----------------------------------------------------------------------------------------------------
test_getBindingSites <- function(viz=FALSE)
{
    message(sprintf("--- test_getBindingSites"))

    rbp.tool <- RnaBindingProtein$new("DDX3X", "KLF1", "K562")
    roi <- rbp.tool$getGeneRegion()
    tbl.hits <- rbp.tool$getBindingSites(roi)
    checkEquals(dim(tbl.hits), c(2, 11))


} # test_getBindingSites
#----------------------------------------------------------------------------------------------------
test_getBindingSites.inGenicRegions <- function()
{
    message(sprintf("--- test_getBindingSites.inGenicRegions"))

    rbp.tool <- RnaBindingProtein$new("DDX3X", "BACH1", "K562")

    x <- rbp.tool$getBindingSites.inGenicRegions()
    checkEquals(lapply(x, dim), list(big=c(82, 14), small=c(7, 9)))
    tbl.small <- x$small   # unique DDX3X hits, sometimes in multiple regions
    tbl.big   <- x$big     # unique DDX3X hit/genic region pairs

    genic.regions.of.interest <-  c("hg38_genes_3UTRs", "hg38_genes_5UTRs", "hg38_genes_cds")

    tbl.big.sub <- subset(tbl.big, type.gre %in% genic.regions.of.interest)
    checkEquals(as.list(table(tbl.big.sub$type.gre)),
                list(hg38_genes_3UTRs=2,
                     hg38_genes_5UTRs=3,
                     hg38_genes_cds=7))

    viz <- function(){   # check these results visually
       igv <- start.igv("BACH1")
       zoomOut(igv)
         # for testing, get the genic regions associated with targetGene, independently, build roi
       tbl.gre <- rbp.tool$getGenicRegions()
       shoulder <- 1000
       roi <- list(chrom=tbl.gre$chrom[1], start=min(tbl.gre$start)-shoulder, end=max(tbl.gre$end)+shoulder)
            #-------------------------
            # the DDX3X hits first
            #-------------------------
       tbl.hits <- rbp.tool$getBindingSites(roi)
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

} # test_getBindingSites.inGenicRegions
#----------------------------------------------------------------------------------------------------
test_getBindingSites.inUTRs.assorted.other.genes <- function()
{
    nonStandardNames <- c("BCL11A_XL_L")
    # eclip.file <- system.file(package="RnaBindingProtein", "extdata", "ENCFF565FNW.bigBed")
    rbp.tool <- RnaBindingProtein$new("DDX3X", nonStandardNames[1],  "K562")
    x <- rbp.tool$getBindingSites.inGenicRegions()

} # test_getBindingSites.inUTRs.assorted.other.genes
#----------------------------------------------------------------------------------------------------
test_writeFastaFile <- function()
{
   message(sprintf("--- test_writeFastaFile"))

   #eclip.file <- system.file(package="RnaBindingProtein", "extdata", "ENCFF565FNW.bigBed")
   rbp.tool <- RnaBindingProtein$new("DDX3X", "BACH1", "K562")
   x <- rbp.tool$getBindingSites.inGenicRegions(intersectionType="within")
   lapply(x, dim)
   tbl <- x$small
   checkTrue(nrow(tbl) > 5)
   checkEquals(4, writeFastaFile(tbl, "tmp.fa", minimum.sequence.length=40))
   checkEquals(7, writeFastaFile(tbl, "tmp.fa", minimum.sequence.length=5))
   checkTrue(file.exists("tmp.fa"))

   stringSet <- readDNAStringSet("tmp.fa")
   checkEquals(width(stringSet), c(35, 42, 34, 67, 143, 36, 41))

} # test_writeFastaFile
#----------------------------------------------------------------------------------------------------
test_runMeme <- function()
{
   message(sprintf("--- test_runMeme"))

   #file.of.previous.results <- system.file(package="RnaBindingProtein",
   #                                        "extdata", "tbl.sites.28goi.RData")
   #if(!file.exists(file.of.previous.results)){

     # extract and combine binding sites in marjorie's rna-seq cluster 3
   goi <- c("CHD1", "CREBBP", "CTCF", "E2F4", "EGR1", "ELF1", "GATAD2A", "GTF2F1",
            "HDAC2", "HMGB3", "JUND", "KLF1", "KLF3", "LDB1", "NELFE", "NFE2",
            "NR3C1", "NRF1", "POLR2A", "POU2F1", "RAD21", "RCOR1", "RFX5", "SMC3",
            "TAL1", "TCF3", "WDR5", "ZBTB7A")

   tbls.sites <- list()
   for(targetGene in goi){
     tryCatch({
        printf("--- %d) %s", grep(targetGene, goi), targetGene)
        rbp.tool <- RnaBindingProtein$new("DDX3X", targetGene, "K562")
        x <- rbp.tool$getBindingSites.inGenicRegions(intersectionType="any")
        printf(" %d sites in genic regions", nrow(x$small))
        tbls.sites[[targetGene]] <- x$small   # duplicates are removed, genic types collapsed
        roi <- rbp.tool$getGeneRegion()
        tbl <- rbp.tool$getBindingSites(roi)
        printf(" %d sites in whole gene", nrow(tbl))
     }, error = function(e){
         printf("failed with %s", targetGene)
         print(e)
         xyz <- 99
         })
      } # for targetGene
   tbl.sites <- do.call(rbind, tbls.sites)
   save(tbl.sites, file="tbl.sites.28goi.RData")
   load("tbl.sites.28goi.RData")
   load("tbl.sites.10goi.RData")
   threshold <- fivenum(tbl.sites$score)[4]
   tbl.sites.filtered <- subset(tbl.sites, score > threshold)
   rbp.tool <- RnaBindingProtein$new("DDX3X", NA, "K562")
   html.path <- rbp.tool$runMeme("goi28_test_runMeme", tbl.sites.filtered)
   checkTrue(file.exists(html.path))
   browseURL(html.path)

} # test_runMeme
#----------------------------------------------------------------------------------------------------
tmp <- function()
{
  library(org.Hs.eg.db); db <- org.Hs.eg.db
  library(TxDb.Hsapiens.UCSC.hg38.knownGene); txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
  goi <- "BACH1"
  goi <- "GRIK1"
  goi <- "KLF1"
  geneID <- select(org.Hs.eg.db, keys=goi, keytype="SYMBOL", columns=c("ENTREZID"))$ENTREZID
  tbl <- select(txdb, keys=geneID, keytype="GENEID", columns=c("TXCHROM", "TXSTART", "TXEND"))
  dim(tbl)

  track <- DataFrameAnnotationTrack(goi, tbl[, 2:4], color="random")
  displayTrack(igv, track)
  showGenomicRegion(igv, goi)


} # tmp
#----------------------------------------------------------------------------------------------------
if(!interactive())
    runTests()
