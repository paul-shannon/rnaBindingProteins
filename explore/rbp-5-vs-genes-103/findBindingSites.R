library(RnaBindingProtein)
library(RUnit)
#----------------------------------------------------------------------------------------------------
tbl.rbp <- data.frame(rbp=c("ddx3x", "fxr2", "eif3d", "igf2bp2", "rps3"),
                      k562=c("ENCFF565FNW-eclip-ddx3x-k562.bigBed",
                             "ENCFF323QTM-eclip-fxr2-k562.bigBed",
                             NA,
                             "ENCFF416JIV-eclip-igf2bp2-k562.bigBed",
                             "ENCFF199OFQ-eclip-rps3-k562.bigBed"),
                      hepg2=c("ENCFF199OFQ-eclip-rps3-k562.bigBed",
                              "ENCFF483HFI-eclip-fxr2-hepg2.bigBed",
                              "ENCFF733BQQ-eclip-eif3d-hepg2.bigBed",
                              NA,
                              "ENCFF848LEJ-eclip-rps3-hepg2.bigBed"),
                      stringsAsFactors=FALSE)
#----------------------------------------------------------------------------------------------------
identifyTargetGenes <- function()
{
    f <- "srm-103-late-discordance-geneSymbolsAdded.tsv"
    data.dir <- "~/github/TrenaProjectErythropoiesis/inst/extdata/geneSets"
    full.path <- file.path(data.dir, f)
    file.exists(full.path)
    tbl.discordances <- read.table(full.path, sep="\t", as.is=TRUE, header=TRUE)
    dim(tbl.discordances)
    head(tbl.discordances)
    goi <- tbl.discordances$geneSymbol
    checkEquals(length(unique(goi)), 103)
    return(goi)

} # identifyTargetGenes
#----------------------------------------------------------------------------------------------------
get.eclip.bigBed.filename <- function(RBP, cell.line)
{
   RBP <- tolower(RBP)
   cell.line <- tolower(cell.line)

   stopifnot(RBP %in% tbl.rbp$rbp)
   filename <- subset(tbl.rbp, rbp==RBP)[, cell.line]
   if(is.na(filename))
       return(NA)
   full.path <- file.path("~/github/rnaBindingProteins/prep/encode-downloads", filename)
   return(full.path)

} # get.eclip.bigBed.filename
#----------------------------------------------------------------------------------------------------
test_get.eclip.bigBed.filename <- function()
{
   rbp <- "DDX3X"
   checkTrue(file.exists(get.eclip.bigBed.filename(rbp, "k562")))
   checkTrue(file.exists(get.eclip.bigBed.filename(rbp, "hepg2")))
   checkTrue(get.eclip.bigBed.filename(rbp, "k562") !=
             get.eclip.bigBed.filename(rbp, "hepg2"))

   rbp <- "FXR2"
   checkTrue(file.exists(get.eclip.bigBed.filename(rbp, "k562")))
   checkTrue(file.exists(get.eclip.bigBed.filename(rbp, "hepg2")))
   checkTrue(get.eclip.bigBed.filename(rbp, "k562") !=
             get.eclip.bigBed.filename(rbp, "hepg2"))

   rbp <- "RPS3"
   checkTrue(file.exists(get.eclip.bigBed.filename(rbp, "k562")))
   checkTrue(file.exists(get.eclip.bigBed.filename(rbp, "hepg2")))
   checkTrue(get.eclip.bigBed.filename(rbp, "k562") !=
             get.eclip.bigBed.filename(rbp, "hepg2"))

   rbp <- "EIF3D"
   checkTrue(is.na(get.eclip.bigBed.filename(rbp, "k562")))
   checkTrue(file.exists(get.eclip.bigBed.filename(rbp, "hepg2")))

   rbp <- "IGF2BP2"
   checkTrue(is.na(get.eclip.bigBed.filename(rbp, "hepg2")))
   checkTrue(file.exists(get.eclip.bigBed.filename(rbp, "k562")))

} # test_get.eclip.bigBed.filename
#----------------------------------------------------------------------------------------------------
runAll <- function()
{
    rbps <- tbl.rbp$rbp
    cell.lines <- c("K562", "HepG2")
    gois <- identifyTargetGenes()

    tbls <- list()

    i <- 0
    for(rbpoi in rbps){
      for(cell.line in cell.lines){
         eclip.file <- get.eclip.bigBed.filename(rbpoi, cell.line)
         if(is.na(eclip.file)) next
         checkTrue(file.exists(eclip.file))
         for(goi in gois){
            rbp <- RnaBindingProtein$new(rbpoi, goi, eclip.file, cell.line)
            x <- rbp$getBindingSites.inGenicRegions()
            tbl <- x$big
            printf("--- %s -> %s, %s: %d", rbpoi, goi, cell.line, nrow(tbl))
            if(nrow(tbl) == 0) next;
            tbl$cell.line <- cell.line
            tbl$rbp <- rbpoi
            coi <- c("chrom","start","end","width","score","start.gre","end.gre","width.gre",
                     "rbp", "cell.line", "geneSymbol", "type.gre", "name")
            i <- i + 1
            tbls[[i]] <- tbl[, coi]
            #browser()
            #xyz <- 99
            } # for goi
         } # for cell.line
       } # for rbpoi

    tbl <- do.call(rbind, tbls)
    rownames(tbl) <- NULL
    save(tbl, file="tbl.rbp.all.RData")

} # runAll
#----------------------------------------------------------------------------------------------------
runMeme <- function()
{
   tbl <- get(load("tbl.rbp.all.RData"))
   tbl <- subset(tbl, width.gre >= 6)
   dim(tbl)

   eclip.file <- system.file(package="RnaBindingProtein", "extdata", "ENCFF565FNW.bigBed")
   rbp <- RnaBindingProtein$new("DDX3X", "BACH1", eclip.file, "K562")

   fasta.filename <- "rbpTargetGenes.10163bindingSites.fa"
   sequence.count <- rbp$writeFastaFile(tbl, fasta.filename)
   checkEquals(nrow(tbl), sequence.count)
   checkTrue(file.exists(fasta.filename))


   cmd <- sprintf("~/meme/bin/meme %s -dna -oc . -nostatus -time 14400 -mod zoops -nmotifs 5 -minw 6 -maxw 50 -objfun  classic -revcomp -markov_order 0", fasta.filename)

   system(cmd)
   system("open meme.html")


} # runMeme
#----------------------------------------------------------------------------------------------------




