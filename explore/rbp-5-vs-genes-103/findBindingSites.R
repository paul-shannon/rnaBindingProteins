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
identifyTargetGenes <- function(discordantOnly=FALSE)
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
    if(discordantOnly)
        goi <- subset(tbl.discordances, RNA.protein.discordance=="Y")$geneSymbol
    printf("--- targetGeneCount: %d", length(goi))
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
    gois <- identifyTargetGenes(discordantOnly=TRUE)

    tbls <- list()

    i <- 0
    for(rbpoi in rbps){
      for(cell.line in cell.lines){
         eclip.file <- get.eclip.bigBed.filename(rbpoi, cell.line)
         if(is.na(eclip.file)) next
         checkTrue(file.exists(eclip.file))
         for(goi in gois){
            rbp <- RnaBindingProtein$new(rbpoi, goi, eclip.file, cell.line)
            x <- rbp$getBindingSites.inGenicRegions(intersectionType="within")
            tbl <- x$big
            printf("--- %s -> %s, %s: %d", rbpoi, goi, cell.line, nrow(tbl))
            if(nrow(tbl) == 0) next;
            tbl$cell.line <- cell.line
            tbl$rbp <- rbpoi
            i <- i + 1
            tbls[[i]] <- tbl
            } # for goi
         } # for cell.line
       } # for rbpoi

    tbl <- do.call(rbind, tbls)
    rownames(tbl) <- NULL
    filename <- sprintf("tbl.rbp.discordant-%s.RData", gsub(" ", ".", Sys.time(), fixed=TRUE))
    printf("--- saving %d binding sites to %s", nrow(tbl), filename)
    save(tbl, file=filename)

} # runAll
#----------------------------------------------------------------------------------------------------
memePrep <- function()
{
   tbl <- get(load("~/github/rnaBindingProteins/explore/rbp-5-vs-genes-103/tbl.rbp.all-2022-01-08.16:17:54.RData"))
   dim(tbl)  # 7395 15
   tbl.3utr <- subset(tbl, type.gre=="hg38_genes_3UTRs")
   tbl.5utr <- subset(tbl, type.gre=="hg38_genes_5UTRs")
   tbl.utr  <- rbind(tbl.3utr, tbl.5utr)
   dim(tbl.3utr)  # 2255 15
   dim(tbl.5utr)  # 1086 15

   table(tbl.utr$rbp)
    #  ddx3x   eif3d    fxr2 igf2bp2    rps3
    #    767     361     823     722     668

   table(tbl.3utr$rbp)
    #  ddx3x   eif3d    fxr2 igf2bp2    rps3
    #    381     204     564     670     436
   table(tbl.5utr$rbp)
    #  ddx3x   eif3d    fxr2 igf2bp2    rps3
    #   386     157     259      52     232

   for(rbp.oi in c("ddx3x", "eif3d", "fxr2", "igf2bp2", "rps3")){
      fasta.filename <- sprintf("%s-utr.fa", rbp.oi)
      tbl.sub <- subset(tbl.utr, rbp==rbp.oi)
      sequence.count <- writeFastaFile(tbl.sub, fasta.filename)
      printf("wrote %d %s sequences to %s", sequence.count, rbp.oi, fasta.filename)
      } # for rbp.oi

     # [1] wrote 406 ddx3x sequences to ddx3x-utr.fa
     # [1] wrote 194 eif3d sequences to eif3d-utr.fa
     # [1] wrote 484 fxr2 sequences to fxr2-utr.fa
     # [1] wrote 347 igf2bp2 sequences to igf2bp2-utr.fa
     # [1] wrote 374 rps3 sequences to rps3-utr.fa

} # memePrep
#----------------------------------------------------------------------------------------------------
# ~/meme/bin/meme targetGenes-97-ddx3x-k562.fa -dna -oc meme.out.targetGenes-97-ddx3x-k562
#   -nostatus -time 14400 -mod zoops -nmotifs 3 -minw 6 -maxw 50 -objfun  classic -revcomp -markov_order 0
runMeme <- function(file.basename, data.dir, meme.results.dir)
{
   fasta.filename <- file.path(data.dir, sprintf("%s.fa", file.basename))
   stopifnot(file.exists(fasta.filename))
   cmd.p1 <- "~/meme/bin/meme "
   cmd.p2 <- "-dna -oc"
   cmd.p3 <- "-nostatus -time 14400 -mod zoops"
   cmd.p4 <- "-nmotifs 6 -minw 5 -maxw 50 -objfun  classic -revcomp -markov_order 0"
   cmd <- sprintf("%s %s %s %s %s %s", cmd.p1, fasta.filename, cmd.p2, meme.results.dir, cmd.p3, cmd.p4)
   system(cmd)

   #system(sprintf("open %s.meme/meme.html", rbp.region))

} # runMeme
#----------------------------------------------------------------------------------------------------
multi.meme <- function()
{
    tbl <- get(load("tbl.rbp.discordant-2022-01-09.21:33:35.RData"))

    dim(tbl)
    colnames(tbl)
    table(tbl$rbp)
      #  ddx3x   eif3d    fxr2 igf2bp2    rps3
      #   2089    1138    3459    1404    2072
    table(tbl$type.gre)
      # hg38_genes_3UTRs hg38_genes_5UTRs   hg38_genes_cds
      #      2783             1687             5692
    table(tbl$cellType)
      # HepG2  K562
      #  4651  5511

    output.dir <- "rbp.5-genes.47.discordant-ov.within"

    for(rbp.oi in sort(unique(tbl$rbp))){
        printf("-----------------------")
        for(cellType.oi in sort(unique(tbl$cellType))){
          for(genomicRegion in c("hg38_genes_3UTRs", "hg38_genes_5UTRs", "hg38_genes_cds", "UTRs", "hg38_genes")){
             tbl.oi <- subset(tbl, rbp==rbp.oi & cellType==cellType.oi & grepl(genomicRegion, type.gre))
             printf("%8s %6s %20s: %d", rbp.oi, cellType.oi, genomicRegion, nrow(tbl.oi))
             if(nrow(tbl.oi) > 0){
                tag <- sprintf("%s-%s-%s-%d-within", rbp.oi, cellType.oi, genomicRegion, nrow(tbl.oi))
                fasta.filename <- sprintf("%s.fa", tag)
                full.path <- file.path(output.dir,fasta.filename)
                sequence.count <- writeFastaFile(tbl.oi, full.path)
                stopifnot(file.exists(full.path))
                #browser()
                #xyz <- 99
                meme.results.dir <- file.path(output.dir, sprintf("%s.meme.out", tag))
                runMeme(tag, data.dir=output.dir, meme.results.dir)
                } # nrow > 0
             } # for genomicRegion
            } # for cellType.oi
        } # for rbp.oi

} # multi.meme
#----------------------------------------------------------------------------------------------------


