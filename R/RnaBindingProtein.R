#' @title RnaBindingProtein
#' @description assess RNA-binding protein binding sites, their structure and motif matches
#' @name RnaBindingProtein

#' @import annotatr
#' @import rtracklayer
#' @import TxDb.Hsapiens.UCSC.hg38.knownGene
#' @import BSgenome.Hsapiens.UCSC.hg38
#'
#' @field rbp the gene symbol name of the protein
#'
#' @examples
#'   rbp <- RnaBindingProtein("DDX3X", "path/to/ENCFF565FNW.bigBed", "K562", "path/to/ddx3x.meme")
#' @export
#----------------------------------------------------------------------------------------------------
RnaBindingProtein = R6Class("RnaBindingProtein",

    #--------------------------------------------------------------------------------
    private = list(rbp=NULL,
                   tbl.encodeMappings=NULL,
                   genome=NULL,
                   cellType=NULL,
                   bigBedFile=NULL,       # likely from encode, eCLIP results
                   tbl.rbpHits=NULL,
                   biomart=NULL,
                   motifs.meme.file=NULL,
                   targetGene=NULL,
                      # originally only 3'UTR and 5'UTR, then added CDS when some DDX3X binding found there
                   tbl.genicRegions=NULL,
                   annotations.fileName=NULL,
                   igv=NULL,
                   importBigBedFile=function(filename){
                       gr <- import(filename)
                       stopifnot("field9" %in% colnames(mcols(gr)))
                       gr$score <- as.numeric(gr$field9)
                       gr <- gr[gr$score > 1]
                       tbl<- as.data.frame(gr)
                       colnames(tbl)[1] <- "chrom"
                       tbl$chrom <- as.character(tbl$chrom)
                       invisible(tbl)
                       } # importBigBedFile
                   ),


    #--------------------------------------------------------------------------------
    public = list(
         #' @description
         #' Creates a new instance of this [R6][R6::R6Class] class.
         #'
         #' @param rpb character, gene symbol identifier for the rbp
         #' @param targetGene character, gene symbol presumably having associated rbp binding sites
         #' @param cellType character, currently "K562" or "HepG2"
         #' @param motifs.meme.file character, full path to an optional motifs file,
         #'        calculated separately, describing motifs for the rbp
        initialize = function(rbp, targetGene, cellType, motifs.meme.file=NA){
            private$cellType <- cellType
            private$targetGene <- targetGene
            private$rbp <- rbp
            private$tbl.encodeMappings <-
                get(load(system.file(package="RnaBindingProtein", "extdata", "vanNostrand-2020",
                                     "encode.ids.targetGenes.cellTypes.RData")))
            tbl.sub <- subset(private$tbl.encodeMappings,
                              cellType==private$cellType & target==private$rbp)
            if(nrow(tbl.sub) == 0){
                message(sprintf("no ENCODE bigBed file for rbp %s in cellType %s", rbp, cellType))
                stop();
                }
            private$bigBedFile <- system.file(package="RnaBindingProtein", "extdata", "vanNostrand-2020",
                                              sprintf("%s.bigBed", tbl.sub$id))
            stopifnot(file.exists(private$bigBedFile))
            private$genome <- "hg38"            # fixed for now
            private$tbl.rbpHits <- private$importBigBedFile(private$bigBedFile)
            private$motifs.meme.file <- motifs.meme.file
            private$annotations.fileName <- system.file(package="RnaBindingProtein", "extdata",
                                                        "tbl.anno-hg38-11-categories.RData")
                                        #"tbl.anno.hg38.utrs.CDS.RData")
            stopifnot(file.exists(private$annotations.fileName))
            }, # ctor

        #------------------------------------------------------------
          #' @description
          #' @return the table of RBPs/cellTypes for which we have bigBed files from ENCODE
        getEncodeCatalog = function(){
            invisible(private$tbl.encodeMappings)
            },


        #------------------------------------------------------------
          #' @description
          #' @return the slightly processed and filtered data.frame version of the ENCODE bigbed
        getBindingTable = function(){
           invisible(private$tbl.rbpHits)
           },

        #------------------------------------------------------------
          #' @description
          #' the bioc annotatr class offers a large set of genic annotations, ultimately
          #' obtained from UCSC via TxDb.Hsapiens.UCSC.hg38.knownGene
          #' @return genic annotation types for hg38
        getAnnotationTypes = function(){
           grep("hg38", builtin_annotations(), v=TRUE)
           },

        #------------------------------------------------------------
          #' @description
          #' start, and retain a reference to, an igvR for the current genome and targetGene
          #' @return igv reference
        start.igv = function(){
           igv <- igvR()
           setGenome(igv, private$genome)
           showGenomicRegion(igv, private$targetGene)
           setBrowserWindowTitle(igv, private$targetGene)
           # setTrackHeight(igv, "Refseq Genes", 120)
           zoomOut(igv)
           private$igv <- igv
           igv
           },

        #------------------------------------------------------------
          #' @description
          #' extract all genic regions, mulitiple types, no longer just
          #' 3' and 5' UTRs and CDS) for the targetGene
          #' @return data.frame
        getGenicRegions = function(){
           if(is.null(private$tbl.genicRegions))
               private$tbl.genicRegions <- get(load(private$annotations.fileName))

           if(!private$targetGene %in% private$tbl.genicRegions$symbol){
              message(sprintf("targetGene %s not annotated wrt region", private$targetGene))
              return(data.frame())
              }
           subset(private$tbl.genicRegions, symbol==private$targetGene)
           },

        #------------------------------------------------------------
          #' @description
          #' access to the full whole genome annotations table
          #' @return data.frame
        getAllGenicAnnotations = function(){
           if(is.null(private$tbl.genicRegions))
                       private$tbl.genicRegions <- get(load(private$annotations.fileName))

           invisible(private$tbl.genicRegions)
           },

        #------------------------------------------------------------
          #' @description
          #' retrieve the rbp bindings sites within the specified region
          #' @param roi list with chrom, start, end fields
          #' @return data.frame
        getBindingSites = function(roi){
          result <- data.frame()
          if(is.list(roi) & all(c("chrom", "start", "end") %in% names(roi))){
             result <- subset(private$tbl.rbpHits, chrom==roi$chrom & start>=roi$start & end <= roi$end)
             }
           result
           }, # getBindingSites

        #------------------------------------------------------------
          #' @description
          #' retrieve the chrom, start and end of the targetGene, if !na
          #' @return list with chrom, start end, geneSymbol fields
        getGeneRegion = function(){
          if(is.na(private$targetGene)) return(list())
          suppressMessages({
             geneID <- select(org.Hs.eg.db, keys=private$targetGene,
                              keytype="SYMBOL", columns=c("ENTREZID"))$ENTREZID
             tbl <- select(TxDb.Hsapiens.UCSC.hg38.knownGene,
                           keys=geneID, keytype="GENEID", columns=c("TXCHROM", "TXSTART", "TXEND"))
             })
          roi <- list(chrom=tbl$TXCHROM[1], start=min(tbl$TXSTART), end=max(tbl$TXEND))
          roi
          }, # getGeneRegion

        #------------------------------------------------------------
          #' @description
          #' retrieve the rbp bindings which intersect with annotated genic regions of interest,
          #' currently UTRs and CDSfor the targetGene
          #' @param intersectionType character, one of "any" or "within"
          #' @return data.frame
        getBindingSites.inGenicRegions = function(intersectionType="within"){
          stopifnot(intersectionType %in% c("any", "within"))
          tbl.genicRegions <- self$getGenicRegions()
          roi <- list(chrom=tbl.genicRegions$chrom[1], start=min(tbl.genicRegions$start),
                      end=max(tbl.genicRegions$end))
          tbl.bindingSites <- self$getBindingSites(roi)
          tbl.big <- data.frame()
          tbl.small <- data.frame()
          tbl.ov <- as.data.frame(findOverlaps(GRanges(tbl.bindingSites), GRanges(tbl.genicRegions),
                                               type=intersectionType))
          if(nrow(tbl.ov) > 0){
             colnames(tbl.ov) <- c("bs", "gr")
             tbl.big <- cbind(tbl.bindingSites[tbl.ov[,1],], tbl.genicRegions[tbl.ov[,2],])
             colnames(tbl.big)[c(12:15, 18, 20, 21)] <- c("chrom.gre", "start.gre", "end.gre", "width.gre", "enst", "geneSymbol", "type.gre")
             coi <- c("chrom", "start", "end", "width", "name", "score", "start.gre", "end.gre", "width.gre", "geneSymbol", "enst", "type.gre")
             tbl.big <- tbl.big[, coi]
             rownames(tbl.big) <- NULL
             tbl.big$rbp <- private$rbp
             tbl.big$cellType <- private$cellType
                # compress tbl.big to tbl.small, just one row for each binding site
             hits <- unique(tbl.ov$bs)
             tbls <- list()
             i <- 1
             for(hit in hits){
                grs <- subset(tbl.ov, bs==hit)$gr
                gr.types <- paste(unique(tbl.genicRegions[grs, "type"]), collapse=";")
                tbl.bs <- tbl.bindingSites[hit,]
                tbl.bs$regionType <- gr.types
                tbls[[i]] <- tbl.bs
                i <- i + 1
                } # for hit

             coi <- c("chrom", "start", "end", "width", "score", "regionType")
             tbl.small <- do.call(rbind, tbls) [, coi]
             tbl.small$targetGene <- private$targetGene
             tbl.small$cellType <- private$cellType
             tbl.small$rbp <- private$rbp
             rownames(tbl.small) <- NULL
             }
          list(big=tbl.big, small=tbl.small)
          }, # getBindingSites

        #------------------------------------------------------------
          #' @description
          #' create fasta file from filtered rbp binding sites, run meme
          #' @param tagName character something memorable and explanatory
          #' @param top.quartiles.to.include numeric, 1 or more
          #' @return character string, the path to the meme html output
        runMeme = function(tag.name, bl.sitesFiltered) {
           base.name <- sprintf("%s.rbp.%s.cellType.%s.sites.count.%d",
                                tag.name, private$rbp, private$cellType, nrow(tbl.sites))

           printf("--- base.name: %s", base.name)
           fasta.filename <- sprintf("%s.fa", base.name)
           writeFastaFile(tbl.sites, fasta.filename, minimum.sequence.length=6, private$rbp,
                          private$targetGene, private$cellType)
           args1 <- sprintf("-dna -oc %s", base.name)
           args2 <- "-nostatus -time 14400 -mod zoops -nmotifs 10 -minw 6 -maxw 100 -objfun  classic -revcomp -markov_order 0"
           cmd <- sprintf("~/meme/bin/meme %s %s %s", fasta.filename, args1, args2)
           system(cmd)
           return(sprintf("%s/meme.html", base.name))
           } # runMeme

        #------------------------------------------------------------
       ) # public
    ) # class
#--------------------------------------------------------------------------------
#' @name writeFastaFile
#' @description
#' retrieve and write to file the DNA sequence for every genomic region in tbl
#' @param tbl data.frame produced by this class
#' @param fasta.filename character, typically with an ".fa" file extension
#' @param minimum.sequence.length integer, default 5
#' @param rbp character the name of the rna binding protein
#' @param targetGene character - if any
#' @param cellType character e.g. K562, HepG2
#' @return integer count of sequences
#' @export
writeFastaFile = function(tbl, fasta.filename, minimum.sequence.length=5, rbp, targetGene, cellType)
{
   tbl <- subset(tbl, width >= minimum.sequence.length)
   stopifnot(nrow(tbl) > 0)

   location.duplicates <- which(duplicated(tbl[, c("chrom", "start", "end")]))
   if(length(location.duplicates) > 0)
       tbl <- tbl[-location.duplicates,]
   dna.string.set <- with(tbl, getSeq(BSgenome.Hsapiens.UCSC.hg38, chrom, start, end))
   names <- with(tbl, sprintf("%s-%s-%s-%s:%d-%d",
                             rbp, targetGene, cellType, chrom, start, end))
   names(dna.string.set) <- names
   writeXStringSet(dna.string.set, fasta.filename)
   length(dna.string.set)

} # writeFastaFile
#--------------------------------------------------------------------------------
