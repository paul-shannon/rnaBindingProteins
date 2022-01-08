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
                   genome=NULL,
                   cellType=NULL,
                   bigBedFile=NULL,       # likely from encode, eCLIP results
                   tbl.rbpHits=NULL,
                   motifs.meme.file=NULL,
                   targetGene=NULL,
                      # originally only 3'UTR and 5'UTR, then added CDS when some DDX3X binding found there
                   tbl.genicRegions=NULL,
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
         #' @param bigBedFile character, full path to the ENCODE file with binding sites
         #' @param cellType character, typically "K562"
         #' @param motifs.meme.file character, full path to an optional motifs file,
         #'        calculated separately, describing motifs for the rbp
        initialize = function(rbp, targetGene, bigBedFile, cellType, motifs.meme.file=NA){
            private$rbp <- rbp
            private$targetGene <- targetGene
            private$genome <- "hg38"            # fixed for now
            stopifnot(file.exists(bigBedFile))
            private$bigBedFile <- bigBedFile
            private$tbl.rbpHits <- private$importBigBedFile(bigBedFile)
            private$cellType <- cellType
            private$motifs.meme.file <- motifs.meme.file
            annotations.file <- system.file(package="RnaBindingProtein", "extdata",
                                            "tbl.anno.hg38.utrs.CDS.RData")
            stopifnot(file.exists(annotations.file))
            private$tbl.genicRegions <- get(load(annotations.file))
            if(!targetGene %in% private$tbl.genicRegions$symbol){
               msg <- sprintf("--- RnaBindingProtein ctor, unrecognized gene symbol '%s'", targetGene)
               stop(msg)
               }
            },

        #------------------------------------------------------------
          #' @description
          #' returns the slightly processed and filtered data.frame version of the ENCODE bigbed
        getBindingTable = function(){
           invisible(private$tbl.rbpHits)
           },

        #------------------------------------------------------------
          #' @description
          #' the bioc annotatr class offers a large set of genic annotations, ultimately
          #' obtained from UCSC via TxDb.Hsapiens.UCSC.hg38.knownGene
          #' return genic annotation types for hg38
        getAnnotationTypes = function(){
           grep("hg38", builtin_annotations(), v=TRUE)
           },

        #------------------------------------------------------------
          #' @description
          #' start, and retain a reference to, an igvR for the current genome and targetGene
          #' return igv reference
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
          #' extract all 3' and 5' UTRs and CDS for the targetGene
          #' return data.frame
        getGenicRegions = function(){
           subset(private$tbl.genicRegions, symbol==private$targetGene)
           },

        #------------------------------------------------------------
          #' @description
          #' access to the full whole genome annotations table
          #' return data.frame
        getAllGenicAnnotations = function(){
           invisible(private$tbl.genicRegions)
           },

        #------------------------------------------------------------
          #' @description
          #' retrieve the rbp bindings sites within the specified region
          #' @param roi list with chrom, start, end fields
          #' return data.frame
        getBindingSites = function(roi){
          result <- data.frame()
          if(is.list(roi) & all(c("chrom", "start", "end") %in% names(roi))){
             result <- subset(private$tbl.rbpHits, chrom==roi$chrom & start>=roi$start & end <= roi$end)
             }
           result
           }, # getBindingSites

        #------------------------------------------------------------
          #' @description
          #' retrieve the rbp bindings which intersect with annotated genic regions of interest,
          #' currently UTRs and CDSfor the targetGene
          #' @param intersectionType character, one of "any" or "within"
          #' return data.frame
        getBindingSites.inGenicRegions = function(intersectionType="any"){
          stopifnot(intersectionType %in% c("any", "within"))
          tbl.genicRegions <- self$getGenicRegions()
          roi <- list(chrom=tbl.genicRegions$chrom[1], start=min(tbl.genicRegions$start),
                      end=max(tbl.genicRegions$end))
          tbl.bindingSites <- self$getBindingSites(roi)
          tbl.big <- data.frame()
          tbl.small <- data.frame()
          tbl.ov <- as.data.frame(findOverlaps(GRanges(tbl.bindingSites), GRanges(tbl.genicRegions)))
          colnames(tbl.ov) <- c("bs", "gr")
          if(nrow(tbl.ov) > 0){
             tbl.big <- cbind(tbl.bindingSites[tbl.ov[,1],], tbl.genicRegions[tbl.ov[,2],])
             colnames(tbl.big)[c(12:15, 18, 20, 21)] <- c("chrom.gre", "start.gre", "end.gre", "width.gre", "enst", "geneSymbol", "type.gre")
             coi <- c("chrom", "start", "end", "width", "name", "score", "start.gre", "end.gre", "width.gre", "geneSymbol", "enst", "type.gre")
             tbl.big <- tbl.big[, coi]
             rownames(tbl.big) <- NULL
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
             tbl.small <- do.call(rbind, tbls)[, coi]
             tbl.small$targetGene <- private$targetGene
             tbl.small$rbp <- private$rbp
             rownames(tbl.small) <- NULL
             }
          list(big=tbl.big, small=tbl.small)
          }, # getBindingSites
        #------------------------------------------------------------
          #' @description
          #' retrieve and write to file the DNA sequence for every genomic region in tbl
          #' @param tbl data.frame produced by this class
          #' @param fasta.filename character, typically with an ".fa" file extension
          #' return integer count of sequences
        writeFastaFile = function(tbl, fasta.filename){
           stopifnot(all(c("chrom", "start", "end", "rbp", "geneSymbol", "cell.line") %in% colnames(tbl)))
           dna.string.set <- with(tbl, getSeq(BSgenome.Hsapiens.UCSC.hg38, chrom, start, end))
           names <- with(tbl, sprintf("%s-%s-%s-%s:%d-%d",
                                      rbp, geneSymbol, cell.line, chrom, start, end))
           names(dna.string.set) <- names
           writeXStringSet(dna.string.set, fasta.filename)
           length(dna.string.set)
           }
       ) # public
    ) # class
#--------------------------------------------------------------------------------
