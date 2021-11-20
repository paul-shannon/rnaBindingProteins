#' @title RnaBindingProtein
#' @description assess RNA-binding protein binding sites, their structure and motif matches
#' @name RnaBindingProtein

#' @import annotatr
#' @import rtracklayer
#' @import TxDb.Hsapiens.UCSC.hg38.knownGene
#' @import igvR
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
                   utrAnnotations=NULL,
                   tbl.3utr=NULL,
                   tbl.5utr=NULL,
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
            annotations.file <- system.file(package="RnaBindingProtein", "extdata", "anno.hg38.utrs.RData")
            stopifnot(file.exists(annotations.file))
            private$utrAnnotations <- get(load(annotations.file))
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
          #' extract all 3' and 5' UTRs for thhe targetGene
          #' return data.frame
        getUTRs = function(){
           genic.anno <- private$utrAnnotations[private$utrAnnotations$symbol %in% private$targetGene]
           tbl.anno <- as.data.frame(genic.anno)
           colnames(tbl.anno)[1] <- "chrom"
           tbl.anno$chrom <- as.character(tbl.anno$chrom)
           dups <- which(duplicated(tbl.anno[, c("start", "end", "type")]))
           if(length(dups) > 0)
               tbl.anno <- tbl.anno[-dups,]
           tbl.3utr <- tbl.anno[grep("3UTR", tbl.anno$type),]
           dim(tbl.3utr)
           tbl.5utr <- tbl.anno[grep("5UTR", tbl.anno$type),]
           dim(tbl.5utr)
           #tbl.track <- tbl.5utr[, c("chrom", "start", "end", "id")]
           #track <- DataFrameAnnotationTrack("5'UTR", tbl.track, color="brown", displayMode="Collapse")
           #displayTrack(private$igv, track)
           #tbl.track <- tbl.3utr[, c("chrom", "start", "end", "id")]
           #track <- DataFrameAnnotationTrack("3'UTR", tbl.track, color="black")
           #displayTrack(private$igv, track)
           private$tbl.3utr  <- tbl.3utr
           private$tbl.5utr  <- tbl.5utr
           browser()
           xyz <- 99
           },

        #------------------------------------------------------------
          #' @description
          #' retrieve the rbp bindings sites within the specified region
          #' @param annotation.types  character, one or more annotatr categoies
          #' @param roi list with chrom, start, end fiels
          #' return data.frame
        getBindingSites = function(annotation.types, roi){
          result <- data.frame()
          if(is.list(roi) & all(c("chrom", "start", "end") %in% names(roi))){
             result <- subset(private$tbl.rbpHits, chrom==roi$chrom & start>=roi$start & end <= roi$end)
             }
           result
           } # getBindingSites

       ) # public
    ) # class
#--------------------------------------------------------------------------------
