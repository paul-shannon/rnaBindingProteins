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
                   tbl.utrs=NULL,
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
            annotations.file <- system.file(package="RnaBindingProtein", "extdata", "UTRS-ucsc.RData")
            stopifnot(file.exists(annotations.file))
            private$tbl.utrs <- get(load(annotations.file))
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
          #' extract all 3' and 5' UTRs for the targetGene
          #' return data.frame
        getUTRs = function(){
           subset(private$tbl.utrs, symbol==private$targetGene)
           },

        #------------------------------------------------------------
          #' @description
          #' retrieve the rbp bindings sites within the specified region
          #' @param roi list with chrom, start, end fields
          #' @param annotation.types  character, one or more annotatr categoies, default is
          #' c("hg38_genes_3UTRs", "hg38_genes_5UTRs")
          #' return data.frame
        getBindingSites = function(roi, annotation.types=c("hg38_genes_3UTRs", "hg38_genes_5UTRs")){
          result <- data.frame()
          if(is.list(roi) & all(c("chrom", "start", "end") %in% names(roi))){
             result <- subset(private$tbl.rbpHits, chrom==roi$chrom & start>=roi$start & end <= roi$end)
             }
           result
           }, # getBindingSites

        #------------------------------------------------------------
          #' @description
          #' retrieve the rbp bindings which intersect with UTRs for the targetGene
          #' @param intersectionType character, one of "any" or "within"
          #' return data.frame
        getBindingSites.inUTRs = function(intersectionType="any"){
          stopifnot(intersectionType %in% c("any", "within"))
          tbl.utrs <- self$getUTRs()
          roi <- list(chrom=tbl.utrs$chrom[1], start=min(tbl.utrs$start), end=max(tbl.utrs$end))
          tbl.bindingSites <- self$getBindingSites(roi)
          tbl.out <- data.frame()
          tbl.ov <- as.data.frame(findOverlaps(GRanges(tbl.bindingSites), GRanges(tbl.utrs)))
          if(nrow(tbl.ov) > 0){
             tbl.out <- cbind(tbl.bindingSites[tbl.ov[,1],], tbl.utrs[tbl.ov[,2],])
             colnames(tbl.out)[c(12:15, 18, 20, 21)] <- c("chrom.utr", "start.utr", "end.utr", "width.utr", "enst", "geneSymbol", "type.utr")
             coi <- c("chrom", "start", "end", "width", "name", "score", "start.utr", "end.utr", "width.utr", "geneSymbol", "enst", "type.utr")
             tbl.out <- tbl.out[, coi]
             colnames(tbl.out) <- NULL
             }
          tbl.out
          } # getBindingSites

       ) # public
    ) # class
#--------------------------------------------------------------------------------
