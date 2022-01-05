library(RnaBindingProtein)
library(BSgenome.Hsapiens.UCSC.hg38)
library(annotatr)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
source("~/github/fimoService/batchMode/fimoBatchTools.R")

#----------------------------------------------------------------------------------------------------
f <- "srm-103-late-discordance-geneSymbolsAdded.tsv"

data.dir <- "~/github/TrenaProjectErythropoiesis/explore/rbp/ddx3x"
full.path <- file.path(data.dir, f)
file.exists(full.path)
tbl.discordances <- read.table(full.path, sep="\t", as.is=TRUE, header=TRUE)
eclip.file <- system.file(package="RnaBindingProtein", "extdata", "ENCFF565FNW.bigBed")

if(!exists("annotations")){
   all.anno <- grep("hg38", builtin_annotations(), v=TRUE)
   #aoi <- all.anno[c(1:11, 16, 19)]
   #aoi <- c("hg38_genes_3UTRs", "hg38_genes_5UTRs", "hg38_genes_cds")
   aoi <- "hg38_genes_5UTRs"
   annotations <- build_annotations(genome="hg38", annotations=aoi)
   tbl.anno <- as.data.frame(annotations)
   }

#----------------------------------------------------------------------------------------------------
find.binding.sites <- function(gene)
{
   rbp <- RnaBindingProtein$new("DDX3X", gene, eclip.file, "K562")

   tbl.roi <- subset(tbl.anno, symbol==gene)
   colnames(tbl.roi)[1] <- "chrom"
   tbl.roi$chrom <- as.character(tbl.roi$chrom)

   #tbl.bindingSites <- rbp$getBindingSites.inGenicRegions()
   #bl.5utr <- subset(tbl.bindingSites$small, regionType=="hg38_genes_5UTRs")

   meme.file <- "ddx3x-ggc.meme"
   tbl.fimo <- fimoBatch(tbl.roi[, c("chrom", "start", "end")],
                         matchThreshold=1e-3, genomeName="hg38", pwmFile=meme.file,
                         expandResultsWithMotifDb=FALSE)
   if(nrow(tbl.fimo) > 0)
      tbl.fimo$gene <- gene

   tbl.fimo

} # find.binding.sites
#----------------------------------------------------------------------------------------------------
goi <- tbl.discordances$geneSymbol
length(goi)
tbls.fimo <- lapply(goi, find.binding.sites)
tbl.all <- do.call(rbind, tbls.fimo)
