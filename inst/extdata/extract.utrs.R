library(GenomicRanges)
library(annotatr)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
hg38.annotation.types <- grep("hg38", builtin_annotations(), v=TRUE)

hg38.utrs <- grep("UTR", hg38.annotation.types, value=TRUE)
hg38.utrs   # "hg38_genes_5UTRs" "hg38_genes_3UTRs"
annotations <- build_annotations(genome="hg38", annotations=c(hg38.utrs))
tbl.anno <- as.data.frame(annotations)
colnames(tbl.anno)[1] <- "chrom"
tbl.anno$chrom <- as.character(tbl.anno$chrom)
nrow(tbl.anno)
roi <- list(chrom="chr21", start=29320940, end=29321745, string="chr21:29,320,940-29,321,745")
subset(tbl.anno, chrom==roi$chrom & start > roi$start & end < roi$end)

dups <- which(duplicated(tbl.anno[, c("chrom", "start", "end", "gene_id", "type")]))
length(dups)
tbl.anno <- tbl.anno[-dups,]
dim(tbl.anno)   # 231324  10
save(tbl.anno, file="UTRS-ucsc.RData")
# check BACH1, first exon, where a strong DDX3X peak occurs, and where I once found a 5'UTR


# length(annotations)  # 325799
# length(unique(annotations$symbol)) # [1] 19024
# tbl.anno <- as.data.frame(annotations)
# dups <- which(duplicated(tbl.anno[, c("chrom", "start", "end", "type")]))
# length(dups)   # 95527
# tbl.anno <- tbl.anno[-dups,]
# dim(tbl.anno)    # 230272 10
#table(tbl.anno$type)   #   hg38_genes_3UTRs hg38_genes_5UTRs
                        #          123965           106307

# tbl.anno <- get(load("UTRS-ucsc.RData"))
# dim(tbl.anno)  # 230272 10
# gr.hg38.utrs <- get(load("anno.hg38.utrs.RData"))
# length(gr.hg38.utrs)   # 325799
# gr.hg38.utrs
# tbl.utrs <- as.data.frame(gr.hg38.utrs)
# dim(tbl.utrs) #
