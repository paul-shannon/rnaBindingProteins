# hg38.utrs <- grep("UTR", hg38.annotation.types, value=TRUE)
# annotations <- build_annotations(genome="hg38", annotations=c(hg38.utrs))
# length(annotations)  # 325799
# length(unique(annotations$symbol)) # [1] 19024
# tbl.anno <- as.data.frame(annotations)
# dups <- which(duplicated(tbl.anno[, c("chrom", "start", "end", "type")]))
# length(dups)   # 95527
# tbl.anno <- tbl.anno[-dups,]
# dim(tbl.anno)    # 230272 10
#table(tbl.anno$type)   #   hg38_genes_3UTRs hg38_genes_5UTRs
                        #          123965           106307

tbl.anno <- get(load("UTRS-ucsc.RData"))
dim(tbl.anno)  # 230272 10
gr.hg38.utrs <- get(load("anno.hg38.utrs.RData"))
length(gr.hg38.utrs)   # 325799
gr.hg38.utrs
tbl.utrs <- as.data.frame(gr.hg38.utrs)
dim(tbl.utrs) #
