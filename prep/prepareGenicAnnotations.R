# the larger goal of this package is to  find all binding sites for DDX3X (and soon, other RBPs also)
# which overlap 5’ or 3’ UTRs. Interesting that the first BACH1 exon is annotated
# as a 5’UTR.  after discussing with Jeff, this resolution proposed:  annotate exon intersections
# (T/F) for each UTR/DDX3x intersection.  should be infrequent.  the import or inconsequentiality
# of this can be figured out later

# here we prepare a data.frame of UTR annotations, with exon intersect noted as an extra column

library(annotatr)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
#------------------------------------------------------------------------------------------------------------------------
# try getting exons and utrs direct from the txdb:
#  bach1 <- "571"
#  transcriptsBy(txdb, by="gene")[[bach1]]
#  tr.ids <- transcriptsBy(txdb, by="gene")[["571"]]$tx_name
#  fiveUTRsByTranscript(txdb, use.names=TRUE)[["ENST00000420190.6"]]
#  lapply(tr.ids, function(tr.id) fiveUTRsByTranscript(txdb, use.names=TRUE)[tr.id])
#  x <- fiveUTRsByTranscript(txdb, use.names=TRUE)    # 86618
#  tr.ids <- intersect(tr.ids, names(x))
#  gr <- unlist(x[tr.ids])
#  tbl.bach1.5utr <- as.data.frame(gr, row.names=NULL)
#  tx.ids <- transcriptsBy(txdb, by="gene")[[bach1]] #$`tx_name`
#  tbl.5utrs.bytx <- as.data.frame(x)
#  status (19 nov, 3pm:  every utr is treated as an exon.  some 5'utrs bear no
#  apparent relation to tss.  next up: try biomart.
#  this may be relevant:  https://www.biostars.org/p/474912/
#     Retrieve 5'UTR sequences for ensembl_transcript_id's with unique start/end positions
#
# see slides: ddx3x 55, 56, with jeff concluding that all the annotatr UTRs are worthwhile

hg38.annotation.types <- grep("hg38", builtin_annotations(), v=TRUE)
hg38.utrs <- grep("UTR", hg38.annotation.types, value=TRUE)
#hg38.exons <- "hg38_genes_exons" # grep("exon", hg38.annotation.types, value=TRUE, ignore.case=TRUE)
#hg38.introns <- "hg38_genes_introns"  #grep("intron", hg38.annotation.types, value=TRUE, ignore.case=TRUE)
annotations <- build_annotations(genome="hg38", annotations=c(hg38.utrs))
length(annotations)  # 325799
length(unique(annotations$symbol)) # [1] 19024

tbl.anno <- as.data.frame(annotations)
colnames(tbl.anno)[1] <- "chrom"
tbl.anno$chrom <- as.character(tbl.anno$chrom)

dim(tbl.anno)   # 325799 10
head(tbl.anno)

dups <- which(duplicated(tbl.anno[, c("chrom", "start", "end", "type")]))
length(dups)   # 95527
tbl.anno <- tbl.anno[-dups,]
dim(tbl.anno)    # 230272 10
table(tbl.anno$type)   #   hg38_genes_3UTRs hg38_genes_5UTRs
                       #          123965           106307

tbl.utrs <- tbl.anno
save(tbl.utrs, file="../inst/extdata/UTRS-ucsc.RData")

    #---------------------------------------------------
    # just exploratory & now abandoned code beyond here
    #---------------------------------------------------

tbl.5utr  <- subset(tbl.anno, type=="hg38_genes_5UTRs")
tbl.3utr  <- subset(tbl.anno, type=="hg38_genes_3UTRs")

    #--------------------------------------------------
    #  take a look
    #--------------------------------------------------

igv <- start.igv("BACH1", "hg38")
zoomOut(igv)
roi <- getGenomicRegion(igv)

    # exons in red

tbl.track <- as.data.frame(exonsBy(txdb, "gene")[["571"]])

#tbl.track <- subset(tbl.exons.txdb, chrom==roi$chrom & start >= roi$start & end <= roi$end)
dim(tbl.track)
track <- DataFrameAnnotationTrack("exons", tbl.track, color="red")
displayTrack(igv, track)

   # 5UTR in blue

tbl.track <- subset(tbl.5utr, chrom==roi$chrom & start >= roi$start & end <= roi$end)
dim(tbl.track)
track <- DataFrameAnnotationTrack("5UTR", tbl.track, color="blue")
displayTrack(igv, track)

   # 3UTR in green

tbl.track <- subset(tbl.3utr, chrom==roi$chrom & start >= roi$start & end <= roi$end)
dim(tbl.track)
track <- DataFrameAnnotationTrack("3UTR", tbl.track, color="darkgreen")
displayTrack(igv, track)




tbl.ov <- as.data.frame(findOverlaps(GRanges(tbl.5utr), GRanges(tbl.exons)))
head(tbl.ov)
length(unique(tbl.ov$queryHits))  # 4606/106307  ~4%

tbl.ov <- as.data.frame(findOverlaps(GRanges(tbl.5utr), GRanges(tbl.3utr)))
head(tbl.ov)
length(unique(tbl.ov$queryHits))  # 4606/106307  ~4%



# examine the puzzling 5'utr in an exon of BACH1
# chr21:29321200-29321600
tbl.small <- subset(tbl.anno, chrom=="chr21" & start >= 29321200 & end <= 29321600)

#         chrom    start      end width strand           id             tx_id gene_id    symbol             type
# 146977  chr21 29321221 29321280    60      +  5UTR:146977 ENST00000548219.5  193629 LINC00189 hg38_genes_5UTRs
# 146981  chr21 29321221 29321280    60      +  5UTR:146981 ENST00000546469.5     571     BACH1 hg38_genes_5UTRs
# 146984  chr21 29321221 29321280    60      +  5UTR:146984 ENST00000547141.5     571     BACH1 hg38_genes_5UTRs
# 146987  chr21 29321221 29321280    60      +  5UTR:146987 ENST00000550131.5     571     BACH1 hg38_genes_5UTRs
# 146989  chr21 29321221 29321280    60      +  5UTR:146989 ENST00000286800.8     571     BACH1 hg38_genes_5UTRs
# 146991  chr21 29321221 29321280    60      +  5UTR:146991 ENST00000399921.5     571     BACH1 hg38_genes_5UTRs
# 146993  chr21 29321221 29321280    60      +  5UTR:146993 ENST00000548467.1     571     BACH1 hg38_genes_5UTRs
# 146995  chr21 29321221 29321280    60      +  5UTR:146995 ENST00000451655.5     571     BACH1 hg38_genes_5UTRs
# 146997  chr21 29321221 29321280    60      +  5UTR:146997 ENST00000447177.5     571     BACH1 hg38_genes_5UTRs
# 146999  chr21 29321221 29321280    60      +  5UTR:146999 ENST00000435072.1     571     BACH1 hg38_genes_5UTRs
# 1492652 chr21 29321221 29321355   135      + exon:1334651 ENST00000548219.5  193629 LINC00189 hg38_genes_exons
# 1492656 chr21 29321221 29321355   135      + exon:1334655 ENST00000546469.5     571     BACH1 hg38_genes_exons
# 1492659 chr21 29321221 29321355   135      + exon:1334658 ENST00000547141.5     571     BACH1 hg38_genes_exons
# 1492662 chr21 29321221 29321355   135      + exon:1334661 ENST00000550131.5     571     BACH1 hg38_genes_exons
# 1492665 chr21 29321221 29321514   294      + exon:1334664 ENST00000286800.8     571     BACH1 hg38_genes_exons
# 1492670 chr21 29321221 29321514   294      + exon:1334669 ENST00000399921.5     571     BACH1 hg38_genes_exons
# 1492675 chr21 29321221 29321355   135      + exon:1334674 ENST00000548467.1     571     BACH1 hg38_genes_exons
# 1492677 chr21 29321221 29321514   294      + exon:1334676 ENST00000451655.5     571     BACH1 hg38_genes_exons
# 1492680 chr21 29321221 29321514   294      + exon:1334679 ENST00000447177.5     571     BACH1 hg38_genes_exons
# 1492683 chr21 29321221 29321514   294      + exon:1334682 ENST00000435072.1     571     BACH1 hg38_genes_exons

igv <- start.igv("BACH1")
showGenomicRegion(igv, "chr21:29321200-29321600")
tbl.track <- unique(subset(tbl.small, type=="hg38_genes_5UTRs")[, c("chrom", "start", "end")])
track <- DataFrameAnnotationTrack("5'utr", tbl.track, color="red")
displayTrack(igv, track)

tbl.track <- unique(subset(tbl.small, type=="hg38_genes_exons")[, c("chrom", "start", "end")])
track <- DataFrameAnnotationTrack("exons", tbl.track, color="blue")
displayTrack(igv, track)

# no promoters in this region, so that cannot be used to identify a TSS

tbl.track <- unique(subset(tbl.small, type=="hg38_genes_promoters")[, c("chrom", "start", "end")])
dim(tbl.track)
track <- DataFrameAnnotationTrack("promoters", tbl.track, color="black")
displayTrack(igv, track)
