library(RnaBindingProtein)
f <- "srm-103-late-discordance-geneSymbolsAdded.tsv"

data.dir <- "~/github/TrenaProjectErythropoiesis/explore/rbp/ddx3x"
full.path <- file.path(data.dir, f)
file.exists(full.path)
tbl.discordances <- read.table(full.path, sep="\t", as.is=TRUE, header=TRUE)
dim(tbl.discordances)
head(tbl.discordances)

table(tbl.discordances$clear.drop.at.day.10_11)   #  ?  N  Y
                                                  #  5  8 90
clear.drop.genes <- subset(tbl.discordances, clear.drop.at.day.10_11=="Y")$geneSymbol
no.drop.genes    <- subset(tbl.discordances, clear.drop.at.day.10_11=="N")$geneSymbol
eclip.file <- system.file(package="RnaBindingProtein", "extdata", "ENCFF565FNW.bigBed")

#----------------------------------------------------------------------------------------------------
run.search <- function(gene)
{
    printf("--- run.search(%s)", gene)

    rbp <- RnaBindingProtein$new("DDX3X", gene, eclip.file, "K562")
    x <- rbp$getBindingSites.inGenicRegions()
    print(lapply(x, dim))
    tbl.small <- x$small   # unique DDX3X hits, sometimes in multiple regions
    tbl.big   <- x$big     # unique DDX3X hit/genic region pairs

    return(tbl.small)

} # run.search
#----------------------------------------------------------------------------------------------------
go <- function()
{

  run.search("BACH1")
  run.search("BCL11A")
  length(clear.drop.genes) # 90
  x.drop <- lapply(clear.drop.genes, run.search)
  tbl.drop <- do.call(rbind,x.drop)
  dim(tbl.drop)  # 386 8
  length(unique(tbl.drop$targetGene))   # 77/90

  length(no.drop.genes)   # 8
  x.noDrop <- lapply(no.drop.genes, run.search)
  tbl.noDrop <- do.call(rbind, x.noDrop)
  dim(tbl.noDrop)  # 10 8

  x.noDrop <- lapply(no.drop.genes, run.search)
  length(unique(tbl.noDrop$targetGene))   # 4/8
  nrow(tbl.noDrop)/length(unique(tbl.noDrop$targetGene)) # 2.5 binding sites per gene
  nrow(tbl.drop)/length(unique(tbl.drop$targetGene))     # 5.0 binding sites per gene

  mean(as.numeric(lapply(unique(tbl.drop$targetGene),
                         function(gene) mean(subset(tbl.drop, targetGene==gene)$score))))
      # [1] 3.848331
  mean(as.numeric(lapply(unique(tbl.noDrop$targetGene),
                         function(gene) mean(subset(tbl.noDrop, targetGene==gene)$score))))
      # [1] 3.583512

  dim(tbl.drop)  # 386 DDX3X/utr-cds hits
  tbl.drop.3utr <- tbl.drop[grep("3UTR", tbl.drop$regionType),]    # 66
  tbl.drop.5utr <- tbl.drop[grep("5UTR", tbl.drop$regionType),]    # 184
  tbl.drop.cds <- subset(tbl.drop, regionType=="hg38_genes_cds")   # 145


  tbl.noDrop.3utr <- tbl.noDrop[grep("3UTR", tbl.noDrop$regionType),]  # 2
  tbl.noDrop.5utr <- tbl.noDrop[grep("5UTR", tbl.noDrop$regionType),]  # 6
  tbl.noDrop.cds <- subset(tbl.noDrop, regionType=="hg38_genes_cds")   # 2

  mean(as.numeric(lapply(unique(tbl.drop.3utr$targetGene),
                         function(gene) mean(subset(tbl.drop.3utr, targetGene==gene)$score))))
      # 2.05
  mean(as.numeric(lapply(unique(tbl.drop.5utr$targetGene),
                         function(gene) mean(subset(tbl.drop.5utr, targetGene==gene)$score))))
      # 5.58

  mean(as.numeric(lapply(unique(tbl.drop.cds$targetGene),
                         function(gene) mean(subset(tbl.drop.cds, targetGene==gene)$score))))
      # 2.74


  mean(as.numeric(lapply(unique(tbl.drop.cds$targetGene),
                         function(gene) mean(subset(tbl.drop.cds, targetGene==gene)$score))))
      # 2.74

  mean(as.numeric(lapply(unique(tbl.noDrop.3utr$targetGene),
                         function(gene) mean(subset(tbl.noDrop.3utr, targetGene==gene)$score))))
      # 1.27
  mean(as.numeric(lapply(unique(tbl.noDrop.5utr$targetGene),
                         function(gene) mean(subset(tbl.noDrop.5utr, targetGene==gene)$score))))
      # 6.55
  mean(as.numeric(lapply(unique(tbl.noDrop.cds$targetGene),
                         function(gene) mean(subset(tbl.noDrop.cds, targetGene==gene)$score))))
      # 1.36


  length(unique(tbl.noDrop$targetGene))       # 4/8    4 had no binding sites
  length(unique(tbl.noDrop.5utr$targetGene))  # 4/8
  length(unique(tbl.noDrop.3utr$targetGene))  # 1/8
  length(unique(tbl.noDrop.cds$targetGene))   # 2/8

  length(unique(tbl.drop$targetGene))       # 77/90    13 had no binding sites
  length(unique(tbl.drop.5utr$targetGene))  # 68/90
  length(unique(tbl.drop.3utr$targetGene))  # 32/90
  length(unique(tbl.drop.cds$targetGene))   # 50/90

} # go
#----------------------------------------------------------------------------------------------------
go.random.1000 <- function()
{

      #--------------------------
      #  get 1000 geneSymbols
      #--------------------------
   gene <- "BACH1"
   rbp <- RnaBindingProtein$new("DDX3X", gene, eclip.file, "K562")
   tbl.anno <- rbp$getGenicRegions()
   dim(tbl.anno)

    x <- rbp$getBindingSites.inGenicRegions()
    print(lapply(x, dim))
    tbl.small <- x$small   # unique DDX3X hits, sometimes in multiple regions
    tbl.big   <- x$big     # unique DDX3X hit/genic region pairs


  run.search("BACH1")
  run.search("BCL11A")
  length(clear.drop.genes) # 90
  x.drop <- lapply(clear.drop.genes, run.search)
  tbl.drop <- do.call(rbind,x.drop)
  dim(tbl.drop)  # 386 8
  length(unique(tbl.drop$targetGene))   # 77/90

  length(no.drop.genes)   # 8
  x.noDrop <- lapply(no.drop.genes, run.search)
  tbl.noDrop <- do.call(rbind, x.noDrop)
  dim(tbl.noDrop)  # 10 8

  x.noDrop <- lapply(no.drop.genes, run.search)
  length(unique(tbl.noDrop$targetGene))   # 4/8
  nrow(tbl.noDrop)/length(unique(tbl.noDrop$targetGene)) # 2.5 binding sites per gene
  nrow(tbl.drop)/length(unique(tbl.drop$targetGene))     # 5.0 binding sites per gene

  mean(as.numeric(lapply(unique(tbl.drop$targetGene),
                         function(gene) mean(subset(tbl.drop, targetGene==gene)$score))))
      # [1] 3.848331
  mean(as.numeric(lapply(unique(tbl.noDrop$targetGene),
                         function(gene) mean(subset(tbl.noDrop, targetGene==gene)$score))))
      # [1] 3.583512

  dim(tbl.drop)  # 386 DDX3X/utr-cds hits
  tbl.drop.3utr <- tbl.drop[grep("3UTR", tbl.drop$regionType),]    # 66
  tbl.drop.5utr <- tbl.drop[grep("5UTR", tbl.drop$regionType),]    # 184
  tbl.drop.cds <- subset(tbl.drop, regionType=="hg38_genes_cds")   # 145


  tbl.noDrop.3utr <- tbl.noDrop[grep("3UTR", tbl.noDrop$regionType),]  # 2
  tbl.noDrop.5utr <- tbl.noDrop[grep("5UTR", tbl.noDrop$regionType),]  # 6
  tbl.noDrop.cds <- subset(tbl.noDrop, regionType=="hg38_genes_cds")   # 2

  mean(as.numeric(lapply(unique(tbl.drop.3utr$targetGene),
                         function(gene) mean(subset(tbl.drop.3utr, targetGene==gene)$score))))
      # 2.05
  mean(as.numeric(lapply(unique(tbl.drop.5utr$targetGene),
                         function(gene) mean(subset(tbl.drop.5utr, targetGene==gene)$score))))
      # 5.58

  mean(as.numeric(lapply(unique(tbl.drop.cds$targetGene),
                         function(gene) mean(subset(tbl.drop.cds, targetGene==gene)$score))))
      # 2.74


  mean(as.numeric(lapply(unique(tbl.drop.cds$targetGene),
                         function(gene) mean(subset(tbl.drop.cds, targetGene==gene)$score))))
      # 2.74

  mean(as.numeric(lapply(unique(tbl.noDrop.3utr$targetGene),
                         function(gene) mean(subset(tbl.noDrop.3utr, targetGene==gene)$score))))
      # 1.27
  mean(as.numeric(lapply(unique(tbl.noDrop.5utr$targetGene),
                         function(gene) mean(subset(tbl.noDrop.5utr, targetGene==gene)$score))))
      # 6.55
  mean(as.numeric(lapply(unique(tbl.noDrop.cds$targetGene),
                         function(gene) mean(subset(tbl.noDrop.cds, targetGene==gene)$score))))
      # 1.36


  length(unique(tbl.noDrop$targetGene))       # 4/8    4 had no binding sites
  length(unique(tbl.noDrop.5utr$targetGene))  # 4/8
  length(unique(tbl.noDrop.3utr$targetGene))  # 1/8
  length(unique(tbl.noDrop.cds$targetGene))   # 2/8

  length(unique(tbl.drop$targetGene))       # 77/90    13 had no binding sites
  length(unique(tbl.drop.5utr$targetGene))  # 68/90
  length(unique(tbl.drop.3utr$targetGene))  # 32/90
  length(unique(tbl.drop.cds$targetGene))   # 50/90

} # go.random.1000
#----------------------------------------------------------------------------------------------------


