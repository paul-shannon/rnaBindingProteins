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
go.drop.discord.vs.noDiscord <- function()
{
  data.dir <- "~/github/rnaBindingProteins/explore/discordant.100/spreadsheetsFromJeff"
  file <- "erytrhoGenesProteinsDropAndDiscordance.RData"
  tbl <- get(load(file.path(data.dir, file)))
  dim(tbl) # 80 5

       #---------------------------------------------------------
       #  examine the 32 "clearDrop" with no protein/rna discord
       #  will we see fewer binding sites here?
       #----------------------------------------------------------

  goi.noDiscord <- subset(tbl, clearProteinDrop & !rnaProteinDiscord)$geneSymbol
  length(goi.noDiscord)  # 32

  tbls.goi.noDiscord <- lapply(goi.noDiscord, run.search)
  tbl.noDiscord <- do.call(rbind, tbls.goi.noDiscord)

  tbl.noDiscord.3utr <- tbl.noDiscord[grep("3UTR", tbl.noDiscord$regionType),]    # 66
  tbl.noDiscord.5utr <- tbl.noDiscord[grep("5UTR", tbl.noDiscord$regionType),]    # 184
  tbl.noDiscord.cds <- subset(tbl.noDiscord, regionType=="hg38_genes_cds")   # 145

  mean(as.numeric(lapply(unique(tbl.noDiscord.3utr$targetGene),
                         function(gene) mean(subset(tbl.noDiscord.3utr, targetGene==gene)$score))))
      # 1.879
  mean(as.numeric(lapply(unique(tbl.noDiscord.5utr$targetGene),
                         function(gene) mean(subset(tbl.noDiscord.5utr, targetGene==gene)$score))))
      # 6.62

  mean(as.numeric(lapply(unique(tbl.noDiscord.cds$targetGene),
                         function(gene) mean(subset(tbl.noDiscord.cds, targetGene==gene)$score))))
      # 2.48

  length(unique(tbl.noDiscord$targetGene))       # 29/32  91%
  length(unique(tbl.noDiscord.5utr$targetGene))  # 26/32  81%
  length(unique(tbl.noDiscord.3utr$targetGene))  # 15/32  47%
  length(unique(tbl.noDiscord.cds$targetGene))   # 20/32  63$

       #---------------------------------------------------------
       #  examine the 42 "clearDrop" with protein/rna discord
       #  will we see more binding sites here?
       #----------------------------------------------------------

  goi.discord <- subset(tbl, clearProteinDrop & rnaProteinDiscord)$geneSymbol
  length(goi.discord)  # 42

  tbls.goi.discord <- lapply(goi.discord, run.search)
  tbl.discord <- do.call(rbind, tbls.goi.discord)

  fivenum(as.numeric(lapply(tbls.goi.discord, nrow)))  # [1] 0  2   3.5   6   23

  tbl.discord.3utr <- tbl.discord[grep("3UTR", tbl.discord$regionType),]    # 66
  tbl.discord.5utr <- tbl.discord[grep("5UTR", tbl.discord$regionType),]    # 184
  tbl.discord.cds <- subset(tbl.discord, regionType=="hg38_genes_cds")   # 145

  mean(as.numeric(lapply(unique(tbl.discord.3utr$targetGene),
                         function(gene) mean(subset(tbl.discord.3utr, targetGene==gene)$score))))
      # 1.4
  mean(as.numeric(lapply(unique(tbl.discord.5utr$targetGene),
                         function(gene) mean(subset(tbl.discord.5utr, targetGene==gene)$score))))
      # 4.93

  mean(as.numeric(lapply(unique(tbl.discord.cds$targetGene),
                         function(gene) mean(subset(tbl.discord.cds, targetGene==gene)$score))))
      # 3.05

  length(unique(tbl.discord$targetGene))       #  39/42  93%
  length(unique(tbl.discord.5utr$targetGene))  #  36/42  86%
  length(unique(tbl.discord.3utr$targetGene))  #  13/42  31%
  length(unique(tbl.discord.cds$targetGene))   #  23/42  55%

} # go.drop.discord.vs.noDiscord
#----------------------------------------------------------------------------------------------------
go.random.1000 <- function()
{

      #--------------------------
      #  get 1000 geneSymbols
      #--------------------------
   gene <- "BACH1"
   rbp <- RnaBindingProtein$new("DDX3X", gene, eclip.file, "K562")
   tbl.anno <- rbp$getAllGenicAnnotations()
   genes <- sort(unique(tbl.anno$symbol))
   deleters <- which(is.na(genes))
   deleters <- c(deleters, grep("^NA$", genes, ignore.case=TRUE))
   if(length(deleters) > 0)
       genes <- genes[-deleters]
   length(genes)

   set.seed(17)
   set.seed(37)
   random.genes <- sort(sample(genes, size=1000))
   tbls.random <- lapply(random.genes, run.search)
   length(tbls.random)
   tbl.random <- do.call(rbind, tbls.random)
   dim(tbl.random)  # 1595 8
   length(unique(tbl.random$targetGene)) # [1] 423, 394

   length(grep("3UTR", tbl.random$regionType)) # 370, 246
   length(grep("5UTR", tbl.random$regionType)) # 779, 696
   nrow(subset(tbl.random, regionType=="hg38_genes_cds"))   # 501, 451

  nrow(tbl.random)/length(unique(tbl.random$targetGene))     # 3.77

  mean(as.numeric(lapply(unique(tbl.random$targetGene),
                         function(gene) mean(subset(tbl.random, targetGene==gene)$score))))
      # [1] 6.73

  dim(tbl.random)  # 1595 DDX3X/utr-cds hits
  tbl.random.3utr <- tbl.random[grep("3UTR", tbl.random$regionType),]    # 66
  tbl.random.5utr <- tbl.random[grep("5UTR", tbl.random$regionType),]    # 184
  tbl.random.cds <- subset(tbl.random, regionType=="hg38_genes_cds")   # 145

  mean(as.numeric(lapply(unique(tbl.random.3utr$targetGene),
                         function(gene) mean(subset(tbl.random.3utr, targetGene==gene)$score))))
      # 2.46, 2.13
  mean(as.numeric(lapply(unique(tbl.random.5utr$targetGene),
                         function(gene) mean(subset(tbl.random.5utr, targetGene==gene)$score))))
      # 11.32, 9.21

  mean(as.numeric(lapply(unique(tbl.random.cds$targetGene),
                         function(gene) mean(subset(tbl.random.cds, targetGene==gene)$score))))
      # 5.87, 7.02

  length(unique(tbl.random$targetGene))       # 77/90    13 had no binding sites
  length(unique(tbl.random.5utr$targetGene))  # 68/90
  length(unique(tbl.random.3utr$targetGene))  # 32/90
  length(unique(tbl.random.cds$targetGene))   # 50/90

} # go.random.1000
#----------------------------------------------------------------------------------------------------


