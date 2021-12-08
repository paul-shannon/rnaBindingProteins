tbl.discordants <-
    read.table("group-1-lateDiscordants.txt", sep="\t", as.is=TRUE, header=TRUE, nrow=-1)

tbl.concordants <-
    read.table("group-2-lateConcordants.txt", sep="\t", as.is=TRUE, header=TRUE, nrow=-1)

dim(tbl.discordants)   # 42 5
dim(tbl.concordants)   # 38 5

    # no geneSymbol shared between sets:
stopifnot(length(intersect(tbl.discordants$geneSymbol, tbl.concordants$geneSymbol)) == 0)

all(colnames(tbl.discordants) == colnames(tbl.concordants))

table(tbl.discordants[, 5])
table(tbl.concordants[, 5])

coi <- c("protein", "geneSymbol", "clearProteinDrop", "tmtAgreement", "rnaProteinDiscord")
colnames(tbl.discordants) <- coi
colnames(tbl.concordants) <- coi

tbl.both <- rbind(tbl.discordants, tbl.concordants)
tbl.both$rnaProteinDiscord[tbl.both$rnaProteinDiscord=="Y"] <- "T"
tbl.both$rnaProteinDiscord[tbl.both$rnaProteinDiscord=="N"] <- "F"

tbl.both$rnaProteinDiscord <- as.logical(tbl.both$rnaProteinDiscord)

table(tbl.both$rnaProteinDiscord)   # FALSE  TRUE

                                    #    38    42

tbl.both$clearDrop[tbl.both$clearDrop=="Y"] <- "T"
tbl.both$clearDrop[tbl.both$clearDrop=="N"] <- "F"
tbl.both$clearDrop <- as.logical(tbl.both$clearDrop)

dim(tbl.both) # 80 5

   #------------------------------------------------------------
   # but 103 genes/proteins were used before.  what of those 23?
   #------------------------------------------------------------

f <- "srm-103-late-discordance-geneSymbolsAdded.tsv"

data.dir <- "~/github/TrenaProjectErythropoiesis/explore/rbp/ddx3x"
full.path <- file.path(data.dir, f)
file.exists(full.path)
tbl.old <- read.table(full.path, sep="\t", header=TRUE, as.is=TRUE)

setdiff(tbl.old$geneSymbol, tbl.both$geneSymbol)

  #  [1] "POU2F1"  "STAT1"   "STAT2"   "ETF1"    "EP300"   "CEBPB"   "CBFA2T3"
  #  [8] "STAT3"   "GATAD2B" "HDAC2"   "JUN"     "JUND"    "DR1"     "RXRB"
  # [15] "SAP130"  "TRIM28"  "UBTF"    "FOXO3"   "NRF1"    "NFKB1"   "TTF2"



with(tbl.both, table(clearDrop, rnaProteinDiscord))
    #          rnaProteinDiscord
    # clearDrop FALSE TRUE
    #     FALSE     6    0
    #     TRUE     32   42



save(tbl.both, file="erytrhoGenesProteinsDropAndDiscordance.RData")
tbl <- get(load("~/github/rnaBindingProteins/explore/discordant.100/spreadsheetsFromJeff/erytrhoGenesProteinsDropAndDiscordance.RData")
# tbl.both <- get(load(

dim(tbl.discordants)          # 42 5
colnames(tbl.concordants)
dim(tbl.concordants)          # 38 5

head(tbl.discordants)


table(tbl.discordants[,3])
table(tbl.discordants[,4])


tbl.dis <- tbl.discordants[, c("Protein", "geneSymbol", "clear.drop.at.day.10_11,

