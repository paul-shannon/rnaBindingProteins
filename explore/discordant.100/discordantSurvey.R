library(RnaBindingProtein)
f <- "tbl.discordances-from-jeff-103x6.RData"
data.dir <- "~/github/TrenaProjectErythropoiesis/explore/rbp/ddx3x"
full.path <- file.path(data.dir, f)
file.exists(full.path)
tbl.discordances <- get(load(full.path))
dim(tbl.discordances)
head(tbl.discordances)

table(tbl.discordances$clear.drop.at.day.10_11)   #  ?  N  Y
                                                  #  5  8 90
clear.drop.genes <- subset(tbl.discordances, clear.drop.at.day.10_11=="Y")$Protein
no.drop.genes    <- subset(tbl.discordances, clear.drop.at.day.10_11=="N")$Protein
eclip.file <- system.file(package="RnaBindingProtein", "extdata", "ENCFF565FNW.bigBed")

#----------------------------------------------------------------------------------------------------
run.search <- function(gene)
{
    rbp <- RnaBindingProtein$new("DDX3X", gene, eclip.file, "K562")

    x <- rbp$getBindingSites.inGenicRegions()
    print(lapply(x, dim))
    tbl.small <- x$small   # unique DDX3X hits, sometimes in multiple regions
    tbl.big   <- x$big     # unique DDX3X hit/genic region pairs

    return(tbl.small)

} # run.search
#----------------------------------------------------------------------------------------------------
run.search("BACH1")

