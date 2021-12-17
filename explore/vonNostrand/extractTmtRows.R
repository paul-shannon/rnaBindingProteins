library(RUnit)
tbl.rbp <-
    read.table("RBPsFromVonNostrand_2020.txt", sep="\t", skip=2, header=FALSE, nrows=-1, fill=TRUE)[, 1:2]
dim(tbl.rbp)  # 356 2
colnames(tbl.rbp) <- c("geneSymbol", "ensg")
length(tbl.rbp$ensg)           # 356
length(unique(tbl.rbp$ensg))   # 356

head(tbl.rbp)
lapply(tbl.rbp, class)

tbl.tmt <-
    read.table("tmt-proteins.txt", sep="\t", skip=2, as.is=TRUE, nrows=-1, fill=TRUE, quote="")
dim(tbl.tmt)  #  5775   36
which(unlist(lapply(1:36, function(col) grepl("ENSG", tbl.tmt[2,col]))))  # 20

head(tbl.tmt[, 20])

tmt.ensg <- tbl.tmt[,20]
length(tmt.ensg)   # 5775
length(unique(tmt.ensg)) # 5326
tmt.ensg.trimmed <- sub("\\..*$", "", tmt.ensg)
length(tmt.ensg.trimmed)            # 5775
length(unique(tmt.ensg.trimmed))    # 5326

ensg.both <- intersect(tbl.rbp$ensg, tmt.ensg.trimmed)
indices <- match(ensg.both, tmt.ensg.trimmed)
stopifnot(all(tmt.ensg.trimmed[indices] %in% tbl.rbp$ensg))

tmt.rows <- unique(sort(unlist(lapply(ensg.both, function(ensg.no.suffix) grep(ensg.no.suffix, tbl.tmt[,20])))))
length(tmt.rows)

tbl.tmt.rbp <- tbl.tmt[tmt.rows,]
dim(tbl.tmt.rbp)   # 327  36

ensg.dups <- tbl.tmt.rbp[which(duplicated(tbl.tmt.rbp[, 20])), 20]
tbl.dups <- subset(tbl.tmt.rbp, ensg %in% ensg.dups)[, c("Accession", "ensg")]
new.order <- order(tbl.dups$ensg)
tbl.dups <- tbl.dups[new.order,]
row.names(tbl.dups) <- NULL

    #----------------------------------------------------------------------------------
    #  verify: are all no-suffix ensgs found in the 20th column of with-suffix ensgs?
    #----------------------------------------------------------------------------------

x <- unlist(lapply(ensg.both, function(ensg) as.logical(grep(ensg, tbl.tmt.rbp[,20]))))
table(x)

    #-------------------------------------------------------------------------------------
    #  verify: with-suffix ensg and protein from each successive row of the new table
    #  make sure that that pairing is also found in the original tbl.tmt
    #-------------------------------------------------------------------------------------

for(row in seq_len(nrow(tbl.tmt.rbp))){
    ensg.suf <- tbl.tmt.rbp[row, 20]
    protein  <- tbl.tmt.rbp[row, 2]
    indices.orig <- grep(ensg.suf, tbl.tmt[,20])
    checkTrue(protein %in% tbl.tmt[indices.orig, 2])
    }


   #-------------------------------------------
   # try to reconstruct colnames from tbl.tmt
   #-------------------------------------------

s <- 'Protein FDR Confidence: Combined	Accession	Description	Exp. q-value: Combined	Sum PEP Score	Coverage [%]	# Peptides	# PSMs	# Unique Peptides	# AAs	MW [kDa]	calc. pI	Score Sequest HT: Sequest HT	# Peptides (by Search Engine): Sequest HT	Biological Process	Cellular Component	Molecular Function	Pfam IDs	Entrez Gene ID	Ensembl Gene ID	Chromosome	WikiPathways	Gene Symbol	KEGG Pathways	Reactome Pathways	# Protein Pathway Groups	# Razor Peptides		"126, Day 2"	"127N, Day 7.5"	"128N, Day 8"	"129C, Day 10"	"130N, Day 12"	"131, Day 14"	# Protein Groups	Modifications'
nchar(s)
coi <- strsplit(s, "\t")[[1]]
length(coi)
dim(tbl.tmt)
dim(tbl.tmt.rbp)

colnames(tbl.tmt.rbp) <- coi
write.table(tbl.tmt.rbp, file="rbp.tmt.tsv", sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)


