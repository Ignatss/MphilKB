# Analysis of circadian RNA-seq data from Marchantia Polymorpha

# DATA PRE-PROCESSING
# Set WD where all the files are

options(stringsAsFactors = FALSE)

# Hmisc - package for Co-Expression Analysis
library('Hmisc') # rcorr, %nin%
library('dplyr')

# Import genes of interest
TF.list   <- read.csv('TF-table.csv')$Gene_ID_New
Meri.list <- read.csv("meristeme-genes.csv", header=T)$id



#==================================================
# Cleaned-up PlantRegMap output (promoter mining)
#==================================================

fimo.counts <- read.csv('fimo-clean-counts.csv')

# My hit-count table, `fimo.counts`, uses two different gene name conventions
# for the TF (X.pattern.name) and the target (sequence.name):
#    TF:      "Mapoly0093s0032"
#    target:  "SCR2"
# neither of which suits us, because we need Tak5.1 notation like "Mp5g1110"

# omit ARF3 which i dropped and GRAS1 which is same as SCR2
# (GRAS1 is SCR2, mistake on my part in pre-prearing data, getting rid of doubles)

fimo.counts <- fimo.counts[
    which(fimo.counts$sequence.name %nin% c('GRAS1', 'ARF3')),
]

# TF.list contains correctly mapped names (e.g. "Mapoly0093s0032" -> "Mp5g11100"),
# so I reuse it for the TF column

name_id_map <- unique(read.csv('TF-table.csv')[,c('Gene_ID', 'Gene_ID_New')])
rownames(name_id_map) <- name_id_map$Gene_ID
fimo.counts$TF_id <- name_id_map[fimo.counts$X.pattern.name, 'Gene_ID_New']

# convert the short names used for the target (e.g. "SCR2" -> "Mp1g20490")

sequence_name_map <- read.csv('meristeme-genes.csv')
rownames(sequence_name_map) <- sequence_name_map$name
fimo.counts$target_id <- sequence_name_map[fimo.counts$sequence.name, 'id']

# some IDs are missing from our name-mapping tables and end up as NA,
# they were deleted as the genome version was updated so i drop them

fimo.counts <- na.omit(fimo.counts)





#============================================
# Create coexpression network based on normalized edge counts
#============================================

# CHECK Co-EXPRESSION PATTERNS AMONG GRAS GENES

# Import Total Data for Gene Counts
normalized.counts <- read.csv("normalized-counts-edgeR.csv", row.names='Genes')

#Get Co-expression Data for Genes of Interest.
#Combine Parametric/Non-Parametric Methods for best results
normalized.counts.rcorr <- rcorr(t(as.matrix(normalized.counts)))


# Flatten a correlation matrix
# (convert from matrix to dataframe)
flatten.corr.matrix <- function(cormat, pmat) {
    ut <- upper.tri(cormat)
    data.frame(
        row = rownames(cormat)[row(cormat)[ut]],
        column = rownames(cormat)[col(cormat)[ut]],
        cor = cormat[ut],
        p = pmat[ut]
    )
}

normalized.counts.rcorr.flat <- flatten.corr.matrix(
    normalized.counts.rcorr$r,
    normalized.counts.rcorr$P
)


Meri.TF.list <- unique(c(Meri.list, TF.list))
filtered.network <- filter(
    normalized.counts.rcorr.flat,
    row %in% Meri.list | column %in% Meri.list,
    row %in% Meri.TF.list, column %in% Meri.TF.list,
    p <= 0.001,
    cor <= -0.8 | cor >= 0.8
)

names(filtered.network)[names(filtered.network) == 'column'] <- 'TF_id'
names(filtered.network)[names(filtered.network) == 'row']    <- 'target_id'




# ===========================================
# Annotate the graph - does the graph edge appear in fimo.counts
# ===========================================

pairs <- paste(fimo.counts$TF_id, fimo.counts$target_id)
filtered.network$in_fimo <- sapply(
    paste(filtered.network$TF_id, filtered.network$target_id),
    function(x) x %in% pairs
)




# ===========================================
# convert the short names used for the target (e.g. "Mp1g20490" -> "SCR2")
# ===========================================


names <- read.csv('gene-names.tsv', sep="\t", na.strings = c(''));
names <- names[order(names$id),]
rownames(names) <- names$id

missing.ids <- rownames(names)[is.na(names$name)]
names[missing.ids, 'name'] <- read.csv('gene-names-missing.tsv', row.names='id', sep='\t')[missing.ids, 'name']

filtered.network$target_name <- names[filtered.network$target_id, 'name']
filtered.network$TF_name     <- names[filtered.network$TF_id,     'name']

# any names missing?
which(is.na(filtered.network$target_name) | is.na(filtered.network$TF_name))


write.table(filtered.network, file="network.tsv", sep="\t", quote = FALSE, row.names = FALSE)


