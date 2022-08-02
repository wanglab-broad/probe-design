# Find shortest isoforms with whole transcriptome input files

rm(list = ls(all.names = TRUE)) #will clear all objects includes hidden objects.
gc() #free up memrory and report the memory usage.

# Loading library
library(seqinr)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(org.Rn.eg.db)
library(stringr)
library(stringi)

# Parameter
# species <- 'H_sapiens'
# species <- 'M_musculus'
species <- 'R_norvegicus'

# Loading input file
file_name <- "GCF_015227675.2_mRatBN7.2_rna.fna"
file_path <- paste('../../Data/Seq', species, file_name, sep = '/')
fa <- read.fasta(file_path, forceDNAtolower = FALSE)

# Filter all records with "X"
filtered_fa <- fa[!grepl('X', names(fa))]

# Change filtered fa name 
#names(filtered_fa) <-  sapply(strsplit(names(filtered_fa), "\\|"), function(v) {return(v[2])})

# Construct seqinr object
seq <- getSequence(filtered_fa)
annot <- getAnnot(filtered_fa)

# Manipulate RNA matrix
df <- data.frame(matrix(unlist(annot), nrow=length(filtered_fa), byrow=T),stringsAsFactors=FALSE)
colnames(df) <- 'old_header'
df$refseq <- sapply(strsplit(df$old_header, "\\ "), function(v) {return(v[1])})
df$refseq <- gsub(">", "", df$refseq)
# df$refseq <- sapply(strsplit(df$old_header, "\\|"), function(v) {return(v[2])})
df$refseq_quary <- gsub("\\..*", "", df$refseq)

# Add seq information 
df$seq <- sapply(seq, c2s)
df$seq_length <- nchar(df$seq)

# Load org.Hs.eg.db
if (species == 'H_sapiens') {
  db <- org.Hs.eg.db
} else if (species == 'M_musculus') {
  db <- org.Mm.eg.db
} else if (species == 'R_norvegicus') {
  db <- org.Rn.eg.db
}

# Use refseq id find gene symbol
temp <- select(db, 
       keys = df$refseq_quary,
       columns = c("SYMBOL", "REFSEQ"),
       keytype = "REFSEQ")

df$gene_db <- temp$SYMBOL

df$gene_annot <- gsub("[\\(\\)]", "", regmatches(df$old_header, gregexpr("\\(.*?\\)", df$old_header)))

# if there's multiple gene symbols brom the annotation, choose the last one without : or space 
print(length(df[startsWith(df$gene_annot, "c\""), 'gene_annot'] ))
df[startsWith(df$gene_annot, "c\""), 'gene_annot'] <- sapply(strsplit(df[startsWith(df$gene_annot, "c\""), 'gene_annot'], '\\"'), function(x) {return(tail(x[!grepl(':| ', x)], n=1))})

# if there's multiple gene symbols brom the annotation, choose the last one
#df[startsWith(df$gene_annot, "c"), 'gene_annot'] <- sapply(strsplit(df[startsWith(df$gene_annot, "c"), 'gene_annot'], '\\"'), function(x) {return(tail(x, n=1))})

# debug
a <- data.frame(sapply(strsplit(df[startsWith(df$gene_annot, "c\""), 'gene_annot'], '\\"'), function(x) {return(tail(x[!grepl(':| ', x)], n=1))}))

# NR or NM
df$refseq_prefix <-  sapply(strsplit(df$refseq_quary, "_"), function(v) {return(v[1])})

# Compute the final gene
# df$gene_db[is.na(df$gene_db)] <- df$gene_annot[is.na(df$gene_db)]

print(sum(is.na(df$gene_db)))
print(sum(is.na(df$gene_annot)))

# Get new headers
df$new_header <- paste(df$refseq, df$gene_annot, sep = '|')

# Sort the df by gene and seq length 
sorted_df <- df[order(df$gene_annot, df$refseq_prefix, df$seq_length), ]
sorted_df <- sorted_df[!duplicated(sorted_df[, c('gene_annot')]),]

# Extracting longest transcript for each gene
# filtered_seq <- tapply(sequences, genes, function(v) {return(v[which(nchar(v)==min(nchar(v)))])})

# Creating an object suitable for the write.fasta function
obj <- tapply(sorted_df$seq, 1:length(sorted_df$seq), s2c)

# Writing output file
current_date <- Sys.Date()
current_date <- format(current_date, format="%m_%d_%y")
output_name <- paste(species, '_rna_shortest_isoforms_', current_date, '.fa', sep = '')
output_path <- paste('../../Data/Seq', species, 'output', output_name, sep = '/')
write.fasta(obj , sorted_df$new_header, file=output_path)
