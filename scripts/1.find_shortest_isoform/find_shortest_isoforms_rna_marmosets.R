# Find shortest isoforms for Marmosets files

rm(list = ls(all.names = TRUE)) #will clear all objects includes hidden objects.
# gc() #free up memrory and report the memory usage.

# Loading library
library(seqinr)
library(stringr)
library(stringi)

# Parameter
species <- 'Marmosets'

# Loading input file
file_path <- paste('../../Data/Seq', species, 'GCF_000004665.1_Callithrix_jacchus-3.2_rna.fna', sep = '/')
fa <- read.fasta(file_path, forceDNAtolower = FALSE)

# Filter all records with "X"
filtered_fa <- fa
#filtered_fa <- fa[!grepl('X', names(fa))]

# Construct seqinr object
seq <- getSequence(filtered_fa)
annot <- getAnnot(filtered_fa)

# Manipulate RNA matrix
df <- data.frame(matrix(unlist(annot), nrow=length(filtered_fa), byrow=T),stringsAsFactors=FALSE)
colnames(df) <- 'old_header'
df$old_header <- gsub('>', '', df$old_header)
df$refseq <- sapply(strsplit(df$old_header, " "), function(v) {return(v[1])})

# Add seq information 
df$seq <- sapply(seq, c2s)
df$seq_length <- nchar(df$seq)

df$gene_annot <- gsub("[\\(\\)]", "", regmatches(df$old_header, gregexpr("\\(.*?\\)", df$old_header)))

# if there's multiple gene symbols brom the annotation, choose the last one without : or space 
print(length(df[startsWith(df$gene_annot, "c\""), 'gene_annot'] ))

# debug
a <- data.frame(sapply(strsplit(df[startsWith(df$gene_annot, "c\""), 'gene_annot'], '\\"'), function(x) {return(tail(x[!grepl(':| ', x)], n=1))}))

df[startsWith(df$gene_annot, "c\""), 'gene_annot'] <- sapply(strsplit(df[startsWith(df$gene_annot, "c\""), 'gene_annot'], '\\"'), function(x) {return(tail(x[!grepl(':| ', x)], n=1))})

# if there's multiple gene symbols brom the annotation, choose the last one
#df[startsWith(df$gene_annot, "c"), 'gene_annot'] <- sapply(strsplit(df[startsWith(df$gene_annot, "c"), 'gene_annot'], '\\"'), function(x) {return(tail(x, n=1))})

# XM or NM
df$refseq_prefix <-  sapply(strsplit(df$refseq, "_"), function(v) {return(v[1])})

# Get new headers
df$new_header <- paste(df$refseq, df$gene_annot, sep = '|')

# Sort the df by gene and seq length 
sorted_df <- df[order(df$gene_annot, df$refseq_prefix, df$seq_length), ]
sorted_df <- sorted_df[!duplicated(sorted_df[,c('gene_annot')]),]

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
