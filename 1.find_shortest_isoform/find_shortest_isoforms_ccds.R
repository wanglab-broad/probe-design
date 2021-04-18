# Find shortest isoforms 

# Loading library
library(seqinr)
library(stringr)
library(stringi)

# Parameter
#species <- 'H_sapiens'
species <- 'M_musculus'

# Loading input file
file_path <- paste('../../Data/Seq', species, 'CCDS_nucleotide.current.fna', sep = '/')
fa <- read.fasta(file_path, forceDNAtolower = FALSE)

annot_path <- paste('../../Data/Seq', species, 'CCDS.current.txt', sep = '/')
annot_table <- read.table(annot_path, sep = '\t')


seq = getSequence(fa)
annot = getAnnot(fa)

# Manipulate RNA matrix
df <- data.frame(matrix(unlist(annot), nrow=length(fa), byrow=T),stringsAsFactors=FALSE)
colnames(df) = 'old_header'
df$old_header <- gsub('>', '', df$old_header)
df$ccds_id <- sapply(strsplit(df$old_header, "\\|"), function(v) {return(v[1])})


# Add seq information 
df$seq <- sapply(seq, c2s)
df$seq_length <- nchar(df$seq)

# join two table
colnames(annot_table)[2:5] <- c('refseq_id', 'gene', 'gene_id', 'ccds_id')
df <- merge(df, annot_table, by = 'ccds_id', all.x = T)


print(length(unique(df$gene)))

# Get new headers
df$new_header <- paste(df$ccds_id, df$gene, sep = '|')

# Sort the df by gene and seq length 
sorted_df <- df[order(df$gene, df$seq_length), ]
sorted_df <- sorted_df[!duplicated(sorted_df[,c('gene')]),]

# Creating an object suitable for the write.fasta function
obj <- tapply(sorted_df$seq, 1:length(sorted_df$seq), s2c)

# Writing output file
current_date <- Sys.Date()
current_date <- format(current_date, format="%m_%d_%y")


output_name <- paste(species, '_ccds_shortest_isoforms_', current_date, '.fa', sep = '')
output_path <- paste('../../Data/Seq', species, 'output', output_name, sep = '/')
write.fasta(obj , sorted_df$new_header, file=output_path)
