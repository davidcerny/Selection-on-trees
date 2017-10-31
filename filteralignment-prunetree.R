library(ape)
goloboff <- read.tree("/Users/David/Grive/CRG_Research/Goloboff_molecules_only_shortest.tre")
plants <- read.table("/Users/David/Grive/CRG_Research/RBCL.fin", stringsAsFactors = FALSE)

plants_new <- data.frame(matrix(vector(), ncol=1, nrow=2*nrow(plants)))
for (i in 1:nrow(plants)) {
  plants_new[2*i-1,] <- paste("> ", gsub("____.*", "", plants[i,1]), sep="")
  plants_new[2*i,] <- plants[i,2]
}
write.table(plants_new,
            "/Users/David/Grive/CRG_Research/RBCL_cleanedup.fin",
            col.names = FALSE,
            row.names = FALSE,
            quote = FALSE)

# Produces error:
# 
# Fasta parsing error, RAxML expects an alignment.
# the sequence before taxon > Abatia_parviflora
# : seems to have a different length
# 
# What is the most common sequence length in the sample?
seqs <- plants_new[c(seq(2, nrow(plants_new), by=2)),]
which.max(table(nchar(seqs)))

# How many sequences have lengths other than 1296?
length(seqs[nchar(seqs) != 1296])

# Remove these
rows_to_remove <- 2*which(nchar(seqs) != 1296)
# Remove the corresponding taxon names as well, using the fact that each taxon name
# is located one row above the corresponding sequence
rows_to_remove <- sort(c(rows_to_remove, rows_to_remove - 1))
plants_aligned <- plants_new[-rows_to_remove,]
write.table(plants_aligned,
            "/Users/David/Grive/CRG_Research/RBCL_aligned.fin",
            col.names = FALSE,
            row.names = FALSE,
            quote = FALSE)

# A test run on the resulting functional FASTA alignment showed that 808 sequences were
# exactly identical to other sequences in the dataset. The warning messages identifying
# these sequences were saved into a text file.
# Delete empty lines: sed '/^\s*$/d' raxml_warning.txt > raxml_warnings.txt
# Read in so that spaces are ignored:
warnings <- read.table("/Users/David/Grive/CRG_Research/raxml_warnings.txt", stringsAsFactors = FALSE, sep = "|")

# Grab the names of the 'duplicate' taxa:
duplicates <- vector()
for (i in 1:nrow(warnings)) {
  duplicates <- append(duplicates, paste("> ", gsub("^.*\\sand\\s*|\\s*are\\s.*$", "", warnings[i,1]), sep = ""))
}
duplicates

# Remove their sequences from the dataset:
plants_frame <- as.data.frame(plants_aligned, stringsAsFactors = FALSE)
taxa_to_remove <- pmatch(duplicates, plants_frame[,])
taxa_to_remove <- sort(c(taxa_to_remove, taxa_to_remove + 1))
plants_unique <- plants_frame[-taxa_to_remove,]
# Check:
# nrow(as.data.frame(plants_unique)) == nrow(plants_frame) - length(taxa_to_remove)
write.table(plants_unique,
            "/Users/David/Grive/CRG_Research/RBCL_unique.fin",
            col.names = FALSE,
            row.names = FALSE,
            quote = FALSE)

# Finally, RAxML automatically generated a reduced version of the same dataset excluding
# those columns that only contained undetermined values.

# Remove the tips that are not represented in the alignment. Start with a vector of the
# names of those taxa that are to be retained:
unique_seqs <- as.data.frame(plants_unique, stringsAsFactors = FALSE)[c(seq(1, length(plants_unique), by=2)),]
unique_seq_names <- gsub(".*>\\s", "", unique_seqs)

# Now, make sure that all of these are present in the tree:
match(unique_seq_names, goloboff$tip.label)

# Remove '.' from every occurrence of 'sp.':
corrected_names <- gsub('[.]', "", unique_seq_names)

# Make sure that this took care of all the mismatches:
which(is.na(match(corrected_names, goloboff$tip.label)))

# Change the names in the alignment, too:
plants_unique[c(seq(1,length(plants_unique),by=2))] <- gsub('[.]', "", plants_unique[c(seq(1,length(plants_unique),by=2))])
write.table(plants_unique,
            "/Users/David/Grive/CRG_Research/RBCL_renamed.fin",
            col.names = FALSE,
            row.names = FALSE,
            quote = FALSE)

# See http://blog.phytools.org/2011/03/prune-tree-to-list-of-taxa.html
pruned_tree <- drop.tip(goloboff, setdiff(goloboff$tip.label, corrected_names))
write.tree(pruned_tree, "/Users/David/Grive/CRG_Research/pruned-goloboff-tree.tre")

