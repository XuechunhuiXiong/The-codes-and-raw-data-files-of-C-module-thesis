library(readr)

library(dplyr)

# Read in the FPKM data file and the transcription factor TSV file

fpkm_data <- read.csv("~/Project/Data/gene_fpkm.csv", header = TRUE)

tf_data <- read_tsv("~/Project/Data/protein_class_Transcription.tsv")


# Extract the gene names from the FPKM data, skipping "-" values

fpkm_genes <- fpkm_data$gene_short_name[fpkm_data$gene_short_name != "-"]

# Extract the gene names and gene synonyms from the transcription factor data, the code was modified because there were multiple names in the column “gene synonyms”

tf_genes <- c(tf_data$Gene, unlist(strsplit(tf_data$`Gene synonym`, ", ")))


# Check the duplication. And remove any leading/trailing whitespace before checking for duplicates for tf_genes

fpkm_genes_duplicates <- fpkm_genes[duplicated(fpkm_genes)]

tf_genes_withoutwhitespace <- trimws(tf_genes)  
tf_genes_duplicates <- tf_genes_withoutwhitespace[duplicated(tf_genes_withoutwhitespace)]

print(fpkm_genes_duplicates)
print(tf_genes_duplicates)

# Remove duplicates

fpkm_genes_uni <- unique(fpkm_genes) 
tf_genes_uni <- unique(tf_genes) 

# Find the intersection of gene names between FPKM data and transcription factor data (exact match)

common_genes <- intersect(fpkm_genes_uni, tf_genes_uni)


# Check the duplication

common_genes_duplicates <- common_genes[duplicated(common_genes)]

print(common_genes_duplicates)


# Subset the FPKM data based on the common genes. This selects only the rows from the FPKM data where the gene names are present in `common_genes`.

fpkm_subset_match <- fpkm_data[match(common_genes, fpkm_data$gene_short_name), ]


# Save the data frame to a CSV file

write.csv(fpkm_subset_match, "~/Project/Data/fpkm_subset_match.csv", row.names = FALSE)

