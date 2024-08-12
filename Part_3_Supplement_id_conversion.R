library(biomaRt)
library(dplyr)
library(tidyr)


# Read the TF filtered 10.csv file
file_path <- "~/Project/Data_New/TF_filtered_10.csv"
tf_filtered_data_new <- read.csv(file_path, stringsAsFactors = FALSE)

# Extract the gene short name column
gene_names_new <- tf_filtered_data_new$gene_short_name

# Connect to the Ensembl database using biomaRt
ensembl <- useEnsembl(biomart = "ensembl", dataset = "ggallus_gene_ensembl")

# Get homologous human gene names
homologs_new <- getBM(attributes = c("external_gene_name", "hsapiens_homolog_associated_gene_name"),
                      filters = "external_gene_name",
                      values = gene_names_new,
                      mart = ensembl)

# Remove the weight and keep a human gene name for each chicken gene name
homologs_unique_new <- homologs_new %>%
  group_by(external_gene_name) %>%
  slice(1) %>%
  ungroup()

# Connect the human genetic database using biomaRt
human_ensembl <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")

# Get the human genetic ID
human_gene_ids_new <- getBM(attributes = c("external_gene_name", "ensembl_gene_id"),
                            filters = "external_gene_name",
                            values = homologs_unique_new$hsapiens_homolog_associated_gene_name,
                            mart = human_ensembl)

# Merge conversion result
final_result_new <- merge(homologs_unique_new, human_gene_ids_new, by.x = "hsapiens_homolog_associated_gene_name", by.y = "external_gene_name")


# Group the final results by human gene names, keeping the first gene ID corresponding to each gene name
final_unique_new <- final_result_new %>%
  group_by(hsapiens_homolog_associated_gene_name) %>%
  slice(1) %>%
  ungroup()

# Save the end result
write.csv(final_unique_new, "~/Project/Data_New/TF_filtered_with_human_ids_new.csv", row.names = FALSE)


merged_data_id_new <- tf_filtered_data_new %>%
  left_join(final_unique_new, by = c("gene_short_name" = "external_gene_name"))

# Save the final merge result
write.csv(merged_data_id_new, "~/Project/Data_New/TF_filtered_with_human_ids_merged_new.csv", row.names = FALSE)
