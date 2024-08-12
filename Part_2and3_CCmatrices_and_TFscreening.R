library(readr)

library(dplyr)

# TF_New.csv is is obtained by selecting one of the columns in the fpkm_subset_match obtained in the previous step

# Read TF data
TF_New <- read_csv("~/Project/Data_New/TF_New.csv", col_types = cols(
  
  gene_short_name = col_character(),
  
  gene_id = col_character(),
  
  PS_cEpiblast = col_double(),
  NP_Epiblast = col_double(),
  Induced_5hrs = col_double(),
  Induced_9hrs = col_double(),
  Induced_12hrs = col_double(),
  
  FC_0hrs = col_double(),
  FC_5hrs = col_double(),
  FC_9hrs = col_double(),
  FC_12hrs = col_double(),
  FC_PS = col_double(),
  FC_NP = col_double()
))


# Discard rows that contain zero values

TF_filtered <- TF_New %>% filter(FC_0hrs != 0 & FC_PS != 0 & FC_5hrs != 0 & FC_9hrs != 0 & FC_12hrs != 0 & FC_NP != 0)

write.csv(TF_filtered, "~/Project/Data/TF_filtered.csv")


# Select genes that have a value greater than 10 in at least one of the six columns

TF_filtered_10 <- TF_filtered %>%
  rowwise() %>%
  filter(any(c(PS_cEpiblast, NP_Epiblast, Induced_5hrs, Induced_9hrs, Induced_12hrs) > 10)) %>%
  ungroup()


# Save the filtered data

write.csv(TF_filtered_10, "~/Project/Data_New/TF_filtered_10.csv")


# Add two columns from Part_3_Supplement_id_conversion, Extract the gene name and gene id
TF_ids_new <-read_csv("~/Project/Data_New/TF_filtered_with_human_ids_merged_new.csv")
TF_info_new <- TF_ids_new %>% select(gene_short_name, ensembl_gene_id)


# Create a TF data matrix (using both gene names and ids)

data_matrix_10 <- as.matrix(TF_filtered_10[, 3:8])
rownames(data_matrix_10) <- TF_filtered_10$gene_short_name
colnames(data_matrix_10) <- c("FC_0hrs", "FC_PS", "FC_5hrs", "FC_9hrs", "FC_12hrs", "FC_NP")

# Calculate corelation coefficient and p-value

n <- nrow(data_matrix_10)
pcc_matrix_10 <- matrix(NA, nrow = n, ncol = n)
p_value_matrix_10 <- matrix(NA, nrow = n, ncol = n)
rownames(pcc_matrix_10) <- rownames(data_matrix_10)
colnames(pcc_matrix_10) <- rownames(data_matrix_10)
rownames(p_value_matrix_10) <- rownames(data_matrix_10)
colnames(p_value_matrix_10) <- rownames(data_matrix_10)

for (i in 1:(n - 1)) {
  for (j in (i + 1):n) {
    test_result_10 <- cor.test(data_matrix_10[i, ], data_matrix_10[j, ], method = "pearson")
    pcc_matrix_10[i, j] <- test_result_10$estimate
    p_value_matrix_10[i, j] <- test_result_10$p.value
    pcc_matrix_10[j, i] <- test_result_10$estimate
    p_value_matrix_10[j, i] <- test_result_10$p.value
  }
}

pcc_matrix_df_10 <- as.data.frame(pcc_matrix_10)
p_value_matrix_df_10 <- as.data.frame(p_value_matrix_10)
write_csv(pcc_matrix_df_10, "~/Project/Data_New/pcc_matrix_New.csv")
write_csv(p_value_matrix_df_10, "~/Project/Data_New/p_value_matrix_New")

# Set thresholds

pcc_threshold_10 <- 0.7


# Pre-screen eligible gene pairs

valid_pairs_new <- which(pcc_matrix_10 >= pcc_threshold, arr.ind = TRUE)
valid_pairs_new <- valid_pairs_new[valid_pairs_new[, 1] != valid_pairs_new[, 2], ]


# Create a new result data frame for saving all results, and sequence gene pairs for further eliminating duplicates

coexpressed_pairs_10 <- data.frame(
  Gene1 = pmin(rownames(pcc_matrix_10)[valid_pairs_new[, 1]], colnames(pcc_matrix_10)[valid_pairs_new[, 2]]),
  Gene2 = pmax(rownames(pcc_matrix_10)[valid_pairs_new[, 1]], colnames(pcc_matrix_10)[valid_pairs_new[, 2]]),
  CC = pcc_matrix_10[valid_pairs_new],
  p_value = p_value_matrix_10[valid_pairs_new],
  
  
  Gene1_ID = TF_info_new$ensembl_gene_id[match(pmin(rownames(pcc_matrix_10)[valid_pairs_new[, 1]], colnames(pcc_matrix_10)[valid_pairs_new[, 2]]), TF_info_new$gene_short_name)],
  Gene2_ID = TF_info_new$ensembl_gene_id[match(pmax(rownames(pcc_matrix_10)[valid_pairs_new[, 1]], colnames(pcc_matrix_10)[valid_pairs_new[, 2]]), TF_info_new$gene_short_name)]
)


# Remove duplicate gene pairs (A-B or B-A)

coexpressed_pairs_uni_10 <- coexpressed_pairs_10 %>% distinct(Gene1, Gene2, .keep_all = TRUE)

write_csv(coexpressed_pairs_uni_10, "~/Project/Data_new/coexpressed_pairs_uni_new.csv")


# Read Y2H data

Y2H_data <- read_tsv("~/Project/Data/HI-union.tsv", col_names = FALSE)


# Rename Y2H data column

colnames(Y2H_data) <- c("Gene1_ID", "Gene2_ID")


# Add an Interaction column for each pair, indicating interaction

Y2H_pairs <- Y2H_data %>% mutate(Interaction = TRUE)


# Merge PCC results and Y2H interaction information

co_expressed_pairs_with_interaction_10 <- coexpressed_pairs_uni_10 %>%
  left_join(Y2H_pairs, by = c("Gene1_ID", "Gene2_ID")) %>%
  mutate(Interaction = ifelse(is.na(Interaction), FALSE, Interaction))


# Gene pairs in reverse order

Y2H_pairs_reversed <- Y2H_data %>%
  rename(Gene1_ID = Gene2_ID, Gene2_ID = Gene1_ID) %>%
  mutate(Interaction = TRUE)

co_expressed_pairs_withinteraction_10 <- co_expressed_pairs_with_interaction_10 %>%
  left_join(Y2H_pairs_reversed, by = c("Gene1_ID", "Gene2_ID"), suffix = c("", ".reversed")) %>%
  mutate(Interaction = ifelse(Interaction | !is.na(Interaction.reversed), TRUE, FALSE)) %>%
  select(-Interaction.reversed)

write.csv(co_expressed_pairs_withinteraction_10, "~/Project/Data_New/co_expressed_pairs_withinteraction_10.csv")


# Screen all gene pairs with Interaction listed as TRUE

co_expressed_pairs_with_interaction_true_10 <- co_expressed_pairs_with_interaction_10 %>%
  filter(Interaction == TRUE)

write.csv(co_expressed_pairs_with_interaction_true_10, "~/Project/Data_New/co_expressed_pairs_with_interaction_true_10.csv")



