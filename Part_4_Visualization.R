
library(igraph)
library(ggraph)
library(tidyverse)
library(readr)
library(dplyr)

## Visualization
# After screening for pairs with p values less than 0.05:co_expressed_pairs_with_interaction_true_pvalue_ture_10.csv
# Read data
interaction_data_new_pvalue <- read_csv("~/Project/Data_New/co_expressed_pairs_with_interaction_true_pvalue_ture_10.csv")

interactions_pairs_new <- interaction_data_new_pvalue %>% select(Gene1, Gene2)

# Convert the data frame into a igraph object
graph_interaction_new <- graph_from_data_frame(interactions_pairs_new, directed = FALSE)

# Use ggraph for visualization
ggraph(graph_interaction_new, layout = "fr") + 
  geom_edge_link(aes(edge_alpha = 0.8), show.legend = FALSE) +
  geom_node_point(color = "steelblue", size = 5) +
  geom_node_text(aes(label = name), repel = TRUE, size = 3) +
  theme_void()



