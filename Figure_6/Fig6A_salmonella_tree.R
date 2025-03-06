## Integration and plotting of sarA blast with Salmonella phylogeny
## Author: Angela Jones
## Date started: April 4, 2023

# load libraries
library(treeio)
library(ggtree)
library(readxl)
library(tidyverse)
library(magrittr)

# import tree data
tree <- read.nexus("Data/worley_mbio2018_salm_newick_tree.txt")
tree_table <- tree %>% as.treedata() %>% as_tibble()

# import strain information
strain_df <- read_xlsx("Data/worley_mbio2018_tables1.xlsx")

# import sarA blast results
blast_df <- read_csv("sara_pblast_allisolates.csv")
blast_df %<>% mutate(sara_containing = ifelse(Alignment.length>=100 & E.value <0.001,"yes","no"))

# summarize sarA presence in different serovars from blast results
count_yes <- blast_df %>% group_by(Serovar_new, sara_containing) %>% count() %>% filter(sara_containing=="yes")
count_total <- blast_df %>% group_by(Serovar_new) %>% count() %>% filter(Serovar_new %in% count_yes$Serovar_new)
count_yes$n_total <- count_total$n
count_yes %<>% mutate(prop_yes = n/n_total) #%>% filter(prop_yes >= 0.1)
blast_wSara <- blast_df %>% filter(Serovar_new %in% count_yes$Serovar_new) %>%
  mutate(sara_containing_log10e=-log10(E.value))
blast_wSara_summary <- blast_wSara %>% group_by(Serovar_new) %>%
  summarize(mean_log10e=mean(sara_containing_log10e),
  mean_alignment_length=mean(Alignment.length))
blast_woSara <- blast_df %>% filter(Alignment.length <100 | E.value>+0.001) %>%
  mutate(sara_containing_log10e=NA)
blast_woSara_summary <- blast_woSara %>% group_by(Serovar_new) %>%
  summarize(mean_log10e=mean(sara_containing_log10e),
            mean_alignment_length=mean(Alignment.length)) %>%
  filter(!Serovar_new %in% blast_wSara_summary$Serovar_new)
blast_summary <- bind_rows(blast_wSara_summary, blast_woSara_summary)

# join blast summary with strain data
strain_blast_df <- left_join(strain_df, blast_summary)
names(strain_blast_df)[1] <- "label"
strain_blast_df_clean <- strain_blast_df %>% drop_na(label)

# join tree table with blast info
tree_table_withinfo <- left_join(tree_table,strain_blast_df_clean, by = "label")

# remove duplicate serovars and filter to remove blast NAs
distinct_serovars <- tree_table_withinfo %>% distinct(Serovar_new, .keep_all = TRUE) %>%
  drop_na(mean_alignment_length)
tips_to_trim <- anti_join(tree_table_withinfo, distinct_serovars, by = "label")

# turn back into tree data
tree_data_withinfo <- tree_table_withinfo %>% as.treedata()

# drop tips
tree_table_withinfo_reduced <- drop.tip(tree_data_withinfo, tips_to_trim$label) %>% as.tibble()

# relabel tips with serovar
tree_data_withinfo_reduced <- tree_table_withinfo_reduced %>% as.treedata()

# reroot to arizonae
trda2 <- root(tree_data_withinfo_reduced, outgroup = "FDA00001960", edgelabel = TRUE) %>% as.tibble() %>% rename("isolate"=label,"label"=Serovar_new) %>% as.treedata()

# plot annotated tree
trda2 %>%
  ggtree(layout = "circular") +
  geom_tippoint(mapping = aes(color = mean_log10e), size=2.5) +
  #geom_label(aes(x=branch, label=Serovar_new)) + 
  geom_tiplab(offset = 0.5) + 
  scale_color_gradientn(colors = rev(RColorBrewer::brewer.pal(10,"Spectral"))) +
  theme(legend.position = "right") +
  labs(color="-Log10[E-Value]")
ggsave("pubmlst_tree_alignlength100_all.pdf", width = 10, height=9)
ggsave("pubmlst_tree_alignlength100_all.png", width =10, height=9)
