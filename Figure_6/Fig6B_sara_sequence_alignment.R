## Multiple sequence alignment and sequence logo of sarA sequences
## Author: Angela Jones
## Date started: April 4, 2023

# load libraries
library(tidyverse)
library(magrittr)
library(seqinr)
library(Biostrings)
library(msa)

# load unalinged sequences in fasta format
sara_seq <- readAAStringSet("Data/sara_filtered_n4554_isolates_translated.fasta")

# conduct multiple sequence alignment
sara_msa <- msa(sara_seq, method = "ClustalW")

# define function that exports msa file as fasta
alignment2Fasta <- function(alignment, filename) {
  sink(filename)
  
  n <- length(rownames(alignment))
  for(i in seq(1, n)) {
    cat(paste0('>', rownames(alignment)[i]))
    cat('\n')
    the.sequence <- toString(unmasked(alignment)[[i]])
    cat(the.sequence)
    cat('\n')  
  }
  
  sink(NULL)
}

# export msa as aligned fasta
alignment2Fasta(test_msa, 'Data/sara_filtered_n4554isolates_msa_out.fasta')

# make sequence motif plot for aligned AAs 190-201
test_seq <- readAAStringSet("Data/sara_filtered_n4554isolates_msa_out.fasta")
pwm <- consensusMatrix(test_seq) # get consensus matrix from aligned fasta with Biostrings
pwm_filtered <- as.matrix(pwm[,190:201]) # subset for AAs of interest
mylabels <- as.character(seq(190,201,1)) # make labels that match the aligned numbers of AAs
ggplot() + 
  geom_logo(pwm_filtered, method="p") +
  theme_logo() +
  scale_x_continuous(breaks=seq(1,12,1),labels=mylabels)
ggsave("sara_yxxq_sequence_logo.png")
ggsave("sara_yxxq_sequence_logo.pdf", width = 7.75, height = 3)
