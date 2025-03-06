## Multiple sequence alignment and sequence logo of GSK3B homologs
## Author: Angela Jones
## Date started: October 16, 2023

# load libraries
library(tidyverse)
library(magrittr)
library(seqinr)
library(Biostrings)
library(msa)
library(ggseqlogo)

# load unalinged sequences in fasta format
gsk3b_seq <- readAAStringSet("Data/Human_GSK3B_orthologues_unaligned.fa")

# conduct multiple sequence alignment
gsk3b_msa <- msa(gsk3b_seq) # default method is ClustalW (Muscle and ClustalOmega also available)

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
alignment2Fasta(gsk3b_msa, 'Data/GSK3B_orthologues_msa_out.fasta')

# make sequence motif plot for aligned AAs 948-959
test_seq <- readAAStringSet("Data/GSK3B_orthologues_msa_out.fasta")
pwm <- consensusMatrix(test_seq) # get consensus matrix from aligned fasta with Biostrings
pwm_filtered <- as.matrix(pwm[,948:959]) # subset for AAs of interest
mylabels <- as.character(seq(948,959,1)) # make labels that match the aligned numbers of AAs
ggplot() + 
  geom_logo(pwm_filtered, method = "p") +
  theme_logo() +
  scale_x_continuous(breaks=seq(1,12,1),labels=mylabels)
ggsave("gsk3b_yxxq_sequence_logo.png")
ggsave("gsk3b_yxxq_sequence_logo.pdf", width = 7.75, height = 3)
