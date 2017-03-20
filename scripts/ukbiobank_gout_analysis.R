# This script is to load the ukbiobank data and process it into smaller Rdata files
# that can then be loaded into the thesis as needed

# libraries ####
library(tidyverse)
library(GenomicRanges)

# helper functions ####

filter_gwas_to_coreExome <- function(){
  #coreExome Markers
  markers <- read.table('/media/xsan/archive/merrimanlab/reference_files/VCF/1000Genomes_vcf_files/Phase3_Sept2014/coreExome_1kgp_markers.txt')
  names(markers) <- c("CHR","BP")
  markersGR <- GRanges(seqnames = paste0('chr',markers$CHR), IRanges(start = markers$BP, width = 1))
}

load_gwas <- function(){
  # takes a directory of plink assoc files that have been converted to .tsv
   temp <- list()
  for( f in list.files("~/Git_repos/UkBioBank/Gout_GWAS",pattern = ".+chr[1-9].+combined.+tsv", full.names = TRUE)){
    temp[[f]] <- data.table::fread(f, header = TRUE, sep = '\t', colClasses = c('integer','integer',rep('chraracter',4),rep('numeric',6), rep('character',3), rep('numeric',3)))
  }
  ukbio_gwas <- bind_rows(temp) %>% arrange(CHR,BP)
}


# all - gwas ####
# initial design is to go through for the all case and then split it into generic functions
# so that each criteria is analysed easily
all_age_sex <- load_gwas()

#find the significant snps
all_sig <- all_age_sex %>% filter(P < 1e-8)

# top 10 SNPs
all_sig %>% arrange(P) %>% head(n=10)

p <- all_age_sex %>%  filter(P < 0.001) %>% ggman::ggmanhattan(.)
ggsave(plot = p, filename = 'images/05_selection_and_association/all_age_sex.png', device ='png')

#find the genes that each signifcant snp is associated with

# hosp - gwas ####

