# This script is to load the ukbiobank data and process it into smaller Rdata files
# that can then be loaded into the thesis as needed

# libraries ####
library(GenomicRanges)
library(VariantAnnotation)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)
library(tidyverse)



# helper functions ####

txdb_gene_annotate <- function(gr){
  # find all unique regions from gr and annotate on gene info
  # take the returned data.frame and merge by c('seqnames','start','end') onto original data
  
  # library(TxDb.Hsapiens.UCSC.hg19.knownGene)
  # library(dplyr)
  # library(org.Hs.eg.db)
  txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
  all_transcripts <- transcripts(txdb) # uses txname as column in GR
  #filter for unique regions in gr
  gr <- gr %>%  data.frame() %>% dplyr::select(seqnames, start, end, strand, marker) %>% distinct() %>% GRanges() 
  # find the matching rows between the 2 GR objects
  overlap_hits <- data.frame(findOverlaps(gr,all_transcripts))
  # pull out the GENEID column for the matching ucsc transcripts
  matching_transcripts <- AnnotationDbi::select(TxDb.Hsapiens.UCSC.hg19.knownGene , 
                                                columns = c("TXID","TXCHROM","TXSTART","TXEND","TXSTRAND", "GENEID"),
                                                keys = as.character(overlap_hits$subjectHits), keytype="TXID")
  # match GENEID to the gene SYMBOL
  gene_list <- AnnotationDbi::select(org.Hs.eg.db, keys =  unique(matching_transcripts$GENEID), columns = c('SYMBOL','ENTREZID'), keytype="ENTREZID")
  # add gene SYMBOL and the row matching info on to the ucsc transcripts using a double merge
  annotated_ucsc_transcripts <- merge(merge(gene_list, matching_transcripts, by.x = "ENTREZID", by.y ='GENEID'), overlap_hits, by.x = "TXID", by.y = "subjectHits")
  # select strand, gene symbol and row matching and remove duplicate rows
  filtered_annotated_ucsc_transcripts <- annotated_ucsc_transcripts %>% dplyr::select(SYMBOL, TXSTRAND, queryHits) %>% filter(!is.na(SYMBOL)) %>% distinct()
  # column to allow merging of original data with the gene info
  gr$id <- 1:NROW(gr)
  # add the gene info onto the original
  merge(gr, filtered_annotated_ucsc_transcripts, by.x = "id", by.y="queryHits", all.x = TRUE) %>% dplyr::select(-id)
}

filter_gwas_to_coreExome <- function(){
  #coreExome Markers
  markers <- read.table('/media/xsan/archive/merrimanlab/reference_files/VCF/1000Genomes_vcf_files/Phase3_Sept2014/coreExome_1kgp_markers.txt')
  names(markers) <- c("CHR","BP")
  markersGR <- GRanges(seqnames = paste0('chr',markers$CHR), IRanges(start = markers$BP, width = 1))
}

load_gwas <- function(files){
  # takes a directory of plink assoc files that have been converted to .tsv
   temp <- list()
  for( f in files){
    temp[[f]] <- data.table::fread(f, header = TRUE, sep = '\t', colClasses = c('integer','integer',rep('chraracter',4),rep('numeric',6), rep('character',3), rep('numeric',3)))
  }
  ukbio_gwas <- bind_rows(temp) %>% arrange(CHR,BP)
}

annotate_snps <- function(gwas_data) {
  target <- with(gwas_data,
                 GRanges( seqnames = Rle(paste0('chr',CHR)),
                          ranges   = IRanges(BP, end= BP + 1),
                          strand   = Rle(strand("*")),
                          marker = SNP
                 ))
  txdb_gene_annotate(target)
}




genome_sig_p <- list()

# all - gwas ####
# initial design is to go through for the all case and then split it into generic functions
# so that each criteria is analysed easily
data <- load_gwas( list.files("~/Git_repos/UkBioBank/Gout_GWAS",pattern = ".+chr[1-9].+combined.+tsv", full.names = TRUE))

#find the significant snps
genome_sig_p[["all_age_sex"]] <-data %>% filter(P < 1e-8)

annotate_snps(genome_sig_p[["all_age_sex"]]) %>% group_by(SYMBOL) %>% tally()
# top 10 SNPs
all_sig %>% arrange(P) %>% head(n=10)

p <- data %>%  filter(P < 0.001) %>% ggman::ggmanhattan(.)
ggsave(plot = p, filename = 'images/05_selection_and_association/all_age_sex.png', device ='png')

# unadjusted all
data <- load_gwas( list.files("/media/xsan/archive/merrimanlab/central_datasets/ukbiobank/gwas/gout/unadjusted/",pattern = ".+chr[1-9].+combined.+tsv", full.names = TRUE))
genome_sig_p[["all_unadjusted"]] <-data %>% filter(P < 1e-8)
tally_and_annotate_snps(genome_sig_p[["all_unadjusted"]])
p <- data %>%  filter(P < 0.001) %>% ggman::ggmanhattan(.)
ggsave(plot = p, filename = 'images/05_selection_and_association/all_unadjusted.png', device ='png')


#find the genes that each signifcant snp is associated with

# hosp - gwas ####
# adjusted age sex
list.files("/media/xsan/scratch/merrimanlab/murray/working_dir/UkBio/GWAS_all_controls/controls/hosp/adjusted/",pattern = ".+age_sex_chr[1-9].+logistic.+tsv", full.names = TRUE)

# unadjusted
list.files("/media/xsan/scratch/merrimanlab/murray/working_dir/UkBio/GWAS_all_controls/controls/hosp/adjusted/",pattern = ".+hosp_chr[1-9].+logistic.+tsv", full.names = TRUE)


