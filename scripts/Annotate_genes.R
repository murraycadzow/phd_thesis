library(dplyr)
library(GenomicRanges)
library(VariantAnnotation)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)

# target <- with(all_age_sex,
#                GRanges( seqnames = Rle(paste0('chr',CHR)),
#                         ranges   = IRanges(BP, end= BP + 1),
#                         strand   = Rle(strand("*")),
#                         marker = SNP
#                ) )
#annotate variants


txdb_gene_annotate <- function(gr){
  # find all unique regions from gr and annotate on gene info
  # take the returned data.frame and merge by c('seqnames','start','end') onto original data
  
  # library(TxDb.Hsapiens.UCSC.hg19.knownGene)
  # library(dplyr)
  # library(org.Hs.eg.db)
  txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
  all_transcripts <- transcripts(txdb) # uses txname as column in GR
  #filter for unique regions in gr
  if("marker" %in% colnames(gr)){
  gr <- gr %>%  data.frame() %>% dplyr::select(seqnames, start, end, strand, marker) %>% distinct() %>% GRanges()
  } else{
    gr <- gr %>%  data.frame() %>% dplyr::select(seqnames, start, end, strand) %>% distinct() %>% GRanges()
  }
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