load("~/Dropbox/db_ihs_nsl_sig_genes.Rdata")

library(dplyr)

library(DBI)
library(GenomicRanges)
library(ggplot2)
con <- dbConnect(RPostgres::Postgres(),dbname = 'selectiondw_test', 
                 host = 'biocmerrimanlab.otago.ac.nz', # i.e. 'ec2-54-83-201-96.compute-1.amazonaws.com'
                 port = 5432, # or any other port specified by your DBA
                 user = 'murray_admin')

dimgene <- dbGetQuery(con, "select * from dimgene inner join dimpos on dimgene.posid = dimpos.posid;") %>% data.frame() %>% dplyr::select(-posid, -posid.1)

markers <- read.table('/media/xsan/archive/merrimanlab/reference_files/VCF/1000Genomes_vcf_files/Phase3_Sept2014/coreExome_1kgp_markers.txt')
markersGR <- GRanges(seqnames = paste0('chr',markers$V1), IRanges(start = markers$V2, width = 1))

dimgeneGR <- GRanges(seqnames = paste0('chr',dimgene$chrom), IRanges(start = dimgene$chrom_start, end= dimgene$chrom_end), gene = dimgene$genename)

ce_dimgene_matches <- as.matrix(findOverlaps(markersGR, dimgeneGR, ignore.strand=TRUE))


geneCounts <- data.frame(table(as.data.frame(dimgeneGR)[ce_dimgene_matches[,2],]$gene))
names(geneCounts)[1] <- 'gene'
nrow(geneCounts)
nrow(bind_rows(ihs_sig_genes))


bind_rows(ihs_sig_genes) %>% inner_join(.,geneCounts, by = 'gene') %>% View()

# count vs proportion
bind_rows(ihs_sig_genes) %>% inner_join(.,geneCounts, by = 'gene') %>% 
  tidyr::gather(.,"pop","count", c(1:37))%>% 
  ggplot(., aes(x = count, y = count/Freq)) + 
  geom_point() + facet_wrap(~pop, scales = "free")

# count versus count weighted on proportion
bind_rows(ihs_sig_genes) %>% inner_join(.,geneCounts, by = 'gene') %>% 
  tidyr::gather(.,"pop","count", c(1:37))%>% 
  ggplot(., aes(x = count, y = count * count/Freq)) + 
  geom_point() + facet_wrap(~pop, scales = "free")


bind_rows(ihs_sig_genes) %>% inner_join(.,geneCounts, by = 'gene') %>% 
  tidyr::gather(.,"pop","count", c(1:37)) %>% mutate(weighted = count*count/Freq) %>% View()


bind_rows(ihs_sig_genes) %>% inner_join(.,geneCounts, by = 'gene') %>% replace(is.na(.),0) %>%
tidyr::gather(.,"pop","count", c(1:37)) %>% mutate(weighted = count*count/Freq) %>% filter(pop %in% c("CIM","NZM","SAM","TON","PAC")) %>% tidyr::spread("pop", count) %>% View()

bind_rows(ihs_sig_genes) %>% inner_join(.,geneCounts, by = 'gene') %>% replace(is.na(.),0) %>%
  tidyr::gather(.,"pop","count", c(1:37)) %>% 
  inner_join(., data.frame(dimgeneGR) %>% dplyr::select(gene,width), "gene") %>%
  mutate(snpsPerbp = Freq/width) %>% 
  mutate( weighted = count/snpsPerbp) %>% 
  filter(pop %in% c("CIM","NZM","SAM","TON","PAC")) %>% 
  tidyr::spread("pop", count) %>% View()



# pathway analysis
library(GSEABase)
library(limma)
library(doParallel)
pathwayList <- GSA::GSA.read.gmt('~/Dropbox/Phase3_analysis/Pathways/PathwayCommons/Pathway Commons.7.All.GSEA.hgnc.gmt')
registerDoParallel(cores = 6)
#larger number = better
ranked_lists <- list()
results<-list()
for(p in unique( colnames(bind_rows(ihs_sig_genes))[!colnames(bind_rows(ihs_sig_genes)) %in% 'gene'] ) ){
  ranked_lists[[p]] <- bind_rows(ihs_sig_genes) %>% inner_join(.,geneCounts, by = 'gene') %>% replace(is.na(.),0) %>%
    tidyr::gather(.,"pop","count", c(1:37)) %>% 
    inner_join(., data.frame(dimgeneGR) %>% dplyr::select(gene,width), "gene") %>%
    mutate(snpsPerbp = Freq/width) %>% 
    mutate( weighted = count/snpsPerbp) %>% 
    filter(pop == p) %>% 
    tidyr::spread("pop", count) %>% arrange(weighted) %>% dplyr::select(gene)
  ranked_lists[[p]]$rank <- as.numeric(row.names(ranked_lists[[p]]))
  results[[p]] <- foreach(path_index=1:length(pathwayList$genesets)) %dopar% {
    limma::geneSetTest(index = ranked_lists[[p]]$gene %in% pathwayList$genesets[[path_index]], statistics = ranked_lists[[p]]$rank, ranks.only = TRUE, alternative = 'mixed')
  }
}
lapply(results,function(x){table(p.adjust(unlist(x), method='bonferroni') < 0.05)})
lapply(results, function(x){which(p.adjust(unlist(x), method='bonferroni') < 0.05)})
limma::barcodeplot(statistics = ranked_lists[["PAC"]]$rank,index = ranked_lists[["PAC"]]$gene %in% pathwayList$genesets[[1863]] )

pathwayList$geneset.names[[1863]]
pathwayList$geneset.names[[3084]]
pathwayList$geneset.names[[3325]]
pathwayList$geneset.names[[4917]]

pathwayList$geneset.descriptions[grep('urate', x = pathwayList$geneset.descriptions)]









