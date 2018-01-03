panel <- read.delim(paste0('~/data/NZ_coreExome_1kgp/nz_1kgp.panel'), stringsAsFactors = FALSE)
markers <- read.table('~/data/NZ_coreExome_1kgp/nz_1kg_markers.txt', stringsAsFactors = FALSE, header = FALSE)
names(markers) <- c("chrom","chrom_start","marker","ref","alt")
gwas_cat <- read.delim('~/data/gwas_catalog//gwas_catalog_v1.0.1-associations_e89_r2017-06-19.tsv', header=TRUE, stringsAsFactors = FALSE, sep='\t')


# load the significant data
## frequency based
lower_sig_stats <- readRDS('~/data/NZ_coreExome_1kgp/100kbWindow_intra/filtered/100kbwindows_lower_sig_stat_genes_filtered_ns3_resized-31-8-2017.RDS')
lower_sig_stats <- bind_rows(lapply(names(lower_sig_stats), function(y){
  lower_sig_stats[[y]] <- bind_rows(
    lapply(names(lower_sig_stats[[y]]), function(x){
      lower_sig_stats[[y]][[x]] %>% mutate(pop = x, statname =y) 
    }
    ))
})) %>% filter(!pop %in% c("NAD","EPN","WPN","POL"))
lower_sig_stats_binned <-readRDS('~/data/NZ_coreExome_1kgp/100kbWindow_intra/filtered/100kbwindows_lower_sig_stat_filtered_ns3_resized_binned-2-11-2017.RDS') %>% mutate(chrom = paste0('chr',chrom)) %>% select('statname' = stat, everything())



upper_sig_stats <-readRDS('~/data/NZ_coreExome_1kgp/100kbWindow_intra/filtered/100kbwindows_upper_sig_stat_genes_filtered_ns3_resized-31-8-2017.RDS')
upper_sig_stats <- bind_rows(lapply(names(upper_sig_stats), function(y){
  upper_sig_stats[[y]] <- bind_rows(
    lapply(names(upper_sig_stats[[y]]), function(x){
      upper_sig_stats[[y]][[x]] %>% mutate(pop = x, statname =y) 
    }
    ))
})) %>% filter(!pop %in% c("NAD","EPN","WPN","POL"))

upper_sig_stats_binned <-readRDS('~/data/NZ_coreExome_1kgp/100kbWindow_intra/filtered/100kbwindows_upper_sig_stat_filtered_ns3_resized_binned-2-11-2017.RDS') %>% mutate(chrom = paste0('chr',chrom)) %>% select('statname' = stat, everything())

## permutation data
# reorder the lower and upper columns so that the numbers make sense in the tables
fdr_low <- readRDS('~/data/NZ_coreExome_1kgp/100kbWindow_intra/filtered/100kbwindows_chr_stats_popgenome_long_ns3_resized_binned_fdr_low.2-11-2017.RDS') %>% select(stat, pop, lower_per, 'lowerCI' = upper, 'FDR' = med, 'upperCI' = lower, lowest_quant, super_pop)

fdr_high <- readRDS('~/data/NZ_coreExome_1kgp/100kbWindow_intra/filtered/100kbwindows_chr_stats_popgenome_long_ns3_resized_binned_fdr_high.2-11-2017.RDS') %>% select(stat, pop, upper_per, 'lowerCI' = lower, 'FDR' = med, 'upperCI' = upper, highest_quant, super_pop)

perm_median <- readRDS('~/data/NZ_coreExome_1kgp/100kbWindow_intra/filtered/permutations_median_thresholds.2-11-2017.RDS')
perm_ci2.5 <- readRDS('~/data/NZ_coreExome_1kgp/100kbWindow_intra/filtered/permutations_ci2.5_thresholds.2-11-2017.RDS')
perm_ci97.5 <- readRDS('~/data/NZ_coreExome_1kgp/100kbWindow_intra/filtered/permutations_ci97.5_thresholds.2-11-2017.RDS')


## load global summary info
global_summary <- readRDS("~/data/NZ_coreExome_1kgp/100kbWindow_intra/filtered/100kbwindows_summary_popgenome_filtered_3ns_resized.31-8-2017.RDS") %>% filter(!pop %in% c("WPN","EPN","POL","NAD"))
global_super_summary <- readRDS('~/data/NZ_coreExome_1kgp/100kbWindow_intra/filtered/100kbwindows_summary_superpop_popgenome_filtered_3ns_resized.6-11-2017.RDS')

create_sel_summary_table <- function(s){
  global_summary %>% ungroup%>% filter(stat == s) %>% left_join(., panel %>% select(pop, super_pop) %>% distinct(), by = 'pop') %>% arrange(super_pop) %>% select(super_pop, pop, mean, sd, min, lower_1,median,upper_99, max)  %>%  data.frame()
}

## load ihs and nsl data
ihs_clus_regions <- readRDS('~/data/NZ_coreExome_1kgp/haplotype/ihs_clus_regions-15-12-2017.RDS')
ihs_clus <- readRDS('~/data/NZ_coreExome_1kgp/haplotype/sig_ihs_clus-15-12-2017.RDS')
sig_ihs_val <- readRDS('~/data/NZ_coreExome_1kgp/haplotype/sig_ihs_values_with_genes_14-11-2017.RDS') %>% filter(abs(statvalue) > 2.6 )

nsl_clus_regions <- readRDS('~/data/NZ_coreExome_1kgp/haplotype/nsl_clus_regions-15-12-2017.RDS')
nsl_clus <- readRDS('~/data/NZ_coreExome_1kgp/haplotype/sig_nsl_clus-15-12-2017.RDS')
sig_nsl_val <- readRDS('~/data/NZ_coreExome_1kgp/haplotype/sig_nsl_values_with_genes_14-11-2017.RDS')  %>% filter(abs(statvalue) > 2.6 )


#nsl_clus <-readRDS('~/data/NZ_coreExome_1kgp/haplotype/nsl_clus_regions-14-7-2017.RDS')

#100kb window 10kb slide pairwise fst
pol_fst <- readRDS('~/data/NZ_coreExome_1kgp/100kbWindow_intra/windowed_poly_chr_fst_popgenome_ns3.14-11-2017.RDS')

# load xpehh
xpehh <- read_delim('~/data/NZ_coreExome_1kgp/haplotype/sig_xpehh_29-10-2017.csv', delim = ',') %>% filter(!pop1 %in% c("POL","EPN","WPN", "NAD"), !pop2 %in% c("POL","EPN","WPN", "NAD")) %>% filter(abs(xpehh_value) > 2.6)

# find all of the genes from POL that have something to suggest they may have been selected
table(unique(lower_sig_stats$SYMBOL) %in% (lower_sig_stats %>% filter(statname %in% c('Tajima.D','Fay.Wu.H','Zeng.E', 'Fu.Li.F', 'Fu.Li.D'), pop %in% c('CIM','NZM','SAM','TON')) %>% .[['SYMBOL']] %>% unique()))

load('~/data/gwas_catalog/diseaseGR-25-7-2017.RData')# brings in objects called {gc_urate_gout,urate_gout,kd,metsyn,obesity,t2d}_GR
load('~/data/gwas_catalog/zhang_and_extra_18-11-2017.RData') # brings in objects called malaria_GR, zhang_immune_GR, neurological_GR, psychiatric_GR

# made from Thesis/selectionpipeline/consec_regions.R
consec_lower_regions <- readRDS('~/data/NZ_coreExome_1kgp/100kbWindow_intra/filtered/100kbwindows_lower_sig_consec_regions_annotated-15-12-2017.RDS') %>% filter(!pop %in% c('WPN','EPN','NAD'))
# this is incorrect
consec_upper_regions <- readRDS('~/data/NZ_coreExome_1kgp/100kbWindow_intra/filtered/100kbwindows_upper_sig_consec_regions_annotated-15-12-2017.RDS') %>% filter(!pop %in% c('WPN','EPN','NAD'))


plot_gene <- function(p, g){
  wh <- genesymbol[g]
  wh <- range(wh, ignore.strand = TRUE)
  wh_df <- data.frame(wh)
  stats_plot<- bind_rows(sig_ihs_val, sig_nsl_val) %>% filter(pop == p, genename == g)%>% mutate(pop = p, statname = case_when(statid == 20 ~ "nSL", statid == 26 ~ "iHS"), chrom = paste0('chr',chrom)) %>%  select(SYMBOL = genename, chrom, chrom_start, chrom_end, statname) %>%  bind_rows(., lower_sig_stats %>% filter(pop == p, SYMBOL == g) %>% select(contains("chrom"), statname, SYMBOL) ) %>% GenomicRanges::GRanges() %>% ggbio::autoplot(., aes(fill = statname, colour = statname)) + theme(legend.position = 'bottom')
  #ggbio::fixed(stats_plot) <- TRUE
  
  marker_plot <- markers %>% mutate(chrom = paste0('chr', chrom), chrom_end = chrom_start +1) %>% filter(chrom %in% wh_df$seqnames, chrom_start >= wh_df$start -5000, chrom_start <= wh_df$end + 5000) %>% GRanges() %>%  ggbio::autoplot(. ) + theme(legend.position = 'none')
  
  gene_plot <- ggbio::autoplot(Homo.sapiens::Homo.sapiens, which = wh, gap.geom = "chevron")
  #ggbio::fixed(gene_plot) <- TRUE
  hasAxis(gene_plot) <- TRUE
  ggbio::tracks(marker_plot, gene_plot, stats_plot ) + ggtitle(g)
}

#load the pathways tables from panther
panther <- list()
for (file in list.files(path = 'data/Panther/')){
  panther[[file]] <- read.delim(paste0('data/Panther/',file), header = FALSE, stringsAsFactors = FALSE) %>% mutate(file  = file, pop = substring(file, 1, 3))
}

# load the kegg pathways analysis
kegg <- list()
list.files('~/Git_repos/bookdown_thesis/data/enrichr_kegg2016/', include.dirs = FALSE, pattern = '*.txt')
for(file in list.files('~/Git_repos/bookdown_thesis/data/enrichr_kegg2016/', include.dirs = FALSE, pattern = '*.txt')){
  kegg[[file]] <- read.delim(paste0('~/Git_repos/bookdown_thesis/data/enrichr_kegg2016/', file), header = TRUE, stringsAsFactors = FALSE)
  
}

kegg_sig <-list()
kegg_sig <- lapply(names(kegg), function(n){kegg[[n]] %>% filter(Adjusted.P.value < 0.05) %>% arrange(Adjusted.P.value) %>% mutate(filename = n)} )
names(kegg_sig) <- names(kegg)
kegg_sig<- bind_rows(kegg_sig) %>% mutate(stat = lapply(filename, function(fn){strsplit(fn, '_')[[1]][2]})) %>% mutate(stat = substring(stat, 1, nchar(stat)-4), pop = substring(filename, 1,3)) %>% select(-filename)

# for use with fst data where the pop is "POP1_POP2"
# will return the populations so that the pops are aligned to the specified pop2
popswitch <- function(p, pop2){
  p1 <- substring(p,1,3)
  p2 <- substring(p,5,8)
  if(p2 == pop2){
    return(paste0(p1,p2))
  } else{
    return(paste0(p2,p1))
  }
}


# for a given population pull out the genes that had haplotypic significant results AND intersected a significant window from the intra-population statistics
extract_sig <- function(popname){
  bind_rows(sig_ihs_val, sig_nsl_val) %>% filter(pop == popname) %>% select(pop, 'SYMBOL' = genename) %>% inner_join(., bind_rows(., lower_sig_stats %>% filter(pop == popname) %>% select(pop, SYMBOL)), by = c("SYMBOL","pop")) %>% filter(pop == popname, !is.na(SYMBOL)) %>% .[['SYMBOL']] %>% unique()
}

# for a given super population pull out the genes that intersected a significant window from the intra-population statistics
# only uses genes that have both haplotypic AND SFS
extract_sig_super <- function(superpopname){
  bind_rows(sig_ihs_val, sig_nsl_val) %>% 
    select(pop, superpop, 'SYMBOL' = genename) %>% 
    filter(superpop == superpopname) %>% 
    inner_join(., lower_sig_stats %>%
                 left_join(., panel %>% select(pop, "superpop" = super_pop) %>% distinct(), by = 'pop')%>% 
                 select(superpop, SYMBOL, pop) %>% filter(superpop == superpopname), by = c('pop','superpop', "SYMBOL")) %>% 
    filter(!is.na(SYMBOL)) %>% .[['SYMBOL']] %>% unique()
}


#Calculate the genome wide DHS regions using a threshold of at least 60 cell lines having a score > 400
dhs_clusters <- readr::read_delim('~/data/UCSC_DNase1HypersensitivityClusters/hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeRegDnaseClustered/wgEncodeRegDnaseClusteredV3.bed.gz', delim = '\t', col_names = c("chrom","chromStart","chromEnd","name","score","sourceCount","sourceIds","sourceScores")) %>% data.frame()
dhs_regions <- dhs_clusters %>% filter(sourceCount > 60 & score > 400) %>% GRanges()

