panel <- read.delim(paste0('~/data/NZ_coreExome_1kgp/nz_1kgp.panel'), stringsAsFactors = FALSE)

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
ihs_clus_regions <- readRDS('~/data/NZ_coreExome_1kgp/haplotype/ihs_clus_regions-14-7-2017.RDS')
#sig_ihs <- readRDS('~/data/NZ_coreExome_1kgp/haplotype/sig_ihs_clus-14-7-2017.RDS')
sig_ihs_val <- readRDS('~/data/NZ_coreExome_1kgp/haplotype/sig_ihs_values_with_genes_14-11-2017.RDS')



ihs_clus <-readRDS('~/data/NZ_coreExome_1kgp/haplotype/ihs_clus_regions-14-7-2017.RDS')
nsl_clus_regions <- readRDS('~/data/NZ_coreExome_1kgp/haplotype/nsl_clus_regions-14-7-2017.RDS')
sig_nsl_val <- readRDS('~/data/NZ_coreExome_1kgp/haplotype/sig_nsl_values_with_genes_14-11-2017.RDS')
#sig_nsl <- readRDS('~/data/NZ_coreExome_1kgp/haplotype/sig_nsl_clus-14-7-2017.RDS')

nsl_clus <-readRDS('~/data/NZ_coreExome_1kgp/haplotype/nsl_clus_regions-14-7-2017.RDS')

#100kb window 10kb slide pairwise fst
pol_fst <- readRDS('~/data/NZ_coreExome_1kgp/100kbWindow_intra/windowed_poly_chr_fst_popgenome_ns3.14-11-2017.RDS')

# load xpehh
xpehh <- read_delim('~/data/NZ_coreExome_1kgp/haplotype/sig_xpehh_29-10-2017.csv', delim = ',') %>% filter(!pop1 %in% c("POL","EPN","WPN", "NAD"), !pop2 %in% c("POL","EPN","WPN", "NAD"))

# find all of the genes from POL that have something to suggest they may have been selected
table(unique(lower_sig_stats$SYMBOL) %in% (lower_sig_stats %>% filter(statname %in% c('Tajima.D','Fay.Wu.H','Zeng.E', 'Fu.Li.F', 'Fu.Li.D'), pop %in% c('CIM','NZM','SAM','TON')) %>% .[['SYMBOL']] %>% unique()))

load('~/data/gwas_catalog/diseaseGR-25-7-2017.RData')# brings in objects called {gc_urate_gout,urate_gout,kd,metsyn,obesity,t2d}_GR

# made from Thesis/selectionpipeline/consec_regions.R
consec_lower_regions <- readRDS('~/data/NZ_coreExome_1kgp/100kbWindow_intra/filtered/100kbwindows_lower_sig_consec_regions_annotated-17-10-2017.RDS')
consec_upper_regions <- readRDS('~/data/NZ_coreExome_1kgp/100kbWindow_intra/filtered/100kbwindows_upper_sig_consec_regions_annotated-17-10-2017.RDS')
