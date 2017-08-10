panel <- read.delim(paste0('~/data/NZ_coreExome_1kgp/nz_1kgp.panel'), stringsAsFactors = FALSE)

# load the significant data
## frequency based
lower_sig_stats <- readRDS('~/data/NZ_coreExome_1kgp/100kbWindow_intra/100kbwindows_lower_sig_stat_genes_filtered_ns3-27-7-2017.RDS')
lower_sig_stats <- bind_rows(lapply(names(lower_sig_stats), function(y){
  lower_sig_stats[[y]] <- bind_rows(
    lapply(names(lower_sig_stats[[y]]), function(x){
      lower_sig_stats[[y]][[x]] %>% mutate(pop = x, statname =y) 
    }
    ))
}))

upper_sig_stats <-readRDS('~/data/NZ_coreExome_1kgp/100kbWindow_intra/100kbwindows_upper_sig_stat_genes_filtered_ns3-27-7-2017.RDS')
upper_sig_stats <- bind_rows(lapply(names(upper_sig_stats), function(y){
  upper_sig_stats[[y]] <- bind_rows(
    lapply(names(upper_sig_stats[[y]]), function(x){
      upper_sig_stats[[y]][[x]] %>% mutate(pop = x, statname =y) 
    }
    ))
}))

## load ihs and nsl data
ihs_clus_regions <- readRDS('~/data/NZ_coreExome_1kgp/haplotype/ihs_clus_regions-14-7-2017.RDS')
nsl_clus_regions <- readRDS('~/data/NZ_coreExome_1kgp/haplotype/nsl_clus_regions-14-7-2017.RDS')

# find all of the genes from POL that have something to suggest they may have been selected
table(unique(lower_sig_stats$SYMBOL) %in% (lower_sig_stats %>% filter(statname %in% c('Tajima.D','Fay.Wu.H','Zeng.E', 'Fu.Li.F', 'Fu.Li.D'), pop %in% c('CIM','NZM','SAM','TON')) %>% .[['SYMBOL']] %>% unique()))

load('~/data/gwas_catalog/diseaseGR-25-7-2017.RData')# brings in objects called {gc_urate_gout,urate_gout,kd,metsyn,obesity,t2d}_GR
