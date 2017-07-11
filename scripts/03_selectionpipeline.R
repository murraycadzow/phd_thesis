
lower_sig_stats <- readRDS('~/data/NZ_coreExome_1kgp/30kbWindow_intra/30kbwindows_lower_sig_stat_genes-27-6-2017.RDS')
lower_sig_stats <- bind_rows(lapply(names(lower_sig_stats), function(y){
  lower_sig_stats[[y]] <- bind_rows(
    lapply(names(lower_sig_stats[[y]]), function(x){
      lower_sig_stats[[y]][[x]] %>% mutate(pop = x, statname =y) 
    }
    ))
}))

upper_sig_stats <-readRDS('~/data/NZ_coreExome_1kgp/30kbWindow_intra/30kbwindows_upper_sig_stat_genes-27-6-2017.RDS')

upper_sig_stats <- bind_rows(lapply(names(upper_sig_stats), function(y){
  upper_sig_stats[[y]] <- bind_rows(
    lapply(names(upper_sig_stats[[y]]), function(x){
      upper_sig_stats[[y]][[x]] %>% mutate(pop = x, statname =y) 
    }
    ))
}))


# find all of the genes from POL that have something to suggest they may have been selected
table(unique(lower_sig_stats$SYMBOL) %in% (lower_sig_stats %>% filter(statname %in% c('Tajima.D','Fay.Wu.H','Zeng.E', 'Fu.Li.F', 'Fu.Li.D'), pop %in% c('CIM','NZM','SAM','TON')) %>% .[['SYMBOL']] %>% unique()))
