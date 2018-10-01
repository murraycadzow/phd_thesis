library(tidyverse)
library(ggdendro)
library(scales)
library(ggpubr)

poly_pop <- c('CIM','NZM','TON','SAM')
panel <- read.delim(paste0('~/data/NZ_coreExome_1kgp/nz_1kgp.panel'), stringsAsFactors = FALSE) %>% filter(!pop %in% c('NAD','WPN','POL','EPN'))
###
# made in clustering/coreExome_1kg_popgenome_clustering_100kb_filtered_3ns.Rmd
mat_d_list <- readRDS('~/data/NZ_coreExome_1kgp/100kbWindow_intra/100kbwindows_filtered_3ns_mat_d_list-31-8-2017.RDS')
mat_d_list["fld_neg"] <- NULL
mat_d_list["fld_pos"] <- NULL
mat_d_list["fld_chr"]<-NULL
prop_list <- readRDS('~/data/NZ_coreExome_1kgp/100kbWindow_intra/100kbwindows_filtered_3ns_prop_list-31-8-2017.RDS')
prop_list["fld_neg"] <- NULL
prop_list["fld_pos"] <- NULL
prop_list["fld_chr"]<-NULL
clus_r_list <- readRDS('~/data/NZ_coreExome_1kgp/100kbWindow_intra/100kbwindows_filtered_3ns_clus_r_list-31-8-2017.RDS')
clus_r_list["fld_neg"] <- NULL
clus_r_list["fld_pos"] <- NULL
clus_r_list["fld_chr"]<-NULL
clus_c_list <- readRDS('~/data/NZ_coreExome_1kgp/100kbWindow_intra/100kbwindows_filtered_3ns_clus_c_list-31-8-2017.RDS')
clus_c_list["fld_neg"] <- NULL
clus_c_list["fld_pos"] <- NULL
clus_c_list["fld_chr"]<-NULL
###

## load global summary info
global_summary <- readRDS("~/data/NZ_coreExome_1kgp/100kbWindow_intra/filtered/100kbwindows_summary_popgenome_filtered_3ns_resized.31-8-2017.RDS") %>% filter(!pop %in% c("WPN","EPN","POL","NAD"), stat != 'Fu.Li.D')
global_super_summary <- readRDS('~/data/NZ_coreExome_1kgp/100kbWindow_intra/filtered/100kbwindows_summary_superpop_popgenome_filtered_3ns_resized.6-11-2017.RDS') %>% filter(stat != 'Fu.Li.D')

## load fst
#100kb window 10kb slide pairwise fst
# FST can't be negative. This keeps NAs
pol_fst <- readRDS('~/data/NZ_coreExome_1kgp/100kbWindow_intra/windowed_poly_chr_fst_popgenome_ns3.14-11-2017.RDS') %>% mutate_at(.vars = vars(-contains('chrom')), function(x){ifelse(x > 0, x, 0)})
pol_fst_99 <- pol_fst %>%  gather('pop', "fst", -contains('chrom')) %>% group_by(pop) %>% filter(fst > quantile(fst, 0.99, na.rm = TRUE))
#source('~/Git_repos/bookdown_thesis/scripts/Annotate_genes.R')
#pol_fst_99genes <- txdb_gene_annotate(pol_fst_99 %>% mutate(chrom = paste0('chr',chrom)) %>% dplyr::select(contains("chrom")) %>% GenomicRanges::GRanges())
pol_fst_99genes <- readRDS("~/data/NZ_coreExome_1kgp/100kbWindow_intra/windowed_poly_chr_fst_popgenome_ns3_99per_genes.25-1-2018.RDS") 

marker_loc <- readRDS('~/data/NZ_coreExome_1kgp/snpEff_Annotated/snpeff_terms.RDS') %>% select(chrom = V1, pos = V2, snp = V3, ref = V4, alt = V5, eff )

# genelists
load('~/data/gwas_catalog/diseaseGR-12-3-2018.RData')
genelists <- data.frame(pheno = c(rep("obesity", length(obesity_GR$SYMBOL)), rep("urate", length( gc_urate_gout_GR$SYMBOL)), rep('t2d', length(t2d_GR$SYMBOL)), rep("kd", length(kd_GR$SYMBOL)), rep("metsyn", length(metsyn_GR$SYMBOL))),  genename = c(obesity_GR$SYMBOL, gc_urate_gout_GR$SYMBOL, t2d_GR$SYMBOL, kd_GR$SYMBOL, metsyn_GR$SYMBOL), stringsAsFactors = FALSE) %>% distinct() %>%  mutate(present = 1) %>% spread(pheno, present, 0)

# ihs and nsl sig markers
sig_ihs_val <- readRDS('~/data/NZ_coreExome_1kgp/haplotype/sig_ihs_values_with_genes_14-11-2017.RDS') %>% filter(abs(statvalue) > 2.6 )
sig_nsl_val <- readRDS('~/data/NZ_coreExome_1kgp/haplotype/sig_nsl_values_with_genes_14-11-2017.RDS') %>% filter(abs(statvalue) > 2.6 )


# markers
markers <- read.table('~/data/NZ_coreExome_1kgp/nz_1kg_markers.txt', stringsAsFactors = FALSE, header = FALSE)
names(markers) <- c("chrom","chrom_start","marker","ref","alt")


create_sel_summary_table <- function(s){
  global_summary %>% ungroup%>% filter(stat == s) %>% left_join(., panel %>% select(pop, super_pop) %>% distinct(), by = 'pop') %>% arrange(super_pop) %>% select(super_pop, pop, mean, sd, min, lower_1,median,upper_99, max)  %>%  data.frame()
}

create_sel_summary_table <- function(s){
  global_summary %>% ungroup%>% filter(stat == s) %>% left_join(., panel %>% select(pop, super_pop) %>% distinct(), by = 'pop') %>% arrange(super_pop) %>% select(super_pop, pop, mean, sd, min, lower_1,median,upper_99, max)  %>%  data.frame()
}

gs_td_table <- create_sel_summary_table('Tajima.D')
gs_fwh_table <- create_sel_summary_table('Fay.Wu.H')
#gs_fld_table <- create_sel_summary_table('Fu.Li.D')
gs_flf_table <- create_sel_summary_table('Fu.Lu.F')
gs_ze_table <- create_sel_summary_table('Zeng.E')

td_super_pop_means <- global_super_summary %>% filter(stat == 'Tajima.D') %>% mutate(range = max - min)
fwh_super_pop_means <- global_super_summary %>% filter(stat == 'Fay.Wu.H') %>% mutate(range = max - min)
#fld_super_pop_means <- global_super_summary %>% filter(stat == 'Fu.Li.D') %>% mutate(range = max - min)
flf_super_pop_means <- global_super_summary %>% filter(stat == 'Fu.Li.F') %>% mutate(range = max - min)
ze_super_pop_means <- global_super_summary %>% filter(stat == 'Zeng.E') %>% mutate(range = max - min)



heatmap_col <- scale_fill_gradient(low = "white", high = "steelblue", guide = 'colourbar', breaks = c(0.0, 0.25, 0.5, 0.75, 1.0 ), limits = c(0,1))

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

super_pop_colours <- cbind(panel %>% select(super_pop) %>% distinct() %>% arrange(super_pop) %>% select('group' = super_pop),colour = gg_color_hue(panel %>% select(super_pop) %>% distinct() %>% tally() %>% .[['n']]))

create_windowTable <- function(){
  bind_rows(lapply(names(mat_d_list)[!grepl('chr',names(mat_d_list))], function(x){mat_d_list[[x]] %>% select(-contains('chrom'), -posid) %>% gather(., pop, value, 1:NCOL(.)) %>% group_by(pop) %>% summarise(total_windows =  sum(value)) %>% mutate(d = x)})) %>% mutate(super = sapply(pop, function(x){strsplit(x, '_')[[1]][1]})) %>%  select(-pop) %>%  group_by(d, super) %>% summarise(min = min(total_windows), mean = mean(total_windows), max = max(total_windows)) %>% filter(!d %in% c("fld_neg","fld_pos")) %>% data.frame()
}

create_windowSummaryTable <- function(){
  bind_rows(lapply(names(mat_d_list)[!grepl('chr',names(mat_d_list))], function(x){mat_d_list[[x]] %>% select(-posid, -contains('chrom')) %>% summarise(total_windows = NROW(.), min = min(rowSums(.)), mean = mean(rowSums(.)), max = max(rowSums(.)), median = median(rowSums(.)), sd = sd(rowSums(.)) ) %>% mutate( stat = x)})) %>% filter(!stat %in% c("fld_neg","fld_pos")) %>% select(stat, total_windows, min, median, max, mean, sd) %>% data.frame()
}


theme_heat <- theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), legend.position = 'none', axis.title.y = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank(), axis.title.x = element_blank(), panel.grid = element_blank())


four_plot <- function(mat_d_name,statname){
  hr_c_neg <- clus_c_list[[paste0(mat_d_name,'_neg')]]
  hr_r_neg <- clus_r_list[[paste0(mat_d_name,'_neg')]]
  mat_d_list[[paste0(mat_d_name,'_neg')]] <- mat_d_list[[paste0(mat_d_name,'_neg')]][hr_r_neg$order,]
  
  ddata_x_neg <- dendro_data(hr_c_neg)
  p1 <- ggplot(segment(ddata_x_neg)) +
    geom_segment(aes(x=x, y=y, xend=xend, yend=yend))
  labs_neg <- label(ddata_x_neg)
  labs_neg$group <- unlist(lapply(as.character(labs_neg$label), function(x){strsplit(x, split = '_')[[1]][1]}))
  p1 <- p1 + geom_text(data=label(ddata_x_neg),
                       aes(label=label, x=x, y=y-0.1, colour=labs_neg$group, angle = 90, hjust = 1.1), size = 3)  + theme_dendro() + theme(legend.position = 'none') + coord_cartesian(ylim = c(ddata_x_neg$segments %>% filter(x == xend) %>% summarise(min = min(yend) - (max(y)), max = max(y)) %>% t()))
  
  
  hr_c_pos <- clus_c_list[[paste0(mat_d_name,'_pos')]]
  hr_r_pos <- clus_r_list[[paste0(mat_d_name,'_pos')]]
  mat_d_list[[paste0(mat_d_name,'_pos')]] <- mat_d_list[[paste0(mat_d_name,'_pos')]][hr_r_pos$order,]
  
  ddata_x_pos <- dendro_data(hr_c_pos)
  p2 <- ggplot(segment(ddata_x_pos)) +
    geom_segment(aes(x=x, y=y, xend=xend, yend=yend))
  labs_pos <- label(ddata_x_pos)
  labs_pos$group <- unlist(lapply(as.character(labs_pos$label), function(x){strsplit(x, split = '_')[[1]][1]}))
  p2 <- p2 + geom_text(data=label(ddata_x_pos),
                       aes(label=label, x=x, y=y-0.1, colour=labs_pos$group, angle = 90, hjust = 1.1), size = 3)  + theme_dendro() + theme(legend.position = 'none') + coord_cartesian(ylim = c(ddata_x_pos$segments %>% filter(x == xend) %>% summarise(min = min(yend) - (max(y)), max = max(y)) %>% t()))
  
  lower_dat <- mat_d_list[[paste0(mat_d_name,'_neg')]] %>% select(-contains('chrom')) %>%  mutate(posid = as.numeric(row.names(.))) %>% gather("pop", "value",2:NCOL(.))%>%  mutate(pop = factor(pop)) %>% mutate(pop = factor(pop, levels(pop)[hr_c_neg$order]))

  
  multiplot(
    p1 + ggtitle("A.", paste(statname,"lower tail")),
    
    
    lower_dat %>% ggplot(., aes(x = pop, y = posid)) + geom_tile(aes(fill = value)) + heatmap_col + theme_bw() + theme_heat + ggtitle("C.") + theme(axis.text.x=element_text(colour = as.character(left_join(labs_neg, super_pop_colours, by = 'group' ) %>% .[['colour']]))) + scale_y_continuous(expand = c(0,0)),
    
    p2 + ggtitle("B.",paste(statname,"upper tail")),
    
    mat_d_list[[paste0(mat_d_name,'_pos')]] %>% select(-contains('chrom')) %>%  mutate(posid = as.numeric(row.names(.))) %>% gather("pop", "value",2:NCOL(.))%>%  mutate(pop = factor(pop)) %>% mutate(pop = factor(pop, levels(pop)[hr_c_pos$order])) %>% ggplot(., aes(x = pop, y = posid)) + geom_tile(aes(fill = value)) + heatmap_col + theme_bw() + theme_heat + ggtitle("D.") + theme(axis.text.x=element_text(colour = as.character(left_join(labs_pos, super_pop_colours, by = 'group' ) %>% .[['colour']]))) + scale_y_continuous(expand = c(0,0))
    , cols = 2) 
}

create_sums <- function(mat_d_name){
  mat_d_list[[mat_d_name]] %>% mutate(POLsum =select(., contains("POL")) %>% rowSums(.), EURsum = select(., contains("EUR")) %>% rowSums(.), EASsum = select(., contains("EAS")) %>% rowSums(.), SASsum = select(., contains("SAS")) %>% rowSums(.), AMRsum = select(., contains("AMR")) %>% rowSums(.), AFRsum = select(., contains("AFR")) %>% rowSums(.), ALLsum = select(., -contains('chrom'), -posid) %>% rowSums(.))
}


create_fst_dendro <- function(){
  #data from rnotebooks/CoreExome/coreExome_fst.Rmd
  #return(ggdendrogram(readRDS("~/data/NZ_coreExome_1kgp/whole_chr_fst_popgenome_clust.3-7-2017.RDS"), method = "complete"))
  fst <- readRDS('~/data/NZ_coreExome_1kgp/whole_chr_fst_popgenome.3-7-2017.RDS')
  a <- bind_rows(lapply(fst, function(x){x %>% mutate(pop1 = as.character(pop1), pop2 = as.character(pop2))%>% select(-p1,-p2) %>% filter(!pop1 %in% c('WPN','EPN','NAD'), !pop2 %in% c('WPN','EPN','NAD')) %>% left_join(., panel %>% mutate(super_pop = paste(super_pop, pop, sep = '_')) %>% select("pop1" = pop, super_pop) %>% distinct(), by = "pop1" ) %>% select(-pop1) %>% spread(., "super_pop", fst )}))
  mat_d <- as.matrix(a[,c(-1,-2)])
  rownames(mat_d) <- paste(a$pop2, a$chrom, sep = '_')
  #hr_d <- hclust(dist(mat_d, method="euclidean"), method = "complete")
  hc_d <- hclust(dist(t(mat_d), method="euclidean"), method = "complete")
  
  
  ddata_x <- dendro_data(hc_d)
  p1 <- ggplot(segment(ddata_x)) +
    geom_segment(aes(x=x, y=y, xend=xend, yend=yend))
  labs <- label(ddata_x)
  labs$group <- unlist(lapply(as.character(labs$label), function(x){strsplit(x, split = '_')[[1]][1]}))
  p1 <- p1 + geom_text(data=label(ddata_x),
                       aes(label=label, x=x, y=y-0.1, colour=labs$group, angle = 90, hjust = 1), size = 3)  + theme_dendro() + theme(legend.position = 'none') + coord_cartesian(ylim = c(ddata_x$segments %>% filter(x == xend) %>% summarise(min = min(yend) - (max(y)), max = max(y)) %>% t()))
  return(p1)
}




ihs_dendro_plot <- function(){
  # data from rnotebooks/CoreExome/coreExome_1kg_ihs_nsl_clustering.Rmd
  mat_ihs <- readRDS('~/data/NZ_coreExome_1kgp/haplotype/ihs_clustered_dist_mat-15-12-2017.RDS')
  clus_c <- hclust(dist(t(mat_ihs), method="euclidean"), method = "complete")
  clus_r <- hclust(dist(mat_ihs, method="euclidean"), method = "complete")
  mat_ihs <- mat_ihs[clus_r$order,clus_c$order]
  
  ddata_x <- dendro_data(clus_c)
  ddata_y <- dendro_data(clus_r)
  
  labs_x <- label(ddata_x)
  labs_y <- label(ddata_y)
  labs_x$group <- unlist(lapply(as.character(labs_x$label), function(x){strsplit(x, split = '_')[[1]][1]}))
  labs_y$group <- unlist(lapply(as.character(labs_y$label), function(x){strsplit(x, split = '_')[[1]][1]}))

  p1 <- ggplot(segment(ddata_x)) +
    geom_segment(aes(x=x, y=y, xend=xend, yend=yend))
  
  p1 <- p1 + geom_text(data=label(ddata_x),
                       aes(label=label, x=x, y=y-0.1, colour=labs_x$group, angle = 90, hjust = 1), size = 3)  + theme_dendro() + theme(legend.position = 'none', plot.margin = unit(c(0.5,0.5,0.1,0.5), "cm"), aspect.ratio = 1/2) + coord_cartesian(ylim = c(ddata_x$segments %>% filter(x == xend) %>% summarise(min = min(yend) - (max(y)), max = max(y)) %>% t()))
  
  p2 <- as.data.frame(mat_ihs) %>% mutate(pop1 = rownames(.)) %>% gather("pop2", "value",1:(NCOL(.)-1)) %>% mutate(pop1 = factor(pop1), pop2 = factor(pop2)) %>% mutate(pop1 = factor(pop1, levels = levels(pop1)[clus_r$order]), pop2 = factor(pop2, levels = levels(pop2)[clus_c$order])) %>%  ggplot(., aes(x = pop2, y = pop1, fill = value)) + geom_tile() + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, colour = as.character(left_join(labs_x, super_pop_colours, by = 'group' ) %>% .[['colour']])), 
          axis.text.y = element_text(hjust = 1, vjust = 0.5, colour = as.character(left_join(labs_y, super_pop_colours, by = 'group' ) %>% .[['colour']])), 
          axis.title.y = element_blank(), 
          axis.title.x = element_blank(),
          axis.ticks = element_blank(),
          legend.position = 'bottom', 
          panel.grid = element_blank(), 
          plot.margin = unit(c(0.1,0.5,0.1,0.5), 'cm'), 
          aspect.ratio = 1, 
          panel.border = element_rect(colour = "black", fill=NA, size= 1)) + 
    heatmap_col + labs(fill = "Proportion") + guides(fill = guide_colorbar(barwidth = 6, barheight = 0.5, label.hjust = 0.1))

  ggarrange(p1, p2, heights = c(1, 2),
            ncol = 1, nrow = 2, align = "v", labels = c("A","B"))
}


nsl_dendro_plot <- function(){
  # data from rnotebooks/CoreExome/coreExome_1kg_ihs_nsl_clustering.Rmd
  mat_nsl <- readRDS('~/data/NZ_coreExome_1kgp/haplotype/nsl_clustered_dist_mat-15-12-2017.RDS')
  clus_c <- hclust(dist(t(mat_nsl), method="euclidean"), method = "complete")
  clus_r <- hclust(dist(mat_nsl, method="euclidean"), method = "complete")
  mat_nsl <- mat_nsl[clus_r$order,clus_c$order]
  
  ddata_x <- dendro_data(clus_c)
  ddata_y <- dendro_data(clus_r)
  
  labs_x <- label(ddata_x)
  labs_y <- label(ddata_y)
  labs_x$group <- unlist(lapply(as.character(labs_x$label), function(x){strsplit(x, split = '_')[[1]][1]}))
  labs_y$group <- unlist(lapply(as.character(labs_y$label), function(x){strsplit(x, split = '_')[[1]][1]}))
  
  p1 <- ggplot(segment(ddata_x)) +
    geom_segment(aes(x=x, y=y, xend=xend, yend=yend))
  
  p1 <- p1 + geom_text(data=label(ddata_x),
                       aes(label=label, x=x, y=y-0.1, colour=labs_x$group, angle = 90, hjust = 1), size = 3)  + theme_dendro() + theme(legend.position = 'none', plot.margin = unit(c(0.5,0.5,0.1,0.5), "cm"), aspect.ratio = 1/2) + coord_cartesian(ylim = c(ddata_x$segments %>% filter(x == xend) %>% summarise(min = min(yend) - (max(y)), max = max(y)) %>% t()))
  
  p2 <- as.data.frame(mat_nsl) %>% mutate(pop1 = rownames(.)) %>% gather("pop2", "value",1:(NCOL(.)-1)) %>% mutate(pop1 = factor(pop1), pop2 = factor(pop2)) %>% mutate(pop1 = factor(pop1, levels = levels(pop1)[clus_r$order]), pop2 = factor(pop2, levels = levels(pop2)[clus_c$order])) %>%  ggplot(., aes(x = pop2, y = pop1, fill = value)) + geom_tile() + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, colour = as.character(left_join(labs_x, super_pop_colours, by = 'group' ) %>% .[['colour']])), 
          axis.text.y = element_text(hjust = 1, vjust = 0.5, colour = as.character(left_join(labs_y, super_pop_colours, by = 'group' ) %>% .[['colour']])), 
          axis.title.y = element_blank(), 
          axis.title.x = element_blank(), 
          axis.ticks = element_blank(),
          legend.position = 'bottom', 
          panel.grid = element_blank(), 
          aspect.ratio = 1, 
          panel.border = element_rect(colour = "black", fill=NA, size= 1)
          ) + 
    heatmap_col + 
    labs(fill = "Proportion")  + guides(fill = guide_colorbar(barwidth = 6, barheight = 0.5, draw.llim = TRUE, draw.ulim = TRUE , label.hjust = 0.1))
  
  ggarrange(p1, p2, heights = c(1, 2),
            ncol = 1, nrow = 2, align = "v", labels = c("A","B"))
}




prop_unique <- function(statname){
  create_sums(statname) %>% 
    filter(ALLsum == 1) %>% 
    select(-posid, -contains('sum')) %>% 
    gather(., "pop", "value", -contains('chrom')) %>% 
    group_by(pop) %>% summarise(unique_windows = sum(value))  %>%  
    left_join(., create_sums(statname) %>% 
                select(-posid, -contains('sum')) %>% 
                gather(., "pop", "value", -contains('chrom')) %>% 
                group_by(pop) %>% 
                summarise(total_windows = sum(value)) , by = "pop") %>%  
    mutate(prop = unique_windows / total_windows)
}

# made from clustering/popgenome_freq_stats.Rmd
lower_sig_stats <- readRDS('~/data/NZ_coreExome_1kgp/100kbWindow_intra/filtered/100kbwindows_lower_sig_stat_genes_filtered_ns3_resized-31-8-2017.RDS')
lower_sig_stats <- bind_rows(lapply(names(lower_sig_stats), function(y){
  lower_sig_stats[[y]] <- bind_rows(
    lapply(names(lower_sig_stats[[y]]), function(x){
      lower_sig_stats[[y]][[x]] %>% mutate(pop = x, statname =y) 
    }
    ))
})) %>% filter(!pop %in% c("NAD","EPN","WPN","POL"), statname != 'Fu.Li.D')

# made from clustering/popgenome_freq_stats.Rmd
upper_sig_stats <-readRDS('~/data/NZ_coreExome_1kgp/100kbWindow_intra/filtered/100kbwindows_upper_sig_stat_genes_filtered_ns3_resized-31-8-2017.RDS')

upper_sig_stats <- bind_rows(lapply(names(upper_sig_stats), function(y){
  upper_sig_stats[[y]] <- bind_rows(
    lapply(names(upper_sig_stats[[y]]), function(x){
      upper_sig_stats[[y]][[x]] %>% mutate(pop = x, statname =y) 
    }
    ))
})) %>% filter(!pop %in% c("NAD","EPN","WPN","POL"), statname != "Fu.Li.D")


gene_dendros <- function(gene_data){
  clus_c <- hclust(dist(t(gene_data), method="euclidean"), method = "complete")
  
  ddata_x <- dendro_data(clus_c)
  p1 <- ggplot(segment(ddata_x)) +
    geom_segment(aes(x=x, y=y, xend=xend, yend=yend))
  labs <- label(ddata_x)
  labs$group <- unlist(lapply(as.character(labs$label), function(x){strsplit(x, split = '_')[[1]][1]}))
  p1 <- p1 + geom_text(data=label(ddata_x),
                       aes(label=label, x=x, y=y, colour=labs$group, angle = 90, vjust = 0.5, hjust = 1.1), size = 3)  + theme_dendro() + theme(legend.position = 'none') + coord_cartesian(ylim = c(ddata_x$segments %>% filter(x == xend) %>% summarise(min = min(yend) - (max(y)), max = max(y)) %>% t()))
  return(p1)
}

# genes that are pretty much only in Poly
pol_tail_filter <- function(mat_d_name, pol_min, all_max){
create_sums(mat_d_name) %>% filter(POLsum >= pol_min & ALLsum <= all_max)%>% mutate(chrom = paste0('chr', chrom)) %>%select(1:3) %>% GenomicRanges::GRanges() %>% GenomicRanges::reduce() %>% txdb_gene_annotate() %>% as_tibble() %>% filter(!is.na(SYMBOL))
}


pickrell2009 <- list(
  biaka = 'COX5B,ACTR1B,ZAP70,KIAA1641,MAN2A1,HIST1H4G,HIST1H3F,HIST1H2BH,HIST1H3G,HIST1H2BI,HIST1H4H,BTN3A2,BTN2A2,BTN3A1,BTN2A3,BTN3A3,BTN2A1,BTN1A1,HMGN4,RBMS3,OR10Q1,OR10W1,OR5B17,OR5B3,OR5B2,OR5B12,OR5B21,CPNE4,TBX18,DCAMKL1,SOHLH2,AFF3,SLC28A3,PRTG,NEDD4,SUHW4,TCF12,TBC1D15,TPH2,TCF12,RTN4,FLJ42562,SGEF,KIAA0350,ITGA6,PDK1,CREB5,NSUN6,ARL5B,NEK10,SLC4A7,SCRN1,FKBP14,PLEKHA8,WIPF3,C7orf10,CSNK1G3,SETBP1,PHACTR3,THADA,PLEKHH2,CTNNA3,FLJ37357,SOX6,GDPD4,PAK1,MYO7A,RIOK3,C18orf8,NPC1,ANKRD29,C18orf45,C3orf48,EFHB,RAB5A,PCAF,SGOL1,PFKP,PITRM1,TRIO,EML4,COX7A2L,DNER,PID1,SPOCK1,EPAS1,ATP6V1E2,RHOQ,LOC388965,INPP4A,MGC26733,CNGA3,OPRM1,DKFZP586P0123,PPME1,P4HA3,PGM2L1,KCNE3,B2M,TRIM69,C15orf43,NRXN1,SOD2,WTAP,ACAT2,TCP1,MRPL18,PNLDC1,MAS1,GNA12,CARD11,TTC23,LRRC28,DMN,MAP2,PCDHGA1,PCDHGA2,PCDHGA3,PCDHGB1,PCDHGA4,PCDHGB2,PCDHGA5,PCDHGB3,PCDHGA6,PCDHGA7,PCDHGB4,PCDHGA8,PCDHGB5,PCDHGA9,PCDHGB6,PDCD7,CLPX,CILP,PARP16,PUNC,LRRN1,FER,EPB41L2,TGFBI,SMAD5,TRPC7,FUT8,C14orf145,TSHR,TSCOT,ZFP37,POMC,DNMT3A,SYN3,TIMP3,C8orf34,NPHP4,KCNAB2,DSCAM,EPB41L4B,C9orf4,C9orf5,CTNNAL1,SCLT1,C4orf33,CSS3,COL21A1,GMDS,TAAR2,TAAR5,TAAR1,VNN1,VNN3,VNN2,C6orf192,RPS12' %>% strsplit(., split = ',') %>% unlist() %>% unique(),
  bantu = 'RYR1,SIPA1L3,DPF1,PPP1R14A,SPINT2,LOC541469,C19orf33,YIF1B,KCNK6,C19orf15,PSMD8,GGN,SPRED3,RASGRP4,FAM98C,FDX1,ARHGAP20,PLXDC2,SBF2,ANP32B,C9orf156,FOXE1,HEMGN,HS3ST2,USP31,KIAA1627,SEC24D,SYNPO2,NDST3,PRSS12,MAN2A1,SUMF1,ITPR1,PPP1R12B,SYT2,JARID1B,ETFB,CLDND2,NKG7,LIM2,SIGLEC12,SIGLEC6,ZNF175,SIGLEC5,SIGLEC10,SIGLEC8,LOC729767,KCNJ3,PREX1,ARFGEF2,SYT1,MCM5,RASD2,MB,LOC284912,APOL6,APOL5,SIM1,SUSD1,ROD1,HSDL2,GAN,CMIP,C12orf26,TMTC2,TRPV1,CTNS,TRPV3,CARKL,TAX1BP3,TMEM93,P2RX5,ITGAE,GSG2,ITGB8,LOC440131,PTHLH,FSTL5,NCAM1,RNF32,LMBR1,NOM1,LINGO2,SETBP1,NUBPL,MAMDC4,LCN10,LCN6,LCN8,UNQ2541,TMEM141,C9orf86,PHPT1,EDF1,TRAF2,FBXW5,C8G,PTGDS,FLJ45224,C9orf142,CLIC3,ABCA2,C9orf139,FUT7,KIAA1984,LCN1,OR52L1,OR56A1,OR56B4,OR52B2,OR52W1,C11orf42,C11orf56,CNGA4,CCKBR,OR56A4,CSMD3,ANKRD19,BICD2,ZNF484,FGD3,SAMD3,KIAA1913,COL8A1,FLJ44076,C3orf26,FILIP1L,RNF165,LOXHD1,RIMS1,EPB41L4B,C9orf4,C9orf5,CTNNAL1,PCDH17,FLJ40296,TSCOT,ZFP37,CCDC66,C3orf63,ARHGEF3,CCDC102B,DLG2,APOL3,APOL4,APOL2,APOL1,MYH9,RAB38,FYCO1,CXCR6,XCR1,CCR1,PPFIA2,HS3ST4,FHIT,HIST1H4G,HIST1H3F,HIST1H2BH,HIST1H3G,HIST1H2BI,HIST1H4H,BTN3A2,BTN2A2,BTN3A1,BTN2A3,BTN3A3,BTN2A1,BTN1A1,HMGN4,AXIIR,ZNF131,MGC42105,CCT8,C21orf7,BACH1,MAS1L,OR2H2,GABBR1,MOG,HLA−F,INPP4A,MGC26733,CNGA3,MGC40405,CDK6,DKFZP564O0523,PEX1,RFXDC2,TEX9,SEMA3E,AGK,FLJ40852,SSBP1,TAS2R3,TAS2R4,TAS2R5,LOC136242,KIAA1147,RGS6,TNN,KIAA0040,TNR,SCLT1,C4orf33,PKN2,PDE9A,WDR4,NDUFV3,WDR40C,C10orf59,LIPL1,LIPF,CPA6,DEPDC2,LMBRD1,COL19A1,TCERG1,GPR151,PPP2R2B,CDC27,MYL4,ITGB3,C17orf57,VPS13C,C21orf29,KRTAP10−1,KRTAP10−2,KRTAP10−3,KRTAP10−4,KRTAP10−5,KRTAP10−6,KRTAP10−7,KRTAP10−8,KRTAP10−9,KRTAP10−10,KRTAP10−11,KRTAP12−4,KR,JMJD1A,VPS24,RNF103,RMND5A,CTNNA3,ISL2,ZNF291,ETFA' %>% strsplit(., split = ',') %>% unlist() %>% unique(),
  europe = 'CMYA3,RGS9,SLC44A1,TMEM30B,PRKCH,HIF1A,DCP2,MCC,RFX3,FLJ35880,NGLY1,OXSM,PDE11A,AGPS,TTC30A,TTC30B,PCMTD1,C8ORFK32,COL22A1,KCND2,TSPAN12,HLA−A,HCG9,ZNRD1,RNF39,TRIM31,TRIM40,TRIM15,TRIM26,FLJ45422,TRIM39,RPP21,PPP1R11,LAMA3,C18orf17,GTF3C3,UNQ6411,PGAP1,ANKRD44,TCERG1,GPR151,PPP2R2B,COL22A1,CDH12,ACVR1C,ACVR1,UPP2,GLCE,PAQR5,PREP,TYRP1,C9orf150,NRG3,HMGCR,COL4A3BP,POLK,MPHOSPH6,SNX25,LRP2BP,ANKRD37,C4orf20,CCDC110,RHOJ,GPHB5,PPP2R5E,DENND1A,FGF1,ARHGAP26,PDZRN4,ELMO1,EDNRA,AK3L1,AK3L2,JAK1,PTPRM,ATF6,OLFML2B,NOS1AP,C1orf111,SH2D1B,GRIK1,CCDC53,NUP37,PMCH,C12orf48,COL6A3,PRLH,RAB17,MLPH,NRG3,FAM109A,SH2B3,ATXN2,CUTL2,SMYD3,NRSN1,DCDC2,TMEM66,LEPROTL1,DCTN6,NTRK2,FERD3L,TMEM100,PCTP,MAF,AUH,NFIL3,SLC8A1,AGBL1,FBXO2,FBXO44,FBXO6,EPB41L4B,C9orf4,C9orf5,CTNNAL1,APBA2,NELL1,CDH10,SYNE2,KLHL15,EIF2S3,ZFX,IER3,PRR3,HLA−E,GNL1,ABCF1,PPP1R10,MRPS18B,C6orf134,C6orf136,DHX16,NRM,MDC1,TUBB,FLOT1,DKFZP686A01247,PHOX2B,TMEM33,RYBP,COL8A1,TUBA3C,USP15,MON2,WDR16,STX8,USP43,LCN2,CIZ1,DNM1,GOLGA2,TRUB2,COQ4,SLC27A4,URM1,CEECAM1,C9orf16,PRB4,PRB1,PRB2,NCOR2,FAM101A,C1GALT1,COL28A1,WWOX,MAPK8IP3,IFT140,CRAMP1L,HN1L,NME3,MRPS34,EME2,SPSB3,NUBP2,IGFALS,FAHD1,C16orf73,HAGH,DRG1,EIF4ENIF1,SFI1,PISD,C22orf30,GRIK2,RPS6KA2,RNASET2,FGFR1OP,COMMD10,SEMA6A,ZPLD1,TIMP2,LGALS3BP,CANT1,C1QTNF1,FLJ21865,DNAH7,STK17B,HECW2' %>% strsplit(., split = ',') %>% unlist() %>% unique(),
  mideast = 'SLC44A1,ERBB4,KCNH5,RHOJ,GPHB5,PPP2R5E,CDH18,C8ORFK32,COL22A1,RGS9,APBA2,KCND2,TSPAN12,NGLY1,OXSM,CCDC102B,NRG3,C20orf19,XRN2,NKX2−2,PDZRN4,CRB2,DENND1A,PRKCH,HIF1A,BACH1,GRIK1,FLJ20699,PPARA,LOC150383,PKDREJ,GTSE1,TRMU,CELSR1,SMYD3,RAVER2,JAK1,COL6A3,PRLH,RAB17,MLPH,TRPS1,COL22A1,TMEM132B,TCERG1,GPR151,PPP2R2B,PRDM10,C11orf37,APLP2,ST14,ZBTB44,SCML4,FLJ10159,ACVR1C,ACVR1,UPP2,FGF1,ARHGAP26,COMMD3,FLJ35880,TSLP,TMPRSS7,PLCXD2,PHLDB2,ABHD10,TAGLN3,CYP26B1,BBS9,GLCE,PAQR5,KIF23,RPLP1,LOC729416,MAPK8IP3,IFT140,CRAMP1L,HN1L,NME3,MRPS34,EME2,SPSB3,NUBP2,IGFALS,FAHD1,C16orf73,HAGH,WWOX,ANKS1B,FLJ46363,CLEC12A,CLEC1B,CLEC12B,CLEC9A,CLEC1A,CLEC7A,GABARAPL1,KLRD1,KLRK1,C12orf59,OLR1,GPR64,MLC1,BRD1,ZBED4,ALG12,CRELD2,PIM3,FLJ41993,LOC164714,FOXA2,ERN1,TEX2,PECAM1,PRICKLE2,SYNE2,KITLG,PARD6B,BCAS4,ADNP,DPM1,MOCS3,PCMTD1,CBR4,SH3RF1,NEK1,DLG2,C14orf124,TGM1,CHMP4A,MDP−1,NEDD8,GMPR2,TINF2,RABGGTA,DHRS1,C14orf21,LTB4R2,LTB4R,ADCY4,NFATC4,KIAA0323,CMA1,RIPK3,CIDEB,CBLN3,MYT1L,TMEM100,PCTP,SMYD2,PTPN14,NRG3,OR2T8,OR2L13,OR2L8,OR2AK2,OR2L2,OR2L3,OR2M5,OR2M2,OR2M3,CDH12,EDNRA,POT1,GPR37,LOC401398,LYPLA1,MRPL15,DOK5,MYOD1,SERGEF,KCNC1,NEB,ARL5A,CACNB4,CDH10,HLA−A,HCG9,ZNRD1,RNF39,TRIM31,TRIM40,TRIM15,TRIM26,FLJ45422,TRIM39,RPP21,PPP1R11,FAM109A,SH2B3,ATXN2,CUTL2,ST6GALNAC1,MXRA7,SFRS2,MFSD11,MGAT5B,PTDSR,LOC124512,STK32B,COMMD10,SEMA6A,CBLN4,CAND1,FAM83B,HTR1F,CGGBP1,GRM8,PHF17,SCLT1,C4orf33,DCP2,MCC,SMAD6,PTPRM,HMGCR,COL4A3BP,POLK,CHRNA6,THAP1,RNF170,CHRNB3,SLC24A5,MYEF2' %>% strsplit(., split = ',') %>% unlist() %>% unique(),
  sasia = 'SLC24A5,MYEF2,SLC12A1,DUT,TMEM100,PCTP,KCND2,TSPAN12,KCNH5,RHOJ,GPHB5,PPP2R5E,C8orf15,MTMR9,AMAC1L2,XKR6,FLJ20699,PPARA,LOC150383,PKDREJ,GTSE1,TRMU,CELSR1,NRG3,NGLY1,OXSM,C6orf148,KCNQ5,C8ORFK32,COL22A1,CRB2,DENND1A,PRR3,RNF39,TRIM31,TRIM40,TRIM15,TRIM26,FLJ45422,TRIM39,RPP21,HLA−E,GNL1,ABCF1,FBXW2,PSMD5,PHF19,TRAF1,C5,CEP110,RAB14,FLJ35880,CXCR7,LGALS8,HEATR1,ACTN2,MTR,SLC44A1,N6AMT1,ZNF294,C20orf133,ACVR1C,ACVR1,UPP2,GPR55,CAB39,ITM2C,PSMD1,LCN2,CIZ1,DNM1,GOLGA2,TRUB2,COQ4,SLC27A4,URM1,CEECAM1,C9orf16,LIN28B,BVES,POPDC3,PREP,APBA2,TCERG1,GPR151,PPP2R2B,RGS9,SLC25A36,SPSB4,ACPL2,DKFZP686A01247,PHOX2B,TMEM33,CCDC102B,DCP2,MCC,GPBAR1,PNKD,ARPC2,TMBIM1,MGC50811,SLC11A1,CTDSP1,VIL1,USP37,AAMP,TMTC2,NUDT18,DOK2,XPO7,NPM2,FGF17,EPB49,RAI16,HR,USP25,C21orf34,KIAA1324L,DMTF1,RARB,TOP2B,ACSL3,KCNE4,ANKS1B,PLEKHG3,SPTB,CHURC1,GPX2,RAB15,FNTB,MAX,PDE11A,AGPS,TTC30A,TTC30B,NTRK2,PCMTD1,EFCBP1,TMEM55A,OTUD6B,PRDM10,C11orf37,APLP2,ST14,ZBTB44,ROR2,SPTLC1,PARD6B,BCAS4,ADNP,DPM1,MOCS3,RP1L1,C8orf74,SOX7,PINX1,XKR6,CHRNA6,THAP1,RNF170,CHRNB3,IPMK,ZCD1,UBE2D1,TFAM,CCT8,C21orf7,BACH1,GRIK1,ELMO1,HTR3A,ZBTB16,ISX,C12orf64,GALNT17,CBR4,SH3RF1,SYNE2,SMAD6,NDUFB6,TAF1L,LOC401498,DNAH7,STK17B,HECW2,PLCB4,C20orf103,PAK7,PLCZ1,CAPZA3,COL22A1,FOXA2,TSPYL5,UBE2N,MRPL42,SOCS2,CRADD,OR2B3,OR5U1,OR5V1,OR12D3,LOC651503,ZNF783,PDIA4,ZNF786,ZNF425,ZNF398,ZNF212,ZNF282,KITLG,CDH13,HSBP1,EFCBP2,MBTPS1,OSGIN1,MLYCD,LOC146167,ADRBK2,NAALADL2,FTMT,DIAPH3' %>% strsplit(., split = ',') %>% unlist() %>% unique(),
  easia = 'RANBP2,FLJ32745,EDAR,ERBB4,HEY2,NCOA7,SLC5A7,SULT1C3,SULT1C1,SULT1C2,GCC2,PCDH15,EXOC6,CYP26C1,CYP26A1,FER1L3,CEP55,APC,SRP19,DCP2,MCC,REEP5,ERC1,MRPL27,EME1,LRRC59,FLJ20920,CHAD,RSAD1,EPN3,SPATA20,CACNA1G,ABCC3,ANKRD40,CROP,XYLT2,MYCBPAP,CBR4,SH3RF1,NEK1,NRG1,MAP3K7IP1,MGAT3,SMCR7L,ATF4,RPS19BP1,CACNA1I,FER,PRKG1,KIF11,HHEX,EXOC6,LOC284757,CRB2,DENND1A,XYLT1,GRM3,KIAA1324L,ITGA6,PDK1,CYP1B1,MGC34824,ARL6IP2,LEPR,PDE4B,CXXC4,DAPK2,FAM96A,SNX1,SNX22,PPIB,CSNK1G1,SDK1,CUGBP2,ADAMTS18,SLC44A5,ACADM,RABGGTB,MSH4,EPB41L4B,C9orf4,C9orf5,CTNNAL1,COL11A1,DYNLT1,VIL2,LOC202459,SYTL3,FOXP1,FLJ34931,CLIP4,ALK,APOL3,APOL4,APOL2,APOL1,MYH9,TXN2,FOXRED2,EIF3S7,USP38,GAB1,ARHGEF7,C13orf16,CADPS,C10orf38,ITGA8,GRM8,ZNF800,GCC1,ARF5,FSCN3,PAX4,NAALADL2,ZNF277P,IFRD1,FLJ39575,NEK10,SLC4A7,UBE2D2,CXXC5,PSD2,NRG2,ZPBP,SERPINI1,GOLPH4,C9orf93,IRF2BP2,TNNI3K,C1orf173,PDE11A,C6orf10,TNXB,CREBL1,FKBPL,PRRT1,PPT2,EGFL8,AGPAT1,RNF5,AGER,GPSM3,NOTCH4,LOC401252,PBX2,TNFSF18,TNFSF4,FAM44B,LIN28B,BVES,POPDC3,PREP,ZC3HAV1L,ZC3HAV1,TM4SF19,FLJ25996,ZDHHC19,OSTalpha,PCYT1A,MGC33212,C8orf15,MTMR9,AMAC1L2,XKR6,ASCC1,C10orf104,DDIT4,DNAJB12,CBARA1,BOC,WDR52,CCDC52,THADA,PLEKHH2' %>% strsplit(., split = ',') %>% unlist() %>% unique(),
  america = 'PRM1,KIAA0350,SOCS1,TNP2,PRM3,PRM2,C16orf75,LOC400499,GPR115,GPR111,CD2AP,C1QDC1,TSPAN11,HSF2,SERINC1,PKIB,HDAC3,C5orf16,FCHSD1,CENTD3,PCDH1,DIAPH1,TMPRSS11F,CHR415SYT,TMPRSS11B,FUT10,RBM13,C8orf41,RNF122,DUSP26,C16orf68,BCAS1,CYP24A1,FLJ16641,NRCAM,PNPLA8,THAP5,DNAJB9,FCGR2A,HSPA6,FCGR3A,FCGR2C,FCGR3B,FCGR2B,FCRLA,FCRLB,DUSP12,ATF6,OLFML2B,C21orf87,BRWD1,HMGN1,WRB,C21orf13,SH3BGR,LOC728776,PKIB,FABP7,SMPDL3A,MUC17,TRIM56,SERPINE1,AP1S1,VGF,FLJ39237,MOGAT3,PLOD3,ZNHIT1,RIMS1,RAD54B,KIAA1429,RBM35A,DPY19L4,ZNF659,FAM83B,PRMT8,EFCAB4B,PARP11,NEK7,SSBP2,B3GALT5,PCP4,LOC150084,GALNT10,SAP30L,HAND1,CCND1,FLJ42258,ORAOV1,FGF19,TSSC1,TTC15,LHX1,AATF,ACACA,C22orf34,BRD1,ZBED4,SPTBN1,RTN4,FLJ42562,KIAA0409,ILK,TAF10,DCHS1,MRPL17,OR2AG2,OR2AG1,OR6A2,OR10A5,OR10A2,TPP1,KCNC2,PLXNA4B,SMC6,FLJ40869,KCNS3,MCMDC1,ASF1A,C6orf60,SSH1,DAO,SVOP,USP30,ALKBH2,UNG,WDR25,KIAA1446,C14orf70,FBXL17,LRRC55,AGTRL1,TNKS1BP1,BHLHB3,SSPN,ITPR2,ITGB6,RBMS1,FLJ21908,P11,RAPGEF3,RNF144,RSAD2,OPRM1,PIP3−E,FBXO47,PLXDC1,CACNB1,RPL19,STAC2,FHOD3,C18orf10,KIAA1328,ZNF438,KRT28,KRT40,KRT25,KRT26,KRT27,KRT10,KRT12,KRT20,KRT23,KRT39,KRTAP3−3,KRTAP3−2,KRTAP3−1,KRTAP1−5,KRTAP1−3,TMEM99,GPR113,SELI,C2orf39,OTOF,HADHB,THRB,VCX−C,KAL1,NOSTRIN,SPBC25,G6PC2,ABCB11,DHRS9,GPR116,GPR110,AQP9,CENPA,KCNK3,DPYSL5,MAPRE3,C2orf18,MAN1A1,GBP6,LRRC8B,APOBEC1,GDF3,DPPA3,CLEC4C,NANOG,SLC13A1,PTPRA,FAM113A,VPS16,GNRH2,MRPS26,OXT,AVP,UBOX5,CD302,LY75,PLA2R1,KIAA1370,ONECUT1,C6orf107,TAF11,ANKS1A,TCP11,SLC16A10,KIAA1919,REV3L,LARP2,PGRMC2,CYP39A1,TDRD6,PLA2G7,SLC25A27,ADAM10,SLTM,FAM63B' %>% strsplit(., split = ',') %>% unlist() %>% unique(),
  oceania = 'FAM77D,SUCLG2,GRB14,COBLL1,FLJ39822,SMOC1,SLC8A3,OR5U1,OR5V1,OR12D3,OR12D2,OR11A1,OR10C1,OR2H1,MAS1L,C3orf57,LOC131149,CRADD,SCN2A,SCN3A,FAM130A2,MTBP,SNTB1,TSSC1,TTC15,ADI1,RNASEH1,RPS7,COLEC11,RGS7BP,P18SRP,C20orf196,CHGB,CGI−09,MCM8,CRLS1,C20orf75,C20orf42,EDEM3,FAM129A,CPEB2,FLJ39743,IGF1R,C13orf21,CTGF,ENPP1,WDR72,C9orf126,PPP6C,RABEPK,HSPA5,DUPD1,DUSP13,SAMD8,VDAC2,MYST4,ANKRD10,ARHGEF7,WWOX,LRRC40,SFRS11,LRRC7,PPM1L,B3GALNT1,NMD3,SGCZ,KIF6,LRCH1,ESD,FAM77D,C4orf16,TIFA,ALPK1,ING3,FLJ21986,WNT16,FAM3C,IGF2R,PNLDC1,MAS1,PDE11A,AGPS,TTC30A,TTC30B,CCDC46,ACTN1,WDR22,ADAM11,DBF4B,HIGD1B,EFTUD2,CCDC103,GFAP,C1QL1,GJA7,LOC146909,ERBB4,FGL1,PCM1,ASAH1,CRB2,STRBP,DENND1A,LACE1,FOXO3A,LRRC18,C10orf72,XYLB,ENDOGL1,ACVR2B,SCN5A,ERBB4,SDCCAG10,ADAMTS6,LRRC3B,TNRC9,FGF14,C1orf21,ECM1,ADAMTSL4,C1orf138,MCL1,ENSA,GOLPH3L,HORMAD1,CTSS,CTSK,ARNT,RBM43,NMI,RASSF4,C10orf10,C10orf25,ZNF22,FRMPD4,ANKRD15,DMRT1,DMRT3,DMRT2,NFAM1,SERHL,CGI−96,FEM1C,TICAM2,TMED7,PTPRT,PHC3,PRKCI,SKIL,CLDN11,RFTN2,MARS2,BOLL,LHFPL3,ORC5L,FERD3L,C11orf69,CD59,FBXO3,LMO2,LCP2,KCNIP1,KCNMB1,ASTN1,FAM5B,ARFGAP3,PACSIN2,TTLL1,BIK,IKZF2,ZNF239,ZNF485,ZNF32,COL21A1,DST,ADAMTS20,COQ5,COX6A1,TRIAP1,15E1.2,SFRS9,DYNLL1,RNF10,POP5,CABP1,KIAA0152,ACADS,MGC5139,ANKRD12,TWSG1,PCSK2,BFSP1,ADK,PDE4D,STK31' %>% strsplit(., split = ',') %>% unlist() %>% unique()
)


# 
# hr_c <- hclust(dist(t(mat_d_list[[paste0('td','_neg')]][,-1:-4]), method="euclidean"), method = "complete")
# hr_r <- hclust(dist((mat_d_list[[paste0('td','_neg')]][,-1:-4]), method="euclidean"), method = "complete")
# mat_d_list[["td_neg"]] <- mat_d_list[["td_neg"]][hr_r$order,]
# 
# mat_d_list[["td_neg"]] %>% select(-contains('chrom')) %>%  mutate(posid = as.numeric(row.names(.))) %>% gather("pop", "value",2:NCOL(.))%>%  mutate(pop = factor(pop)) %>% mutate(pop = factor(pop, levels(pop)[hr_c$order])) %>% ggplot(., aes(x = pop, y = posid)) + geom_tile(aes(fill = value)) + heatmap_col + theme_bw() + theme_heat + ggtitle("D.")

#  mat_nsl <- readRDS('~/data/NZ_coreExome_1kgp/haplotype/nsl_clustered_dist_mat-6-7-2017.RDS')
# clus <- hclust(dist(t(mat_nsl), method="euclidean"), method = "complete")
# ddata_x <- dendro_data(clus)
# ddata_x$segments <- ddata_x$segments %>% mutate(yend = ifelse(yend == 0, 0, yend))
# p2 <- ggplot(segment(ddata_x)) +
#   geom_segment(aes(x=x, y=y, xend=xend, yend=yend))
# labs <- label(ddata_x)
# labs$group <- unlist(lapply(as.character(labs$label), function(x){strsplit(x, split = '_')[[1]][1]}))
# p2 + geom_text(data=label(ddata_x),
#                aes(label=label, x=x, y=y-0.1, colour=labs$group, angle = 90))  + theme_dendro() + theme(legend.position = 'none')
# 
# # ggdendrogram( nsl_clus,theme_dendro = FALSE, rotate = TRUE,color = 'tomato') + ggtitle("A")
