library(tidyverse)
library(ggdendro)
library(scales)
library(ggpubr)



###
# made in clustering/coreExome_1kg_popgenome_clustering_100kb_filtered_3ns.Rmd
mat_d_list <- readRDS('~/data/NZ_coreExome_1kgp/100kbWindow_intra/100kbwindows_filtered_3ns_mat_d_list-31-8-2017.RDS')
prop_list <- readRDS('~/data/NZ_coreExome_1kgp/100kbWindow_intra/100kbwindows_filtered_3ns_prop_list-31-8-2017.RDS')
clus_r_list <- readRDS('~/data/NZ_coreExome_1kgp/100kbWindow_intra/100kbwindows_filtered_3ns_clus_r_list-31-8-2017.RDS')
clus_c_list <- readRDS('~/data/NZ_coreExome_1kgp/100kbWindow_intra/100kbwindows_filtered_3ns_clus_c_list-31-8-2017.RDS')
###

panel <- read.delim(paste0('~/data/NZ_coreExome_1kgp/nz_1kgp.panel'), stringsAsFactors = FALSE)

heatmap_col <- scale_fill_gradient(low = "white", high = "steelblue", guide = 'colourbar', breaks = c(0.0, 0.25, 0.5, 0.75, 1.0 ), limits = c(0,1))

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

super_pop_colours <- cbind(panel %>% select(super_pop) %>% distinct() %>% arrange(super_pop) %>% select('group' = super_pop),colour = gg_color_hue(panel %>% select(super_pop) %>% distinct() %>% tally() %>% .[['n']]))

create_windowTable <- function(){
  bind_rows(lapply(names(mat_d_list)[!grepl('chr',names(mat_d_list))], function(x){mat_d_list[[x]] %>% select(-contains('chrom'), -posid) %>% gather(., pop, value, 1:NCOL(.)) %>% group_by(pop) %>% summarise(total_windows =  sum(value)) %>% mutate(d = x)})) %>% mutate(super = sapply(pop, function(x){strsplit(x, '_')[[1]][1]})) %>% select(-pop) %>% group_by(d, super) %>% summarise(min = min(total_windows), mean = mean(total_windows), max = max(total_windows)) %>% data.frame()
}

create_windowSummaryTable <- function(){
  bind_rows(lapply(names(mat_d_list)[!grepl('chr',names(mat_d_list))], function(x){mat_d_list[[x]] %>% select(-posid, -contains('chrom')) %>% summarise(total_windows = NROW(.), min = min(rowSums(.)), mean = mean(rowSums(.)), max = max(rowSums(.)), median = median(rowSums(.)), sd = sd(rowSums(.)) ) %>% mutate( stat = x)})) %>% select(stat, total_windows, min, median, max, mean, sd) %>% data.frame()
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
  mat_ihs <- readRDS('~/data/NZ_coreExome_1kgp/haplotype/ihs_clustered_dist_mat-6-7-2017.RDS')
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
  mat_nsl <- readRDS('~/data/NZ_coreExome_1kgp/haplotype/nsl_clustered_dist_mat-6-7-2017.RDS')
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
})) %>% filter(!pop %in% c("NAD","EPN","WPN","POL"))

# made from clustering/popgenome_freq_stats.Rmd
upper_sig_stats <-readRDS('~/data/NZ_coreExome_1kgp/100kbWindow_intra/filtered/100kbwindows_upper_sig_stat_genes_filtered_ns3_resized-31-8-2017.RDS')

upper_sig_stats <- bind_rows(lapply(names(upper_sig_stats), function(y){
  upper_sig_stats[[y]] <- bind_rows(
    lapply(names(upper_sig_stats[[y]]), function(x){
      upper_sig_stats[[y]][[x]] %>% mutate(pop = x, statname =y) 
    }
    ))
})) %>% filter(!pop %in% c("NAD","EPN","WPN","POL"))


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
