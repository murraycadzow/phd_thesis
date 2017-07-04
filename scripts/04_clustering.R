library(tidyverse)
library(ggdendro)


mat_d_list <- readRDS('~/data/NZ_coreExome_1kgp/30kbWindow_intra/30kbwindows_mat_d_list-28-6-2017.RDS')

create_windowTable <- function(){
  bind_rows(lapply(names(mat_d_list), function(x){mat_d_list[[x]] %>% select(-contains('chrom'), -posid) %>% gather(., pop, value, 1:NCOL(.)) %>% group_by(pop) %>% summarise(total_windows =  sum(value)) %>% mutate(d = x)})) %>% mutate(super = sapply(pop, function(x){strsplit(x, '_')[[1]][1]})) %>% select(-pop) %>% group_by(d, super) %>% summarise(min = min(total_windows), mean = mean(total_windows), max = max(total_windows)) %>% data.frame()
}

create_windowSummaryTable <- function(){
  bind_rows(lapply(names(mat_d_list), function(x){mat_d_list[[x]] %>% select(-posid, -contains('chrom')) %>% summarise(total_windows = NROW(.), min = min(rowSums(.)), mean = mean(rowSums(.)), max = max(rowSums(.)), median = median(rowSums(.)), sd = sd(rowSums(.)) ) %>% mutate( stat = x)})) %>% select(stat, total_windows, min, median, max, mean, sd) %>% data.frame()
}

theme_heat <- theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), legend.position = 'none', axis.title.y = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank(), axis.title.x = element_blank())
four_plot <- function(mat_d_name,statname){
  multiplot(
    ggdendrogram(hclust(dist(t(mat_d_list[[paste0(mat_d_name,'_neg')]][,-1:-4]), method="euclidean"), method = "complete")) + ggtitle("A.", paste(statname,"lower tail")),
    
    mat_d_list[[paste0(mat_d_name,'_neg')]] %>% select(-contains('chrom')) %>% gather("pop", "value",2:NCOL(.))%>%  ggplot(., aes(x = pop, y = posid)) + geom_tile(aes(fill = value)) +scale_fill_gradient2()+ theme_bw() + theme_heat + ggtitle("C."),
    
    ggdendrogram(hclust(dist(t(mat_d_list[[paste0(mat_d_name,'_pos')]][,-1:-4]), method="euclidean"), method = "complete"), theme_dendro = TRUE) + ggtitle("B.",paste(statname,"upper tail")),
    
    mat_d_list[[paste0(mat_d_name,'_pos')]] %>% select(-contains('chrom')) %>% gather("pop", "value",2:NCOL(.))%>%  ggplot(., aes(x = pop, y = posid)) + geom_tile(aes(fill = value)) +scale_fill_gradient2()+ theme_bw() + theme_heat + ggtitle("D.")
    , cols = 2) 
}

create_sums <- function(mat_d_name){
  mat_d_list[[mat_d_name]] %>% mutate(POLsum =select(., contains("POL")) %>% rowSums(.), EURsum = select(., contains("EUR")) %>% rowSums(.), EASsum = select(., contains("EAS")) %>% rowSums(.), SASsum = select(., contains("SAS")) %>% rowSums(.), AMRsum = select(., contains("AMR")) %>% rowSums(.), AFRsum = select(., contains("AFR")) %>% rowSums(.), ALLsum = select(., -contains('chrom'), -posid) %>% rowSums(.))
}

create_fst_dendro <- function(){
    return(ggdendrogram(readRDS("~/data/NZ_coreExome_1kgp/whole_chr_fst_popgenome_clust.3-7-2017.RDS"), method = "complete"))
}

create_gwas_cat_table <-function(){
  load('~/data/gwas_catalog/disease_ref_table-26-6-2017.RData')
  return(gwas_int_ref_table %>% select(-pmid, -PUBMEDID))
}