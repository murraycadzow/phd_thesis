library(tidyverse)
library(RColorBrewer)
library(scales)

set.seed(713)
scratch_dir <- '/media/xsan/scratch/merrimanlab/murray/working_dir/coreExome_selection/NZ_coreExome_1kgp/data/NZ_coreExome/'
archive_dir <- '/media/xsan/archive/merrimanlab/murray/Bioinformatics/Projects/coreExome_selection/NZ_coreExome_1kgp/data/NZ_coreExome/'


cim <- read.table(paste0(archive_dir,'CIM.panel'), stringsAsFactors = FALSE) %>% select('SUBJECT' = V1, 'pop' = V2, 'super_pop' = V3)
nzm <- read.table(paste0(archive_dir,'NZM.panel'), stringsAsFactors = FALSE)%>% select('SUBJECT' = V1, 'pop' = V2, 'super_pop' = V3)
sam <- read.table(paste0(archive_dir,'SAM.panel'), stringsAsFactors = FALSE)%>% select('SUBJECT' = V1, 'pop' = V2, 'super_pop' = V3)
ton <- read.table(paste0(archive_dir,'TON.panel'), stringsAsFactors = FALSE)%>% select('SUBJECT' = V1, 'pop' = V2, 'super_pop' = V3)
nzc <- read.table(paste0(archive_dir,'NZC.panel'), stringsAsFactors = FALSE)%>% select('SUBJECT' = V1, 'pop' = V2, 'super_pop' = V3)

exclude_peeps <- read.delim(paste0(archive_dir,'src_data/QC1_7-BlanketExclusions.txt'),header=FALSE, stringsAsFactors = FALSE)

euro_samples <- read.delim(paste0(archive_dir,'src_data/QC1_7-FullMerge1.0and1.1_ExDup_MattBrown_Czech-plus_correctAff.fam'), stringsAsFactors = FALSE, header=FALSE, sep=' ')
names(euro_samples) <- c('FID','IID','PID','MID','SEX','AFF')

poly_samples <- read.delim(paste0(archive_dir,'src_data/QC1_7-FullMerge1.0and1.1_ExDup_MattBrown_Czech-plus_polygoutAff.fam'), stringsAsFactors = FALSE, header=FALSE)
names(poly_samples) <- c('FID','IID','PID','MID','SEX','AFF')

pca_ancestry <- read.delim(paste0(archive_dir,"src_data/PCA_afterRuthsApp_210317.txt"), stringsAsFactors = FALSE, header =TRUE, sep ='\t')

table(duplicated(euro_samples$IID))
table(euro_samples$AFF, exclude = NULL)
table(poly_samples$AFF, exclude = NULL)

euro_pca <- pca_ancestry %>% inner_join(., euro_samples %>% select(SUBJECT=IID, AFF), by = 'SUBJECT') %>% filter(AFF != -9 & !is.na(AFF) & !Ethnicity.Matching.PCA  %in% c("No","unknown") & STUDYCOHORTNAME == "Gout" & !SUBJECT %in% exclude_peeps$V2)


poly_pca <- pca_ancestry %>% inner_join(., poly_samples %>% select(SUBJECT=IID, AFF), by = 'SUBJECT') %>% filter(AFF != -9 & !is.na(AFF) & !Ethnicity.Matching.PCA  %in% c("No","unknown") & STUDYCOHORTNAME == "Gout" & !SUBJECT %in% exclude_peeps$V2)


ggsave(rbind(euro_pca,poly_pca) %>% inner_join(., bind_rows(cim,nzm,sam,ton,nzc), by = 'SUBJECT') %>% filter(!SRSpecific %in% c("","African","Iberian","East Asian", "Native American","South Asian", "Melanesian","Middle Eastern"))%>%  slice(grep("Maori|European|Samoan|Tongan", SRSpecific)) %>% group_by(SRSpecific) %>% mutate(Population = ifelse(SRSpecific == 'CI Maori/Tahitian', "CI Maori",SRSpecific )) %>% ggplot(., aes(x = PCA2, y =PCA4, colour = Population)) + geom_point(aes(alpha = 0.3)) + theme_bw() + scale_alpha(guide = FALSE) + xlab("PC2") + ylab("PC4"), width = 6.5, height = 6.5, filename = '~/Git_repos/bookdown_thesis/images/02_methods/pca_plot.png')

#%>% mutate(pca2mean = mean(PCA2), pca2sd = sd(PCA2), pca4mean = mean(PCA4), pca4sd = sd(PCA4)) %>%  filter(filter_pca(PCA2, 1) & filter_pca(PCA4, 1))
