#devtools::install_github("hms-dbmi/UpSetR")
library(UpSetR)
library(tidyverse)

# load('/media/xsan/staff_groups/merrimanlab/Merriman_Documents/Murray/ukbiobank_util/data/ukbiobank_genotyped2016-04-26.RData')
# all_fam <- read.table("~/Git_repos/bookdown_thesis/data/Ukbiobank/chr1impv1.fam", header=FALSE)
# colnames(all_fam) <- c("FID","IID","PID","MID", "SEX","AFF")
# 
# gout_patients <- genotyped %>% select(f.eid, goutaff, goutult, goutself, gouthosp, goutall, gout_exclude_drug, goutwinnard, renaldisease, f.21000.0.0)  %>% filter(goutaff == 1)
# 
# 
# 
# # A
# only_self <- gout_patients  %>% filter(goutself == 1 & is.na(goutult) & is.na(gouthosp)) %>% select(f.eid) %>% tally()
# 
# # B
# only_hosp <- gout_patients  %>% filter(gouthosp == 1 & is.na(goutult) & is.na(goutself)) %>% select(f.eid) %>% tally()
# 
# # C
# only_ult <- gout_patients  %>% filter(goutult == 1 & is.na(gouthosp) & is.na(goutself)) %>% select(f.eid) %>% tally()
# 
# # D
# union_hosp_ult <- gout_patients  %>% filter(gouthosp == 1 & goutult == 1 & is.na(goutself)) %>% select(f.eid) %>% tally()
# 
# # E
# union_hosp_self <- gout_patients  %>% filter(gouthosp == 1 & goutself == 1 & is.na(goutult)) %>% select(f.eid) %>% tally()
# 
# # F
# union_ult_self <- gout_patients  %>% filter(goutself == 1 & goutult == 1 & is.na(gouthosp)) %>% select(f.eid) %>% tally()
# 
# # G
# union_hosp_ult_self <- gout_patients  %>% filter(gouthosp == 1 & goutult == 1 & goutself == 1) %>% select(f.eid) %>% tally()
# 
# 
# 
# 
# gout_patients$f.eid <-  as.numeric(gout_patients$f.eid)
# ult <- gout_patients  %>%  filter(goutult == 1 & f.eid %in% all_fam$IID & f.21000.0.0 %in% c(1001,1002,1003) )  %>% select(f.eid)
# self <- gout_patients  %>%  filter(goutself == 1 & f.eid %in% all_fam$IID & f.21000.0.0 %in% c(1001,1002,1003))  %>% select(f.eid)
# hosp <- gout_patients  %>%  filter(gouthosp == 1 & f.eid %in% all_fam$IID & f.21000.0.0 %in% c(1001,1002,1003))  %>% select(f.eid)
# winnard <- gout_patients  %>%  filter(goutwinnard == 1 & f.eid %in% all_fam$IID & f.21000.0.0 %in% c(1001,1002,1003))  %>% select(f.eid)

listInput <- list(Hospital=hosp[,1], "Self-Reported" = self[,1], Winnard = winnard[,1], ULT = ult[,1])
#saveRDS(listInput, file = '~/Git_repos/bookdown_thesis/data/Ukbiobank/upsetPlotList.RDS')
upset(fromList(listInput), order.by = "freq", text.scale = 1.5 )


#venn.diagram(list(Hospital=hosp[,1], "Self-Reported" = self[,1], Winnard = winnard[,1], ULT = ult[,1]),fill = c('#FFFF80','#AAFFAA','#FFAAAA','#AAAAFF') ,filename = "genotyped_venn.tiff")