library(tidyverse)

# the filtered files were filtered for P < 5e-5
p_thres <- 5e-5

self_as <-bind_rows(lapply(list.files('~/data/ukbiobank_gwas/filtered/', pattern = 'controlsself_age_sex_chr*', full.names = TRUE), function(f){
  read_delim(f, col_names = TRUE, delim = '\t', col_types = list(CHR = col_integer(), SNP = col_character(),  BP = col_integer(), A1 = col_character(), TEST = col_character(), NMISS = col_integer(),OR = col_double(), SE = col_double(), L95 = col_double(), U95 = col_double(), STAT = col_double(), P = col_double())) %>% filter(P < p_thres)})) 

self_asb <- bind_rows(lapply(list.files('~/data/ukbiobank_gwas/filtered/', pattern = 'controlsself_age_sex_bmi_chr*', full.names = TRUE), function(f){read_delim(f,col_names = TRUE, delim = '\t', col_types = list(CHR = col_integer(), SNP = col_character(),  BP = col_integer(), A1 = col_character(), TEST = col_character(), NMISS = col_integer(),OR = col_double(), SE = col_double(), L95 = col_double(), U95 = col_double(), STAT = col_double(), P = col_double())) %>% filter(P < p_thres)}))


winnard_as <- bind_rows(lapply(list.files('~/data/ukbiobank_gwas/filtered/', pattern = 'controlswinnard_age_sex_chr*', full.names = TRUE), function(f){read_delim(f,col_names = TRUE, delim = '\t', col_types = list(CHR = col_integer(), SNP = col_character(),  BP = col_integer(), A1 = col_character(), TEST = col_character(), NMISS = col_integer(),OR = col_double(), SE = col_double(), L95 = col_double(), U95 = col_double(), STAT = col_double(), P = col_double())) %>% filter(P < p_thres)}))

winnard_asb <- bind_rows(lapply(list.files('~/data/ukbiobank_gwas/filtered/', pattern = 'controlswinnard_age_sex_bmi_chr*', full.names = TRUE),function(f){read_delim(f,col_names = TRUE, delim = '\t', col_types = list(CHR = col_integer(), SNP = col_character(),  BP = col_integer(), A1 = col_character(), TEST = col_character(), NMISS = col_integer(),OR = col_double(), SE = col_double(), L95 = col_double(), U95 = col_double(), STAT = col_double(), P = col_double())) %>% filter(P < p_thres)}))

selfult_as <- bind_rows(lapply(list.files('~/data/ukbiobank_gwas/filtered/', pattern = 'controlsself_ult_age_sex_chr*', full.names = TRUE),function(f){read_delim(f,col_names = TRUE, delim = '\t', col_types = list(CHR = col_integer(), SNP = col_character(),  BP = col_integer(), A1 = col_character(), TEST = col_character(), NMISS = col_integer(),OR = col_double(), SE = col_double(), L95 = col_double(), U95 = col_double(), STAT = col_double(), P = col_double())) %>% filter(P < p_thres)}))

selfult_asb <- bind_rows(lapply(list.files('~/data/ukbiobank_gwas/filtered/', pattern = 'controlsself_ult_age_sex_bmi_chr*', full.names = TRUE),function(f){read_delim(f,col_names = TRUE, delim = '\t', col_types = list(CHR = col_integer(), SNP = col_character(),  BP = col_integer(), A1 = col_character(), TEST = col_character(), NMISS = col_integer(),OR = col_double(), SE = col_double(), L95 = col_double(), U95 = col_double(), STAT = col_double(), P = col_double())) %>% filter(P < p_thres)}))

hosp_as <- bind_rows(lapply(list.files('~/data/ukbiobank_gwas/filtered/', pattern = 'controlshosp_age_sex_chr*', full.names = TRUE),function(f){read_delim(f,col_names = TRUE, delim = '\t', col_types = list(CHR = col_integer(), SNP = col_character(),  BP = col_integer(), A1 = col_character(), TEST = col_character(), NMISS = col_integer(),OR = col_double(), SE = col_double(), L95 = col_double(), U95 = col_double(), STAT = col_double(), P = col_double())) %>% filter(P < p_thres)}))

hosp_asb <- bind_rows(lapply(list.files('~/data/ukbiobank_gwas/filtered/', pattern = 'controlshosp_age_sex_bmi_chr*', full.names = TRUE),function(f){read_delim(f,col_names = TRUE, delim = '\t', col_types = list(CHR = col_integer(), SNP = col_character(),  BP = col_integer(), A1 = col_character(), TEST = col_character(), NMISS = col_integer(),OR = col_double(), SE = col_double(), L95 = col_double(), U95 = col_double(), STAT = col_double(), P = col_double())) %>% filter(P < p_thres)}))

all_as <- bind_rows(lapply(list.files('~/data/ukbiobank_gwas/filtered/', pattern = 'controlsall_age_sex_chr*', full.names = TRUE),function(f){read_delim(f,col_names = TRUE, delim = '\t', col_types = list(CHR = col_integer(), SNP = col_character(),  BP = col_integer(), A1 = col_character(), TEST = col_character(), NMISS = col_integer(),OR = col_double(), SE = col_double(), L95 = col_double(), U95 = col_double(), STAT = col_double(), P = col_double())) %>% filter(P < p_thres)}))

all_asb <- bind_rows(lapply(list.files('~/data/ukbiobank_gwas/filtered/', pattern = 'controlsall_age_sex_bmi_chr*', full.names = TRUE),function(f){read_delim(f,col_names = TRUE, delim = '\t', col_types = list(CHR = col_integer(), SNP = col_character(),  BP = col_integer(), A1 = col_character(), TEST = col_character(), NMISS = col_integer(),OR = col_double(), SE = col_double(), L95 = col_double(), U95 = col_double(), STAT = col_double(), P = col_double())) %>% filter(P < p_thres)}))

ult_as <- bind_rows(lapply(list.files('~/data/ukbiobank_gwas/filtered/', pattern = 'controlsult_age_sex_chr*', full.names = TRUE),function(f){read_delim(f,col_names = TRUE, delim = '\t', col_types = list(CHR = col_integer(), SNP = col_character(),  BP = col_integer(), A1 = col_character(), TEST = col_character(), NMISS = col_integer(),OR = col_double(), SE = col_double(), L95 = col_double(), U95 = col_double(), STAT = col_double(), P = col_double())) %>% filter(P < p_thres)}))

ult_asb <- bind_rows(lapply(list.files('~/data/ukbiobank_gwas/filtered/', pattern = 'controlsult_age_sex_bmi_chr*', full.names = TRUE),function(f){read_delim(f,col_names = TRUE, delim = '\t', col_types = list(CHR = col_integer(), SNP = col_character(),  BP = col_integer(), A1 = col_character(), TEST = col_character(), NMISS = col_integer(),OR = col_double(), SE = col_double(), L95 = col_double(), U95 = col_double(), STAT = col_double(), P = col_double())) %>% filter(P < p_thres)}))


#hosp_asb[!hosp_asb$SNP %in% hosp_as$SNP ,]

#self_asb[!self_asb$SNP %in% self_as$SNP,]

