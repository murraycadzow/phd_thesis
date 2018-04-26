library(tidyverse)

# the filtered files were filtered for P < 1e-5
p_thres <- 1e-3

# self_as <-bind_rows(lapply(list.files('~/data/ukbiobank_gwas/filtered/', pattern = 'controlsself_age_sex_chr*', full.names = TRUE), function(f){
#   read_delim(f, col_names = TRUE, delim = '\t', col_types = list(CHR = col_integer(), SNP = col_character(),  BP = col_integer(), A1 = col_character(), TEST = col_character(), NMISS = col_integer(),OR = col_double(), SE = col_double(), L95 = col_double(), U95 = col_double(), STAT = col_double(), P = col_double())) %>% filter(P < p_thres)})) 

self_asb <- bind_rows(lapply(list.files('~/data/ukbiobank_gwas/filtered_p0.05/', pattern = 'controlsself_age_sex_bmi_chr*', full.names = TRUE), function(f){read_delim(f,col_names = TRUE, delim = '\t', col_types = list(CHR = col_integer(), SNP = col_character(),  BP = col_integer(), A1 = col_character(), TEST = col_character(), NMISS = col_integer(),OR = col_double(), SE = col_double(), L95 = col_double(), U95 = col_double(), STAT = col_double(), P = col_double())) %>% filter(P < p_thres, nchar(A1) ==1 , TEST == 'ADD')}))


# winnard_as <- bind_rows(lapply(list.files('~/data/ukbiobank_gwas/filtered/', pattern = 'controlswinnard_age_sex_chr*', full.names = TRUE), function(f){read_delim(f,col_names = TRUE, delim = '\t', col_types = list(CHR = col_integer(), SNP = col_character(),  BP = col_integer(), A1 = col_character(), TEST = col_character(), NMISS = col_integer(),OR = col_double(), SE = col_double(), L95 = col_double(), U95 = col_double(), STAT = col_double(), P = col_double())) %>% filter(P < p_thres, nchar(A1) ==1 )}))

winnard_asb <- bind_rows(lapply(list.files('~/data/ukbiobank_gwas/filtered_p0.05/', pattern = 'controlswinnard_age_sex_bmi_chr*', full.names = TRUE),function(f){read_delim(f,col_names = TRUE, delim = '\t', col_types = list(CHR = col_integer(), SNP = col_character(),  BP = col_integer(), A1 = col_character(), TEST = col_character(), NMISS = col_integer(),OR = col_double(), SE = col_double(), L95 = col_double(), U95 = col_double(), STAT = col_double(), P = col_double())) %>% filter(P < p_thres, nchar(A1) ==1 )}))

# selfult_as <- bind_rows(lapply(list.files('~/data/ukbiobank_gwas/filtered/', pattern = 'controlsself_ult_age_sex_chr*', full.names = TRUE),function(f){read_delim(f,col_names = TRUE, delim = '\t', col_types = list(CHR = col_integer(), SNP = col_character(),  BP = col_integer(), A1 = col_character(), TEST = col_character(), NMISS = col_integer(),OR = col_double(), SE = col_double(), L95 = col_double(), U95 = col_double(), STAT = col_double(), P = col_double())) %>% filter(P < p_thres, nchar(A1) ==1 )}))

selfult_asb <- bind_rows(lapply(list.files('~/data/ukbiobank_gwas/filtered_p0.05/', pattern = 'controlsself_ult_age_sex_bmi_chr*', full.names = TRUE),function(f){read_delim(f,col_names = TRUE, delim = '\t', col_types = list(CHR = col_integer(), SNP = col_character(),  BP = col_integer(), A1 = col_character(), TEST = col_character(), NMISS = col_integer(),OR = col_double(), SE = col_double(), L95 = col_double(), U95 = col_double(), STAT = col_double(), P = col_double())) %>% filter(P < p_thres, nchar(A1) ==1 )}))

# hosp_as <- bind_rows(lapply(list.files('~/data/ukbiobank_gwas/filtered/', pattern = 'controlshosp_age_sex_chr*', full.names = TRUE),function(f){read_delim(f,col_names = TRUE, delim = '\t', col_types = list(CHR = col_integer(), SNP = col_character(),  BP = col_integer(), A1 = col_character(), TEST = col_character(), NMISS = col_integer(),OR = col_double(), SE = col_double(), L95 = col_double(), U95 = col_double(), STAT = col_double(), P = col_double())) %>% filter(P < p_thres, nchar(A1) ==1 )}))

hosp_asb <- bind_rows(lapply(list.files('~/data/ukbiobank_gwas/filtered_p0.05/', pattern = 'controlshosp_age_sex_bmi_chr*', full.names = TRUE),function(f){read_delim(f,col_names = TRUE, delim = '\t', col_types = list(CHR = col_integer(), SNP = col_character(),  BP = col_integer(), A1 = col_character(), TEST = col_character(), NMISS = col_integer(),OR = col_double(), SE = col_double(), L95 = col_double(), U95 = col_double(), STAT = col_double(), P = col_double())) %>% filter(P < p_thres, nchar(A1) ==1,TEST == 'ADD' )}))

# all_as <- bind_rows(lapply(list.files('~/data/ukbiobank_gwas/filtered/', pattern = 'controlsall_age_sex_chr*', full.names = TRUE),function(f){read_delim(f,col_names = TRUE, delim = '\t', col_types = list(CHR = col_integer(), SNP = col_character(),  BP = col_integer(), A1 = col_character(), TEST = col_character(), NMISS = col_integer(),OR = col_double(), SE = col_double(), L95 = col_double(), U95 = col_double(), STAT = col_double(), P = col_double())) %>% filter(P < p_thres, nchar(A1) ==1 )}))

all_asb <- bind_rows(lapply(list.files('~/data/ukbiobank_gwas/filtered/', pattern = 'controlsall_age_sex_bmi_chr*', full.names = TRUE),function(f){read_delim(f,col_names = TRUE, delim = '\t', col_types = list(CHR = col_integer(), SNP = col_character(),  BP = col_integer(), A1 = col_character(), TEST = col_character(), NMISS = col_integer(),OR = col_double(), SE = col_double(), L95 = col_double(), U95 = col_double(), STAT = col_double(), P = col_double())) %>% filter(P < p_thres, nchar(A1) ==1,TEST == 'ADD' )}))

# ult_as <- bind_rows(lapply(list.files('~/data/ukbiobank_gwas/filtered/', pattern = 'controlsult_age_sex_chr*', full.names = TRUE),function(f){read_delim(f,col_names = TRUE, delim = '\t', col_types = list(CHR = col_integer(), SNP = col_character(),  BP = col_integer(), A1 = col_character(), TEST = col_character(), NMISS = col_integer(),OR = col_double(), SE = col_double(), L95 = col_double(), U95 = col_double(), STAT = col_double(), P = col_double())) %>% filter(P < p_thres, nchar(A1) ==1 )}))

ult_asb <- bind_rows(lapply(list.files('~/data/ukbiobank_gwas/filtered_p0.05/', pattern = 'controlsult_age_sex_bmi_chr*', full.names = TRUE),function(f){read_delim(f,col_names = TRUE, delim = '\t', col_types = list(CHR = col_integer(), SNP = col_character(),  BP = col_integer(), A1 = col_character(), TEST = col_character(), NMISS = col_integer(),OR = col_double(), SE = col_double(), L95 = col_double(), U95 = col_double(), STAT = col_double(), P = col_double())) %>% filter(P < p_thres, nchar(A1) ==1 , TEST == 'ADD')}))


#hosp_asb[!hosp_asb$SNP %in% hosp_as$SNP ,]

#self_asb[!self_asb$SNP %in% self_as$SNP,]

