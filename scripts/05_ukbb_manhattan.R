library(tidyverse)
library(ggman)
library(qqman)
for(gout in c("all","winnard","hosp","ult","self","self_ult")){
l <-list()
for(f in list.files(path = '~/data/ukbiobank_gwas/',pattern = paste0("controls",gout,"_age*.+bmi*.+"), full.names = TRUE)){
l[[f]] <- read_delim(file = f, delim = '\t', col_types = list(
  CHR = col_integer(),
  SNP = col_character(),
  BP = col_integer(),
  A1 = col_character(),
  TEST = col_character(),
  NMISS = col_integer(),
  OR = col_double(),
  SE = col_double(),
  L95 = col_double(),
  U95 = col_double(),
  STAT = col_double(),
  P = col_double()
)) %>% filter(P < 0.1, TEST == "ADD") %>% select(-STAT, -L95,-U95)
}


dat <- do.call(rbind,l)

ggsave(ggmanhattan(dat), file = paste0('~/Git_repos/bookdown_thesis/images/05_selection_and_association/',gout,'_age_sex_bmi.png'), width = 7,scale = 2)




for(f in list.files(path = '~/data/ukbiobank_gwas/',pattern = paste0("controls",gout,"_age.+sex_chr.+"), full.names = TRUE)){
  l[[f]] <- read_delim(file = f, delim = '\t', col_types = list(
    CHR = col_integer(),
    SNP = col_character(),
    BP = col_integer(),
    A1 = col_character(),
    TEST = col_character(),
    NMISS = col_integer(),
    OR = col_double(),
    SE = col_double(),
    L95 = col_double(),
    U95 = col_double(),
    STAT = col_double(),
    P = col_double()
  )) %>% filter(P < 0.1, TEST == "ADD") %>% select(-STAT, -L95,-U95)
}


dat <- do.call(rbind,l)

ggsave(ggmanhattan(dat), file = paste0('~/Git_repos/bookdown_thesis/images/05_selection_and_association/',gout,'_age_sex.png'), width = 7, scale = 2)
}
