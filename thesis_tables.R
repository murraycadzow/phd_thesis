## ------------------------------------------------------------------------
library(tidyverse)
library(pander)
gout_phenotypes <- read_delim('~/gout_phenotypes.txt', delim='\t')


## ---- echo = FALSE, results='asis'---------------------------------------
tab <- gout_phenotypes %>% select(SUBJECT, BMICALC,BMI, GOUT, ACRGOUTAFFSTAT, DIABETES, AGECOL,SEX, HEART,KIDNEY, FATTYLIVER,WEIGHT,HEIGHT,WAIST,SYSTOLIC,DIASTOLIC,ETHCLASS) %>% filter(!is.na(GOUT), !is.na(SEX), ETHCLASS %in% c(1,2,3)) %>% 
  group_by(ETHCLASS,GOUT) %>% 
  summarise(Age = paste0(round(mean(AGECOL, na.rm=TRUE),2)," (", round(sd(AGECOL, na.rm=TRUE),2),")"), 
            BMI = paste0(round(mean(BMI, na.rm=TRUE),2), " (", round(sd(BMI, na.rm = TRUE),2),")"),
            Waist = paste0(round(mean(WAIST, na.rm=TRUE),2), " (", round(sd(WAIST, na.rm=TRUE),2), ")"),
            n_diabetes = sum(DIABETES == 2, na.rm=TRUE),
            n_gout= NROW(GOUT),
            n_kidney = sum(KIDNEY==2, na.rm=TRUE),
            n_Fatty_liver = sum(FATTYLIVER, na.rm=TRUE),
            #percentMale = signif(sum(SEX == 1)/NROW(SEX),digits = 2) * 100,
            meanBP = paste0(round(mean(SYSTOLIC,na.rm=TRUE), 1),'/',round(mean(DIASTOLIC, na.rm=TRUE),1))) %>% ungroup() %>% 
  mutate(GOUT = ifelse(GOUT == 2, "Gout", "Control"), ETHCLASS = unlist(lapply(ETHCLASS, function(x){switch(x, "1"={return("EP")}, "2" ={return("WP")}, "3" = {return("CAU")})}))) %>% t() 
colnames(tab) <- c(1,2,3,4,5,6)
knitr::kable(tab, longtable = TRUE, format = 'latex', booktabs = TRUE )


