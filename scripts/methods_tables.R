library(tidyverse)
library(formattable)
library(pander)

ce_clinical_table <- function() {
  gout_phenotypes <- read_delim('data/gout_phenotypes.txt', delim='\t')
  tab <- gout_phenotypes %>% 
    select(SUBJECT, BMICALC,BMI, GOUT, ACRGOUTAFFSTAT, DIABETES, AGECOL,SEX, HEART,KIDNEY, FATTYLIVER,WEIGHT,HEIGHT,WAIST,SYSTOLIC,DIASTOLIC,ETHCLASS) %>% 
    filter(!is.na(GOUT), !is.na(SEX), ETHCLASS %in% c(1,2,3)) %>% # should I ensure Gout and sex are complete?
    group_by(ETHCLASS,GOUT) %>% 
    summarise(Age = paste0(round(mean(AGECOL, na.rm=TRUE),2)," (", round(sd(AGECOL, na.rm=TRUE),2),")"),
              Sex = round(sum(SEX ==1, na.rm = TRUE)/ sum(SEX %in% c(1,2), na.rm =TRUE),2),
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
 return(tab)
}


ce_populations_table <- function() {
  # return table with:
  # n
  # % male
  # % gout
  # mean BMI (sd)
  # % diabetes
  load('data/ce_selected_ids_11-4-2017.Rdata')
  gout_phenotypes <- read_delim('data/gout_phenotypes.txt', delim='\t')
  si <- bind_rows(lapply(names(selected_ids), function(x){selected_ids[[x]]$pop <- x; return(selected_ids[[x]])}))
  tab <- gout_phenotypes %>% 
    select(SUBJECT, BMICALC,BMI, GOUT, ACRGOUTAFFSTAT, DIABETES, AGECOL,SEX, HEART,KIDNEY, FATTYLIVER,WEIGHT,HEIGHT,WAIST,SYSTOLIC,DIASTOLIC,ETHCLASS) %>% 
    inner_join(., si, "SUBJECT") %>% arrange(pop) %>% 
    group_by(pop) %>% 
    summarise( n = NROW(pop),
      Age = paste0(round(mean(AGECOL, na.rm=TRUE),2)," (", round(sd(AGECOL, na.rm=TRUE),2),")"),
      Sex = round(sum(SEX == 1, na.rm = TRUE) / sum(SEX %in% c(1,2), na.rm =TRUE), 2),
            BMI = paste0(round(mean(BMI, na.rm=TRUE),2), " (", round(sd(BMI, na.rm = TRUE),2),")"),
            Waist = paste0(round(mean(WAIST, na.rm=TRUE),2), " (", round(sd(WAIST, na.rm=TRUE),2), ")"),
            per_diabetes = sum(DIABETES == 2, na.rm=TRUE) / NROW(DIABETES),
            per_gout= sum(GOUT == 2, na.rm=TRUE) / NROW(GOUT)
    )
  return(tab)
  
}