library(tidyverse)
library(formattable)
library(pander)

low_cov_seq <- function(){
  gout_phenotypes <- read_delim('data/gout_phenotypes.txt', delim='\t')
  si <- data.frame(SUBJECT= "AT0115")
  tab <- gout_phenotypes %>% 
    select(SUBJECT, BMICALC,BMI, GOUT, ACRGOUTAFFSTAT, DIABETES, AGECOL,SEX, HEART,KIDNEY, FATTYLIVER,WEIGHT,HEIGHT,WAIST,SYSTOLIC,DIASTOLIC,ETHCLASS) %>% 
    inner_join(., si, "SUBJECT") %>% arrange(pop) %>% 
    group_by(pop) %>% 
    summarise( n = NROW(pop),
               Age = paste0(formattable(mean(AGECOL, na.rm=TRUE), digits = 2, format = 'f')," (", formattable(sd(AGECOL, na.rm=TRUE),digits = 2, format = 'f'),")"),
               `Sex (% Male)` = formattable(sum(SEX == 1, na.rm = TRUE) / sum(SEX %in% c(1,2), na.rm =TRUE) * 100, digits = 2, format = 'f'),
               BMI = paste0(formattable(mean(BMI, na.rm=TRUE),digits = 2, format = 'f'), " (", formattable(sd(BMI, na.rm = TRUE),digits = 2, format = 'f'),")"),
               Waist = paste0(formattable(mean(WAIST, na.rm=TRUE),digits = 2, format = 'f'), " (", formattable(sd(WAIST, na.rm=TRUE),digits = 2, format ='f'), ")"),
               `Diabetes (%)` = formattable(sum(DIABETES == 2, na.rm=TRUE) / NROW(DIABETES) * 100, digits = 2, format = 'f') ,
               `Gout (%)` = formattable(sum(GOUT == 2, na.rm=TRUE) / NROW(GOUT)  * 100, digits = 2, format = 'f')
    )
  return(tab)
}

ce_clinical_table <- function() {
  gout_phenotypes <- read_delim('data/gout_phenotypes.txt', delim='\t')
  tab <- gout_phenotypes %>% 
    select(SUBJECT, BMICALC,BMI, GOUT, ACRGOUTAFFSTAT, DIABETES, AGECOL,SEX, HEART,KIDNEY, FATTYLIVER,WEIGHT,HEIGHT,WAIST,SYSTOLIC,DIASTOLIC,ETHCLASS) %>% 
    filter(!is.na(GOUT), !is.na(SEX), ETHCLASS %in% c(1,2,3)) %>% # should I ensure Gout and sex are complete?
    group_by(ETHCLASS,GOUT) %>% 
    summarise(Age = paste0(round(mean(AGECOL, na.rm=TRUE),2)," (", round(sd(AGECOL, na.rm=TRUE),2),")"),
              Sex = paste0(sum(SEX ==1, na.rm = TRUE), " (",round(sum(SEX ==1, na.rm = TRUE)/ sum(SEX %in% c(1,2), na.rm =TRUE) * 100,2),")"),
              BMI = paste0(round(mean(BMI, na.rm=TRUE),2), " (", round(sd(BMI, na.rm = TRUE),2),")"),
              Waist = paste0(round(mean(WAIST, na.rm=TRUE),2), " (", round(sd(WAIST, na.rm=TRUE),2), ")"),
              n_diabetes = paste0(sum(DIABETES == 2, na.rm=TRUE)," (",round(sum(DIABETES == 2, na.rm=TRUE)/NROW(GOUT)*100,2),")"),
              n_gout= NROW(GOUT),
              n_kidney = paste0(sum(KIDNEY==2, na.rm=TRUE)," (", round(sum(KIDNEY==2, na.rm=TRUE)/NROW(GOUT)*100, 2),")"),
              n_heart = paste0(sum(HEART == 2, na.rm = TRUE), " (",round(sum(HEART == 2, na.rm = TRUE)/NROW(GOUT) * 100,2),")"),
              n_Fatty_liver = paste0(sum(FATTYLIVER==2, na.rm=TRUE)," (",round(sum(FATTYLIVER == 2, na.rm = TRUE)/NROW(GOUT)*100,2),")"),
              #percentMale = signif(sum(SEX == 1)/NROW(SEX),digits = 2) * 100,
              meanSystolic = paste0(round(mean(SYSTOLIC,na.rm=TRUE), 1), " (",round(sd(SYSTOLIC,na.rm=TRUE), 1) ,")"),
              meanDiastolic = paste0(round(mean(DIASTOLIC, na.rm=TRUE),1), " (", round(sd(DIASTOLIC,na.rm=TRUE), 1),")")
              ) %>% ungroup() %>% 
    mutate(GOUT = ifelse(GOUT == 2, "Gout", "Control"), ETHCLASS = unlist(lapply(ETHCLASS, function(x){switch(x, "1"={return("EP")}, "2" ={return("WP")}, "3" = {return("CAU")})}))) %>% unite(category, c("ETHCLASS","GOUT"), sep = '_') %>% gather("var","value", -category) %>% spread(category, value) %>% mutate(var = c('Mean age, years (SD)', "Mean BMI (SD)", "Mean Systolic BP (SD)", "Mean Diastolic BP (SD)","Diabetes, n (%)", "Fatty liver, n (%)", "n", "Heart disease, n (%)", "Kidney disease, n (%)", "Sex, n male (% male)", "Mean Waist cm (SD)"))
 return(tab[c(6,9,4,5,7,8,1:3,10),])
}


ce_populations_table <- function() {
  # return table with:
  # n
  # % male
  # % gout
  # mean BMI (sd)
  # % diabetes
  load('data/ce_selected_ids_27-4-2017.Rdata')
  gout_phenotypes <- read_delim('data/gout_phenotypes.txt', delim='\t')
  si <- bind_rows(lapply(names(selected_ids), function(x){selected_ids[[x]]$pop <- x; return(selected_ids[[x]])})) %>% filter(!pop %in% c('WPN','EPN','NAD'))
  tab <- gout_phenotypes %>% 
    select(SUBJECT, BMICALC,BMI, GOUT, ACRGOUTAFFSTAT, DIABETES, AGECOL,SEX, HEART,KIDNEY, FATTYLIVER,WEIGHT,HEIGHT,WAIST,SYSTOLIC,DIASTOLIC,ETHCLASS) %>% 
    inner_join(., si, "SUBJECT") %>% arrange(pop) %>% 
    group_by(pop) %>% 
    summarise( n = NROW(pop),
      Age = paste0(formattable(mean(AGECOL, na.rm=TRUE), digits = 2, format = 'f')," (", formattable(sd(AGECOL, na.rm=TRUE),digits = 2, format = 'f'),")"),
      `Sex (% Male)` = round(sum(SEX == 1, na.rm = TRUE) / sum(SEX %in% c(1,2), na.rm =TRUE) * 100, digits = 0),
            BMI = paste0(formattable(mean(BMI, na.rm=TRUE),digits = 2, format = 'f'), " (", formattable(sd(BMI, na.rm = TRUE),digits = 2, format = 'f'),")"),
            Waist = paste0(formattable(mean(WAIST, na.rm=TRUE),digits = 2, format = 'f'), " (", formattable(sd(WAIST, na.rm=TRUE),digits = 2, format ='f'), ")"),
            `Diabetes (%)` = round(sum(DIABETES == 2, na.rm=TRUE) / NROW(DIABETES) * 100, digits = 0) ,
            `Gout (%)` = round(sum(GOUT == 2, na.rm=TRUE) / NROW(GOUT)  * 100, digits = 0),
      `Kidney (%)` = round(sum(KIDNEY == 2, na.rm=TRUE) / NROW(KIDNEY)  * 100, digits = 0),
      `Heart (%)` = round(sum(HEART == 2, na.rm=TRUE) / NROW(HEART)  * 100, digits = 0)
    )
  return(tab)
  
}
