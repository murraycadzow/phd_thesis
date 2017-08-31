library(tidyverse)
load('~/data/ukbiobank_gout_bd_refined2016-04-26.RData')


# f.31.0.0 = sex, 1 = male
# f.48.0.0 = waist
# f.50.0.0 = height
# f.21001.0.0 = BMI
# f.21000.0.0 = ethnicity
# f.20002.* =  self report illness

ukbio_gout <- bd_refined_gout %>% filter() %>% select(goutaff, 'age' = f.21003.0.0,"bmi" = f.21001.0.0, 'sex' = f.31.0.0, 'height' = f.50.0.0, 'waist' = f.48.0.0, 'hip' = f.49.0.0, 'eth' = f.21000.0.0 ) %>% mutate(waist_height_ratio = waist / height, sex_factor = factor(if_else(sex == 0, "female", "male", missing = NULL)), cc = if_else(goutaff == 1, "case", "control", missing = NULL), waist_hip_ratio = waist / hip) %>% filter(!is.na(goutaff) & eth %in% c(1001,1002,1003)) %>% group_by(cc) 



ukbio_gout %>% summarise(per_sex = sum(sex)/NROW(.), mean_age = mean(age, na.rm = TRUE), mean_bmi = mean(bmi, na.rm = TRUE), mean_waistheight_ratio = mean(waist_height_ratio, na.rm = TRUE), mean_waisthip_ratio = mean(waist_hip_ratio, na.rm = TRUE), n = n())

table(ukbio_gout$eth)

models <- list(
  list(adjuster = "age", model = glm( goutaff ~  age, family = 'binomial',data = ukbio_gout)),
  list(adjuster = "sex", model = glm( goutaff ~ sex_factor, family = 'binomial',data = ukbio_gout)),
  list(adjuster = "age,sex", model = glm( goutaff ~ age + sex_factor, family = 'binomial',data = ukbio_gout)),
  list(adjuster = "bmi", model = glm( goutaff ~ bmi, family = 'binomial',data = ukbio_gout)),
  list(adjuster = "age,sex,bmi", model = glm( goutaff ~ age + sex_factor + bmi, family = 'binomial',data = ukbio_gout)),
  list(adjuster = "waist", model = glm( goutaff ~ waist, family = 'binomial',data = ukbio_gout)),
  list(adjuster = "waist_hip_ratio", model = glm( goutaff ~ waist_hip_ratio, family = 'binomial',data = ukbio_gout)),
  list(adjuster = "waist_height_ratio", model = glm( goutaff ~ waist_height_ratio, family = 'binomial',data = ukbio_gout)),
  list(adjuster = "sex,waist_height_ratio", model = glm( goutaff ~ sex_factor + waist_height_ratio, family = 'binomial',data = ukbio_gout)),
  list(adjuster = "age,sex,waist_height_ratio", model = glm( goutaff ~ age + sex_factor + waist_height_ratio, family = 'binomial',data = ukbio_gout)),
  list(adjuster = "age,sex,waist,waist_height_ratio", model = glm( goutaff ~ age + sex_factor + waist + waist_height_ratio, family = 'binomial',data = ukbio_gout)),
  list(adjuster = "age,sex,bmi,waist,waist_height_ratio", model = glm( goutaff ~ age + sex_factor + bmi + waist + waist_height_ratio, family = 'binomial',data = ukbio_gout)),
  list(adjuster = "eth", model = glm( goutaff ~  factor(eth), family = 'binomial',data = ukbio_gout)) 
)


bind_rows(lapply(models, function(x){cbind(adjuster = x$adjuster,broom::glance(x$model))})) %>% arrange(AIC)

lapply(models, function(x){broom::tidy(x$model)})

