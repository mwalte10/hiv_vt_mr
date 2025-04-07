# install.packages('meta')
#install.packages('ggforestplot')
rm(list = ls())
library(data.table)
library(ggplot2)
library(stringr)
library(tidyverse)
library(tidybayes)
library(glmmTMB)
library(ggpubr)


setwd('C:/Users/mwalters/OneDrive - Imperial College London/Papers/VT estimation/paper/')
table_dir <- './tables/'
figure_dir <- './figures/'
input_data_dir <- './data/'
output_data_dir <- './data/model_data/'
results_dir <- './results/'

data <- fread(paste0(input_data_dir,'/public_data.csv'))
data[,studlab := paste0(author, '_', study_year)]
data[,id := 1:nrow(data)]
data[,NEW_2024 := ifelse(NEW_2024 == 'Y', 'Yes', 'No')]
data[inferred_ART_weeks < 5, model2_pvt_sc := 'less_four']
data[inferred_ART_weeks >4 & inferred_ART_weeks < 40, model2_pvt_sc := 'more_four']
data[inferred_ART_weeks == 40, model2_pvt_sc := 'pre_conception']
data[,weighted_prop := weighted.mean(event.e/ n.e, n.e), by = c('tt', 'model2_pvt_sc')]
data[,cd4_mid := (cd4_mid - 500) / 100]
data[,time_on_art_median := inferred_ART_weeks - 20]
data[,late_start := ifelse(inferred_ART_weeks < 5, T, F)]
data[,pvt_type := paste0(model2_pvt_sc, '_', tt)]
data[,on_art := ifelse(inferred_ART_weeks == 40, T, F)]

get_predictions <- function(model_name, study_data){
  if(model_name == 'model1'){
    model <- readRDS(paste0('./public_results/model_output/model_1.RDS'))
  }
  if(model_name == 'model2'){
    model <- readRDS(paste0('./public_results/model_output/model_2.RDS'))
  }
  if(model_name == 'model3'){
    model <- readRDS(paste0('./public_results/model_output/model_3.RDS'))
  }
  if(model_name == 'model4'){
    model <- readRDS(paste0('./public_results/model_output/model_4.RDS'))
  }
  
  stud_cov = predict(model, newdata = study_data, cov.fit = T)$cov.fit
  stud_cov <- mvtnorm::rmvnorm(3000, mean = predict(model, study_data), sigma = stud_cov)
  stud_cov <- data.table(stud_cov)[,draw := 1:3000]
  stud_cov <- melt(stud_cov, id.vars = 'draw')
  stud_cov[,value := plogis(value)]
  stud_cov <- stud_cov[,.(predictions = median(value), lower = quantile(value, 0.025), upper = quantile(value, 0.975)), by = 'variable']
  stud_cov <- stud_cov[,.(predictions, lower, upper, id = gsub(pattern = 'V', replacement = '', variable))]
  stud_cov[,id := study_data$id]
  
  preds <- merge(stud_cov, study_data, by = 'id')
  return(preds)
  
}

mod1 <- get_predictions(model_name = 'model1', study_data = data[model == 'model1'])
mod2 <- get_predictions(model_name = 'model2', study_data = data[model == 'model2'])
mod3 <- get_predictions(model_name = 'model3', study_data = data[model == 'model3'])
mod4 <- get_predictions(model_name = 'model4', study_data = data[model == 'model4'])
dt <- rbind(mod1, mod2, mod3, mod4)
dt[!is.na(cd4_mid), cd4_mid := (cd4_mid + 5) * 100]
dt[!is.na(time_on_art_median), time_on_art_median := (time_on_art_median + 20)]


stud_include <- fread(paste0(input_data_dir, '/studies_included_by_analysis_public.csv'))
stud_include <- stud_include[,.(studlab = paste0(author, '_', study_year), 
                                cd4_factor, tt = ifelse(type == 1, 'peri', 'bf'),
                                type = paste0(cd4_factor, '_', type),
                                study_year,
                                `Default`  = default, 
                                `2019 SR` = sr_2019,
                                `2024\nmethods update` = update_2024,
                                `2024\nSR update` = SR_2024,
                                reported_prop, 
                                lower_SI = lower, 
                                upper_SI = upper,
                                N_SI = n.e,
                                cd4_mid_SI = cd4_mid)]
stud_include <- unique(stud_include)
stud_include[`2024\nmethods update` != 'Included' & Default == '' & `2019 SR` == '', remove := T]
stud_include <- stud_include[is.na(remove),]
stud_include[,remove := NULL]
stud_include[lower_SI < 0, lower_SI := 0]
stud_include[upper_SI <0 , upper_SI := 0]
stud_include[type %in% c('A_1', 'B_1', 'C_1', 'A_2', 'B_2', 'C_2'), model := 'model1']
stud_include[type %in% c('D_1', 'D_2', 'H_1', 'H_2',
                         'I_1', 'I_2', 'J_1', 'J_2',
                         'K_1', 'K_2', 'L_2', 'M_2'), model := 'model2']
stud_include[type %in% c('E_1', 'F_1', 'G_1'), model := 'model3']
stud_include[type %in% c('E_2', 'F_2', 'G_2'), model := 'model4']


##Switch to a shorter word
# stud_include[,Default := ifelse(Default == 'Included', 'X','')]
# stud_include[,`2019 SR` := ifelse(`2019 SR` == 'Included', 'X','')]
# stud_include[,`2024\nmethods update` := ifelse(`2024\nmethods update` == 'Included', 'X','')]
# stud_include[,`2024\nSR update` := ifelse(`2024\nSR update` == 'Included', 'X','')]

dt <- merge(dt, stud_include, by = c('studlab', 'tt', 'cd4_factor', 'study_year', 'model'), all.x = T, all.y =T )
dt[`2024\nmethods update` != 'Included', predictions := reported_prop]
dt[`2024\nmethods update` != 'Included', n.e := N_SI]
dt[`2024\nmethods update` != 'Included', cd4_mid := cd4_mid_SI]
dt[`2024\nmethods update` != 'Included', lower := lower_SI]
dt[`2024\nmethods update` != 'Included', upper := upper_SI]
dt[,type := paste0(cd4_factor, ifelse(tt == 'peri','_1', '_2'))]

dt <- dt[,.(id, studlab, author, study_year, model, type, cd4_factor, tt, cd4_mid, time_on_art_median, on_art, predictions, lower, upper, n.e, `Default`, `2019 SR`, `2024\nmethods update`, `2024\nSR update`, NEW_2024)]

get_model_predictions <- function(ext = ''){
  est <- readRDS('./public_results/model_output/spectrum_estimates_CI.RDS')
  est <- data.table(melt(est))
  est[,tt := type]
  map <- fread('data/public_data.csv')
  map <- map[model2_pvt_sc != '',.(model2_pvt_sc, cd4_factor)] %>% unique()
  est[variable == 100, type := ifelse(type == 'peri', 'A_1', 'A_2')]
  est[variable == 275, type := ifelse(type == 'peri', 'B_1', 'B_2')]
  est[variable == 500, type := ifelse(type == 'peri', 'C_1', 'C_2')]
  
  est <- merge(est, map, by.x = 'variable', by.y = 'model2_pvt_sc', all.x = T)
  est[!is.na(cd4_factor), type := ifelse(type == 'peri', paste0(cd4_factor, '_', 1), paste0(cd4_factor, '_', 2))]
  est[variable == 2, type := 'E_1']
  est[variable == 20, type := 'F_1']
  est[variable == 40, type := 'G_1']
  est[variable == 'onart', type := 'G_2']
  est[variable == 'startart', type := 'F_2']
  est <- dcast(est[,.(variable, type, tt, model, variable.1, value)], 
               variable + type + tt +  model ~ variable.1, value.var = 'value')
  mod <- est[,.(variable, type, tt, predictions = (median), 
                lower = (lower), upper = (upper), model)]

  if(ext == ''){
    mod[,studlab := 'Modelled']
  }else{
    mod[,studlab := 'Default']
  }

  mod[,label := studlab]
  mod[model == 'model1', cd4_mid := variable]
  mod[model == 'model3', time_on_art_median := variable]
  mod[model == 'model4' & variable == 'onart', time_on_art_median := 40]
  mod[model == 'model4' & variable == 'startart', time_on_art_median := 20]
  mod[,variable := NULL]
  
  
  return(mod)
}
mod <- get_model_predictions()
pdt <- rbind(dt, mod, fill = T)

pdt[grep('D_', type),level_1 := 'Infection']
pdt[grep('H_', type),level_1 := 'Option B']
pdt[grep('I_', type),level_1 := 'Option A']
pdt[grep('J_', type),level_1 := 'SDNVP']
pdt[grep('K_', type),level_1 := 'Dual ARV']
pdt[grep('L_', type),level_1 := 'SDNVP, >350']
pdt[grep('M_', type),level_1 := 'SDNVP, <350']


##model1 levels
pdt[cd4_mid < 200 & model == 'model1', level_1 := 'CD4 midpoint [0-200)']
pdt[type %in% c('A_1', 'A_2'), level_1 := 'CD4 midpoint [0-200)']
pdt[cd4_mid > 199 & cd4_mid < 350& model == 'model1', level_1 := 'CD4 midpoint [200-350]']
pdt[type %in% c('B_1', 'B_2'), level_1 := 'CD4 midpoint [200-350]']
pdt[cd4_mid >= 350& model == 'model1', level_1 := 'CD4 midpoint >350']
pdt[type %in% c('C_1', 'C_2')& model == 'model1', level_1 := 'CD4 midpoint >350']

##model3 levels
pdt[time_on_art_median < 5, level_1 := 'ART initiated in last month']
pdt[time_on_art_median < 40 & time_on_art_median > 4, level_1 := 'ART initiated before last month']
pdt[type == 'F_1', level_1 :='ART initiated before last month']
pdt[time_on_art_median == 40, level_1 := 'ART initiated preconception']
pdt[type == 'G_1', level_1 := 'ART initiated preconception']
pdt[type == 'E_1' & is.na(level_1), level_1 := 'ART initiated in last month']

##model4 levels
pdt[on_art == TRUE, level_1 := 'ART initiated preconception']
pdt[type == 'G_2', level_1 := 'ART initiated preconception']
pdt[type == 'F_2', level_1 := 'ART initiated during pregnancy']

missing_authors <- pdt[is.na(author) & studlab != 'Modelled',studlab]
missing_years <- unlist(lapply(strsplit(missing_authors, split = '_'), '[[', 2))
missing_authors <- unlist(lapply(strsplit(missing_authors, split = '_'), '[[', 1))

pdt[is.na(author) & studlab != 'Modelled', author := missing_authors]
pdt[is.na(author) & studlab != 'Modelled', study_year := missing_years]

pdt[is.na(label), label := paste0(stringr::str_to_title(author), ', ', study_year)]
setnames(pdt, 'label', 'study_label')

pdt[,predictions := predictions * 100]
pdt[,lower := lower * 100]
pdt[,upper := upper * 100]

pdt[model == 'model1' & tt == 'peri',formatted_pred := paste0(sprintf("%.1f",predictions), ' (', 
                                                sprintf("%.1f",lower ), ', ',
                                                sprintf("%.1f",upper), ')')]

pdt[model == 'model1' & tt == 'bf',formatted_pred := paste0(sprintf("%.2f",predictions), ' (', 
                                                               sprintf("%.2f",lower), ', ',
                                                               sprintf("%.2f",upper), ')')]
#pdt[model == 'model1' & tt == 'bf' & is.na(n.e),formatted_pred := paste0(sprintf("%.2f",predictions))]

pdt[model == 'model2' & level_1 == 'Infection',formatted_pred := paste0(sprintf("%.1f",predictions), ' (', 
                                                               sprintf("%.1f",lower ), ', ',
                                                               sprintf("%.1f",upper), ')')]

pdt[model == 'model2' & tt == 'peri' & level_1 != 'Infection',formatted_pred := paste0(sprintf("%.1f",predictions), ' (', 
                                                sprintf("%.1f",lower), ', ',
                                                sprintf("%.1f",upper), ')')]

pdt[model == 'model2' & tt == 'bf'& level_1 != 'Infection',formatted_pred := paste0(sprintf("%.2f",predictions), ' (', 
                                                               sprintf("%.2f",lower), ', ',
                                                               sprintf("%.2f",upper), ')')]


pdt[model == 'model3',formatted_pred := paste0(sprintf("%.2f",predictions), ' (', 
                                                               sprintf("%.2f",lower), ', ',
                                                               sprintf("%.2f",upper), ')')]

pdt[model == 'model4',formatted_pred := paste0(sprintf("%.2f",predictions), ' (', 
                                                sprintf("%.2f",lower), ', ',
                                                sprintf("%.2f",upper), ')')]



pdt[`2024\nmethods update`  == 'Included',color := 'black']
pdt[`2024\nmethods update`  == '',color := 'grey']
pdt[NEW_2024 == 'Yes',color := 'red']
pdt[studlab == 'Modelled', color := 'black']
pdt[,fontface := ifelse(study_label == 'Modelled', 'bold', 'plain')]
pdt[model == 'model1',left_lab := cd4_mid]
pdt[model == 'model3',left_lab := time_on_art_median]

pdt <- pdt[,.(level_0 = model, level_1 = tt, level_2 = level_1, 
              color, fontface, 
              study_label, n.e, left_lab, 
              est = predictions, lower, upper, formatted_pred, 
              sr_default = Default, sr_2019 = `2019 SR`, sr_2024 = `2024\nmethods update`)]

pdt[is.na(sr_default),sr_default := '']
pdt[is.na(sr_2019),sr_2019 := '']
pdt[is.na(sr_2024),sr_2024 := '']

pdt[level_0 == 'model1', level_0 := 'Model 1']
pdt[level_0 == 'model2', level_0 := 'Model 2']
pdt[level_0 == 'model3', level_0 := 'Model 3']
pdt[level_0 == 'model4', level_0 := 'Model 4']

write.csv(pdt, paste0('./public_results/input_for_forest_plots_public.csv'), row.names = F)


##TO DO:
###Ensure that the grey studies that are in the supplement for model three actually should be included
##Ensure that the grey studies in this pdt for model 4 should actually be excluded
##Add on the part of the data frame that is column headers


