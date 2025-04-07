# install.packages('meta')
#install.packages('ggforestplot')
rm(list = ls())
library(data.table)
library(ggplot2)
library(stringr)
library(tidyverse)
library(tidybayes)
library(glmmTMB)
setwd('C:/Users/mwalters/OneDrive - Imperial College London/Papers/VT estimation/paper/')
table_dir <- './tables/'
figure_dir <- './figures/'
input_data_dir <- './data/'
output_data_dir <- './data/model_data/'
results_dir <- './results/'
source('./functions.R')

set.seed(925)

data <- fread(paste0(input_data_dir,'/public_data.csv'))
data[,studlab := paste0(author, '_', study_year)]
data[,id := 1:nrow(data)]

data <- data[model %in% c('model3')]

data[,LPV := ifelse(grepl('LPV', data$primary_treatment_type), 1, 0)]
data[,NVP := ifelse(grepl('NVP', data$primary_treatment_type), 1, 0)]
data[,RAL := ifelse(grepl('RAL', data$primary_treatment_type), 1, 0)]
data[,ATV := ifelse(grepl('ATV', data$primary_treatment_type), 1, 0)]
data[,EFV := ifelse(grepl('EFV', data$primary_treatment_type), 1, 0)]
data[,DTG:= ifelse(grepl('DTG', data$primary_treatment_type), 1, 0)]

data[,treatment_info := DTG+ LPV + EFV + NVP + RAL + ATV]
data[,misc_reg := ifelse(treatment_info == 0 , 1, 0)]
data[,time_on_art_median := inferred_ART_weeks]

data <- data[,.(studlab, type = tt, location, study_year, year_start, year_end, inferred_primary, event.e, n.e, time_on_art_median, vls_threshold, vls_prop, vl_time, id, DTG, LPV, ATV, EFV, NVP, RAL, misc_reg, NEW_2024)]
data <- melt(data, id.vars = c('studlab', 'type','inferred_primary' , 'location', 'study_year','year_start', 'year_end',  'event.e', 'n.e', 'time_on_art_median', 'vls_threshold', 'vls_prop',  'vl_time','id', 'NEW_2024'))

data <- data[value == 1,]
data <- data[,.(studlab, type, study_year, year_start, year_end,  inferred_primary, location, event.e, n.e, time_on_art_median, vls_threshold, vls_prop, vl_time, id, regimen = variable, NEW_2024)]

data[regimen %in% c('ATV', 'LPV'),class:='PI']
data[regimen %in% c('NVP','EFV'), class := 'NNRTI']
data[regimen %in% c('RAL', 'DTG'), class := 'INSTI']
data[regimen %in% c('misc_reg'), class := 'misc_reg']


region_map <- fread('./public_results/data/region_map.csv')

data$location <- tools::toTitleCase(data$location)
data[,location := gsub('Ivory Coast', "Côte d'Ivoire", data$location)]
data <- merge(data, region_map, by.x = 'location', by.y = 'country', all.x = T, all.y = F)
data <- unique(data)
data[location %in% c('India', 'Korea'), Region := "Asia and the Pacific"]
data[location %in% c( "Malawi and Mozambique"), Region := "Eastern and southern Africa"]
data[location %in% c("Argentina, Brazil, South Africa, Tanzania, Thailand, US",
                   "Uganda, Zimbabwe", "South Africa and Uganda" ), Region := 'Multiple regions']
data[location %in% c("Europe", 'Uk', 'Uk and Ireland', 'United States', 'Usa'), Region := "Western and central Europe and North America"]
data[Region %in% c( "Eastern and southern Africa","Western and central Africa"   ),loc := 'SSA']
data[Region %in% c('Latin America', 'Asia and the Pacific', 'Western and central Europe and North America'),loc := 'non SSA']
data[Region == "Multiple regions" ,loc := 'mixed regions']
data$loc <- factor(data$loc, levels = c('SSA', 'non SSA', 'mixed regions'))

data[,time_on_art_median := round(time_on_art_median)]
data[time_on_art_median < 7, tri := '3rd trimester']
data[time_on_art_median %in% 7:20, tri := '2nd trimester']
data[time_on_art_median %in% 21:39, tri := '1st trimester']
data[time_on_art_median == 40, tri := 'Preconception']

write.csv(data, './public_results/data/public_data_assigned_trt.csv', row.names = F)
