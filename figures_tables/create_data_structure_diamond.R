get_diamond <- function(){
  dt <- readRDS('./public_results/model_output/estimate_draws.RDS')
  dt <- dt[variable == 100,cd4_factor := 'A']
  dt <- dt[variable == 275,cd4_factor := 'B']
  dt <- dt[variable == 500,cd4_factor := 'C']
  dt <- dt[variable == 'mat_sero',cd4_factor := 'D']
  dt <- dt[variable == 'dual_arv',cd4_factor := 'K']
  dt <- dt[variable == 'opt_a',cd4_factor := 'I']
  dt <- dt[variable == 'opt_b',cd4_factor := 'H']
  dt <- dt[variable == 'sdnvp',cd4_factor := 'J']
  dt <- dt[variable == 'sdnvp_gte350',cd4_factor := 'L']
  dt <- dt[variable == 'sdnvp_lte350',cd4_factor := 'M']
  dt <- dt[variable == 2 & model == 'model3',cd4_factor := 'E']
  dt <- dt[variable == 20& model == 'model3',cd4_factor := 'F']
  dt <- dt[variable == 40& model == 'model3',cd4_factor := 'G']
  dt <- dt[variable == 'onart',cd4_factor := 'G']
  dt <- dt[variable == 'startart',cd4_factor := 'F']
  setnames(dt, 'type', 'tt')
  dt[,study_label := 'Modelled: Updated']
  dt[,type := paste0(cd4_factor, '_', ifelse(tt == 'peri', 1, 2))]
  dt <- dt[!is.na(cd4_factor)]
  dt[,median := quantile(value, 0.5), by = 'type']
  dt[,lower := quantile(value, 0.025), by = 'type']
  dt[,upper := quantile(value, 0.975), by = 'type']
  dt <- dt[!(model == 'model1' & cd4_factor == 'E')]
  dt <- dt[!(model == 'model1' & cd4_factor == 'F')]
  dt <- dt[!(model == 'model1' & cd4_factor == 'G')]
  

  data <- fread(paste0(input_data_dir,'/public_data.csv'))
  data[,type := paste0(cd4_factor, '_', ifelse(tt == 'peri', 1, 2))]
  data_wa <- data[,.(model, cd4_factor, tt, type, med = weighted.mean(prop_infected, n.e), sum_n = sum(n.e)), by = 'type']
  data_wa <- unique(data_wa)
  confint <- DescTools::BinomCI(data_wa$med * data_wa$sum_n, n = data_wa$sum_n,  method = 'wilson')
  data_wa[,lower := as.vector(confint[,2])]
  data_wa[,upper:=  as.vector(confint[,3])]
  data_wa[lower < 0 ,lower := 0]
  data_wa[,study_label := 'Weighted average: Updated']
  
  data_default <- fread(paste0(input_data_dir, '/default_values_public.csv'))
  data_default[,tt := ifelse(tt == 'peri', 1, 2)]
  data_default <- data_default[,.(model,
                              med = weighted.mean(event.e / n.e, n.e), 
                              sum_n = sum(n.e), type = paste0(cd4_factor, '_', tt)), by = c('tt', 'cd4_factor')]
  data_default <- unique(data_default[!is.na(tt),.(model, tt, type, cd4_factor, med, sum_n)])
  data_default[,study_label := 'Weighted average: Former']
  confint <- DescTools::BinomCI(data_default$med * data_default$sum_n, n = data_default$sum_n,  method = 'wilson')
  data_default[,lower := as.vector(confint[,2])]
  data_default[,upper:= as.vector(confint[,3])]
  data_default_wa <- data_default
  
  ##for model 1 we do median default rather than the modelled default
  data_default <- fread(paste0(input_data_dir, '/default_values_public.csv'))
  data_default[,tt := ifelse(tt == 'peri', 1, 2)]
  data_default <- data_default[model== 'model1']
  data_default[,prop := event.e / n.e]
  data_default[is.na(n.e), prop := event.e] ## fix for c2
  data_default <- unique(data_default[,.(model, median = median(prop)), by = c('tt', 'cd4_factor')])
  data_default[,type := paste0(cd4_factor, '_', tt)]
  
  modelled_default <- readRDS('./public_results/model_output/old_estimate_draws.RDS')
  modelled_default <- modelled_default[variable == 100,cd4_factor := 'A']
  modelled_default <- modelled_default[variable == 275,cd4_factor := 'B']
  modelled_default <- modelled_default[variable == 500,cd4_factor := 'C']
  modelled_default <- modelled_default[variable == 'mat_sero',cd4_factor := 'D']
  modelled_default <- modelled_default[variable == 'dual_arv',cd4_factor := 'K']
  modelled_default <- modelled_default[variable == 'opt_a',cd4_factor := 'I']
  modelled_default <- modelled_default[variable == 'opt_b',cd4_factor := 'H']
  modelled_default <- modelled_default[variable == 'sdnvp',cd4_factor := 'J']
  modelled_default <- modelled_default[variable == 'sdnvp_gte350',cd4_factor := 'L']
  modelled_default <- modelled_default[variable == 'sdnvp_lte350',cd4_factor := 'M']
  modelled_default <- modelled_default[variable == 2 & model == 'model3',cd4_factor := 'E']
  modelled_default <- modelled_default[variable == 20 & model == 'model3',cd4_factor := 'F']
  modelled_default <- modelled_default[variable == 40 & model == 'model3',cd4_factor := 'G']
  modelled_default <- modelled_default[variable == 'onart',cd4_factor := 'G']
  modelled_default <- modelled_default[variable == 'startart',cd4_factor := 'F']
  setnames(modelled_default, 'type', 'tt')
  modelled_default <- modelled_default[!is.na(cd4_factor),]
  modelled_default[,study_label := 'Modelled: Former']
  modelled_default[,type := paste0(cd4_factor, '_', ifelse(tt == 'peri', 1, 2))]
  modelled_default[,median := quantile(value, 0.5), by = 'type']
  modelled_default[,lower := quantile(value, 0.025), by = 'type']
  modelled_default[,upper := quantile(value, 0.975), by = 'type']
  modelled_default <- modelled_default[!(model == 'model1' & cd4_factor == 'E')]
  modelled_default <- modelled_default[!(model == 'model1' & cd4_factor == 'F')]
  modelled_default <- modelled_default[!(model == 'model1' & cd4_factor == 'G')]
  ##we don't put in these old modelled estimates because the data hadn't been extracted this way
  modelled_default <- modelled_default[model %in%  c('model2', 'model4'),]
  
  dt <- rbind(dt[,.(model, tt, cd4_factor, med = plogis(median), 
                    lower = plogis(lower), upper = plogis(upper), study_label, type)], 
              data_default_wa[,.(model, tt, cd4_factor, med, lower, upper, study_label, type)], 
              data_default[,.(model = 'model1', tt, cd4_factor, med = median, 
                              lower = NA, upper = NA, study_label = 'Median: Former', type)], 
              data_wa[,.(model, tt, cd4_factor,  med, upper, lower, study_label, type)],
              modelled_default[,.(model, tt, cd4_factor,  med = plogis(median), 
                                  lower = plogis(lower), upper = plogis(upper),  study_label, type)])
  dt <- dt[!is.na(med)]
  dt <- unique(dt)
  
  diamond_test <- dt
  diamond_test[tt == 1, tt:='peri']
  diamond_test[tt == 2, tt := 'bf']
  diamond_test <- diamond_test[,.(level_0 = model, level_1 = tt, level_2 = cd4_factor, study_label, med, lower, upper)]
  map <- fread('./public_results/naming_map.csv')
  map[,level_2 := unlist(lapply(unlist(lapply(map$type, strsplit, split = '_'), recursive = F), '[[', 1))]
  diamond <- merge(map, diamond_test, by.x = c('tt', 'level_2'), by.y = c('level_1', 'level_2'), all.y = T)
  
  diamond <- diamond[,.(level_0, level_1 = tt, level_2 = name, study_label,
                        med = med * 100, lower = lower * 100, upper = upper* 100)]
  diamond[lower < 0, lower := 0]
  diamond[level_2 == 'Infection',formatted_pred := paste0(sprintf("%.1f",med), ' (', 
                                                          sprintf("%.1f",lower ), ', ',
                                                          sprintf("%.1f",upper), ')')]
  
  diamond[ level_1 == 'peri' & level_2 != 'Infection',formatted_pred := paste0(sprintf("%.1f",med), ' (', 
                                                                               sprintf("%.1f",lower), ', ',
                                                                               sprintf("%.1f",upper), ')')]
  
  diamond[level_1 == 'bf'& level_2 != 'Infection',formatted_pred := paste0(sprintf("%.2f",med), ' (', 
                                                                           sprintf("%.2f",lower), ', ',
                                                                           sprintf("%.2f",upper), ')')]
  diamond[level_0 == 'model1' & is.na(lower)  & level_1 =='bf', formatted_pred := sprintf("%.2f",med)]
  diamond[level_0 == 'model1' & is.na(lower) & level_1 =='peri', formatted_pred := sprintf("%.1f",med)]
  
  
  diamond <- rbind(diamond, data.table(
    level_0 = 'model1', level_1 = 'bf', 
    level_2 = c('CD4 midpoint >350', "CD4 midpoint [200-350]", 
                "CD4 midpoint [200-350]", "CD4 midpoint [200-350]"),
    study_label = c('Weighted average: Former',
                    'Weighted average: Former',
                    'Median: Former',
                    'Weighted average: Updated'), 
    med = '', lower  = '', upper = '',
    formatted_pred = c('Missing N',rep('No studies',3))))
  
  diamond$study_label <- factor(diamond$study_label, levels = c('Median: Former',
                                                                'Weighted average: Former', 
                                                                'Modelled: Former', 
                                                                'Weighted average: Updated', 
                                                                'Modelled: Updated'))
  
  diamond[,y := as.integer(fct_rev(study_label))]
  
  expand_dim <- data.table(name = c('lwr', 'prop', 'upr', 'prop2'), 
                           y_scalar = c(0,0.15,0,-0.15),
                           y = rep(seq(1:5), each = 4))
  
  
  dia <- merge(diamond, expand_dim, by = c('y'), allow.cartesian = T)
  dia[,y := y + y_scalar]
  dia[name %in% c('prop', 'prop2'),x := med]
  dia[name %in% c('lwr'),x := lower]
  dia[name %in% c('upr'),x := upper]
  dia <- dia[,.(level_0,level_1, level_2, study_label, left_lab = '', name, y, x, formatted_pred, sr_default = '',
                sr_2019 = '', sr_2024 = '')]
  dia[study_label %in% c('Modelled: Former', 'Weighted average: Former','Median: Former'), sr_default := 'Estimate']
  dia[study_label %in% c('Modelled: Updated', 'Weighted average: Updated'), sr_2024 := 'Estimate']
  
  dia[level_2 == 'CD4 midpoint [0-200)' & study_label == 'Modelled', left_lab := 100]
  dia[level_2 == 'CD4 midpoint [200-350]'& study_label == 'Modelled', left_lab := 275]
  dia[level_2 == 'CD4 midpoint >350'& study_label == 'Modelled', left_lab := 500]
  
  dia[,id := paste0(level_0, level_1, level_2)]
  
  saveRDS(dia, paste0('./public_results/figures_tables/appendix/diamonds.RDS'))
  
}


