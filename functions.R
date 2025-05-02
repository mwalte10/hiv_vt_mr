get_posterior <- function(model){
  set.seed(925)
  cov <- vcov(model)$cond
  for(i in 1:nrow(cov)){
    for(j in 1:ncol(cov)){
      cov[i,j] <- round(cov[i,j], 8)
    }
  }
  posterior_draws <- mvtnorm::rmvnorm(3000, mean = fixef(model)$cond, sigma = cov)
  posterior_draws <- data.table(posterior_draws)
  posterior_draws[,draw := 1:3000]
  posterior_draws <- melt(posterior_draws, id.vars = c('draw'))
  posterior_draws <- dcast(posterior_draws, draw ~ variable)
  return(posterior_draws)
}

get_peri_tr <- function(x, cd4){
  cd4 <- cd4 / 100 - 5
  cd4_mid <- x[['cd4_mid']]
  ttperi <- x[['(Intercept)']]
  int <- 0
  
  
  tr <- cd4_mid * cd4 + ttperi + int * cd4
  return(tr)
}

get_bf_tr <- function(x, cd4){
  cd4 <- cd4 / 100 - 5
  cd4_mid <- x[['cd4_mid']]
  ttbf <- x[['ttbf']] + x[['(Intercept)']]
  int <- x[['cd4_mid:ttbf']]
  
  tr <- cd4_mid * cd4 + ttbf  + int * cd4 
  return(tr)
}

get_model1_estimates <- function(posterior, cd4){
  peri <- get_peri_tr(x = posterior, cd4 = cd4)
  bf <- get_bf_tr(x = posterior, cd4 = cd4)
  
  
  return(data.table(cd4 = cd4, peri = peri, bf = bf, draw = 1:length(peri)))
}


get_model3_estimates <- function(posterior, week){
  week_val = week - 20 
  if(length(week) > 1){
    add <- c(rep(posterior[['late_start']],4), rep(0,37))
  }else if(week < 5){
    add <- posterior[['late_start']]
  }else{
    add = 0
  }
  peri <- posterior[['time_on_art_median']] * week_val + posterior[['(Intercept)']] + add
  
  dt <- data.table(value = peri, week = week, draw = 1:length(peri))
  
  return(dt)
}


get_model4_estimates <- function(posterior){
  
  dt <- data.table(onart = posterior[['(Intercept)']],
                   startart = posterior[['(Intercept)']] + posterior[['on_artFALSE']])
  dt[,draw := 1:nrow(dt)]
  
  return(dt)
}

prep_fp_diamonds_mod1 <- function(){
  ####model label
  mod <- readRDS('./model_output/spectrum_estimates_CI.RDS')
  mod <- mod[model == 'model1']
  mod[,studlab := 'Modelled: 2024']
  mod[,label := studlab]
  mod[variable == 100, level_2 := 'CD4 midpoint [0-200)']
  mod[variable == 275, level_2 := 'CD4 midpoint [200-350]']
  mod[variable == 500, level_2 := 'CD4 midpoint >350']
  setnames(mod, 'type', 'tt')
  
  mod <- mod[,.(level_0 = 'Model 1', level_1 = tt, level_2, study_label = studlab,
                est = median, lower, upper, formatted_pred, sr_default = '', sr_2019 = '', sr_2024 = 'Estimate')]
  dia <- rbind(dia[,.(level_0, level_1, level_2, study_label, est, lower, upper, formatted_pred, sr_default = 'Estimate', sr_2019 = '', sr_2024 = '')],
               mod[,.(level_0, level_1, level_2, study_label, est, lower, upper, formatted_pred, sr_default, sr_2019, sr_2024)])
  
  model_label = 'Model 1'
  data <- fread(paste0(input_data_dir,'/public_data.csv'))
  data[,studlab := paste0(author, '_', study_year)]
  data[,id := 1:nrow(data)]
  data[, type := paste0(cd4_factor, ifelse(tt == 'peri', '_1', '_2'))]
  
  mod_data <- data[model == 'model1']
  mod_data[,cd4_mid :=(cd4_mid - 500)/100]
  mod_data$tt <- factor(mod_data$tt, levels = c('peri', 'bf'))
  mod_data[,est := weighted.mean(event.e/n.e, n.e), by = 'type']
  mod_data[,sum_n := sum(n.e), by = 'type']
  mod_data <- unique(mod_data[,.(cd4_factor, tt, type, est, sum_n)])
  mod_data <- mod_data[,.(level_0 = 'Model 1', level_1 = tt, level_2 = cd4_factor, est = est * 100, sum_n)]
  mod_data[level_2 == 'A', level_2 := 'CD4 midpoint [0-200)']
  mod_data[level_2 == 'B', level_2 := 'CD4 midpoint [200-350]']
  mod_data[level_2 == 'C', level_2 := 'CD4 midpoint >350']
  mod_data[,lower := (est/100 - ((1.96)/sqrt(sum_n)) * sqrt(est/100 * (1-est/100)))*100]
  mod_data[,upper := (est/100 + ((1.96)/sqrt(sum_n)) * sqrt(est/100 * (1-est/100)))*100]
  mod_data[lower < 0 , lower:=0]
  mod_data[level_1 == 'peri',formatted_pred := paste0(sprintf("%.1f",est), ' (', 
                                                      sprintf("%.1f",lower ), ', ',
                                                      sprintf("%.1f",upper), ')')]
  
  mod_data[level_1 == 'bf',formatted_pred := paste0(sprintf("%.2f",est), ' (', 
                                                    sprintf("%.2f",lower), ', ',
                                                    sprintf("%.2f",upper), ')')]
  mod_data[,study_label := 'Weighted average: 2024']
  mod_data <- mod_data[,.(level_0, level_1, level_2, study_label, est, lower, upper, formatted_pred, sr_default = '', sr_2019 = '', sr_2024 = 'Estimate')]
  mod_data <- rbind(mod_data, data.table(level_0 = model_label, level_1 = 'bf', level_2 = 'CD4 midpoint [200-350]',
                                         study_label = 'Weighted average: 2024', est = '', lower = '', upper = '',
                                         formatted_pred = 'No studies', sr_default = '', sr_2019 = '', sr_2024 = ''))
  
  dia <- rbind(dia, dt)
  
  dia$study_label <- factor(dia$study_label, levels = c('Median: default', 
                                                        'Weighted average: default',
                                                        'Weighted average: 2024',
                                                        'Modelled: 2024'))
  dia[,y := as.integer(fct_rev(study_label))]
  
  expand_dim <- data.table(name = c('lwr', 'prop', 'upr', 'prop2'), 
                           y_scalar = c(0,0.2,0,-0.2),
                           y = rep(seq(1:4), each = 4))
  dia <- merge(dia, expand_dim, by = c('y'), allow.cartesian = T)
  dia[,y := y + y_scalar]
  dia[name %in% c('prop', 'prop2'),x := est]
  dia[name %in% c('lwr'),x := lower]
  dia[name %in% c('upr'),x := upper]
  
  dia <- dia[,.(level_0,level_1, level_2, study_label, left_lab = '', name, y, x, formatted_pred, sr_default,
                sr_2019, sr_2024)]
  dia[level_2 == 'CD4 midpoint [0-200)' & study_label == 'Modelled', left_lab := 100]
  dia[level_2 == 'CD4 midpoint [200-350]'& study_label == 'Modelled', left_lab := 275]
  dia[level_2 == 'CD4 midpoint >350'& study_label == 'Modelled', left_lab := 500]
}

prep_files <- function(pjnz){
  source('C:/Users/mwalters/frogger//scripts/read_spectrum.R')
  dp <- eppasm::read_dp(pjnz)
  dp.x = dp
  demp <- prepare_leapfrog_demp(pjnz)
  proj <- prepare_leapfrog_projp(pjnz)
  proj <- prepare_hc_leapfrog_projp(pjnz, proj, dp.x = dp, demp = demp)
  country <- eppasm::read_country(pjnz)
  
  demp$netmigr <- leapfrog:::read_netmigr(pjnz, adjust_u5mig = FALSE)
  demp$netmigr_adj <- leapfrog:::adjust_spectrum_netmigr(demp$netmigr)
  
  ## projection parameters
  dpfile <- grep(".DP$", utils::unzip(pjnz, list=TRUE)$Name, value=TRUE)
  
  yr_start <- as.integer(dpsub(dp = dp.x, "<FirstYear MV2>",2,4))
  yr_end <- as.integer(dpsub(dp = dp.x, "<FinalYear MV2>",2,4))
  proj.years <- yr_start:yr_end
  timedat.idx <- 4+1:length(proj.years)-1
  year.idx <- 1:length(proj.years)
  
  proj$ctx_effect <- 0.33
  proj$laf <- 1
  proj$paed_art_elig_age <- as.integer(proj$paed_art_elig_age)
  proj$mat_prev_input <- rep(TRUE, length(year.idx))
  pmtct_new <- array(0, dim = c(7, length(year.idx)), dimnames = list(pmtct = c("Option A", "Option B", "SDNVP", "Dual ARV", "Option B+: before pregnancy", "Option B+: >4 weeks", "Option B+: <4 weeks")))
  ## pick out which ones were inserted as numbers
  pmtct_new[, which(colSums(proj$pmtct)[, 1] > 0)] <- proj$pmtct[, (which(colSums(proj$pmtct)[, 1] > 0)), 1]
  ## pick out which ones were inserted as percent
  pmtct_new[, which(colSums(proj$pmtct)[, 1] == 0)] <- proj$pmtct[, which(colSums(proj$pmtct)[, 1] == 0), 2]
  proj$pmtct <- pmtct_new
  
  saveRDS(proj, paste0(spec_file_dir, '/proj/', eppasm::read_country(pjnz), '.RDS'))
  saveRDS(demp, paste0(spec_file_dir, '/demp/', eppasm::read_country(pjnz), '.RDS'))
  
}

prep_spectrum_format <- function(){
  spec_estimates <- readRDS('./model_output/spectrum_estimates.RDS')
  spec_estimates <- spec_estimates[,.(type, variable, value = plogis(median))]
  
  mtct <- array(data = NA, dim = c(8,2), dimnames = list(cd4 = c('>500', '350-500', '250-349', '200-249', '100-199', '50-99', '<50', 'INFECTION'),
                                                         trans_type = c('peri', 'bf')))
  pmtct_mtct <- array(data= NA, dim = c(7,7,2), dimnames = list(cd4 = c('>500', '350-500', '250-349', '200-249', '100-199', '50-99', '<50'),
                                                                pmtct_reg = c('option A', 'option B', 'single dose nevirapine', 'WHO 2006 dual ARV regimen',
                                                                              'ART before pregnancy', 'ART >4 weeks before delivery', 'ART <4 weeks before delivery'),
                                                                transission_type = c('peri', 'bf')))
  
  mtct[c('100-199', '50-99', '<50'),'peri'] <- spec_estimates[type == 'peri' & variable == 100, value] 
  mtct[c('200-249', '250-349'),'peri'] <- spec_estimates[type == 'peri' & variable == 275, value] 
  mtct[c('350-500', '>500'),'peri'] <- spec_estimates[type == 'peri' & variable == 500, value] 
  
  mtct[c('100-199', '50-99', '<50'),'bf'] <- spec_estimates[type == 'bf' & variable == 100, value] 
  mtct[c('200-249', '250-349'),'bf'] <- spec_estimates[type == 'bf' & variable == 275, value] 
  mtct[c('350-500', '>500'),'bf'] <- spec_estimates[type == 'bf' & variable == 500, value]
  
  mtct['INFECTION', 'peri'] <- spec_estimates[type == 'peri' & variable == 'mat_sero', value] 
  mtct['INFECTION', 'bf'] <- spec_estimates[type == 'bf' & variable == 'mat_sero', value] 
  
  pmtct_mtct[,'option A','peri'] <- spec_estimates[type == 'peri' & variable == 'opt_a', value] 
  pmtct_mtct[,'option A','bf'] <- spec_estimates[type == 'bf' & variable == 'opt_a', value] 
  
  pmtct_mtct[,'option B','peri'] <- spec_estimates[type == 'peri' & variable == 'opt_b', value] 
  pmtct_mtct[,'option B','bf'] <- spec_estimates[type == 'bf' & variable == 'opt_b', value] 
  
  pmtct_mtct[,'single dose nevirapine','peri'] <- spec_estimates[type == 'peri' & variable == 'sdnvp', value] 
  pmtct_mtct[c('>500', '350-500'),'single dose nevirapine','bf'] <- spec_estimates[type == 'bf' & variable == 'sdnvp_gte350', value] 
  pmtct_mtct[c('250-349', '200-249', '100-199', '50-99', '<50'),'single dose nevirapine','bf'] <- spec_estimates[type == 'bf' & variable == 'sdnvp_lte350', value] 
  
  pmtct_mtct[,'WHO 2006 dual ARV regimen','peri'] <- spec_estimates[type == 'peri' & variable == 'dual_arv', value]
  pmtct_mtct[,'WHO 2006 dual ARV regimen','bf'] <- spec_estimates[type == 'bf' & variable == 'dual_arv', value]
  
  pmtct_mtct[,'ART before pregnancy','peri'] <- spec_estimates[type == 'peri' & variable == 40,value]
  pmtct_mtct[,'ART >4 weeks before delivery','peri'] <- spec_estimates[type == 'peri' & variable == 20,value]
  pmtct_mtct[,'ART <4 weeks before delivery','peri'] <- spec_estimates[type == 'peri' & variable == 2,value]
  
  pmtct_mtct[,'ART before pregnancy','bf'] <- spec_estimates[type == 'bf' & variable == 'onart',value]
  pmtct_mtct[,'ART >4 weeks before delivery','bf'] <- spec_estimates[type == 'bf' & variable == 'startart',value]
  pmtct_mtct[,'ART <4 weeks before delivery','bf'] <- spec_estimates[type == 'bf' & variable == 'startart',value]
  
  saveRDS(list(pmtct_mtct = pmtct_mtct, mtct = mtct), './spectrum_files/updated_parms.RDS')
}

get_stacked_bar_outputs <- function(pjnz){
  parameters <- readRDS(paste0('./spectrum_files//proj/', eppasm::read_country(pjnz), '.RDS'))
  demp <- readRDS(paste0('./spectrum_files/demp/', eppasm::read_country(pjnz), '.RDS'))
  country <- eppasm::read_country(pjnz)
  
  out <- run_model(demp, parameters, 1970:2023, NULL, run_child_model = TRUE)
  
  strat <- out$paed_inf_strat
  dimnames(strat) <- list(
    pmtct = c('On ART pre-conception', 'Started ART during pregnancy',
              'On ART, dropped out', 'Started ART, dropped out',
              'PVT, non-ART', 
              'No PVT',
              'Maternal seroconversion'),
    tt = c('perinatal', 'bf0-6', 'bf7-12', 'bf12+', 'bf24+'),
    year = 1970:2023)
  
  an_inf <- apply(out$p_infections[1:10,,], 3, sum)
  an_inf <- data.table(year = 1970:2023, inf_total = an_inf)
  
  strat <- data.table(melt(strat))
  cat_map <- data.table(pmtct = c('On ART pre-conception', 'Started ART during pregnancy',
                                  'On ART, dropped out', 'Started ART, dropped out',
                                  'PVT, non-ART', 
                                  'No PVT',
                                  'Maternal seroconversion'),
                        broad_pmtct = c('On ART pre-conception', 'Started ART during pregnancy',
                                        'On ART, dropped out', 'Started ART, dropped out',
                                        'PVT, non-ART', 
                                        'No PVT',
                                        'Maternal seroconversion'))
  
  strat <- merge(strat, cat_map, by = 'pmtct')
  strat[,tt := ifelse(tt == 'perinatal', 'perinatal', 'bf')]
  strat <- strat[,.(value = sum(value)), by = c('tt', 'year', 'broad_pmtct')]
  strat[,total_inf := sum(value), by = 'year']
  strat$broad_pmtct <- factor(strat$broad_pmtct, levels = c('Maternal seroconversion', 'No PVT', 'Started ART, dropped out',
                                                            'On ART, dropped out', 'PVT, non-ART',
                                                            'Started ART during pregnancy','On ART pre-conception'))
  strat[,run := 'Former VT']
  default <- copy(strat)
  
  new_parms <- readRDS(paste0('./spectrum_files/updated_parms.RDS'))
  parameters$mtct <- new_parms$mtct
  parameters$pmtct_mtct <-  new_parms$pmtct_mtct
  
  out <- run_model(demp, parameters, 1970:2023, NULL, run_child_model = TRUE)
  strat <- out$paed_inf_strat
  dimnames(strat) <- list(
    pmtct = c('On ART pre-conception', 'Started ART during pregnancy',
              'On ART, dropped out', 'Started ART, dropped out',
              'PVT, non-ART', 
              'No PVT',
              'Maternal seroconversion'),
    tt = c('perinatal', 'bf0-6', 'bf7-12', 'bf12+', 'bf24+'),
    year = 1970:2023)
  
  strat <- data.table(melt(strat))
  
  strat <- merge(strat, cat_map, by = 'pmtct')
  strat[,tt := ifelse(tt == 'perinatal', 'perinatal', 'bf')]
  strat <- strat[,.(value = sum(value)), by = c('tt', 'year', 'broad_pmtct')]
  strat[,total_inf := sum(value), by = 'year']
  strat$broad_pmtct <- factor(strat$broad_pmtct, levels = c('Maternal seroconversion', 'No PVT', 'Started ART, dropped out',
                                                            'On ART, dropped out', 'PVT, non-ART',
                                                            'Started ART during pregnancy','On ART pre-conception'))
  strat[,run := 'Updated VT']
  
  strat <- rbind(strat, default)
  strat <- merge(strat, an_inf, by = 'year')
  
  write.csv(strat, paste0('./spectrum_files/spectrum_format/inf_diff/', country,'.csv'),row.names = F)
  
}

plot_country_results <- function(country, year.x){
  strat <- fread(paste0('./spectrum_files/spectrum_format/inf_diff/', country,'.csv'))
  
  dt <- strat[year == year.x]
  dt[value < 0 ,value := 0]
  dt[,prop := value / total_inf]
  
  dt[, tt := ifelse(tt == 'perinatal', 'Perinatal', 'Breastfeeding')]
  dt$tt <- factor(dt$tt, levels = c('Perinatal', 'Breastfeeding'))
  
  dt[,max := sum(value), by = c('tt', 'run')]
  
  pct_diff <- unique(dt[,.(tt, run, max)])
  pct_diff <- dcast(pct_diff, tt ~ run)
  copy_pd <- copy(pct_diff)
  copy_pd <- copy_pd[,.(`Former VT` = sum(`Former VT`), `Updated VT` = sum(`Updated VT`))]
  copy_pd[,tt := 'total']
  pct_diff <- rbind(copy_pd, pct_diff)
  pct_diff[,pct_diff := ((`Updated VT` - `Former VT`) / `Former VT`) * 100]
  pct_diff[,color := ifelse(pct_diff < 0, 'Red', 'Green')]
  pct_diff[,pct_diff := ifelse(pct_diff <0, sprintf("%.1f",pct_diff), paste0('+',sprintf("%.1f",pct_diff)))]
  pct_diff[,x := 'Updated VT']
  pct_diff <- merge(pct_diff, unique(dt[run == 'Updated VT',.(tt, max)]), by = 'tt', all.x = T)
  pct_diff[,y := ifelse(max < 1000, max*1.08 ,max *1.05)]
  
  color_map <- unique(dt[,.(tt, broad_pmtct)])
  color_map[broad_pmtct == 'No PVT',color := 'deepskyblue4']
  color_map[broad_pmtct == 'Started ART, dropped out',color := 'deepskyblue']
  color_map[broad_pmtct == 'On ART, dropped out',color := 'darkslategray1']
  color_map[broad_pmtct == 'Maternal seroconversion',color := 'darkorange']
  color_map[broad_pmtct == 'PVT, non-ART',color := 'deeppink4']
  color_map[broad_pmtct == 'Started ART during pregnancy',color := 'darkolivegreen4']
  color_map[broad_pmtct == 'On ART pre-conception',color := 'darkgreen']
  
  color_map$broad_pmtct <- factor(color_map$broad_pmtct, levels = c('Maternal seroconversion', 'No PVT',
                                                                    'Started ART, dropped out',
                                                                    'On ART, dropped out',
                                                                    'PVT, non-ART',
                                                                    'Started ART during pregnancy', 'On ART pre-conception'))
  
  color_map <- color_map[order(broad_pmtct),]
  
  
  gg <- ggplot()  +
    geom_bar(data = dt, aes(run, value, fill = broad_pmtct),position="stack", stat="identity") +
    scale_fill_manual(values =color_map$color,
                      breaks = color_map$broad_pmtct) +
    theme_bw(base_size = 8)+
    theme(
      strip.text = element_text(size = rel(1.3)),
      axis.text = element_text(size = rel(1.2)),
      axis.title = element_text(size = rel(1.3)),
      plot.caption = element_text(size = rel(1.1), color = ifelse(pct_diff[tt == 'total',color] == 'Green', 'darkgreen', 'red'))) +
    labs(fill = 'Transmission category', title = paste0(country, ', ', year.x)) +
    guides(fill = guide_legend(title.position = 'top')) +
    labs(x = NULL, y = 'Paediatric infections',
         caption = paste0('Total number of vertical transmissions differs by ', pct_diff[tt == 'total',pct_diff], '%')) +
    facet_wrap(~tt) + coord_cartesian(ylim = c(0,max(dt$max) * 1.1)) +
    geom_text(data = pct_diff[tt != 'total'], aes(x = x, y = y, label = paste0(pct_diff, '%'), col = color), show.legend = F,
              size = 5) +
    scale_color_manual(values = c('darkgreen', 'red'), breaks = c('Green', 'Red'))
  gg
  
  
  
  png(paste0('./figures_tables/appendix/spec_stacked_bar/', country, '_', year.x, '.png'), width = 6, height = 3.5, units = 'in', res = 900)
  print(gg)
  dev.off()
}

format_spectrum_mtct <- function(spec_estimates){
  mtct_format <- array(data = NA, dim = c(8,2), dimnames = list(cd4 = c(">500", "350-500", 
                                                                        "250-349", "200-249",
                                                                        "100-199",   "50-99",    
                                                                        "<50", "INFECTION"),
                                                                trans_type = c('perinatal', 'bf')))
  
  mtct_format[c('>500', '350-500'),'perinatal'] <- spec_estimates[type == 'peri' & variable == 500, median]
  mtct_format[c('250-349', '200-249'),'perinatal'] <- spec_estimates[type == 'peri' & variable == 275, median]
  mtct_format[c('100-199', '50-99', '<50'),'perinatal'] <- spec_estimates[type == 'peri' & variable == 100, median]
  mtct_format[c('INFECTION'),'perinatal'] <- spec_estimates[type == 'peri' & variable == 'mat_sero', median]
  
  mtct_format[c('>500', '350-500'),'bf'] <- spec_estimates[type != 'peri' & variable == 500, median]
  mtct_format[c('250-349', '200-249'),'bf'] <- spec_estimates[type != 'peri' & variable == 275, median]
  mtct_format[c('100-199', '50-99', '<50'),'bf'] <- spec_estimates[type != 'peri' & variable == 100, median]
  mtct_format[c('INFECTION'),'bf'] <- spec_estimates[type != 'peri' & variable == 'mat_sero', median]
  
  pmtct_format <- array(data = NA, dim = c(7,7,2), dimnames = list(cd4 = c(">500", "350-500", 
                                                                           "250-349", "200-249",
                                                                           "100-199",   "50-99",    
                                                                           "<50"),
                                                                   pmtct_reg = c("option A", "option B",                    
                                                                                 "single dose nevirapine",
                                                                                 "WHO 2006 dual ARV regimen",
                                                                                 "ART before pregnancy",
                                                                                 "ART >4 weeks before delivery",
                                                                                 "ART <4 weeks before delivery"),
                                                                   transmission_type = c('perinatal', 'bf')))
  pmtct_format[,'option A','perinatal'] <- spec_estimates[type == 'peri' & variable == 'opt_a',median]
  pmtct_format[,'option A','bf'] <- spec_estimates[type != 'peri' & variable == 'opt_a',median]
  pmtct_format[,'option B','perinatal'] <- spec_estimates[type == 'peri' & variable == 'opt_b',median]
  pmtct_format[,'option B','bf'] <- spec_estimates[type != 'peri' & variable == 'opt_b',median]
  pmtct_format[,'WHO 2006 dual ARV regimen','perinatal'] <- spec_estimates[type == 'peri' & variable == 'dual_arv',median]
  pmtct_format[,'WHO 2006 dual ARV regimen','bf'] <- spec_estimates[type != 'peri' & variable == 'dual_arv',median]
  pmtct_format[,'single dose nevirapine','perinatal'] <- spec_estimates[type == 'peri' & variable == 'sdnvp',median]
  pmtct_format[c('>500', '350-500'),'single dose nevirapine','bf'] <- spec_estimates[type != 'peri' & variable == 'sdnvp_gte350',median]
  pmtct_format[c('250-349', '200-249', '100-199', '50-99', '<50'),'single dose nevirapine','bf'] <- spec_estimates[type != 'peri' & variable == 'sdnvp_lte350',median]
  pmtct_format[,'ART before pregnancy','perinatal'] <- spec_estimates[type == 'peri' & variable == '40',median]
  pmtct_format[,'ART before pregnancy','bf'] <- spec_estimates[type != 'peri' & variable == 'onart',median]
  pmtct_format[,'ART >4 weeks before delivery','perinatal'] <- spec_estimates[type == 'peri' & variable == '20',median]
  pmtct_format[,'ART >4 weeks before delivery','bf'] <- spec_estimates[type != 'peri' & variable == 'startart',median]
  pmtct_format[,'ART <4 weeks before delivery','perinatal'] <- spec_estimates[type == 'peri' & variable == '2',median]
  pmtct_format[,'ART <4 weeks before delivery','bf'] <- spec_estimates[type != 'peri' & variable == 'startart',median]
  
  return(list(mtct = mtct_format, pmtct_mtct = pmtct_format))
  
}

get_diamond <- function(){
  dt <- readRDS('./model_output/estimate_draws.RDS')
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
  
  modelled_default <- readRDS('./model_output/old_estimate_draws.RDS')
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
  map <- fread('./naming_map.csv')
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
  
  saveRDS(dia, paste0('./figures_tables/appendix/diamonds.RDS'))
  
}

get_model_predictions <- function(ext = ''){
  est <- readRDS('./model_output/spectrum_estimates_CI.RDS')
  est <- data.table(melt(est))
  est[,tt := type]
  map <- fread('./data/public_data.csv')
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

get_predictions <- function(model_name, study_data){
  if(model_name == 'model1'){
    model <- readRDS(paste0('./model_output/model_1.RDS'))
  }
  if(model_name == 'model2'){
    model <- readRDS(paste0('./model_output/model_2.RDS'))
  }
  if(model_name == 'model3'){
    model <- readRDS(paste0('./model_output/model_3.RDS'))
  }
  if(model_name == 'model4'){
    model <- readRDS(paste0('./model_output/model_4.RDS'))
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

prep_forest_plot_data <- function(){
  data <- fread(paste0('./data/public_data.csv'))
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
  
  mod1 <- get_predictions(model_name = 'model1', study_data = data[model == 'model1'])
  mod2 <- get_predictions(model_name = 'model2', study_data = data[model == 'model2'])
  mod3 <- get_predictions(model_name = 'model3', study_data = data[model == 'model3'])
  mod4 <- get_predictions(model_name = 'model4', study_data = data[model == 'model4'])
  dt <- rbind(mod1, mod2, mod3, mod4)
  dt[!is.na(cd4_mid), cd4_mid := (cd4_mid + 5) * 100]
  dt[!is.na(time_on_art_median), time_on_art_median := (time_on_art_median + 20)]
  
  stud_include <- fread(paste0('./data/studies_included_by_analysis_public.csv'))
  stud_include <- stud_include[,.(studlab = paste0(author, '_', study_year), 
                                  app_ref = `Appendix reference`,
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
  
  dt <- dt[,.(id, studlab, author, study_year, app_ref, model, type, cd4_factor, tt, cd4_mid, time_on_art_median, on_art, predictions, lower, upper, n.e, `Default`, `2019 SR`, `2024\nmethods update`, `2024\nSR update`, NEW_2024)]
  
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
  
  pdt[is.na(label), label := paste0("'",stringr::str_to_title(author), ', ', study_year, "'^", app_ref)]
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
                study_label, n.e = round(n.e), left_lab, 
                est = predictions, lower, upper, formatted_pred, 
                sr_default = Default, sr_2019 = `2019 SR`, sr_2024 = `2024\nmethods update`)]
  
  pdt[is.na(sr_default),sr_default := '']
  pdt[is.na(sr_2019),sr_2019 := '']
  pdt[is.na(sr_2024),sr_2024 := '']
  
  pdt[level_0 == 'model1', level_0 := 'Model 1']
  pdt[level_0 == 'model2', level_0 := 'Model 2']
  pdt[level_0 == 'model3', level_0 := 'Model 3']
  pdt[level_0 == 'model4', level_0 := 'Model 4']
  
  write.csv(pdt, paste0('./input_for_forest_plots_public.csv'), row.names = F)
  
}

get_level_2_plot <- function(dt_in, level_2_in, diamond = diamond_in){
  dt2 <- dt_in[level_2 == level_2_in]
  if(level_2_in == 'CD4 midpoint [200-350]' & nrow(dt2) == 0){
    tt <- 'bf'
    model <- 'Model 1'
  }else{
    tt <- unique(dt2[,level_1])
    model <-unique(dt2[,level_0]) 
  }
  
  #fp_parms <- fread(paste0(figure_dir, '/forest_plot_space_parms.csv'))
  space_parms <- fp_parms[level_0 == model & level_1 == tt & level_2 == level_2_in,]
  
  
  header <- data.table(study_label = 'Study', sr_default = 'Default',
                       sr_2019 = '2019',
                       sr_2024 = '2024',
                       formatted_pred = 'VT probability (%)',
                       n.e = 'N',
                       color = 'black', 
                       fontface = 'bold')
  
  if(model == 'Model 1'){
    header[,left_lab := 'CD4']
  }
  if(model == 'Model 3'){
    header[,left_lab := 'Weeks']
  }
  dt2 <- rbind(header, dt2, fill = T)
  
  diamond <- diamond[level_2 == level_2_in & level_1 == tt,]
  diamond[,x := as.numeric(x)]
  diamond[study_label == 'Modelled: 2024',formatted_pred := dt2[study_label == 'Modelled', formatted_pred]]
  diamond[study_label == 'Modelled: 2024' & name == 'upr',x := dt2[study_label == 'Modelled', upper]]
  diamond[study_label == 'Modelled: 2024' & name == 'lwr',x := dt2[study_label == 'Modelled', lower]]
  diamond[study_label == 'Modelled: 2024' & name %in% c('prop', 'prop2'),x := dt2[study_label == 'Modelled', est]]
  
  
  if(model %in% c('Model 1', 'Model 3', 'Model 4')){
    order <- setdiff(dt2[order(left_lab), study_label], c('Modelled', 'Study'))
    order <- c('Study',order, 'Modelled')
  }else if(model == 'Model 2'){
    order <- setdiff(dt[order(study_label), study_label], c('Modelled', 'Study'))
    order <- c('Study',order, 'Modelled')
  } 
  dt2[,save := study_label]
  
  if(model == 'Model 3' & level_2_in == 'ART initiated before last month'){
    dt2[,study_label := paste0(study_label, '_', left_lab, '_', 1:nrow(dt2))]
    order <- setdiff(dt2[order(as.integer(left_lab)), study_label], c('Modelled_20_71', 'Study_Weeks_1'))
    order <- c('Study_Weeks_1', order, 'Modelled_20_71')
  }
  
  dt2$study_label <- factor(dt2$study_label, levels = order)
  dt2[save != 'Study', left_lab := round(as.integer(left_lab))]
  
  if(all(sort(unique(dt2$color)) == c('black', 'grey', 'red'))){
    cols <- c('black', 'grey', 'red')
  }else if(all(sort(unique(dt2$color)) == c('black', 'grey'))){
    cols <- c('black', 'grey')
  } else if(all(sort(unique(dt2$color)) == c('black', 'red'))){
    cols <- c('black', 'red')
  } else if(all(sort(unique(dt2$color)) == c('black'))){
    cols <- 'black'
  }
  
  dt2$color <- factor(dt2$color, cols)
  
  
  dt2 <- dt2[study_label != 'Modelled',]
  dt2 <- dt2[save != 'Modelled',]
  
  #####Spacing parameters
  xmax <- space_parms$xmax
  pc_left_x <- c(0,xmax)
  pc_right_x <- c(0,xmax)
  pc_mid_x <- c(0,xmax)
  pc_mid_y <- c(1,nrow(dt2))
  left_x <- c(0, 1.6)
  right_x <- c(0, 2.3)
  pc_height = 0.1 *  nrow(dt2)
  blh <- 0.05/pc_height
  title_h <- 0.1/pc_height
  glh <- 0.1/pc_height
  mh <- 0.2/pc_height
  fontsize <- 7.5
  if(model == 'Model 2'){
    text_size = 7.5 / .pt
    right_text_size = 5/.pt
  }else{
    text_size = 7.5 / .pt
    right_text_size = 5/.pt 
  }
  if(model == 'Model 3'){
    text_size = 6 / .pt
    right_text_size = 4.5/.pt
  }
  
  left_lab_x <- space_parms$left_lab_scalar * xmax
  n_x <- space_parms$n_scalar * xmax
  left_lab_x <- 1.5
  n_x <- 1.25
  ds <- space_parms$ds
  aspect <- diff(pc_left_x) / nrow(dt2)
  
  
  main <- ggplot(data = dt2, aes(y = fct_rev(study_label))) +
    geom_linerange(aes(xmin = lower, xmax = upper, col = color)) +
    geom_point(aes(x =est, col = color), size = 2, shape = 15, show.legend = F) +
    coord_cartesian(xlim = pc_mid_x) +
    theme_void() +
    theme(plot.margin = unit(c(t=0,r=0,b=0.2,l=0), "cm"),
          legend.position = 'none') + labs(x = 'VT probability (%)') +
    scale_color_manual(values = cols) 
  
  left <- ggplot(data = dt2, aes(y = fct_rev(study_label))) +
    geom_text(aes(x = 0,  label = save, col = color), hjust = 0, fontface ='bold', show.legend = F, size = text_size, parse = T) +
    geom_text(aes(x = n_x,  label = n.e, fontface = fontface, col = color), hjust = 0.5, show.legend = F, size = text_size) +
    geom_text(aes(x = left_lab_x,  label = left_lab, fontface = fontface, col = color), hjust = 0.5, show.legend = F, size = text_size) +
    scale_color_manual(values = cols)  +
    theme_void(base_size = fontsize) +
    theme(plot.margin = unit(c(t=0,r=0,b=0.2,l=0), "cm")) +
    coord_cartesian(xlim =left_x)
  
  right <-  ggplot(data = dt2, aes(y = fct_rev(study_label))) +
    geom_text(aes(x = 0 , label = formatted_pred, fontface = fontface, col = color), hjust = 0, show.legend = F, size = text_size) +
    geom_text(aes(x = 1.4, label = sr_default, fontface = fontface, col = color), hjust = 0.5, show.legend = F, size = right_text_size) +
    geom_text(aes(x = 1.8 , label = sr_2019, fontface = fontface, col = color), hjust = 0.5, show.legend = F, size = right_text_size) +
    geom_text(aes(x = 2.2 , label = sr_2024 , fontface = fontface, col = color), hjust = 0.5, show.legend = F, size = right_text_size) +
    theme_void() + coord_cartesian(xlim = right_x) +
    scale_color_manual(values = cols) +
    theme(plot.margin = unit(c(t=0,r=0,b=0.2,l=0), "cm"))
  
  study_plot <- (left + main + right) +
    plot_layout(widths = c(1.3, 1, 1))  
  
  diamond_left <- ggplot(data = unique(diamond[,.(study_label, left_lab)]), aes(y = fct_rev(study_label))) +
    geom_text(aes(x = 0 ,  label = study_label), hjust = 0, fontface ='bold', show.legend = F, size = text_size) +
    geom_text(aes(x = left_lab_x , label = left_lab), hjust = 0.5, show.legend = F, size = text_size)  +
    theme_void(base_size = fontsize) +
    theme(plot.margin = unit(c(t=0,r=0,b=0.2,l=0), "cm")) +
    coord_cartesian(xlim = left_x)
  
  
  if(model %in% c('Model 3')){
    breaks_dia = c(1,2)
    ylim_dia =  c(0.55,3.45)
  }
  if(model%in% c('Model 1', 'Model 2', 'Model 4')){
    breaks_dia = c(1,3)
    ylim_dia =  c(0.55,4.45)
  }
  
  if(any(diamond$y > 5)){
    diamond[y >2.15, y:=y -1]
  }
  needs_point <- diamond[study_label %in% c(diamond[is.na(x),study_label])]
  
  if(max(diamond$y)-0.15 > length(unique(diamond$study_label))){
    diamond[y >2.15, y:=y -1]
  }
  
  diamond_main <-  ggplot() +
    geom_polygon(data = diamond, aes(x = x, y = y, group = study_label), 
                 fill = 'gray', color = 'black', size = 0.4) +
    geom_point(data = needs_point[name == 'prop'], aes(x = x, y = y - 0.15), size = 2, pch = 18) + 
    coord_cartesian(xlim =pc_mid_x, ylim = ylim_dia) +
    scale_y_continuous(breaks = breaks_dia) +
    theme_classic(base_size = fontsize) +
    theme(axis.line.y = element_blank(),
          axis.ticks.y= element_blank(),
          axis.text.y= element_blank(),
          axis.title.y= element_blank(),
          axis.text.x = element_text(size = rel(1.1)),
          axis.title.x = element_text(size = rel(1.1), margin = margin(b = 5)), 
          plot.margin = unit(c(t=0,r=0,b=.2,l=0), "cm")) + labs(x = 'VT probability (%)')  
  
  diamond_right <- ggplot(data = unique(diamond[,.(study_label, formatted_pred, sr_default, sr_2019, sr_2024)]), aes(y = fct_rev(study_label))) +
    geom_text(aes(x = 0,label = formatted_pred), fontface = 'bold', hjust = 0, show.legend = F, size = text_size) +
    geom_text(aes(x = 1.4, label = sr_default), hjust = 0.5, show.legend = F, size = right_text_size) +
    geom_text(aes(x = 1.8 , label = sr_2019), hjust = 0.5, show.legend = F, size = right_text_size) +
    geom_text(aes(x = 2.2 , label = sr_2024), hjust = 0.5, show.legend = F, size = right_text_size) +
    theme_classic(base_size = fontsize) + coord_cartesian(xlim = right_x) +
    theme(plot.margin = unit(c(t=0,r=0.3,b=.2,l=0), "cm"),
          axis.text = element_blank(),
          axis.line = element_blank(),
          axis.title = element_blank(),
          axis.ticks = element_blank()) 
  
  diamond_plot  <- (diamond_left + diamond_main + diamond_right) +
    plot_layout(widths = c(1.3, 1, 1))  
  
  
  line <- ggplot() +
    geom_segment(aes(x = 0, xend = 5, y = 1.5, yend = 1.5), 
                 linetype = "solid", colour = "gray45", size = 0.15) +
    theme_void() +
    theme(plot.margin = unit(c(t=0,r=0,b=0,l=0), "cm")) +
    coord_cartesian(ylim = c(0, 2), xlim = c(0.14, 4.78))
  
  if(model == 'Model 2'){
    level_2_in = ''
  }
  if(model == 'Model 1') {
    if(level_2_in == 'CD4 midpoint [0-200)'){
      level_2_in <- paste0('A. ', 'CD4 midpoint [0-200)')
    }
    if(level_2_in == 'CD4 midpoint [200-350]'){
      level_2_in <- paste0('B. ', 'CD4 midpoint [200-350]')
    }
    if(level_2_in == 'CD4 midpoint >350'){
      level_2_in <- paste0('C. ', 'CD4 midpoint >350')
    }
  }
  
  if(model == 'Model 3') {
    if(level_2_in == 'ART initiated in last month'){
      level_2_in <- paste0('A. ', 'ART initiated in last month')
    }
    if(level_2_in == 'ART initiated before last month'){
      level_2_in <- paste0('B. ', 'ART initiated before last month')
    }
    if(level_2_in == 'ART initiated preconception'){
      level_2_in <- paste0('C. ', 'ART initiated preconception')
    }
  }
  
  if(model == 'Model 4') {
    if(level_2_in ==  "ART initiated during pregnancy"){
      level_2_in <- paste0('A. ',  "ART initiated during pregnancy")
    }
    
    if(level_2_in == 'ART initiated preconception'){
      level_2_in <- paste0('B. ', 'ART initiated preconception')
    }
  }
  height = nrow(dt2) #+ length(unique(diamond$study_label)) 
  
  plot_level2 <-(study_plot / line / diamond_plot) +
    plot_annotation(title = level_2_in) +
    plot_layout(heights = c(height / 3, 0.1, 1)) &
    theme(plot.margin = unit(c(t=0,r=0,b=0,l=5), "pt"),
          plot.title = element_text(size = fontsize * 1.1, face="bold", hjust = 0, #-0.014
                                    margin = margin(t = 0, r = 0, b = 0, l = 0)),
          plot.title.position = "plot")
  
  plot_level2 
  return(plot_level2)
}

plot_fp <- function(pdt, level_0_in, level_1_in){
  dt <- pdt[level_0 == level_0_in & level_1 == level_1_in, ]
  dt_in = dt
  
  if(level_0_in == 'Model 1'){
   # fp_parms <- fread(paste0(figure_dir, '/forest_plot_space_parms.csv'))
    space_parms <- fp_parms[level_0 == level_0_in  & level_1 == level_1_in,]
    
    cd4_lte200 <- get_level_2_plot(dt_in = dt, level_2_in = 'CD4 midpoint [0-200)')
    cd4_mid <- get_level_2_plot(dt_in = dt, level_2_in = 'CD4 midpoint [200-350]') ##Peri is missing all 2019
    cd4_gte350 <- get_level_2_plot(dt_in = dt, level_2_in = 'CD4 midpoint >350')
    
    heights_tt = c(nrow(dt[level_2 == 'CD4 midpoint [0-200)'])+6,
                   nrow(dt[level_2 == 'CD4 midpoint [200-350]'])+6,
                   nrow(dt[level_2 == 'CD4 midpoint >350'])+6)
    heights_tt[heights_tt == 5] <- 6
    height_plot <- 7.5 * (sum(heights_tt) / 46)
    
    # unique(space_parms$height_plot)
    
    heights_tt = heights_tt / sum(heights_tt)
    full_title <- unique(space_parms$plot_title)
    
    
    plot <- ggarrange(cd4_lte200,cd4_mid, cd4_gte350,  ncol = 1,
                      heights = heights_tt)
    # plot
    
    title  <- ggplot() +
      labs(title = full_title) + 
      theme_void() +
      theme(plot.margin = unit(c(t=0,r=0,b=0.15,l=0.35), "cm"),
            plot.title = element_text(face= 'bold', size = 12))
    
    plot_full <- ggarrange(title, plot,  ncol = 1,
                           heights = c(0.04,1))
    
    png(paste0("./figures_tables/appendix/forest_plots/", level_0_in, '_', level_1_in, '.png'),
        height = height_plot, width = 7.5, res = 900, units = 'in')
    print(plot_full)
    dev.off()
    
    
  }
  
  if(level_0_in == 'Model 2'){
    for(level_2 in unique(dt$level_2)){
      #fp_parms <- fread(paste0(figure_dir, '/forest_plot_space_parms.csv'))
      level_2_in = level_2
      space_parms <- fp_parms[level_0 == level_0_in  & level_1 == level_1_in & level_2 == level_2_in,]
      
      plot  <- get_level_2_plot(dt_in = dt, level_2_in = level_2)
      lev_2 = level_2
      height_plot <- 7.5 * ((nrow(dt[level_2 == lev_2])+5) / 46)
      if(level_2 %in% c("SDNVP, >350", "SDNVP, <350", "Dual ARV", 'Option A')){
        height_plot = height_plot * 1.1
      }
      title_fac <- 7.5 / height_plot
      
      full_title <- unique(space_parms$plot_title)
      # plot
      
      title  <- ggplot() +
        labs(title = full_title) + 
        theme_void() +
        theme(plot.margin = unit(c(t=0,r=0,b=0,l=.2), "cm"),
              plot.title = element_text(face= 'bold', size = 12))
      
      plot_full <- ggarrange(title, plot,  ncol = 1,
                             heights = c(0.04 * title_fac,1))
      
      if(level_2_in == "SDNVP, >350"){
        level_2_in = 'sdnvp_high'
      }
      if(level_2_in == "SDNVP, <350"){
        level_2_in = 'sdnvp_low'
      }
      
      png(paste0("./figures_tables/appendix/forest_plots/", level_0_in, '_', level_1_in, '_', level_2_in,'.png'),
          height = height_plot, width = 7.5, res = 900, units = 'in')
      print(plot_full)
      dev.off()
    }
  }
  
  if(level_0_in == 'Model 3'){
   # fp_parms <- fread(paste0(figure_dir, '/forest_plot_space_parms.csv'))
    space_parms <- fp_parms[level_0 == level_0_in  & level_1 == level_1_in,]
    
    ft <- get_level_2_plot(dt_in = dt, level_2_in = 'ART initiated in last month')
    st <- get_level_2_plot(dt_in = dt, level_2_in = 'ART initiated before last month') ##Peri is missing all 2019
    pc <- get_level_2_plot(dt_in = dt, level_2_in = 'ART initiated preconception')
    
    heights_tt = c(nrow(dt[level_2 == 'ART initiated in last month'])+6,
                   nrow(dt[level_2 == 'ART initiated before last month'])+6,
                   nrow(dt[level_2 == 'ART initiated preconception'])+6)
    heights_tt[heights_tt == 5] <- 6
    height_plot <- 7.5 * (sum(heights_tt) / 46)
    full_title <- unique(space_parms$plot_title)    
    heights_tt = heights_tt / sum(heights_tt)
    heights_tt[1] <- heights_tt[1] *1.1
    
    plot <- ggarrange(ft, st, pc,  ncol = 1,
                      heights = heights_tt)
    # plot
    
    title  <- ggplot() +
      labs(title = full_title) + 
      theme_void() +
      theme(plot.margin = unit(c(t=0,r=0,b=0,l=.2), "cm"),
            plot.title = element_text(face= 'bold', size = 12))
    
    plot_full <- ggarrange(title, plot,  ncol = 1,
                           heights = c(0.02,1))
    
    if(height_plot > 11.5){
      height_plot = 11.5
    }
    
    png(paste0("./figures_tables/appendix/forest_plots/", level_0_in, '_', level_1_in, '.png'),
        height = height_plot , width = 7.5, res = 900, units = 'in')
    print(plot_full)
    dev.off()
    
    
  }
  
  if(level_0_in == 'Model 4'){
   # fp_parms <- fread(paste0(figure_dir, '/forest_plot_space_parms.csv'))
    space_parms <- fp_parms[level_0 == level_0_in  & level_1 == level_1_in,]
    
    startart <- get_level_2_plot(dt_in = dt, level_2_in = "ART initiated during pregnancy" )
    pc <- get_level_2_plot(dt_in = dt, level_2_in = 'ART initiated preconception')
    
    heights_tt = c(nrow(dt[level_2 == "ART initiated during pregnancy" ])+6,
                   nrow(dt[level_2 == 'ART initiated preconception'])+6)
    heights_tt[heights_tt == 5] <- 6
    height_plot <- 7.5 * (sum(heights_tt) / 46)
    full_title <- unique(space_parms$plot_title)    
    heights_tt = heights_tt / sum(heights_tt)
    heights_tt[1] <- heights_tt[1] *1.1
    
    plot <- ggarrange(startart, pc,  ncol = 1,
                      heights = heights_tt)
    # plot
    
    title  <- ggplot() +
      labs(title = full_title) + 
      theme_void() +
      theme(plot.margin = unit(c(t=0,r=0,b=0,l=.2), "cm"),
            plot.title = element_text(face= 'bold', size = 12))
    
    plot_full <- ggarrange(title, plot,  ncol = 1,
                           heights = c(0.05,1))
    
    png(paste0("./figures_tables/appendix/forest_plots/", level_0_in, '_', level_1_in, '.png'),
        height = height_plot, width = 7.5, res = 900, units = 'in')
    print(plot_full)
    dev.off()
    
    
  }
  
}
