# install.packages('meta')
#install.packages('ggforestplot')
rm(list = ls())
library(data.table)
library(ggplot2)
library(stringr)
library(tidyverse)
library(tidybayes)
library(glmmTMB)
input_data_dir <- './data/'
source('./functions.R')

set.seed(925)

data <- fread(paste0(input_data_dir,'/public_data.csv'))
data[,studlab := paste0(author, '_', study_year)]
data[,id := 1:nrow(data)]

###############################################################################
##Models, new data
###############################################################################
mod1_data <- data[model == 'model1']
mod1_data[,cd4_mid :=(cd4_mid - 500)/100]
mod1_data$tt <- factor(mod1_data$tt, levels = c('peri', 'bf'))
form_binom_mod1 <- cbind(event.e,(n.e-event.e)) ~ cd4_mid  + tt + cd4_mid * tt + (1 | studlab) + (1| id)
model1  <- glmmTMB(form_binom_mod1, data = mod1_data, family = binomial())
saveRDS(model1, './model_output/model_1.RDS')
summary(model1)
mod1_posterior <- get_posterior(model = model1)
mod1_result_draws <- rbindlist(lapply(c(100,275,500), get_model1_estimates, posterior = mod1_posterior))
mod1_result_draws_spec <- mod1_result_draws 
mod1_result_summary <- melt(mod1_result_draws, id.vars = c('cd4', 'draw'))
setnames(mod1_result_summary, c('variable', 'cd4'), c('type', 'variable'))
mod1_result_summary <- mod1_result_summary[,.(median = quantile(value, 0.5), 
                                              lower = quantile(value, 0.025),
                                              upper = quantile(value, 0.975)), by = c('type', 'variable')]
mod1_result_draws <- rbindlist(lapply(0:700, get_model1_estimates, posterior = mod1_posterior))
mod1_result_draws <- melt(mod1_result_draws, id.vars = c('cd4', 'draw'))
setnames(mod1_result_draws, c('cd4', 'variable'), c('variable', 'type'))



mod2_data <- data[model == 'model2']
mod2_data[,pvt_type := paste0(model2_pvt_sc, '_', tt)]
form_binom_mod2 <- cbind(event.e,(n.e-event.e)) ~ -1 + pvt_type + (1 | studlab)  + (1|id) 
model2 <- glmmTMB(form_binom_mod2, data = mod2_data, family = binomial())
saveRDS(model2, './model_output/model_2.RDS')
summary(model2)
mod2_posterior <- get_posterior(model = model2)
mod2_result_draws <- melt(mod2_posterior, id.vars = c('draw'))
mod2_result_draws[,type := ifelse(grepl('bf', mod2_result_draws$variable), 'bf', 'peri')]
mod2_result_draws[,variable := gsub('pvt_type', '', variable)]
mod2_result_draws[,variable := gsub('_bf', '', variable)]
mod2_result_draws[,variable := gsub('_peri', '', variable)]
mod2_result_draws_spec <- mod2_result_draws
mod2_result_summary <- mod2_result_draws[,.(median = quantile(value, 0.5), 
                                            lower = quantile(value, 0.025),
                                            upper = quantile(value, 0.975)), by = c('type', 'variable')]

mod3_data <- data[model == 'model3']
mod3_data[,time_on_art_median := inferred_ART_weeks - 20]
mod3_data[,late_start := ifelse(inferred_ART_weeks < 4, 1, 0)]
form_binom_mod3 <- cbind(event.e,(n.e-event.e)) ~  time_on_art_median + late_start + (1 | studlab) + (1 | id)
model3 <- glmmTMB(form_binom_mod3, data = mod3_data, family = binomial())
saveRDS(model3, './model_output/model_3.RDS')
summary(model3)
mod3_posterior <- get_posterior(model = model3)
mod3_result_draws <- rbindlist(lapply(c(2,20,40), get_model3_estimates, posterior = mod3_posterior))
mod3_result_draws <- mod3_result_draws[,.(variable = week, type = 'peri', draw, value)]
mod3_result_draws_spec <- mod3_result_draws
mod3_result_summary <- mod3_result_draws[,.(median = quantile(value, 0.5), 
                                            lower = quantile(value, 0.025),
                                            upper = quantile(value, 0.975)), by = c('type', 'variable')]
mod3_result_draws <- rbindlist(lapply(0:40, get_model3_estimates, posterior = mod3_posterior))
mod3_result_draws <- mod3_result_draws[,.(draw, type = 'peri', variable = week, value)]


mod4_data <- data[model == 'model4']
mod4_data[,event.e := event.e *10]
mod4_data[,n.e := n.e *10]
mod4_data[,on_art := ifelse(inferred_ART_weeks == 40,T,F)]
mod4_data$on_art <- factor(mod4_data$on_art, levels = c(T, F))
form_binom_mod4 <- cbind(event.e,(n.e-event.e)) ~ on_art + (1 | id) #+ (1 | studlab)
model4 <- glmmTMB(form_binom_mod4, data = mod4_data, family = binomial)
saveRDS(model4, './model_output/model_4.RDS')
summary(model4)
mod4_posterior <- get_posterior(model = model4)
mod4_result_draws <- get_model4_estimates(posterior = mod4_posterior)
mod4_result_draws <- melt(mod4_result_draws, id.vars = 'draw')
mod4_result_draws[,type := 'bf']
mod4_result_draws_spec <- mod4_result_draws
mod4_result_summary <- mod4_result_draws[,.(median = quantile(value, 0.5), 
                                            lower = quantile(value, 0.025),
                                            upper = quantile(value, 0.975)), by = c('type', 'variable')]

###############################################################################
##Models, pre 2018 data
###############################################################################
old_data <- data[NEW_2024 == 'N',]
mod1_old_data <- old_data[model == 'model1']
mod1_old_data[,cd4_mid :=(cd4_mid - 500)/100]
mod1_old_data$tt <- factor(mod1_old_data$tt, levels = c('peri', 'bf'))
model1_old  <- glmmTMB(form_binom_mod1, data = mod1_old_data, family = binomial())
summary(model1_old)
mod1_old_posterior <- get_posterior(model = model1_old)
mod1_old_result_draws <- rbindlist(lapply(0:700, get_model1_estimates, posterior = mod1_old_posterior))
mod1_old_result_draws <- melt(mod1_old_result_draws, id.vars = c('cd4', 'draw'))
setnames(mod1_old_result_draws, c('cd4', 'variable'), c('variable', 'type'))

mod2_old_data <- old_data[model == 'model2']
mod2_old_data[,pvt_type := paste0(model2_pvt_sc, '_', tt)]
mod2_old_data[,pvt_type := paste0(model2_pvt_sc, '_', tt)]
model2_old <- glmmTMB(form_binom_mod2, data = mod2_old_data, family = binomial())
summary(model2_old)
mod2_old_posterior <- get_posterior(model = model2_old )
mod2_old_result_draws <- melt(mod2_old_posterior, id.vars = c('draw'))
mod2_old_result_draws[,type := ifelse(grepl('bf', mod2_old_result_draws$variable), 'bf', 'peri')]
mod2_old_result_draws[,variable := gsub('pvt_type', '', variable)]
mod2_old_result_draws[,variable := gsub('_bf', '', variable)]
mod2_old_result_draws[,variable := gsub('_peri', '', variable)]

mod3_old_data <- old_data[model == 'model3']
mod3_old_data[,time_on_art_median := inferred_ART_weeks - 20]
mod3_old_data[,late_start := ifelse(inferred_ART_weeks < 4, 1, 0)]
model3_old <- glmmTMB(form_binom_mod3, data = mod3_old_data, family = binomial())
summary(model3_old)
mod3_old_posterior <- get_posterior(model = model3_old)
mod3_old_result_draws <- rbindlist(lapply(0:40, get_model3_estimates, posterior = mod3_old_posterior))
mod3_old_result_draws <- mod3_old_result_draws[,.(draw, type = 'peri', variable = week, value)]

mod4_old_data <- old_data[model == 'model4']
mod4_old_data[,event.e := event.e *10]
mod4_old_data[,n.e := n.e *10]
mod4_old_data[,on_art := ifelse(inferred_ART_weeks == 40,T,F)]
mod4_old_data$on_art <- factor(mod4_old_data$on_art, levels = c(T, F))
model4_old <- glmmTMB(form_binom_mod4, data = mod4_old_data, family = binomial)
summary(model4_old)
mod4_old_posterior <- get_posterior(model = model4_old)
mod4_old_result_draws <- get_model4_estimates(posterior = mod4_old_posterior)
mod4_old_result_draws <- melt(mod4_old_result_draws, id.vars = 'draw')
mod4_old_result_draws[,type := 'bf']


###############################################################################
##Save results
###############################################################################
draws <- rbind(mod1_result_draws[,model := 'model1'],
               mod2_result_draws[,model := 'model2'],
               mod3_result_draws[,model := 'model3'],
               mod4_result_draws[,model := 'model4'])
saveRDS(draws, './model_output/estimate_draws.RDS')

draws <- rbind(mod1_old_result_draws[,model := 'model1'],
               mod2_old_result_draws[,model := 'model2'],
               mod3_old_result_draws[,model := 'model3'],
               mod4_old_result_draws[,model := 'model4'])
saveRDS(draws, './model_output/old_estimate_draws.RDS')

vt_parms_spec_draws <- rbind(draws[model == 'model1' & variable %in% c(100,275,500),],
                             draws[model == 'model2',],
                             draws[model == 'model3' & variable %in% c(2,20,40),],
                             draws[model == 'model4',])
vt_parms_spec_draws[,value := plogis(value)]
saveRDS(vt_parms_spec_draws, './model_output/spec_draws.RDS')

spec_estimates <- rbind(mod1_result_summary[,model := 'model1'],
                        mod2_result_summary[,model := 'model2'],
                        mod3_result_summary[,model := 'model3'],
                        mod4_result_summary[,model := 'model4'])
spec_estimates <- spec_estimates[,.(type, variable, model, median = plogis(median), 
                                    lower = plogis(lower), 
                                    upper = plogis(upper))]
spec_estimates[type == 'peri',formatted_pred := paste0(round(median*100, 1), ' (', round(lower*100 , 1), ', ',round(upper*100,1), ')')]
spec_estimates[type != 'peri',formatted_pred := paste0(round(median*100, 2), ' (', round(lower*100 , 2), ', ',round(upper*100,2), ')')]
saveRDS(spec_estimates, './model_output/spectrum_estimates_CI.RDS')
spec_estimates <- format_spectrum_mtct(spec_estimates)
saveRDS(spec_estimates, './model_output/spectrum_estimates.RDS')

posterior_estimates <- rbind(melt(mod1_posterior, id.vars = 'draw')[,model := 'model1'],
                             melt(mod2_posterior, id.vars = 'draw')[,model := 'model2'],
                             melt(mod3_posterior, id.vars = 'draw')[,model := 'model3'],
                             melt(mod4_posterior, id.vars = 'draw')[,model := 'model4'])
saveRDS(posterior_estimates, './model_output/posterior_estimates.RDS')


model_table <- posterior_estimates
model_table <- model_table[,.(median = median(value),
                              lower = quantile(value, 0.025),
                              upper = quantile(value, 0.975)), by = c('variable', 'model')]
model_table <- melt(model_table, id.vars = c('variable', 'model'))
model_table[,odds := exp(value)]
model_table[,value := round(value, 2)]
model_table[,odds := round(odds, 3)]
model_table <- melt(model_table, id.vars = c('variable', 'model', 'variable.1'))
model_table <- dcast(model_table, variable + model + variable.2 ~ variable.1, value.var = 'value')
model_table <- model_table[,value := paste0(median, ' (', lower, ', ', upper, ')')]
model_table <- dcast(model_table[,.(variable, model, variable.2, value)], model + variable ~ variable.2, value.var = 'value')
saveRDS(model_table, './model_output/regression_table.RDS')

################################################################################
##Treatment type analysis
################################################################################
data <- fread('./data/public_data_assigned_trt.csv')
data[,time_on_art_median := time_on_art_median - 20]
data[,late_start := ifelse((time_on_art_median + 20) < 4, 1, 0)]
data$class <- factor(data$class, levels = c('NNRTI', 'INSTI', 'PI', 'misc_reg'))
data$loc <- factor(data$loc, levels = c('SSA', 'non SSA', 'mixed regions'))
data$tri <- factor(data$tri, levels = c('Preconception', '1st trimester', '2nd trimester', '3rd trimester'))

form_binom_mod3_art <- cbind(event.e,(n.e-event.e)) ~  time_on_art_median + late_start + class + (1 | studlab) + (1 | id)
model3_art <- glmmTMB(form_binom_mod3_art, data = data, family = binomial())
summary(model3_art)
posterior_mod3_art <- get_posterior(model3_art)
saveRDS(posterior_mod3_art, './model_output/posterior_draws_mod3_art.RDS')

form_binom_mod3_art_region <- cbind(event.e,(n.e-event.e)) ~  time_on_art_median + late_start + class + loc + (1 | studlab) + (1 | id)
model3_art_region <- glmmTMB(form_binom_mod3_art_region, data = data, family = binomial())
summary(model3_art_region)
posterior_mod3_art_region <- get_posterior(model3_art_region)
saveRDS(posterior_mod3_art_region, './model_output/posterior_draws_mod3_art_region.RDS')

################################################################################
##Viral load suppression analysis
################################################################################
data <- fread('./data/public_data_assigned_trt.csv')
data[,time_on_art_median := time_on_art_median - 20]
data[,late_start := ifelse((time_on_art_median + 20) < 4, 1, 0)]
data$class <- factor(data$class, levels = c('NNRTI', 'INSTI', 'PI', 'misc_reg'))
data$loc <- factor(data$loc, levels = c('SSA', 'non SSA', 'mixed regions'))
data$tri <- factor(data$tri, levels = c('Preconception', '1st trimester', '2nd trimester', '3rd trimester'))
data_vls <- data[vls_threshold == 50 & vl_time %in% c('delivery', 'third trimester') & type == 'peri']
data_vls[,n_vls := (vls_prop * n.e)]
data_vls[,time := ifelse(tri %in% c('Preconception', '1st trimester'), 'early', 'late')]
data_vls$time <- factor(data_vls$time, levels = c('early', 'late'))

form_binom_vls <- cbind(n_vls,(n.e-n_vls)) ~  class + time + class * time + (1 | studlab) + (1 | id)
model_time_vls <- glmmTMB(form_binom_vls, data = data_vls, family = binomial())
summary(model_time_vls)
posterior_vls <- get_posterior(model_time_vls)
saveRDS(posterior_vls, './model_output/posterior_draws_vls.RDS')
