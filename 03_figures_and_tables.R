# install.packages('meta')
#install.packages('ggforestplot')
rm(list = ls())
library(data.table)
library(ggplot2)
library(stringr)
library(tidyverse)
library(tidybayes)
library(glmmTMB)
library(patchwork)
library(ggpubr)
##all files paths should be directed towards the public_results directory
input_data_dir <- './data/'
output_data_dir <- './data/model_data/'
figure_dir <- './figures_tables/appendix/forest_plots/'
source('./functions.R')
set.seed(925)

###############################################################################
##Tables 1 and 2 
###############################################################################
{
  spec_estimates <- readRDS('./model_output/spectrum_estimates.RDS')
  spec_estimates <- melt(spec_estimates, id.vars = c('type', 'model', 'variable'))
  spec_estimates[,value := plogis(value)]
  spec_estimates[type == 'peri' & variable != 40, value_round := sprintf("%.1f", value*100)]
  spec_estimates[type == 'peri' & variable == 40, value_round := sprintf("%.2f", value*100)]
  spec_estimates[type == 'bf', value_round := sprintf("%.2f", value*100)]
  spec_estimates[,value := NULL]
  spec_estimates <- dcast(spec_estimates, type + model + variable ~ variable.1, value.var = 'value_round')
  spec_estimates[,lab := paste0(median, ' (', lower, ', ', upper, ')')]  
  spec_estimates <- spec_estimates[,.(type, model, variable, lab)]
  write.csv(spec_estimates, './figures_tables/tab_1_2.csv', row.names = F)
}

###############################################################################
##Figure 2 (data and model results)
###############################################################################
{
  draws <- readRDS('./model_output/estimate_draws.RDS')
  data <- fread(paste0(input_data_dir,'/public_data.csv'))
  size_limits = c(0,5000)
  range_limits = c(0,5)
  text_size = 9
  pt_alpha_level = 0.5
  wrap_length_strip = 2
  
  data <- data[,.(model, cd4_mid, time_art = inferred_ART_weeks, n.e, prop_infected, NEW_2024, model2_pvt_sc, type = tt)]
  data[,NEW_2024 := ifelse(NEW_2024 == 'Y', 'Yes', 'No')]
  draws[,value := plogis(value)]
  draws <- draws[,.(value, 
                    median = median(value),
                    lower = quantile(value, 0.025),
                    upper = quantile(value, 0.975)), by = c('variable', 'type', 'model')]
  
  model_map <- data.table(model = c('model1', 'model1', 'model2', 'model2', 'model3', 'model4'), 
                          type = rep(c('peri', 'bf'), 3),
                          model_name = c('Model 1: No PVT', 'Model 1: No PVT, per month',
                                         'Model 2: Short-course PVT', 'Model 2: Short-course PVT',
                                         'Model 3: Perinatal, lifelong ART', 'Model 4: BF, lifelong ART, per month'))
  draws <- merge(draws, model_map, by = c('model', 'type'))
  data <- merge(data, model_map, by = c('model', 'type'))
  data[model == 'model4',model2_pvt_sc := ifelse(time_art == 40, 'onart', 'startart')]
  draws[variable == 'mat_sero', model_name := 'Model 2: Maternal seroconversion']
  data[model2_pvt_sc == 'mat_sero', model_name := 'Model 2: Maternal seroconversion']
  
  pvt_map <- data.table(model2_pvt_sc = c('mat_sero', 'sdnvp','sdnvp_lte350','sdnvp_gte350',  'dual_arv', 'opt_a', 'opt_b', 'startart', 'onart'), 
                        variable = c('mat_sero', 'sdnvp','sdnvp_lte350','sdnvp_gte350',  'dual_arv', 'opt_a', 'opt_b', 'startart', 'onart'), 
                        xlab = c('Maternal\nseroconversion', 'sdNVP', 'sdNVP, <350', 'sdNVP, >350', 'Dual ARV', 'Option A', 'Option B',
                                 'Initiated ART during pregnancy',
                                 'Initiated ART preconception'))
  draws <- merge(draws, pvt_map, all.x = T, by = c('variable'))
  data <- merge(data, pvt_map, all.x = T, by = c('model2_pvt_sc'))
  
  draws$xlab <- factor(draws$xlab, levels = c('Maternal\nseroconversion', 'sdNVP', 'sdNVP, <350', 'sdNVP, >350', 'Dual ARV', 'Option A', 'Option B',
                                              'Initiated ART preconception',
                                              'Initiated ART during pregnancy'))
  data$xlab <- factor(data$xlab, levels = c('Maternal\nseroconversion', 'sdNVP', 'sdNVP, <350', 'sdNVP, >350', 'Dual ARV', 'Option A', 'Option B',
                                            'Initiated ART preconception',
                                            'Initiated ART during pregnancy'
  ))
  
  
  
  mod1_gg_peri <- ggplot() +
    geom_point(data = data[model == 'model1' & type == 'peri'], aes((cd4_mid), prop_infected*100, size = n.e, col = NEW_2024), alpha = pt_alpha_level) +
    geom_line(data = unique(draws[model == 'model1' & type == 'peri',.(variable, median)]), aes(as.integer(variable), median*100), col = 'darkcyan') + 
    geom_ribbon(data = unique(draws[model == 'model1' & type == 'peri',.(variable, lower, upper)]), aes(as.integer(variable), ymin = lower*100, ymax = upper*100), alpha = 0.2, fill = 'darkcyan', col = NA) + 
    theme_bw(base_size = text_size) + 
    labs(x = 'CD4 midpoint (per cubic millimetre)', y = 'Perinatal transmission probability (%)', color = 'Study identified in 2024 review', size = 'Study size') + 
    ylim(0, 0.5*100) + 
    facet_wrap(~model_name, labeller = labeller(groupwrap = label_wrap_gen(wrap_length_strip))) + scale_size_continuous(limits = size_limits, range = range_limits) + 
    scale_color_manual(values = c('black', 'red'))  + 
    theme(legend.position = 'bottom', axis.title.x = element_text(margin = margin(t = 10)))
  mod1_gg_peri
  
  mod1_gg_bf <- ggplot() +
    geom_point(data = data[model == 'model1' & type == 'bf'], aes((cd4_mid), prop_infected*100, size = n.e, col = NEW_2024), alpha = pt_alpha_level) +
    geom_line(data = unique(draws[model == 'model1' & type == 'bf',.(variable, median)]), aes(as.integer(variable), median*100), col = 'darkcyan') + 
    geom_ribbon(data = unique(draws[model == 'model1' & type == 'bf',.(variable, lower, upper)]), aes(as.integer(variable), ymin = lower*100, ymax = upper*100), alpha = 0.2, fill = 'darkcyan', col = NA) + 
    theme_bw(base_size = text_size) + 
    labs(x = 'CD4 midpoint (per cubic millimetre)', y = 'Breastfeeding transmission probability (%)') + 
    ylim(0, 0.06*100) + 
    facet_wrap(~model_name, labeller = labeller(groupwrap = label_wrap_gen(wrap_length_strip))) + scale_size_continuous(limits = size_limits, range = range_limits) + 
    scale_color_manual(values = c('black', 'red')) +
    theme(    axis.title.x = element_text(margin = margin(t = 10)))  # Adds 10 points of margin above x-axis label
  mod1_gg_bf
  
  mod2_gg_peri <- ggplot() + 
    geom_violin(data = draws[xlab != 'Maternal\nseroconversion' & type == 'peri'], aes(xlab, value*100),  col = 'darkcyan', fill = 'darkcyan', alpha = 0.2) +
    geom_point(data = data[xlab != 'Maternal\nseroconversion' & type == 'peri'], aes(xlab, (prop_infected)*100, size = n.e, col = NEW_2024), alpha = pt_alpha_level) + 
    facet_wrap(~model_name, labeller = labeller(groupwrap = label_wrap_gen(wrap_length_strip))) + scale_size_continuous(limits = size_limits, range = range_limits) + 
    ylim(0, 0.2*100) +
    theme_bw(base_size = text_size) + 
    labs(x = 'Short-course PVT', y = NULL) + 
    scale_color_manual(values = c('black', 'red')) + 
    scale_x_discrete(labels = function(x) str_wrap(x, width = 5))
  mod2_gg_peri
  
  mod2_gg_peri_inf <- ggplot() + 
    geom_violin(data = draws[xlab == 'Maternal\nseroconversion' & type == 'peri'], aes(xlab, value*100),  col = 'darkcyan', fill = 'darkcyan', alpha = 0.2) +
    geom_point(data = data[xlab == 'Maternal\nseroconversion' & type == 'peri'], aes(xlab, (prop_infected)*100, size = n.e, col = NEW_2024), alpha = pt_alpha_level) + 
    facet_wrap(~model_name, labeller = label_wrap_gen()) + scale_size_continuous(limits = size_limits, range = range_limits) + 
    theme_bw(base_size = text_size)  + 
    ylim(0, 0.5*100) + 
    labs(x = 'Acute infection', y = NULL) +   scale_color_manual(values = c('black', 'red')) + 
    scale_x_discrete(labels = function(x) str_wrap(x, width = 5))
  
  mod2_gg_peri_inf
  
  mod2_gg_bf <- ggplot() + 
    geom_violin(data = draws[xlab != 'Maternal\nseroconversion' & type != 'peri' & model == 'model2'], aes(xlab, value*100),  col = 'darkcyan', fill = 'darkcyan', alpha = 0.2) +
    geom_point(data = data[xlab != 'Maternal\nseroconversion' & type != 'peri'& model == 'model2'], aes(xlab, (prop_infected)*100, size = n.e, col = NEW_2024), alpha = pt_alpha_level) + 
    facet_wrap(~model_name) + scale_size_continuous(limits = size_limits, range = range_limits) + 
    ylim(0, 0.03*100) +
    theme_bw(base_size = text_size) + 
    labs(x = 'Short-course PVT', y = NULL) + 
    scale_color_manual(values = c('black', 'red'))  + 
    scale_x_discrete(labels = function(x) str_wrap(x, width = 5))
  mod2_gg_bf
  
  mod2_gg_bf_inf <- ggplot() + 
    geom_violin(data = draws[xlab == 'Maternal\nseroconversion' & type != 'peri'], aes(xlab, value*100),  col = 'darkcyan', fill = 'darkcyan', alpha = 0.2) +
    geom_point(data = data[xlab == 'Maternal\nseroconversion' & type != 'peri'], aes(xlab, (prop_infected)*100, size = n.e, col = NEW_2024), alpha = pt_alpha_level) + 
    facet_wrap(~model_name, labeller = label_wrap_gen()) + scale_size_continuous(limits = size_limits, range = range_limits) + 
    theme_bw(base_size = text_size)  + 
    ylim(0, 0.5*100) + 
    labs(x = 'Acute infection', y = NULL) +
    scale_color_manual(values = c('black', 'red'))  + 
    scale_x_discrete(labels = function(x) str_wrap(x, width = 5)) 
  mod2_gg_bf_inf
  
  mod3_gg <- ggplot() + 
    geom_line(data = unique(draws[model == 'model3' & type == 'peri',.(variable, median)]), aes(-as.numeric(variable), median*100), col = 'darkcyan') + 
    geom_ribbon(data = unique(draws[model == 'model3' & type == 'peri',.(variable, lower, upper)]), aes(-as.numeric(variable), ymin = lower*100, ymax = upper*100), fill = 'darkcyan', col = NA, alpha = 0.2) +
    geom_point(data = data[model == 'model3' & type == 'peri'], aes(-time_art, (prop_infected)*100, size = n.e, col = NEW_2024), alpha = pt_alpha_level) + 
    facet_wrap(~model_name, labeller = labeller(groupwrap = label_wrap_gen(wrap_length_strip))) + scale_size_continuous(limits = size_limits, range = range_limits) + 
    theme_bw(base_size = text_size) + 
    ylim(0,0.2*100) + 
    scale_x_continuous(breaks = -seq(0,40, by= 10), labels = ( c('Delivery', paste0(seq(10,30, by = 10), '\nweeks'), 'Pre-\nconception')) ) +
    labs(x = 'Timing of ART intiation',y = NULL) + 
    scale_color_manual(values = c('black', 'red')) +
    theme(plot.margin = margin(5.5,13,5.5,5.5))
  mod3_gg
  
  mod4_gg <- ggplot() + 
    geom_violin(data = draws[model == 'model4'], aes(xlab, value*100), col = 'darkcyan', fill = 'darkcyan', alpha = 0.2) + 
    geom_point(data = data[model == 'model4'], aes(xlab, (prop_infected)*100, size = n.e, col = NEW_2024), alpha = pt_alpha_level) + 
    facet_wrap(~model_name, labeller = labeller(groupwrap = label_wrap_gen(wrap_length_strip))) + scale_size_continuous(limits = size_limits, range = range_limits) + 
    theme_bw(base_size = text_size) + 
    labs(x = 'Timing of ART intiation', y = NULL) + 
    scale_color_manual(values = c('black', 'red'))  + 
    scale_x_discrete(labels = function(x) str_wrap(x, width = 20))  +
    theme(plot.margin = margin(5.5,13,5.5,5.5))
  
  mod4_gg
  
  legend <- ggpubr::get_legend(mod1_gg_peri)
  mod1_gg_peri <- mod1_gg_peri + theme(legend.position = 'none')
  mod2_gg_peri_inf <- mod2_gg_peri_inf + theme(legend.position = 'none')
  mod2_gg_peri <- mod2_gg_peri + theme(legend.position = 'none')
  mod3_gg <- mod3_gg + theme(legend.position = 'none')
  
  peri <- grid.arrange(mod1_gg_peri, mod2_gg_peri_inf, mod2_gg_peri, mod3_gg, nrow = 1,
                       widths = c(1,.7,1,1))
  
  peri <- grid.arrange(mod1_gg_peri, mod2_gg_peri_inf, mod2_gg_peri, mod3_gg, nrow = 1,
                       widths = c(1,.7,1,1))
  
  mod1_gg_bf <- mod1_gg_bf + theme(legend.position = 'none')
  mod2_gg_bf_inf <- mod2_gg_bf_inf + theme(legend.position = 'none')
  mod2_gg_bf <- mod2_gg_bf + theme(legend.position = 'none')
  mod4_gg <- mod4_gg + theme(legend.position = 'none')
  
  bf <- grid.arrange(mod1_gg_bf, mod2_gg_bf_inf, mod2_gg_bf, mod4_gg, nrow = 1,
                     widths = c(1,.7,1,1))
  
  grid.arrange(peri, bf, legend, nrow = 3, heights = c(1,1, 0.2))
  w = 8.24
  h = 6
  
  png(paste0('./figures_tables/figure_2.png'), width = w, height =h ,res = 900, units = 'in')
  grid.arrange(peri, bf, legend, nrow = 3, heights = c(1,1, 0.2))
  dev.off()
  
  pdf(paste0('./figures_tables/figure_2.pdf'), width = w, height =h )
  grid.arrange(peri, bf, legend, nrow = 3, heights = c(1,1, 0.2))
  dev.off()
}

###############################################################################
##Figure 3 (stacked bar)
###############################################################################
{
  mwi <- fread(paste0('./data/spectrum_files/spectrum_format/inf_diff/Malawi.csv'))
  drc <- fread(paste0('./data/spectrum_files/spectrum_format/inf_diff/Democratic Republic of the Congo.csv'))
  bfa <- fread(paste0('./data/spectrum_files/spectrum_format/inf_diff/Burkina Faso.csv'))
  rwa <- fread(paste0('./data/spectrum_files/spectrum_format/inf_diff/Rwanda.csv'))
  
  dt <- rbind(mwi[,loc := 'Malawi'],
              drc[,loc := 'Democratic Republic of the Congo'],
              bfa[,loc := 'Burkina Faso'],
              rwa[,loc := 'Rwanda'])
  
  save <- copy(dt)
  dt <- dt[year == 2023]
  
  pct_diff <- unique(dt[,.(tt, value, run, loc)])
  pct_diff <- pct_diff[,.(value = sum(value)), by = c('tt', 'run', 'loc')]
  pct_diff_total <- unique(pct_diff[,.(tt = 'total', value = sum(value)), by = c('run', 'loc')])
  pct_diff <- rbind(pct_diff, pct_diff_total, fill = T)
  pct_diff <- dcast(pct_diff, tt + loc ~ run)
  pct_diff[,pct_diff := round((`Updated VT` - `Former VT`) / (`Former VT`) * 100)]
  pct_diff[,color := ifelse(pct_diff <0 , 'Red', 'Green')]
  pct_diff[,pct_diff := ifelse(pct_diff < 0, paste0(pct_diff, '%'), paste0('+', pct_diff, '%'))]
  pct_diff[,max := ifelse(`Updated VT` > `Former VT`, `Updated VT`, `Former VT`)]
  pct_diff[,y := ifelse(max < 1000, max*1.08 ,max *1.05)]
  pct_diff[loc == 'Burkina Faso' & tt == 'perinatal',y := 545 + 15]
  pct_diff[loc == 'Democratic Republic of the Congo' & tt == 'perinatal', y:= 3875+ 40]
  pct_diff[loc == 'Malawi' & tt == 'perinatal', y:= 950+ 20]
  pct_diff[loc == 'Rwanda' & tt == 'perinatal', y:= 180+ 5]
  
  pct_diff[loc == 'Burkina Faso' & tt != 'perinatal',y := 440 +37]
  pct_diff[loc == 'Democratic Republic of the Congo' & tt != 'perinatal', y:= 3100 + 250]
  pct_diff[loc == 'Malawi' & tt != 'perinatal', y:= 1370 + 68]
  pct_diff[loc == 'Rwanda' & tt != 'perinatal', y:= 282 + 8]
  pct_diff[,x := 'Updated VT']
  pct_diff_total <- pct_diff[tt == 'total',]
  pct_diff <- pct_diff[tt != 'total',]
  
  
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
  dt$broad_pmtct <- factor(dt$broad_pmtct, levels= c('Maternal seroconversion', 'No PVT',
                                                     'Started ART, dropped out',
                                                     'On ART, dropped out',
                                                     'PVT, non-ART',
                                                     'Started ART during pregnancy', 'On ART pre-conception'))
  
  color_map <- color_map[order(broad_pmtct),]
  dt$loc <- factor(dt$loc, levels = c('Rwanda', 'Burkina Faso', 'Malawi', 'Democratic Republic of the Congo'))
  
  
  dt[,max := sum(value), by = c('tt', 'run', 'loc')]
  
  get_max <- unique(dt[,.(tt, loc, max)])
  get_max[,max_val := max(max), by = c('loc')]
  get_max <- get_max[max == max_val,.(loc, max)]
  pct_diff[,max := NULL]
  pct_diff <- merge(pct_diff, get_max, by = c('loc'))
  pct_diff[,new_y := (`Updated VT` / max * 1.1) * max * 0.95]
  pct_diff[,y := new_y]
  
  
  perinatal <- ggplot()  +
    geom_bar(data = dt[tt == 'perinatal'], aes(run, value, fill = broad_pmtct),position="stack", stat="identity") +
    scale_fill_manual(values =color_map$color,
                      breaks = color_map$broad_pmtct,
                      labels = function(x) str_wrap(x, width = 12)) +
    geom_point(data = pct_diff, aes(x, new_y), shape = NA) + 
    theme_bw(base_size = 8)+
    theme(
      strip.text = element_text(size = rel(1)),
      axis.text = element_text(size = rel(1)),
      axis.title = element_text(size = rel(1.3))) +
    labs(fill = 'Transmission category') +
    guides(fill = guide_legend(title.position = 'top')) +
    labs(x = NULL, y = 'Perinatal infections\n(number)') +
    facet_wrap(~loc, nrow = 1, scales = 'free_y', 
               labeller = label_wrap_gen()) + 
    geom_text(data = pct_diff[tt == 'perinatal'], aes(x = x, y = y, label = pct_diff, col = color), show.legend = F,
              size = 4) +
    scale_color_manual(values = c('darkgreen', 'red'), breaks = c('Green', 'Red')) + 
    scale_y_continuous(expand = expansion(mult = c(0.05,0.05)))
  
  perinatal
  
  bf <- ggplot()  +
    geom_bar(data = dt[tt != 'perinatal'], aes(run, value, fill = broad_pmtct),position="stack", stat="identity") +
    scale_fill_manual(values =color_map$color,
                      breaks = color_map$broad_pmtct) +
    geom_point(data = pct_diff, aes(x, new_y), shape = NA) + 
    theme_bw(base_size = 8)+
    theme(strip.text = element_text(size = rel(1)),
          axis.text = element_text(size = rel(1)),
          axis.title = element_text(size = rel(1.3))) +
    labs(fill = 'Transmission category') +
    guides(fill = guide_legend(title.position = 'top')) +
    labs(y = 'Infections during breastfeeding\n(number)', x = NULL) +
    facet_wrap(~loc, nrow = 1, scales = 'free_y', 
               labeller = label_wrap_gen()) +
    geom_text(data = pct_diff[tt != 'perinatal'], aes(x = x, y = y, label = pct_diff, col = color), show.legend = F,
              size = 4) +
    scale_color_manual(values = c('darkgreen', 'red'), breaks = c('Green', 'Red'))
  bf
  
  
  xaxis <- ggplot() +
    labs(x = 'Vertical transmission probability') +
    theme_void() + # Remove all default plot elements
    theme(
      axis.title.x = element_text(size = 8, margin = margin(b = 0)),
      axis.line.x = element_blank(), # Optional: Add a visible x-axis line
      axis.ticks.x = element_blank(), # Remove ticks
      axis.text.x = element_blank()   # Remove tick labels
    ) +
    coord_cartesian(clip = "off") # Ensures space for the label if needed
  
  
  
  library(ggplot2)
  library(gridExtra)
  library(cowplot)
  
  
  legend <- ggpubr::get_legend(perinatal)
  perinatal <- perinatal + theme(legend.position = 'none')
  bf <- bf + theme(legend.position = 'none')
  
  pct_diff_total$loc <- factor(pct_diff_total$loc, levels = c('Rwanda', 'Burkina Faso', 'Malawi', 'Democratic Republic of the Congo'))
  pct_diff_total[,run := 'Updated VT']
  copy <- copy(pct_diff_total)
  pct_diff_total <- rbind(pct_diff_total, copy[,run := 'Former VT'])
  pct_diff_total[run == 'Former VT', pct_diff := 'Percent difference in\ntotal number of\n infections:']
  pct_diff_total$run <- factor(pct_diff_total$run, levels = c('Former VT', 'Updated VT'))

  total_pct_diff <- ggplot() +
    labs(y = ' \n ') +
    geom_point(data = pct_diff_total, aes(x = run, y = max), col = NA) + 
    geom_text(data = pct_diff_total[run == 'Former VT' & loc == 'Rwanda',], aes(x = run, y = max , label = pct_diff), col = 'black', size = 2.1) +
    geom_text(data = pct_diff_total[run != 'Former VT',], aes(x = run, y = max , label = pct_diff), col = 'darkgreen', size = 4, fontface = 'bold') +
    facet_wrap(~loc, nrow = 1, scales = 'free') + 
    theme_bw()  + theme(panel.grid.major = element_blank(), 
                        panel.grid.minor = element_blank(),
                        axis.line=element_blank(),
                        axis.text.x=element_blank(),
                        axis.text.y=element_blank(),axis.ticks=element_blank(),
                        axis.title.x=element_blank(),
                        #axis.title.y=element_blank(),
                        legend.position="none",
                        panel.border = element_blank(),
                        strip.text = element_blank() )

  total_pct_diff
  
  plots <- plot_grid(perinatal, bf, total_pct_diff, nrow = 3, rel_heights = c(1,1,0.3))
  plots
  
  final_plot <- plot_grid(
    plots, legend, ncol = 2, rel_widths = c(1,0.2)  # Adjust height of legend
  )
  h = 5.9
  w = 8.27
  dev.new(noRStudioGD = T, width = w, height = h)
  png(filename =paste0('./figures_tables/figure_3.png'), height = h, width = w, units = 'in', res = 900)
  final_plot
  dev.off()
  
  pdf(paste0('./figures_tables/figure_3.pdf'), height = h, width = w)
  final_plot
  dev.off()
}

###############################################################################
##Appendix Figures 3.2.1-3.2.15 (forest plots) 
###############################################################################
{
  prep_forest_plot_data()
  get_diamond()
  pdt <- fread(paste0('./input_for_forest_plots_public.csv'))
  diamond_in <- readRDS(paste0('./figures_tables/appendix/diamonds.RDS'))
  fp_parms <- fread(paste0('./figures_tables/forest_plot_space_parms_public.csv'))
  plot_fp(pdt = pdt, level_0_in = 'Model 1', level_1_in = 'peri')
  plot_fp(pdt = pdt, level_0_in = 'Model 1', level_1_in = 'bf')
  
  ##this needs to be ran interactively
  plot_fp(pdt = pdt, level_0_in = 'Model 2', level_1_in = 'peri')
  plot_fp(pdt = pdt, level_0_in = 'Model 2', level_1_in = 'bf')
  
  plot_fp(pdt = pdt, level_0_in = 'Model 3', level_1_in = 'peri')
  plot_fp(pdt = pdt, level_0_in = 'Model 4', level_1_in = 'bf')
}

###############################################################################
##Appendix table 4.1.1
###############################################################################
{
  posterior_estimates <- readRDS('./model_output/posterior_estimates.RDS')
  posterior_estimates <- posterior_estimates[,.(median = quantile(value, 0.5), 
                                                lower = quantile(value, 0.025),
                                                upper = quantile(value, 0.975)), by = c('variable', 'model')]
  posterior_estimates <- melt(posterior_estimates, id.vars = c('variable', 'model'))
  posterior_estimates[,or := exp(value)]
  ##model 2-4 odds ratio is three decimals
  posterior_estimates[, or_round := ifelse(model != 'model1',sprintf("%.3f", or), sprintf("%.2f", or))]
  posterior_estimates[, value_round := sprintf("%.2f", value)]
  posterior_estimates <- posterior_estimates[,.(variable, model, variable.1, or_round, value_round)]
  posterior_estimates <- melt(posterior_estimates, id.vars = c('variable', 'model', 'variable.1'))
  posterior_estimates <- dcast(posterior_estimates, variable + model + variable.2 ~ variable.1, value.var = 'value')
  posterior_estimates[, value := paste0(median, ' (', lower, ', ', upper, ')')]
  posterior_estimates <- dcast(posterior_estimates[,.(variable, model, variable.2, value)], variable + model ~ variable.2, value.var = 'value')
  posterior_estimates <- posterior_estimates[order(model),.(model, variable, value_round, or_round)]
  write.csv(posterior_estimates, './figures_tables/appendix/tab_4_1_1.csv', row.names = F) 
}

###############################################################################
##Appendix table 4.1.2
###############################################################################
{
  posterior_mod3_art <- readRDS('./model_output/posterior_draws_mod3_art.RDS')
  posterior_mod3_art <- melt(posterior_mod3_art, id.vars = 'draw')
  mod3_art_posterior_draws <- posterior_mod3_art[,.(median = quantile(value, 0.5), 
                                                    lower = quantile(value, 0.025),
                                                    upper = quantile(value, 0.975)), by = c('variable')]
  mod3_art_posterior_draws <- melt(mod3_art_posterior_draws, id.vars = 'variable')
  mod3_art_posterior_draws[,odds := round(exp(value),2)]
  mod3_art_posterior_draws[,value := round(value,2)]
  mod3_art_posterior_draws <- melt(mod3_art_posterior_draws, id.vars = c('variable', 'variable.1'))
  mod3_art_posterior_draws <- dcast(mod3_art_posterior_draws, variable + variable.2 ~ variable.1, value.var = 'value')
  mod3_art_posterior_draws[,value := paste0(median, ' (', lower, ', ', upper, ')')]
  mod3_art_posterior_draws <- mod3_art_posterior_draws[,.(variable, variable.2, value)]
  mod3_art_posterior_draws <- dcast(mod3_art_posterior_draws, variable ~ variable.2, value.var = 'value')
  write.csv(mod3_art_posterior_draws, './figures_tables/appendix/tab_4_1_2.csv', row.names = F)
  
}

###############################################################################
##Appendix table 4.1.3
###############################################################################
{
  posterior_mod3_art_region <- readRDS('./model_output/posterior_draws_mod3_art_region.RDS')
  posterior_mod3_art_region <- melt(posterior_mod3_art_region, id.vars = 'draw')
  mod3_art_region_posterior_draws <- posterior_mod3_art_region[,.(median = quantile(value, 0.5), 
                                                    lower = quantile(value, 0.025),
                                                    upper = quantile(value, 0.975)), by = c('variable')]
  mod3_art_region_posterior_draws <- melt(mod3_art_region_posterior_draws, id.vars = 'variable')
  mod3_art_region_posterior_draws[,odds := round(exp(value),3)]
  mod3_art_region_posterior_draws[,value := round(value,2)]
  mod3_art_region_posterior_draws <- melt(mod3_art_region_posterior_draws, id.vars = c('variable', 'variable.1'))
  mod3_art_region_posterior_draws <- dcast(mod3_art_region_posterior_draws, variable + variable.2 ~ variable.1, value.var = 'value')
  mod3_art_region_posterior_draws[,value := paste0(median, ' (', lower, ', ', upper, ')')]
  mod3_art_region_posterior_draws <- mod3_art_region_posterior_draws[,.(variable, variable.2, value)]
  mod3_art_region_posterior_draws <- dcast(mod3_art_region_posterior_draws, variable ~ variable.2, value.var = 'value')
  write.csv(mod3_art_region_posterior_draws, './figures_tables/appendix/tab_4_1_3.csv', row.names = F)
}

###############################################################################
##Appendix table 4.1.4
###############################################################################
{
  posterior_mod_vls <- readRDS('./model_output/posterior_draws_vls.RDS')
  posterior_mod_vls <- melt(posterior_mod_vls, id.vars = 'draw')
  mod_vls_posterior_draws <- posterior_mod_vls[,.(median = quantile(value, 0.5), 
                                                                  lower = quantile(value, 0.025),
                                                                  upper = quantile(value, 0.975)), by = c('variable')]
  mod_vls_posterior_draws <- melt(mod_vls_posterior_draws, id.vars = 'variable')
  mod_vls_posterior_draws[,odds := round(exp(value),2)]
  mod_vls_posterior_draws[,value := round(value,2)]
  mod_vls_posterior_draws <- melt(mod_vls_posterior_draws, id.vars = c('variable', 'variable.1'))
  mod_vls_posterior_draws <- dcast(mod_vls_posterior_draws, variable + variable.2 ~ variable.1, value.var = 'value')
  mod_vls_posterior_draws[,value := paste0(median, ' (', lower, ', ', upper, ')')]
  mod_vls_posterior_draws <- mod_vls_posterior_draws[,.(variable, variable.2, value)]
  mod_vls_posterior_draws <- dcast(mod_vls_posterior_draws, variable ~ variable.2, value.var = 'value')
  write.csv(mod_vls_posterior_draws, './figures_tables/appendix/tab_4_1_4.csv', row.names = F)
}


###############################################################################
##Appendix figure 5.2.1.2
###############################################################################
{
  data <- fread(paste0(input_data_dir,'/public_data.csv'))
  data[,studlab := paste0(author, '_', study_year)]
  data[,id := 1:nrow(data)]
  
  mod3_data <- data[model == 'model3']
  mod3_data[,time_on_art_median := inferred_ART_weeks - 20]
  mod3_data[,late_start := ifelse(inferred_ART_weeks < 4, 1, 0)]
  form_binom_mod3 <- cbind(event.e,(n.e-event.e)) ~  time_on_art_median + late_start + (1 | studlab) + (1 | id)
  model3 <- glmmTMB(form_binom_mod3, data = mod3_data, family = binomial())
  model_range <- glmmTMB(form_binom_mod3, data = mod3_data[code %in% c(3,4)], family = binomial())
  model_median <- glmmTMB(form_binom_mod3, data = mod3_data[code %in% c(1,2,4) | late_start == 1], family = binomial())
  
  model3 <- get_posterior(model = model3)
  model_range_post <- get_posterior(model = model_range)
  model_median_post <- get_posterior(model = model_median)
 
  tab <- rbind(model3[,mod := 'All data'],
               model_range_post[,mod := 'Range'],
               model_median_post[,mod := 'Median']) 
  tab <- melt(tab, id.vars = c('mod', 'draw'))
  tab <- tab[,.(median = quantile(value, 0.5),
                lower = quantile(value, 0.025),
                upper = quantile(value, 0.975)), by = c('mod', 'variable')]
  
  tab <- melt(tab, id.vars = c('mod', 'variable'))
  tab[,value := round(value,2)]
  tab <- dcast(tab, mod + variable ~ variable.1, value.var = 'value')
  tab <- tab[,.(mod, variable, value = paste0(median, ' (', lower, ', ', upper, ')'))]
  tab <- dcast(tab, mod ~ variable, value.var = 'value')
 
}


###############################################################################
##Appendix figure 5.2.1
###############################################################################
{
  data <- fread(paste0(input_data_dir,'/public_data.csv'))
  data[,studlab := paste0(author, '_', study_year)]
  data[,id := 1:nrow(data)]
  
  mod3_data <- data[model == 'model3']
  mod3_data[,time_on_art_median := inferred_ART_weeks - 20]
  mod3_data[,late_start := ifelse(inferred_ART_weeks < 4, 1, 0)]
  form_binom_mod3 <- cbind(event.e,(n.e-event.e)) ~  time_on_art_median + late_start + (1 | studlab) + (1 | id)
  model_range <- glmmTMB(form_binom_mod3, data = mod3_data[code %in% c(3,4)], family = binomial())
  model_median <- glmmTMB(form_binom_mod3, data = mod3_data[code %in% c(1,2,4) | late_start == 1], family = binomial())
  
  model_range_post <- get_posterior(model = model_range)
  model_median_post <- get_posterior(model = model_median)
  
  dt <- rbind(mod3_data[code %in% c(1,2,4) | late_start == 1,.(studlab, prop_infected, n.e, inferred_ART_weeks, mod = 'Median weeks reported')],
              mod3_data[code %in% c(3,4),.(studlab, prop_infected, n.e, inferred_ART_weeks, mod = 'Range of weeks reported')])
  
  mod3_range_draws <- rbindlist(lapply(0:40, get_model3_estimates, posterior = model_range_post))
  mod3_median_draws <- rbindlist(lapply(0:40, get_model3_estimates, posterior = model_median_post))
  mod3_sens_draw <- rbind(mod3_range_draws[,mod := 'Range of weeks reported'],
                          mod3_median_draws[,mod := 'Median weeks reported'])
  mod3_sens_draw[,value := plogis(value)]
  
  
  mod3_result_summary <- mod3_sens_draw[,.(median = quantile(value, 0.5), 
                                           lower = quantile(value, 0.025),
                                           upper = quantile(value, 0.975)), by = c('week', 'mod')]
  
  gg <- ggplot() + geom_point(data = dt, aes(-inferred_ART_weeks, prop_infected, size = n.e)) +
    geom_line(data = mod3_result_summary, aes(-week, median))  +
    geom_ribbon(data = mod3_result_summary, aes(-week, ymin = lower, ymax = upper), alpha = 0.2) +
    facet_wrap(~mod) +
    scale_x_continuous(breaks = seq(-40, 0, by = 10), labels = c('Pre-\nconception', 30, 20, 10, 'Delivery')) +
    labs(x = 'Weeks on ART', size = 'Study size', y = 'Proportion infected') +
    theme_bw(base_size = 11) + theme(legend.position = 'bottom',
                                     legend.margin = margin(0, 0, 0, 0),        # Remove inner padding within the legend box
                                     legend.box.margin = margin(-10, 0, 0, 0),
                                     axis.text.x = element_text(size = rel(0.8)))
  
  
  h = 3.5
  w = 6.24
  
  png(paste0('./figures_tables/appendix/5_2_1.png'), height = h, width = w, res = 900, units = 'in')
  print(gg)
  dev.off()
}

###############################################################################
##Appendix figures 6.1-6.16
###############################################################################
spec_file_dir <- './spectrum_files/'
spec_files <- list.files(paste0(spec_file_dir, '/pjnz/'), full.names = T)
mapply(plot_country_results, rep(unlist(lapply(spec_files, eppasm::read_country)),4), year.x = rep(c(2000,2010,2015, 2023),each = 4))

###############################################################################
##Appendix figure 6.17
###############################################################################
{
  dt <- save
  dt <- dt[,.(inf = sum(value)), by = c('year', 'tt', 'run', 'loc')]
  dt <- dcast(dt, year + tt + loc ~ run, value.var = 'inf')
  dt[,pct_diff := 100 * (`Updated VT` - `Former VT`)  / (`Former VT`)]
  dt[tt == 'perinatal' & year %in% c(2000,2010,2015,2023),pct_diff] %>% abs() %>% max()
  dt[tt == 'perinatal' & year %in% c(2000,2010,2015,2023),pct_diff] %>% mean()
  dt[tt != 'perinatal' & year %in% c(2000,2010,2015,2023),pct_diff] %>% abs() %>% range()
  dt[tt != 'perinatal' & year %in% c(2000,2010,2015,2023),pct_diff] %>% mean()
  dt <- dt[year%in% c(2000,2010,2015,2023)]
  
  dt <- dt
  dt[,tt := ifelse(tt == 'bf', 'Breastfeeding', 'Perinatal')]
  
  dt[,pos := pct_diff * 1.15]
  dt[pct_diff < 18, pos := pct_diff * 1.2]
  dt[tt == 'Perinatal', pos := pct_diff - 3]
  dt[tt == 'Breastfeeding' & pct_diff < 6, pos := pct_diff + 3]
  dt[,lab := paste0(round(pct_diff, digits = 1), '%')]
  dt[pct_diff > 0, lab := paste0('+', lab)]
  dt$tt <- factor(dt$tt, levels = c('Perinatal', 'Breastfeeding'))
  
  tornado_plot <- ggplot(dt, aes(y = loc, x = pct_diff, fill = tt)) +
    geom_bar(stat = "identity", width = 0.7) +
    theme_minimal() +
    labs(
      x = "Percent difference from former VT probabilities",
      y = "",
      fill = 'Infection type'
    ) +
    scale_y_discrete(labels = scales::label_wrap(20)) + 
    scale_fill_manual(values  =c('Perinatal' = '#0f86b6', 'Breastfeeding' =  wesanderson::wes_palette("Zissou1")[4])) +
    theme_bw() + 
    theme(legend.position = 'bottom') + 
    facet_wrap(~year, ncol = 1) + 
    xlim(-15,30) +
    geom_vline(xintercept  = 0)
  
  h = 6
  w = 4.7
  png(filename = paste0('./figures_tables/appendix/6_17.png'), width = w, height = h, units = 'in', res = 900)
  tornado_plot  
  dev.off()
}

