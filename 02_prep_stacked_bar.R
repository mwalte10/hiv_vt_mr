# install.packages('meta')
#install.packages('ggforestplot')
rm(list = ls())
library(data.table)
library(ggplot2)
library(stringr)
library(tidyverse)
library(tidybayes)
library(ggpubr)

library(data.table)
library(ggplot2)
library(tidyverse)
library(grid)
devtools::load_all('C:/Users/mwalters/frogger/') ## Must be on mkw_spec_2024
source('C:/Users/mwalters/frogger//scripts/spectrum_inputs_paeds.R')
source('C:/Users/mwalters/frogger//scripts/read_spectrum.R')
source('C:/Users/mwalters/frogger/tests/testthat/helper-child-model.R')
source('./public_results/functions.R')

spec_file_dir <- './public_results/data/spectrum_files/'
spec_files <- list.files(paste0(spec_file_dir, '/pjnz/'), full.names = T)

###############################################################################
##Prepare Spectrum files
###############################################################################
##these take a while to run, so if they already exist don't rerun
if(!file.exists('./public_results/data/spectrum_files/proj/Rwanda.RDS')){
  lapply(spec_files, prep_files)
}

###############################################################################
##Prep Spectrum parameters
###############################################################################
prep_spectrum_format()

###############################################################################
##Country level stacked bars for appendix
###############################################################################
lapply(spec_files, get_stacked_bar_outputs)
