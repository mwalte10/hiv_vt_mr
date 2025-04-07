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

##frogger should be loaded from https://doi.org/10.5281/zenodo.15166687
library(frogger)
source('./functions.R')

spec_file_dir <- './data/spectrum_files/'
spec_files <- list.files(paste0(spec_file_dir, '/pjnz/'), full.names = T)

###############################################################################
##Prepare Spectrum files
###############################################################################
##these take a while to run, so if they already exist don't rerun
if(!file.exists('./data/spectrum_files/proj/Rwanda.RDS')){
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
