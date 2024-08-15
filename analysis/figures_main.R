# filename:    figures_main.R
# created:     06 March 2023
# updated:     24 July 2024
# author:      S.C. McClelland
# description: This file creates main text figures for manuscript.
#-------------------------------------------------------------------------------
library(data.table)
library(forcats)
library(ggplot2)
library(ggpubr)
library(ggridges)
library(patchwork)
library(RColorBrewer)
library(rstudioapi)
library(sf)
library(terra)
#-------------------------------------------------------------------------------
source('results_functions.R')
source('plot_discrete_color_bar.R')
#-------------------------------------------------------------------------------
options(scipen = 999, digits = 4)
options(rgl.printRglwidget = TRUE)
base.path = dirname(getActiveDocumentContext()$path)
data.path = paste(base.path, 'manuscript-data', sep = '/')
fig.path  = paste(base.path, 'manuscript-figures', sep = '/')
lu.path   = paste(base.path, 'gis', sep = '/')
raster    = 'msw-cropland-rf-ir-area.tif'
raster_cc = 'msw-masked-cropland-rf-ir-area.tif'
shp       = 'WB_countries_Admin0_10m.shp'
#-------------------------------------------------------------------------------
# FIGURE 1 | FULL ADOPTION, YIELD CONSTRAINED, REVENUE POTENTIAL
#-------------------------------------------------------------------------------
# GLOBAL
unconstr_2050     = fread(paste(data.path, 'global-unconstrained-GHG-mitigation-potential-2050.csv', sep = '/'))
yield_cst_2050   = fread(paste(data.path, 'global-yield-constrained-GHG-mitigation-potential-2050.csv', sep = '/'))
yield_cst_2050_l = fread(paste(data.path, 'global-yield-constrained-low-price-GHG-mitigation-potential-2050.csv', sep = '/'))
yield_cst_2050_h = fread(paste(data.path, 'global-yield-constrained-high-price-GHG-mitigation-potential-2050.csv', sep = '/'))
# REGIONAL
ipcc_unconstr_2050     = fread(paste(data.path, 'regional-unconstrained-GHG-mitigation-potential-2050.csv', sep = '/'))
ipcc_yield_cst_2050   = fread(paste(data.path, 'regional-yield-constrained-GHG-mitigation-potential-2050.csv', sep = '/'))
ipcc_yield_cst_2050_l = fread(paste(data.path, 'regional-yield-constrained-low-price-GHG-mitigation-potential-2050.csv', sep = '/'))
ipcc_yield_cst_2050_h = fread(paste(data.path, 'regional-yield-constrained-high-price-GHG-mitigation-potential-2050.csv', sep = '/'))

# background map
map_background = IPCC_map(lu.path, shp, raster)
ggsave(paste0(fig.path,'/IPCC-map-background.pdf'), map_background$IPCC, units = 'mm', width = 180, height = 185, device='pdf', dpi=300)
# global
global_126_ptl = gl_bar_potential(unconstr_2050, yield_cst_2050, yield_cst_2050_l, yield_cst_2050_h,'ssp126')
ggsave(paste0(fig.path,'/ssp126-2050-global-ghg-mitigation-potential.pdf'), global_126_ptl$GHG, bg = 'transparent',units = 'mm', width = 58, height = 58, device='pdf', dpi=300)
# regional
ADP_126_ptl    = rg_bar_potential(ipcc_unconstr_2050, ipcc_yield_cst_2050, ipcc_yield_cst_2050_l, ipcc_yield_cst_2050_h, 'ssp126', 'ADP')
ggsave(paste0(fig.path,'/ssp126-2050-ADP-ghg-mitigation-potential.pdf'), ADP_126_ptl$GHG, bg = 'transparent',units = 'mm', width = 58, height = 58, device='pdf', dpi=300)
AME_126_ptl    = rg_bar_potential(ipcc_unconstr_2050, ipcc_yield_cst_2050, ipcc_yield_cst_2050_l, ipcc_yield_cst_2050_h, 'ssp126', 'AME')
ggsave(paste0(fig.path,'/ssp126-2050-AME-ghg-mitigation-potential.pdf'), AME_126_ptl$GHG, bg = 'transparent',units = 'mm', width = 58, height = 58, device='pdf', dpi=300)
DEV_126_ptl    = rg_bar_potential(ipcc_unconstr_2050, ipcc_yield_cst_2050, ipcc_yield_cst_2050_l, ipcc_yield_cst_2050_h, 'ssp126', 'DEV')
ggsave(paste0(fig.path,'/ssp126-2050-DEV-ghg-mitigation-potential.pdf'), DEV_126_ptl$GHG, bg = 'transparent',units = 'mm', width = 58, height = 58, device='pdf', dpi=300)
EEWCA_126_ptl  = rg_bar_potential(ipcc_unconstr_2050, ipcc_yield_cst_2050, ipcc_yield_cst_2050_l, ipcc_yield_cst_2050_h, 'ssp126', 'EEWCA')
ggsave(paste0(fig.path,'/ssp126-2050-EEWCA-ghg-mitigation-potential.pdf'), EEWCA_126_ptl$GHG, bg = 'transparent',units = 'mm', width = 58, height = 58, device='pdf', dpi=300)
LAC_126_ptl    = rg_bar_potential(ipcc_unconstr_2050, ipcc_yield_cst_2050, ipcc_yield_cst_2050_l, ipcc_yield_cst_2050_h, 'ssp126', 'LAC')
ggsave(paste0(fig.path,'/ssp126-2050-LAC-ghg-mitigation-potential.pdf'), LAC_126_ptl$GHG, bg = 'transparent',units = 'mm', width = 58, height = 58, device='pdf', dpi=300)
#-------------------------------------------------------------------------------
# FIGURE 2 | MAPS
#-------------------------------------------------------------------------------
Mg_ha = 100L
kg_ha = 10L
C_bio = 0.42 # Ma et al. 2018 | value for 'reproductive organs'
# LOAD DATA
# NTILL-RES
ntill_res_dt = readRDS(paste(data.path, 'imputed-ensemble-relative-responses-weighted-mean-ntill-res.rds', sep = '/'))
ntill_res_dt = dcast(ntill_res_dt,
                     cell + x + y + y_block + ssp + total_crop_area_ha ~ variable,
                     value.var = 'value')
ntill_res_dt[, m_cr_grain := (m_cr_grain*kg_ha)/C_bio]
ntill_res_dt[, s_cr_grain := (s_cr_grain/Mg_ha)/C_bio]

# CCG-NTILL-RES
ccg_ntill_dt = readRDS(paste(data.path, 'imputed-ensemble-relative-responses-weighted-mean-ccg-ntill.rds', sep = '/'))
ccg_ntill_dt = dcast(ccg_ntill_dt,
                     cell + x + y + y_block + ssp + total_crop_area_ha ~ variable,
                     value.var = 'value')
ccg_ntill_dt[, m_cr_grain := (m_cr_grain*kg_ha)/C_bio]
ccg_ntill_dt[, s_cr_grain := (s_cr_grain/Mg_ha)/C_bio]
# CCL-NTILL-RES
ccl_ntill_dt = readRDS(paste(data.path, 'imputed-ensemble-relative-responses-weighted-mean-ccl-ntill.rds', sep = '/'))
ccl_ntill_dt = dcast(ccl_ntill_dt,
                     cell + x + y + y_block + ssp + total_crop_area_ha ~ variable,
                     value.var = 'value')
ccl_ntill_dt[, m_cr_grain := (m_cr_grain*kg_ha)/C_bio]
ccl_ntill_dt[, s_cr_grain := (s_cr_grain/Mg_ha)/C_bio]

# MAPS | SSP1-2.6, 2050
ntill_res_ssp126_2050 = annual_map(ntill_res_dt, 2050, 'ssp126')
ccg_ntill_ssp126_2050 = annual_map(ccg_ntill_dt, 2050, 'ssp126')
ccl_ntill_ssp126_2050 = annual_map(ccl_ntill_dt, 2050, 'ssp126')

ssp126_2050_maps      = ntill_res_ssp126_2050$GHG + ntill_res_ssp126_2050$Yield +
  ccg_ntill_ssp126_2050$GHG + ccg_ntill_ssp126_2050$Yield +
  ccl_ntill_ssp126_2050$GHG + ccl_ntill_ssp126_2050$Yield +
  plot_layout(ncol = 2, nrow = 3, guides = 'collect') &
  theme(legend.position = 'none')
ssp126_2050_maps     = ssp126_2050_maps + plot_annotation(tag_levels = list(c('a', 'd', 'b', 'e', 'c', 'f'))) &
  theme(plot.tag = element_text(size = 7))

ggsave(paste0(fig.path,'/', 'ssp126-annual-2050-ghg-yield.pdf'),     ssp126_2050_maps,  units = 'mm', width = 180, height = 185, device='pdf', dpi=300)
ggsave(paste0(fig.path,'/', 'ssp126-annual-2050-ghg-legend.pdf'),    ntill_res_ssp126_2050$legend1,  units = 'mm', width = 180, height = 50, device='pdf', dpi=300)
ggsave(paste0(fig.path,'/', 'ssp126-annual-2050-yield-legend.pdf'),  ntill_res_ssp126_2050$legend2,  units = 'mm', width = 180, height = 50, device='pdf', dpi=300)
#-------------------------------------------------------------------------------
# FIGURE 3 | BEST MANAGEMENT PRACTICE
#-------------------------------------------------------------------------------
IPCC_dt = ipcc_name(lu.path, shp, raster)
ghg_bmp_dt = fread(paste(data.path, 'best-management-practice-GHG-mitigation-2050.csv', sep = '/'))
ghg_bmp_dt = ghg_bmp_dt[IPCC_dt[, .(cell, IPCC_NAME, WB_NAME)], 
                        on = .(cell = cell) ]
ghg_bmp_dt = ghg_bmp_dt[!is.na(ssp)]

ym_bmp_dt = fread(paste(data.path, 'best-management-practice-GHG-mitigation-yield-maximum-2050.csv', sep = '/'))
ym_bmp_dt = ym_bmp_dt[IPCC_dt[, .(cell, IPCC_NAME, WB_NAME)], 
                      on = .(cell = cell) ]
ym_bmp_dt = ym_bmp_dt[!is.na(ssp)]

yc_bmp_dt  = fread(paste(data.path, 'best-management-practice-GHG-mitigation-yield-constrained-2050.csv', sep = '/'))
yc_bmp_dt  = yc_bmp_dt[IPCC_dt[, .(cell, IPCC_NAME, WB_NAME)], 
                       on = .(cell = cell) ]
yc_bmp_dt  = yc_bmp_dt[!is.na(ssp)]

# GHG | SSP1-2.6
ghg_126_bmp = bmp_map(ghg_bmp_dt, 'ssp126')

# YM | SSP1-2.6
ym_126_bmp = bmp_map(ym_bmp_dt, 'ssp126')

# YC  | SSP1-2.6
yc_126_bmp = bmp_map(yc_bmp_dt, 'ssp126')

# Proportions by region, YC | SSP1-2.6
yc_126_bmp_bar = bmp_bplot(yc_bmp_dt, 'ssp126')

bmp_plots = ghg_126_bmp$bmp + yc_126_bmp_bar$freq + ym_126_bmp$bmp + yc_126_bmp_bar$ghg + yc_126_bmp$bmp +
  plot_layout(ncol = 2, nrow = 3, guides = 'collect') &
  theme(legend.position = 'none')

bmp_plots = bmp_plots + plot_annotation(tag_levels = list(c('a', 'd', 'b', 'e', 'c'))) &
  theme(plot.tag = element_text(size = 7))

# save 
ggsave(paste0(fig.path,'/', 'ssp126-bmp-plot.pdf'),bmp_plots,  units = 'mm', width = 180, height = 185, device='pdf', dpi=300)
#-------------------------------------------------------------------------------
# FIGURE 4 | GLOBAL HECTARE RESPONSE
#-------------------------------------------------------------------------------
gl_hectare_gcm  = fread(paste(data.path, 'gcm-decadal-global-crop-soil-responses.csv', sep = '/'))
rg_hectare_gcm  = fread(paste(data.path, 'gcm-decadal-IPCC-region-crop-soil-responses.csv', sep = '/'))

# SSP1-2.6 | 2050
# global
soc_gl_126_2050 = soc_gl_ha(gl_hectare_gcm, 'ssp126', 2050)
ggsave(soc_gl_126_2050, file = paste(fig.path, 'gcm-soc-seq-rate-global-ssp126-2050.pdf', sep = '/'), units = 'mm', width = 88, height = 88, device='pdf', dpi=300)
#-------------------------------------------------------------------------------
# Extended Data 1-4 | SUPERVISED CLUSTERING
#-------------------------------------------------------------------------------
# Fig. 1 SC Features, SOC and Yield, 1G-0L-0T-1R and 0G-1L-0T-1R
# SOC
load(paste(data.path, "ccl-syn-soc_SHAP.Rdata", sep = '/'))
load(paste(data.path, "ccg-syn-soc_SHAP.Rdata", sep = '/'))
# Yield
load(paste(data.path, "ccl-ntill-yield_SHAP.Rdata", sep = '/'))
load(paste(data.path, "ccg-ntill-yield_SHAP.Rdata", sep = '/'))

colors = c('Climate' = "#1B9E77", 'Site' = "#D95F02", 'Management' = "#7570B3", 'State' = "#E7298A")
# CCL-NTILL | SOC
ccl_f_names    = c('Soil Decomposition Multiplier', 'Surface Decomposition Multiplier','MAT',
                   'ET','Sand Fraction',parse(text = 'Delta~Grain~Yield'),
                   parse(text = 'Delta~N~Leach'), parse(text = 'Delta~Soil~Decomposition~Multiplier'), parse(text = 'Delta~N~Fixation'),
                   parse(text = 'Delta~ET'),'N Fixation','N Inputs',
                   parse(text = 'Delta~Residue~Biomass'), 'Initial Soil C Stock', 'Cover Crop Biomass') # reverse order
ccl_f_types    = c('State', 'Site', 'State', 'Management', 'State', 'State', 'State',
                   'State', 'State', 'State', 'Site', 'State', 'Climate', 'State', 'State') # actual order
ccl_feature    = feature_p(ccl_f_names, ccl_SOC$gl_gg$data$feature, 
                           ccl_SOC$gl_gg$data$value, colors, ccl_f_types)

# CCL-NTILL | Yield
ccl_ntill_y_f_names    = c('Nitrate Fertilizer Fraction', parse(text = 'Delta~Surface~Decomposition~Multiplier'),'Residue Return Fraction',
                           'MAT', 'Cover Crop Biomass', 'Clay Fraction',
                           'Soil Bulk Density', 'Sand Fraction', parse(text = 'Delta~Soil~Organic~Carbon'),
                           parse(text = 'Delta~N~Leached'), parse(text = 'Delta~Gross~Nitrification'), 'N Inputs',
                           parse(text='Delta~ET'), 'MAP', parse(text = 'Delta~N~Fixation')) # reverse order
ccl_ntill_y_f_types    = c('State', 'Climate', 'State', 'Management', 'State', 'State', 'State',
                           'Site', 'Site', 'Site', 'State', 'Climate', 'Management', 'State', 'Management') # actual order
ccl_ntill_y_feature    = feature_yield_p(ccl_ntill_y_f_names, ccl_ntill_YIELD$gl_gg$data$feature,
                                         ccl_ntill_YIELD$gl_gg$data$value, colors, ccl_ntill_y_f_types)
ccl_ntill_y_feature    = ccl_ntill_y_feature + coord_cartesian(xlim =c(0, 12.5))


# CCG-NTILL | SOC
ccg_f_names    = c(parse(text = 'Delta~Cover~Crop~Biomass'), parse(text = 'Delta~ET'),'Urea Fertilizer Fraction',
                   'MAT',parse(text = 'Delta~Grain~Yield'),parse(text = 'Delta~Soil~Decomposition~Multiplier'),
                   'Initial Residue Return Fraction', 'Sand Fraction', parse(text = 'Delta~Gross~Nitrification'),
                   parse(text = 'Delta~Residue~Biomass'),'Gross Nitrification','N Inputs',
                   parse(text = 'Delta~N~Leach'), 'Initial Soil C Stock', 'Cover Crop Biomass') # reverse order
ccg_f_types    = c('State', 'Site', 'State', 'Management', 'State', 'State', 'State',
                   'Site', 'Management', 'State', 'State', 'Climate', 'Management', 'State', 'State') # actual order
ccg_feature    = feature_p(ccg_f_names, ccg_SOC$gl_gg$data$feature, 
                           ccg_SOC$gl_gg$data$value, colors, ccg_f_types)

# CCG-NTILL | Yield
ccg_ntill_y_f_names    = c(parse(text = 'Delta~Surface~Decomposition~Multiplier'), 'Initial SOC stock', parse(text = 'Delta~Soil~Decomposition~Multiplier'),
                           'N Inputs', 'MAT', 'Clay Fraction',
                           'Soil Bulk Density', parse(text = 'Delta~ET'), 'Residue Return Fraction',
                           'Sand Fraction', 'MAP', parse(text = 'Delta~Gross~Nitrification'),
                           'Cover Crop Biomass',parse(text = 'Delta~N~Leached'), parse(text = 'Delta~Soil~Organic~Carbon')) # reverse order
ccg_ntill_y_f_types    = c('State', 'State', 'State', 'State', 'Climate', 'Site', 'Management',
                           'State', 'Site', 'Site', 'Climate', 'Management', 'State', 'Site', 'State') # actual order
ccg_ntill_y_feature    = feature_yield_p(ccg_ntill_y_f_names, ccg_ntill_YIELD$gl_gg$data$feature,
                                         ccg_ntill_YIELD$gl_gg$data$value, colors, ccg_ntill_y_f_types)
ccg_ntill_y_feature    = ccg_ntill_y_feature + coord_cartesian(xlim =c(0, 4))

# COMBINE
Feature_Rank   = ccl_feature + ccl_ntill_y_feature + ccg_feature + ccg_ntill_y_feature +
  plot_layout(ncol = 2, nrow = 2, guides = 'collect') &
  theme(legend.position = 'bottom')
Feature_Rank = Feature_Rank + plot_annotation(tag_levels = list(c('a', 'c', 'b', 'd'))) &
  theme(plot.tag = element_text(size = 7))
ggsave(paste0(fig.path,'/ext-data/SOC-Yield-SHAP-features-all.pdf'), Feature_Rank, units = 'mm', width = 180, height = 185, device='pdf', dpi=300)

# Fig. 2 SC Main Features, SOC, 1G-0L-0T-1R and 0G-1L-0T-1R
# CCL-NTILL
ccl_map     = cluster_map(ccl_k_SOC$kmeans_dt, 'ssp126')
ccl_density = SHAP_density(copy(ccl_k_SOC$kmeans_dt), 'ssp126', 
                           c('s_SOC', 's_cc_biomass', 'd_s_cr_residC', 'N.amt'),
                           c('s_SOC'         = "Soil Organic Carbon Difference\n(Mg CO\u2082-eq ha\u207b\u00b9 yr\u207b\u00b9)",
                             's_cc_biomass'  = "Cover Crop Biomass\n(Mg DM ha\u207b\u00b9 yr\u207b\u00b9)",
                             'd_s_cr_residC' = "Residue Biomass Difference\n(Mg DM ha\u207b\u00b9 yr\u207b\u00b9)",
                             'N.amt'         = "Nitrogen Inputs\n(g N m\u207b\u00b2 yr\u207b\u00b9)"))

# CCG-NTILL
ccg_map     = cluster_map(ccg_k_SOC$kmeans_dt, 'ssp126')
ccg_density = SHAP_density(copy(ccg_k_SOC$kmeans_dt), 'ssp126', 
                           c('s_SOC', 's_cc_biomass', 'd_s_N_leach', 'N.amt'),
                           c('s_SOC'         = "Soil Organic Carbon Difference\n(Mg CO\u2082-eq ha\u207b\u00b9 yr\u207b\u00b9)",
                             's_cc_biomass'  = "Cover Crop Biomass\n(Mg DM ha\u207b\u00b9 yr\u207b\u00b9)",
                             'd_s_N_leach'   = "Leached Nitrogen Difference\n(g N m\u207b\u00b2 yr\u207b\u00b9)",
                             'N.amt'         = "Nitrogen Inputs\n(g N m\u207b\u00b2 yr\u207b\u00b9)"))
SOC_main_features = ccl_density + ccg_density + ccl_map + ccg_map +
  plot_layout(ncol = 2, nrow = 2, guides = 'collect') &
  theme(legend.position = 'none')
SOC_main_features = SOC_main_features + plot_annotation(tag_levels = list(c('a', 'c', 'b', 'd'))) &
  theme(plot.tag = element_text(size = 7))
# add legend manually
ggsave(paste0(fig.path,'/ext-data/SOC-main-features-all.tiff'), SOC_main_features, units = 'mm', width = 180, height = 160, device='tiff', dpi=300)
# Fig. 3 SC Main Features, Yield, 1G-0L-0T-1R and 0G-1L-0T-1R
# CCL-NTILL
ccl_ntill_y_map     = cluster_map(ccl_ntill_k_YIELD$kmeans_dt, 'ssp126')
ccl_ntill_y_density = SHAP_y_density(copy(ccl_ntill_k_YIELD$kmeans_dt), 'ssp126', 
                                        c('s_cr_grain', 'm_nfix', 'hist_bio12', 's_annet'),
                                        c('s_cr_grain'    = "Grain Yield Difference\n(Mg ha\u207b\u00b9 yr\u207b\u00b9)",
                                          'm_nfix'        = 'N Fixation Difference\n(g N m\u207b\u00b2 yr\u207b\u00b9)',
                                          'hist_bio12'    = 'Mean Annual Precipitation\n(mm yr\u207b\u00b9)',
                                          's_annet'       = 'ET Difference\n(cm H\u2082O yr\u207b\u00b9)'))

# CCG-NTILL
ccg_ntill_y_map     = cluster_map(ccg_ntill_k_YIELD$kmeans_dt, 'ssp126')
ccg_ntill_y_density = SHAP_y_density(copy(ccg_ntill_k_YIELD$kmeans_dt), 'ssp126', 
                                     c('s_cr_grain', 's_SOC', 's_N_leach', 's_cc_biomass'),
                                     c('s_cr_grain'    = "Grain Yield Difference\n(Mg ha\u207b\u00b9 yr\u207b\u00b9)",
                                       's_SOC'         = "Soil Organic Carbon Difference\n(Mg CO\u2082-eq ha\u207b\u00b9 yr\u207b\u00b9)",
                                       's_N_leach'     = "Leached Nitrogen Difference\n(g N m\u207b\u00b2 yr\u207b\u00b9)",
                                       's_cc_biomass'  = "Cover Crop Biomass\n(Mg DM ha\u207b\u00b9 yr\u207b\u00b9)"))
Yield_main_features = ccl_ntill_y_density + ccg_ntill_y_density + ccl_ntill_y_map + ccg_ntill_y_map +
  plot_layout(ncol = 2, nrow = 2, guides = 'collect') &
  theme(legend.position = 'none')
Yield_main_features = Yield_main_features + plot_annotation(tag_levels = list(c('a', 'c', 'b', 'd'))) &
  theme(plot.tag = element_text(size = 7))
# add legend manually
ggsave(paste0(fig.path,'/ext-data/yield-main-features-all.tiff'), Yield_main_features, units = 'mm', width = 180, height = 160, device='tiff', dpi=300)
# Fig. 4 SC Feature Rank, other practices
load(paste(data.path, "res-yield_SHAP.Rdata", sep = '/'))
load(paste(data.path, "ntill-res-yield_SHAP.Rdata", sep = '/'))
load(paste(data.path, "ccl-res-yield_SHAP.Rdata", sep = '/'))
load(paste(data.path, "ccg-res-yield_SHAP.Rdata", sep = '/'))

# RES
res_y_f_names    = c(parse(text = 'Delta~Soil~Decomposition~Multiplier'), 'Ammonium Fertilizer Fraction',parse(text = 'Delta~N~Fixation'),
                     parse(text = 'Delta~Surface~Decomposition~Multiplier'),'MAP','MAT',
                     'Sand Fraction', 'Nitrate Fertilizer Fraction', 'Urea Fertilizer Fraction',
                     'N Inputs', 'Residue Return Fraction', parse(text = 'Delta~N~Leached'),
                     parse(text='Delta~Gross~Nitrification'),parse(text = 'Delta~ET'), parse(text = 'Delta~Soil~Organic~Carbon')) # reverse order
res_y_f_types    = c('State', 'State', 'State', 'State', 'Management', 'Management', 'Management',
                     'Management', 'Site', 'Climate', 'Climate', 'State', 'State', 'Management', 'State') # actual order
res_y_feature    = feature_yield_p(res_y_f_names, res_YIELD$gl_gg$data$feature,
                                   res_YIELD$gl_gg$data$value, colors, res_y_f_types)
# NTILL-RES
ntill_res_y_f_names    = c('Initial SOC Stock', parse(text = 'Delta~Surface~Decomposition~Multiplier'),'Ammonium Fertilizer Fraction',
                           parse(text = 'Delta~Soil~Decomposition~Multiplier'),'N Inputs','Urea Fertilizer Fraction',
                           'MAP', 'Residue Return Fraction', 'Nitrate Fertilizer Fraction',
                           parse(text = 'Delta~Soil~Organic~Carbon'), 'MAT', 'Sand Fraction',
                           parse(text='Delta~ET'),parse(text = 'Delta~Gross~Nitrification'), parse(text = 'Delta~N~Leached')) # reverse order
ntill_res_y_f_types    = c('State', 'State', 'State', 'Site', 'Climate', 'State', 'Management',
                           'Management', 'Climate', 'Management', 'Management', 'State', 'Management', 'State', 'Site') # actual order
ntill_res_y_feature    = feature_yield_p(ntill_res_y_f_names, ntill_res_YIELD$gl_gg$data$feature,
                                         ntill_res_YIELD$gl_gg$data$value, colors, ntill_res_y_f_types)
# CCL-RES
ccl_res_y_f_names    = c('Urea Fertilizer Fraction', 'MAT','Nitrate Fertilizer Fraction',
                         parse(text = 'Delta~Surface~Decomposition~Multiplier'),parse(text = 'Delta~Gross~Nitrification'), 'Cover Crop Biomass',
                         'Clay Fraction', 'Sand Fraction', 'Soil Bulk Density',
                         'N Inputs', 'MAP', parse(text = 'Delta~N~Leached'),
                         parse(text='Delta~ET'),parse(text = 'Delta~Soil~Organic~Carbon'), parse(text = 'Delta~N~Fixation')) # reverse order
ccl_res_y_f_types    = c('State', 'State', 'State', 'State', 'Climate', 'Management', 'Site',
                         'Site', 'Site', 'State', 'State', 'State', 'Management', 'Climate', 'Management') # actual order
ccl_res_y_feature    = feature_yield_p(ccl_res_y_f_names, ccl_res_YIELD$gl_gg$data$feature,
                                       ccl_res_YIELD$gl_gg$data$value, colors, ccl_res_y_f_types)
ccl_res_y_feature    = ccl_res_y_feature + coord_cartesian(xlim =c(0, 12.5))
# CCG-RES
ccg_res_y_f_names    = c(parse(text = 'Delta~Surface~Decomposition~Multiplier'), 'Soil pH',parse(text = 'Delta~Gross~Nitrification'),
                         'MAT','N Inputs', parse(text = 'Delta~Soil~Organic~Carbon'),
                         parse(text = 'Delta~Soil~Decomposition~Multiplier'), 'Residue Return Fraction', 'Clay Fraction',
                         'Soil Bulk Density', 'Sand Fraction', parse(text = 'Delta~N~Leached'),
                         parse(text='Delta~ET'),'MAP', 'Cover Crop Biomass') # reverse order
ccg_res_y_f_types    = c('State', 'Climate', 'State', 'State', 'Site', 'Site', 'Site',
                         'Management', 'State', 'State', 'Management', 'Climate', 'State', 'Site', 'State') # actual order
ccg_res_y_feature    = feature_yield_p(ccg_res_y_f_names, ccg_res_YIELD$gl_gg$data$feature,
                                       ccg_res_YIELD$gl_gg$data$value, colors, ccg_res_y_f_types)
# COMBINE
Yield_other_feature   = res_y_feature + ntill_res_y_feature + ccl_res_y_feature + ccg_res_y_feature +
  plot_layout(ncol = 2, nrow = 2, guides = 'collect') &
  theme(legend.position = 'bottom')
Yield_other_feature = Yield_other_feature + plot_annotation(tag_levels = list(c('a', 'c', 'b', 'd'))) &
  theme(plot.tag = element_text(size = 7))
ggsave(paste0(fig.path,'/ext-data/Yield-SHAP-features-other.pdf'), Yield_other_feature, units = 'mm', width = 180, height = 185, device='pdf', dpi=300)
#-------------------------------------------------------------------------------
# OLD - SAVE FOR NOW
#-------------------------------------------------------------------------------
# Figure(s) X (Global & Regional Cumulative Totals)
#-------------------------------------------------------------------------------
# GLOBAL
  # RESIDUE
res_gl_mCI   = fread(paste(data.path, 'residue-global-decadal-cumulative-responses.csv', sep = '/'))
res_gl_sum   = fread(paste(data.path, 'residue-global-decadal-cumulative-sums.csv', sep = '/'))
res_global   = cumul_line(res_gl_mCI, res_gl_sum)
ggsave(paste0(fig.path,'/', 'cumul-time-series-ghg-residue.tiff'),  res_global$GHG,  units = 'in', width = 8, height = 5, device='tiff', dpi=300)
  # NTILL
ntill_gl_mCI = fread(paste(data.path, 'ntill-global-decadal-cumulative-responses.csv', sep = '/'))
ntill_gl_sum = fread(paste(data.path, 'ntill-global-decadal-cumulative-sums.csv', sep = '/'))
ntill_global = cumul_line(ntill_gl_mCI, ntill_gl_sum)
ggsave(paste0(fig.path,'/', 'cumul-time-series-ghg-ntill.tiff'),  ntill_global$GHG,  units = 'in', width = 8, height = 5, device='tiff', dpi=300)
  # CCG
ccg_gl_mCI = fread(paste(data.path, 'ccg-global-decadal-cumulative-responses.csv', sep = '/'))
ccg_gl_sum = fread(paste(data.path, 'ccg-global-decadal-cumulative-sums.csv', sep = '/'))
ccg_global = cumul_line(ccg_gl_mCI, ccg_gl_sum)
ggsave(paste0(fig.path,'/', 'cumul-time-series-ghg-ccg.tiff'),  ccg_global$GHG,  units = 'in', width = 8, height = 5, device='tiff', dpi=300)
  # CCL
ccl_gl_mCI = fread(paste(data.path, 'ccl-global-decadal-cumulative-responses.csv', sep = '/'))
ccl_gl_sum = fread(paste(data.path, 'ccl-global-decadal-cumulative-sums.csv', sep = '/'))
ccl_global = cumul_line(ccl_gl_mCI, ccl_gl_sum)
ggsave(paste0(fig.path,'/', 'cumul-time-series-ghg-ccl.tiff'),  ccl_global$GHG,  units = 'in', width = 8, height = 5, device='tiff', dpi=300)
  # CCG-NTILL
ccg_ntill_gl_mCI = fread(paste(data.path, 'ccg_ntill-global-decadal-cumulative-responses.csv', sep = '/'))
ccg_ntill_gl_sum = fread(paste(data.path, 'ccg_ntill-global-decadal-cumulative-sums.csv', sep = '/'))
ccg_ntill_global = cumul_line_st(ccg_ntill_gl_mCI, ccg_ntill_gl_sum)
ggsave(paste0(fig.path,'/', 'cumul-time-series-ghg-ccg-ntill.tiff'),  ccg_ntill_global$GHG,  units = 'in', width = 8, height = 5, device='tiff', dpi=300)

ccl_ntill_gl_mCI = fread(paste(data.path, 'ccl_ntill-global-decadal-cumulative-responses.csv', sep = '/'))
ccl_ntill_gl_sum = fread(paste(data.path, 'ccl_ntill-global-decadal-cumulative-sums.csv', sep = '/'))
ccl_ntill_global = cumul_line_st(ccl_ntill_gl_mCI, ccl_ntill_gl_sum)
ggsave(paste0(fig.path,'/', 'cumul-time-series-ghg-ccl-ntill.tiff'),  ccl_ntill_global$GHG,  units = 'in', width = 8, height = 5, device='tiff', dpi=300)

# REGIONAL
res_rg_sum = fread(paste(data.path, 'residue-regional-decadal-cumulative-sums.csv', sep = '/'))
res_rg_mCI = fread(paste(data.path, 'residue-regional-decadal-cumulative-responses.csv', sep = '/'))
ntill_rg_sum = fread(paste(data.path, 'ntill-regional-decadal-cumulative-sums.csv', sep = '/'))
ntill_rg_mCI = fread(paste(data.path, 'ntill-regional-decadal-cumulative-responses.csv', sep = '/'))
ccg_rg_sum = fread(paste(data.path, 'ccg-regional-decadal-cumulative-sums.csv', sep = '/'))
ccg_rg_mCI = fread(paste(data.path, 'ccg-regional-decadal-cumulative-responses.csv', sep = '/'))
ccl_rg_sum = fread(paste(data.path, 'ccl-regional-decadal-cumulative-sums.csv', sep = '/'))
ccl_rg_mCI = fread(paste(data.path, 'ccl-regional-decadal-cumulative-responses.csv', sep = '/'))
ccg_ntill_rg_sum = fread(paste(data.path, 'ccg_ntill-regional-decadal-cumulative-sums.csv', sep = '/'))
ccg_ntill_rg_mCI = fread(paste(data.path, 'ccg_ntill-regional-decadal-cumulative-responses.csv', sep = '/'))
ccl_ntill_rg_sum = fread(paste(data.path, 'ccl_ntill-regional-decadal-cumulative-sums.csv', sep = '/'))
ccl_ntill_rg_mCI = fread(paste(data.path, 'ccl_ntill-regional-decadal-cumulative-responses.csv', sep = '/'))
  # ADP
res_region_ADP   = cumul_rg_line(res_rg_mCI[IPCC_NAME %in% 'ADP'], res_rg_sum[IPCC_NAME %in% 'ADP'])
ntill_region_ADP = cumul_rg_line(ntill_rg_mCI[IPCC_NAME %in% 'ADP'], ntill_rg_sum[IPCC_NAME %in% 'ADP'])
ccg_region_ADP   = cumul_rg_line(ccg_rg_mCI[IPCC_NAME %in% 'ADP'], ccg_rg_sum[IPCC_NAME %in% 'ADP'])
ccl_region_ADP   = cumul_rg_line(ccl_rg_mCI[IPCC_NAME %in% 'ADP'], ccl_rg_sum[IPCC_NAME %in% 'ADP'])
ccg_ntill_region_ADP = cumul_rg_line(ccg_ntill_rg_mCI[IPCC_NAME %in% 'ADP'], ccg_ntill_rg_sum[IPCC_NAME %in% 'ADP'])
ccl_ntill_region_ADP = cumul_rg_line(ccl_ntill_rg_mCI[IPCC_NAME %in% 'ADP'], ccl_ntill_rg_sum[IPCC_NAME %in% 'ADP'])
ggsave(paste0(fig.path,'/', 'cumul-time-series-ghg-residue-ADP.tiff'),  res_region_ADP$GHG,  units = 'in', width = 8, height = 5, device='tiff', dpi=300)
ggsave(paste0(fig.path,'/', 'cumul-time-series-ghg-ntill-ADP.tiff'),  ntill_region_ADP$GHG,  units = 'in', width = 8, height = 5, device='tiff', dpi=300)
ggsave(paste0(fig.path,'/', 'cumul-time-series-ghg-ccg-ADP.tiff'),  ccg_region_ADP$GHG,  units = 'in', width = 8, height = 5, device='tiff', dpi=300)
ggsave(paste0(fig.path,'/', 'cumul-time-series-ghg-ccl-ADP.tiff'),  ccl_region_ADP$GHG,  units = 'in', width = 8, height = 5, device='tiff', dpi=300)
ggsave(paste0(fig.path,'/', 'cumul-time-series-ghg-ccg-ntill-ADP.tiff'),  ccg_ntill_region_ADP$GHG,  units = 'in', width = 8, height = 5, device='tiff', dpi=300)
ggsave(paste0(fig.path,'/', 'cumul-time-series-ghg-ccl-ntill-ADP.tiff'),  ccl_ntill_region_ADP$GHG,  units = 'in', width = 8, height = 5, device='tiff', dpi=300)
  # AME
res_region_AME   = cumul_rg_line(res_rg_mCI[IPCC_NAME %in% 'AME'], res_rg_sum[IPCC_NAME %in% 'AME'])
ntill_region_AME = cumul_rg_line(ntill_rg_mCI[IPCC_NAME %in% 'AME'], ntill_rg_sum[IPCC_NAME %in% 'AME'])
ccg_region_AME   = cumul_rg_line(ccg_rg_mCI[IPCC_NAME %in% 'AME'], ccg_rg_sum[IPCC_NAME %in% 'AME'])
ccl_region_AME   = cumul_rg_line(ccl_rg_mCI[IPCC_NAME %in% 'AME'], ccl_rg_sum[IPCC_NAME %in% 'AME'])
ccg_ntill_region_AME = cumul_rg_line(ccg_ntill_rg_mCI[IPCC_NAME %in% 'AME'], ccg_ntill_rg_sum[IPCC_NAME %in% 'AME'])
ccl_ntill_region_AME = cumul_rg_line(ccl_ntill_rg_mCI[IPCC_NAME %in% 'AME'], ccl_ntill_rg_sum[IPCC_NAME %in% 'AME'])
ggsave(paste0(fig.path,'/', 'cumul-time-series-ghg-residue-AME.tiff'),  res_region_AME$GHG,  units = 'in', width = 8, height = 5, device='tiff', dpi=300)
ggsave(paste0(fig.path,'/', 'cumul-time-series-ghg-ntill-AME.tiff'),  ntill_region_AME$GHG,  units = 'in', width = 8, height = 5, device='tiff', dpi=300)
ggsave(paste0(fig.path,'/', 'cumul-time-series-ghg-ccg-AME.tiff'),  ccg_region_AME$GHG,  units = 'in', width = 8, height = 5, device='tiff', dpi=300)
ggsave(paste0(fig.path,'/', 'cumul-time-series-ghg-ccl-AME.tiff'),  ccl_region_AME$GHG,  units = 'in', width = 8, height = 5, device='tiff', dpi=300)
ggsave(paste0(fig.path,'/', 'cumul-time-series-ghg-ccg-ntill-AME.tiff'),  ccg_ntill_region_AME$GHG,  units = 'in', width = 8, height = 5, device='tiff', dpi=300)
ggsave(paste0(fig.path,'/', 'cumul-time-series-ghg-ccl-ntill-AME.tiff'),  ccl_ntill_region_AME$GHG,  units = 'in', width = 8, height = 5, device='tiff', dpi=300)
  # DEV
res_region_DEV   = cumul_rg_line(res_rg_mCI[IPCC_NAME %in% 'DEV'], res_rg_sum[IPCC_NAME %in% 'DEV'])
ntill_region_DEV = cumul_rg_line(ntill_rg_mCI[IPCC_NAME %in% 'DEV'], ntill_rg_sum[IPCC_NAME %in% 'DEV'])
ccg_region_DEV   = cumul_rg_line(ccg_rg_mCI[IPCC_NAME %in% 'DEV'], ccg_rg_sum[IPCC_NAME %in% 'DEV'])
ccl_region_DEV   = cumul_rg_line(ccl_rg_mCI[IPCC_NAME %in% 'DEV'], ccl_rg_sum[IPCC_NAME %in% 'DEV'])
ccg_ntill_region_DEV = cumul_rg_line(ccg_ntill_rg_mCI[IPCC_NAME %in% 'DEV'], ccg_ntill_rg_sum[IPCC_NAME %in% 'DEV'])
ccl_ntill_region_DEV = cumul_rg_line(ccl_ntill_rg_mCI[IPCC_NAME %in% 'DEV'], ccl_ntill_rg_sum[IPCC_NAME %in% 'DEV'])
ggsave(paste0(fig.path,'/', 'cumul-time-series-ghg-residue-DEV.tiff'),  res_region_DEV$GHG,  units = 'in', width = 8, height = 5, device='tiff', dpi=300)
ggsave(paste0(fig.path,'/', 'cumul-time-series-ghg-ntill-DEV.tiff'),  ntill_region_DEV$GHG,  units = 'in', width = 8, height = 5, device='tiff', dpi=300)
ggsave(paste0(fig.path,'/', 'cumul-time-series-ghg-ccg-DEV.tiff'),  ccg_region_DEV$GHG,  units = 'in', width = 8, height = 5, device='tiff', dpi=300)
ggsave(paste0(fig.path,'/', 'cumul-time-series-ghg-ccl-DEV.tiff'),  ccl_region_DEV$GHG,  units = 'in', width = 8, height = 5, device='tiff', dpi=300)
ggsave(paste0(fig.path,'/', 'cumul-time-series-ghg-ccg-ntill-DEV.tiff'),  ccg_ntill_region_DEV$GHG,  units = 'in', width = 8, height = 5, device='tiff', dpi=300)
ggsave(paste0(fig.path,'/', 'cumul-time-series-ghg-ccl-ntill-DEV.tiff'),  ccl_ntill_region_DEV$GHG,  units = 'in', width = 8, height = 5, device='tiff', dpi=300)
  # EEWCA
res_region_EEWCA   = cumul_rg_line(res_rg_mCI[IPCC_NAME %in% 'EEWCA'], res_rg_sum[IPCC_NAME %in% 'EEWCA'])
ntill_region_EEWCA = cumul_rg_line(ntill_rg_mCI[IPCC_NAME %in% 'EEWCA'], ntill_rg_sum[IPCC_NAME %in% 'EEWCA'])
ccg_region_EEWCA   = cumul_rg_line(ccg_rg_mCI[IPCC_NAME %in% 'EEWCA'], ccg_rg_sum[IPCC_NAME %in% 'EEWCA'])
ccl_region_EEWCA   = cumul_rg_line(ccl_rg_mCI[IPCC_NAME %in% 'EEWCA'], ccl_rg_sum[IPCC_NAME %in% 'EEWCA'])
ccg_ntill_region_EEWCA = cumul_rg_line(ccg_ntill_rg_mCI[IPCC_NAME %in% 'EEWCA'], ccg_ntill_rg_sum[IPCC_NAME %in% 'EEWCA'])
ccl_ntill_region_EEWCA = cumul_rg_line(ccl_ntill_rg_mCI[IPCC_NAME %in% 'EEWCA'], ccl_ntill_rg_sum[IPCC_NAME %in% 'EEWCA'])
ggsave(paste0(fig.path,'/', 'cumul-time-series-ghg-residue-EEWCA.tiff'),  res_region_EEWCA$GHG,  units = 'in', width = 8, height = 5, device='tiff', dpi=300)
ggsave(paste0(fig.path,'/', 'cumul-time-series-ghg-ntill-EEWCA.tiff'),  ntill_region_EEWCA$GHG,  units = 'in', width = 8, height = 5, device='tiff', dpi=300)
ggsave(paste0(fig.path,'/', 'cumul-time-series-ghg-ccg-EEWCA.tiff'),  ccg_region_EEWCA$GHG,  units = 'in', width = 8, height = 5, device='tiff', dpi=300)
ggsave(paste0(fig.path,'/', 'cumul-time-series-ghg-ccl-EEWCA.tiff'),  ccl_region_EEWCA$GHG,  units = 'in', width = 8, height = 5, device='tiff', dpi=300)
ggsave(paste0(fig.path,'/', 'cumul-time-series-ghg-ccg-ntill-EEWCA.tiff'),  ccg_ntill_region_EEWCA$GHG,  units = 'in', width = 8, height = 5, device='tiff', dpi=300)
ggsave(paste0(fig.path,'/', 'cumul-time-series-ghg-ccl-ntill-EEWCA.tiff'),  ccl_ntill_region_EEWCA$GHG,  units = 'in', width = 8, height = 5, device='tiff', dpi=300)
  # LAC
res_region_LAC   = cumul_rg_line(res_rg_mCI[IPCC_NAME %in% 'LAC'], res_rg_sum[IPCC_NAME %in% 'LAC'])
ntill_region_LAC = cumul_rg_line(ntill_rg_mCI[IPCC_NAME %in% 'LAC'], ntill_rg_sum[IPCC_NAME %in% 'LAC'])
ccg_region_LAC   = cumul_rg_line(ccg_rg_mCI[IPCC_NAME %in% 'LAC'], ccg_rg_sum[IPCC_NAME %in% 'LAC'])
ccl_region_LAC   = cumul_rg_line(ccl_rg_mCI[IPCC_NAME %in% 'LAC'], ccl_rg_sum[IPCC_NAME %in% 'LAC'])
ccg_ntill_region_LAC = cumul_rg_line(ccg_ntill_rg_mCI[IPCC_NAME %in% 'LAC'], ccg_ntill_rg_sum[IPCC_NAME %in% 'LAC'])
ccl_ntill_region_LAC = cumul_rg_line(ccl_ntill_rg_mCI[IPCC_NAME %in% 'LAC'], ccl_ntill_rg_sum[IPCC_NAME %in% 'LAC'])
ggsave(paste0(fig.path,'/', 'cumul-time-series-ghg-residue-LAC.tiff'),  res_region_LAC$GHG,  units = 'in', width = 8, height = 5, device='tiff', dpi=300)
ggsave(paste0(fig.path,'/', 'cumul-time-series-ghg-ntill-LAC.tiff'),  ntill_region_LAC$GHG,  units = 'in', width = 8, height = 5, device='tiff', dpi=300)
ggsave(paste0(fig.path,'/', 'cumul-time-series-ghg-ccg-LAC.tiff'),  ccg_region_LAC$GHG,  units = 'in', width = 8, height = 5, device='tiff', dpi=300)
ggsave(paste0(fig.path,'/', 'cumul-time-series-ghg-ccl-LAC.tiff'),  ccl_region_LAC$GHG,  units = 'in', width = 8, height = 5, device='tiff', dpi=300)
ggsave(paste0(fig.path,'/', 'cumul-time-series-ghg-ccg-ntill-LAC.tiff'),  ccg_ntill_region_LAC$GHG,  units = 'in', width = 8, height = 5, device='tiff', dpi=300)
ggsave(paste0(fig.path,'/', 'cumul-time-series-ghg-ccl-ntill-LAC.tiff'),  ccl_ntill_region_LAC$GHG,  units = 'in', width = 8, height = 5, device='tiff', dpi=300)
#-------------------------------------------------------------------------------
# Figure(s) X (Global & Regional Annual Totals)
#-------------------------------------------------------------------------------
# GLOBAL
gl_sum       = fread(paste(data.path, 'all-global-decadal-cumulative-sums.csv', sep = '/'))
ann_p_gl_370 = annual_total_gl_bp(gl_sum, 'ssp370')
ggsave(paste0(fig.path,'/', 'annual-ssp370-N2O-boxplot-global.tiff'),  ann_p_gl_370$N2O,  
       units = 'in', width = 10, height = 10, device='tiff', dpi=300)
ggsave(paste0(fig.path,'/', 'annual-ssp370-SOC-boxplot-global.tiff'),  ann_p_gl_370$SOC,  
       units = 'in', width = 10, height = 10, device='tiff', dpi=300)
ggsave(paste0(fig.path,'/', 'annual-ssp370-GHG-boxplot-global.tiff'),  ann_p_gl_370$GHG,  
       units = 'in', width = 10, height = 10, device='tiff', dpi=300)
ggsave(paste0(fig.path,'/', 'annual-ssp370-yield-boxplot-global.tiff'),  ann_p_gl_370$Yield,  
       units = 'in', width = 10, height = 10, device='tiff', dpi=300)
ggsave(paste0(fig.path,'/', 'annual-ssp370-residue-boxplot-global.tiff'),  ann_p_gl_370$Resid,  
       units = 'in', width = 10, height = 10, device='tiff', dpi=300)
ggsave(paste0(fig.path,'/', 'annual-ssp370-ccbiomass-boxplot-global.tiff'),  ann_p_gl_370$CC_bio,  
       units = 'in', width = 10, height = 10, device='tiff', dpi=300)

ann_p_gl_126 = annual_total_gl_bp(gl_sum, 'ssp126')
ggsave(paste0(fig.path,'/', 'annual-ssp126-N2O-boxplot-global.tiff'),  ann_p_gl_126$N2O,  
       units = 'in', width = 10, height = 10, device='tiff', dpi=300)
ggsave(paste0(fig.path,'/', 'annual-ssp126-SOC-boxplot-global.tiff'),  ann_p_gl_126$SOC,  
       units = 'in', width = 10, height = 10, device='tiff', dpi=300)
ggsave(paste0(fig.path,'/', 'annual-ssp126-GHG-boxplot-global.tiff'),  ann_p_gl_126$GHG,  
       units = 'in', width = 10, height = 10, device='tiff', dpi=300)
ggsave(paste0(fig.path,'/', 'annual-ssp126-yield-boxplot-global.tiff'),  ann_p_gl_126$Yield,  
       units = 'in', width = 10, height = 10, device='tiff', dpi=300)
ggsave(paste0(fig.path,'/', 'annual-ssp126-residue-boxplot-global.tiff'),  ann_p_gl_126$Resid,  
       units = 'in', width = 10, height = 10, device='tiff', dpi=300)
ggsave(paste0(fig.path,'/', 'annual-ssp126-ccbiomass-boxplot-global.tiff'),  ann_p_gl_126$CC_bio,  
       units = 'in', width = 10, height = 10, device='tiff', dpi=300)
# REGIONAL
rg_sum        = fread(paste(data.path, 'all-regional-decadal-cumulative-sums.csv', sep = '/'))
  # ADP
ADP_ann_p_370 = annual_total_rg_bp(rg_sum, 'ssp370', 'ADP')
ggsave(paste0(fig.path,'/', 'annual-ssp370-N2O-boxplot-ADP.tiff'),  ADP_ann_p_370$N2O,  
       units = 'in', width = 10, height = 10, device='tiff', dpi=300)
ggsave(paste0(fig.path,'/', 'annual-ssp370-SOC-boxplot-ADP.tiff'),  ADP_ann_p_370$SOC,  
       units = 'in', width = 10, height = 10, device='tiff', dpi=300)
ggsave(paste0(fig.path,'/', 'annual-ssp370-GHG-boxplot-ADP.tiff'),  ADP_ann_p_370$GHG,  
       units = 'in', width = 10, height = 10, device='tiff', dpi=300)
ggsave(paste0(fig.path,'/', 'annual-ssp370-yield-boxplot-ADP.tiff'),  ADP_ann_p_370$Yield,  
       units = 'in', width = 10, height = 10, device='tiff', dpi=300)
ggsave(paste0(fig.path,'/', 'annual-ssp370-residue-boxplot-ADP.tiff'),  ADP_ann_p_370$Resid,  
       units = 'in', width = 10, height = 10, device='tiff', dpi=300)
ggsave(paste0(fig.path,'/', 'annual-ssp370-ccbiomass-boxplot-ADP.tiff'),  ADP_ann_p_370$CC_bio,  
       units = 'in', width = 10, height = 10, device='tiff', dpi=300)
ADP_ann_p_126 = annual_total_rg_bp(rg_sum, 'ssp126', 'ADP')
ggsave(paste0(fig.path,'/', 'annual-ssp126-N2O-boxplot-ADP.tiff'),  ADP_ann_p_126$N2O,  
       units = 'in', width = 10, height = 10, device='tiff', dpi=300)
ggsave(paste0(fig.path,'/', 'annual-ssp126-SOC-boxplot-ADP.tiff'),  ADP_ann_p_126$SOC,  
       units = 'in', width = 10, height = 10, device='tiff', dpi=300)
ggsave(paste0(fig.path,'/', 'annual-ssp126-GHG-boxplot-ADP.tiff'),  ADP_ann_p_126$GHG,  
       units = 'in', width = 10, height = 10, device='tiff', dpi=300)
ggsave(paste0(fig.path,'/', 'annual-ssp126-yield-boxplot-ADP.tiff'),  ADP_ann_p_126$Yield,  
       units = 'in', width = 10, height = 10, device='tiff', dpi=300)
ggsave(paste0(fig.path,'/', 'annual-ssp126-residue-boxplot-ADP.tiff'),  ADP_ann_p_126$Resid,  
       units = 'in', width = 10, height = 10, device='tiff', dpi=300)
ggsave(paste0(fig.path,'/', 'annual-ssp126-ccbiomass-boxplot-ADP.tiff'),  ADP_ann_p_126$CC_bio,  
       units = 'in', width = 10, height = 10, device='tiff', dpi=300)
  # AME
AME_ann_p_370 = annual_total_rg_bp(rg_sum, 'ssp370', 'AME')
ggsave(paste0(fig.path,'/', 'annual-ssp370-N2O-boxplot-AME.tiff'),  AME_ann_p_370$N2O,  
       units = 'in', width = 10, height = 10, device='tiff', dpi=300)
ggsave(paste0(fig.path,'/', 'annual-ssp370-SOC-boxplot-AME.tiff'),  AME_ann_p_370$SOC,  
       units = 'in', width = 10, height = 10, device='tiff', dpi=300)
ggsave(paste0(fig.path,'/', 'annual-ssp370-GHG-boxplot-AME.tiff'),  AME_ann_p_370$GHG,  
       units = 'in', width = 10, height = 10, device='tiff', dpi=300)
ggsave(paste0(fig.path,'/', 'annual-ssp370-yield-boxplot-AME.tiff'),  AME_ann_p_370$Yield,  
       units = 'in', width = 10, height = 10, device='tiff', dpi=300)
ggsave(paste0(fig.path,'/', 'annual-ssp370-residue-boxplot-AME.tiff'),  AME_ann_p_370$Resid,  
       units = 'in', width = 10, height = 10, device='tiff', dpi=300)
ggsave(paste0(fig.path,'/', 'annual-ssp370-ccbiomass-boxplot-AME.tiff'),  AME_ann_p_370$CC_bio,  
       units = 'in', width = 10, height = 10, device='tiff', dpi=300)
AME_ann_p_126 = annual_total_rg_bp(rg_sum, 'ssp126', 'AME')
ggsave(paste0(fig.path,'/', 'annual-ssp126-N2O-boxplot-AME.tiff'),  AME_ann_p_126$N2O,  
       units = 'in', width = 10, height = 10, device='tiff', dpi=300)
ggsave(paste0(fig.path,'/', 'annual-ssp126-SOC-boxplot-AME.tiff'),  AME_ann_p_126$SOC,  
       units = 'in', width = 10, height = 10, device='tiff', dpi=300)
ggsave(paste0(fig.path,'/', 'annual-ssp126-GHG-boxplot-AME.tiff'),  AME_ann_p_126$GHG,  
       units = 'in', width = 10, height = 10, device='tiff', dpi=300)
ggsave(paste0(fig.path,'/', 'annual-ssp126-yield-boxplot-AME.tiff'),  AME_ann_p_126$Yield,  
       units = 'in', width = 10, height = 10, device='tiff', dpi=300)
ggsave(paste0(fig.path,'/', 'annual-ssp126-residue-boxplot-AME.tiff'),  AME_ann_p_126$Resid,  
       units = 'in', width = 10, height = 10, device='tiff', dpi=300)
ggsave(paste0(fig.path,'/', 'annual-ssp126-ccbiomass-boxplot-AME.tiff'),  AME_ann_p_126$CC_bio,  
       units = 'in', width = 10, height = 10, device='tiff', dpi=300)
  # DEV
DEV_ann_p_370 = annual_total_rg_bp(rg_sum, 'ssp370', 'DEV')
ggsave(paste0(fig.path,'/', 'annual-ssp370-N2O-boxplot-DEV.tiff'),  DEV_ann_p_370$N2O,  
       units = 'in', width = 10, height = 10, device='tiff', dpi=300)
ggsave(paste0(fig.path,'/', 'annual-ssp370-SOC-boxplot-DEV.tiff'),  DEV_ann_p_370$SOC,  
       units = 'in', width = 10, height = 10, device='tiff', dpi=300)
ggsave(paste0(fig.path,'/', 'annual-ssp370-GHG-boxplot-DEV.tiff'),  DEV_ann_p_370$GHG,  
       units = 'in', width = 10, height = 10, device='tiff', dpi=300)
ggsave(paste0(fig.path,'/', 'annual-ssp370-yield-boxplot-DEV.tiff'),  DEV_ann_p_370$Yield,  
       units = 'in', width = 10, height = 10, device='tiff', dpi=300)
ggsave(paste0(fig.path,'/', 'annual-ssp370-residue-boxplot-DEV.tiff'),  DEV_ann_p_370$Resid,  
       units = 'in', width = 10, height = 10, device='tiff', dpi=300)
ggsave(paste0(fig.path,'/', 'annual-ssp370-ccbiomass-boxplot-DEV.tiff'),  DEV_ann_p_370$CC_bio,  
       units = 'in', width = 10, height = 10, device='tiff', dpi=300)
DEV_ann_p_126 = annual_total_rg_bp(rg_sum, 'ssp126', 'DEV')
ggsave(paste0(fig.path,'/', 'annual-ssp126-N2O-boxplot-DEV.tiff'),  DEV_ann_p_126$N2O,  
       units = 'in', width = 10, height = 10, device='tiff', dpi=300)
ggsave(paste0(fig.path,'/', 'annual-ssp126-SOC-boxplot-DEV.tiff'),  DEV_ann_p_126$SOC,  
       units = 'in', width = 10, height = 10, device='tiff', dpi=300)
ggsave(paste0(fig.path,'/', 'annual-ssp126-GHG-boxplot-DEV.tiff'),  DEV_ann_p_126$GHG,  
       units = 'in', width = 10, height = 10, device='tiff', dpi=300)
ggsave(paste0(fig.path,'/', 'annual-ssp126-yield-boxplot-DEV.tiff'),  DEV_ann_p_126$Yield,  
       units = 'in', width = 10, height = 10, device='tiff', dpi=300)
ggsave(paste0(fig.path,'/', 'annual-ssp126-residue-boxplot-DEV.tiff'),  DEV_ann_p_126$Resid,  
       units = 'in', width = 10, height = 10, device='tiff', dpi=300)
ggsave(paste0(fig.path,'/', 'annual-ssp126-ccbiomass-boxplot-DEV.tiff'),  DEV_ann_p_126$CC_bio,  
       units = 'in', width = 10, height = 10, device='tiff', dpi=300)
  # EEWCA
EEWCA_ann_p_370 = annual_total_rg_bp(rg_sum, 'ssp370', 'EEWCA')
ggsave(paste0(fig.path,'/', 'annual-ssp370-N2O-boxplot-EEWCA.tiff'),  EEWCA_ann_p_370$N2O,  
       units = 'in', width = 10, height = 10, device='tiff', dpi=300)
ggsave(paste0(fig.path,'/', 'annual-ssp370-SOC-boxplot-EEWCA.tiff'),  EEWCA_ann_p_370$SOC,  
       units = 'in', width = 10, height = 10, device='tiff', dpi=300)
ggsave(paste0(fig.path,'/', 'annual-ssp370-GHG-boxplot-EEWCA.tiff'),  EEWCA_ann_p_370$GHG,  
       units = 'in', width = 10, height = 10, device='tiff', dpi=300)
ggsave(paste0(fig.path,'/', 'annual-ssp370-yield-boxplot-EEWCA.tiff'),  EEWCA_ann_p_370$Yield,  
       units = 'in', width = 10, height = 10, device='tiff', dpi=300)
ggsave(paste0(fig.path,'/', 'annual-ssp370-residue-boxplot-EEWCA.tiff'),  EEWCA_ann_p_370$Resid,  
       units = 'in', width = 10, height = 10, device='tiff', dpi=300)
ggsave(paste0(fig.path,'/', 'annual-ssp370-ccbiomass-boxplot-EEWCA.tiff'),  EEWCA_ann_p_370$CC_bio,  
       units = 'in', width = 10, height = 10, device='tiff', dpi=300)
EEWCA_ann_p_126 = annual_total_rg_bp(rg_sum, 'ssp126', 'EEWCA')
ggsave(paste0(fig.path,'/', 'annual-ssp126-N2O-boxplot-EEWCA.tiff'),  EEWCA_ann_p_126$N2O,  
       units = 'in', width = 10, height = 10, device='tiff', dpi=300)
ggsave(paste0(fig.path,'/', 'annual-ssp126-SOC-boxplot-EEWCA.tiff'),  EEWCA_ann_p_126$SOC,  
       units = 'in', width = 10, height = 10, device='tiff', dpi=300)
ggsave(paste0(fig.path,'/', 'annual-ssp126-GHG-boxplot-EEWCA.tiff'),  EEWCA_ann_p_126$GHG,  
       units = 'in', width = 10, height = 10, device='tiff', dpi=300)
ggsave(paste0(fig.path,'/', 'annual-ssp126-yield-boxplot-EEWCA.tiff'),  EEWCA_ann_p_126$Yield,  
       units = 'in', width = 10, height = 10, device='tiff', dpi=300)
ggsave(paste0(fig.path,'/', 'annual-ssp126-residue-boxplot-EEWCA.tiff'),  EEWCA_ann_p_126$Resid,  
       units = 'in', width = 10, height = 10, device='tiff', dpi=300)
ggsave(paste0(fig.path,'/', 'annual-ssp126-ccbiomass-boxplot-EEWCA.tiff'),  EEWCA_ann_p_126$CC_bio,  
       units = 'in', width = 10, height = 10, device='tiff', dpi=300)
  # LAC
LAC_ann_p_370 = annual_total_rg_bp(rg_sum, 'ssp370', 'LAC')
ggsave(paste0(fig.path,'/', 'annual-ssp370-N2O-boxplot-LAC.tiff'),  LAC_ann_p_370$N2O,  
       units = 'in', width = 10, height = 10, device='tiff', dpi=300)
ggsave(paste0(fig.path,'/', 'annual-ssp370-SOC-boxplot-LAC.tiff'),  LAC_ann_p_370$SOC,  
       units = 'in', width = 10, height = 10, device='tiff', dpi=300)
ggsave(paste0(fig.path,'/', 'annual-ssp370-GHG-boxplot-LAC.tiff'),  LAC_ann_p_370$GHG,  
       units = 'in', width = 10, height = 10, device='tiff', dpi=300)
ggsave(paste0(fig.path,'/', 'annual-ssp370-yield-boxplot-LAC.tiff'),  LAC_ann_p_370$Yield,  
       units = 'in', width = 10, height = 10, device='tiff', dpi=300)
ggsave(paste0(fig.path,'/', 'annual-ssp370-residue-boxplot-LAC.tiff'),  LAC_ann_p_370$Resid,  
       units = 'in', width = 10, height = 10, device='tiff', dpi=300)
ggsave(paste0(fig.path,'/', 'annual-ssp370-ccbiomass-boxplot-LAC.tiff'),  LAC_ann_p_370$CC_bio,  
       units = 'in', width = 10, height = 10, device='tiff', dpi=300)
LAC_ann_p_126 = annual_total_rg_bp(rg_sum, 'ssp126', 'LAC')
ggsave(paste0(fig.path,'/', 'annual-ssp126-N2O-boxplot-LAC.tiff'),  LAC_ann_p_126$N2O,  
       units = 'in', width = 10, height = 10, device='tiff', dpi=300)
ggsave(paste0(fig.path,'/', 'annual-ssp126-SOC-boxplot-LAC.tiff'),  LAC_ann_p_126$SOC,  
       units = 'in', width = 10, height = 10, device='tiff', dpi=300)
ggsave(paste0(fig.path,'/', 'annual-ssp126-GHG-boxplot-LAC.tiff'),  LAC_ann_p_126$GHG,  
       units = 'in', width = 10, height = 10, device='tiff', dpi=300)
ggsave(paste0(fig.path,'/', 'annual-ssp126-yield-boxplot-LAC.tiff'),  LAC_ann_p_126$Yield,  
       units = 'in', width = 10, height = 10, device='tiff', dpi=300)
ggsave(paste0(fig.path,'/', 'annual-ssp126-residue-boxplot-LAC.tiff'),  LAC_ann_p_126$Resid,  
       units = 'in', width = 10, height = 10, device='tiff', dpi=300)
ggsave(paste0(fig.path,'/', 'annual-ssp126-ccbiomass-boxplot-LAC.tiff'),  LAC_ann_p_126$CC_bio,  
       units = 'in', width = 10, height = 10, device='tiff', dpi=300)
#-------------------------------------------------------------------------------
# Figure(s) X (Explanatory Variable Maps)
#-------------------------------------------------------------------------------
# 2100
ntill_ann_2100     = ann_var_map(rel_ntill_flux_dt, '2100', 'ssp370', c('m_sC.N', 'm_sldcmp', 'm_sfdcmp', 'm_annet', 'm_cr_NPP'))
ntill_ann_2100$NPP
ntill_ann_2100$annet
ntill_ann_2100$sldcmp
# grass cover crop
ccg_ann_2100       = ann_var_map(rel_ccg_flux_dt, '2100', 'ssp370', c('m_sC.N', 'm_sldcmp', 'm_sfdcmp', 'm_annet', 'm_cr_NPP'))
ccg_ann_2100$NPP
ccg_ann_2100$sldcmp
ccg_ann_2100$sfdcmp
ccg_ann_2100$annet
ccg_ann_2100$annet_legend
ccg_ann_2100$dcmp_legend
ccg_ann_2100$NPP_legend
ccg_ntill_ann_2100 = ann_var_map(rel_ccg_ntill_flux_dt, '2100', 'ssp370', c('m_sC.N', 'm_sldcmp', 'm_sfdcmp', 'm_annet', 'm_cr_NPP'))
ccg_ntill_ann_2100$NPP
ccg_ntill_ann_2100$sldcmp
ccg_ntill_ann_2100$sfdcmp
ccg_ntill_ann_2100$annet
# legume cover crop
ccl_ann_2100       = ann_var_map(rel_ccl_flux_dt, '2100', 'ssp370', c('m_sC.N', 'm_sldcmp', 'm_sfdcmp', 'm_annet', 'm_cr_NPP'))
ccl_ann_2100$NPP
ccl_ann_2100$sldcmp
ccl_ann_2100$sfdcmp
ccl_ann_2100$annet
ccl_ntill_ann_2100 = ann_var_map(rel_ccl_ntill_flux_dt, '2100', 'ssp370', c('m_sC.N', 'm_sldcmp', 'm_sfdcmp', 'm_annet', 'm_cr_NPP'))
ccl_ntill_ann_2100$NPP
ccl_ntill_ann_2100$sldcmp
ccl_ntill_ann_2100$sfdcmp
ccl_ntill_ann_2100$annet
# 2050
ntill_ann_2050     = ann_var_map(rel_ntill_flux_dt, '2050', 'ssp370', c('m_sC.N', 'm_sldcmp', 'm_sfdcmp', 'm_annet', 'm_cr_NPP'))
ntill_ann_2050$NPP
ntill_ann_2050$annet
# grass cover crop
ccg_ann_2050       = ann_var_map(rel_ccg_flux_dt, '2050', 'ssp370', c('m_sC.N', 'm_sldcmp', 'm_sfdcmp', 'm_annet', 'm_cr_NPP'))
ccg_ann_2050$NPP
ccg_ann_2050$sldcmp
ccg_ntill_ann_2050 = ann_var_map(rel_ccg_ntill_flux_dt, '2050', 'ssp370', c('m_sC.N', 'm_sldcmp', 'm_sfdcmp', 'm_annet', 'm_cr_NPP'))
ccg_ntill_ann_2050$NPP
ccg_ntill_ann_2050$sldcmp
# legume cover crop
ccl_ann_2050       = ann_var_map(rel_ccl_flux_dt, '2050', 'ssp370', c('m_sC.N', 'm_sldcmp', 'm_sfdcmp', 'm_annet', 'm_cr_NPP'))
ccl_ann_2050$NPP
ccl_ann_2050$sldcmp
ccl_ann_2050$annet
ccl_ntill_ann_2050 = ann_var_map(rel_ccl_ntill_flux_dt, '2050', 'ssp370', c('m_sC.N', 'm_sldcmp', 'm_sfdcmp', 'm_annet', 'm_cr_NPP'))
ccl_ntill_ann_2050$NPP
ccl_ntill_ann_2050$sldcmp
ccl_ntill_ann_2050$annet
#-------------------------------------------------------------------------------
# LOAD Global, Regional, Country Means
#-------------------------------------------------------------------------------
d_r_mSE # load here
#-------------------------------------------------------------------------------
# Figure(s) X (Regional Emissions + Yield Responses)
#-------------------------------------------------------------------------------
# BASE MAP
map_base        = IPCC_map(crop_area_dt, lu.path)
ggsave(paste0(fig.path,'/IPCC-map.tiff'),    map_base$IPCC,  units = 'in', width = 9, 
       height = 5, device='tiff', dpi=300)

# SINGLE PRACTICE
ADP_sgl_370_mSE = rg_annual_sgl_mSE(d_r_mSE, 'ssp370', c('res','ntill','ccg','ccl'), 
                                    c('ADP', 'AME', 'DEV', 'EEWCA', 'LAC'),
                                    'ADP')
ggsave(paste0(fig.path,'/ADP-ssp370-annual-GHG-response-2100.png'),    
       ADP_sgl_370_mSE$GHG, bg = 'transparent', units = 'in', width = 11, height = 4, 
       device='png', dpi=300)
ggsave(paste0(fig.path,'/ADP-ssp370-annual-yield-response-2100.png'),    
       ADP_sgl_370_mSE$Yield, bg = 'transparent', units = 'in', width = 11, height = 4, 
       device='png', dpi=300)
ADP_sgl_126_mSE = rg_annual_sgl_mSE(d_r_mSE, 'ssp126', c('res','ntill','ccg','ccl'), 
                                    c('ADP', 'AME', 'DEV', 'EEWCA', 'LAC'),
                                    'ADP')
ggsave(paste0(fig.path,'/ADP-ssp126-annual-GHG-response-2100.png'),    
       ADP_sgl_126_mSE$GHG, bg = 'transparent', units = 'in', width = 11, height = 4, 
       device='png', dpi=300)
ggsave(paste0(fig.path,'/ADP-ssp126-annual-yield-response-2100.png'),    
       ADP_sgl_126_mSE$Yield, bg = 'transparent', units = 'in', width = 11, height = 4, 
       device='png', dpi=300)
AME_sgl_370_mSE = rg_annual_sgl_mSE(d_r_mSE, 'ssp370', c('res','ntill','ccg','ccl'), 
                                    c('ADP', 'AME', 'DEV', 'EEWCA', 'LAC'),
                                    'AME')
ggsave(paste0(fig.path,'/AME-ssp370-annual-GHG-response-2100.png'),    
       AME_sgl_370_mSE$GHG, bg = 'transparent', units = 'in', width = 11, height = 4, 
       device='png', dpi=300)
ggsave(paste0(fig.path,'/AME-ssp370-annual-yield-response-2100.png'),    
       AME_sgl_370_mSE$Yield, bg = 'transparent', units = 'in', width = 11, height = 4, 
       device='png', dpi=300)
AME_sgl_126_mSE = rg_annual_sgl_mSE(d_r_mSE, 'ssp126', c('res','ntill','ccg','ccl'), 
                                    c('ADP', 'AME', 'DEV', 'EEWCA', 'LAC'),
                                    'AME')
ggsave(paste0(fig.path,'/AME-ssp126-annual-GHG-response-2100.png'),    
       AME_sgl_126_mSE$GHG, bg = 'transparent', units = 'in', width = 11, height = 4, 
       device='png', dpi=300)
ggsave(paste0(fig.path,'/AME-ssp126-annual-yield-response-2100.png'),    
       AME_sgl_126_mSE$Yield, bg = 'transparent', units = 'in', width = 11, height = 4, 
       device='png', dpi=300)
DEV_sgl_370_mSE = rg_annual_sgl_mSE(d_r_mSE, 'ssp370', c('res','ntill','ccg','ccl'), 
                                    c('ADP', 'AME', 'DEV', 'EEWCA', 'LAC'),
                                    'DEV')
ggsave(paste0(fig.path,'/DEV-ssp370-annual-GHG-response-2100.png'),    
       DEV_sgl_370_mSE$GHG, bg = 'transparent', units = 'in', width = 11, height = 4, 
       device='png', dpi=300)
ggsave(paste0(fig.path,'/DEV-ssp370-annual-yield-response-2100.png'),    
       DEV_sgl_370_mSE$Yield, bg = 'transparent', units = 'in', width = 11, height = 4, 
       device='png', dpi=300)
DEV_sgl_126_mSE = rg_annual_sgl_mSE(d_r_mSE, 'ssp126', c('res','ntill','ccg','ccl'), 
                                    c('ADP', 'AME', 'DEV', 'EEWCA', 'LAC'),
                                    'DEV')
ggsave(paste0(fig.path,'/DEV-ssp126-annual-GHG-response-2100.png'),    
       DEV_sgl_126_mSE$GHG, bg = 'transparent', units = 'in', width = 11, height = 4, 
       device='png', dpi=300)
ggsave(paste0(fig.path,'/DEV-ssp126-annual-yield-response-2100.png'),    
       DEV_sgl_126_mSE$Yield, bg = 'transparent', units = 'in', width = 11, height = 4, 
       device='png', dpi=300)
EEWCA_sgl_370_mSE = rg_annual_sgl_mSE(d_r_mSE, 'ssp370', c('res','ntill','ccg','ccl'), 
                                    c('ADP', 'AME', 'DEV', 'EEWCA', 'LAC'),
                                    'EEWCA')
ggsave(paste0(fig.path,'/EEWCA-ssp370-annual-GHG-response-2100.png'),    
       EEWCA_sgl_370_mSE$GHG, bg = 'transparent', units = 'in', width = 11, height = 4, 
       device='png', dpi=300)
ggsave(paste0(fig.path,'/EEWCA-ssp370-annual-yield-response-2100.png'),    
       EEWCA_sgl_370_mSE$Yield, bg = 'transparent', units = 'in', width = 11, height = 4, 
       device='png', dpi=300)
EEWCA_sgl_126_mSE = rg_annual_sgl_mSE(d_r_mSE, 'ssp126', c('res','ntill','ccg','ccl'), 
                                    c('ADP', 'AME', 'DEV', 'EEWCA', 'LAC'),
                                    'EEWCA')
ggsave(paste0(fig.path,'/EEWCA-ssp126-annual-GHG-response-2100.png'),    
       EEWCA_sgl_126_mSE$GHG, bg = 'transparent', units = 'in', width = 11, height = 4, 
       device='png', dpi=300)
ggsave(paste0(fig.path,'/EEWCA-ssp126-annual-yield-response-2100.png'),    
       EEWCA_sgl_126_mSE$Yield, bg = 'transparent', units = 'in', width = 11, height = 4, 
       device='png', dpi=300)
LAC_sgl_370_mSE = rg_annual_sgl_mSE(d_r_mSE, 'ssp370', c('res','ntill','ccg','ccl'), 
                                    c('ADP', 'AME', 'DEV', 'EEWCA', 'LAC'),
                                    'LAC')
ggsave(paste0(fig.path,'/LAC-ssp370-annual-GHG-response-2100.png'),    
       LAC_sgl_370_mSE$GHG, bg = 'transparent', units = 'in', width = 11, height = 4, 
       device='png', dpi=300)
ggsave(paste0(fig.path,'/LAC-ssp370-annual-yield-response-2100.png'),    
       LAC_sgl_370_mSE$Yield, bg = 'transparent', units = 'in', width = 11, height = 4, 
       device='png', dpi=300)
LAC_sgl_126_mSE = rg_annual_sgl_mSE(d_r_mSE, 'ssp126', c('res','ntill','ccg','ccl'), 
                                    c('ADP', 'AME', 'DEV', 'EEWCA', 'LAC'),
                                    'LAC')
ggsave(paste0(fig.path,'/LAC-ssp126-annual-GHG-response-2100.png'),    
       LAC_sgl_126_mSE$GHG, bg = 'transparent', units = 'in', width = 11, height = 4, 
       device='png', dpi=300)
ggsave(paste0(fig.path,'/LAC-ssp126-annual-yield-response-2100.png'),    
       LAC_sgl_126_mSE$Yield, bg = 'transparent', units = 'in', width = 11, height = 4, 
       device='png', dpi=300)
# COMBINED PRACTICE
ADP_cmb_370_mSE = rg_annual_cmb_mSE(d_r_mSE, 'ssp370', c('ccg-ntill','ccl-ntill'), 
                                    c('ADP', 'AME', 'DEV', 'EEWCA', 'LAC'),
                                    'ADP')
ggsave(paste0(fig.path,'/ADP-ssp370-annual-GHG-response-2100-combined.png'),    
       ADP_cmb_370_mSE$GHG, bg = 'transparent', units = 'in', width = 8, height = 4, 
       device='png', dpi=300)
ggsave(paste0(fig.path,'/ADP-ssp370-annual-yield-response-2100-combined.png'),    
       ADP_cmb_370_mSE$Yield, bg = 'transparent', units = 'in', width = 8, height = 4, 
       device='png', dpi=300)
ADP_cmb_126_mSE = rg_annual_cmb_mSE(d_r_mSE, 'ssp126', c('ccg-ntill','ccl-ntill'), 
                                    c('ADP', 'AME', 'DEV', 'EEWCA', 'LAC'),
                                    'ADP')
ggsave(paste0(fig.path,'/ADP-ssp126-annual-GHG-response-2100-combined.png'),    
       ADP_cmb_126_mSE$GHG, bg = 'transparent', units = 'in', width = 8, height = 4, 
       device='png', dpi=300)
ggsave(paste0(fig.path,'/ADP-ssp126-annual-yield-response-2100-combined.png'),    
       ADP_cmb_126_mSE$Yield, bg = 'transparent', units = 'in', width = 8, height = 4, 
       device='png', dpi=300)
AME_cmb_370_mSE = rg_annual_cmb_mSE(d_r_mSE, 'ssp370', c('ccg-ntill','ccl-ntill'), 
                                    c('ADP', 'AME', 'DEV', 'EEWCA', 'LAC'),
                                    'AME')
ggsave(paste0(fig.path,'/AME-ssp370-annual-GHG-response-2100-combined.png'),    
       AME_cmb_370_mSE$GHG, bg = 'transparent', units = 'in', width = 8, height = 4, 
       device='png', dpi=300)
ggsave(paste0(fig.path,'/AME-ssp370-annual-yield-response-2100-combined.png'),    
       AME_cmb_370_mSE$Yield, bg = 'transparent', units = 'in', width = 8, height = 4, 
       device='png', dpi=300)
AME_cmb_126_mSE = rg_annual_cmb_mSE(d_r_mSE, 'ssp126', c('ccg-ntill','ccl-ntill'), 
                                    c('ADP', 'AME', 'DEV', 'EEWCA', 'LAC'),
                                    'AME')
ggsave(paste0(fig.path,'/AME-ssp126-annual-GHG-response-2100-combined.png'),    
       AME_cmb_126_mSE$GHG, bg = 'transparent', units = 'in', width = 8, height = 4, 
       device='png', dpi=300)
ggsave(paste0(fig.path,'/AME-ssp126-annual-yield-response-2100-combined.png'),    
       AME_cmb_126_mSE$Yield, bg = 'transparent', units = 'in', width = 8, height = 4, 
       device='png', dpi=300)
DEV_cmb_370_mSE = rg_annual_cmb_mSE(d_r_mSE, 'ssp370', c('ccg-ntill','ccl-ntill'), 
                                    c('ADP', 'AME', 'DEV', 'EEWCA', 'LAC'),
                                    'DEV')
ggsave(paste0(fig.path,'/DEV-ssp370-annual-GHG-response-2100-combined.png'),    
       DEV_cmb_370_mSE$GHG, bg = 'transparent', units = 'in', width = 8, height = 4, 
       device='png', dpi=300)
ggsave(paste0(fig.path,'/DEV-ssp370-annual-yield-response-2100-combined.png'),    
       DEV_cmb_370_mSE$Yield, bg = 'transparent', units = 'in', width = 8, height = 4, 
       device='png', dpi=300)
DEV_cmb_126_mSE = rg_annual_cmb_mSE(d_r_mSE, 'ssp126', c('ccg-ntill','ccl-ntill'), 
                                    c('ADP', 'AME', 'DEV', 'EEWCA', 'LAC'),
                                    'DEV')
ggsave(paste0(fig.path,'/DEV-ssp126-annual-GHG-response-2100-combined.png'),    
       DEV_cmb_126_mSE$GHG, bg = 'transparent', units = 'in', width = 8, height = 4, 
       device='png', dpi=300)
ggsave(paste0(fig.path,'/DEV-ssp126-annual-yield-response-2100-combined.png'),    
       DEV_cmb_126_mSE$Yield, bg = 'transparent', units = 'in', width = 8, height = 4, 
       device='png', dpi=300)
EEWCA_cmb_370_mSE = rg_annual_cmb_mSE(d_r_mSE, 'ssp370', c('ccg-ntill','ccl-ntill'), 
                                    c('ADP', 'AME', 'DEV', 'EEWCA', 'LAC'),
                                    'EEWCA')
ggsave(paste0(fig.path,'/EEWCA-ssp370-annual-GHG-response-2100-combined.png'),    
       EEWCA_cmb_370_mSE$GHG, bg = 'transparent', units = 'in', width = 8, height = 4, 
       device='png', dpi=300)
ggsave(paste0(fig.path,'/EEWCA-ssp370-annual-yield-response-2100-combined.png'),    
       EEWCA_cmb_370_mSE$Yield, bg = 'transparent', units = 'in', width = 8, height = 4, 
       device='png', dpi=300)
EEWCA_cmb_126_mSE = rg_annual_cmb_mSE(d_r_mSE, 'ssp126', c('ccg-ntill','ccl-ntill'), 
                                    c('ADP', 'AME', 'DEV', 'EEWCA', 'LAC'),
                                    'EEWCA')
ggsave(paste0(fig.path,'/EEWCA-ssp126-annual-GHG-response-2100-combined.png'),    
       EEWCA_cmb_126_mSE$GHG, bg = 'transparent', units = 'in', width = 8, height = 4, 
       device='png', dpi=300)
ggsave(paste0(fig.path,'/EEWCA-ssp126-annual-yield-response-2100-combined.png'),    
       EEWCA_cmb_126_mSE$Yield, bg = 'transparent', units = 'in', width = 8, height = 4, 
       device='png', dpi=300)
LAC_cmb_370_mSE = rg_annual_cmb_mSE(d_r_mSE, 'ssp370', c('ccg-ntill','ccl-ntill'), 
                                    c('ADP', 'AME', 'DEV', 'EEWCA', 'LAC'),
                                    'LAC')
ggsave(paste0(fig.path,'/LAC-ssp370-annual-GHG-response-2100-combined.png'),    
       LAC_cmb_370_mSE$GHG, bg = 'transparent', units = 'in', width = 8, height = 4, 
       device='png', dpi=300)
ggsave(paste0(fig.path,'/LAC-ssp370-annual-yield-response-2100-combined.png'),    
       LAC_cmb_370_mSE$Yield, bg = 'transparent', units = 'in', width = 8, height = 4, 
       device='png', dpi=300)
LAC_cmb_126_mSE = rg_annual_cmb_mSE(d_r_mSE, 'ssp126', c('ccg-ntill','ccl-ntill'), 
                                    c('ADP', 'AME', 'DEV', 'EEWCA', 'LAC'),
                                    'LAC')
ggsave(paste0(fig.path,'/LAC-ssp126-annual-GHG-response-2100-combined.png'),    
       LAC_cmb_126_mSE$GHG, bg = 'transparent', units = 'in', width = 8, height = 4, 
       device='png', dpi=300)
ggsave(paste0(fig.path,'/LAC-ssp126-annual-yield-response-2100-combined.png'),    
       LAC_cmb_126_mSE$Yield, bg = 'transparent', units = 'in', width = 8, height = 4, 
       device='png', dpi=300)