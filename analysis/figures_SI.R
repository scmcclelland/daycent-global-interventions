# filename:    figures_main.R
# created:     06 March 2023
# updated:     23 July 2024
# author:      S.C. McClelland
# description: This file creates SI text figures for manuscript.
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
# FIGURE S1 | UNCONSTRAINED & YIELD CONSTRAINED POTENTIAL
#-------------------------------------------------------------------------------
# GLOBAL
biophys_2050     = fread(paste(data.path, 'global-unconstrained-GHG-mitigation-potential-2050.csv', sep = '/'))
yield_cst_2050   = fread(paste(data.path, 'global-yield-constrained-GHG-mitigation-potential-2050.csv', sep = '/'))
yield_cst_2050_l = fread(paste(data.path, 'global-yield-constrained-low-price-GHG-mitigation-potential-2050.csv', sep = '/'))
yield_cst_2050_h = fread(paste(data.path, 'global-yield-constrained-high-price-GHG-mitigation-potential-2050.csv', sep = '/'))
# REGIONAL
ipcc_biophys_2050     = fread(paste(data.path, 'regional-unconstrained-GHG-mitigation-potential-2050.csv', sep = '/'))
ipcc_yield_cst_2050   = fread(paste(data.path, 'regional-yield-constrained-GHG-mitigation-potential-2050.csv', sep = '/'))
ipcc_yield_cst_2050_l = fread(paste(data.path, 'regional-yield-constrained-low-price-GHG-mitigation-potential-2050.csv', sep = '/'))
ipcc_yield_cst_2050_h = fread(paste(data.path, 'regional-yield-constrained-high-price-GHG-mitigation-potential-2050.csv', sep = '/'))

# background map
map_background = IPCC_map(lu.path, shp, raster)
ggsave(paste0(fig.path,'/si/IPCC-map-background.pdf'), map_background$IPCC, units = 'mm', width = 180, height = 185, device='pdf', dpi=300)
# global
global_370_ptl = gl_bar_potential(biophys_2050, yield_cst_2050, yield_cst_2050_l, yield_cst_2050_h,'ssp370')
ggsave(paste0(fig.path,'/si/ssp370-2050-global-ghg-mitigation-potential.pdf'), global_370_ptl$GHG, bg = 'transparent',units = 'mm', width = 58, height = 58,  device='pdf', dpi=300)
# regional
ADP_370_ptl    = rg_bar_potential(ipcc_biophys_2050, ipcc_yield_cst_2050, ipcc_yield_cst_2050_l, ipcc_yield_cst_2050_h, 'ssp370', 'ADP')
ggsave(paste0(fig.path,'/si/ssp370-2050-ADP-ghg-mitigation-potential.pdf'), ADP_370_ptl$GHG, bg = 'transparent',units = 'mm', width = 58, height = 58,  device='pdf', dpi=300)
AME_370_ptl    = rg_bar_potential(ipcc_biophys_2050, ipcc_yield_cst_2050, ipcc_yield_cst_2050_l, ipcc_yield_cst_2050_h, 'ssp370', 'AME')
ggsave(paste0(fig.path,'/si/ssp370-2050-AME-ghg-mitigation-potential.pdf'), AME_370_ptl$GHG, bg = 'transparent',units = 'mm', width = 58, height = 58,  device='pdf', dpi=300)
DEV_370_ptl    = rg_bar_potential(ipcc_biophys_2050, ipcc_yield_cst_2050, ipcc_yield_cst_2050_l, ipcc_yield_cst_2050_h, 'ssp370', 'DEV')
ggsave(paste0(fig.path,'/si/ssp370-2050-DEV-ghg-mitigation-potential.pdf'), DEV_370_ptl$GHG, bg = 'transparent',units = 'mm', width = 58, height = 58,  device='pdf', dpi=300)
EEWCA_370_ptl  = rg_bar_potential(ipcc_biophys_2050, ipcc_yield_cst_2050, ipcc_yield_cst_2050_l, ipcc_yield_cst_2050_h, 'ssp370', 'EEWCA')
ggsave(paste0(fig.path,'/si/ssp370-2050-EEWCA-ghg-mitigation-potential.pdf'), EEWCA_370_ptl$GHG, bg = 'transparent',units = 'mm', width = 58, height = 58,  device='pdf', dpi=300)
LAC_370_ptl    = rg_bar_potential(ipcc_biophys_2050, ipcc_yield_cst_2050, ipcc_yield_cst_2050_l, ipcc_yield_cst_2050_h, 'ssp370', 'LAC')
ggsave(paste0(fig.path,'/si/ssp370-2050-LAC-ghg-mitigation-potential.pdf'), LAC_370_ptl$GHG, bg = 'transparent',units = 'mm', width = 58, height = 58,  device='pdf', dpi=300)
#-------------------------------------------------------------------------------
# FIGURE S2-4 | MAPS
#-------------------------------------------------------------------------------
Mg_ha = 100L
kg_ha = 10L
C_bio = 0.42 # Ma et al. 2018 | value for 'reproductive organs'
# LOAD DATA
# RES
res_dt = readRDS(paste(data.path, 'imputed-ensemble-relative-responses-weighted-mean-res.rds', sep = '/'))
res_dt = dcast(res_dt,
               cell + x + y + y_block + ssp + total_crop_area_ha ~ variable,
               value.var = 'value')
res_dt[, m_cr_grain := (m_cr_grain*kg_ha)/C_bio]
res_dt[, s_cr_grain := (s_cr_grain/Mg_ha)/C_bio]
# NTILL-RES
ntill_res_dt = readRDS(paste(data.path, 'imputed-ensemble-relative-responses-weighted-mean-ntill-res.rds', sep = '/'))
ntill_res_dt = dcast(ntill_res_dt,
                     cell + x + y + y_block + ssp + total_crop_area_ha ~ variable,
                     value.var = 'value')
ntill_res_dt[, m_cr_grain := (m_cr_grain*kg_ha)/C_bio]
ntill_res_dt[, s_cr_grain := (s_cr_grain/Mg_ha)/C_bio]
# CCG-RES
ccg_res_dt = readRDS(paste(data.path, 'imputed-ensemble-relative-responses-weighted-mean-ccg-res.rds', sep = '/'))
ccg_res_dt = dcast(ccg_res_dt,
                   cell + x + y + y_block + ssp + total_crop_area_ha ~ variable,
                   value.var = 'value')
ccg_res_dt[, m_cr_grain := (m_cr_grain*kg_ha)/C_bio]
ccg_res_dt[, s_cr_grain := (s_cr_grain/Mg_ha)/C_bio]
# CCL-RES
ccl_res_dt = readRDS(paste(data.path, 'imputed-ensemble-relative-responses-weighted-mean-ccl-res.rds', sep = '/'))
ccl_res_dt = dcast(ccl_res_dt,
                   cell + x + y + y_block + ssp + total_crop_area_ha ~ variable,
                   value.var = 'value')
ccl_res_dt[, m_cr_grain := (m_cr_grain*kg_ha)/C_bio]
ccl_res_dt[, s_cr_grain := (s_cr_grain/Mg_ha)/C_bio]
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

# MAPS | SSP3-7.0, 2050 
ntill_res_ssp370_2050 = annual_map(ntill_res_dt, 2050, 'ssp370')
ccg_ntill_ssp370_2050 = annual_map(ccg_ntill_dt, 2050, 'ssp370')
ccl_ntill_ssp370_2050 = annual_map(ccl_ntill_dt, 2050, 'ssp370')

ssp370_2050_maps      = ntill_res_ssp370_2050$GHG + ntill_res_ssp370_2050$Yield +
  ccg_ntill_ssp370_2050$GHG + ccg_ntill_ssp370_2050$Yield +
  ccl_ntill_ssp370_2050$GHG + ccl_ntill_ssp370_2050$Yield +
  plot_layout(ncol = 2, nrow = 3, guides = 'collect') &
  theme(legend.position = 'none')
ssp370_2050_maps     = ssp370_2050_maps + plot_annotation(tag_levels = list(c('a', 'd', 'b', 'e', 'c', 'f'))) &
  theme(plot.tag = element_text(size = 7))
ggsave(paste0(fig.path,'/', 'si/ssp370-annual-2050-ghg-yield.pdf'),     ssp370_2050_maps,  units = 'mm', width = 180, height = 185, device='pdf', dpi=300)

res_ssp370_2050     = annual_map(res_dt, 2050,     'ssp370')
ccg_res_ssp370_2050 = annual_map(ccg_res_dt, 2050, 'ssp370')
ccl_res_ssp370_2050 = annual_map(ccl_res_dt, 2050, 'ssp370')

ssp370_2050_maps2      = res_ssp370_2050$GHG + res_ssp370_2050$Yield +
  ccg_res_ssp370_2050$GHG + ccg_res_ssp370_2050$Yield +
  ccl_res_ssp370_2050$GHG + ccl_res_ssp370_2050$Yield +
  plot_layout(ncol = 2, nrow = 3, guides = 'collect') &
  theme(legend.position = 'none')
ssp370_2050_maps2     = ssp370_2050_maps2 + plot_annotation(tag_levels = list(c('a', 'd', 'b', 'e', 'c', 'f'))) &
  theme(plot.tag = element_text(size = 7))
ggsave(paste0(fig.path,'/', 'si/ssp370-annual-2050-ghg-yield2.pdf'),     ssp370_2050_maps2,  units = 'mm', width = 180, height = 185, device='pdf', dpi=300)

# Maps | SSP1-2.6, 2050
res_ssp126_2050     = annual_map(res_dt, 2050,     'ssp126')
ccg_res_ssp126_2050 = annual_map(ccg_res_dt, 2050, 'ssp126')
ccl_res_ssp126_2050 = annual_map(ccl_res_dt, 2050, 'ssp126')

ssp126_2050_maps2      = res_ssp126_2050$GHG + res_ssp126_2050$Yield +
  ccg_res_ssp126_2050$GHG + ccg_res_ssp126_2050$Yield +
  ccl_res_ssp126_2050$GHG + ccl_res_ssp126_2050$Yield +
  plot_layout(ncol = 2, nrow = 3, guides = 'collect') &
  theme(legend.position = 'none')
ssp126_2050_maps2     = ssp126_2050_maps2 + plot_annotation(tag_levels = list(c('a', 'd', 'b', 'e', 'c', 'f'))) &
  theme(plot.tag = element_text(size = 7))
ggsave(paste0(fig.path,'/', 'si/ssp126-annual-2050-ghg-yield2.pdf'),     ssp126_2050_maps2,  units = 'mm', width = 180, height = 185, device='pdf', dpi=300)
#-------------------------------------------------------------------------------
# FIGURE S5-S16, S18-S21 | GLOBAL HECTARE RESPONSES
#-------------------------------------------------------------------------------
gl_hectare_gcm  = fread(paste(data.path, 'gcm-decadal-global-crop-soil-responses.csv', sep = '/'))
rg_hectare_gcm  = fread(paste(data.path, 'gcm-decadal-IPCC-region-crop-soil-responses.csv', sep = '/'))

# Grain 
# SSP1-2.6 | 2030
grain_gl_rg_126_2030 = grain_gl_rg_ha(gl_hectare_gcm, rg_hectare_gcm, 'ssp126', 2030)
ggsave(grain_gl_rg_126_2030, file = paste(fig.path, '/si/gcm-grain-diff-global-regional-ssp126-2030.pdf', sep = '/'), units = 'mm', width = 180, height = 185, device='pdf', dpi=300)

# SSP1-2.6 | 2050
grain_gl_rg_126_2050 = grain_gl_rg_ha(gl_hectare_gcm, rg_hectare_gcm, 'ssp126', 2050)
ggsave(grain_gl_rg_126_2050, file = paste(fig.path, '/si/gcm-grain-diff-global-regional-ssp126-2050.pdf', sep = '/'), units = 'mm', width = 180, height = 185, device='pdf', dpi=300)

# SSP3-7.0 | 2030
grain_gl_rg_370_2030 = grain_gl_rg_ha(gl_hectare_gcm, rg_hectare_gcm, 'ssp370', 2030)
ggsave(grain_gl_rg_370_2030, file = paste(fig.path, '/si/gcm-grain-diff-global-regional-ssp370-2030.pdf', sep = '/'), units = 'mm', width = 180, height = 185, device='pdf', dpi=300)

# SSP3-7.0 | 2050
grain_gl_rg_370_2050 = grain_gl_rg_ha(gl_hectare_gcm, rg_hectare_gcm, 'ssp370', 2050)
ggsave(grain_gl_rg_370_2050, file = paste(fig.path, '/si/gcm-grain-diff-global-regional-ssp370-2050.pdf', sep = '/'), units = 'mm', width = 180, height = 185, device='pdf', dpi=300)

# SOC
# SSP1-2.6 | 2050 (regions-only)
soc_rg_126_2050 = soc_rg_ha(rg_hectare_gcm, 'ssp126', 2050)
ggsave(soc_rg_126_2050, file = paste(fig.path, '/si/gcm-soc-seq-rate-IPCC-region-ssp126-2050.pdf', sep = '/'), units = 'mm', width = 180, height = 185, device='pdf', dpi=300)

# SSP1-2.6 | 2030
soc_gl_rg_126_2030 = soc_gl_rg_ha(gl_hectare_gcm, rg_hectare_gcm, 'ssp126', 2030)
ggsave(soc_gl_rg_126_2030, file = paste(fig.path, '/si/gcm-soc-seq-rate-global-regional-ssp126-2030.pdf', sep = '/'), units = 'mm', width = 180, height = 185, device='pdf', dpi=300)

# SSP3-7.0 | 2050
soc_gl_rg_370_2050 = soc_gl_rg_ha(gl_hectare_gcm, rg_hectare_gcm, 'ssp370', 2050)
ggsave(soc_gl_rg_370_2050, file = paste(fig.path, '/si/gcm-soc-seq-rate-global-regional-ssp370-2050.pdf', sep = '/'), units = 'mm', width = 180, height = 185, device='pdf', dpi=300)

# SSP3-7.0 | 2030
soc_gl_rg_370_2030 = soc_gl_rg_ha(gl_hectare_gcm, rg_hectare_gcm, 'ssp370', 2030)
ggsave(soc_gl_rg_370_2030, file = paste(fig.path, '/si/gcm-soc-seq-rate-global-regional-ssp370-2030.pdf', sep = '/'), units = 'mm', width = 180, height = 185, device='pdf', dpi=300)

# N2O
# SSP1-2.6 | 2030
n2o_gl_rg_126_2030 = n2o_gl_rg_ha(gl_hectare_gcm, rg_hectare_gcm, 'ssp126', 2030)
ggsave(n2o_gl_rg_126_2030, file = paste(fig.path, '/si/gcm-n2o-rate-global-regional-ssp126-2030.pdf', sep = '/'), units = 'mm', width = 180, height = 185, device='pdf', dpi=300)

# SSP1-2.6 | 2050
n2o_gl_rg_126_2050 = n2o_gl_rg_ha(gl_hectare_gcm, rg_hectare_gcm, 'ssp126', 2050)
ggsave(n2o_gl_rg_126_2050, file = paste(fig.path, '/si/gcm-n2o-rate-global-regional-ssp126-2050.pdf', sep = '/'), units = 'mm', width = 180, height = 185, device='pdf', dpi=300)

# SSP3-7.0 | 2030
n2o_gl_rg_370_2030 = n2o_gl_rg_ha(gl_hectare_gcm, rg_hectare_gcm, 'ssp370', 2030)
ggsave(n2o_gl_rg_370_2030, file = paste(fig.path, '/si/gcm-n2o-rate-global-regional-ssp370-2030.pdf', sep = '/'), units = 'mm', width = 180, height = 185, device='pdf', dpi=300)

# SSP3-7.0 | 2050
n2o_gl_rg_370_2050 = n2o_gl_rg_ha(gl_hectare_gcm, rg_hectare_gcm, 'ssp370', 2050)
ggsave(n2o_gl_rg_370_2050, file = paste(fig.path, '/si/gcm-n2o-rate-global-regional-ssp370-2050.pdf', sep = '/'), units = 'mm', width = 180, height = 185, device='pdf', dpi=300)

# GHG
# SSP1-2.6 | 2030
ghg_gl_rg_126_2030 = ghg_gl_rg_ha(gl_hectare_gcm, rg_hectare_gcm, 'ssp126', 2030)
ggsave(ghg_gl_rg_126_2030, file = paste(fig.path, '/si/gcm-ghg-seq-rate-global-regional-ssp126-2030.pdf', sep = '/'), units = 'mm', width = 180, height = 185, device='pdf', dpi=300)

# SSP1-2.6 | 2050
ghg_gl_rg_126_2050 = ghg_gl_rg_ha(gl_hectare_gcm, rg_hectare_gcm, 'ssp126', 2050)
ggsave(ghg_gl_rg_126_2050, file = paste(fig.path, '/si/gcm-ghg-seq-rate-global-regional-ssp126-2050.pdf', sep = '/'), units = 'mm', width = 180, height = 185, device='pdf', dpi=300)

# SSP3-7.0 | 2030
ghg_gl_rg_370_2030 = ghg_gl_rg_ha(gl_hectare_gcm, rg_hectare_gcm, 'ssp370', 2030)
ggsave(ghg_gl_rg_370_2030, file = paste(fig.path, '/si/gcm-ghg-seq-rate-global-regional-ssp370-2030.pdf', sep = '/'), units = 'mm', width = 180, height = 185, device='pdf', dpi=300)

# SSP3-7.0 | 2050
ghg_gl_rg_370_2050 = ghg_gl_rg_ha(gl_hectare_gcm, rg_hectare_gcm, 'ssp370', 2050)
ggsave(ghg_gl_rg_370_2050, file = paste(fig.path, '/si/gcm-ghg-seq-rate-global-regional-ssp370-2050.pdf', sep = '/'), units = 'mm', width = 180, height = 185, device='pdf', dpi=300)

# COMPARISONS & PERCENT CHANGE ESTIMATES FOR TEXT
# GHG 
rg_GHG = rg_hectare_gcm[, lapply(.SD, mean), .SDcols = 's_GHG', by = .(ssp, scenario, y_block, IPCC_NAME)]
rg_GHG = rg_GHG[y_block == 2050, ]
rg_GHG = rg_GHG[, s_GHG := s_GHG/35L] # annual
rg_GHG[ssp %in% 'ssp126']

gl_GHG = gl_hectare_gcm[, lapply(.SD, mean), .SDcols = 's_GHG', by = .(ssp, scenario, y_block)]
gl_GHG = gl_GHG[y_block == 2050, ]
gl_GHG = gl_GHG[, s_GHG := s_GHG/35L] # annual
gl_GHG[ssp %in% 'ssp126']
# N2O
rg_N2O = rg_hectare_gcm[, lapply(.SD, mean), .SDcols = 's_N2O', by = .(ssp, scenario, y_block, IPCC_NAME)]
rg_N2O = rg_N2O[y_block == 2050, ]
rg_N2O = rg_N2O[, s_N2O := s_N2O/35L] # annual
rg_N2O[ssp %in% 'ssp126']
# SOC
rg_SOC = rg_hectare_gcm[, lapply(.SD, mean), .SDcols = 's_SOC', by = .(ssp, scenario, y_block, IPCC_NAME)]
rg_SOC = rg_SOC[y_block == 2050, ]
rg_SOC = rg_SOC[, s_SOC := s_SOC/35L] # annual
rg_SOC[ssp %in% 'ssp126']

gl_SOC = gl_hectare_gcm[, lapply(.SD, mean), .SDcols = 's_SOC', by = .(ssp, scenario, y_block, IPCC_NAME)]
gl_SOC = gl_SOC[y_block == 2050, ]
gl_SOC = gl_SOC[, s_SOC := s_SOC/35L] # annual
gl_SOC[ssp %in% 'ssp126']

# COMPARISONS
# Grain 
gl_grain = gl_hectare_gcm[, lapply(.SD, mean), .SDcols = 's_cr_grain', by = .(ssp, scenario, y_block)]
gl_grain = gl_grain[y_block == 2050, ]
gl_grain = gl_grain[, s_cr_grain := s_cr_grain/35L] # annual
gl_grain[ssp %in% 'ssp126']

rg_grain = rg_hectare_gcm[, lapply(.SD, mean), .SDcols = 's_cr_grain', by = .(ssp, scenario, y_block, IPCC_NAME)]
rg_grain = rg_grain[y_block == 2050, ]
rg_grain = rg_grain[, s_cr_grain := s_cr_grain/35L] # annual
rg_grain[ssp %in% 'ssp126']

# load reference case for percent change estimates
# values below are absolute emissions and yield
gl_hectare_conv  = fread(paste(data.path, 'ensemble-decadal-global-crop-soil-responses-reference-case.csv', sep = '/'))
rg_hectare_conv  = fread(paste(data.path, 'ensemble-decadal-IPCC-region-crop-soil-responses-reference-case.csv', sep = '/'))
# global
gl_hectare_conv = gl_hectare_conv[y_block == 2050]
gl_hectare_conv = gl_hectare_conv[, s_grain_m := s_grain_m/35L] # annual
gl_hectare_conv[ssp %in% 'ssp126']
print(gl_hectare_conv[, c('IPCC_NAME', 'ssp', 'scenario', 'y_block', 's_grain_m')])
# regional
rg_hectare_conv = rg_hectare_conv[y_block == 2050]
rg_hectare_conv = rg_hectare_conv[, s_grain_m := s_grain_m/35L] # annual
rg_hectare_conv[ssp %in% 'ssp126']

print(rg_hectare_conv[, c('IPCC_NAME', 'ssp', 'scenario', 'y_block', 's_grain_m')])
#-------------------------------------------------------------------------------
# FIGURE S17 | BEST MANAGEMENT PRACTICE
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
ghg_370_bmp = bmp_map(ghg_bmp_dt, 'ssp370')

# YM | SSP1-2.6
ym_370_bmp = bmp_map(ym_bmp_dt, 'ssp370')

# YC  | SSP1-2.6
yc_370_bmp = bmp_map(yc_bmp_dt, 'ssp370')

# Proportions by region, YC | SSP1-2.6
yc_370_bmp_bar = bmp_bplot(yc_bmp_dt, 'ssp370')

bmp_plots = ghg_370_bmp$bmp + yc_370_bmp_bar$freq + ym_370_bmp$bmp + yc_370_bmp_bar$ghg + yc_370_bmp$bmp +
  plot_layout(ncol = 2, nrow = 3, guides = 'collect') &
  theme(legend.position = 'none')

bmp_plots = bmp_plots + plot_annotation(tag_levels = list(c('a', 'd', 'b', 'e', 'c'))) &
  theme(plot.tag = element_text(size = 7))

# save 
ggsave(paste0(fig.path,'/', 'si/ssp370-bmp-plot.pdf'),bmp_plots,  units = 'mm', width = 180, height = 185, device='pdf', dpi=300)

#-------------------------------------------------------------------------------
# FIGURE S22-S28 | SUPERVISED CLUSTERING
#-------------------------------------------------------------------------------
# Legume cover crop, residue, no-tillage | SOC
load(paste(data.path, "ccl-syn-soc_SHAP.Rdata", sep = '/'))
ccl_density_SI = SHAP_density_SI(copy(ccl_k_SOC$kmeans_dt), 'ssp126', 
                                 c('SOMC_sum_', 'm_nfix', 'd_s_annet', 'd_m_nfix',
                                   'd_m_sldcmp', 'd_s_N_leach', 'd_s_cr_grain', 'SLSAND'),
                                 c('SOMC_sum_'    = 'Initial SOC Stock\n(g C m\u207b\u00b2)',
                                   'm_nfix'       = 'N Fixation\n(g N m\u207b\u00b2 yr\u207b\u00b9)',
                                   'd_s_annet'    = 'ET Difference\n(cm H\u2082O yr\u207b\u00b9)',
                                   'd_m_nfix'     = 'N Fixation Difference\n(g N m\u207b\u00b2 yr\u207b\u00b9)',
                                   'd_m_sldcmp'   = 'Soil Decomposition\nMultiplier Difference',
                                   'd_s_N_leach'  = 'Leached Nitrogen Difference\n(g N m\u207b\u00b2 yr\u207b\u00b9)',
                                   'd_s_cr_grain' = 'Grain Yield Difference\n(Mg ha\u207b\u00b9 yr\u207b\u00b9)',
                                   'SLSAND'       = 'Sand Fraction'))
ggsave(ccl_density_SI, file = paste(fig.path, 'si/ccl-syn-SOC-SI-ssp126-2050.tiff', sep = '/'), units = 'mm', width = 180, height = 185, device='tiff', dpi=300)

# Grass cover crop, residue, no-tillage | SOC
load(paste(data.path, "ccg-syn-soc_SHAP.Rdata", sep = '/'))
ccg_density_SI = SHAP_density_SI(copy(ccg_k_SOC$kmeans_dt), 'ssp126', 
                                 c('SOMC_sum_', 's_gr_nit', 'd_s_cr_residC', 'd_s_gr_nit',
                                   'SLSAND', 'res.rtrn.amt', 'd_m_sldcmp', 'd_s_cr_grain'),
                                 c('SOMC_sum_'    = 'Initial SOC Stock\n(g C m\u207b\u00b2)',
                                   's_gr_nit'       = 'Gross Nitrification\n(g N m\u207b\u00b2 yr\u207b\u00b9)',
                                   'd_s_cr_residC' = 'Crop Residue Difference\n(Mg DM ha\u207b\u00b9 yr\u207b\u00b9)',
                                   'd_s_gr_nit'    = 'Gross Nitrification Difference\n(g N m\u207b\u00b2 yr\u207b\u00b9)',
                                   'SLSAND'        = 'Sand Fraction',
                                   'res.rtrn.amt'  = 'Initial Residue\nReturn Fraction',
                                   'd_m_sldcmp'    = 'Soil Decomposition\nMultiplier Difference',
                                   'd_s_cr_grain'  = 'Grain Yield Difference\n(Mg ha\u207b\u00b9 yr\u207b\u00b9)'))
ggsave(ccg_density_SI, file = paste(fig.path, 'si/ccg-syn-SOC-SI-ssp126-2050.tiff', sep = '/'), units = 'mm', width = 180, height = 185, device='tiff', dpi=300)

# Grass cover crop, residue, no-tillage | Yield
load(paste(data.path, "ccg-ntill-yield_SHAP.Rdata", sep = '/'))

ccg_ntill_y_density_SI = SHAP_y_density_SI(copy(ccg_ntill_k_YIELD$kmeans_dt), 'ssp126', 
                                 c('s_gr_nit', 'hist_bio12', 'SLSAND', 'res.rtrn.amt',
                                   's_annet', 'SLBLKD', 'SLCLAY', 'hist_bio1'),
                                 c('s_gr_nit'    = 'Gross Nitrification\n(g N m\u207b\u00b2 yr\u207b\u00b9)', 
                                   'hist_bio12'  = 'Mean Annual Precipitation\n(mm yr\u207b\u00b9)',
                                   'SLSAND'      = 'Sand Fraction',
                                   'res.rtrn.amt'= 'Initial Residue\nReturn Fraction',
                                   's_annet'     = 'ET Difference\n(cm H\u2082O yr\u207b\u00b9)',
                                   'SLBLKD'      = 'Soil Bulk Density\n(g m\u207b\u00b3',
                                   'SLCLAY'      = 'Clay Fraction',
                                   'hist_bio1'   = 'Mean Annual Temperature\n(\u2103)'))
ggsave(ccg_ntill_y_density_SI, file = paste(fig.path, 'si/ccg-ntill-yield-SI-ssp126-2050.tiff', sep = '/'), units = 'mm', width = 180, height = 185, device='tiff', dpi=300)

# Legume cover crop, residue, no-tillage | Yield
load(paste(data.path, "ccl-ntill-yield_SHAP.Rdata", sep = '/'))
ccl_ntill_y_density_SI = SHAP_y_density_SI(copy(ccl_ntill_k_YIELD$kmeans_dt), 'ssp126', 
                                           c('N.amt', 's_gr_nit', 's_N_leach', 's_SOC',
                                             'SLSAND', 'SLBLKD', 'SLCLAY', 's_cc_biomass'),
                                           c('N.amt'       = 'Nitrogen Inputs\n(g N m\u207b\u00b2 yr\u207b\u00b9)', 
                                             's_gr_nit'    = 'Gross Nitrification Difference\n(g N m\u207b\u00b2 yr\u207b\u00b9)',
                                             's_N_leach'   = "Leached Nitrogen Difference\n(g N m\u207b\u00b2 yr\u207b\u00b9)",
                                             's_SOC'       = "Soil Organic Carbon Difference\n(Mg CO\u2082-eq ha\u207b\u00b9 yr\u207b\u00b9)",
                                             'SLSAND'      = 'Sand Fraction',
                                             'SLBLKD'      = 'Soil Bulk Density\n(g m\u207b\u00b3',
                                             'SLCLAY'      = 'Clay Fraction',
                                             's_cc_biomass'= "Cover Crop Biomass\n(Mg DM ha\u207b\u00b9 yr\u207b\u00b9)"))
ggsave(ccl_ntill_y_density_SI, file = paste(fig.path, 'si/ccl-ntill-yield-SI-ssp126-2050.tiff', sep = '/'), units = 'mm', width = 180, height = 185, device='tiff', dpi=300)

# No-tillage, residue | Yield
load(paste(data.path, "ntill-res-yield_SHAP.Rdata", sep = '/'))
ntill_res_y_map     = cluster_map(ntill_res_k_YIELD$kmeans_dt, 'ssp126')
ntill_res_y_density   = SHAP_y_density2(copy(ntill_res_k_YIELD$kmeans_dt), 'ssp126', 
                           c('s_cr_grain', 's_N_leach', 's_gr_nit', 's_annet'),
                           c('s_cr_grain' = "Grain Yield Difference\n(Mg ha\u207b\u00b9 yr\u207b\u00b9)",
                             's_N_leach'  = "Leached Nitrogen Difference\n(g N m\u207b\u00b2 yr\u207b\u00b9)",
                             's_gr_nit'   = "Gross Nitrification Difference\n(g N m\u207b\u00b2 yr\u207b\u00b9)",
                             's_annet'    = "ET Difference\n(cm H\u2082O yr\u207b\u00b9)"))
yield_ntill_res_features = ntill_res_y_density + ntill_res_y_map  +
  plot_layout(ncol = 1, nrow = 2, guides = 'collect') &
  theme(legend.position = 'none')
yield_ntill_res_features = yield_ntill_res_features + plot_annotation(tag_levels = 'a') &
  theme(plot.tag = element_text(size = 7))
# add legend manually
ggsave(paste0(fig.path,'/si/yield-ntill-res-features-main.tiff'), yield_ntill_res_features, units = 'mm', width = 90, height = 160, device='tiff', dpi=300)
# other features
ntill_res_y_density_SI = SHAP_y_density_SI2(copy(ccl_ntill_k_YIELD$kmeans_dt), 'ssp126', 
                                            c('SLSAND', 'hist_bio1', 's_SOC', 'frac_NO3',
                                              'res.rtrn.amt', 'hist_bio12', 'frac_Urea', 'N.amt'),
                                            c('SLSAND'      = 'Sand Fraction', 
                                              'hist_bio1'   = 'Mean Annual Temperature\n(\u2103)',
                                              's_SOC'       = "Soil Organic Carbon Difference\n(Mg CO\u2082-eq ha\u207b\u00b9 yr\u207b\u00b9)",
                                              'frac_NO3'    = "Nitrate Fertilizer Fraction",
                                              'res.rtrn.amt'= 'Initial Residue\nReturn Fraction',
                                              'hist_bio12'  = 'Mean Annual Precipitation\n(mm yr\u207b\u00b9)',
                                              'frac_Urea'   = 'Urea Fertilizer Fraction',
                                              'N.amt'       = "Nitrogen Inputs\n(g N m\u207b\u00b2 yr\u207b\u00b9)"))
ggsave(ntill_res_y_density_SI, file = paste(fig.path, 'si/ntill-res-yield-SI-ssp126-2050.tiff', sep = '/'), units = 'mm', width = 180, height = 185, device='tiff', dpi=300)


#-------------------------------------------------------------------------------
# OLD#
#-------------------------------------------------------------------------------
# CLIMATE SENSITIVITY
#-------------------------------------------------------------------------------
# global sensitivity
gl_sens_dt = fread(paste(data.path, 'ensemble-decadal-global-crop-soil-responses.csv', sep = '/'))
gl_cols    = c('ssp','scenario','y_block','s_GHG_m', 's_grain_m')
gl_sens_dt = gl_sens_dt[, ..gl_cols]

# separate historical
gl_sens_hist_dt = gl_sens_dt[ssp %in% 'historical',]
gl_sens_dt      = gl_sens_dt[!ssp %in% 'historical',]
# update historical col names
setnames(gl_sens_hist_dt, 's_GHG_m','hist_s_GHG_m')
setnames(gl_sens_hist_dt, 's_grain_m','hist_s_grain_m')
# join to ssps
gl_sens_dt      = gl_sens_dt[gl_sens_hist_dt, on = .(scenario = scenario,
                                                     y_block  = y_block)]
gl_sens_dt[, i.ssp := NULL]
# compute % change | 2030, 2050, 2100
gl_sens_dt      = gl_sens_dt[y_block %in% c(2030, 2050, 2100)]
gl_sens_dt[, ghg_diff      := hist_s_GHG_m - s_GHG_m]
gl_sens_dt[, ghg_perc_diff := (ghg_diff/hist_s_GHG_m)*100]
gl_sens_dt[, ghg_perc_diff := round(ghg_perc_diff, digits = 0)]

gl_sens_dt[, gr_diff      := abs(hist_s_grain_m) - abs(s_grain_m)] # for interpretation
gl_sens_dt[, gr_perc_diff := (gr_diff/hist_s_grain_m)*100]
gl_sens_dt[, gr_perc_diff := round(gr_perc_diff, digits = 0)]

setcolorder(gl_sens_dt, c('ssp', 'scenario', 'y_block', 'gr_perc_diff', 'ghg_perc_diff'))

# regional sensitivity
rg_sens_dt = fread(paste(data.path, 'ensemble-decadal-IPCC-region-crop-soil-responses.csv', sep = '/'))
rg_cols    = c('ssp','scenario','y_block','IPCC_NAME','s_GHG_m', 's_grain_m')
rg_sens_dt = rg_sens_dt[, ..rg_cols]

# separate historical
rg_sens_hist_dt = rg_sens_dt[ssp %in% 'historical',]
rg_sens_dt      = rg_sens_dt[!ssp %in% 'historical',]
# update historical col names
setnames(rg_sens_hist_dt, 's_GHG_m','hist_s_GHG_m')
setnames(rg_sens_hist_dt, 's_grain_m','hist_s_grain_m')
# join to ssps
rg_sens_dt      = rg_sens_dt[rg_sens_hist_dt, on = .(scenario = scenario,
                                                     y_block  = y_block,
                                                     IPCC_NAME = IPCC_NAME)]
rg_sens_dt[, i.ssp := NULL]
# compute % change | 2030, 2050, 2100
rg_sens_dt      = rg_sens_dt[y_block %in% c(2030, 2050, 2100)]
rg_sens_dt[, ghg_diff      := hist_s_GHG_m - s_GHG_m]
rg_sens_dt[, ghg_perc_diff := (ghg_diff/hist_s_GHG_m)*100]
rg_sens_dt[, ghg_perc_diff := round(ghg_perc_diff, digits = 0)]

rg_sens_dt[, gr_diff      := abs(hist_s_grain_m) - abs(s_grain_m)] # for interpretation
rg_sens_dt[, gr_perc_diff := (gr_diff/hist_s_grain_m)*100]
rg_sens_dt[, gr_perc_diff := round(gr_perc_diff, digits = 0)]

setcolorder(rg_sens_dt, c('ssp', 'scenario', 'y_block', 'IPCC_NAME','gr_perc_diff', 'ghg_perc_diff'))
