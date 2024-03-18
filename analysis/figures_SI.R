# filename:    figures_main.R
# created:     06 March 2023
# updated:     12 March 2024
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
# FIGURE S1 | BIOPHYSICAL & YIELD CONSTRAINED POTENTIAL
#-------------------------------------------------------------------------------
# GLOBAL
biophys_2050     = fread(paste(data.path, 'global-biophysical-GHG-mitigation-potential-2050.csv', sep = '/'))
yield_cst_2050   = fread(paste(data.path, 'global-yield-constrained-GHG-mitigation-potential-2050.csv', sep = '/'))
yield_cst_2050_l = fread(paste(data.path, 'global-yield-constrained-low-price-GHG-mitigation-potential-2050.csv', sep = '/'))
yield_cst_2050_h = fread(paste(data.path, 'global-yield-constrained-high-price-GHG-mitigation-potential-2050.csv', sep = '/'))
# REGIONAL
ipcc_biophys_2050     = fread(paste(data.path, 'regional-biophysical-GHG-mitigation-potential-2050.csv', sep = '/'))
ipcc_yield_cst_2050   = fread(paste(data.path, 'regional-yield-constrained-GHG-mitigation-potential-2050.csv', sep = '/'))
ipcc_yield_cst_2050_l = fread(paste(data.path, 'regional-yield-constrained-low-price-GHG-mitigation-potential-2050.csv', sep = '/'))
ipcc_yield_cst_2050_h = fread(paste(data.path, 'regional-yield-constrained-high-price-GHG-mitigation-potential-2050.csv', sep = '/'))

# background map
map_background = IPCC_map(lu.path, shp, raster)
ggsave(paste0(fig.path,'/IPCC-map-background.tiff'), map_background$IPCC, units = 'in', width = 9, height = 5, device='tiff', dpi=300)
# global
global_370_ptl = gl_bar_potential(biophys_2050, yield_cst_2050, yield_cst_2050_l, yield_cst_2050_h,'ssp370')
ggsave(paste0(fig.path,'/si/ssp370-2050-global-ghg-mitigation-potential.tiff'), global_370_ptl$GHG, bg = 'transparent',units = 'in', width = 8, height = 6, device='tiff', dpi=300)
# regional
ADP_370_ptl    = rg_bar_potential(ipcc_biophys_2050, ipcc_yield_cst_2050, ipcc_yield_cst_2050_l, ipcc_yield_cst_2050_h, 'ssp370', 'ADP')
ggsave(paste0(fig.path,'/si/ssp370-2050-ADP-ghg-mitigation-potential.tiff'), ADP_370_ptl$GHG, bg = 'transparent',units = 'in', width = 8, height = 6, device='tiff', dpi=300)
AME_370_ptl    = rg_bar_potential(ipcc_biophys_2050, ipcc_yield_cst_2050, ipcc_yield_cst_2050_l, ipcc_yield_cst_2050_h, 'ssp370', 'AME')
ggsave(paste0(fig.path,'/si/ssp370-2050-AME-ghg-mitigation-potential.tiff'), AME_370_ptl$GHG, bg = 'transparent',units = 'in', width = 8, height = 6, device='tiff', dpi=300)
DEV_370_ptl    = rg_bar_potential(ipcc_biophys_2050, ipcc_yield_cst_2050, ipcc_yield_cst_2050_l, ipcc_yield_cst_2050_h, 'ssp370', 'DEV')
ggsave(paste0(fig.path,'/si/ssp370-2050-DEV-ghg-mitigation-potential.tiff'), DEV_370_ptl$GHG, bg = 'transparent',units = 'in', width = 8, height = 6, device='tiff', dpi=300)
EEWCA_370_ptl  = rg_bar_potential(ipcc_biophys_2050, ipcc_yield_cst_2050, ipcc_yield_cst_2050_l, ipcc_yield_cst_2050_h, 'ssp370', 'EEWCA')
ggsave(paste0(fig.path,'/si/ssp370-2050-EEWCA-ghg-mitigation-potential.tiff'), EEWCA_370_ptl$GHG, bg = 'transparent',units = 'in', width = 8, height = 6, device='tiff', dpi=300)
LAC_370_ptl    = rg_bar_potential(ipcc_biophys_2050, ipcc_yield_cst_2050, ipcc_yield_cst_2050_l, ipcc_yield_cst_2050_h, 'ssp370', 'LAC')
ggsave(paste0(fig.path,'/si/ssp370-2050-LAC-ghg-mitigation-potential.tiff'), LAC_370_ptl$GHG, bg = 'transparent',units = 'in', width = 8, height = 6, device='tiff', dpi=300)
#-------------------------------------------------------------------------------
# FIGURE SX-X | GLOBAL HECTARE RESPONSE
#-------------------------------------------------------------------------------
gl_hectare_gcm  = fread(paste(data.path, 'gcm-decadal-global-crop-soil-responses.csv', sep = '/'))
rg_hectare_gcm  = fread(paste(data.path, 'gcm-decadal-IPCC-region-crop-soil-responses.csv', sep = '/'))

  # SOC
# SSP3-7.0 | 2050
soc_gl_rg_370_2050 = soc_gl_rg_ha(gl_hectare_gcm, rg_hectare_gcm, 'ssp370', 2050)
ggsave(soc_gl_rg_370_2050, file = paste(fig.path, '/si/gcm-soc-seq-rate-global-regional-ssp370-2050.tiff', sep = '/'), units = 'in', width = 9, height = 10, device='tiff', dpi=300)

# SSP1-2.6 | 2030
soc_gl_rg_126_2030 = soc_gl_rg_ha(gl_hectare_gcm, rg_hectare_gcm, 'ssp126', 2030)
ggsave(soc_gl_rg_126_2030, file = paste(fig.path, '/si/gcm-soc-seq-rate-global-regional-ssp126-2030.tiff', sep = '/'), units = 'in', width = 9, height = 10, device='tiff', dpi=300)

# SSP3-7.0 | 2030
soc_gl_rg_370_2030 = soc_gl_rg_ha(gl_hectare_gcm, rg_hectare_gcm, 'ssp370', 2030)
ggsave(soc_gl_rg_370_2030, file = paste(fig.path, '/si/gcm-soc-seq-rate-global-regional-ssp370-2030.tiff', sep = '/'), units = 'in', width = 9, height = 10, device='tiff', dpi=300)

  # GHG
# SSP1-2.6 | 2030
ghg_gl_rg_126_2030 = ghg_gl_rg_ha(gl_hectare_gcm, rg_hectare_gcm, 'ssp126', 2030)
ggsave(ghg_gl_rg_126_2030, file = paste(fig.path, '/si/gcm-ghg-seq-rate-global-regional-ssp126-2030.tiff', sep = '/'), units = 'in', width = 9, height = 10, device='tiff', dpi=300)

# SSP1-2.6 | 2050
ghg_gl_rg_126_2050 = ghg_gl_rg_ha(gl_hectare_gcm, rg_hectare_gcm, 'ssp126', 2050)
ggsave(ghg_gl_rg_126_2050, file = paste(fig.path, '/si/gcm-ghg-seq-rate-global-regional-ssp126-2050.tiff', sep = '/'), units = 'in', width = 9, height = 10, device='tiff', dpi=300)

# SSP3-7.0 | 2030
ghg_gl_rg_370_2030 = ghg_gl_rg_ha(gl_hectare_gcm, rg_hectare_gcm, 'ssp370', 2030)
ggsave(ghg_gl_rg_370_2030, file = paste(fig.path, '/si/gcm-ghg-seq-rate-global-regional-ssp370-2030.tiff', sep = '/'), units = 'in', width = 9, height = 10, device='tiff', dpi=300)

# SSP3-7.0 | 2050
ghg_gl_rg_370_2050 = ghg_gl_rg_ha(gl_hectare_gcm, rg_hectare_gcm, 'ssp370', 2050)
ggsave(ghg_gl_rg_370_2050, file = paste(fig.path, '/si/gcm-ghg-seq-rate-global-regional-ssp370-2050.tiff', sep = '/'), units = 'in', width = 9, height = 10, device='tiff', dpi=300)

# Grain 
  # SSP1-2.6 | 2030
grain_gl_rg_126_2030 = grain_gl_rg_ha(gl_hectare_gcm, rg_hectare_gcm, 'ssp126', 2030)
gl_grain = gl_hectare_gcm[, lapply(.SD, mean), .SDcols = 's_cr_grain', by = .(ssp, scenario, y_block)]
gl_grain = gl_grain[y_block == 2030, ]
gl_grain = gl_grain[, s_cr_grain := s_cr_grain/15L]
gl_grain[ssp %in% 'ssp126']

rg_grain = rg_hectare_gcm[, lapply(.SD, mean), .SDcols = 's_cr_grain', by = .(ssp, scenario, y_block, IPCC_NAME)]
rg_grain = rg_grain[y_block == 2030, ]
rg_grain = rg_grain[, s_cr_grain := s_cr_grain/15L]
rg_grain[ssp %in% 'ssp126']

  # SSP1-2.6 | 2050
grain_gl_rg_126_2050 = grain_gl_rg_ha(gl_hectare_gcm, rg_hectare_gcm, 'ssp126', 2050)
gl_grain = gl_hectare_gcm[, lapply(.SD, mean), .SDcols = 's_cr_grain', by = .(ssp, scenario, y_block)]
gl_grain = gl_grain[y_block == 2050, ]
gl_grain = gl_grain[, s_cr_grain := s_cr_grain/35L]
gl_grain[ssp %in% 'ssp126']

rg_grain = rg_hectare_gcm[, lapply(.SD, mean), .SDcols = 's_cr_grain', by = .(ssp, scenario, y_block, IPCC_NAME)]
rg_grain = rg_grain[y_block == 2050, ]
rg_grain = rg_grain[, s_cr_grain := s_cr_grain/35L]
rg_grain[ssp %in% 'ssp126']
#-------------------------------------------------------------------------------
# FIGURE SX-X | MAPS
#-------------------------------------------------------------------------------
Mg_ha = 100L
kg_ha = 10L
C_bio = 0.45
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

  # MAPS | SSP3-7.0, 2030
ntill_res_ssp370_2030 = annual_map(ntill_res_dt, 2030, 'ssp370')
ccg_ntill_ssp370_2030 = annual_map(ccg_ntill_dt, 2030, 'ssp370')
ccl_ntill_ssp370_2030 = annual_map(ccl_ntill_dt, 2030, 'ssp370')

ssp370_2030_maps      = ntill_res_ssp370_2030$GHG + ntill_res_ssp370_2030$Yield +
  ccg_ntill_ssp370_2030$GHG + ccg_ntill_ssp370_2030$Yield +
  ccl_ntill_ssp370_2030$GHG + ccl_ntill_ssp370_2030$Yield +
  plot_layout(ncol = 2, nrow = 3, guides = 'collect') &
  theme(legend.position = 'none')
ssp370_2030_maps     = ssp370_2030_maps + plot_annotation(tag_levels = 'a')
ggsave(paste0(fig.path,'/', 'si/ssp370-annual-2030-ghg-yield.tiff'),     ssp370_2030_maps,  units = 'in', width = 10, height = 14, device='tiff', dpi=300)

res_ssp370_2030     = annual_map(res_dt, 2030,     'ssp370')
ccg_res_ssp370_2030 = annual_map(ccg_res_dt, 2030, 'ssp370')
ccl_res_ssp370_2030 = annual_map(ccl_res_dt, 2030, 'ssp370')

ssp370_2030_maps2      = res_ssp370_2030$GHG + res_ssp370_2030$Yield +
  ccg_res_ssp370_2030$GHG + ccg_res_ssp370_2030$Yield +
  ccl_res_ssp370_2030$GHG + ccl_res_ssp370_2030$Yield +
  plot_layout(ncol = 2, nrow = 3, guides = 'collect') &
  theme(legend.position = 'none')
ssp370_2030_maps2     = ssp370_2030_maps2 + plot_annotation(tag_levels = 'a')
ggsave(paste0(fig.path,'/', 'si/ssp370-annual-2030-ghg-yield2.tiff'),     ssp370_2030_maps2,  units = 'in', width = 10, height = 14, device='tiff', dpi=300)

  # MAPS | SSP3-7.0, 2050 
ntill_res_ssp370_2050 = annual_map(ntill_res_dt, 2050, 'ssp370')
ccg_ntill_ssp370_2050 = annual_map(ccg_ntill_dt, 2050, 'ssp370')
ccl_ntill_ssp370_2050 = annual_map(ccl_ntill_dt, 2050, 'ssp370')

ssp370_2050_maps      = ntill_res_ssp370_2050$GHG + ntill_res_ssp370_2050$Yield +
  ccg_ntill_ssp370_2050$GHG + ccg_ntill_ssp370_2050$Yield +
  ccl_ntill_ssp370_2050$GHG + ccl_ntill_ssp370_2050$Yield +
  plot_layout(ncol = 2, nrow = 3, guides = 'collect') &
  theme(legend.position = 'none')
ssp370_2050_maps     = ssp370_2050_maps + plot_annotation(tag_levels = 'a')
ggsave(paste0(fig.path,'/', 'si/ssp370-annual-2050-ghg-yield.tiff'),     ssp370_2050_maps,  units = 'in', width = 10, height = 14, device='tiff', dpi=300)

res_ssp370_2050     = annual_map(res_dt, 2050,     'ssp370')
ccg_res_ssp370_2050 = annual_map(ccg_res_dt, 2050, 'ssp370')
ccl_res_ssp370_2050 = annual_map(ccl_res_dt, 2050, 'ssp370')

ssp370_2050_maps2      = res_ssp370_2050$GHG + res_ssp370_2050$Yield +
  ccg_res_ssp370_2050$GHG + ccg_res_ssp370_2050$Yield +
  ccl_res_ssp370_2050$GHG + ccl_res_ssp370_2050$Yield +
  plot_layout(ncol = 2, nrow = 3, guides = 'collect') &
  theme(legend.position = 'none')
ssp370_2050_maps2     = ssp370_2050_maps2 + plot_annotation(tag_levels = 'a')
ggsave(paste0(fig.path,'/', 'si/ssp370-annual-2050-ghg-yield2.tiff'),     ssp370_2050_maps2,  units = 'in', width = 10, height = 14, device='tiff', dpi=300)


  # MAPS | SSP1-2.6, 2030
ntill_res_ssp126_2030 = annual_map(ntill_res_dt, 2030, 'ssp126')
ccg_ntill_ssp126_2030 = annual_map(ccg_ntill_dt, 2030, 'ssp126')
ccl_ntill_ssp126_2030 = annual_map(ccl_ntill_dt, 2030, 'ssp126')

ssp126_2030_maps      = ntill_res_ssp126_2030$GHG + ntill_res_ssp126_2030$Yield +
  ccg_ntill_ssp126_2030$GHG + ccg_ntill_ssp126_2030$Yield +
  ccl_ntill_ssp126_2030$GHG + ccl_ntill_ssp126_2030$Yield +
  plot_layout(ncol = 2, nrow = 3, guides = 'collect') &
  theme(legend.position = 'none')
ssp126_2030_maps     = ssp126_2030_maps + plot_annotation(tag_levels = 'a')
ggsave(paste0(fig.path,'/', 'si/ssp126-annual-2030-ghg-yield.tiff'),     ssp126_2030_maps,  units = 'in', width = 10, height = 14, device='tiff', dpi=300)

res_ssp126_2030     = annual_map(res_dt, 2030,     'ssp126')
ccg_res_ssp126_2030 = annual_map(ccg_res_dt, 2030, 'ssp126')
ccl_res_ssp126_2030 = annual_map(ccl_res_dt, 2030, 'ssp126')

ssp126_2030_maps2      = res_ssp126_2030$GHG + res_ssp126_2030$Yield +
  ccg_res_ssp126_2030$GHG + ccg_res_ssp126_2030$Yield +
  ccl_res_ssp126_2030$GHG + ccl_res_ssp126_2030$Yield +
  plot_layout(ncol = 2, nrow = 3, guides = 'collect') &
  theme(legend.position = 'none')
ssp126_2030_maps2     = ssp126_2030_maps2 + plot_annotation(tag_levels = 'a')
ggsave(paste0(fig.path,'/', 'si/ssp126-annual-2030-ghg-yield2.tiff'),     ssp126_2030_maps2,  units = 'in', width = 10, height = 14, device='tiff', dpi=300)

res_ssp126_2050     = annual_map(res_dt, 2050,     'ssp126')
ccg_res_ssp126_2050 = annual_map(ccg_res_dt, 2050, 'ssp126')
ccl_res_ssp126_2050 = annual_map(ccl_res_dt, 2050, 'ssp126')

ssp126_2050_maps2      = res_ssp126_2050$GHG + res_ssp126_2050$Yield +
  ccg_res_ssp126_2050$GHG + ccg_res_ssp126_2050$Yield +
  ccl_res_ssp126_2050$GHG + ccl_res_ssp126_2050$Yield +
  plot_layout(ncol = 2, nrow = 3, guides = 'collect') &
  theme(legend.position = 'none')
ssp126_2050_maps2     = ssp126_2050_maps2 + plot_annotation(tag_levels = 'a')
ggsave(paste0(fig.path,'/', 'si/ssp126-annual-2050-ghg-yield2.tiff'),     ssp126_2050_maps2,  units = 'in', width = 10, height = 14, device='tiff', dpi=300)
#-------------------------------------------------------------------------------
# FIGURE SX-X | CLIMATE SENSITIVITY
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
# FIGURE SX-SX | BEST MANAGEMENT PRACTICE
#-------------------------------------------------------------------------------
IPCC_dt = ipcc_name(lu.path, shp, raster)
ghg_bmp_dt = fread(paste(data.path, 'best-management-practice-GHG-mitigation-2050.csv', sep = '/'))
ghg_bmp_dt = ghg_bmp_dt[IPCC_dt[, .(cell, IPCC_NAME, WB_NAME)], 
           on = .(cell = cell) ]
ghg_bmp_dt = ghg_bmp_dt[!is.na(ssp)]

yc_bmp_dt  = fread(paste(data.path, 'best-management-practice-GHG-mitigation-yield-constrained-2050.csv', sep = '/'))
yc_bmp_dt  = yc_bmp_dt[IPCC_dt[, .(cell, IPCC_NAME, WB_NAME)], 
                        on = .(cell = cell) ]
yc_bmp_dt  = ghg_bmp_dt[!is.na(ssp)]

# GHG | SSP3-7.0
ghg_370_bmp = bmp_map(ghg_bmp_dt, 'ssp370')
ggsave(paste0(fig.path,'/ghg-bmp-ssp370-2050.tiff'), ghg_370_bmp$bmp, units = 'in', width = 9, height = 5, device='tiff', dpi=300)

# YC  | SSP3-7.0
yc_370_bmp = bmp_map(yc_bmp_dt, 'ssp370')
ggsave(paste0(fig.path,'/yield-constrained-bmp-ssp370-2050.tiff'), yc_370_bmp$bmp, units = 'in', width = 9, height = 5, device='tiff', dpi=300)
