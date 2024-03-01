# filename:    figures_main.R
# created:     06 March 2023
# updated:     26 February 2024
# author:      S.C. McClelland
# description: This file creates main text figures for manuscript.
#-------------------------------------------------------------------------------
library(data.table)
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
# FIGURE 1 | BIOPHYSICAL & YIELD CONSTRAINED POTENTIAL
#-------------------------------------------------------------------------------
# GLOBAL
biophys_2050   = fread(paste(data.path, 'global-biophysical-GHG-mitigation-potential-2050.csv', sep = '/'))
yield_cst_2050 = fread(paste(data.path, 'global-yield-constrained-GHG-mitigation-potential-2050.csv', sep = '/'))

# REGIONAL
ipcc_biophys_2050   = fread(paste(data.path, 'regional-biophysical-GHG-mitigation-potential-2050.csv', sep = '/'))
ipcc_yield_cst_2050 = fread(paste(data.path, 'regional-yield-constrained-GHG-mitigation-potential-2050.csv', sep = '/'))

# background map
map_background = IPCC_map(lu.path, shp, raster)
ggsave(paste0(fig.path,'/IPCC-map-background.tiff'), map_background$IPCC, units = 'in', width = 9, height = 5, device='tiff', dpi=300)
# global
global_126_ptl = gl_bar_potential(biophys_2050, yield_cst_2050, 'ssp126')
ggsave(paste0(fig.path,'/ssp126-2050global-ghg-mitigation-potential.tiff'), global_126_ptl$GHG, bg = 'transparent',units = 'in', width = 8, height = 6, device='tiff', dpi=300)
# regional
ADP_126_ptl    = rg_bar_potential(ipcc_biophys_2050, ipcc_yield_cst_2050, 'ssp126', 'ADP')
ggsave(paste0(fig.path,'/ssp126-2050-ADP-ghg-mitigation-potential.tiff'), ADP_126_ptl$GHG, bg = 'transparent',units = 'in', width = 8, height = 6, device='tiff', dpi=300)
AME_126_ptl    = rg_bar_potential(ipcc_biophys_2050, ipcc_yield_cst_2050, 'ssp126', 'AME')
ggsave(paste0(fig.path,'/ssp126-2050-AME-ghg-mitigation-potential.tiff'), AME_126_ptl$GHG, bg = 'transparent',units = 'in', width = 8, height = 6, device='tiff', dpi=300)
DEV_126_ptl    = rg_bar_potential(ipcc_biophys_2050, ipcc_yield_cst_2050, 'ssp126', 'DEV')
ggsave(paste0(fig.path,'/ssp126-2050-DEV-ghg-mitigation-potential.tiff'), DEV_126_ptl$GHG, bg = 'transparent',units = 'in', width = 8, height = 6, device='tiff', dpi=300)
EEWCA_126_ptl  = rg_bar_potential(ipcc_biophys_2050, ipcc_yield_cst_2050, 'ssp126', 'EEWCA')
ggsave(paste0(fig.path,'/ssp126-2050-EEWCA-ghg-mitigation-potential.tiff'), EEWCA_126_ptl$GHG, bg = 'transparent',units = 'in', width = 8, height = 6, device='tiff', dpi=300)
LAC_126_ptl    = rg_bar_potential(ipcc_biophys_2050, ipcc_yield_cst_2050, 'ssp126', 'LAC')
ggsave(paste0(fig.path,'/ssp126-2050-LAC-ghg-mitigation-potential.tiff'), LAC_126_ptl$GHG, bg = 'transparent',units = 'in', width = 8, height = 6, device='tiff', dpi=300)
#-------------------------------------------------------------------------------
# FIGURE 2 | MAPS
#-------------------------------------------------------------------------------
Mg_ha = 100L
kg_ha = 10L
C_bio = 0.45
  # LOAD DATA
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

  # MAPS | SSP1-2.6, 2050
ntill_res_ssp126_2050 = annual_map(ntill_res_dt, 2050, 'ssp126')
ccg_res_ssp126_2050   = annual_map(ccg_res_dt, 2050, 'ssp126')
ccl_res_ssp126_2050   = annual_map(ccl_res_dt, 2050, 'ssp126')
ccg_ntill_ssp126_2050 = annual_map(ccg_ntill_dt, 2050, 'ssp126')
ccl_ntill_ssp126_2050 = annual_map(ccl_ntill_dt, 2050, 'ssp126')

ssp126_2050_maps      = ntill_res_ssp126_2050$GHG + ntill_res_ssp126_2050$Yield +
  ccg_res_ssp126_2050$GHG   + ccg_res_ssp126_2050$Yield +
  ccl_res_ssp126_2050$GHG   + ccl_res_ssp126_2050$Yield +
  ccg_ntill_ssp126_2050$GHG + ccg_ntill_ssp126_2050$Yield +
  ccl_ntill_ssp126_2050$GHG + ccl_ntill_ssp126_2050$Yield +
  plot_layout(ncol = 2, nrow = 5, guides = 'collect') &
  theme(legend.position = 'none')
ssp126_2050_maps     = ssp126_2050_maps + plot_annotation(tag_levels = 'a')
ggsave(paste0(fig.path,'/', 'ssp126-annual-2050-ghg-yield.tiff'),     ssp126_2050_maps,  units = 'in', width = 10, height = 14, device='tiff', dpi=300)
ggsave(paste0(fig.path,'/', 'ssp126-annual-2050-ghg-legend.tiff'),    ntill_res_ssp126_2050$legend1,  units = 'in', width = 9, height = 5, device='tiff', dpi=300)
ggsave(paste0(fig.path,'/', 'ssp126-annual-2050-yield-legend.tiff'),  ntill_res_ssp126_2050$legend2,  units = 'in', width = 9, height = 5, device='tiff', dpi=300)

# MAPS | SSP1-2.6, 2030
ntill_res_ssp126_2030 = annual_map(ntill_res_dt, 2030, 'ssp126')
ccg_res_ssp126_2030   = annual_map(ccg_res_dt, 2030, 'ssp126')
ccl_res_ssp126_2030   = annual_map(ccl_res_dt, 2030, 'ssp126')
ccg_ntill_ssp126_2030 = annual_map(ccg_ntill_dt, 2030, 'ssp126')
ccl_ntill_ssp126_2030 = annual_map(ccl_ntill_dt, 2030, 'ssp126')

ssp126_2030_maps      = ntill_res_ssp126_2030$GHG + ntill_res_ssp126_2030$Yield +
  ccg_res_ssp126_2030$GHG   + ccg_res_ssp126_2030$Yield +
  ccl_res_ssp126_2030$GHG   + ccl_res_ssp126_2030$Yield +
  ccg_ntill_ssp126_2030$GHG + ccg_ntill_ssp126_2030$Yield +
  ccl_ntill_ssp126_2030$GHG + ccl_ntill_ssp126_2030$Yield +
  plot_layout(ncol = 2, nrow = 5, guides = 'collect') &
  theme(legend.position = 'none')
ssp126_2030_maps     = ssp126_2030_maps + plot_annotation(tag_levels = 'a')
ggsave(paste0(fig.path,'/', 'ssp126-annual-2030-ghg-yield.tiff'),     ssp126_2030_maps,  units = 'in', width = 10, height = 14, device='tiff', dpi=300)
ggsave(paste0(fig.path,'/', 'ssp126-annual-2030-ghg-legend.tiff'),    ntill_res_ssp126_2030$legend1,  units = 'in', width = 9, height = 5, device='tiff', dpi=300)
ggsave(paste0(fig.path,'/', 'ssp126-annual-2030-yield-legend.tiff'),  ntill_res_ssp126_2030$legend2,  units = 'in', width = 9, height = 5, device='tiff', dpi=300)
#-------------------------------------------------------------------------------
# FIGURE 3 | REGIONAL HECTARE RESPONSE
#-------------------------------------------------------------------------------
rg_hectare_gcm  = fread(paste(data.path, 'gcm-decadal-IPCC-region-crop-soil-responses.csv', sep = '/'))

# SSP1-2.6 | 2050
soc_rg_126_2050 = soc_rg_ha(rg_hectare_gcm, 'ssp126', 2050)
ggsave(soc_rg_126_2050, file = paste(fig.path, 'gcm-soc-seq-rate-IPCC-region-ssp126-2050.tiff', sep = '/'), units = 'in', width = 9, height = 10, device='tiff', dpi=300)
#-------------------------------------------------------------------------------
# FIGURE X | BEST MANAGEMENT PRACTICE
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

# GHG | SSP1-2.6
ghg_126_bmp = bmp_map(ghg_bmp_dt, 'ssp126')
ggsave(paste0(fig.path,'/ghg-bmp-ssp126-2050.tiff'), ghg_126_bmp$bmp, units = 'in', width = 9, height = 5, device='tiff', dpi=300)

# YC  | SSP1-2.6
yc_126_bmp = bmp_map(yc_bmp_dt, 'ssp126')
ggsave(paste0(fig.path,'/yield-constrained-bmp-ssp126-2050.tiff'), yc_126_bmp$bmp, units = 'in', width = 9, height = 5, device='tiff', dpi=300)



#-------------------------------------------------------------------------------
# OLD - SAVE FOR NOW
#-------------------------------------------------------------------------------

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