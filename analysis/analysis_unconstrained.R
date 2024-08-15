# filename:    analysis_unconstrained.R
# created:     23 February 2024
# updated:     15 August 2024
# author:      S.C. McClelland
# description: This file analyzes GHG mitigation potential from DayCent simulations
#              for global cropland over time for different interventions and SSPs. 
#              GHG mitigation potential is unconstrained potential.
# note:        All variables expressed as difference between practice(s) and REF
#-------------------------------------------------------------------------------
# Column Headers

# cell: cell number
# x: longitude
# y: latitude
# y_block: decade (mean annual + or - 5 years) or year (cumulative)
# ssp: historical (no climate change), SSP1-2.6, SSP3-7.0
# gcm: general circultion model from CMIP6
# total_crop_area_ha: total crop area (maize, soybean, wheat) in grid cell

# s_GHG: cumulative SOC + N2O differences
# m_GHG: decadal mean annual SOC + N2O differences
# s_SOC: cumulative SOC differences
# m_SOC: decadal mean annual SOC differences
# s_N2O: cumulative N2O (direct and indirect) differences
# m_N2O: decadal mean annual N2O (direct and indirect) differences
# s_dN2O: cumulative direct N2O differences
# m_dN2O: decadal mean annual direct N2O differences
# s_iN2O: cumulative indirect N2O differences
# m_iN2O: decadal mean annual indirect N2O differences
# s_cr_grain: cumulative crop yield differences
# m_cr_grain: decadal mean annual crop yield differences

# GHG units:
# s_: Mg CO2-eq ha-1
# m_: Mg CO2-eq ha-1 yr-1

# grain units:
# s_: g C m-2
# m_: g C m-2 yr-1

# Practices
# res or residue (0G-0L-1T-1R)
# ntill or no-tillage (0G-0L-0T-0R)
# ntill-res or no-tillage + residue (0G-0L-0T-1R)
# ccg or grass cover crop (1G-0L-1T-0R)
# ccg-res or grass cover crop + residue (1G-0L-0T-1R)
# ccl or legume cover crop (0G-1L-1T-0R)
# ccl-res or legume cover crop + residue (0G-1L-1T-1R)
# ccg-ntill or grass cover crop + no-tillage + residue (1G-0L-0T-1R)
# ccl-ntill or legume cover crop + no-tillage + residue (0G-1L-0T-1R)
#-------------------------------------------------------------------------------
library(data.table)
library(factoextra)
library(rstudioapi)
library(scales)
library(sf)
library(terra)
#-------------------------------------------------------------------------------
source('results_functions.R')
#-------------------------------------------------------------------------------
options(scipen = 999, digits = 4)
options(rgl.printRglwidget = TRUE)
base.path = dirname(getActiveDocumentContext()$path)
data.path = paste(base.path, 'manuscript-data', sep = '/')
lu.path   = paste(base.path, 'gis', sep = '/')
raster    = 'msw-cropland-rf-ir-area.tif'
raster_cc = 'msw-masked-cropland-rf-ir-area.tif'
shp       = 'WB_countries_Admin0_10m.shp'
#-------------------------------------------------------------------------------
# LOAD Data
#-------------------------------------------------------------------------------
# Long to wide format

# RESIDUE
residue_dt = readRDS(paste(data.path, 'imputed-gcm-relative-responses-weighted-mean-res.rds', sep = '/'))
residue_dt = dcast(residue_dt,
                   cell + x + y + y_block + ssp + gcm + total_crop_area_ha ~ variable,
                   value.var = 'value')
# NTILL
ntill_dt = readRDS(paste(data.path, 'imputed-gcm-relative-responses-weighted-mean-ntill.rds', sep = '/'))
ntill_dt = dcast(ntill_dt,
                 cell + x + y + y_block + ssp + gcm + total_crop_area_ha ~ variable,
                 value.var = 'value')
# NTILL-RES
ntill_res_dt = readRDS(paste(data.path, 'imputed-gcm-relative-responses-weighted-mean-ntill-res.rds', sep = '/'))
ntill_res_dt = dcast(ntill_res_dt,
                     cell + x + y + y_block + ssp + gcm + total_crop_area_ha ~ variable,
                     value.var = 'value')
# CCG
ccg_dt = readRDS(paste(data.path, 'imputed-gcm-relative-responses-weighted-mean-ccg.rds', sep = '/'))
ccg_dt = dcast(ccg_dt,
               cell + x + y + y_block + ssp + gcm + total_crop_area_ha ~ variable,
               value.var = 'value')
# CCG-RES
ccg_res_dt = readRDS(paste(data.path, 'imputed-gcm-relative-responses-weighted-mean-ccg-res.rds', sep = '/'))
ccg_res_dt = dcast(ccg_res_dt,
                   cell + x + y + y_block + ssp + gcm + total_crop_area_ha ~ variable,
                   value.var = 'value')
# CCL
ccl_dt = readRDS(paste(data.path, 'imputed-gcm-relative-responses-weighted-mean-ccl.rds', sep = '/'))
ccl_dt = dcast(ccl_dt,
               cell + x + y + y_block + ssp + gcm + total_crop_area_ha ~ variable,
               value.var = 'value')
# CCL-RES
ccl_res_dt = readRDS(paste(data.path, 'imputed-gcm-relative-responses-weighted-mean-ccl-res.rds', sep = '/'))
ccl_res_dt = dcast(ccl_res_dt,
                   cell + x + y + y_block + ssp + gcm + total_crop_area_ha ~ variable,
                   value.var = 'value')
# CCG-NTILL-RES
ccg_ntill_dt = readRDS(paste(data.path, 'imputed-gcm-relative-responses-weighted-mean-ccg-ntill.rds', sep = '/'))
ccg_ntill_dt = dcast(ccg_ntill_dt,
                     cell + x + y + y_block + ssp + gcm + total_crop_area_ha ~ variable,
                     value.var = 'value')
# CCL-NTILL-RES
ccl_ntill_dt = readRDS(paste(data.path, 'imputed-gcm-relative-responses-weighted-mean-ccl-ntill.rds', sep = '/'))
ccl_ntill_dt = dcast(ccl_ntill_dt,
                     cell + x + y + y_block + ssp + gcm + total_crop_area_ha ~ variable,
                     value.var = 'value')
#-------------------------------------------------------------------------------
# ADD IPCC Region Names
#-------------------------------------------------------------------------------
IPCC_dt = ipcc_name(lu.path, shp, raster)
# RESIDUE
residue_dt  = residue_dt[IPCC_dt[, .(cell, IPCC_NAME, WB_NAME)], 
                         on = .(cell = cell) ]
residue_dt  = residue_dt[!is.na(gcm)]
# NTILL
ntill_dt    = ntill_dt[IPCC_dt[, .(cell, IPCC_NAME, WB_NAME)], 
                       on = .(cell = cell) ]
ntill_dt    = ntill_dt[!is.na(gcm)]
# NTILL-RES
ntill_res_dt= ntill_res_dt[IPCC_dt[, .(cell, IPCC_NAME, WB_NAME)], 
                           on = .(cell = cell) ]
ntill_res_dt= ntill_res_dt[!is.na(gcm)]
# CCG
ccg_dt      = ccg_dt[IPCC_dt[, .(cell, IPCC_NAME, WB_NAME)], 
                     on = .(cell = cell) ]
ccg_dt      = ccg_dt[!is.na(gcm)]
# CCG-RES
ccg_res_dt  = ccg_res_dt[IPCC_dt[, .(cell, IPCC_NAME, WB_NAME)], 
                         on = .(cell = cell) ]
ccg_res_dt  = ccg_res_dt[!is.na(gcm)]
# CCL
ccl_dt      = ccl_dt[IPCC_dt[, .(cell, IPCC_NAME, WB_NAME)], 
                     on = .(cell = cell) ]
ccl_dt      = ccl_dt[!is.na(gcm)]
# CCL-RES
ccl_res_dt  = ccl_res_dt[IPCC_dt[, .(cell, IPCC_NAME, WB_NAME)], 
                         on = .(cell = cell) ]
ccl_res_dt  = ccl_res_dt[!is.na(gcm)]
# CCG-NTILL
ccg_ntill_dt  = ccg_ntill_dt[IPCC_dt[, .(cell, IPCC_NAME, WB_NAME)], 
                             on = .(cell = cell) ]
ccg_ntill_dt  = ccg_ntill_dt[!is.na(gcm)]
# CCL-NTILL
ccl_ntill_dt  = ccl_ntill_dt[IPCC_dt[, .(cell, IPCC_NAME, WB_NAME)], 
                             on = .(cell = cell) ]
ccl_ntill_dt  = ccl_ntill_dt[!is.na(gcm)]
#-------------------------------------------------------------------------------
# CONVERT Biomass Units: g C m-2 to kg ha-1 yr-1 or Mg ha-1
#-------------------------------------------------------------------------------
Mg_ha = 100L
kg_ha = 10L
C_gr  = 0.42 # Ma et al. 2018 | value for 'reproductive organs'
# RESIDUE
residue_dt[, m_cr_grain := (m_cr_grain*kg_ha)/C_gr]
residue_dt[, s_cr_grain := (s_cr_grain/Mg_ha)/C_gr]
# NTILL
ntill_dt[, m_cr_grain := (m_cr_grain*kg_ha)/C_gr]
ntill_dt[, s_cr_grain := (s_cr_grain/Mg_ha)/C_gr]
# NTILL-RES
ntill_res_dt[, m_cr_grain := (m_cr_grain*kg_ha)/C_gr]
ntill_res_dt[, s_cr_grain := (s_cr_grain/Mg_ha)/C_gr]
# CCG
ccg_dt[, m_cr_grain := (m_cr_grain*kg_ha)/C_gr]
ccg_dt[, s_cr_grain := (s_cr_grain/Mg_ha)/C_gr]
# CCG-RES
ccg_res_dt[, m_cr_grain := (m_cr_grain*kg_ha)/C_gr]
ccg_res_dt[, s_cr_grain := (s_cr_grain/Mg_ha)/C_gr]
# CCL
ccl_dt[, m_cr_grain := (m_cr_grain*kg_ha)/C_gr]
ccl_dt[, s_cr_grain := (s_cr_grain/Mg_ha)/C_gr]
# CCL-RES
ccl_res_dt[, m_cr_grain := (m_cr_grain*kg_ha)/C_gr]
ccl_res_dt[, s_cr_grain := (s_cr_grain/Mg_ha)/C_gr]
# CCG-NTILL-RES
ccg_ntill_dt[, m_cr_grain := (m_cr_grain*kg_ha)/C_gr]
ccg_ntill_dt[, s_cr_grain := (s_cr_grain/Mg_ha)/C_gr]
# CCL-NTILL-RES
ccl_ntill_dt[, m_cr_grain := (m_cr_grain*kg_ha)/C_gr]
ccl_ntill_dt[, s_cr_grain := (s_cr_grain/Mg_ha)/C_gr]
#-------------------------------------------------------------------------------
# UNCONSTRAINED MITIGATION POTENTIAL | GLOBAL
#-------------------------------------------------------------------------------
k_cols = c('cell','ssp', 'gcm','y_block','total_crop_area_ha','GHG_area', 'GRAIN_area')
  # RESIDUE
residue_cst_dt = residue_dt[, GHG_area   := s_GHG*total_crop_area_ha]
residue_cst_dt = residue_dt[, GRAIN_area := s_cr_grain*total_crop_area_ha]
# constrained | 2100
residue_2100_dt  = residue_cst_dt[y_block == 2100,]
residue_2100_dt  = residue_2100_dt[, ..k_cols]
residue_2100_gcm = residue_2100_dt[, lapply(.SD, sum), .SDcols = c('GHG_area','GRAIN_area','total_crop_area_ha'), by = .(ssp, gcm, y_block)]
residue_2100_gcm[, scenario := 'res']
# potential | 2100
residue_2100_gl  = global_mean_CI(residue_2100_gcm, '2100')
residue_2100_gl[, scenario := 'res']
# constrained | 2050
residue_2050_dt  = residue_cst_dt[y_block == 2050,]
residue_2050_dt  = residue_2050_dt[, ..k_cols]
residue_2050_gcm = residue_2050_dt[, lapply(.SD, sum), .SDcols = c('GHG_area','GRAIN_area','total_crop_area_ha'), by = .(ssp, gcm, y_block)]
residue_2050_gcm[, scenario := 'res']
# potential | 2050
residue_2050_gl  = global_mean_CI(residue_2050_gcm, '2050')
residue_2050_gl[, scenario := 'res']

# constrained | 2030
residue_2030_dt  = residue_cst_dt[y_block == 2030,]
residue_2030_dt  = residue_2030_dt[, ..k_cols]
residue_2030_gcm = residue_2030_dt[, lapply(.SD, sum), .SDcols = c('GHG_area','GRAIN_area','total_crop_area_ha'), by = .(ssp, gcm, y_block)]
residue_2030_gcm[, scenario := 'res']
# potential | 2030
residue_2030_gl  = global_mean_CI(residue_2030_gcm, '2030')
residue_2030_gl[, scenario := 'res']

# NTILL
ntill_cst_dt = ntill_dt[, GHG_area   := s_GHG*total_crop_area_ha]
ntill_cst_dt = ntill_dt[, GRAIN_area := s_cr_grain*total_crop_area_ha]
# constrained | 2100
ntill_2100_dt  = ntill_cst_dt[y_block == 2100,]
ntill_2100_dt  = ntill_2100_dt[, ..k_cols]
ntill_2100_gcm = ntill_2100_dt[, lapply(.SD, sum), .SDcols = c('GHG_area','GRAIN_area','total_crop_area_ha'), by = .(ssp, gcm, y_block)]
ntill_2100_gcm[, scenario := 'ntill']
# potential | 2100
ntill_2100_gl  = global_mean_CI(ntill_2100_gcm, '2100')
ntill_2100_gl[, scenario := 'ntill']
# constrained | 2050
ntill_2050_dt  = ntill_cst_dt[y_block == 2050,]
ntill_2050_dt  = ntill_2050_dt[, ..k_cols]
ntill_2050_gcm = ntill_2050_dt[, lapply(.SD, sum), .SDcols = c('GHG_area','GRAIN_area','total_crop_area_ha'), by = .(ssp, gcm, y_block)]
ntill_2050_gcm[, scenario := 'ntill']
# potential | 2050
ntill_2050_gl  = global_mean_CI(ntill_2050_gcm, '2050')
ntill_2050_gl[, scenario := 'ntill']
# constrained | 2030
ntill_2030_dt  = ntill_cst_dt[y_block == 2030,]
ntill_2030_dt  = ntill_2030_dt[, ..k_cols]
ntill_2030_gcm = ntill_2030_dt[, lapply(.SD, sum), .SDcols = c('GHG_area','GRAIN_area','total_crop_area_ha'), by = .(ssp, gcm, y_block)]
ntill_2030_gcm[, scenario := 'ntill']
# potential | 2030
ntill_2030_gl  = global_mean_CI(ntill_2030_gcm, '2030')
ntill_2030_gl[, scenario := 'ntill']

# NTILL-RES
ntill_res_cst_dt = ntill_res_dt[, GHG_area   := s_GHG*total_crop_area_ha]
ntill_res_cst_dt = ntill_res_dt[, GRAIN_area := s_cr_grain*total_crop_area_ha]
# constrained | 2100
ntill_res_2100_dt  = ntill_res_cst_dt[y_block == 2100,]
ntill_res_2100_dt  = ntill_res_2100_dt[, ..k_cols]
ntill_res_2100_gcm = ntill_res_2100_dt[, lapply(.SD, sum), .SDcols = c('GHG_area','GRAIN_area','total_crop_area_ha'), by = .(ssp, gcm, y_block)]
ntill_res_2100_gcm[, scenario := 'ntill-res']
# potential | 2100
ntill_res_2100_gl  = global_mean_CI(ntill_res_2100_gcm, '2100')
ntill_res_2100_gl[, scenario := 'ntill-res']
# constrained | 2050
ntill_res_2050_dt  = ntill_res_cst_dt[y_block == 2050,]
ntill_res_2050_dt  = ntill_res_2050_dt[, ..k_cols]
ntill_res_2050_gcm = ntill_res_2050_dt[, lapply(.SD, sum), .SDcols = c('GHG_area','GRAIN_area','total_crop_area_ha'), by = .(ssp, gcm, y_block)]
ntill_res_2050_gcm[, scenario := 'ntill-res']
# potential | 2050
ntill_res_2050_gl  = global_mean_CI(ntill_res_2050_gcm, '2050')
ntill_res_2050_gl[, scenario := 'ntill-res']
# constrained | 2030
ntill_res_2030_dt  = ntill_res_cst_dt[y_block == 2030,]
ntill_res_2030_dt  = ntill_res_2030_dt[, ..k_cols]
ntill_res_2030_gcm = ntill_res_2030_dt[, lapply(.SD, sum), .SDcols = c('GHG_area','GRAIN_area','total_crop_area_ha'), by = .(ssp, gcm, y_block)]
ntill_res_2030_gcm[, scenario := 'ntill-res']
# potential | 2030
ntill_res_2030_gl  = global_mean_CI(ntill_res_2030_gcm, '2030')
ntill_res_2030_gl[, scenario := 'ntill-res']

# CCG
ccg_cst_dt = ccg_dt[, GHG_area   := s_GHG*total_crop_area_ha]
ccg_cst_dt = ccg_dt[, GRAIN_area := s_cr_grain*total_crop_area_ha]
# constrained | 2100
ccg_2100_dt  = ccg_cst_dt[y_block == 2100,]
ccg_2100_dt  = ccg_2100_dt[, ..k_cols]
ccg_2100_gcm = ccg_2100_dt[, lapply(.SD, sum), .SDcols = c('GHG_area','GRAIN_area','total_crop_area_ha'), by = .(ssp, gcm, y_block)]
ccg_2100_gcm[, scenario := 'ccg']
# potential | 2100
ccg_2100_gl  = global_mean_CI(ccg_2100_gcm, '2100')
ccg_2100_gl[, scenario := 'ccg']
# constrained | 2050
ccg_2050_dt  = ccg_cst_dt[y_block == 2050,]
ccg_2050_dt  = ccg_2050_dt[, ..k_cols]
ccg_2050_gcm = ccg_2050_dt[, lapply(.SD, sum), .SDcols = c('GHG_area','GRAIN_area','total_crop_area_ha'), by = .(ssp, gcm, y_block)]
ccg_2050_gcm[, scenario := 'ccg']
# potential | 2050
ccg_2050_gl  = global_mean_CI(ccg_2050_gcm, '2050')
ccg_2050_gl[, scenario := 'ccg']
# constrained | 2030
ccg_2030_dt  = ccg_cst_dt[y_block == 2030,]
ccg_2030_dt  = ccg_2030_dt[, ..k_cols]
ccg_2030_gcm = ccg_2030_dt[, lapply(.SD, sum), .SDcols = c('GHG_area','GRAIN_area','total_crop_area_ha'), by = .(ssp, gcm, y_block)]
ccg_2030_gcm[, scenario := 'ccg']
# potential | 2030
ccg_2030_gl  = global_mean_CI(ccg_2030_gcm, '2030')
ccg_2030_gl[, scenario := 'ccg']

# CCG-RES
ccg_res_cst_dt = ccg_res_dt[, GHG_area   := s_GHG*total_crop_area_ha]
ccg_res_cst_dt = ccg_res_dt[, GRAIN_area := s_cr_grain*total_crop_area_ha]
# constrained | 2100
ccg_res_2100_dt  = ccg_res_cst_dt[y_block == 2100,]
ccg_res_2100_dt  = ccg_res_2100_dt[, ..k_cols]
ccg_res_2100_gcm = ccg_res_2100_dt[, lapply(.SD, sum), .SDcols = c('GHG_area','GRAIN_area','total_crop_area_ha'), by = .(ssp, gcm, y_block)]
ccg_res_2100_gcm[, scenario := 'ccg-res']
# potential | 2100
ccg_res_2100_gl  = global_mean_CI(ccg_res_2100_gcm, '2100')
ccg_res_2100_gl[, scenario := 'ccg-res']
# constrained | 2050
ccg_res_2050_dt  = ccg_res_cst_dt[y_block == 2050,]
ccg_res_2050_dt  = ccg_res_2050_dt[, ..k_cols]
ccg_res_2050_gcm = ccg_res_2050_dt[, lapply(.SD, sum), .SDcols = c('GHG_area','GRAIN_area','total_crop_area_ha'), by = .(ssp, gcm, y_block)]
ccg_res_2050_gcm[, scenario := 'ccg-res']
# potential | 2050
ccg_res_2050_gl  = global_mean_CI(ccg_res_2050_gcm, '2050')
ccg_res_2050_gl[, scenario := 'ccg-res']

# SUM HIGHER / LOWER CROP YIELD AREA ONLY
ccg_res_2050_H_area = ccg_res_2050_dt[GRAIN_area >= 0, lapply(.SD, sum), .SDcols = c('total_crop_area_ha'), by = .(ssp, gcm, y_block)]
mean(ccg_res_2050_H_area[ssp %in% 'ssp126', total_crop_area_ha])
# 154,101,013
ccg_res_2050_L_area = ccg_res_2050_dt[GRAIN_area < 0, lapply(.SD, sum), .SDcols = c('total_crop_area_ha'), by = .(ssp, gcm, y_block)]
mean(ccg_res_2050_L_area[ssp %in% 'ssp126', total_crop_area_ha])
# 251,450,122

# constrained | 2030
ccg_res_2030_dt  = ccg_res_cst_dt[y_block == 2030,]
ccg_res_2030_dt  = ccg_res_2030_dt[, ..k_cols]
ccg_res_2030_gcm = ccg_res_2030_dt[, lapply(.SD, sum), .SDcols = c('GHG_area','GRAIN_area','total_crop_area_ha'), by = .(ssp, gcm, y_block)]
ccg_res_2030_gcm[, scenario := 'ccg-res']
# potential | 2030
ccg_res_2030_gl  = global_mean_CI(ccg_res_2030_gcm, '2030')
ccg_res_2030_gl[, scenario := 'ccg-res']

# CCL
ccl_cst_dt = ccl_dt[, GHG_area   := s_GHG*total_crop_area_ha]
ccl_cst_dt = ccl_dt[, GRAIN_area := s_cr_grain*total_crop_area_ha]
# constrained | 2100
ccl_2100_dt  = ccl_cst_dt[y_block == 2100,]
ccl_2100_dt  = ccl_2100_dt[, ..k_cols]
ccl_2100_gcm = ccl_2100_dt[, lapply(.SD, sum), .SDcols = c('GHG_area','GRAIN_area','total_crop_area_ha'), by = .(ssp, gcm, y_block)]
ccl_2100_gcm[, scenario := 'ccl']
# potential | 2100
ccl_2100_gl  = global_mean_CI(ccl_2100_gcm, '2100')
ccl_2100_gl[, scenario := 'ccl']
# constrained | 2050
ccl_2050_dt  = ccl_cst_dt[y_block == 2050,]
ccl_2050_dt  = ccl_2050_dt[, ..k_cols]
ccl_2050_gcm = ccl_2050_dt[, lapply(.SD, sum), .SDcols = c('GHG_area','GRAIN_area','total_crop_area_ha'), by = .(ssp, gcm, y_block)]
ccl_2050_gcm[, scenario := 'ccl']
# potential | 2050
ccl_2050_gl  = global_mean_CI(ccl_2050_gcm, '2050')
ccl_2050_gl[, scenario := 'ccl']
# constrained | 2030
ccl_2030_dt  = ccl_cst_dt[y_block == 2030,]
ccl_2030_dt  = ccl_2030_dt[, ..k_cols]
ccl_2030_gcm = ccl_2030_dt[, lapply(.SD, sum), .SDcols = c('GHG_area','GRAIN_area','total_crop_area_ha'), by = .(ssp, gcm, y_block)]
ccl_2030_gcm[, scenario := 'ccl']
# potential | 2030
ccl_2030_gl  = global_mean_CI(ccl_2030_gcm, '2030')
ccl_2030_gl[, scenario := 'ccl']

# CCL-RES
ccl_res_cst_dt = ccl_res_dt[, GHG_area   := s_GHG*total_crop_area_ha]
ccl_res_cst_dt = ccl_res_dt[, GRAIN_area := s_cr_grain*total_crop_area_ha]
# constrained | 2100
ccl_res_2100_dt  = ccl_res_cst_dt[y_block == 2100,]
ccl_res_2100_dt  = ccl_res_2100_dt[, ..k_cols]
ccl_res_2100_gcm = ccl_res_2100_dt[, lapply(.SD, sum), .SDcols = c('GHG_area','GRAIN_area','total_crop_area_ha'), by = .(ssp, gcm, y_block)]
ccl_res_2100_gcm[, scenario := 'ccl-res']
# potential | 2100
ccl_res_2100_gl  = global_mean_CI(ccl_res_2100_gcm, '2100')
ccl_res_2100_gl[, scenario := 'ccl-res']
# constrained | 2050
ccl_res_2050_dt  = ccl_res_cst_dt[y_block == 2050,]
ccl_res_2050_dt  = ccl_res_2050_dt[, ..k_cols]
ccl_res_2050_gcm = ccl_res_2050_dt[, lapply(.SD, sum), .SDcols = c('GHG_area','GRAIN_area','total_crop_area_ha'), by = .(ssp, gcm, y_block)]
ccl_res_2050_gcm[, scenario := 'ccl-res']
# potential | 2050
ccl_res_2050_gl  = global_mean_CI(ccl_res_2050_gcm, '2050')
ccl_res_2050_gl[, scenario := 'ccl-res']

# SUM HIGHER / LOWER CROP YIELD AREA ONLY
ccl_res_2050_H_area = ccl_res_2050_dt[GRAIN_area >= 0, lapply(.SD, sum), .SDcols = c('total_crop_area_ha'), by = .(ssp, gcm, y_block)]
mean(ccl_res_2050_H_area[ssp %in% 'ssp126', total_crop_area_ha])
# 276,790,093
ccl_res_2050_L_area = ccl_res_2050_dt[GRAIN_area < 0, lapply(.SD, sum), .SDcols = c('total_crop_area_ha'), by = .(ssp, gcm, y_block)]
mean(ccl_res_2050_L_area[ssp %in% 'ssp126', total_crop_area_ha])
# 128,761,042

# constrained | 2030
ccl_res_2030_dt  = ccl_res_cst_dt[y_block == 2030,]
ccl_res_2030_dt  = ccl_res_2030_dt[, ..k_cols]
ccl_res_2030_gcm = ccl_res_2030_dt[, lapply(.SD, sum), .SDcols = c('GHG_area','GRAIN_area','total_crop_area_ha'), by = .(ssp, gcm, y_block)]
ccl_res_2030_gcm[, scenario := 'ccl-res']
# potential | 2030
ccl_res_2030_gl  = global_mean_CI(ccl_res_2030_gcm, '2030')
ccl_res_2030_gl[, scenario := 'ccl-res']

# CCG-NTILL-RES
ccg_ntill_cst_dt = ccg_ntill_dt[, GHG_area   := s_GHG*total_crop_area_ha]
ccg_ntill_cst_dt = ccg_ntill_dt[, GRAIN_area := s_cr_grain*total_crop_area_ha]
# constrained | 2100
ccg_ntill_2100_dt  = ccg_ntill_cst_dt[y_block == 2100,]
ccg_ntill_2100_dt  = ccg_ntill_2100_dt[, ..k_cols]
ccg_ntill_2100_gcm = ccg_ntill_2100_dt[, lapply(.SD, sum), .SDcols = c('GHG_area','GRAIN_area','total_crop_area_ha'), by = .(ssp, gcm, y_block)]
ccg_ntill_2100_gcm[, scenario := 'ccg-ntill']
# potential | 2100
ccg_ntill_2100_gl  = global_mean_CI(ccg_ntill_2100_gcm, '2100')
ccg_ntill_2100_gl[, scenario := 'ccg-ntill']
# constrained | 2050
ccg_ntill_2050_dt  = ccg_ntill_cst_dt[y_block == 2050,]
ccg_ntill_2050_dt  = ccg_ntill_2050_dt[, ..k_cols]
ccg_ntill_2050_gcm = ccg_ntill_2050_dt[, lapply(.SD, sum), .SDcols = c('GHG_area','GRAIN_area','total_crop_area_ha'), by = .(ssp, gcm, y_block)]
ccg_ntill_2050_gcm[, scenario := 'ccg-ntill']
# potential | 2050
ccg_ntill_2050_gl  = global_mean_CI(ccg_ntill_2050_gcm, '2050')
ccg_ntill_2050_gl[, scenario := 'ccg-ntill']

# SUM HIGHER / LOWER CROP YIELD AREA ONLY
ccg_ntill_2050_H_area = ccg_ntill_2050_dt[GRAIN_area >= 0, lapply(.SD, sum), .SDcols = c('total_crop_area_ha'), by = .(ssp, gcm, y_block)]
mean(ccg_ntill_2050_H_area[ssp %in% 'ssp126', total_crop_area_ha])
# 85,408,566
ccg_ntill_2050_L_area = ccg_ntill_2050_dt[GRAIN_area < 0, lapply(.SD, sum), .SDcols = c('total_crop_area_ha'), by = .(ssp, gcm, y_block)]
mean(ccg_ntill_2050_L_area[ssp %in% 'ssp126', total_crop_area_ha])
# 320,144,911

# constrained | 2030
ccg_ntill_2030_dt  = ccg_ntill_cst_dt[y_block == 2030,]
ccg_ntill_2030_dt  = ccg_ntill_2030_dt[, ..k_cols]
ccg_ntill_2030_gcm = ccg_ntill_2030_dt[, lapply(.SD, sum), .SDcols = c('GHG_area','GRAIN_area','total_crop_area_ha'), by = .(ssp, gcm, y_block)]
ccg_ntill_2030_gcm[, scenario := 'ccg-ntill']
# potential | 2030
ccg_ntill_2030_gl  = global_mean_CI(ccg_ntill_2030_gcm, '2030')
ccg_ntill_2030_gl[, scenario := 'ccg-ntill']

# CCL-NTILL-RES
ccl_ntill_cst_dt = ccl_ntill_dt[, GHG_area   := s_GHG*total_crop_area_ha]
ccl_ntill_cst_dt = ccl_ntill_dt[, GRAIN_area := s_cr_grain*total_crop_area_ha]
# constrained | 2100
ccl_ntill_2100_dt  = ccl_ntill_cst_dt[y_block == 2100,]
ccl_ntill_2100_dt  = ccl_ntill_2100_dt[, ..k_cols]
ccl_ntill_2100_gcm = ccl_ntill_2100_dt[, lapply(.SD, sum), .SDcols = c('GHG_area','GRAIN_area','total_crop_area_ha'), by = .(ssp, gcm, y_block)]
ccl_ntill_2100_gcm[, scenario := 'ccl-ntill']
# potential | 2100
ccl_ntill_2100_gl  = global_mean_CI(ccl_ntill_2100_gcm, '2100')
ccl_ntill_2100_gl[, scenario := 'ccl-ntill']
# constrained | 2050
ccl_ntill_2050_dt  = ccl_ntill_cst_dt[y_block == 2050,]
ccl_ntill_2050_dt  = ccl_ntill_2050_dt[, ..k_cols]
ccl_ntill_2050_gcm = ccl_ntill_2050_dt[, lapply(.SD, sum), .SDcols = c('GHG_area','GRAIN_area','total_crop_area_ha'), by = .(ssp, gcm, y_block)]
ccl_ntill_2050_gcm[, scenario := 'ccl-ntill']
# potential | 2050
ccl_ntill_2050_gl  = global_mean_CI(ccl_ntill_2050_gcm, '2050')
ccl_ntill_2050_gl[, scenario := 'ccl-ntill']

# SUM HIGHER / LOWER CROP YIELD AREA ONLY
ccl_ntill_2050_H_area = ccl_ntill_2050_dt[GRAIN_area >= 0, lapply(.SD, sum), .SDcols = c('total_crop_area_ha'), by = .(ssp, gcm, y_block)]
mean(ccl_ntill_2050_H_area[ssp %in% 'ssp126', total_crop_area_ha])
# 238,095,561
ccl_ntill_2050_L_area = ccl_ntill_2050_dt[GRAIN_area < 0, lapply(.SD, sum), .SDcols = c('total_crop_area_ha'), by = .(ssp, gcm, y_block)]
mean(ccl_ntill_2050_L_area[ssp %in% 'ssp126', total_crop_area_ha])
# 167,457,916

# constrained | 2030
ccl_ntill_2030_dt  = ccl_ntill_cst_dt[y_block == 2030,]
ccl_ntill_2030_dt  = ccl_ntill_2030_dt[, ..k_cols]
ccl_ntill_2030_gcm = ccl_ntill_2030_dt[, lapply(.SD, sum), .SDcols = c('GHG_area','GRAIN_area','total_crop_area_ha'), by = .(ssp, gcm, y_block)]
ccl_ntill_2030_gcm[, scenario := 'ccl-ntill']
# potential | 2030
ccl_ntill_2030_gl  = global_mean_CI(ccl_ntill_2030_gcm, '2030')
ccl_ntill_2030_gl[, scenario := 'ccl-ntill']

# COMBINE & SAVE DT | ENSEMBLE
# 2030
unconstr_2030 = rbind(residue_2030_gl, ntill_2030_gl, ntill_res_2030_gl,
                       ccg_2030_gl, ccg_res_2030_gl, ccl_2030_gl, ccl_res_2030_gl,
                       ccg_ntill_2030_gl, ccl_ntill_2030_gl)
setcolorder(unconstr_2030, c('ssp', 'y_block', 'scenario'))
fwrite(unconstr_2030, file = paste(data.path, 'global-unconstrained-GHG-mitigation-potential-2030.csv', sep = '/'))
# 2050
unconstr_2050 = rbind(residue_2050_gl, ntill_2050_gl, ntill_res_2050_gl,
                       ccg_2050_gl, ccg_res_2050_gl, ccl_2050_gl, ccl_res_2050_gl,
                       ccg_ntill_2050_gl, ccl_ntill_2050_gl)
setcolorder(unconstr_2050, c('ssp', 'y_block', 'scenario'))
fwrite(unconstr_2050, file = paste(data.path, 'global-unconstrained-GHG-mitigation-potential-2050.csv', sep = '/'))

# 2100
unconstr_2100 = rbind(residue_2100_gl, ntill_2100_gl, ntill_res_2100_gl,
                      ccg_2100_gl, ccg_res_2100_gl, ccl_2100_gl, ccl_res_2100_gl,
                      ccg_ntill_2100_gl, ccl_ntill_2100_gl)
setcolorder(unconstr_2100, c('ssp', 'y_block', 'scenario'))
fwrite(unconstr_2100, file = paste(data.path, 'global-unconstrained-GHG-mitigation-potential-2100.csv', sep = '/'))

# COMBINE & SAVE DT | GCM
# 2030
unconstr_gcm_2030 = rbind(residue_2030_gcm, ntill_2030_gcm, ntill_res_2030_gcm,
                     ccg_2030_gcm, ccg_res_2030_gcm, ccl_2030_gcm, ccl_res_2030_gcm,
                     ccg_ntill_2030_gcm, ccl_ntill_2030_gcm)
setcolorder(unconstr_gcm_2030, c('ssp', 'y_block', 'scenario'))
fwrite(unconstr_2030, file = paste(data.path, 'global-unconstrained-GHG-mitigation-potential-gcm-2030.csv', sep = '/'))
# 2050
unconstr_gcm_2050 = rbind(residue_2050_gcm, ntill_2050_gcm, ntill_res_2050_gcm,
                     ccg_2050_gcm, ccg_res_2050_gcm, ccl_2050_gcm, ccl_res_2050_gcm,
                     ccg_ntill_2050_gcm, ccl_ntill_2050_gcm)
setcolorder(unconstr_gcm_2050, c('ssp', 'y_block', 'scenario'))
fwrite(unconstr_gcm_2050, file = paste(data.path, 'global-unconstrained-GHG-mitigation-potential-gcm-2050.csv', sep = '/'))
# 2100
unconstr_gcm_2100 = rbind(residue_2100_gcm, ntill_2100_gcm, ntill_res_2100_gcm,
                          ccg_2100_gcm, ccg_res_2100_gcm, ccl_2100_gcm, ccl_res_2100_gcm,
                          ccg_ntill_2100_gcm, ccl_ntill_2100_gcm)
setcolorder(unconstr_gcm_2100, c('ssp', 'y_block', 'scenario'))
fwrite(unconstr_gcm_2100, file = paste(data.path, 'global-unconstrained-GHG-mitigation-potential-gcm-2100.csv', sep = '/'))
#-------------------------------------------------------------------------------
# UNCONSTRAINED MITIGATION POTENTIAL | REGION
#-------------------------------------------------------------------------------
k_cols = c('cell','ssp', 'gcm','IPCC_NAME','y_block','total_crop_area_ha','GHG_area', 'GRAIN_area')
# RESIDUE
residue_cst_dt = residue_dt[, GHG_area   := s_GHG*total_crop_area_ha]
residue_cst_dt = residue_dt[, GRAIN_area := s_cr_grain*total_crop_area_ha]
# constrained | 2100
residue_2100_dt    = residue_cst_dt[y_block == 2100,]
residue_2100_dt    = residue_2100_dt[, ..k_cols]
residue_2100_r_gcm = residue_2100_dt[, lapply(.SD, sum), .SDcols = c('GHG_area','GRAIN_area','total_crop_area_ha'), by = .(ssp, gcm, y_block, IPCC_NAME)]
residue_2100_r_gcm[, scenario := 'res']
# potential | 2100
residue_2100_rg    = regional_mean_CI(residue_2100_r_gcm, '2100')
residue_2100_rg[, scenario := 'res']
# constrained | 2050
residue_2050_dt    = residue_cst_dt[y_block == 2050,]
residue_2050_dt    = residue_2050_dt[, ..k_cols]
residue_2050_r_gcm = residue_2050_dt[, lapply(.SD, sum), .SDcols = c('GHG_area','GRAIN_area','total_crop_area_ha'), by = .(ssp, gcm, y_block, IPCC_NAME)]
residue_2050_r_gcm[, scenario := 'res']
# potential | 2050
residue_2050_rg    = regional_mean_CI(residue_2050_r_gcm, '2050')
residue_2050_rg[, scenario := 'res']
# constrained | 2030
residue_2030_dt    = residue_cst_dt[y_block == 2030,]
residue_2030_dt    = residue_2030_dt[, ..k_cols]
residue_2030_r_gcm = residue_2030_dt[, lapply(.SD, sum), .SDcols = c('GHG_area','GRAIN_area','total_crop_area_ha'), by = .(ssp, gcm, y_block, IPCC_NAME)]
residue_2030_r_gcm[, scenario := 'res']
# potential | 2030
residue_2030_rg    = regional_mean_CI(residue_2030_r_gcm, '2030')
residue_2030_rg[, scenario := 'res']

# NTILL
ntill_cst_dt = ntill_dt[, GHG_area   := s_GHG*total_crop_area_ha]
ntill_cst_dt = ntill_dt[, GRAIN_area := s_cr_grain*total_crop_area_ha]
# constrained | 2100
ntill_2100_dt    = ntill_cst_dt[y_block == 2100,]
ntill_2100_dt    = ntill_2100_dt[, ..k_cols]
ntill_2100_r_gcm = ntill_2100_dt[, lapply(.SD, sum), .SDcols = c('GHG_area','GRAIN_area','total_crop_area_ha'), by = .(ssp, gcm, y_block, IPCC_NAME)]
ntill_2100_r_gcm[, scenario := 'ntill']
# potential | 2100
ntill_2100_rg    = regional_mean_CI(ntill_2100_r_gcm, '2100')
ntill_2100_rg[, scenario := 'ntill']
# constrained | 2050
ntill_2050_dt    = ntill_cst_dt[y_block == 2050,]
ntill_2050_dt    = ntill_2050_dt[, ..k_cols]
ntill_2050_r_gcm = ntill_2050_dt[, lapply(.SD, sum), .SDcols = c('GHG_area','GRAIN_area','total_crop_area_ha'), by = .(ssp, gcm, y_block, IPCC_NAME)]
ntill_2050_r_gcm[, scenario := 'ntill']
# potential | 2050
ntill_2050_rg    = regional_mean_CI(ntill_2050_r_gcm, '2050')
ntill_2050_rg[, scenario := 'ntill']
# constrained | 2030
ntill_2030_dt    = ntill_cst_dt[y_block == 2030,]
ntill_2030_dt    = ntill_2030_dt[, ..k_cols]
ntill_2030_r_gcm = ntill_2030_dt[, lapply(.SD, sum), .SDcols = c('GHG_area','GRAIN_area','total_crop_area_ha'), by = .(ssp, gcm, y_block, IPCC_NAME)]
ntill_2030_r_gcm[, scenario := 'ntill']
# potential | 2030
ntill_2030_rg    = regional_mean_CI(ntill_2030_r_gcm, '2030')
ntill_2030_rg[, scenario := 'ntill']

# NTILL-RES
ntill_res_cst_dt = ntill_res_dt[, GHG_area   := s_GHG*total_crop_area_ha]
ntill_res_cst_dt = ntill_res_dt[, GRAIN_area := s_cr_grain*total_crop_area_ha]
# constrained | 2100
ntill_res_2100_dt    = ntill_res_cst_dt[y_block == 2100,]
ntill_res_2100_dt    = ntill_res_2100_dt[, ..k_cols]
ntill_res_2100_r_gcm = ntill_res_2100_dt[, lapply(.SD, sum), .SDcols = c('GHG_area','GRAIN_area','total_crop_area_ha'), by = .(ssp, gcm, y_block, IPCC_NAME)]
ntill_res_2100_r_gcm[, scenario := 'ntill-res']
# potential | 2100
ntill_res_2100_rg  = regional_mean_CI(ntill_res_2100_r_gcm, '2100')
ntill_res_2100_rg[, scenario := 'ntill-res']
# constrained | 2050
ntill_res_2050_dt    = ntill_res_cst_dt[y_block == 2050,]
ntill_res_2050_dt    = ntill_res_2050_dt[, ..k_cols]
ntill_res_2050_r_gcm = ntill_res_2050_dt[, lapply(.SD, sum), .SDcols = c('GHG_area','GRAIN_area','total_crop_area_ha'), by = .(ssp, gcm, y_block, IPCC_NAME)]
ntill_res_2050_r_gcm[, scenario := 'ntill-res']
# potential | 2050
ntill_res_2050_rg  = regional_mean_CI(ntill_res_2050_r_gcm, '2050')
ntill_res_2050_rg[, scenario := 'ntill-res']
# constrained | 2030
ntill_res_2030_dt    = ntill_res_cst_dt[y_block == 2030,]
ntill_res_2030_dt    = ntill_res_2030_dt[, ..k_cols]
ntill_res_2030_r_gcm = ntill_res_2030_dt[, lapply(.SD, sum), .SDcols = c('GHG_area','GRAIN_area','total_crop_area_ha'), by = .(ssp, gcm, y_block, IPCC_NAME)]
ntill_res_2030_r_gcm[, scenario := 'ntill-res']
# potential | 2030
ntill_res_2030_rg    = regional_mean_CI(ntill_res_2030_r_gcm, '2030')
ntill_res_2030_rg[, scenario := 'ntill-res']

# CCG
ccg_cst_dt = ccg_dt[, GHG_area   := s_GHG*total_crop_area_ha]
ccg_cst_dt = ccg_dt[, GRAIN_area := s_cr_grain*total_crop_area_ha]
# constrained | 2100
ccg_2100_dt    = ccg_cst_dt[y_block == 2100,]
ccg_2100_dt    = ccg_2100_dt[, ..k_cols]
ccg_2100_r_gcm = ccg_2100_dt[, lapply(.SD, sum), .SDcols = c('GHG_area','GRAIN_area','total_crop_area_ha'), by = .(ssp, gcm, y_block, IPCC_NAME)]
ccg_2100_r_gcm[, scenario := 'ccg']
# potential | 2100
ccg_2100_rg    = regional_mean_CI(ccg_2100_r_gcm, '2100')
ccg_2100_rg[, scenario := 'ccg']
# constrained | 2050
ccg_2050_dt    = ccg_cst_dt[y_block == 2050,]
ccg_2050_dt    = ccg_2050_dt[, ..k_cols]
ccg_2050_r_gcm = ccg_2050_dt[, lapply(.SD, sum), .SDcols = c('GHG_area','GRAIN_area','total_crop_area_ha'), by = .(ssp, gcm, y_block, IPCC_NAME)]
ccg_2050_r_gcm[, scenario := 'ccg']
# potential | 2050
ccg_2050_rg    = regional_mean_CI(ccg_2050_r_gcm, '2050')
ccg_2050_rg[, scenario := 'ccg']
# constrained | 2030
ccg_2030_dt    = ccg_cst_dt[y_block == 2030,]
ccg_2030_dt    = ccg_2030_dt[, ..k_cols]
ccg_2030_r_gcm = ccg_2030_dt[, lapply(.SD, sum), .SDcols = c('GHG_area','GRAIN_area','total_crop_area_ha'), by = .(ssp, gcm, y_block, IPCC_NAME)]
ccg_2030_r_gcm[, scenario := 'ccg']
# potential | 2030
ccg_2030_rg    = regional_mean_CI(ccg_2030_r_gcm, '2030')
ccg_2030_rg[, scenario := 'ccg']

# CCG-RES
ccg_res_cst_dt = ccg_res_dt[, GHG_area   := s_GHG*total_crop_area_ha]
ccg_res_cst_dt = ccg_res_dt[, GRAIN_area := s_cr_grain*total_crop_area_ha]
# constrained | 2100
ccg_res_2100_dt    = ccg_res_cst_dt[y_block == 2100,]
ccg_res_2100_dt    = ccg_res_2100_dt[, ..k_cols]
ccg_res_2100_r_gcm = ccg_res_2100_dt[, lapply(.SD, sum), .SDcols = c('GHG_area','GRAIN_area','total_crop_area_ha'), by = .(ssp, gcm, y_block, IPCC_NAME)]
ccg_res_2100_r_gcm[, scenario := 'ccg-res']
# potential | 2100
ccg_res_2100_rg  = regional_mean_CI(ccg_res_2100_r_gcm, '2100')
ccg_res_2100_rg[, scenario := 'ccg-res']
# constrained | 2050
ccg_res_2050_dt    = ccg_res_cst_dt[y_block == 2050,]
ccg_res_2050_dt    = ccg_res_2050_dt[, ..k_cols]
ccg_res_2050_r_gcm = ccg_res_2050_dt[, lapply(.SD, sum), .SDcols = c('GHG_area','GRAIN_area','total_crop_area_ha'), by = .(ssp, gcm, y_block, IPCC_NAME)]
ccg_res_2050_r_gcm[, scenario := 'ccg-res']
# potential | 2050
ccg_res_2050_rg  = regional_mean_CI(ccg_res_2050_r_gcm, '2050')
ccg_res_2050_rg[, scenario := 'ccg-res']
# constrained | 2030
ccg_res_2030_dt    = ccg_res_cst_dt[y_block == 2030,]
ccg_res_2030_dt    = ccg_res_2030_dt[, ..k_cols]
ccg_res_2030_r_gcm = ccg_res_2030_dt[, lapply(.SD, sum), .SDcols = c('GHG_area','GRAIN_area','total_crop_area_ha'), by = .(ssp, gcm, y_block, IPCC_NAME)]
ccg_res_2030_r_gcm[, scenario := 'ccg-res']
# potential | 2030
ccg_res_2030_rg    = regional_mean_CI(ccg_res_2030_r_gcm, '2030')
ccg_res_2030_rg[, scenario := 'ccg-res']

# CCL
ccl_cst_dt = ccl_dt[, GHG_area   := s_GHG*total_crop_area_ha]
ccl_cst_dt = ccl_dt[, GRAIN_area := s_cr_grain*total_crop_area_ha]
# constrained | 2100
ccl_2100_dt    = ccl_cst_dt[y_block == 2100,]
ccl_2100_dt    = ccl_2100_dt[, ..k_cols]
ccl_2100_r_gcm = ccl_2100_dt[, lapply(.SD, sum), .SDcols = c('GHG_area','GRAIN_area','total_crop_area_ha'), by = .(ssp, gcm, y_block, IPCC_NAME)]
ccl_2100_r_gcm[, scenario := 'ccl']
# potential | 2100
ccl_2100_rg    = regional_mean_CI(ccl_2100_r_gcm, '2100')
ccl_2100_rg[, scenario := 'ccl']
# constrained | 2050
ccl_2050_dt    = ccl_cst_dt[y_block == 2050,]
ccl_2050_dt    = ccl_2050_dt[, ..k_cols]
ccl_2050_r_gcm = ccl_2050_dt[, lapply(.SD, sum), .SDcols = c('GHG_area','GRAIN_area','total_crop_area_ha'), by = .(ssp, gcm, y_block, IPCC_NAME)]
ccl_2050_r_gcm[, scenario := 'ccl']
# potential | 2050
ccl_2050_rg    = regional_mean_CI(ccl_2050_r_gcm, '2050')
ccl_2050_rg[, scenario := 'ccl']
# constrained | 2030
ccl_2030_dt    = ccl_cst_dt[y_block == 2030,]
ccl_2030_dt    = ccl_2030_dt[, ..k_cols]
ccl_2030_r_gcm = ccl_2030_dt[, lapply(.SD, sum), .SDcols = c('GHG_area','GRAIN_area','total_crop_area_ha'), by = .(ssp, gcm, y_block, IPCC_NAME)]
ccl_2030_r_gcm[, scenario := 'ccl']
# potential | 2030
ccl_2030_rg    = regional_mean_CI(ccl_2030_r_gcm, '2030')
ccl_2030_rg[, scenario := 'ccl']

# CCL-RES
ccl_res_cst_dt = ccl_res_dt[, GHG_area   := s_GHG*total_crop_area_ha]
ccl_res_cst_dt = ccl_res_dt[, GRAIN_area := s_cr_grain*total_crop_area_ha]
# constrained | 2100
ccl_res_2100_dt    = ccl_res_cst_dt[y_block == 2100,]
ccl_res_2100_dt    = ccl_res_2100_dt[, ..k_cols]
ccl_res_2100_r_gcm = ccl_res_2100_dt[, lapply(.SD, sum), .SDcols = c('GHG_area','GRAIN_area','total_crop_area_ha'), by = .(ssp, gcm, y_block, IPCC_NAME)]
ccl_res_2100_r_gcm[, scenario := 'ccl-res']
# potential | 2100
ccl_res_2100_rg  = regional_mean_CI(ccl_res_2100_r_gcm, '2100')
ccl_res_2100_rg[, scenario := 'ccl-res']
# constrained | 2050
ccl_res_2050_dt    = ccl_res_cst_dt[y_block == 2050,]
ccl_res_2050_dt    = ccl_res_2050_dt[, ..k_cols]
ccl_res_2050_r_gcm = ccl_res_2050_dt[, lapply(.SD, sum), .SDcols = c('GHG_area','GRAIN_area','total_crop_area_ha'), by = .(ssp, gcm, y_block, IPCC_NAME)]
ccl_res_2050_r_gcm[, scenario := 'ccl-res']
# potential | 2050
ccl_res_2050_rg  = regional_mean_CI(ccl_res_2050_r_gcm, '2050')
ccl_res_2050_rg[, scenario := 'ccl-res']
# constrained | 2030
ccl_res_2030_dt    = ccl_res_cst_dt[y_block == 2030,]
ccl_res_2030_dt    = ccl_res_2030_dt[, ..k_cols]
ccl_res_2030_r_gcm = ccl_res_2030_dt[, lapply(.SD, sum), .SDcols = c('GHG_area','GRAIN_area','total_crop_area_ha'), by = .(ssp, gcm, y_block, IPCC_NAME)]
ccl_res_2030_r_gcm[, scenario := 'ccl-res']
# potential | 2030
ccl_res_2030_rg    = regional_mean_CI(ccl_res_2030_r_gcm, '2030')
ccl_res_2030_rg[, scenario := 'ccl-res']

# CCG-NTILL-RES
ccg_ntill_cst_dt = ccg_ntill_dt[, GHG_area   := s_GHG*total_crop_area_ha]
ccg_ntill_cst_dt = ccg_ntill_dt[, GRAIN_area := s_cr_grain*total_crop_area_ha]
# constrained | 2100
ccg_ntill_2100_dt    = ccg_ntill_cst_dt[y_block == 2100,]
ccg_ntill_2100_dt    = ccg_ntill_2100_dt[, ..k_cols]
ccg_ntill_2100_r_gcm = ccg_ntill_2100_dt[, lapply(.SD, sum), .SDcols = c('GHG_area','GRAIN_area','total_crop_area_ha'), by = .(ssp, gcm, y_block, IPCC_NAME)]
ccg_ntill_2100_r_gcm[, scenario := 'ccg-ntill']
# potential | 2100
ccg_ntill_2100_rg  = regional_mean_CI(ccg_ntill_2100_r_gcm, '2100')
ccg_ntill_2100_rg[, scenario := 'ccg-ntill']
# constrained | 2050
ccg_ntill_2050_dt    = ccg_ntill_cst_dt[y_block == 2050,]
ccg_ntill_2050_dt    = ccg_ntill_2050_dt[, ..k_cols]
ccg_ntill_2050_r_gcm = ccg_ntill_2050_dt[, lapply(.SD, sum), .SDcols = c('GHG_area','GRAIN_area','total_crop_area_ha'), by = .(ssp, gcm, y_block, IPCC_NAME)]
ccg_ntill_2050_r_gcm[, scenario := 'ccg-ntill']
# potential | 2050
ccg_ntill_2050_rg  = regional_mean_CI(ccg_ntill_2050_r_gcm, '2050')
ccg_ntill_2050_rg[, scenario := 'ccg-ntill']
# constrained | 2030
ccg_ntill_2030_dt    = ccg_ntill_cst_dt[y_block == 2030,]
ccg_ntill_2030_dt    = ccg_ntill_2030_dt[, ..k_cols]
ccg_ntill_2030_r_gcm = ccg_ntill_2030_dt[, lapply(.SD, sum), .SDcols = c('GHG_area','GRAIN_area','total_crop_area_ha'), by = .(ssp, gcm, y_block, IPCC_NAME)]
ccg_ntill_2030_r_gcm[, scenario := 'ccg-ntill']
# potential | 2030
ccg_ntill_2030_rg    = regional_mean_CI(ccg_ntill_2030_r_gcm, '2030')
ccg_ntill_2030_rg[, scenario := 'ccg-ntill']

# CCL-NTILL-RES
ccl_ntill_cst_dt = ccl_ntill_dt[, GHG_area   := s_GHG*total_crop_area_ha]
ccl_ntill_cst_dt = ccl_ntill_dt[, GRAIN_area := s_cr_grain*total_crop_area_ha]
# constrained | 2100
ccl_ntill_2100_dt    = ccl_ntill_cst_dt[y_block == 2100,]
ccl_ntill_2100_dt    = ccl_ntill_2100_dt[, ..k_cols]
ccl_ntill_2100_r_gcm = ccl_ntill_2100_dt[, lapply(.SD, sum), .SDcols = c('GHG_area','GRAIN_area','total_crop_area_ha'), by = .(ssp, gcm, y_block, IPCC_NAME)]
ccl_ntill_2100_r_gcm[, scenario := 'ccl-ntill']
# potential | 2100
ccl_ntill_2100_rg  = regional_mean_CI(ccl_ntill_2100_r_gcm, '2100')
ccl_ntill_2100_rg[, scenario := 'ccl-ntill']
# constrained | 2050
ccl_ntill_2050_dt    = ccl_ntill_cst_dt[y_block == 2050,]
ccl_ntill_2050_dt    = ccl_ntill_2050_dt[, ..k_cols]
ccl_ntill_2050_r_gcm = ccl_ntill_2050_dt[, lapply(.SD, sum), .SDcols = c('GHG_area','GRAIN_area','total_crop_area_ha'), by = .(ssp, gcm, y_block, IPCC_NAME)]
ccl_ntill_2050_r_gcm[, scenario := 'ccl-ntill']
# potential | 2050
ccl_ntill_2050_rg  = regional_mean_CI(ccl_ntill_2050_r_gcm, '2050')
ccl_ntill_2050_rg[, scenario := 'ccl-ntill']
# constrained | 2030
ccl_ntill_2030_dt    = ccl_ntill_cst_dt[y_block == 2030,]
ccl_ntill_2030_dt    = ccl_ntill_2030_dt[, ..k_cols]
ccl_ntill_2030_r_gcm = ccl_ntill_2030_dt[, lapply(.SD, sum), .SDcols = c('GHG_area','GRAIN_area','total_crop_area_ha'), by = .(ssp, gcm, y_block, IPCC_NAME)]
ccl_ntill_2030_r_gcm[, scenario := 'ccl-ntill']
# potential | 2030
ccl_ntill_2030_rg    = regional_mean_CI(ccl_ntill_2030_r_gcm, '2030')
ccl_ntill_2030_rg[, scenario := 'ccl-ntill']

# COMBINE & SAVE DT | ENSEMBLE
# 2030
ipcc_unconstr_2030 = rbind(residue_2030_rg, ntill_2030_rg, ntill_res_2030_rg,
                            ccg_2030_rg, ccg_res_2030_rg, ccl_2030_rg, ccl_res_2030_rg,
                            ccg_ntill_2030_rg, ccl_ntill_2030_rg)
setcolorder(ipcc_unconstr_2030, c('IPCC_NAME','ssp', 'y_block', 'scenario'))
fwrite(ipcc_unconstr_2030, file = paste(data.path, 'regional-unconstrained-GHG-mitigation-potential-2030.csv', sep = '/'))
# 2050
ipcc_unconstr_2050 = rbind(residue_2050_rg, ntill_2050_rg, ntill_res_2050_rg,
                            ccg_2050_rg, ccg_res_2050_rg, ccl_2050_rg, ccl_res_2050_rg,
                            ccg_ntill_2050_rg, ccl_ntill_2050_rg)
setcolorder(ipcc_unconstr_2050, c('IPCC_NAME','ssp', 'y_block', 'scenario'))
fwrite(ipcc_unconstr_2050, file = paste(data.path, 'regional-unconstrained-GHG-mitigation-potential-2050.csv', sep = '/'))
# 2100
ipcc_unconstr_2100 = rbind(residue_2100_rg, ntill_2100_rg, ntill_res_2100_rg,
                           ccg_2100_rg, ccg_res_2100_rg, ccl_2100_rg, ccl_res_2100_rg,
                           ccg_ntill_2100_rg, ccl_ntill_2100_rg)
setcolorder(ipcc_unconstr_2100, c('IPCC_NAME','ssp', 'y_block', 'scenario'))
fwrite(ipcc_unconstr_2100, file = paste(data.path, 'regional-unconstrained-GHG-mitigation-potential-2100.csv', sep = '/'))

# COMBINE & SAVE DT | GCM
# 2030
ipcc_unconstr_gcm_2030 = rbind(residue_2030_r_gcm, ntill_2030_r_gcm, ntill_res_2030_r_gcm,
                           ccg_2030_r_gcm, ccg_res_2030_r_gcm, ccl_2030_r_gcm, ccl_res_2030_r_gcm,
                           ccg_ntill_2030_r_gcm, ccl_ntill_2030_r_gcm)
setcolorder(ipcc_unconstr_gcm_2030, c('IPCC_NAME','ssp', 'y_block', 'scenario'))
fwrite(ipcc_unconstr_gcm_2030, file = paste(data.path, 'regional-unconstrained-GHG-mitigation-potential-gcm-2030.csv', sep = '/'))
# 2050
ipcc_unconstr_gcm_2050 = rbind(residue_2050_r_gcm, ntill_2050_r_gcm, ntill_res_2050_r_gcm,
                           ccg_2050_r_gcm, ccg_res_2050_r_gcm, ccl_2050_r_gcm, ccl_res_2050_r_gcm,
                           ccg_ntill_2050_r_gcm, ccl_ntill_2050_r_gcm)
setcolorder(ipcc_unconstr_gcm_2050, c('IPCC_NAME','ssp', 'y_block', 'scenario'))
fwrite(ipcc_unconstr_gcm_2050, file = paste(data.path, 'regional-unconstrained-GHG-mitigation-potential-gcm-2050.csv', sep = '/'))
# 2100
ipcc_unconstr_gcm_2100 = rbind(residue_2100_r_gcm, ntill_2100_r_gcm, ntill_res_2100_r_gcm,
                               ccg_2100_r_gcm, ccg_res_2100_r_gcm, ccl_2100_r_gcm, ccl_res_2100_r_gcm,
                               ccg_ntill_2100_r_gcm, ccl_ntill_2100_r_gcm)
setcolorder(ipcc_unconstr_gcm_2100, c('IPCC_NAME','ssp', 'y_block', 'scenario'))
fwrite(ipcc_unconstr_gcm_2100, file = paste(data.path, 'regional-unconstrained-GHG-mitigation-potential-gcm-2100.csv', sep = '/'))
