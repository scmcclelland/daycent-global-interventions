# filename:    analysis_biophysical.R
# created:     23 February 2024
# updated:     26 February 2024
# author:      S.C. McClelland
# description: This file analyzes GHG mitigation potential from DayCent simulations
#              for global cropland over time for different interventions and SSPs. 
#              GHG mitigation potential is biophysical potential.
#-------------------------------------------------------------------------------
library(data.table)
library(factoextra)
library(patchwork)
library(RColorBrewer)
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
# CCG-NTILL-RES
ccg_ntill_dt  = ccg_ntill_dt[IPCC_dt[, .(cell, IPCC_NAME, WB_NAME)], 
                             on = .(cell = cell) ]
ccg_ntill_dt  = ccg_ntill_dt[!is.na(gcm)]
# CCL-NTILL-RES
ccl_ntill_dt  = ccl_ntill_dt[IPCC_dt[, .(cell, IPCC_NAME, WB_NAME)], 
                             on = .(cell = cell) ]
ccl_ntill_dt  = ccl_ntill_dt[!is.na(gcm)]
#-------------------------------------------------------------------------------
# CONVERT Biomass Units
#-------------------------------------------------------------------------------
Mg_ha = 100L
kg_ha = 10L
C_bio = 0.45
# RESIDUE
residue_dt[, m_cr_grain := (m_cr_grain*kg_ha)/C_bio]
residue_dt[, s_cr_grain := (s_cr_grain/Mg_ha)/C_bio]
# NTILL
ntill_dt[, m_cr_grain := (m_cr_grain*kg_ha)/C_bio]
ntill_dt[, s_cr_grain := (s_cr_grain/Mg_ha)/C_bio]
# NTILL-RES
ntill_res_dt[, m_cr_grain := (m_cr_grain*kg_ha)/C_bio]
ntill_res_dt[, s_cr_grain := (s_cr_grain/Mg_ha)/C_bio]
# CCG
ccg_dt[, m_cr_grain := (m_cr_grain*kg_ha)/C_bio]
ccg_dt[, s_cr_grain := (s_cr_grain/Mg_ha)/C_bio]
# CCG-RES
ccg_res_dt[, m_cr_grain := (m_cr_grain*kg_ha)/C_bio]
ccg_res_dt[, s_cr_grain := (s_cr_grain/Mg_ha)/C_bio]
# CCL
ccl_dt[, m_cr_grain := (m_cr_grain*kg_ha)/C_bio]
ccl_dt[, s_cr_grain := (s_cr_grain/Mg_ha)/C_bio]
# CCL-RES
ccl_res_dt[, m_cr_grain := (m_cr_grain*kg_ha)/C_bio]
ccl_res_dt[, s_cr_grain := (s_cr_grain/Mg_ha)/C_bio]
# CCG-NTILL-RES
ccg_ntill_dt[, m_cr_grain := (m_cr_grain*kg_ha)/C_bio]
ccg_ntill_dt[, s_cr_grain := (s_cr_grain/Mg_ha)/C_bio]
# CCG-NTILL-RES
ccl_ntill_dt[, m_cr_grain := (m_cr_grain*kg_ha)/C_bio]
ccl_ntill_dt[, s_cr_grain := (s_cr_grain/Mg_ha)/C_bio]
#-------------------------------------------------------------------------------
# BIOPHYSICAL MITIGATION POTENTIAL | GLOBAL
#-------------------------------------------------------------------------------
k_cols = c('cell','ssp', 'gcm','y_block','total_crop_area_ha','GHG_area', 'GRAIN_area')
  # RESIDUE
residue_cst_dt = residue_dt[, GHG_area   := s_GHG*total_crop_area_ha]
residue_cst_dt = residue_dt[, GRAIN_area := s_cr_grain*total_crop_area_ha]
# constrained | 2050
residue_2050_dt  = residue_cst_dt[y_block == 2050,]
residue_2050_dt  = residue_2050_dt[, ..k_cols]
residue_2050_gcm = residue_2050_dt[, lapply(.SD, sum), .SDcols = c('GHG_area','GRAIN_area','total_crop_area_ha'), by = .(ssp, gcm, y_block)]
# potential | 2050
residue_2050_gl  = global_mean_CI(residue_2050_gcm, '2050')
residue_2050_gl[, scenario := 'res']

# constrained | 2030
residue_2030_dt  = residue_cst_dt[y_block == 2030,]
residue_2030_dt  = residue_2030_dt[, ..k_cols]
residue_2030_gcm = residue_2030_dt[, lapply(.SD, sum), .SDcols = c('GHG_area','GRAIN_area','total_crop_area_ha'), by = .(ssp, gcm, y_block)]
# potential | 2030
residue_2030_gl  = global_mean_CI(residue_2030_gcm, '2030')
residue_2030_gl[, scenario := 'res']

# NTILL
ntill_cst_dt = ntill_dt[, GHG_area   := s_GHG*total_crop_area_ha]
ntill_cst_dt = ntill_dt[, GRAIN_area := s_cr_grain*total_crop_area_ha]
# constrained | 2050
ntill_2050_dt  = ntill_cst_dt[y_block == 2050,]
ntill_2050_dt  = ntill_2050_dt[, ..k_cols]
ntill_2050_gcm = ntill_2050_dt[, lapply(.SD, sum), .SDcols = c('GHG_area','GRAIN_area','total_crop_area_ha'), by = .(ssp, gcm, y_block)]
# potential | 2050
ntill_2050_gl  = global_mean_CI(ntill_2050_gcm, '2050')
ntill_2050_gl[, scenario := 'ntill']

# constrained | 2030
ntill_2030_dt  = ntill_cst_dt[y_block == 2030,]
ntill_2030_dt  = ntill_2030_dt[, ..k_cols]
ntill_2030_gcm = ntill_2030_dt[, lapply(.SD, sum), .SDcols = c('GHG_area','GRAIN_area','total_crop_area_ha'), by = .(ssp, gcm, y_block)]
# potential | 2030
ntill_2030_gl  = global_mean_CI(ntill_2030_gcm, '2030')
ntill_2030_gl[, scenario := 'ntill']

# NTILL-RES
ntill_res_cst_dt = ntill_res_dt[, GHG_area   := s_GHG*total_crop_area_ha]
ntill_res_cst_dt = ntill_res_dt[, GRAIN_area := s_cr_grain*total_crop_area_ha]
# constrained | 2050
ntill_res_2050_dt  = ntill_res_cst_dt[y_block == 2050,]
ntill_res_2050_dt  = ntill_res_2050_dt[, ..k_cols]
ntill_res_2050_gcm = ntill_res_2050_dt[, lapply(.SD, sum), .SDcols = c('GHG_area','GRAIN_area','total_crop_area_ha'), by = .(ssp, gcm, y_block)]
# potential | 2050
ntill_res_2050_gl  = global_mean_CI(ntill_res_2050_gcm, '2050')
ntill_res_2050_gl[, scenario := 'ntill-res']

# constrained | 2030
ntill_res_2030_dt  = ntill_res_cst_dt[y_block == 2030,]
ntill_res_2030_dt  = ntill_res_2030_dt[, ..k_cols]
ntill_res_2030_gcm = ntill_res_2030_dt[, lapply(.SD, sum), .SDcols = c('GHG_area','GRAIN_area','total_crop_area_ha'), by = .(ssp, gcm, y_block)]
# potential | 2030
ntill_res_2030_gl  = global_mean_CI(ntill_res_2030_gcm, '2030')
ntill_res_2030_gl[, scenario := 'ntill-res']

# CCG
ccg_cst_dt = ccg_dt[, GHG_area   := s_GHG*total_crop_area_ha]
ccg_cst_dt = ccg_dt[, GRAIN_area := s_cr_grain*total_crop_area_ha]
# constrained | 2050
ccg_2050_dt  = ccg_cst_dt[y_block == 2050,]
ccg_2050_dt  = ccg_2050_dt[, ..k_cols]
ccg_2050_gcm = ccg_2050_dt[, lapply(.SD, sum), .SDcols = c('GHG_area','GRAIN_area','total_crop_area_ha'), by = .(ssp, gcm, y_block)]
# potential | 2050
ccg_2050_gl  = global_mean_CI(ccg_2050_gcm, '2050')
ccg_2050_gl[, scenario := 'ccg']

# constrained | 2030
ccg_2030_dt  = ccg_cst_dt[y_block == 2030,]
ccg_2030_dt  = ccg_2030_dt[, ..k_cols]
ccg_2030_gcm = ccg_2030_dt[, lapply(.SD, sum), .SDcols = c('GHG_area','GRAIN_area','total_crop_area_ha'), by = .(ssp, gcm, y_block)]
# potential | 2030
ccg_2030_gl  = global_mean_CI(ccg_2030_gcm, '2030')
ccg_2030_gl[, scenario := 'ccg']

# CCG-RES
ccg_res_cst_dt = ccg_res_dt[, GHG_area   := s_GHG*total_crop_area_ha]
ccg_res_cst_dt = ccg_res_dt[, GRAIN_area := s_cr_grain*total_crop_area_ha]
# constrained | 2050
ccg_res_2050_dt  = ccg_res_cst_dt[y_block == 2050,]
ccg_res_2050_dt  = ccg_res_2050_dt[, ..k_cols]
ccg_res_2050_gcm = ccg_res_2050_dt[, lapply(.SD, sum), .SDcols = c('GHG_area','GRAIN_area','total_crop_area_ha'), by = .(ssp, gcm, y_block)]
# potential | 2050
ccg_res_2050_gl  = global_mean_CI(ccg_res_2050_gcm, '2050')
ccg_res_2050_gl[, scenario := 'ccg-res']

# constrained | 2030
ccg_res_2030_dt  = ccg_res_cst_dt[y_block == 2030,]
ccg_res_2030_dt  = ccg_res_2030_dt[, ..k_cols]
ccg_res_2030_gcm = ccg_res_2030_dt[, lapply(.SD, sum), .SDcols = c('GHG_area','GRAIN_area','total_crop_area_ha'), by = .(ssp, gcm, y_block)]
# potential | 2030
ccg_res_2030_gl  = global_mean_CI(ccg_res_2030_gcm, '2030')
ccg_res_2030_gl[, scenario := 'ccg-res']

# CCL
ccl_cst_dt = ccl_dt[, GHG_area   := s_GHG*total_crop_area_ha]
ccl_cst_dt = ccl_dt[, GRAIN_area := s_cr_grain*total_crop_area_ha]
# constrained | 2050
ccl_2050_dt  = ccl_cst_dt[y_block == 2050,]
ccl_2050_dt  = ccl_2050_dt[, ..k_cols]
ccl_2050_gcm = ccl_2050_dt[, lapply(.SD, sum), .SDcols = c('GHG_area','GRAIN_area','total_crop_area_ha'), by = .(ssp, gcm, y_block)]
# potential | 2050
ccl_2050_gl  = global_mean_CI(ccl_2050_gcm, '2050')
ccl_2050_gl[, scenario := 'ccl']

# constrained | 2030
ccl_2030_dt  = ccl_cst_dt[y_block == 2030,]
ccl_2030_dt  = ccl_2030_dt[, ..k_cols]
ccl_2030_gcm = ccl_2030_dt[, lapply(.SD, sum), .SDcols = c('GHG_area','GRAIN_area','total_crop_area_ha'), by = .(ssp, gcm, y_block)]
# potential | 2030
ccl_2030_gl  = global_mean_CI(ccl_2030_gcm, '2030')
ccl_2030_gl[, scenario := 'ccl']

# CCL-RES
ccl_res_cst_dt = ccl_res_dt[, GHG_area   := s_GHG*total_crop_area_ha]
ccl_res_cst_dt = ccl_res_dt[, GRAIN_area := s_cr_grain*total_crop_area_ha]
# constrained | 2050
ccl_res_2050_dt  = ccl_res_cst_dt[y_block == 2050,]
ccl_res_2050_dt  = ccl_res_2050_dt[, ..k_cols]
ccl_res_2050_gcm = ccl_res_2050_dt[, lapply(.SD, sum), .SDcols = c('GHG_area','GRAIN_area','total_crop_area_ha'), by = .(ssp, gcm, y_block)]
# potential | 2050
ccl_res_2050_gl  = global_mean_CI(ccl_res_2050_gcm, '2050')
ccl_res_2050_gl[, scenario := 'ccl-res']

# constrained | 2030
ccl_res_2030_dt  = ccl_res_cst_dt[y_block == 2030,]
ccl_res_2030_dt  = ccl_res_2030_dt[, ..k_cols]
ccl_res_2030_gcm = ccl_res_2030_dt[, lapply(.SD, sum), .SDcols = c('GHG_area','GRAIN_area','total_crop_area_ha'), by = .(ssp, gcm, y_block)]
# potential | 2030
ccl_res_2030_gl  = global_mean_CI(ccl_res_2030_gcm, '2030')
ccl_res_2030_gl[, scenario := 'ccl-res']

# CCG-NTILL-RES
ccg_ntill_cst_dt = ccg_ntill_dt[, GHG_area   := s_GHG*total_crop_area_ha]
ccg_ntill_cst_dt = ccg_ntill_dt[, GRAIN_area := s_cr_grain*total_crop_area_ha]
# constrained | 2050
ccg_ntill_2050_dt  = ccg_ntill_cst_dt[y_block == 2050,]
ccg_ntill_2050_dt  = ccg_ntill_2050_dt[, ..k_cols]
ccg_ntill_2050_gcm = ccg_ntill_2050_dt[, lapply(.SD, sum), .SDcols = c('GHG_area','GRAIN_area','total_crop_area_ha'), by = .(ssp, gcm, y_block)]
# potential | 2050
ccg_ntill_2050_gl  = global_mean_CI(ccg_ntill_2050_gcm, '2050')
ccg_ntill_2050_gl[, scenario := 'ccg-ntill']

# constrained | 2030
ccg_ntill_2030_dt  = ccg_ntill_cst_dt[y_block == 2030,]
ccg_ntill_2030_dt  = ccg_ntill_2030_dt[, ..k_cols]
ccg_ntill_2030_gcm = ccg_ntill_2030_dt[, lapply(.SD, sum), .SDcols = c('GHG_area','GRAIN_area','total_crop_area_ha'), by = .(ssp, gcm, y_block)]
# potential | 2030
ccg_ntill_2030_gl  = global_mean_CI(ccg_ntill_2030_gcm, '2030')
ccg_ntill_2030_gl[, scenario := 'ccg-ntill']

# CCL-NTILL-RES
ccl_ntill_cst_dt = ccl_ntill_dt[, GHG_area   := s_GHG*total_crop_area_ha]
ccl_ntill_cst_dt = ccl_ntill_dt[, GRAIN_area := s_cr_grain*total_crop_area_ha]
# constrained | 2050
ccl_ntill_2050_dt  = ccl_ntill_cst_dt[y_block == 2050,]
ccl_ntill_2050_dt  = ccl_ntill_2050_dt[, ..k_cols]
ccl_ntill_2050_gcm = ccl_ntill_2050_dt[, lapply(.SD, sum), .SDcols = c('GHG_area','GRAIN_area','total_crop_area_ha'), by = .(ssp, gcm, y_block)]
# potential | 2050
ccl_ntill_2050_gl  = global_mean_CI(ccl_ntill_2050_gcm, '2050')
ccl_ntill_2050_gl[, scenario := 'ccl-ntill']

# constrained | 2030
ccl_ntill_2030_dt  = ccl_ntill_cst_dt[y_block == 2030,]
ccl_ntill_2030_dt  = ccl_ntill_2030_dt[, ..k_cols]
ccl_ntill_2030_gcm = ccl_ntill_2030_dt[, lapply(.SD, sum), .SDcols = c('GHG_area','GRAIN_area','total_crop_area_ha'), by = .(ssp, gcm, y_block)]
# potential | 2030
ccl_ntill_2030_gl  = global_mean_CI(ccl_ntill_2030_gcm, '2030')
ccl_ntill_2030_gl[, scenario := 'ccl-ntill']

# COMBINE & SAVE DT
# 2030
biophys_2030 = rbind(residue_2030_gl, ntill_2030_gl, ntill_res_2030_gl,
                       ccg_2030_gl, ccg_res_2030_gl, ccl_2030_gl, ccl_res_2030_gl,
                       ccg_ntill_2030_gl, ccl_ntill_2030_gl)
setcolorder(biophys_2030, c('ssp', 'y_block', 'scenario'))
fwrite(biophys_2030, file = paste(data.path, 'global-biophysical-GHG-mitigation-potential-2030.csv', sep = '/'))
# 2050
biophys_2050 = rbind(residue_2050_gl, ntill_2050_gl, ntill_res_2050_gl,
                       ccg_2050_gl, ccg_res_2050_gl, ccl_2050_gl, ccl_res_2050_gl,
                       ccg_ntill_2050_gl, ccl_ntill_2050_gl)
setcolorder(biophys_2050, c('ssp', 'y_block', 'scenario'))
fwrite(biophys_2050, file = paste(data.path, 'global-biophysical-GHG-mitigation-potential-2050.csv', sep = '/'))
#-------------------------------------------------------------------------------
# BIOPHYSICAL MITIGATION POTENTIAL | REGION
#-------------------------------------------------------------------------------
k_cols = c('cell','ssp', 'gcm','IPCC_NAME','y_block','total_crop_area_ha','GHG_area', 'GRAIN_area')
# RESIDUE
residue_cst_dt = residue_dt[, GHG_area   := s_GHG*total_crop_area_ha]
residue_cst_dt = residue_dt[, GRAIN_area := s_cr_grain*total_crop_area_ha]
# constrained | 2050
residue_2050_dt    = residue_cst_dt[y_block == 2050,]
residue_2050_dt    = residue_2050_dt[, ..k_cols]
residue_2050_r_gcm = residue_2050_dt[, lapply(.SD, sum), .SDcols = c('GHG_area','GRAIN_area','total_crop_area_ha'), by = .(ssp, gcm, y_block, IPCC_NAME)]
# potential | 2050
residue_2050_rg    = regional_mean_CI(residue_2050_r_gcm, '2050')
residue_2050_rg[, scenario := 'res']

# constrained | 2030
residue_2030_dt    = residue_cst_dt[y_block == 2030,]
residue_2030_dt    = residue_2030_dt[, ..k_cols]
residue_2030_r_gcm = residue_2030_dt[, lapply(.SD, sum), .SDcols = c('GHG_area','GRAIN_area','total_crop_area_ha'), by = .(ssp, gcm, y_block, IPCC_NAME)]
# potential | 2030
residue_2030_rg    = regional_mean_CI(residue_2030_r_gcm, '2030')
residue_2030_rg[, scenario := 'res']

# NTILL
ntill_cst_dt = ntill_dt[, GHG_area   := s_GHG*total_crop_area_ha]
ntill_cst_dt = ntill_dt[, GRAIN_area := s_cr_grain*total_crop_area_ha]
# constrained | 2050
ntill_2050_dt    = ntill_cst_dt[y_block == 2050,]
ntill_2050_dt    = ntill_2050_dt[, ..k_cols]
ntill_2050_r_gcm = ntill_2050_dt[, lapply(.SD, sum), .SDcols = c('GHG_area','GRAIN_area','total_crop_area_ha'), by = .(ssp, gcm, y_block, IPCC_NAME)]
# potential | 2050
ntill_2050_rg    = regional_mean_CI(ntill_2050_r_gcm, '2050')
ntill_2050_rg[, scenario := 'ntill']

# constrained | 2030
ntill_2030_dt    = ntill_cst_dt[y_block == 2030,]
ntill_2030_dt    = ntill_2030_dt[, ..k_cols]
ntill_2030_r_gcm = ntill_2030_dt[, lapply(.SD, sum), .SDcols = c('GHG_area','GRAIN_area','total_crop_area_ha'), by = .(ssp, gcm, y_block, IPCC_NAME)]
# potential | 2030
ntill_2030_rg    = regional_mean_CI(ntill_2030_r_gcm, '2030')
ntill_2030_rg[, scenario := 'ntill']

# NTILL-RES
ntill_res_cst_dt = ntill_res_dt[, GHG_area   := s_GHG*total_crop_area_ha]
ntill_res_cst_dt = ntill_res_dt[, GRAIN_area := s_cr_grain*total_crop_area_ha]
# constrained | 2050
ntill_res_2050_dt    = ntill_res_cst_dt[y_block == 2050,]
ntill_res_2050_dt    = ntill_res_2050_dt[, ..k_cols]
ntill_res_2050_r_gcm = ntill_res_2050_dt[, lapply(.SD, sum), .SDcols = c('GHG_area','GRAIN_area','total_crop_area_ha'), by = .(ssp, gcm, y_block, IPCC_NAME)]
# potential | 2050
ntill_res_2050_rg  = regional_mean_CI(ntill_res_2050_r_gcm, '2050')
ntill_res_2050_rg[, scenario := 'ntill-res']

# constrained | 2030
ntill_res_2030_dt    = ntill_res_cst_dt[y_block == 2030,]
ntill_res_2030_dt    = ntill_res_2030_dt[, ..k_cols]
ntill_res_2030_r_gcm = ntill_res_2030_dt[, lapply(.SD, sum), .SDcols = c('GHG_area','GRAIN_area','total_crop_area_ha'), by = .(ssp, gcm, y_block, IPCC_NAME)]
# potential | 2030
ntill_res_2030_rg    = regional_mean_CI(ntill_res_2030_r_gcm, '2030')
ntill_res_2030_rg[, scenario := 'ntill-res']

# CCG
ccg_cst_dt = ccg_dt[, GHG_area   := s_GHG*total_crop_area_ha]
ccg_cst_dt = ccg_dt[, GRAIN_area := s_cr_grain*total_crop_area_ha]
# constrained | 2050
ccg_2050_dt    = ccg_cst_dt[y_block == 2050,]
ccg_2050_dt    = ccg_2050_dt[, ..k_cols]
ccg_2050_r_gcm = ccg_2050_dt[, lapply(.SD, sum), .SDcols = c('GHG_area','GRAIN_area','total_crop_area_ha'), by = .(ssp, gcm, y_block, IPCC_NAME)]
# potential | 2050
ccg_2050_rg    = regional_mean_CI(ccg_2050_r_gcm, '2050')
ccg_2050_rg[, scenario := 'ccg']

# constrained | 2030
ccg_2030_dt    = ccg_cst_dt[y_block == 2030,]
ccg_2030_dt    = ccg_2030_dt[, ..k_cols]
ccg_2030_r_gcm = ccg_2030_dt[, lapply(.SD, sum), .SDcols = c('GHG_area','GRAIN_area','total_crop_area_ha'), by = .(ssp, gcm, y_block, IPCC_NAME)]
# potential | 2030
ccg_2030_rg    = regional_mean_CI(ccg_2030_r_gcm, '2030')
ccg_2030_rg[, scenario := 'ccg']

# CCG-RES
ccg_res_cst_dt = ccg_res_dt[, GHG_area   := s_GHG*total_crop_area_ha]
ccg_res_cst_dt = ccg_res_dt[, GRAIN_area := s_cr_grain*total_crop_area_ha]
# constrained | 2050
ccg_res_2050_dt    = ccg_res_cst_dt[y_block == 2050,]
ccg_res_2050_dt    = ccg_res_2050_dt[, ..k_cols]
ccg_res_2050_r_gcm = ccg_res_2050_dt[, lapply(.SD, sum), .SDcols = c('GHG_area','GRAIN_area','total_crop_area_ha'), by = .(ssp, gcm, y_block, IPCC_NAME)]
# potential | 2050
ccg_res_2050_rg  = regional_mean_CI(ccg_res_2050_r_gcm, '2050')
ccg_res_2050_rg[, scenario := 'ccg-res']

# constrained | 2030
ccg_res_2030_dt    = ccg_res_cst_dt[y_block == 2030,]
ccg_res_2030_dt    = ccg_res_2030_dt[, ..k_cols]
ccg_res_2030_r_gcm = ccg_res_2030_dt[, lapply(.SD, sum), .SDcols = c('GHG_area','GRAIN_area','total_crop_area_ha'), by = .(ssp, gcm, y_block, IPCC_NAME)]
# potential | 2030
ccg_res_2030_rg    = regional_mean_CI(ccg_res_2030_r_gcm, '2030')
ccg_res_2030_rg[, scenario := 'ccg-res']

# CCL
ccl_cst_dt = ccl_dt[, GHG_area   := s_GHG*total_crop_area_ha]
ccl_cst_dt = ccl_dt[, GRAIN_area := s_cr_grain*total_crop_area_ha]
# constrained | 2050
ccl_2050_dt    = ccl_cst_dt[y_block == 2050,]
ccl_2050_dt    = ccl_2050_dt[, ..k_cols]
ccl_2050_r_gcm = ccl_2050_dt[, lapply(.SD, sum), .SDcols = c('GHG_area','GRAIN_area','total_crop_area_ha'), by = .(ssp, gcm, y_block, IPCC_NAME)]
# potential | 2050
ccl_2050_rg    = regional_mean_CI(ccl_2050_r_gcm, '2050')
ccl_2050_rg[, scenario := 'ccl']

# constrained | 2030
ccl_2030_dt    = ccl_cst_dt[y_block == 2030,]
ccl_2030_dt    = ccl_2030_dt[, ..k_cols]
ccl_2030_r_gcm = ccl_2030_dt[, lapply(.SD, sum), .SDcols = c('GHG_area','GRAIN_area','total_crop_area_ha'), by = .(ssp, gcm, y_block, IPCC_NAME)]
# potential | 2030
ccl_2030_rg    = regional_mean_CI(ccl_2030_r_gcm, '2030')
ccl_2030_rg[, scenario := 'ccl']

# CCL-RES
ccl_res_cst_dt = ccl_res_dt[, GHG_area   := s_GHG*total_crop_area_ha]
ccl_res_cst_dt = ccl_res_dt[, GRAIN_area := s_cr_grain*total_crop_area_ha]
# constrained | 2050
ccl_res_2050_dt    = ccl_res_cst_dt[y_block == 2050,]
ccl_res_2050_dt    = ccl_res_2050_dt[, ..k_cols]
ccl_res_2050_r_gcm = ccl_res_2050_dt[, lapply(.SD, sum), .SDcols = c('GHG_area','GRAIN_area','total_crop_area_ha'), by = .(ssp, gcm, y_block, IPCC_NAME)]
# potential | 2050
ccl_res_2050_rg  = regional_mean_CI(ccl_res_2050_r_gcm, '2050')
ccl_res_2050_rg[, scenario := 'ccl-res']

# constrained | 2030
ccl_res_2030_dt    = ccl_res_cst_dt[y_block == 2030,]
ccl_res_2030_dt    = ccl_res_2030_dt[, ..k_cols]
ccl_res_2030_r_gcm = ccl_res_2030_dt[, lapply(.SD, sum), .SDcols = c('GHG_area','GRAIN_area','total_crop_area_ha'), by = .(ssp, gcm, y_block, IPCC_NAME)]
# potential | 2030
ccl_res_2030_rg    = regional_mean_CI(ccl_res_2030_r_gcm, '2030')
ccl_res_2030_rg[, scenario := 'ccl-res']

# CCG-NTILL-RES
ccg_ntill_cst_dt = ccg_ntill_dt[, GHG_area   := s_GHG*total_crop_area_ha]
ccg_ntill_cst_dt = ccg_ntill_dt[, GRAIN_area := s_cr_grain*total_crop_area_ha]
# constrained | 2050
ccg_ntill_2050_dt    = ccg_ntill_cst_dt[y_block == 2050,]
ccg_ntill_2050_dt    = ccg_ntill_2050_dt[, ..k_cols]
ccg_ntill_2050_r_gcm = ccg_ntill_2050_dt[, lapply(.SD, sum), .SDcols = c('GHG_area','GRAIN_area','total_crop_area_ha'), by = .(ssp, gcm, y_block, IPCC_NAME)]
# potential | 2050
ccg_ntill_2050_rg  = regional_mean_CI(ccg_ntill_2050_r_gcm, '2050')
ccg_ntill_2050_rg[, scenario := 'ccg-ntill']

# constrained | 2030
ccg_ntill_2030_dt    = ccg_ntill_cst_dt[y_block == 2030,]
ccg_ntill_2030_dt    = ccg_ntill_2030_dt[, ..k_cols]
ccg_ntill_2030_r_gcm = ccg_ntill_2030_dt[, lapply(.SD, sum), .SDcols = c('GHG_area','GRAIN_area','total_crop_area_ha'), by = .(ssp, gcm, y_block, IPCC_NAME)]
# potential | 2030
ccg_ntill_2030_rg    = regional_mean_CI(ccg_ntill_2030_r_gcm, '2030')
ccg_ntill_2030_rg[, scenario := 'ccg-ntill']

# CCL-NTILL-RES
ccl_ntill_cst_dt = ccl_ntill_dt[, GHG_area   := s_GHG*total_crop_area_ha]
ccl_ntill_cst_dt = ccl_ntill_dt[, GRAIN_area := s_cr_grain*total_crop_area_ha]
# constrained | 2050
ccl_ntill_2050_dt    = ccl_ntill_cst_dt[y_block == 2050,]
ccl_ntill_2050_dt    = ccl_ntill_2050_dt[, ..k_cols]
ccl_ntill_2050_r_gcm = ccl_ntill_2050_dt[, lapply(.SD, sum), .SDcols = c('GHG_area','GRAIN_area','total_crop_area_ha'), by = .(ssp, gcm, y_block, IPCC_NAME)]
# potential | 2050
ccl_ntill_2050_rg  = regional_mean_CI(ccl_ntill_2050_r_gcm, '2050')
ccl_ntill_2050_rg[, scenario := 'ccl-ntill']

# constrained | 2030
ccl_ntill_2030_dt    = ccl_ntill_cst_dt[y_block == 2030,]
ccl_ntill_2030_dt    = ccl_ntill_2030_dt[, ..k_cols]
ccl_ntill_2030_r_gcm = ccl_ntill_2030_dt[, lapply(.SD, sum), .SDcols = c('GHG_area','GRAIN_area','total_crop_area_ha'), by = .(ssp, gcm, y_block, IPCC_NAME)]
# potential | 2030
ccl_ntill_2030_rg    = regional_mean_CI(ccl_ntill_2030_r_gcm, '2030')
ccl_ntill_2030_rg[, scenario := 'ccl-ntill']

# COMBINE & SAVE DT
# 2030
ipcc_biophys_2030 = rbind(residue_2030_rg, ntill_2030_rg, ntill_res_2030_rg,
                            ccg_2030_rg, ccg_res_2030_rg, ccl_2030_rg, ccl_res_2030_rg,
                            ccg_ntill_2030_rg, ccl_ntill_2030_rg)
setcolorder(ipcc_biophys_2030, c('IPCC_NAME','ssp', 'y_block', 'scenario'))
fwrite(ipcc_biophys_2030, file = paste(data.path, 'regional-biophysical-GHG-mitigation-potential-2030.csv', sep = '/'))
# 2050
ipcc_biophys_2050 = rbind(residue_2050_rg, ntill_2050_rg, ntill_res_2050_rg,
                            ccg_2050_rg, ccg_res_2050_rg, ccl_2050_rg, ccl_res_2050_rg,
                            ccg_ntill_2050_rg, ccl_ntill_2050_rg)
setcolorder(ipcc_biophys_2050, c('IPCC_NAME','ssp', 'y_block', 'scenario'))
fwrite(ipcc_biophys_2050, file = paste(data.path, 'regional-biophysical-GHG-mitigation-potential-2050.csv', sep = '/'))
