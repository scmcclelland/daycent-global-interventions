# filename:    analysis_bestpractice.R
# created:     28 February 2024
# updated:     28 February 2024
# author:      S.C. McClelland
# description: This file analyzes GHG mitigation potential from DayCent simulations
#              for global cropland over time for different interventions and SSPs. 
#              Best practice is either (1) max GHG mitigation potential 
#              or (2) max GHG mitigation potential when constrained by yield.
# definition:  Yield constrained mitigation potential is defined as cropland area where 
#              'the difference in ensemble mean cumulative crop-area weighted grain yield 
#              is >= 0 in 2030 or 2050.' The cumulative GHG mitigation potential is then
#              divided by the number of years since the start of the simulation (15 or 35 years).
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
residue_dt = readRDS(paste(data.path, 'imputed-ensemble-relative-responses-weighted-mean-res.rds', sep = '/'))
residue_dt = dcast(residue_dt,
                   cell + x + y + y_block + ssp + total_crop_area_ha ~ variable,
                   value.var = 'value')
residue_dt[, scenario := 'res']
# NTILL
ntill_dt = readRDS(paste(data.path, 'imputed-ensemble-relative-responses-weighted-mean-ntill.rds', sep = '/'))
ntill_dt = dcast(ntill_dt,
                   cell + x + y + y_block + ssp + total_crop_area_ha ~ variable,
                   value.var = 'value')
ntill_dt[, scenario := 'ntill']
# NTILL-RES
ntill_res_dt = readRDS(paste(data.path, 'imputed-ensemble-relative-responses-weighted-mean-ntill-res.rds', sep = '/'))
ntill_res_dt = dcast(ntill_res_dt,
                 cell + x + y + y_block + ssp + total_crop_area_ha ~ variable,
                 value.var = 'value')
ntill_res_dt[, scenario := 'ntill-res']
# CCG
ccg_dt = readRDS(paste(data.path, 'imputed-ensemble-relative-responses-weighted-mean-ccg.rds', sep = '/'))
ccg_dt = dcast(ccg_dt,
                 cell + x + y + y_block + ssp + total_crop_area_ha ~ variable,
                 value.var = 'value')
ccg_dt[, scenario := 'ccg']
# CCG-RES
ccg_res_dt = readRDS(paste(data.path, 'imputed-ensemble-relative-responses-weighted-mean-ccg-res.rds', sep = '/'))
ccg_res_dt = dcast(ccg_res_dt,
               cell + x + y + y_block + ssp + total_crop_area_ha ~ variable,
               value.var = 'value')
ccg_res_dt[, scenario := 'ccg-res']
# CCL
ccl_dt = readRDS(paste(data.path, 'imputed-ensemble-relative-responses-weighted-mean-ccl.rds', sep = '/'))
ccl_dt = dcast(ccl_dt,
               cell + x + y + y_block + ssp + total_crop_area_ha ~ variable,
               value.var = 'value')
ccl_dt[, scenario := 'ccl']
# CCL-RES
ccl_res_dt = readRDS(paste(data.path, 'imputed-ensemble-relative-responses-weighted-mean-ccl-res.rds', sep = '/'))
ccl_res_dt = dcast(ccl_res_dt,
                   cell + x + y + y_block + ssp + total_crop_area_ha ~ variable,
                   value.var = 'value')
ccl_res_dt[, scenario := 'ccl-res']
# CCG-NTILL-RES
ccg_ntill_dt = readRDS(paste(data.path, 'imputed-ensemble-relative-responses-weighted-mean-ccg-ntill.rds', sep = '/'))
ccg_ntill_dt = dcast(ccg_ntill_dt,
                   cell + x + y + y_block + ssp + total_crop_area_ha ~ variable,
                   value.var = 'value')
ccg_ntill_dt[, scenario := 'ccg-ntill']
# CCL-NTILL-RES
ccl_ntill_dt = readRDS(paste(data.path, 'imputed-ensemble-relative-responses-weighted-mean-ccl-ntill.rds', sep = '/'))
ccl_ntill_dt = dcast(ccl_ntill_dt,
                     cell + x + y + y_block + ssp + total_crop_area_ha ~ variable,
                     value.var = 'value')
ccl_ntill_dt[, scenario := 'ccl-ntill']
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
# BEST MANAGEMENT PRACTICE THRU 2050 | NO YIELD CONSTRAINT, AREA-WEIGHTED
#-------------------------------------------------------------------------------
# combine dt
ensemble_dt      = rbind(residue_dt, ntill_dt, ntill_res_dt, ccg_dt, ccg_res_dt,
                    ccl_dt, ccl_res_dt, ccg_ntill_dt, ccl_ntill_dt)
# multiply s_GHG by area | N.B. may not want to embed area into calculation
ensemble_dt[, GHG_area   := s_GHG*total_crop_area_ha]
ensemble_dt[, GRAIN_area := s_cr_grain*total_crop_area_ha]

ensemble_2050_dt = ensemble_dt[y_block == 2050,]
setorder(ensemble_2050_dt, cell, scenario)
setcolorder(ensemble_2050_dt, c('cell', 'x', 'y', 'y_block', 'ssp', 'scenario'))

cols             = c('cell', 'x', 'y', 'ssp', 'scenario', 'GHG_area')
ensemble_2050_dt = ensemble_2050_dt[, ..cols]
# melt
ensemble_2050_dt = melt(ensemble_2050_dt,
             id.vars = c("cell", "x", "y", "ssp", "scenario"),
             measure.vars = "GHG_area")
# get min value
bmp_ensemble_2050_dt = ensemble_2050_dt[, .SD[which.min(value)], by = .(cell, x, y, ssp)]
# dcast
bmp_ensemble_2050_dt = dcast(bmp_ensemble_2050_dt,
                             cell + x + y + ssp + scenario ~ variable,
                             value.var = 'value')
# join with original dt
ensemble_gr_2050_dt  = ensemble_dt[y_block == 2050, c('cell', 'x', 'y', 'ssp',
                                                      'scenario', 'GHG_area',
                                                      'GRAIN_area')]
bmp_ensemble_2050_dt = bmp_ensemble_2050_dt[ensemble_gr_2050_dt, on = .(cell = cell,
                                                                        x = x,
                                                                        y = y,
                                                                        ssp = ssp,
                                                                        scenario = scenario)]
bmp_ensemble_2050_dt = bmp_ensemble_2050_dt[!is.na(GHG_area)]
setorder(bmp_ensemble_2050_dt, cell)
NROW(bmp_ensemble_2050_dt) # 123108
# recode positive GHG, GRAIN cells as 'no change'
bmp_ensemble_2050_dt[GHG_area > 0, scenario := 'BAU']
bmp_ensemble_2050_dt[scenario %in% 'BAU', GHG_area := 0]
bmp_ensemble_2050_dt[scenario %in% 'BAU', GRAIN_area := 0]
bmp_ensemble_2050_dt[, i.GHG_area := NULL]
# check
sum(bmp_ensemble_2050_dt[ssp %in% 'historical', GHG_area])
sum(bmp_ensemble_2050_dt[ssp %in% 'historical', GRAIN_area])
# save 
fwrite(bmp_ensemble_2050_dt, file = paste(data.path, 'best-management-practice-GHG-mitigation-2050.csv', sep = '/'))
#-------------------------------------------------------------------------------
# BEST MANAGEMENT PRACTICE THRU 2050 | YIELD CONSTRAINT
#-------------------------------------------------------------------------------
ensemble_dt         = rbind(residue_dt, ntill_dt, ntill_res_dt, ccg_dt, ccg_res_dt,
                         ccl_dt, ccl_res_dt, ccg_ntill_dt, ccl_ntill_dt)
# multiply s_GHG by area | N.B. may not want to embed area into calculation
ensemble_dt[, GHG_area   := s_GHG*total_crop_area_ha]
ensemble_dt[, GRAIN_area := s_cr_grain*total_crop_area_ha]

ensemble_2050_yc_dt = ensemble_dt[y_block == 2050,]
setorder(ensemble_2050_yc_dt, cell, scenario)
setcolorder(ensemble_2050_yc_dt, c('cell', 'x', 'y', 'y_block', 'ssp', 'scenario'))

cols                = c('cell', 'x', 'y', 'ssp', 'scenario', 'GHG_area', 'GRAIN_area')
ensemble_2050_yc_dt = ensemble_2050_yc_dt[, ..cols]
# filter yield constrained
ensemble_2050_yc_dt = ensemble_2050_yc_dt[GRAIN_area < 0,  YC := 'yes']
ensemble_2050_yc_dt = ensemble_2050_yc_dt[GRAIN_area >= 0, YC := 'no']
# melt
ensemble_2050_yc_dt = melt(ensemble_2050_yc_dt,
                        id.vars = c("cell", "x", "y", "ssp", "scenario", "YC", "GRAIN_area"),
                        measure.vars = "GHG_area")
# get min value for cell + scenario by yield constraint
bmp_ensemble_2050_yc_dt = ensemble_2050_yc_dt[, .SD[which.min(value)], by = .(cell, x, y, ssp, YC)]
# recode positive GHG cells as 'no change'
bmp_ensemble_2050_yc_dt = bmp_ensemble_2050_yc_dt[value > 0, scenario := 'BAU']
# update yield constraint cells
bmp_ensemble_2050_yc_dt = bmp_ensemble_2050_yc_dt[YC %in% 'yes', value := 0]
# refilter for scenarios
bmp_ensemble_2050_yc_dt = bmp_ensemble_2050_yc_dt[, .SD[which.min(value)], by = .(cell, x, y, ssp)]
# update YC selected scenarios to BAU scenario
bmp_ensemble_2050_yc_dt = bmp_ensemble_2050_yc_dt[YC %in% 'yes', scenario := 'BAU']
bmp_ensemble_2050_yc_dt[scenario %in% 'BAU', value      := 0]
bmp_ensemble_2050_yc_dt[scenario %in% 'BAU', GRAIN_area := 0]
bmp_ensemble_2050_yc_dt[, YC := NULL]

# dcast
bmp_ensemble_2050_yc_dt = dcast(bmp_ensemble_2050_yc_dt,
                             cell + x + y + ssp + scenario + GRAIN_area ~ variable,
                             value.var = 'value')
NROW(bmp_ensemble_2050_yc_dt) # 123108
# check
sum(bmp_ensemble_2050_yc_dt[ssp %in% 'historical', GHG_area])
sum(bmp_ensemble_2050_yc_dt[ssp %in% 'historical', GRAIN_area])
# save 
fwrite(bmp_ensemble_2050_yc_dt, file = paste(data.path, 'best-management-practice-GHG-mitigation-yield-constrained-2050.csv', sep = '/'))
