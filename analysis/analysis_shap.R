# filename:    analysis_shap.R
# created:     28 February 2024
# updated:     15 August 2024
# author:      S.C. McClelland
# description: This file contains the supervised clustering analysis
#-------------------------------------------------------------------------------
# NOTE: not all variables described here | to be added at future date

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
library(caret)
library(data.table)
library(factoextra)
library(fastshap)
library(ranger)
library(rstudioapi)
library(scales)
library(sf)
library(shapviz)
library(terra)
library(umap)
#-------------------------------------------------------------------------------
source('results_functions.R')
#-------------------------------------------------------------------------------
options(scipen = 999, digits = 4)
options(rgl.printRglwidget = TRUE)
base.path = dirname(getActiveDocumentContext()$path)
data.path = paste(base.path, 'manuscript-data', sep = '/')
fig.path  = paste(base.path, 'manuscript-figures', sep = '/')
#-------------------------------------------------------------------------------
# LOAD ABSOLUTE DATA (actual responses, not relative to REF)
#-------------------------------------------------------------------------------
# NTILL-RES
ntill_res_flux_dt = readRDS(paste(data.path, 'ensemble-absolute-responses-weighted-mean-ntill-res.rds', sep = '/'))
# CCG
ccg_flux_dt       = readRDS(paste(data.path, 'ensemble-absolute-responses-weighted-mean-ccg.rds', sep = '/'))
# CCL
ccl_flux_dt       = readRDS(paste(data.path, 'ensemble-absolute-responses-weighted-mean-ccl.rds', sep = '/'))
# CCG-NTILL-RES
ccg_ntill_flux_dt = readRDS(paste(data.path, 'ensemble-absolute-responses-weighted-mean-ccg-ntill.rds', sep = '/'))
# CCL-NTILL-RES
ccl_ntill_flux_dt = readRDS(paste(data.path, 'ensemble-absolute-responses-weighted-mean-ccl-ntill.rds', sep = '/'))
#-------------------------------------------------------------------------------
# LOAD RELATIVE DATA (response differences, relative to REF)
#-------------------------------------------------------------------------------
# RESIDUE
rel_res_flux_dt       = readRDS(paste(data.path, 'ensemble-relative-responses-weighted-mean-res.rds', sep = '/'))
# NO-TILLAGE
rel_ntill_flux_dt     = readRDS(paste(data.path, 'ensemble-relative-responses-weighted-mean-ntill.rds', sep = '/'))
# NO-TILLAGE, RESIDUE
rel_ntill_res_flux_dt = readRDS(paste(data.path, 'ensemble-relative-responses-weighted-mean-ntill-res.rds', sep = '/'))
# GRASS COVER CROP
rel_ccg_flux_dt       = readRDS(paste(data.path, 'ensemble-relative-responses-weighted-mean-ccg.rds', sep = '/'))
# GRASS COVER CROP, RESIDUE
rel_ccg_res_flux_dt   = readRDS(paste(data.path, 'ensemble-relative-responses-weighted-mean-ccg-res.rds', sep = '/'))
# LEGUME COVER CROP
rel_ccl_flux_dt       = readRDS(paste(data.path, 'ensemble-relative-responses-weighted-mean-ccl.rds', sep = '/'))
# LEGUME COVER CROP, RESIDUE
rel_ccl_res_flux_dt   = readRDS(paste(data.path, 'ensemble-relative-responses-weighted-mean-ccl-res.rds', sep = '/'))
# GRASS COVER CROP, NO-TILLAGE, RESIDUE
rel_ccg_ntill_flux_dt = readRDS(paste(data.path, 'ensemble-relative-responses-weighted-mean-ccg-ntill.rds', sep = '/'))
# LEGUME COVER CROP, NO-TILLAGE, RESIDUE
rel_ccl_ntill_flux_dt = readRDS(paste(data.path, 'ensemble-relative-responses-weighted-mean-ccl-ntill.rds', sep = '/'))
#-------------------------------------------------------------------------------
# SOIL & CLIMATE DATA
#------------------------------------------------------------------------------
site_wth          = readRDS(paste(data.path, 'ensemble_site_climate_data_decadal.rds', sep = '/'))
# CONSTRAIN TO 2050
site_wth          = site_wth[y_block <= 2050,]
# REDUCE VARIABLES
s_w_cols          = c('gridid', 'y_block', 'ssp', 'bio1', 'hist_bio1', 'bio12', 'hist_bio12',
                      'ELEV', 'SOMC_sum_', 'SLBLKD', 'SLCLAY', 'SLSAND', 'SLPH')
site_wth          = site_wth[, ..s_w_cols]
# ENSEMBLE MEAN OVER TIME
site_wth          = site_wth[, lapply(.SD, mean), .SDcols = s_w_cols[4:13], by = .(gridid, ssp)]
#-------------------------------------------------------------------------------
# MANAGEMENT DATA
#------------------------------------------------------------------------------
load(paste(data.path, 'input_table_by_gridid_crop_irr.RData', sep = '/'))
input_cols        = c('gridid', 'crop', 'irr','fertN.amt', 'orgN.amt', 'res.rtrn.amt','frac_NH4', 'frac_NO3', 'frac_Urea')
main_table        = main_table[, ..input_cols]
main_table[, N.amt := fertN.amt + orgN.amt]
in_mn_cols        = c('N.amt','fertN.amt', 'orgN.amt', 'res.rtrn.amt','frac_NH4', 'frac_NO3', 'frac_Urea')
# gridid mean
main_table        = main_table[, lapply(.SD, mean), .SDcols = in_mn_cols, by = .(gridid)]
#-------------------------------------------------------------------------------
# CONVERT Biomass Units: g C m-2 to kg ha-1 yr-1 or Mg ha-1
#-------------------------------------------------------------------------------
Mg_ha   = 100L
C_gr    = 0.42 # Ma et al. 2018 | value for crop 'reproductive organs'
C_res   = 0.49 # Phyllis2 database | value for all crops
C_shoot = 0.43 # Ma et al. 2018 | value for crop 'stem'
C_root  = 0.38 # Ma et al. 2018 | value for crop 'root'

# ABSOLUTE
# NTILL-RES
ntill_res_flux_dt[, s_cr_grain  := (s_cr_grain/Mg_ha)/C_gr]
ntill_res_flux_dt[, s_cr_residC := (s_cr_residC/Mg_ha)/C_res]
# CCG
ccg_flux_dt[, s_cr_grain  := (s_cr_grain/Mg_ha)/C_gr]
ccg_flux_dt[, s_cr_residC := (s_cr_residC/Mg_ha)/C_res]
ccg_flux_dt[, s_cc_shootC := (s_cc_shootC/Mg_ha)/C_shoot]
ccg_flux_dt[, s_cc_rootC  := (s_cc_rootC/Mg_ha)/C_root]
# CCL
ccl_flux_dt[, s_cr_grain  := (s_cr_grain/Mg_ha)/C_gr]
ccl_flux_dt[, s_cr_residC := (s_cr_residC/Mg_ha)/C_res]
ccl_flux_dt[, s_cc_shootC := (s_cc_shootC/Mg_ha)/C_shoot]
ccl_flux_dt[, s_cc_rootC  := (s_cc_rootC/Mg_ha)/C_root]
# CCG-NTILL-RES
ccg_ntill_flux_dt[, s_cr_grain  := (s_cr_grain/Mg_ha)/C_gr]
ccg_ntill_flux_dt[, s_cr_residC := (s_cr_residC/Mg_ha)/C_res]
ccg_ntill_flux_dt[, s_cc_shootC := (s_cc_shootC/Mg_ha)/C_shoot]
ccg_ntill_flux_dt[, s_cc_rootC  := (s_cc_rootC/Mg_ha)/C_root]
# CCL-NTILL-RES
ccl_ntill_flux_dt[, s_cr_grain  := (s_cr_grain/Mg_ha)/C_gr]
ccl_ntill_flux_dt[, s_cr_residC := (s_cr_residC/Mg_ha)/C_res]
ccl_ntill_flux_dt[, s_cc_shootC := (s_cc_shootC/Mg_ha)/C_shoot]
ccl_ntill_flux_dt[, s_cc_rootC  := (s_cc_rootC/Mg_ha)/C_root]

# RELATIVE
# RES
rel_res_flux_dt[, s_cr_grain  := (s_cr_grain/Mg_ha)/C_gr]
rel_res_flux_dt[, s_cr_residC := (s_cr_residC/Mg_ha)/C_res]
# NTILL
rel_ntill_flux_dt[, s_cr_grain  := (s_cr_grain/Mg_ha)/C_gr]
rel_ntill_flux_dt[, s_cr_residC := (s_cr_residC/Mg_ha)/C_res]
# NTILL-RES
rel_ntill_res_flux_dt[, s_cr_grain  := (s_cr_grain/Mg_ha)/C_gr]
rel_ntill_res_flux_dt[, s_cr_residC := (s_cr_residC/Mg_ha)/C_res]
# CCG
rel_ccg_flux_dt[, s_cr_grain  := (s_cr_grain/Mg_ha)/C_gr]
rel_ccg_flux_dt[, s_cr_residC := (s_cr_residC/Mg_ha)/C_res]
rel_ccg_flux_dt[, s_cc_shootC := (s_cc_shootC/Mg_ha)/C_shoot]
rel_ccg_flux_dt[, s_cc_rootC  := (s_cc_rootC/Mg_ha)/C_root]
# CCG-RES
rel_ccg_res_flux_dt[, s_cr_grain  := (s_cr_grain/Mg_ha)/C_gr]
rel_ccg_res_flux_dt[, s_cr_residC := (s_cr_residC/Mg_ha)/C_res]
rel_ccg_res_flux_dt[, s_cc_shootC := (s_cc_shootC/Mg_ha)/C_shoot]
rel_ccg_res_flux_dt[, s_cc_rootC  := (s_cc_rootC/Mg_ha)/C_root]
# CCL
rel_ccl_flux_dt[, s_cr_grain  := (s_cr_grain/Mg_ha)/C_gr]
rel_ccl_flux_dt[, s_cr_residC := (s_cr_residC/Mg_ha)/C_res]
rel_ccl_flux_dt[, s_cc_shootC := (s_cc_shootC/Mg_ha)/C_shoot]
rel_ccl_flux_dt[, s_cc_rootC  := (s_cc_rootC/Mg_ha)/C_root]
# CCL-RES
rel_ccl_res_flux_dt[, s_cr_grain  := (s_cr_grain/Mg_ha)/C_gr]
rel_ccl_res_flux_dt[, s_cr_residC := (s_cr_residC/Mg_ha)/C_res]
rel_ccl_res_flux_dt[, s_cc_shootC := (s_cc_shootC/Mg_ha)/C_shoot]
rel_ccl_res_flux_dt[, s_cc_rootC  := (s_cc_rootC/Mg_ha)/C_root]
# CCG-NTILL-RES
rel_ccg_ntill_flux_dt[, s_cr_grain  := (s_cr_grain/Mg_ha)/C_gr]
rel_ccg_ntill_flux_dt[, s_cr_residC := (s_cr_residC/Mg_ha)/C_res]
rel_ccg_ntill_flux_dt[, s_cc_shootC := (s_cc_shootC/Mg_ha)/C_shoot]
rel_ccg_ntill_flux_dt[, s_cc_rootC  := (s_cc_rootC/Mg_ha)/C_root]
# CCL-NTILL-RES
rel_ccl_ntill_flux_dt[, s_cr_grain  := (s_cr_grain/Mg_ha)/C_gr]
rel_ccl_ntill_flux_dt[, s_cr_residC := (s_cr_residC/Mg_ha)/C_res]
rel_ccl_ntill_flux_dt[, s_cc_shootC := (s_cc_shootC/Mg_ha)/C_shoot]
rel_ccl_ntill_flux_dt[, s_cc_rootC  := (s_cc_rootC/Mg_ha)/C_root]
#-------------------------------------------------------------------------------
# FACTORIAL CONTRIBUTIONS TO COMBINED RESPONSE | Legume Cover Crop
#-------------------------------------------------------------------------------
# INTERPRETATION NOTE: The random forest response (s_SOC) is the difference between
#                      the combined practice scenario, cover crop + no-tillage +
#                      full residue retention, and the 'additive model' which is 
#                      no-tillage + residue added to cover crop under BAU residue management.
# Additive model: ntill_res, ccl
#-------------------------------------------------------------------------------
# TIDY CCL-NTILL RESPONSE DT
  # filter to 2050
ccl_ntill_flux_dt = ccl_ntill_flux_dt[y_block <= 2050,]
  # mean response for nfix, soil C:N, sfdcmp, slcmp
ccl_ntill_flux_dt[, m_nfix       := mean(m_nfix), by = .(gridid, ssp)]
ccl_ntill_flux_dt[, m_sfdcmp     := mean(m_sfdcmp), by = .(gridid, ssp)]
ccl_ntill_flux_dt[, m_sldcmp     := mean(m_sldcmp), by = .(gridid, ssp)]
# sum cover crop biomass
ccl_ntill_flux_dt[, s_cc_biomass := s_cc_shootC + s_cc_rootC]
  # reduce columns
drop_cols = c('m_cr_rootC','m_cr_shootC','m_cr_shootN','m_cr_grain','m_cr_grainN','m_cr_NPP',    
  'm_cr_residC','m_cr_residN','s_cr_rootC','s_cr_shootC','s_cr_shootN', 
  's_cr_grainN','s_cr_residN','m_cc_rootC','m_cc_shootC','m_cc_shootN',
  's_cc_shootN','m_SOC','m_N2O','m_iN2O','m_dN2O','m_GHG', 'm_annet','m_sC.N','m_cr_irr',    
  's_cr_irr','s_N2O','s_iN2O','s_dN2O','s_GHG','s_cc_shootC','s_cc_rootC', 's_cr_NPP')
ccl_ntill_flux_dt = ccl_ntill_flux_dt[y_block == 2050, -..drop_cols]

# EXPECTED ADDITIVE RESPONSE DT
  # update column names
names = c('s_SOC', 's_cr_grain', 's_cr_residC','s_cc_rootC', 's_cc_shootC', 'm_sfdcmp',
          'm_sldcmp','m_nfix', 's_annet', 's_gr_nit', 's_N_leach')
new_n = as.vector(outer('ccl_', colnames(ccl_flux_dt[,..names]), paste0))
old_n = colnames(ccl_flux_dt[, ..names])
setnames(ccl_flux_dt, old = c(old_n), new = c(new_n))
cols  = c('gridid', 'x', 'y', 'y_block', 'gcm', 'ssp', 'WB_NAME', 'WB_REGION', new_n)
ccl_flux_dt = ccl_flux_dt[, ..cols]
  # reduce cols of ntill_res
cols  = c('gridid', 'x', 'y', 'y_block', 'gcm', 'ssp', 'WB_NAME', 'WB_REGION', old_n)
ntill_res_flux_dt = ntill_res_flux_dt[, ..cols]
  # join ntill_res and cover crop table
additive_dt = ntill_res_flux_dt[ccl_flux_dt, on = .(gridid = gridid,
                                                  x = x,
                                                  y = y,
                                                  y_block = y_block,
                                                  ssp = ssp,
                                                  gcm = gcm,
                                                  WB_NAME = WB_NAME,
                                                  WB_REGION = WB_REGION)]
  # filter to 2050
additive_dt = additive_dt[y_block <= 2050,]
  # mean response for nfix, soil C:N, sfdcmp, slcmp
additive_dt[, m_nfix       := mean(m_nfix),     by = .(gridid, ssp, gcm)]
additive_dt[, ccl_m_nfix   := mean(ccl_m_nfix), by = .(gridid, ssp, gcm)]
additive_dt[, m_sfdcmp     := mean(m_sfdcmp),   by = .(gridid, ssp, gcm)]
additive_dt[, ccl_m_sfdcmp := mean(ccl_m_sfdcmp), by = .(gridid, ssp, gcm)]
additive_dt[, m_sldcmp     := mean(m_sldcmp),     by = .(gridid, ssp, gcm)]
additive_dt[, ccl_m_sldcmp := mean(ccl_m_sldcmp), by = .(gridid, ssp, gcm)]
  # sum cover crop biomass
additive_dt[, s_cc_biomass     := s_cc_shootC + s_cc_rootC]
additive_dt[, ccl_s_cc_biomass := ccl_s_cc_shootC + ccl_s_cc_rootC]
  # reduce columns
drop_cols = c('s_cc_rootC', 's_cc_shootC', 'ccl_s_cc_rootC','ccl_s_cc_shootC')
additive_dt[, scenario := 'additive-model']
additive_dt = additive_dt[, -..drop_cols]
  # additive dt, mean of features | create by hand
additive_dt[, s_SOC        := s_SOC + ccl_s_SOC] # additive response
additive_dt[, s_cr_grain   := (s_cr_grain + ccl_s_cr_grain)/2]
additive_dt[, s_cr_grain   := (s_cr_grain + ccl_s_cr_grain)/2] 
additive_dt[, s_cr_residC  := (s_cr_residC + ccl_s_cr_residC)/2]
additive_dt[, s_cc_biomass := ccl_s_cc_biomass] # only cover crop
additive_dt[, m_sfdcmp     := (m_sfdcmp + ccl_m_sfdcmp)/2] 
additive_dt[, m_sldcmp     := (m_sldcmp + ccl_m_sldcmp)/2] 
additive_dt[, m_nfix       := ccl_m_nfix] # only cover crop
additive_dt[, s_annet      := (s_annet + ccl_s_annet)/2]   
additive_dt[, s_gr_nit     := (s_gr_nit + ccl_s_gr_nit)/2]
additive_dt[, s_N_leach    := (s_N_leach + ccl_s_N_leach)/2]
  # drop ccl columns
drop_ccl_cols = c(new_n, 'scenario', 'ccl_s_cc_biomass')
additive_dt   = additive_dt[y_block == 2050, -..drop_ccl_cols]
  # update column names
names = c('s_SOC', 's_cr_grain', 's_cr_residC','s_cc_biomass', 'm_sfdcmp',
          'm_sldcmp', 'm_nfix', 's_annet', 's_gr_nit', 's_N_leach')
new_n = as.vector(outer('add_', colnames(additive_dt[,..names]), paste0))
old_n = colnames(additive_dt[, ..names])
setnames(additive_dt, old = c(old_n), new = c(new_n))

  # JOIN DT
ranger_dt     = ccl_ntill_flux_dt[additive_dt, on = .(gridid = gridid,
                                                                     x = x,
                                                                     y = y,
                                                                     y_block = y_block,
                                                                     ssp = ssp,
                                                                     gcm = gcm,
                                                                     WB_NAME = WB_NAME,
                                                                     WB_REGION = WB_REGION)]

# DIFFERENCE DT | BY HAND
ranger_dt[, s_SOC          := s_SOC - add_s_SOC]
ranger_dt[, d_s_cr_grain   := s_cr_grain - add_s_cr_grain]
ranger_dt[, d_s_cr_residC  := s_cr_residC - add_s_cr_residC]
ranger_dt[, d_s_cc_biomass := s_cc_biomass - add_s_cc_biomass]
ranger_dt[, d_m_sfdcmp     := m_sfdcmp - add_m_sfdcmp]
ranger_dt[, d_m_sldcmp     := m_sldcmp - add_m_sldcmp]
ranger_dt[, d_m_nfix       := m_nfix - add_m_nfix]
ranger_dt[, d_s_annet      := s_annet - add_s_annet]
ranger_dt[, d_s_gr_nit     := s_gr_nit - add_s_gr_nit]
ranger_dt[, d_s_N_leach    := s_N_leach - add_s_N_leach]
# DROP COLUMNS
ranger_dt     = ranger_dt[, -..new_n]
ranger_dt     = ranger_dt[!is.na(scenario)]
# JOIN SITE_WTH 
ranger_dt     = ranger_dt[site_wth, on = .(gridid = gridid, ssp = ssp)]
ranger_dt     = ranger_dt[!is.na(scenario)]
# JOIN INPUT_DT
ranger_dt     = ranger_dt[main_table, on = .(gridid = gridid)]
ranger_dt     = ranger_dt[!is.na(scenario)]

# RANDOM FOREST
set.seed(1234)
ranger_dt[, fertN.amt := NULL]
ranger_dt[, orgN.amt  := NULL]
ccl_SOC   = soc_supervised_clustering(ranger_dt)
ccl_k_SOC = optimal_k(ccl_SOC$UMAP$layout, 5L, ccl_SOC$full_ranger_dt, ccl_SOC$reduced_ranger_dt)

# checks
ccl_SOC$gl_gg
cluster_bplot(ccl_k_SOC$kmeans_dt, 's_SOC')
cluster_map(ccl_k_SOC$kmeans_dt, 'ssp126')
cluster_bplot(ccl_k_SOC$kmeans_dt, 'd_s_cr_residC')

sv_dependence(ccl_SOC$viz, 'd_m_nfix', color_var = 'd_s_cr_residC')
save(ccl_SOC,  ccl_k_SOC,  file = paste(data.path, "ccl-syn-soc_SHAP.Rdata", sep = '/'))

# YIELD
# ranger_dt   = ranger_dt[, s_cr_residC := NULL]
# ranger_dt   = ranger_dt[, d_s_cr_residC := NULL]
# ccl_YIELD   = yield_supervised_clustering(ranger_dt)
# ccl_k_YIELD = optimal_k(ccl_YIELD$UMAP$layout, 6L, ccl_YIELD$full_ranger_dt, ccl_YIELD$reduced_ranger_dt)
# 
# ccl_YIELD$gl_gg
# cluster_bplot(ccl_k_YIELD$kmeans_dt, 'd_s_cr_grain')
# cluster_map(ccl_k_YIELD$kmeans_dt, 'ssp126')
# cluster_bplot(ccl_k_YIELD$kmeans_dt, 'm_nfix')
# save(ccl_YIELD,  ccl_k_YIELD,  file = paste(data.path, "ccl-syn-yield_SHAP.Rdata", sep = '/'))
#-------------------------------------------------------------------------------
# FACTORIAL CONTRIBUTIONS TO COMBINED RESPONSE | Grass Cover Crop
#-------------------------------------------------------------------------------
# INTERPRETATION NOTE: The random forest response (s_SOC) is the difference between
#                      the combined practice scenario, cover crop + no-tillage +
#                      full residue retention, and the 'additive model' which is 
#                      no-tillage + residue added to cover crop under BAU residue management.
# Additive model: ntill_res, ccg
#-------------------------------------------------------------------------------
# TIDY ccg-NTILL RESPONSE DT
# filter to 2050
ccg_ntill_flux_dt = ccg_ntill_flux_dt[y_block <= 2050,]
# mean response for nfix, soil C:N, sfdcmp, slcmp
ccg_ntill_flux_dt[, m_nfix       := mean(m_nfix), by = .(gridid, ssp)]
ccg_ntill_flux_dt[, m_sfdcmp     := mean(m_sfdcmp), by = .(gridid, ssp)]
ccg_ntill_flux_dt[, m_sldcmp     := mean(m_sldcmp), by = .(gridid, ssp)]
# sum cover crop biomass
ccg_ntill_flux_dt[, s_cc_biomass := s_cc_shootC + s_cc_rootC]
# reduce columns
drop_cols = c('m_cr_rootC','m_cr_shootC','m_cr_shootN','m_cr_grain','m_cr_grainN','m_cr_NPP',    
              'm_cr_residC','m_cr_residN','s_cr_rootC','s_cr_shootC','s_cr_shootN', 
              's_cr_grainN','s_cr_residN','m_cc_rootC','m_cc_shootC','m_cc_shootN',
              's_cc_shootN','m_SOC','m_N2O','m_iN2O','m_dN2O','m_GHG', 'm_annet','m_sC.N','m_cr_irr',    
              's_cr_irr','s_N2O','s_iN2O','s_dN2O','s_GHG','s_cc_shootC','s_cc_rootC', 's_cr_NPP')
ccg_ntill_flux_dt = ccg_ntill_flux_dt[y_block == 2050, -..drop_cols]

# EXPECTED ADDITIVE RESPONSE DT
# update column names
names = c('s_SOC', 's_cr_grain', 's_cr_residC','s_cc_rootC', 's_cc_shootC', 'm_sfdcmp',
          'm_sldcmp','m_nfix', 's_annet', 's_gr_nit', 's_N_leach')
new_n = as.vector(outer('ccg_', colnames(ccg_flux_dt[,..names]), paste0))
old_n = colnames(ccg_flux_dt[, ..names])
setnames(ccg_flux_dt, old = c(old_n), new = c(new_n))
cols  = c('gridid', 'x', 'y', 'y_block', 'gcm', 'ssp', 'WB_NAME', 'WB_REGION', new_n)
ccg_flux_dt = ccg_flux_dt[, ..cols]
# reduce cols of ntill_res
cols  = c('gridid', 'x', 'y', 'y_block', 'gcm', 'ssp', 'WB_NAME', 'WB_REGION', old_n)
ntill_res_flux_dt = ntill_res_flux_dt[, ..cols]
# join ntill_res and cover crop table
additive_dt = ntill_res_flux_dt[ccg_flux_dt, on = .(gridid = gridid,
                                                    x = x,
                                                    y = y,
                                                    y_block = y_block,
                                                    ssp = ssp,
                                                    gcm = gcm,
                                                    WB_NAME = WB_NAME,
                                                    WB_REGION = WB_REGION)]
# filter to 2050
additive_dt = additive_dt[y_block <= 2050,]
# mean response for nfix, soil C:N, sfdcmp, slcmp
additive_dt[, m_nfix       := mean(m_nfix),     by = .(gridid, ssp, gcm)]
additive_dt[, ccg_m_nfix   := mean(ccg_m_nfix), by = .(gridid, ssp, gcm)]
additive_dt[, m_sfdcmp     := mean(m_sfdcmp),   by = .(gridid, ssp, gcm)]
additive_dt[, ccg_m_sfdcmp := mean(ccg_m_sfdcmp), by = .(gridid, ssp, gcm)]
additive_dt[, m_sldcmp     := mean(m_sldcmp),     by = .(gridid, ssp, gcm)]
additive_dt[, ccg_m_sldcmp := mean(ccg_m_sldcmp), by = .(gridid, ssp, gcm)]
# sum cover crop biomass
additive_dt[, s_cc_biomass     := s_cc_shootC + s_cc_rootC]
additive_dt[, ccg_s_cc_biomass := ccg_s_cc_shootC + ccg_s_cc_rootC]
# reduce columns
drop_cols = c('s_cc_rootC', 's_cc_shootC', 'ccg_s_cc_rootC','ccg_s_cc_shootC')
additive_dt[, scenario := 'additive-model']
additive_dt = additive_dt[, -..drop_cols]
# additive dt, mean of features | create by hand
additive_dt[, s_SOC        := s_SOC + ccg_s_SOC] # additive response
additive_dt[, s_cr_grain   := (s_cr_grain + ccg_s_cr_grain)/2]
additive_dt[, s_cr_grain   := (s_cr_grain + ccg_s_cr_grain)/2] 
additive_dt[, s_cr_residC  := (s_cr_residC + ccg_s_cr_residC)/2]
additive_dt[, s_cc_biomass := ccg_s_cc_biomass] # only cover crop
additive_dt[, m_sfdcmp     := (m_sfdcmp + ccg_m_sfdcmp)/2] 
additive_dt[, m_sldcmp     := (m_sldcmp + ccg_m_sldcmp)/2] 
additive_dt[, m_nfix       := (m_nfix + ccg_m_nfix)/2] # since non-legume using mean
additive_dt[, s_annet      := (s_annet + ccg_s_annet)/2]   
additive_dt[, s_gr_nit     := (s_gr_nit + ccg_s_gr_nit)/2]
additive_dt[, s_N_leach    := (s_N_leach + ccg_s_N_leach)/2]
# drop ccg columns
drop_ccg_cols = c(new_n, 'scenario', 'ccg_s_cc_biomass')
additive_dt   = additive_dt[y_block == 2050, -..drop_ccg_cols]
# update column names
names = c('s_SOC', 's_cr_grain', 's_cr_residC','s_cc_biomass', 'm_sfdcmp',
          'm_sldcmp', 'm_nfix', 's_annet', 's_gr_nit', 's_N_leach')
new_n = as.vector(outer('add_', colnames(additive_dt[,..names]), paste0))
old_n = colnames(additive_dt[, ..names])
setnames(additive_dt, old = c(old_n), new = c(new_n))

# JOIN DT
ranger_dt     = ccg_ntill_flux_dt[additive_dt, on = .(gridid = gridid,
                                                      x = x,
                                                      y = y,
                                                      y_block = y_block,
                                                      ssp = ssp,
                                                      gcm = gcm,
                                                      WB_NAME = WB_NAME,
                                                      WB_REGION = WB_REGION)]

# DIFFERENCE DT | BY HAND
ranger_dt[, s_SOC        := s_SOC - add_s_SOC]
ranger_dt[, d_s_cr_grain   := s_cr_grain - add_s_cr_grain]
ranger_dt[, d_s_cr_residC  := s_cr_residC - add_s_cr_residC]
ranger_dt[, d_s_cc_biomass := s_cc_biomass - add_s_cc_biomass]
ranger_dt[, d_m_sfdcmp     := m_sfdcmp - add_m_sfdcmp]
ranger_dt[, d_m_sldcmp     := m_sldcmp - add_m_sldcmp]
ranger_dt[, d_m_nfix       := m_nfix - add_m_nfix]
ranger_dt[, d_s_annet      := s_annet - add_s_annet]
ranger_dt[, d_s_gr_nit     := s_gr_nit - add_s_gr_nit]
ranger_dt[, d_s_N_leach    := s_N_leach - add_s_N_leach]
# DROP COLUMNS
ranger_dt     = ranger_dt[, -..new_n]
ranger_dt     = ranger_dt[!is.na(scenario)]
# JOIN SITE_WTH 
ranger_dt     = ranger_dt[site_wth, on = .(gridid = gridid, ssp = ssp)]
ranger_dt     = ranger_dt[!is.na(scenario)]
# JOIN INPUT_DT
ranger_dt     = ranger_dt[main_table, on = .(gridid = gridid)]
ranger_dt     = ranger_dt[!is.na(scenario)]

# RANDOM FOREST
  # SOC
set.seed(1234)
ranger_dt[, fertN.amt := NULL]
ranger_dt[, orgN.amt  := NULL]
ccg_SOC   = soc_supervised_clustering(ranger_dt)
ccg_k_SOC = optimal_k(ccg_SOC$UMAP$layout, 5L, ccg_SOC$full_ranger_dt, ccg_SOC$reduced_ranger_dt)

ccg_SOC$gl_gg
cluster_bplot(ccg_k_SOC$kmeans_dt, 's_SOC')
cluster_map(ccg_k_SOC$kmeans_dt, 'ssp126')
cluster_bplot(ccg_k_SOC$kmeans_dt, 'd_s_N_leach')

sv_dependence(ccg_SOC$viz, 'SOMC_sum_', color_var = 'd_s_cc_biomass')
save(ccg_SOC,  ccg_k_SOC,  file = paste(data.path, "ccg-syn-soc_SHAP.Rdata", sep = '/'))

  # YIELD
# ranger_dt   = ranger_dt[, s_cr_residC := NULL]
# ranger_dt   = ranger_dt[, d_s_cr_residC := NULL]
# ccg_YIELD   = yield_supervised_clustering(ranger_dt)
# ccg_k_YIELD = optimal_k(ccg_YIELD$UMAP$layout, 6L, ccg_YIELD$full_ranger_dt, ccg_YIELD$reduced_ranger_dt)
# 
# ccg_YIELD$gl_gg
# cluster_bplot(ccg_k_YIELD$kmeans_dt, 'd_s_cr_grain')
# cluster_map(ccg_k_YIELD$kmeans_dt, 'ssp126')
# cluster_bplot(ccg_k_YIELD$kmeans_dt, 'hist_bio12')
# save(ccg_YIELD,  ccg_k_YIELD,  file = paste(data.path, "ccg-syn-yield_SHAP.Rdata", sep = '/'))
#-------------------------------------------------------------------------------
# SUPERVISED CLUSTERING YIELD | Residue
#-------------------------------------------------------------------------------
# filter to 2050
rel_res_flux_dt = rel_res_flux_dt[y_block <= 2050,]
# mean response for nfix, soil C:N, sfdcmp, slcmp
rel_res_flux_dt[, m_nfix       := mean(m_nfix), by = .(gridid, ssp)]
rel_res_flux_dt[, m_sfdcmp     := mean(m_sfdcmp), by = .(gridid, ssp)]
rel_res_flux_dt[, m_sldcmp     := mean(m_sldcmp), by = .(gridid, ssp)]
# reduce columns
drop_cols = c('m_cr_rootC','m_cr_shootC','m_cr_shootN','m_cr_grain','m_cr_grainN','m_cr_NPP',    
              'm_cr_residC','m_cr_residN','s_cr_rootC','s_cr_shootC','s_cr_shootN', 
              's_cr_grainN','s_cr_residN','m_cc_rootC','m_cc_shootC','m_cc_shootN',
              's_cc_shootN','m_SOC','m_N2O','m_iN2O','m_dN2O','m_GHG', 'm_annet','m_sC.N','m_cr_irr',    
              's_cr_irr','s_N2O','s_iN2O','s_dN2O','s_GHG','s_cc_shootC','s_cc_rootC', 's_cr_NPP')
rel_res_flux_dt = rel_res_flux_dt[y_block == 2050, -..drop_cols]

# JOIN SITE_WTH 
ranger_dt     = rel_res_flux_dt[site_wth, on = .(gridid = gridid, ssp = ssp)]
ranger_dt     = ranger_dt[!is.na(scenario)]
# JOIN INPUT_DT
ranger_dt     = ranger_dt[main_table, on = .(gridid = gridid)]
ranger_dt     = ranger_dt[!is.na(scenario)]

# RANDOM FOREST
set.seed(1234)
ranger_dt[, fertN.amt := NULL]
ranger_dt[, orgN.amt  := NULL]
ranger_dt[, s_cr_residC := NULL]

res_YIELD   = yield_noCC_supervised_clustering(ranger_dt)
res_k_YIELD = optimal_k(res_YIELD$UMAP$layout, 4L, res_YIELD$full_ranger_dt, res_YIELD$reduced_ranger_dt)

# checks
res_YIELD$gl_gg
res_YIELD$bee_gg

save(res_YIELD,  res_k_YIELD,  file = paste(data.path, "res-yield_SHAP.Rdata", sep = '/'))
#-------------------------------------------------------------------------------
# SUPERVISED CLUSTERING YIELD | No-tillage
#-------------------------------------------------------------------------------
# filter to 2050
rel_ntill_flux_dt = rel_ntill_flux_dt[y_block <= 2050,]
# mean response for nfix, soil C:N, sfdcmp, slcmp
rel_ntill_flux_dt[, m_nfix       := mean(m_nfix), by = .(gridid, ssp)]
rel_ntill_flux_dt[, m_sfdcmp     := mean(m_sfdcmp), by = .(gridid, ssp)]
rel_ntill_flux_dt[, m_sldcmp     := mean(m_sldcmp), by = .(gridid, ssp)]
# reduce columns
drop_cols = c('m_cr_rootC','m_cr_shootC','m_cr_shootN','m_cr_grain','m_cr_grainN','m_cr_NPP',    
              'm_cr_residC','m_cr_residN','s_cr_rootC','s_cr_shootC','s_cr_shootN', 
              's_cr_grainN','s_cr_residN','m_cc_rootC','m_cc_shootC','m_cc_shootN',
              's_cc_shootN','m_SOC','m_N2O','m_iN2O','m_dN2O','m_GHG', 'm_annet','m_sC.N','m_cr_irr',    
              's_cr_irr','s_N2O','s_iN2O','s_dN2O','s_GHG','s_cc_shootC','s_cc_rootC', 's_cr_NPP')
rel_ntill_flux_dt = rel_ntill_flux_dt[y_block == 2050, -..drop_cols]

# JOIN SITE_WTH 
ranger_dt     = rel_ntill_flux_dt[site_wth, on = .(gridid = gridid, ssp = ssp)]
ranger_dt     = ranger_dt[!is.na(scenario)]
# JOIN INPUT_DT
ranger_dt     = ranger_dt[main_table, on = .(gridid = gridid)]
ranger_dt     = ranger_dt[!is.na(scenario)]

# RANDOM FOREST
set.seed(1234)
ranger_dt[, fertN.amt := NULL]
ranger_dt[, orgN.amt  := NULL]
ranger_dt[, s_cr_residC := NULL]

ntill_YIELD   = yield_noCC_supervised_clustering(ranger_dt)
ntill_k_YIELD = optimal_k(ntill_YIELD$UMAP$layout, 4L, ntill_YIELD$full_ranger_dt, ntill_YIELD$reduced_ranger_dt)

# checks
ntill_YIELD$gl_gg
ntill_YIELD$bee_gg

save(ntill_YIELD,  ntill_k_YIELD,  file = paste(data.path, "ntill-yield_SHAP.Rdata", sep = '/'))
#-------------------------------------------------------------------------------
# SUPERVISED CLUSTERING YIELD | No-tillage, residue
#-------------------------------------------------------------------------------
# filter to 2050
rel_ntill_res_flux_dt = rel_ntill_res_flux_dt[y_block <= 2050,]
# mean response for nfix, soil C:N, sfdcmp, slcmp
rel_ntill_res_flux_dt[, m_nfix       := mean(m_nfix), by = .(gridid, ssp)]
rel_ntill_res_flux_dt[, m_sfdcmp     := mean(m_sfdcmp), by = .(gridid, ssp)]
rel_ntill_res_flux_dt[, m_sldcmp     := mean(m_sldcmp), by = .(gridid, ssp)]
# reduce columns
drop_cols = c('m_cr_rootC','m_cr_shootC','m_cr_shootN','m_cr_grain','m_cr_grainN','m_cr_NPP',    
              'm_cr_residC','m_cr_residN','s_cr_rootC','s_cr_shootC','s_cr_shootN', 
              's_cr_grainN','s_cr_residN','m_cc_rootC','m_cc_shootC','m_cc_shootN',
              's_cc_shootN','m_SOC','m_N2O','m_iN2O','m_dN2O','m_GHG', 'm_annet','m_sC.N','m_cr_irr',    
              's_cr_irr','s_N2O','s_iN2O','s_dN2O','s_GHG','s_cc_shootC','s_cc_rootC', 's_cr_NPP')
rel_ntill_res_flux_dt = rel_ntill_res_flux_dt[y_block == 2050, -..drop_cols]

# JOIN SITE_WTH 
ranger_dt     = rel_ntill_res_flux_dt[site_wth, on = .(gridid = gridid, ssp = ssp)]
ranger_dt     = ranger_dt[!is.na(scenario)]
# JOIN INPUT_DT
ranger_dt     = ranger_dt[main_table, on = .(gridid = gridid)]
ranger_dt     = ranger_dt[!is.na(scenario)]

# RANDOM FOREST
set.seed(1234)
ranger_dt[, fertN.amt := NULL]
ranger_dt[, orgN.amt  := NULL]
ranger_dt[, s_cr_residC := NULL]

ntill_res_YIELD   = yield_noCC_supervised_clustering(ranger_dt)
ntill_res_k_YIELD = optimal_k(ntill_res_YIELD$UMAP$layout, 4L, ntill_res_YIELD$full_ranger_dt, ntill_res_YIELD$reduced_ranger_dt)

# checks
ntill_res_YIELD$gl_gg
ntill_res_YIELD$bee_gg

save(ntill_res_YIELD,  ntill_res_k_YIELD,  file = paste(data.path, "ntill-res-yield_SHAP.Rdata", sep = '/'))
#-------------------------------------------------------------------------------
# SUPERVISED CLUSTERING YIELD | Grass cover crop
#-------------------------------------------------------------------------------
# filter to 2050
rel_ccg_flux_dt = rel_ccg_flux_dt[y_block <= 2050,]
# mean response for nfix, soil C:N, sfdcmp, slcmp
rel_ccg_flux_dt[, m_nfix       := mean(m_nfix), by = .(gridid, ssp)]
rel_ccg_flux_dt[, m_sfdcmp     := mean(m_sfdcmp), by = .(gridid, ssp)]
rel_ccg_flux_dt[, m_sldcmp     := mean(m_sldcmp), by = .(gridid, ssp)]
# add cover crop responses
rel_ccg_flux_dt[, s_cc_biomass := s_cc_shootC + s_cc_rootC]
# reduce columns
drop_cols = c('m_cr_rootC','m_cr_shootC','m_cr_shootN','m_cr_grain','m_cr_grainN','m_cr_NPP',    
              'm_cr_residC','m_cr_residN','s_cr_rootC','s_cr_shootC','s_cr_shootN', 
              's_cr_grainN','s_cr_residN','m_cc_rootC','m_cc_shootC','m_cc_shootN',
              's_cc_shootN','m_SOC','m_N2O','m_iN2O','m_dN2O','m_GHG', 'm_annet','m_sC.N','m_cr_irr',    
              's_cr_irr','s_N2O','s_iN2O','s_dN2O','s_GHG','s_cc_shootC','s_cc_rootC', 's_cr_NPP')
rel_ccg_flux_dt = rel_ccg_flux_dt[y_block == 2050, -..drop_cols]

# JOIN SITE_WTH 
ranger_dt     = rel_ccg_flux_dt[site_wth, on = .(gridid = gridid, ssp = ssp)]
ranger_dt     = ranger_dt[!is.na(scenario)]
# JOIN INPUT_DT
ranger_dt     = ranger_dt[main_table, on = .(gridid = gridid)]
ranger_dt     = ranger_dt[!is.na(scenario)]

# RANDOM FOREST
set.seed(1234)
ranger_dt[, fertN.amt := NULL]
ranger_dt[, orgN.amt  := NULL]
ranger_dt[, s_cr_residC := NULL] # testing

ccg_YIELD   = yield_CC_supervised_clustering(ranger_dt)
ccg_k_YIELD = optimal_k(ccg_YIELD$UMAP$layout, 4L, ccg_YIELD$full_ranger_dt, ccg_YIELD$reduced_ranger_dt)

# checks
ccg_YIELD$gl_gg
ccg_YIELD$bee_gg

save(ccg_YIELD,  ccg_k_YIELD,  file = paste(data.path, "ccg-yield_SHAP.Rdata", sep = '/'))
#-------------------------------------------------------------------------------
# SUPERVISED CLUSTERING YIELD | Grass cover crop, residue
#-------------------------------------------------------------------------------
# filter to 2050
rel_ccg_res_flux_dt = rel_ccg_res_flux_dt[y_block <= 2050,]
# mean response for nfix, soil C:N, sfdcmp, slcmp
rel_ccg_res_flux_dt[, m_nfix       := mean(m_nfix), by = .(gridid, ssp)]
rel_ccg_res_flux_dt[, m_sfdcmp     := mean(m_sfdcmp), by = .(gridid, ssp)]
rel_ccg_res_flux_dt[, m_sldcmp     := mean(m_sldcmp), by = .(gridid, ssp)]
# add cover crop responses
rel_ccg_res_flux_dt[, s_cc_biomass := s_cc_shootC + s_cc_rootC]
# reduce columns
drop_cols = c('m_cr_rootC','m_cr_shootC','m_cr_shootN','m_cr_grain','m_cr_grainN','m_cr_NPP',    
              'm_cr_residC','m_cr_residN','s_cr_rootC','s_cr_shootC','s_cr_shootN', 
              's_cr_grainN','s_cr_residN','m_cc_rootC','m_cc_shootC','m_cc_shootN',
              's_cc_shootN','m_SOC','m_N2O','m_iN2O','m_dN2O','m_GHG', 'm_annet','m_sC.N','m_cr_irr',    
              's_cr_irr','s_N2O','s_iN2O','s_dN2O','s_GHG','s_cc_shootC','s_cc_rootC', 's_cr_NPP')
rel_ccg_res_flux_dt = rel_ccg_res_flux_dt[y_block == 2050, -..drop_cols]

# JOIN SITE_WTH 
ranger_dt     = rel_ccg_res_flux_dt[site_wth, on = .(gridid = gridid, ssp = ssp)]
ranger_dt     = ranger_dt[!is.na(scenario)]
# JOIN INPUT_DT
ranger_dt     = ranger_dt[main_table, on = .(gridid = gridid)]
ranger_dt     = ranger_dt[!is.na(scenario)]

# RANDOM FOREST
set.seed(1234)
ranger_dt[, fertN.amt := NULL]
ranger_dt[, orgN.amt  := NULL]
ranger_dt[, s_cr_residC := NULL]

ccg_res_YIELD   = yield_CC_supervised_clustering(ranger_dt)
ccg_res_k_YIELD = optimal_k(ccg_res_YIELD$UMAP$layout, 4L, ccg_res_YIELD$full_ranger_dt, ccg_res_YIELD$reduced_ranger_dt)

# checks
ccg_res_YIELD$gl_gg
ccg_res_YIELD$bee_gg

save(ccg_res_YIELD,  ccg_res_k_YIELD,  file = paste(data.path, "ccg-res-yield_SHAP.Rdata", sep = '/'))
#-------------------------------------------------------------------------------
# SUPERVISED CLUSTERING YIELD | Legume cover crop
#-------------------------------------------------------------------------------
# filter to 2050
rel_ccl_flux_dt = rel_ccl_flux_dt[y_block <= 2050,]
# mean response for nfix, soil C:N, sfdcmp, slcmp
rel_ccl_flux_dt[, m_nfix       := mean(m_nfix), by = .(gridid, ssp)]
rel_ccl_flux_dt[, m_sfdcmp     := mean(m_sfdcmp), by = .(gridid, ssp)]
rel_ccl_flux_dt[, m_sldcmp     := mean(m_sldcmp), by = .(gridid, ssp)]
# add cover crop responses
rel_ccl_flux_dt[, s_cc_biomass := s_cc_shootC + s_cc_rootC]
# reduce columns
drop_cols = c('m_cr_rootC','m_cr_shootC','m_cr_shootN','m_cr_grain','m_cr_grainN','m_cr_NPP',    
              'm_cr_residC','m_cr_residN','s_cr_rootC','s_cr_shootC','s_cr_shootN', 
              's_cr_grainN','s_cr_residN','m_cc_rootC','m_cc_shootC','m_cc_shootN',
              's_cc_shootN','m_SOC','m_N2O','m_iN2O','m_dN2O','m_GHG', 'm_annet','m_sC.N','m_cr_irr',    
              's_cr_irr','s_N2O','s_iN2O','s_dN2O','s_GHG','s_cc_shootC','s_cc_rootC', 's_cr_NPP')
rel_ccl_flux_dt = rel_ccl_flux_dt[y_block == 2050, -..drop_cols]

# JOIN SITE_WTH 
ranger_dt     = rel_ccl_flux_dt[site_wth, on = .(gridid = gridid, ssp = ssp)]
ranger_dt     = ranger_dt[!is.na(scenario)]
# JOIN INPUT_DT
ranger_dt     = ranger_dt[main_table, on = .(gridid = gridid)]
ranger_dt     = ranger_dt[!is.na(scenario)]

# RANDOM FOREST
set.seed(1234)
ranger_dt[, fertN.amt := NULL]
ranger_dt[, orgN.amt  := NULL]
ranger_dt[, s_cr_residC := NULL] 

ccl_YIELD   = yield_CC_supervised_clustering(ranger_dt)
ccl_k_YIELD = optimal_k(ccl_YIELD$UMAP$layout, 4L, ccl_YIELD$full_ranger_dt, ccl_YIELD$reduced_ranger_dt)

# checks
ccl_YIELD$gl_gg
ccl_YIELD$bee_gg

save(ccl_YIELD,  ccl_k_YIELD,  file = paste(data.path, "ccl-yield_SHAP.Rdata", sep = '/'))
#-------------------------------------------------------------------------------
# SUPERVISED CLUSTERING YIELD | Legume cover crop, residue
#-------------------------------------------------------------------------------
# filter to 2050
rel_ccl_res_flux_dt = rel_ccl_res_flux_dt[y_block <= 2050,]
# mean response for nfix, soil C:N, sfdcmp, slcmp
rel_ccl_res_flux_dt[, m_nfix       := mean(m_nfix), by = .(gridid, ssp)]
rel_ccl_res_flux_dt[, m_sfdcmp     := mean(m_sfdcmp), by = .(gridid, ssp)]
rel_ccl_res_flux_dt[, m_sldcmp     := mean(m_sldcmp), by = .(gridid, ssp)]
# add cover crop responses
rel_ccl_res_flux_dt[, s_cc_biomass := s_cc_shootC + s_cc_rootC]
# reduce columns
drop_cols = c('m_cr_rootC','m_cr_shootC','m_cr_shootN','m_cr_grain','m_cr_grainN','m_cr_NPP',    
              'm_cr_residC','m_cr_residN','s_cr_rootC','s_cr_shootC','s_cr_shootN', 
              's_cr_grainN','s_cr_residN','m_cc_rootC','m_cc_shootC','m_cc_shootN',
              's_cc_shootN','m_SOC','m_N2O','m_iN2O','m_dN2O','m_GHG', 'm_annet','m_sC.N','m_cr_irr',    
              's_cr_irr','s_N2O','s_iN2O','s_dN2O','s_GHG','s_cc_shootC','s_cc_rootC', 's_cr_NPP')
rel_ccl_res_flux_dt = rel_ccl_res_flux_dt[y_block == 2050, -..drop_cols]

# JOIN SITE_WTH 
ranger_dt     = rel_ccl_res_flux_dt[site_wth, on = .(gridid = gridid, ssp = ssp)]
ranger_dt     = ranger_dt[!is.na(scenario)]
# JOIN INPUT_DT
ranger_dt     = ranger_dt[main_table, on = .(gridid = gridid)]
ranger_dt     = ranger_dt[!is.na(scenario)]

# RANDOM FOREST
set.seed(1234)
ranger_dt[, fertN.amt := NULL]
ranger_dt[, orgN.amt  := NULL]
ranger_dt[, s_cr_residC := NULL]

ccl_res_YIELD   = yield_CC_supervised_clustering(ranger_dt)
ccl_res_k_YIELD = optimal_k(ccl_res_YIELD$UMAP$layout, 4L, ccl_res_YIELD$full_ranger_dt, ccl_res_YIELD$reduced_ranger_dt)

# checks
ccl_res_YIELD$gl_gg
ccl_res_YIELD$bee_gg

save(ccl_res_YIELD,  ccl_res_k_YIELD,  file = paste(data.path, "ccl-res-yield_SHAP.Rdata", sep = '/'))
#-------------------------------------------------------------------------------
# SUPERVISED CLUSTERING YIELD | Grass cover crop, no-tillage, residue
#-------------------------------------------------------------------------------
# filter to 2050
rel_ccg_ntill_flux_dt = rel_ccg_ntill_flux_dt[y_block <= 2050,]
# mean response for nfix, soil C:N, sfdcmp, slcmp
rel_ccg_ntill_flux_dt[, m_nfix       := mean(m_nfix), by = .(gridid, ssp)]
rel_ccg_ntill_flux_dt[, m_sfdcmp     := mean(m_sfdcmp), by = .(gridid, ssp)]
rel_ccg_ntill_flux_dt[, m_sldcmp     := mean(m_sldcmp), by = .(gridid, ssp)]
# add cover crop responses
rel_ccg_ntill_flux_dt[, s_cc_biomass := s_cc_shootC + s_cc_rootC]
# reduce columns
drop_cols = c('m_cr_rootC','m_cr_shootC','m_cr_shootN','m_cr_grain','m_cr_grainN','m_cr_NPP',    
              'm_cr_residC','m_cr_residN','s_cr_rootC','s_cr_shootC','s_cr_shootN', 
              's_cr_grainN','s_cr_residN','m_cc_rootC','m_cc_shootC','m_cc_shootN',
              's_cc_shootN','m_SOC','m_N2O','m_iN2O','m_dN2O','m_GHG', 'm_annet','m_sC.N','m_cr_irr',    
              's_cr_irr','s_N2O','s_iN2O','s_dN2O','s_GHG','s_cc_shootC','s_cc_rootC', 's_cr_NPP')
rel_ccg_ntill_flux_dt = rel_ccg_ntill_flux_dt[y_block == 2050, -..drop_cols]

# JOIN SITE_WTH 
ranger_dt     = rel_ccg_ntill_flux_dt[site_wth, on = .(gridid = gridid, ssp = ssp)]
ranger_dt     = ranger_dt[!is.na(scenario)]
# JOIN INPUT_DT
ranger_dt     = ranger_dt[main_table, on = .(gridid = gridid)]
ranger_dt     = ranger_dt[!is.na(scenario)]

# RANDOM FOREST
set.seed(1234)
ranger_dt[, fertN.amt := NULL]
ranger_dt[, orgN.amt  := NULL]
ranger_dt[, s_cr_residC := NULL]

ccg_ntill_YIELD   = yield_CC_supervised_clustering(ranger_dt)
ccg_ntill_k_YIELD = optimal_k(ccg_ntill_YIELD$UMAP$layout, 4L, ccg_ntill_YIELD$full_ranger_dt, ccg_ntill_YIELD$reduced_ranger_dt)

# checks
ccg_ntill_YIELD$gl_gg
ccg_ntill_YIELD$bee_gg

save(ccg_ntill_YIELD,  ccg_ntill_k_YIELD,  file = paste(data.path, "ccg-ntill-yield_SHAP.Rdata", sep = '/'))
#-------------------------------------------------------------------------------
# SUPERVISED CLUSTERING YIELD | Legume cover crop, no-tillage, residue
#-------------------------------------------------------------------------------
# filter to 2050
rel_ccl_ntill_flux_dt = rel_ccl_ntill_flux_dt[y_block <= 2050,]
# mean response for nfix, soil C:N, sfdcmp, slcmp
rel_ccl_ntill_flux_dt[, m_nfix       := mean(m_nfix), by = .(gridid, ssp)]
rel_ccl_ntill_flux_dt[, m_sfdcmp     := mean(m_sfdcmp), by = .(gridid, ssp)]
rel_ccl_ntill_flux_dt[, m_sldcmp     := mean(m_sldcmp), by = .(gridid, ssp)]
# add cover crop responses
rel_ccl_ntill_flux_dt[, s_cc_biomass := s_cc_shootC + s_cc_rootC]
# reduce columns
drop_cols = c('m_cr_rootC','m_cr_shootC','m_cr_shootN','m_cr_grain','m_cr_grainN','m_cr_NPP',    
              'm_cr_residC','m_cr_residN','s_cr_rootC','s_cr_shootC','s_cr_shootN', 
              's_cr_grainN','s_cr_residN','m_cc_rootC','m_cc_shootC','m_cc_shootN',
              's_cc_shootN','m_SOC','m_N2O','m_iN2O','m_dN2O','m_GHG', 'm_annet','m_sC.N','m_cr_irr',    
              's_cr_irr','s_N2O','s_iN2O','s_dN2O','s_GHG','s_cc_shootC','s_cc_rootC', 's_cr_NPP')
rel_ccl_ntill_flux_dt = rel_ccl_ntill_flux_dt[y_block == 2050, -..drop_cols]

# JOIN SITE_WTH 
ranger_dt     = rel_ccl_ntill_flux_dt[site_wth, on = .(gridid = gridid, ssp = ssp)]
ranger_dt     = ranger_dt[!is.na(scenario)]
# JOIN INPUT_DT
ranger_dt     = ranger_dt[main_table, on = .(gridid = gridid)]
ranger_dt     = ranger_dt[!is.na(scenario)]

# RANDOM FOREST
set.seed(1234)
ranger_dt[, fertN.amt := NULL]
ranger_dt[, orgN.amt  := NULL]
ranger_dt[, s_cr_residC := NULL]

ccl_ntill_YIELD   = yield_CC_supervised_clustering(ranger_dt)
ccl_ntill_k_YIELD = optimal_k(ccl_ntill_YIELD$UMAP$layout, 4L, ccl_ntill_YIELD$full_ranger_dt, ccl_ntill_YIELD$reduced_ranger_dt)

# checks
ccl_ntill_YIELD$gl_gg
ccl_ntill_YIELD$bee_gg

save(ccl_ntill_YIELD,  ccl_ntill_k_YIELD,  file = paste(data.path, "ccl-ntill-yield_SHAP.Rdata", sep = '/'))
