# filename:    analysis_shap.R
# created:     28 February 2024
# updated:     28 February 2024
# author:      S.C. McClelland
# description: This file analyzes data from DayCent simulations
#              for global cropland soil N2O, CO2 over time for different 
#              interventions and SSPs.
#-------------------------------------------------------------------------------
library(caret)
library(data.table)
library(factoextra)
library(fastshap)
library(patchwork)
library(ranger)
library(rstudioapi)
library(scales)
library(sf)
library(shapviz)
library(terra)
#-------------------------------------------------------------------------------
source('results_functions.R')
source('shap.R')
#-------------------------------------------------------------------------------
options(scipen = 999, digits = 4)
options(rgl.printRglwidget = TRUE)
base.path = dirname(getActiveDocumentContext()$path)
data.path = paste(base.path, 'manuscript-data', sep = '/')
fig.path  = paste(base.path, 'manuscript-figures', sep = '/')
#-------------------------------------------------------------------------------
# LOAD ENSEMBLE DATA
#-------------------------------------------------------------------------------
# NTILL
ntill_flux_dt     = readRDS(paste(data.path, 'ensemble-relative-responses-weighted-mean-ntill.rds', sep = '/'))
# NTILL-RES
ntill_res_flux_dt = readRDS(paste(data.path, 'ensemble-relative-responses-weighted-mean-ntill-res.rds', sep = '/'))
# CCG
ccg_flux_dt       = readRDS(paste(data.path, 'ensemble-relative-responses-weighted-mean-ccg.rds', sep = '/'))
# CCG-RES
ccg_res_flux_dt   = readRDS(paste(data.path, 'ensemble-relative-responses-weighted-mean-ccg-res.rds', sep = '/'))
# CCL
ccl_flux_dt       = readRDS(paste(data.path, 'ensemble-relative-responses-weighted-mean-ccl.rds', sep = '/'))
# CCL-RES
ccl_res_flux_dt   = readRDS(paste(data.path, 'ensemble-relative-responses-weighted-mean-ccl-res.rds', sep = '/'))
# CCG-NTILL-RES
ccg_ntill_flux_dt = readRDS(paste(data.path, 'ensemble-relative-responses-weighted-mean-ccg-ntill.rds', sep = '/'))
# CCL-NTILL-RES
ccl_ntill_flux_dt = readRDS(paste(data.path, 'ensemble-relative-responses-weighted-mean-ccl-ntill.rds', sep = '/'))
#-------------------------------------------------------------------------------
# CONVERT Biomass Units
#-------------------------------------------------------------------------------
Mg_ha = 100L
C_bio = 0.45
# NTILL
ntill_flux_dt[, s_cr_grain := (s_cr_grain/Mg_ha)/C_bio]
ntill_flux_dt[, s_cr_residC := (s_cr_residC/Mg_ha)/C_bio]
ntill_flux_dt[, s_cr_NPP := (s_cr_NPP/Mg_ha)/C_bio]
# NTILL-RES
ntill_res_flux_dt[, s_cr_grain := (s_cr_grain/Mg_ha)/C_bio]
ntill_res_flux_dt[, s_cr_residC := (s_cr_residC/Mg_ha)/C_bio]
ntill_res_flux_dt[, s_cr_NPP := (s_cr_NPP/Mg_ha)/C_bio]
# CCG
ccg_flux_dt[, s_cr_grain := (s_cr_grain/Mg_ha)/C_bio]
ccg_flux_dt[, s_cr_residC := (s_cr_residC/Mg_ha)/C_bio]
ccg_flux_dt[, s_cr_NPP := (s_cr_NPP/Mg_ha)/C_bio]
# CCG-RES
ccg_res_flux_dt[, s_cr_grain := (s_cr_grain/Mg_ha)/C_bio]
ccg_res_flux_dt[, s_cr_residC := (s_cr_residC/Mg_ha)/C_bio]
ccg_res_flux_dt[, s_cr_NPP := (s_cr_NPP/Mg_ha)/C_bio]
# CCL
ccl_flux_dt[, s_cr_grain := (s_cr_grain/Mg_ha)/C_bio]
ccl_flux_dt[, s_cr_residC := (s_cr_residC/Mg_ha)/C_bio]
ccl_flux_dt[, s_cr_NPP := (s_cr_NPP/Mg_ha)/C_bio]
# CCL-RES
ccl_res_flux_dt[, s_cr_grain := (s_cr_grain/Mg_ha)/C_bio]
ccl_res_flux_dt[, s_cr_residC := (s_cr_residC/Mg_ha)/C_bio]
ccl_res_flux_dt[, s_cr_NPP := (s_cr_NPP/Mg_ha)/C_bio]
# CCG-NTILL-RES
ccg_ntill_flux_dt[, s_cr_grain := (s_cr_grain/Mg_ha)/C_bio]
ccg_ntill_flux_dt[, s_cr_residC := (s_cr_residC/Mg_ha)/C_bio]
ccg_ntill_flux_dt[, s_cr_NPP := (s_cr_NPP/Mg_ha)/C_bio]
# CCL-NTILL-RES
ccl_ntill_flux_dt[, s_cr_grain := (s_cr_grain/Mg_ha)/C_bio]
ccl_ntill_flux_dt[, s_cr_residC := (s_cr_residC/Mg_ha)/C_bio]
ccl_ntill_flux_dt[, s_cr_NPP := (s_cr_NPP/Mg_ha)/C_bio]
#-------------------------------------------------------------------------------
# FACTORIAL CONTRIBUTIONS TO COMBINED RESPONSE | Legume Cover Crop
#-------------------------------------------------------------------------------
# INTERPRETATION NOTE: The random forest response (s_SOC) is the difference between
#                      the combined practice scenario, cover crop + no-tillage +
#                      full residue retention, and the 'additive model' which is 
#                      no-tillage + residue added to cover crop and BAU residue management.
#                      Negative values mean there is more SOC sequestration with the combined 
#                      practice while positive values indicate less SOC sequestration.
#                      Features are also the difference between the combined practice and
#                      the 'additive model'.
# Additive model: ntill_res, ccl
#-------------------------------------------------------------------------------
# TIDY CCL-NTILL RESPONSE DT
# filter to 2050
ccl_ntill_flux_dt = ccl_ntill_flux_dt[y_block <= 2050,]
# mean response for nfix, soil C:N, sfdcmp, slcmp
ccl_ntill_flux_dt[, m_nfix       := mean(m_nfix), by = .(gridid, ssp)]
ccl_ntill_flux_dt[, m_sC.N       := mean(m_sC.N), by = .(gridid, ssp)]
ccl_ntill_flux_dt[, m_sfdcmp     := mean(m_sfdcmp), by = .(gridid, ssp)]
ccl_ntill_flux_dt[, m_sldcmp     := mean(m_sldcmp), by = .(gridid, ssp)]
# sum cover crop biomass
ccl_ntill_flux_dt[, s_cc_biomass := s_cc_shootC + s_cc_rootC]
# reduce columns
drop_cols = c('m_cr_rootC','m_cr_shootC','m_cr_shootN','m_cr_grain','m_cr_grainN','m_cr_NPP',    
  'm_cr_residC','m_cr_residN','s_cr_rootC','s_cr_shootC','s_cr_shootN', 
  's_cr_grainN','s_cr_residN','m_cc_rootC','m_cc_shootC','m_cc_shootN',
  's_cc_shootN','m_SOC','m_N2O','m_iN2O','m_dN2O','m_GHG', 'm_annet','m_cr_irr',    
  's_cr_irr','s_N2O','s_iN2O','s_dN2O','s_GHG','s_cc_shootC','s_cc_rootC')
ccl_ntill_flux_dt = ccl_ntill_flux_dt[y_block == 2050, -..drop_cols]

# EXPECTED ADDITIVE RESPONSE DT
  # update column names
names = c('s_SOC', 's_cr_NPP', 's_cr_grain', 's_cr_residC','s_cc_rootC', 's_cc_shootC', 'm_sfdcmp',
          'm_sldcmp', 'm_sC.N', 'm_nfix', 's_annet', 's_gr_nit')
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
additive_dt[, m_nfix       := mean(m_nfix), by = .(gridid, ssp, gcm)]
additive_dt[, ccl_m_nfix   := mean(ccl_m_nfix), by = .(gridid, ssp, gcm)]
additive_dt[, m_sC.N       := mean(m_sC.N), by = .(gridid, ssp, gcm)]
additive_dt[, ccl_m_sC.N   := mean(ccl_m_sC.N), by = .(gridid, ssp, gcm)]
additive_dt[, m_sfdcmp     := mean(m_sfdcmp), by = .(gridid, ssp, gcm)]
additive_dt[, ccl_m_sfdcmp := mean(ccl_m_sfdcmp), by = .(gridid, ssp, gcm)]
additive_dt[, m_sldcmp     := mean(m_sldcmp), by = .(gridid, ssp, gcm)]
additive_dt[, ccl_m_sldcmp := mean(ccl_m_sldcmp), by = .(gridid, ssp, gcm)]
  # sum cover crop biomass
additive_dt[, s_cc_biomass     := s_cc_shootC + s_cc_rootC]
additive_dt[, ccl_s_cc_biomass := ccl_s_cc_shootC + ccl_s_cc_rootC]
  # reduce columns
drop_cols = c('s_cc_rootC', 's_cc_shootC', 'ccl_s_cc_rootC','ccl_s_cc_shootC')
additive_dt[, scenario := 'additive-model']
additive_dt = additive_dt[, -..drop_cols]
  # difference dt | by hand
additive_dt[, s_SOC        := s_SOC + ccl_s_SOC]
additive_dt[, s_cr_grain   := s_cr_grain + ccl_s_cr_grain]
additive_dt[, s_cr_NPP     := s_cr_NPP + ccl_s_cr_NPP]
additive_dt[, s_cr_residC  := s_cr_residC + ccl_s_cr_residC]
additive_dt[, s_cc_biomass := s_cc_biomass + ccl_s_cc_biomass]
additive_dt[, m_sfdcmp     := m_sfdcmp + ccl_m_sfdcmp]
additive_dt[, m_sldcmp     := m_sldcmp + ccl_m_sldcmp]
additive_dt[, m_sC.N       := m_sC.N + ccl_m_sC.N]
additive_dt[, m_nfix       := m_nfix + ccl_m_nfix]
additive_dt[, s_annet      := s_annet + ccl_s_annet]
additive_dt[, s_gr_nit     := s_gr_nit + ccl_s_gr_nit]
  # drop ccl columns
drop_ccl_cols = c(new_n, 'scenario', 'ccl_s_cc_biomass')
additive_dt   = additive_dt[y_block == 2050, -..drop_ccl_cols]
  # update column names
names = c('s_SOC', 's_cr_NPP', 's_cr_grain', 's_cr_residC','s_cc_biomass', 'm_sfdcmp',
          'm_sldcmp', 'm_sC.N', 'm_nfix', 's_annet', 's_gr_nit')
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
ranger_dt[, s_SOC        := s_SOC - add_s_SOC]
ranger_dt[, s_cr_grain   := s_cr_grain - add_s_cr_grain]
ranger_dt[, s_cr_NPP     := s_cr_NPP - add_s_cr_NPP]
ranger_dt[, s_cr_residC  := s_cr_residC - add_s_cr_residC]
ranger_dt[, s_cc_biomass := s_cc_biomass - add_s_cc_biomass]
ranger_dt[, m_sfdcmp     := m_sfdcmp - add_m_sfdcmp]
ranger_dt[, m_sldcmp     := m_sldcmp - add_m_sldcmp]
ranger_dt[, m_sC.N       := m_sC.N - add_m_sC.N]
ranger_dt[, m_nfix       := m_nfix - add_m_nfix]
ranger_dt[, s_annet      := s_annet - add_s_annet]
ranger_dt[, s_gr_nit     := s_gr_nit - add_s_gr_nit]
# DROP COLUMNS
ranger_dt     = ranger_dt[, -..new_n]
ranger_dt     = ranger_dt[!is.na(scenario)]

# RANDOM FOREST
pred_fun                 = function(object, newdata) {
  predict(object, data = newdata)$predictions
}
  # test table
# s             = sample(ranger_dt[,gridid],5000)
# ranger_dt     = ranger_dt[gridid %in% s]
drop_ranger   = c('gridid', 'x', 'y', 'scenario', 'y_block', 'gcm', 'ssp', 'WB_NAME', 'WB_REGION')
r_ranger_dt   = ranger_dt[, -..drop_ranger]
  # correlated variables
correlationMatrix = cor(r_ranger_dt[,c(1:8,10:11)])
highlyCorrelated  = findCorrelation(correlationMatrix, cutoff=0.75, names = TRUE)
print(highlyCorrelated)
  # remove highly correlated
r_ranger_dt[, s_cr_NPP := NULL]
r_ranger_dt[, m_sldcmp := NULL]

  # remove outliers

print('Running ranger.')
soc_ranger       = ranger(dependent.variable.name = 's_SOC', data = r_ranger_dt,
                          importance = 'permutation', keep.inbag = TRUE, seed = 1234)
soc_ranger

soc_features = r_ranger_dt[, -c('s_SOC')]
# SHAP COMPUTATION - N.B. TAKES A LONG TIME
print('Computing SHAP')
soc_SHAP       = explain(soc_ranger, X = soc_features, pred_wrapper = pred_fun, nsim = 10, adjust = TRUE,
                         shap_only = FALSE)
shv.soc_gl     = shapviz(soc_SHAP)
# initial visualization
soc_SHAP_gl_gg = sv_importance(shv.soc_gl)
soc_SHAP_be_gg = sv_importance(shv.soc_gl, kind = "beeswarm", show_numbers = FALSE, bee_width = 0.2)
sv_dependence(shv.soc_gl, 's_gr_nit', color_var = NULL)
sv_dependence(shv.soc_gl, 'm_sC.N', color_var = NULL)
sv_dependence(shv.soc_gl, 's_cr_residC', color_var = NULL)
sv_dependence(shv.soc_gl, 'm_nfix', color_var = NULL)

sv_dependence(shv.soc_gl, 's_cr_residC', color_var = 'm_sC.N')
#-------------------------------------------------------------------------------
# FACTORIAL CONTRIBUTIONS TO COMBINED RESPONSE | Legume Cover Crop 
#-------------------------------------------------------------------------------
# INTERPRETATION NOTE: The random forest response (s_SOC) is the difference between
#                      the combined practice scenario, cover crop + no-tillage +
#                      full residue retention, and the 'additive model' which is 
#                      no-tillage + residue added to cover crop and BAU residue management.
#                      Negative values mean there is more SOC sequestration with the combined 
#                      practice while positive values indicate less SOC sequestration.
#                      Features are also the difference between the combined practice and
#                      the 'additive model'.
# Additive model: ntill, ccl_res
#-------------------------------------------------------------------------------
# TIDY CCL-NTILL RESPONSE DT
# filter to 2050
ccl_ntill_flux_dt = ccl_ntill_flux_dt[y_block <= 2050,]
# mean response for nfix, soil C:N, sfdcmp, slcmp
ccl_ntill_flux_dt[, m_nfix       := mean(m_nfix), by = .(gridid, ssp)]
ccl_ntill_flux_dt[, m_sC.N       := mean(m_sC.N), by = .(gridid, ssp)]
ccl_ntill_flux_dt[, m_sfdcmp     := mean(m_sfdcmp), by = .(gridid, ssp)]
ccl_ntill_flux_dt[, m_sldcmp     := mean(m_sldcmp), by = .(gridid, ssp)]
# sum cover crop biomass
ccl_ntill_flux_dt[, s_cc_biomass := s_cc_shootC + s_cc_rootC]
# reduce columns
drop_cols = c('m_cr_rootC','m_cr_shootC','m_cr_shootN','m_cr_grain','m_cr_grainN','m_cr_NPP',    
              'm_cr_residC','m_cr_residN','s_cr_rootC','s_cr_shootC','s_cr_shootN', 
              's_cr_grainN','s_cr_residN','m_cc_rootC','m_cc_shootC','m_cc_shootN',
              's_cc_shootN','m_SOC','m_N2O','m_iN2O','m_dN2O','m_GHG', 'm_annet','m_cr_irr',    
              's_cr_irr','s_N2O','s_iN2O','s_dN2O','s_GHG','s_cc_shootC','s_cc_rootC')
ccl_ntill_flux_dt = ccl_ntill_flux_dt[y_block == 2050, -..drop_cols]

# EXPECTED ADDITIVE RESPONSE DT
# update column names
names = c('s_SOC', 's_cr_NPP', 's_cr_grain', 's_cr_residC','s_cc_rootC', 's_cc_shootC', 'm_sfdcmp',
          'm_sldcmp', 'm_sC.N', 'm_nfix', 's_annet', 's_gr_nit')
new_n = as.vector(outer('ccl_', colnames(ccl_res_flux_dt[,..names]), paste0))
old_n = colnames(ccl_res_flux_dt[, ..names])
setnames(ccl_res_flux_dt, old = c(old_n), new = c(new_n))
cols  = c('gridid', 'x', 'y', 'y_block', 'gcm', 'ssp', 'WB_NAME', 'WB_REGION', new_n)
ccl_res_flux_dt = ccl_res_flux_dt[, ..cols]
# reduce cols of ntill_res
cols  = c('gridid', 'x', 'y', 'y_block', 'gcm', 'ssp', 'WB_NAME', 'WB_REGION', old_n)
ntill_flux_dt = ntill_flux_dt[, ..cols]
# join ntill_res and cover crop table
additive_dt = ntill_flux_dt[ccl_res_flux_dt, on = .(gridid = gridid,
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
additive_dt[, m_nfix       := mean(m_nfix), by = .(gridid, ssp, gcm)]
additive_dt[, ccl_m_nfix   := mean(ccl_m_nfix), by = .(gridid, ssp, gcm)]
additive_dt[, m_sC.N       := mean(m_sC.N), by = .(gridid, ssp, gcm)]
additive_dt[, ccl_m_sC.N   := mean(ccl_m_sC.N), by = .(gridid, ssp, gcm)]
additive_dt[, m_sfdcmp     := mean(m_sfdcmp), by = .(gridid, ssp, gcm)]
additive_dt[, ccl_m_sfdcmp := mean(ccl_m_sfdcmp), by = .(gridid, ssp, gcm)]
additive_dt[, m_sldcmp     := mean(m_sldcmp), by = .(gridid, ssp, gcm)]
additive_dt[, ccl_m_sldcmp := mean(ccl_m_sldcmp), by = .(gridid, ssp, gcm)]
# sum cover crop biomass
additive_dt[, s_cc_biomass     := s_cc_shootC + s_cc_rootC]
additive_dt[, ccl_s_cc_biomass := ccl_s_cc_shootC + ccl_s_cc_rootC]
# reduce columns
drop_cols = c('s_cc_rootC', 's_cc_shootC', 'ccl_s_cc_rootC','ccl_s_cc_shootC')
additive_dt[, scenario := 'additive-model']
additive_dt = additive_dt[, -..drop_cols]
# difference dt | by hand
additive_dt[, s_SOC        := s_SOC + ccl_s_SOC]
additive_dt[, s_cr_grain   := s_cr_grain + ccl_s_cr_grain]
additive_dt[, s_cr_NPP     := s_cr_NPP + ccl_s_cr_NPP]
additive_dt[, s_cr_residC  := s_cr_residC + ccl_s_cr_residC]
additive_dt[, s_cc_biomass := s_cc_biomass + ccl_s_cc_biomass]
additive_dt[, m_sfdcmp     := m_sfdcmp + ccl_m_sfdcmp]
additive_dt[, m_sldcmp     := m_sldcmp + ccl_m_sldcmp]
additive_dt[, m_sC.N       := m_sC.N + ccl_m_sC.N]
additive_dt[, m_nfix       := m_nfix + ccl_m_nfix]
additive_dt[, s_annet      := s_annet + ccl_s_annet]
additive_dt[, s_gr_nit     := s_gr_nit + ccl_s_gr_nit]
# drop ccl columns
drop_ccl_cols = c(new_n, 'scenario', 'ccl_s_cc_biomass')
additive_dt   = additive_dt[y_block == 2050, -..drop_ccl_cols]
# update column names
names = c('s_SOC', 's_cr_NPP', 's_cr_grain', 's_cr_residC','s_cc_biomass', 'm_sfdcmp',
          'm_sldcmp', 'm_sC.N', 'm_nfix', 's_annet', 's_gr_nit')
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
ranger_dt[, s_SOC        := s_SOC - add_s_SOC]
ranger_dt[, s_cr_grain   := s_cr_grain - add_s_cr_grain]
ranger_dt[, s_cr_NPP     := s_cr_NPP - add_s_cr_NPP]
ranger_dt[, s_cr_residC  := s_cr_residC - add_s_cr_residC]
ranger_dt[, s_cc_biomass := s_cc_biomass - add_s_cc_biomass]
ranger_dt[, m_sfdcmp     := m_sfdcmp - add_m_sfdcmp]
ranger_dt[, m_sldcmp     := m_sldcmp - add_m_sldcmp]
ranger_dt[, m_sC.N       := m_sC.N - add_m_sC.N]
ranger_dt[, m_nfix       := m_nfix - add_m_nfix]
ranger_dt[, s_annet      := s_annet - add_s_annet]
ranger_dt[, s_gr_nit     := s_gr_nit - add_s_gr_nit]
# DROP COLUMNS
ranger_dt     = ranger_dt[, -..new_n]
ranger_dt     = ranger_dt[!is.na(scenario)]

# RANDOM FOREST
pred_fun                 = function(object, newdata) {
  predict(object, data = newdata)$predictions
}
# test table
# s             = sample(ranger_dt[,gridid],5000)
# ranger_dt     = ranger_dt[gridid %in% s]
drop_ranger   = c('gridid', 'x', 'y', 'scenario', 'y_block', 'gcm', 'ssp', 'WB_NAME', 'WB_REGION')
r_ranger_dt   = ranger_dt[, -..drop_ranger]
# correlated variables
correlationMatrix = cor(r_ranger_dt[,c(1:8,10:11)])
highlyCorrelated  = findCorrelation(correlationMatrix, cutoff=0.75, names = TRUE)
print(highlyCorrelated)
# remove highly correlated
r_ranger_dt[, s_cr_NPP   := NULL]
r_ranger_dt[, s_cr_grain := NULL]
r_ranger_dt[, m_sldcmp   := NULL]

# remove outliers

print('Running ranger.')
soc_ranger       = ranger(dependent.variable.name = 's_SOC', data = r_ranger_dt,
                          importance = 'permutation', keep.inbag = TRUE, seed = 1234)
soc_ranger

soc_features = r_ranger_dt[, -c('s_SOC')]
# SHAP COMPUTATION - N.B. TAKES A LONG TIME
print('Computing SHAP')
soc_SHAP       = explain(soc_ranger, X = soc_features, pred_wrapper = pred_fun, nsim = 10, adjust = TRUE,
                         shap_only = FALSE)
shv.soc_gl     = shapviz(soc_SHAP)
# initial visualization
soc_SHAP_gl_gg = sv_importance(shv.soc_gl)
soc_SHAP_be_gg = sv_importance(shv.soc_gl, kind = "beeswarm", show_numbers = FALSE, bee_width = 0.2)
sv_dependence(shv.soc_gl, 's_gr_nit', color_var = NULL)
sv_dependence(shv.soc_gl, 's_annet', color_var = NULL)
sv_dependence(shv.soc_gl, 'm_sC.N', color_var = NULL)
sv_dependence(shv.soc_gl, 'm_nfix', color_var = NULL)

