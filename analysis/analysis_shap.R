# filename:    analysis_shap.R
# created:     28 February 2024
# updated:     15 March 2024
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
library(umap)
#-------------------------------------------------------------------------------
source('results_functions.R')
# source('shap.R')
#-------------------------------------------------------------------------------
options(scipen = 999, digits = 4)
options(rgl.printRglwidget = TRUE)
base.path = dirname(getActiveDocumentContext()$path)
data.path = paste(base.path, 'manuscript-data', sep = '/')
fig.path  = paste(base.path, 'manuscript-figures', sep = '/')
#-------------------------------------------------------------------------------
# LOAD ABSOLUTE DATA
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
# SOIL & CLIMATE DATA
#------------------------------------------------------------------------------
site_wth          = readRDS(paste(data.path, 'ensemble_site_climate_data_decadal.rds', sep = '/'))
# CONSTRAIN TO 2050
site_wth          = site_wth[y_block <= 2050,]
# REDUCE VARS
s_w_cols          = c('gridid', 'y_block', 'ssp', 'bio1', 'hist_bio1', 'bio12', 'hist_bio12',
                      'ELEV', 'SOMC_sum_', 'SLBLKD', 'SLCLAY', 'SLSAND', 'SLPH')
site_wth          = site_wth[, ..s_w_cols]
# ENSEMBLE MEAN OVER TIME
site_wth          = site_wth[, lapply(.SD, mean), .SDcols = s_w_cols[4:13], by = .(gridid, ssp)]
#-------------------------------------------------------------------------------
# CONVERT Biomass Units
#-------------------------------------------------------------------------------
Mg_ha = 100L
C_bio = 0.45
# NTILL-RES
ntill_res_flux_dt[, s_cr_grain  := (s_cr_grain/Mg_ha)/C_bio]
ntill_res_flux_dt[, s_cr_residC := (s_cr_residC/Mg_ha)/C_bio]
# CCG
ccg_flux_dt[, s_cr_grain  := (s_cr_grain/Mg_ha)/C_bio]
ccg_flux_dt[, s_cr_residC := (s_cr_residC/Mg_ha)/C_bio]
ccg_flux_dt[, s_cc_shootC := (s_cc_shootC/Mg_ha)/C_bio]
ccg_flux_dt[, s_cc_rootC  := (s_cc_rootC/Mg_ha)/C_bio]
# CCL
ccl_flux_dt[, s_cr_grain  := (s_cr_grain/Mg_ha)/C_bio]
ccl_flux_dt[, s_cr_residC := (s_cr_residC/Mg_ha)/C_bio]
ccl_flux_dt[, s_cc_shootC := (s_cc_shootC/Mg_ha)/C_bio]
ccl_flux_dt[, s_cc_rootC  := (s_cc_rootC/Mg_ha)/C_bio]
# CCG-NTILL-RES
ccg_ntill_flux_dt[, s_cr_grain  := (s_cr_grain/Mg_ha)/C_bio]
ccg_ntill_flux_dt[, s_cr_residC := (s_cr_residC/Mg_ha)/C_bio]
ccg_ntill_flux_dt[, s_cc_shootC := (s_cc_shootC/Mg_ha)/C_bio]
ccg_ntill_flux_dt[, s_cc_rootC  := (s_cc_rootC/Mg_ha)/C_bio]
# CCL-NTILL-RES
ccl_ntill_flux_dt[, s_cr_grain  := (s_cr_grain/Mg_ha)/C_bio]
ccl_ntill_flux_dt[, s_cr_residC := (s_cr_residC/Mg_ha)/C_bio]
ccl_ntill_flux_dt[, s_cc_shootC := (s_cc_shootC/Mg_ha)/C_bio]
ccl_ntill_flux_dt[, s_cc_rootC  := (s_cc_rootC/Mg_ha)/C_bio]
#-------------------------------------------------------------------------------
# FACTORIAL CONTRIBUTIONS TO COMBINED RESPONSE | Legume Cover Crop
#-------------------------------------------------------------------------------
# INTERPRETATION NOTE: The random forest response (s_SOC) is the difference between
#                      the combined practice scenario, cover crop + no-tillage +
#                      full residue retention, and the 'additive model' which is 
#                      no-tillage + residue added to cover crop under BAU residue management.
#                      Features are XXX.
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
  's_cr_irr','s_N2O','s_iN2O','s_dN2O','s_GHG','s_cc_shootC','s_cc_rootC')
ccl_ntill_flux_dt = ccl_ntill_flux_dt[y_block == 2050, -..drop_cols]

# EXPECTED ADDITIVE RESPONSE DT
  # update column names
names = c('s_SOC', 's_cr_NPP', 's_cr_grain', 's_cr_residC','s_cc_rootC', 's_cc_shootC', 'm_sfdcmp',
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
additive_dt[, s_cr_NPP     := s_cr_NPP + ccl_s_cr_NPP] # additive
additive_dt[, s_cr_residC  := (s_cr_residC + ccl_s_cr_residC)/2]
additive_dt[, s_cc_biomass := s_cc_biomass + ccl_s_cc_biomass] # additive
additive_dt[, m_sfdcmp     := (m_sfdcmp + ccl_m_sfdcmp)/2] 
additive_dt[, m_sldcmp     := (m_sldcmp + ccl_m_sldcmp)/2] 
additive_dt[, m_nfix       := (m_nfix + ccl_m_nfix)/2]
additive_dt[, s_annet      := (s_annet + ccl_s_annet)/2]   
additive_dt[, s_gr_nit     := (s_gr_nit + ccl_s_gr_nit)/2]
additive_dt[, s_N_leach    := (s_N_leach + ccl_s_N_leach)/2]
  # drop ccl columns
drop_ccl_cols = c(new_n, 'scenario', 'ccl_s_cc_biomass')
additive_dt   = additive_dt[y_block == 2050, -..drop_ccl_cols]
  # update column names
names = c('s_SOC', 's_cr_NPP', 's_cr_grain', 's_cr_residC','s_cc_biomass', 'm_sfdcmp',
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
# ranger_dt[, s_SOC        := s_SOC - add_s_SOC]
# ranger_dt[, s_cr_grain   := s_cr_grain - add_s_cr_grain]
# ranger_dt[, s_cr_NPP     := s_cr_NPP - add_s_cr_NPP]
# ranger_dt[, s_cr_residC  := s_cr_residC - add_s_cr_residC]
# ranger_dt[, s_cc_biomass := s_cc_biomass - add_s_cc_biomass]
# ranger_dt[, m_sfdcmp     := m_sfdcmp - add_m_sfdcmp]
# ranger_dt[, m_sldcmp     := m_sldcmp - add_m_sldcmp]
# ranger_dt[, m_nfix       := m_nfix - add_m_nfix]
# ranger_dt[, s_annet      := s_annet - add_s_annet]
# ranger_dt[, s_gr_nit     := s_gr_nit - add_s_gr_nit]
# ranger_dt[, s_N_leach    := s_N_leach - add_s_N_leach]
ranger_dt[, s_SOC          := s_SOC - add_s_SOC]
ranger_dt[, d_s_cr_grain   := s_cr_grain - add_s_cr_grain]
ranger_dt[, d_s_cr_NPP     := s_cr_NPP - add_s_cr_NPP]
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

# RANDOM FOREST
  # REMOVE NPP, includes actual and difference, estimated additive as mean responses
ranger_dt[, s_cr_NPP := NULL]
set.seed(1234)
ccl_SOC   = soc_supervised_clustering(ranger_dt)
ccl_k_SOC = optimal_k(ccl_SOC$UMAP$layout, 3L, ccl_SOC$full_ranger_dt, ccl_SOC$reduced_ranger_dt)
# plotting
  # clustering is hard to tell where synergy didn't happen (though seems like c3)

# REMOVE NPP, only difference + cc bio as additive not mean, estimated additive as mean responses
ranger_dt[, s_cr_NPP := NULL]
set.seed(1234)
ccl_SOC   = soc_supervised_clustering(ranger_dt)
ccl_k_SOC = optimal_k(ccl_SOC$UMAP$layout, 5L, ccl_SOC$full_ranger_dt, ccl_SOC$reduced_ranger_dt)

cluster_bplot(ccl_k_SOC$kmeans_dt, 's_SOC')
cluster_map(ccl_k_SOC$kmeans_dt, 'ssp126')

# NPP, includes actual and difference, ...this one is getting closer but still not what i'd like to see with cluster...
ranger_dt[, s_cr_residC   := NULL]
ranger_dt[, d_s_cr_residC := NULL]
ranger_dt[, s_cr_grain    := NULL]
ranger_dt[, d_s_cr_grain  := NULL]
ranger_dt[, s_cc_biomass  := NULL]
ranger_dt[, d_s_cc_biomass := NULL]
ccl_SOC   = soc_supervised_clustering(ranger_dt)
ccl_k_SOC = optimal_k(ccl_SOC$UMAP$layout, 6L, ccl_SOC$full_ranger_dt, ccl_SOC$reduced_ranger_dt)

cluster_bplot(ccl_k_SOC$kmeans_dt, 's_SOC') # better to present as ES with 95% CI?
cluster_map(ccl_k_SOC$kmeans_dt, 'ssp126')
cluster_bplot(ccl_k_SOC$kmeans_dt, 'd_s_annet')

#-------------------------------------------------------------------------------
# OLD #
#-------------------------------------------------------------------------------
# pred_fun                 = function(object, newdata) {
#   predict(object, data = newdata)$predictions
# }
  # test table
# s             = sample(ranger_dt[,gridid],5000)
# ranger_dt     = ranger_dt[gridid %in% s]
# drop_ranger   = c('gridid', 'x', 'y', 'scenario', 'y_block', 'gcm', 'ssp', 'WB_NAME', 'WB_REGION')
# r_ranger_dt   = ranger_dt[, -..drop_ranger]
#   # correlated variables
# correlationMatrix = cor(r_ranger_dt[,c(1:5,7:8,10:33)])
# highlyCorrelated  = findCorrelation(correlationMatrix, cutoff=0.75, names = TRUE)
# print(highlyCorrelated)
#   # remove highly correlated
# # r_ranger_dt[, s_cr_NPP    := NULL]
# # r_ranger_dt[, m_sfdcmp    := NULL]
# # r_ranger_dt[, m_sldcmp    := NULL]
# # r_ranger_dt[, s_cr_grain  := NULL]
# # r_ranger_dt[, SLBLKD      := NULL]
# # r_ranger_dt[, SLSAND      := NULL]
# # remove soil C:N
# r_ranger_dt[, m_sC.N := NULL]
# r_ranger_dt[, d_m_sC.N := NULL]
# 
# quant     = quantile(r_ranger_dt[,s_SOC], seq(0,1,.01))
# print(quant)
# r_ranger_dt = r_ranger_dt[s_SOC > quant[[2]],]
# r_ranger_dt = r_ranger_dt[s_SOC < quant[[100]],]
# 
# print('Running ranger.')
# soc_ranger       = ranger(dependent.variable.name = 's_SOC', data = r_ranger_dt,
#                           importance = 'permutation', keep.inbag = TRUE, seed = 1234)
# soc_ranger # R2 = 0.93 | diff only
#            # R2 = 0.94 | diff and absolute
# 
# soc_features = r_ranger_dt[, -c('s_SOC')]
# # SHAP COMPUTATION - N.B. TAKES A LONG TIME
# print('Computing SHAP')
# soc_SHAP       = explain(soc_ranger, X = soc_features, pred_wrapper = pred_fun, nsim = 10, adjust = TRUE,
#                          shap_only = FALSE)
# shv.soc_gl     = shapviz(soc_SHAP)
# # initial visualization
# soc_SHAP_gl_gg = sv_importance(shv.soc_gl)
# soc_SHAP_be_gg = sv_importance(shv.soc_gl, kind = "beeswarm", show_numbers = FALSE, bee_width = 0.2)
# 
# # test of diff and absolute vars from combo scenario, with all cor vars, some mean vars, N leach
# sv_dependence(shv.soc_gl, 's_cc_biomass', color_var = NULL) # actual cover crop biomass
# sv_dependence(shv.soc_gl, 's_cr_NPP',     color_var = NULL)     # actual NPP
# sv_dependence(shv.soc_gl, 'SOMC_sum_',    color_var = NULL)    # initial C stock
# sv_dependence(shv.soc_gl, 'd_s_cr_residC',color_var = NULL)
# sv_dependence(shv.soc_gl, 'd_s_annet',    color_var = NULL)    # strong interaction with cc biomass
# sv_dependence(shv.soc_gl, 'd_m_sldcmp',   color_var = NULL)   # decreasing soil decomp
# sv_dependence(shv.soc_gl, 'd_s_cr_NPP',   color_var = NULL)
# 
# # test with all cor vars, some mean vars, N leach
# sv_dependence(shv.soc_gl, 's_annet',    color_var = NULL)  # suggests water use efficiency
# sv_dependence(shv.soc_gl, 's_annet',    color_var = 's_cr_residC')
# sv_dependence(shv.soc_gl, 'hist_bio1',  color_var = NULL) 
# sv_dependence(shv.soc_gl, 's_cr_residC',  color_var = NULL) 
# sv_dependence(shv.soc_gl, 'hist_bio12',  color_var = NULL) 
# sv_dependence(shv.soc_gl, 'SOMC_sum_',  color_var = NULL) 
# sv_dependence(shv.soc_gl, 'm_sldcmp',  color_var = NULL) 
# sv_dependence(shv.soc_gl, 'SLSAND',  color_var = NULL)
# try umap and clustering
# make nstart 100
# looks ok, 72%, with three clusters (1 and 2 most different)
#-------------------------------------------------------------------------------
# FACTORIAL CONTRIBUTIONS TO COMBINED RESPONSE | Grass Cover Crop
#-------------------------------------------------------------------------------
# INTERPRETATION NOTE: The random forest response (s_SOC) is the difference between
#                      the combined practice scenario, cover crop + no-tillage +
#                      full residue retention, and the 'additive model' which is 
#                      no-tillage + residue added to cover crop under BAU residue management.
#                      Features are XXX.
# Additive model: ntill_res, ccg
#-------------------------------------------------------------------------------
# TIDY CCG-NTILL RESPONSE DT
# filter to 2050
ccg_ntill_flux_dt = ccg_ntill_flux_dt[y_block <= 2050,]
# mean response for nfix, soil C:N, sfdcmp, slcmp
ccg_ntill_flux_dt[, m_nfix       := mean(m_nfix), by = .(gridid, ssp)]
ccg_ntill_flux_dt[, m_sC.N       := mean(m_sC.N), by = .(gridid, ssp)]
ccg_ntill_flux_dt[, m_sfdcmp     := mean(m_sfdcmp), by = .(gridid, ssp)]
ccg_ntill_flux_dt[, m_sldcmp     := mean(m_sldcmp), by = .(gridid, ssp)]
# sum cover crop biomass
ccg_ntill_flux_dt[, s_cc_biomass := s_cc_shootC + s_cc_rootC]
# reduce columns
drop_cols = c('m_cr_rootC','m_cr_shootC','m_cr_shootN','m_cr_grain','m_cr_grainN','m_cr_NPP',    
              'm_cr_residC','m_cr_residN','s_cr_rootC','s_cr_shootC','s_cr_shootN', 
              's_cr_grainN','s_cr_residN','m_cc_rootC','m_cc_shootC','m_cc_shootN',
              's_cc_shootN','m_SOC','m_N2O','m_iN2O','m_dN2O','m_GHG', 'm_annet','m_cr_irr',    
              's_cr_irr','s_N2O','s_iN2O','s_dN2O','s_GHG','s_cc_shootC','s_cc_rootC')
ccg_ntill_flux_dt = ccg_ntill_flux_dt[y_block == 2050, -..drop_cols]

# EXPECTED ADDITIVE RESPONSE DT
# update column names
names = c('s_SOC', 's_cr_NPP', 's_cr_grain', 's_cr_residC','s_cc_rootC', 's_cc_shootC', 'm_sfdcmp',
          'm_sldcmp', 'm_sC.N', 'm_nfix', 's_annet', 's_gr_nit')
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
additive_dt[, m_nfix       := mean(m_nfix), by = .(gridid, ssp, gcm)]
additive_dt[, ccg_m_nfix   := mean(ccg_m_nfix), by = .(gridid, ssp, gcm)]
additive_dt[, m_sC.N       := mean(m_sC.N), by = .(gridid, ssp, gcm)]
additive_dt[, ccg_m_sC.N   := mean(ccg_m_sC.N), by = .(gridid, ssp, gcm)]
additive_dt[, m_sfdcmp     := mean(m_sfdcmp), by = .(gridid, ssp, gcm)]
additive_dt[, ccg_m_sfdcmp := mean(ccg_m_sfdcmp), by = .(gridid, ssp, gcm)]
additive_dt[, m_sldcmp     := mean(m_sldcmp), by = .(gridid, ssp, gcm)]
additive_dt[, ccg_m_sldcmp := mean(ccg_m_sldcmp), by = .(gridid, ssp, gcm)]
# sum cover crop biomass
additive_dt[, s_cc_biomass     := s_cc_shootC + s_cc_rootC]
additive_dt[, ccg_s_cc_biomass := ccg_s_cc_shootC + ccg_s_cc_rootC]
# reduce columns
drop_cols = c('s_cc_rootC', 's_cc_shootC', 'ccg_s_cc_rootC','ccg_s_cc_shootC')
additive_dt[, scenario := 'additive-model']
additive_dt = additive_dt[, -..drop_cols]
# difference dt | by hand
additive_dt[, s_SOC        := s_SOC + ccg_s_SOC]
additive_dt[, s_cr_grain   := s_cr_grain + ccg_s_cr_grain]
additive_dt[, s_cr_NPP     := s_cr_NPP + ccg_s_cr_NPP]
additive_dt[, s_cr_residC  := s_cr_residC + ccg_s_cr_residC]
additive_dt[, s_cc_biomass := s_cc_biomass + ccg_s_cc_biomass]
additive_dt[, m_sfdcmp     := m_sfdcmp + ccg_m_sfdcmp]
additive_dt[, m_sldcmp     := m_sldcmp + ccg_m_sldcmp]
additive_dt[, m_sC.N       := m_sC.N + ccg_m_sC.N]
additive_dt[, m_nfix       := m_nfix + ccg_m_nfix]
additive_dt[, s_annet      := s_annet + ccg_s_annet]
additive_dt[, s_gr_nit     := s_gr_nit + ccg_s_gr_nit]
# drop ccl columns
drop_ccg_cols = c(new_n, 'scenario', 'ccg_s_cc_biomass')
additive_dt   = additive_dt[y_block == 2050, -..drop_ccg_cols]
# update column names
names = c('s_SOC', 's_cr_NPP', 's_cr_grain', 's_cr_residC','s_cc_biomass', 'm_sfdcmp',
          'm_sldcmp', 'm_sC.N', 'm_nfix', 's_annet', 's_gr_nit')
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

# ranger_dt     = ccg_ntill_flux_dt[additive_dt[,.(gridid, x, y,
#                                                  y_block, gcm, ssp,
#                                                  WB_NAME, WB_REGION,
#                                                  add_s_SOC)], on = .(gridid = gridid,
#                                                      x = x,
#                                                      y = y,
#                                                      y_block = y_block,
#                                                      ssp = ssp,
#                                                      gcm = gcm,
#                                                      WB_NAME = WB_NAME,
#                                                      WB_REGION = WB_REGION)]
# DIFFERENCE DT | BY HAND
ranger_dt[, s_SOC        := s_SOC - add_s_SOC]
# ranger_dt[, add_s_SOC    := NULL]
# DEPRECATED FOR NOW
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
# # DROP COLUMNS
ranger_dt     = ranger_dt[, -..new_n]
ranger_dt     = ranger_dt[!is.na(scenario)]
# JOIN SITE_WTH | N.B. this removes historical observations...
ranger_dt     = ranger_dt[site_wth, on = .(gridid = gridid, ssp = ssp)]
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
correlationMatrix = cor(r_ranger_dt[,c(1:5,7:8,10:21)])
highlyCorrelated  = findCorrelation(correlationMatrix, cutoff=0.75, names = TRUE)
print(highlyCorrelated)
# remove highly correlated
r_ranger_dt[, s_cr_NPP := NULL]
r_ranger_dt[, m_sfdcmp := NULL]
r_ranger_dt[, m_sldcmp := NULL]
r_ranger_dt[, s_cr_grain := NULL]
r_ranger_dt[, SLCLAY     := NULL]
r_ranger_dt[, SLSAND     := NULL]
# remove soil C:N
r_ranger_dt[, m_sC.N := NULL]

quant     = quantile(r_ranger_dt[,s_SOC], seq(0,1,.01))
print(quant)
r_ranger_dt = r_ranger_dt[s_SOC > quant[[2]],]
r_ranger_dt = r_ranger_dt[s_SOC < quant[[100]],]


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

# all cases
sv_dependence(shv.soc_gl, 'SOMC_sum_',    color_var = NULL) # lower soil C stocks benefit then tapers
sv_dependence(shv.soc_gl, 's_gr_nit',     color_var = NULL) # unclear + outliers
sv_dependence(shv.soc_gl, 's_cr_residC', color_var = NULL)  # when resid decreased, not synergistic?
sv_dependence(shv.soc_gl, 's_cc_biomass', color_var = NULL) # looks like biomass increased
sv_dependence(shv.soc_gl, 'hist_bio1', color_var = NULL)    # benefit at 10-20C?
