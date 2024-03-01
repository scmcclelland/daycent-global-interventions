# filename:    analysis_aggregated.R
# created:     06 March 2023
# updated:     28 February 2024
# author:      S.C. McClelland
# description: This file analyzes data from DayCent simulations
#              for global cropland soil N2O, CO2 over time for different 
#              interventions and SSPs.
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
fig.path  = paste(base.path, 'manuscript-figures', sep = '/')
#-------------------------------------------------------------------------------
# LOAD GCM DATA
#-------------------------------------------------------------------------------
# RESIDUE
residue_gcm_flux_dt   = readRDS(paste(data.path, 'gcm-relative-responses-weighted-mean-res.rds', sep = '/'))
# NTILL
ntill_gcm_flux_dt     = readRDS(paste(data.path, 'gcm-relative-responses-weighted-mean-ntill.rds', sep = '/'))
# NTILL-RES
ntill_res_gcm_flux_dt = readRDS(paste(data.path, 'gcm-relative-responses-weighted-mean-ntill-res.rds', sep = '/'))
# CCG
ccg_gcm_flux_dt       = readRDS(paste(data.path, 'gcm-relative-responses-weighted-mean-ccg.rds', sep = '/'))
# CCG-RES
ccg_res_gcm_flux_dt   = readRDS(paste(data.path, 'gcm-relative-responses-weighted-mean-ccg-res.rds', sep = '/'))
# CCL
ccl_gcm_flux_dt       = readRDS(paste(data.path, 'gcm-relative-responses-weighted-mean-ccl.rds', sep = '/'))
# CCL-RES
ccl_res_gcm_flux_dt   = readRDS(paste(data.path, 'gcm-relative-responses-weighted-mean-ccl-res.rds', sep = '/'))
# CCG-NTILL-RES
ccg_ntill_gcm_flux_dt = readRDS(paste(data.path, 'gcm-relative-responses-weighted-mean-ccg-ntill.rds', sep = '/'))
# CCL-NTILL-RES
ccl_ntill_gcm_flux_dt = readRDS(paste(data.path, 'gcm-relative-responses-weighted-mean-ccl-ntill.rds', sep = '/'))
#-------------------------------------------------------------------------------
# UPDATE Region Names
#-------------------------------------------------------------------------------
crop_area_dt         = readRDS(paste(data.path, 'crop-area-country-ipcc-region.rds', sep = '/'))
cc_crop_area_dt      = readRDS(paste(data.path, 'cover-crop-crop-area-country-ipcc-region.rds', sep = '/'))
# RESIDUE
residue_gcm_flux_dt  = residue_gcm_flux_dt[crop_area_dt[, .(WB_NAME, IPCC_NAME)], 
                                            on = .(WB_NAME = WB_NAME) ]
residue_gcm_flux_dt  = residue_gcm_flux_dt[!is.na(gcm)]
# NTILL
ntill_gcm_flux_dt    = ntill_gcm_flux_dt[crop_area_dt[, .(WB_NAME, IPCC_NAME)], 
                                           on = .(WB_NAME = WB_NAME) ]
ntill_gcm_flux_dt    = ntill_gcm_flux_dt[!is.na(gcm)]
# NTILL-RES
ntill_res_gcm_flux_dt= ntill_res_gcm_flux_dt[crop_area_dt[, .(WB_NAME, IPCC_NAME)], 
                                         on = .(WB_NAME = WB_NAME) ]
ntill_res_gcm_flux_dt= ntill_res_gcm_flux_dt[!is.na(gcm)]
# CCG
ccg_gcm_flux_dt      = ccg_gcm_flux_dt[cc_crop_area_dt[, .(WB_NAME, IPCC_NAME)], 
                                           on = .(WB_NAME = WB_NAME) ]
ccg_gcm_flux_dt      = ccg_gcm_flux_dt[!is.na(gcm)]
# CCG-RES
ccg_res_gcm_flux_dt  = ccg_res_gcm_flux_dt[cc_crop_area_dt[, .(WB_NAME, IPCC_NAME)], 
                                       on = .(WB_NAME = WB_NAME) ]
ccg_res_gcm_flux_dt  = ccg_res_gcm_flux_dt[!is.na(gcm)]
# CCL
ccl_gcm_flux_dt      = ccl_gcm_flux_dt[cc_crop_area_dt[, .(WB_NAME, IPCC_NAME)], 
                                           on = .(WB_NAME = WB_NAME) ]
ccl_gcm_flux_dt      = ccl_gcm_flux_dt[!is.na(gcm)]
# CCL-RES
ccl_res_gcm_flux_dt  = ccl_res_gcm_flux_dt[cc_crop_area_dt[, .(WB_NAME, IPCC_NAME)], 
                                       on = .(WB_NAME = WB_NAME) ]
ccl_res_gcm_flux_dt  = ccl_res_gcm_flux_dt[!is.na(gcm)]
# CCG-NTILL
ccg_ntill_gcm_flux_dt= ccg_ntill_gcm_flux_dt[cc_crop_area_dt[, .(WB_NAME, IPCC_NAME)], 
                                           on = .(WB_NAME = WB_NAME) ]
ccg_ntill_gcm_flux_dt= ccg_ntill_gcm_flux_dt[!is.na(gcm)]
# CCL-NTILL
ccl_ntill_gcm_flux_dt= ccl_ntill_gcm_flux_dt[cc_crop_area_dt[, .(WB_NAME, IPCC_NAME)], 
                                           on = .(WB_NAME = WB_NAME) ]
ccl_ntill_gcm_flux_dt= ccl_ntill_gcm_flux_dt[!is.na(gcm)]
#-------------------------------------------------------------------------------
# CONVERT Biomass Units
#-------------------------------------------------------------------------------
Mg_ha = 100L
kg_ha = 10L
C_bio = 0.45
# RESIDUE
residue_gcm_flux_dt[, m_cr_grain := (m_cr_grain*kg_ha)/C_bio]
residue_gcm_flux_dt[, s_cr_grain := (s_cr_grain/Mg_ha)/C_bio]
# NTILL
ntill_gcm_flux_dt[, m_cr_grain := (m_cr_grain*kg_ha)/C_bio]
ntill_gcm_flux_dt[, s_cr_grain := (s_cr_grain/Mg_ha)/C_bio]
# NTILL-RES
ntill_res_gcm_flux_dt[, m_cr_grain := (m_cr_grain*kg_ha)/C_bio]
ntill_res_gcm_flux_dt[, s_cr_grain := (s_cr_grain/Mg_ha)/C_bio]
# CCG
ccg_gcm_flux_dt[, m_cr_grain := (m_cr_grain*kg_ha)/C_bio]
ccg_gcm_flux_dt[, s_cr_grain := (s_cr_grain/Mg_ha)/C_bio]
# CCG-RES
ccg_res_gcm_flux_dt[, m_cr_grain := (m_cr_grain*kg_ha)/C_bio]
ccg_res_gcm_flux_dt[, s_cr_grain := (s_cr_grain/Mg_ha)/C_bio]
# CCL
ccl_gcm_flux_dt[, m_cr_grain := (m_cr_grain*kg_ha)/C_bio]
ccl_gcm_flux_dt[, s_cr_grain := (s_cr_grain/Mg_ha)/C_bio]
# CCL-RES
ccl_res_gcm_flux_dt[, m_cr_grain := (m_cr_grain*kg_ha)/C_bio]
ccl_res_gcm_flux_dt[, s_cr_grain := (s_cr_grain/Mg_ha)/C_bio]
# CCG-NTILL-RES
ccg_ntill_gcm_flux_dt[, m_cr_grain := (m_cr_grain*kg_ha)/C_bio]
ccg_ntill_gcm_flux_dt[, s_cr_grain := (s_cr_grain/Mg_ha)/C_bio]
# CCL-NTILL-RES
ccl_ntill_gcm_flux_dt[, m_cr_grain := (m_cr_grain*kg_ha)/C_bio]
ccl_ntill_gcm_flux_dt[, s_cr_grain := (s_cr_grain/Mg_ha)/C_bio]
#-------------------------------------------------------------------------------
# MEAN GRID CELL RESPONSE BY DECADE | GLOBAL
#-------------------------------------------------------------------------------
# GLOBAL
  # RESIDUE
d_res_gl_mn  = global_decadal_mean(residue_gcm_flux_dt)
d_res_gl_mCI = global_decadal_mCI(d_res_gl_mn)
d_res_gl_mn[, scenario := 'res']
d_res_gl_mCI[, scenario := 'res']
  # NTILL
d_ntill_gl_mn  = global_decadal_mean(ntill_gcm_flux_dt)
d_ntill_gl_mCI = global_decadal_mCI(d_ntill_gl_mn)
d_ntill_gl_mn[, scenario := 'ntill']
d_ntill_gl_mCI[, scenario := 'ntill']
  # NTILL-RES
d_ntill_res_gl_mn  = global_decadal_mean(ntill_res_gcm_flux_dt)
d_ntill_res_gl_mCI = global_decadal_mCI(d_ntill_res_gl_mn)
d_ntill_res_gl_mn[, scenario := 'ntill-res']
d_ntill_res_gl_mCI[, scenario := 'ntill-res']
  # CCG
d_ccg_gl_mn  = global_decadal_mean(ccg_gcm_flux_dt)
d_ccg_gl_mCI = global_decadal_mCI(d_ccg_gl_mn)
d_ccg_gl_mn[, scenario := 'ccg']
d_ccg_gl_mCI[, scenario := 'ccg']
  # CCG-RES
d_ccg_res_gl_mn  = global_decadal_mean(ccg_res_gcm_flux_dt)
d_ccg_res_gl_mCI = global_decadal_mCI(d_ccg_res_gl_mn)
d_ccg_res_gl_mn[, scenario := 'ccg-res']
d_ccg_res_gl_mCI[, scenario := 'ccg-res']
  # CCL
d_ccl_gl_mn  = global_decadal_mean(ccl_gcm_flux_dt)
d_ccl_gl_mCI = global_decadal_mCI(d_ccl_gl_mn)
d_ccl_gl_mn[, scenario := 'ccl']
d_ccl_gl_mCI[, scenario := 'ccl']
  # CCL-RES
d_ccl_res_gl_mn  = global_decadal_mean(ccl_res_gcm_flux_dt)
d_ccl_res_gl_mCI = global_decadal_mCI(d_ccl_res_gl_mn)
d_ccl_res_gl_mn[, scenario := 'ccl-res']
d_ccl_res_gl_mCI[, scenario := 'ccl-res']
  # CCG-NTILL-RES
d_ccg_ntill_gl_mn  = global_decadal_mean(ccg_ntill_gcm_flux_dt)
d_ccg_ntill_gl_mCI = global_decadal_mCI(d_ccg_ntill_gl_mn)
d_ccg_ntill_gl_mn[, scenario := 'ccg-ntill']
d_ccg_ntill_gl_mCI[, scenario := 'ccg-ntill']
  # CCL-NTILL-RES
d_ccl_ntill_gl_mn  = global_decadal_mean(ccl_ntill_gcm_flux_dt)
d_ccl_ntill_gl_mCI = global_decadal_mCI(d_ccl_ntill_gl_mn)
d_ccl_ntill_gl_mn[, scenario := 'ccl-ntill']
d_ccl_ntill_gl_mCI[, scenario := 'ccl-ntill']

# global response across gcm
d_g_mn = rbind(d_res_gl_mn, d_ntill_gl_mn, d_ntill_res_gl_mn, d_ccg_gl_mn, d_ccg_res_gl_mn,
               d_ccl_gl_mn, d_ccl_res_gl_mn,d_ccg_ntill_gl_mn, d_ccl_ntill_gl_mn)
d_g_mn[, IPCC_NAME := 'GLOBE']
d_g_mn[, m_GHG_units := 'Mg CO2-eq ha-1 yr-1']
d_g_mn[, s_GHG_units := 'Mg CO2-eq ha-1']
d_g_mn[, m_biomass_units := 'kg C ha-1 yr-1']
d_g_mn[, s_biomass_units := 'Mg C ha-1']
setcolorder(d_g_mn, c('IPCC_NAME','gcm', 'ssp', 'scenario', 'y_block'))
fwrite(d_g_mn, paste(data.path, 'gcm-decadal-global-crop-soil-responses.csv', sep = '/'))

# global response mean, CI
d_g_mCI = rbind(d_res_gl_mCI, d_ntill_gl_mCI, d_ntill_res_gl_mCI, d_ccg_gl_mCI, d_ccg_res_gl_mCI,
                d_ccl_gl_mCI, d_ccl_res_gl_mCI,d_ccg_ntill_gl_mCI, d_ccl_ntill_gl_mCI)
d_g_mCI[, IPCC_NAME := 'GLOBE']
d_g_mCI[, m_GHG_units := 'Mg CO2-eq ha-1 yr-1']
d_g_mCI[, s_GHG_units := 'Mg CO2-eq ha-1']
d_g_mCI[, m_biomass_units := 'kg C ha-1 yr-1']
d_g_mCI[, s_biomass_units := 'Mg C ha-1']
setcolorder(d_g_mCI, c('IPCC_NAME','ssp', 'scenario', 'y_block'))
fwrite(d_g_mCI, paste(data.path, 'ensemble-decadal-global-crop-soil-responses.csv', sep = '/'))

# REGION CI
  # RESIDUE
d_res_r_mn  = regional_decadal_mean(residue_gcm_flux_dt)
d_res_r_mCI = regional_decadal_mCI(d_res_r_mn)
d_res_r_mn[, scenario := 'res']
d_res_r_mCI[, scenario := 'res']
  # NTILL
d_ntill_r_mn  = regional_decadal_mean(ntill_gcm_flux_dt)
d_ntill_r_mCI = regional_decadal_mCI(d_ntill_r_mn)
d_ntill_r_mn[, scenario := 'ntill']
d_ntill_r_mCI[, scenario := 'ntill']
  # NTILL-RES
d_ntill_res_r_mn  = regional_decadal_mean(ntill_res_gcm_flux_dt)
d_ntill_res_r_mCI = regional_decadal_mCI(d_ntill_res_r_mn)
d_ntill_res_r_mn[, scenario := 'ntill-res']
d_ntill_res_r_mCI[, scenario := 'ntill-res']
  # CCG
d_ccg_r_mn  = regional_decadal_mean(ccg_gcm_flux_dt)
d_ccg_r_mCI = regional_decadal_mCI(d_ccg_r_mn)
d_ccg_r_mn[, scenario := 'ccg']
d_ccg_r_mCI[, scenario := 'ccg']
  # CCG-RES
d_ccg_res_r_mn  = regional_decadal_mean(ccg_res_gcm_flux_dt)
d_ccg_res_r_mCI = regional_decadal_mCI(d_ccg_res_r_mn)
d_ccg_res_r_mn[, scenario := 'ccg-res']
d_ccg_res_r_mCI[, scenario := 'ccg-res']
  # CCL
d_ccl_r_mn  = regional_decadal_mean(ccl_gcm_flux_dt)
d_ccl_r_mCI = regional_decadal_mCI(d_ccl_r_mn)
d_ccl_r_mn[, scenario := 'ccl']
d_ccl_r_mCI[, scenario := 'ccl']
  # CCL-RES
d_ccl_res_r_mn  = regional_decadal_mean(ccl_res_gcm_flux_dt)
d_ccl_res_r_mCI = regional_decadal_mCI(d_ccl_res_r_mn)
d_ccl_res_r_mn[, scenario := 'ccl-res']
d_ccl_res_r_mCI[, scenario := 'ccl-res']
  # CCG-NTILL-RES
d_ccg_ntill_r_mn  = regional_decadal_mean(ccg_ntill_gcm_flux_dt)
d_ccg_ntill_r_mCI = regional_decadal_mCI(d_ccg_ntill_r_mn)
d_ccg_ntill_r_mn[, scenario := 'ccg-ntill']
d_ccg_ntill_r_mCI[, scenario := 'ccg-ntill']
  # CCL-NTILL-RES
d_ccl_ntill_r_mn  = regional_decadal_mean(ccl_ntill_gcm_flux_dt)
d_ccl_ntill_r_mCI = regional_decadal_mCI(d_ccl_ntill_r_mn)
d_ccl_ntill_r_mn[, scenario := 'ccl-ntill']
d_ccl_ntill_r_mCI[, scenario := 'ccl-ntill']

# regional response across gcm
d_r_mn = rbind(d_res_r_mn, d_ntill_r_mn, d_ntill_res_r_mn, d_ccg_r_mn, d_ccg_res_r_mn,
                d_ccl_r_mn, d_ccl_res_r_mn, d_ccg_ntill_r_mn, d_ccl_ntill_r_mn)
d_r_mn[, m_GHG_units := 'Mg CO2-eq ha-1 yr-1']
d_r_mn[, s_GHG_units := 'Mg CO2-eq ha-1']
d_r_mn[, m_biomass_units := 'kg C ha-1 yr-1']
d_r_mn[, s_biomass_units := 'Mg C ha-1']
setcolorder(d_r_mn, c('IPCC_NAME','ssp', 'gcm','scenario', 'y_block'))
fwrite(d_r_mn, paste(data.path, 'gcm-decadal-IPCC-region-crop-soil-responses.csv', sep = '/'))

# regional response mean, CI
d_r_mCI = rbind(d_res_r_mCI, d_ntill_r_mCI, d_ntill_res_r_mCI, d_ccg_r_mCI, d_ccg_res_r_mCI,
                d_ccl_r_mCI, d_ccl_res_r_mCI,d_ccg_ntill_r_mCI, d_ccl_ntill_r_mCI)
d_r_mCI[, m_GHG_units := 'Mg CO2-eq ha-1 yr-1']
d_r_mCI[, s_GHG_units := 'Mg CO2-eq ha-1']
d_r_mCI[, m_biomass_units := 'kg C ha-1 yr-1']
d_r_mCI[, s_biomass_units := 'Mg C ha-1']
setcolorder(d_r_mCI, c('IPCC_NAME','ssp', 'scenario', 'y_block'))
fwrite(d_r_mCI, paste(data.path, 'ensemble-decadal-IPCC-region-crop-soil-responses.csv', sep = '/'))
#-------------------------------------------------------------------------------
# FACTORIAL CONTRIBUTIONS TO COMBINED RESPONSE | Legume Cover Crop
#-------------------------------------------------------------------------------
# LOAD VARIABLE DATA
  # NTILL-RES
ntill_res_gcm_var_flux_dt = readRDS(paste(data.path, 'gcm-relative-driver-responses-weighted-mean-ntill-res.rds', sep = '/'))
  # CCL
ccl_gcm_var_flux_dt       = readRDS(paste(data.path, 'gcm-relative-driver-responses-weighted-mean-ccl.rds', sep = '/'))
  # CCL-NTILL-RES
ccl_ntill_gcm_var_flux_dt = readRDS(paste(data.path, 'gcm-relative-driver-responses-weighted-mean-ccl-ntill.rds', sep = '/'))
# BIND DT
ntill_res_gcm_dt = ntill_res_gcm_flux_dt[ntill_res_gcm_var_flux_dt, on = .(gridid = gridid,
                                                                          x = x,
                                                                          y = y,
                                                                          scenario = scenario,
                                                                          y_block = y_block,
                                                                          ssp = ssp,
                                                                          gcm = gcm,
                                                                          WB_NAME = WB_NAME,
                                                                          WB_REGION = WB_REGION,
                                                                          s_cr_grain = s_cr_grain)]
rm(ntill_res_gcm_flux_dt, ntill_res_gcm_var_flux_dt)
gc()
ccl_gcm_dt       = ccl_gcm_flux_dt[ccl_gcm_var_flux_dt, on = .(gridid = gridid,
                                                                             x = x,
                                                                             y = y,
                                                                             scenario = scenario,
                                                                             y_block = y_block,
                                                                             ssp = ssp,
                                                                             gcm = gcm,
                                                                             WB_NAME = WB_NAME,
                                                                             WB_REGION = WB_REGION,
                                                                             s_cr_grain = s_cr_grain)]
rm(ccl_gcm_flux_dt, ccl_gcm_var_flux_dt)
gc()
ccl_ntill_gcm_dt  = ccl_ntill_gcm_flux_dt[ccl_ntill_gcm_var_flux_dt, on = .(gridid = gridid,
                                                                  x = x,
                                                                  y = y,
                                                                  scenario = scenario,
                                                                  y_block = y_block,
                                                                  ssp = ssp,
                                                                  gcm = gcm,
                                                                  WB_NAME = WB_NAME,
                                                                  WB_REGION = WB_REGION,
                                                                  s_cr_grain = s_cr_grain)]
rm(ccl_ntill_gcm_flux_dt, ccl_ntill_gcm_var_flux_dt)
gc()
# TIDY CCL-NTILL RESPONSE DT
# filter to 2050
ccl_ntill_gcm_dt = ccl_ntill_gcm_dt[y_block <= 2050,]
# mean response for nfix, soil C:N, sfdcmp, slcmp
ccl_ntill_gcm_dt[, m_nfix       := mean(m_nfix), by = .(gridid, ssp, gcm)]
ccl_ntill_gcm_dt[, m_sC.N       := mean(m_sC.N), by = .(gridid, ssp, gcm)]
ccl_ntill_gcm_dt[, m_sfdcmp     := mean(m_sfdcmp), by = .(gridid, ssp, gcm)]
ccl_ntill_gcm_dt[, m_sldcmp     := mean(m_sldcmp), by = .(gridid, ssp, gcm)]
# sum cover crop biomass
ccl_ntill_gcm_dt[, s_cc_biomass := s_cc_shootC + s_cc_rootC]
# reduce columns
drop_cols = c('m_cr_grain', 'm_SOC', 'm_N2O', 'm_iN2O', 'm_dN2O', 'm_GHG', 's_N2O',
              's_iN2O', 's_dN2O', 's_GHG', 's_cc_rootC', 's_cc_shootC')
ccl_ntill_gcm_dt = ccl_ntill_gcm_dt[y_block == 2050, -..drop_cols]
# EXPECTED ADDITIVE RESPONSE DT
  # update column names
names = c('s_SOC', 's_cr_NPP', 's_cr_grain', 's_cr_residC','s_cc_rootC', 's_cc_shootC', 'm_sfdcmp',
          'm_sldcmp', 'm_sC.N', 'm_nfix', 's_annet', 's_gr_nit')
new_n = as.vector(outer('ccl_', colnames(ccl_gcm_dt[,..names]), paste0))
old_n = colnames(ccl_gcm_dt[, ..names])
setnames(ccl_gcm_dt, old = c(old_n), new = c(new_n))
cols  = c('gridid', 'x', 'y', 'y_block', 'gcm', 'ssp', 'WB_NAME', 'WB_REGION', new_n)
ccl_gcm_dt = ccl_gcm_dt[, ..cols]
  # join ntill_res and cover crop table
additive_dt = ntill_res_gcm_dt[ccl_gcm_dt, on = .(gridid = gridid,
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
drop_cols = c('m_cr_grain', 'm_SOC', 'm_N2O', 'm_iN2O', 'm_dN2O', 'm_GHG', 's_N2O',
              's_iN2O', 's_dN2O', 's_GHG', 's_cc_rootC', 's_cc_shootC', 'ccl_s_cc_rootC',
              'ccl_s_cc_shootC')
additive_dt[, scenario := 'additive-model']
additive_dt = additive_dt[, -..drop_cols]
  # difference dt | by hand
additive_dt[, s_SOC      := s_SOC + ccl_s_SOC]
additive_dt[, s_cr_grain := s_cr_grain + ccl_s_cr_grain]
additive_dt[, s_cr_NPP   := s_cr_NPP + ccl_s_cr_NPP]
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
ranger_dt     = ccl_ntill_gcm_dt[additive_dt, on = .(gridid = gridid,
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
ranger_dt = ranger_dt[, -..new_n]
# RANDOM FOREST

library(ranger)
library(fastshap)
library(shapviz)
pred_fun                 = function(object, newdata) {
  predict(object, data = newdata)$predictions
}
# diff between ccg and ccg_ntill
s             = sample(ranger_dt[,gridid],1000)
ranger_dt     = ranger_dt[gridid %in% s]
ranger_dt     = ranger_dt[!is.na(scenario)]
ranger_dt     = ranger_dt[complete.cases(ranger_dt)] # remove this!
drop_ranger   = c('gridid', 'x', 'y', 'scenario', 'y_block', 'gcm', 'ssp', 'WB_NAME', 'WB_REGION')
red_ranger_dt = ranger_dt[, -..drop_ranger]
library(caret)
correlationMatrix = cor(red_ranger_dt[,c(1,3:11)])
highlyCorrelated  = findCorrelation(correlationMatrix, cutoff=0.75, names = TRUE)
print(highlyCorrelated)
red_ranger_dt[, s_cr_NPP := NULL]

print('Running ranger.')
soc_ranger       = ranger(dependent.variable.name = 's_SOC', data = red_ranger_dt,
                          importance = 'permutation', keep.inbag = TRUE, seed = 1234)
soc_ranger

soc_features = red_ranger_dt[, -c('s_SOC')]
# SHAP COMPUTATION - N.B. TAKES A LONG TIME
print('Computing SHAP')
soc_SHAP       = explain(soc_ranger, X = soc_features, pred_wrapper = pred_fun, nsim = 10, adjust = TRUE,
                         shap_only = FALSE)
shv.soc_gl     = shapviz(soc_SHAP)
# initial visualization
soc_SHAP_gl_gg = sv_importance(shv.soc_gl)
soc_SHAP_be_gg = sv_importance(shv.soc_gl, kind = "beeswarm", show_numbers = FALSE, bee_width = 0.2)
# sv_dependence(shv.soc_gl, 'm_cr_NPP', color_var = NULL)
# sv_dependence(shv.soc_gl, 'm_sldcmp', color_var = NULL)
sv_dependence(shv.soc_gl, 's_gr_nit', color_var = NULL)
sv_dependence(shv.soc_gl, 'm_sldcmp', color_var = NULL)

# OLD #

# # REGION SE
#   # RESIDUE
# d_res_r_mSE = regional_decadal_mSE(d_res_r_mn)
# d_res_r_mSE[, scenario := 'res']
#   # NTILL
# d_ntill_r_mSE = regional_decadal_mSE(d_ntill_r_mn)
# d_ntill_r_mSE[, scenario := 'ntill']
#   # CCG
# d_ccg_r_mSE = regional_decadal_mSE(d_ccg_r_mn)
# d_ccg_r_mSE[, scenario := 'ccg']
#   # CCL
# d_ccl_r_mSE = regional_decadal_mSE(d_ccl_r_mn)
# d_ccl_r_mSE[, scenario := 'ccl']
#   # CCG-NTILL
# d_ccg_ntill_r_mSE = regional_decadal_mSE(d_ccg_ntill_r_mn)
# d_ccg_ntill_r_mSE[, scenario := 'ccg-ntill']
#   # CCL-NTILL
# d_ccl_ntill_r_mSE = regional_decadal_mSE(d_ccl_ntill_r_mn)
# d_ccl_ntill_r_mSE[, scenario := 'ccl-ntill']
# 
# d_r_mSE = rbind(d_res_r_mSE, d_ntill_r_mSE, d_ccg_r_mSE, d_ccl_r_mSE, d_ccg_ntill_r_mSE, d_ccl_ntill_r_mSE)
# d_r_mSE[, m_GHG_units := 'Mg CO2-eq ha-1 yr-1']
# d_r_mSE[, s_GHG_units := 'Mg CO2-eq ha-1']
# d_r_mSE[, m_biomass_units := 'kg C ha-1 yr-1']
# d_r_mSE[, s_biomass_units := 'Mg C ha-1']
# fwrite(d_r_mSE, paste(data.path, 'decadal-IPCC-region-crop-soil-responses-se.csv', sep = '/'))

# # COUNTRY
#   # RESIDUE
# d_res_c_mn  = country_decadal_mean(rel_gcm_res_flux_dt)
# d_res_c_mCI = country_decadal_mCI(d_res_c_mn)
#   # NTILL
# d_ntill_c_mn  = country_decadal_mean(rel_gcm_ntill_flux_dt)
# d_ntill_c_mCI = country_decadal_mCI(d_ntill_c_mn)
#   # CCG
# d_ccg_c_mn  = country_decadal_mean(rel_gcm_ccg_flux_dt)
# d_ccg_c_mCI = country_decadal_mCI(d_ccg_c_mn)
#   # CCL
# d_ccl_c_mn  = country_decadal_mean(rel_gcm_ccl_flux_dt)
# d_ccl_c_mCI = country_decadal_mCI(d_ccl_c_mn)
#   # CCG-NTILL
# d_ccg_ntill_c_mn  = country_decadal_mean(rel_gcm_ccg_ntill_flux_dt)
# d_ccg_ntill_c_mCI = country_decadal_mCI(d_ccg_ntill_c_mn)
#   # CCL-NTILL
# d_ccl_ntill_c_mn  = country_decadal_mean(rel_gcm_ccl_ntill_flux_dt)
# d_ccl_ntill_c_mCI = country_decadal_mCI(d_ccl_ntill_c_mn)
# 
# d_c_mCI = rbind(d_res_c_mCI, d_ntill_c_mCI, d_ccg_c_mCI, d_ccl_c_mCI, d_ccg_ntill_c_mCI, d_ccl_ntill_c_mCI)
# d_c_mCI[, m_GHG_units := 'Mg CO2-eq ha-1 yr-1']
# d_c_mCI[, s_GHG_units := 'Mg CO2-eq ha-1']
# d_c_mCI[, m_biomass_units := 'kg C ha-1 yr-1']
# d_c_mCI[, s_biomass_units := 'Mg C ha-1']
# fwrite(d_c_mCI, paste(data.path, 'decadal-country-crop-soil-responses.csv', sep = '/'))
#-------------------------------------------------------------------------------
# TEST CODE: EXPLANATORY VARIABLES
#-------------------------------------------------------------------------------
# COMPARE & PLOT CHANGES
var = c('m_gr_nit','s_gr_nit','m_annet','s_annet','m_nfix',
        'm_sfdcmp', 'm_sldcmp', 'm_sC.N', 'm_cr_NPP', 's_cr_NPP', 'm_cc_shootC', 
        's_cc_shootC','m_cc_rootC','s_cc_rootC')
ntill_var = rel_ntill_flux_dt[WB_REGION %in% 'Other', lapply(.SD, mean), .SDcols = var, by = .(gcm, ssp, y_block)]
ccg_var   = rel_ccg_flux_dt[WB_REGION %in% 'Other', lapply(.SD, mean), .SDcols = var, by = .(gcm, ssp, y_block)]
ccl_var   = rel_ccl_flux_dt[WB_REGION %in% 'Other', lapply(.SD, mean), .SDcols = var, by = .(gcm, ssp, y_block)]

ccg_ntill_var = rel_ccg_ntill_flux_dt[WB_REGION %in% 'Other', lapply(.SD, mean), .SDcols = var, by = .(gcm, ssp, y_block)]
ccl_ntill_var = rel_ccl_ntill_flux_dt[WB_REGION %in% 'Other', lapply(.SD, mean), .SDcols = var, by = .(gcm, ssp, y_block)]
#-------------------------------------------------------------------------------
# TEST CODE: EXPLANATORY VARIABLES | Supervised Clustering
#-------------------------------------------------------------------------------
var = c('m_SOC','m_annet','m_sfdcmp', 'm_sldcmp', 'm_sC.N', 'm_cr_NPP', 'm_cc_shootC')
ntill_var = rel_ntill_flux_dt[, lapply(.SD, mean), .SDcols = var, by = .(gcm, ssp, gridid)]
ntill_var[, m_cc_shootC := NULL]

ccl_var       = rel_ccl_flux_dt[, lapply(.SD, mean), .SDcols = var, by = .(gcm, ssp, gridid)]
ccl_ntill_var = rel_ccl_ntill_flux_dt[, lapply(.SD, mean), .SDcols = var, by = .(gcm, ssp, gridid)]

ccg_var       = rel_ccg_flux_dt[, lapply(.SD, mean), .SDcols = var, by = .(gcm, ssp, gridid)]
ccg_ntill_var = rel_ccg_ntill_flux_dt[, lapply(.SD, mean), .SDcols = var, by = .(gcm, ssp, gridid)]

# SOMETHING LIKE BELOW MIGHT WORK - NEED TO DISCUSS 

#soc RF and SHAP
library(ranger)
library(fastshap)
library(shapviz)
# diff between ccg and ccg_ntill
s = sample(ccl_ntill_var[,gridid],10000)
ccl_ntill_var_r = ccl_ntill_var[ssp %in% 'ssp370' & gridid %in% s, c(3:9)]
ccl_var_r       = ccl_var[ssp %in% 'ssp370' & gridid %in% s, c(3:9)]
ntill_var_r     = ntill_var[ssp %in% 'ssp370' & gridid %in% ccl_ntill_var_r$gridid, c(3:9)]

# diff_ccl = ccl_ntill_var_r[, c(2:7)] - ntill_var_r[, c(2:7)]
diff_ccl = ccl_ntill_var_r[, c(2:7)] - ccl_var_r[, c(2:7)]

print('Running ranger.')
soc_ranger       = ranger(dependent.variable.name = 'm_SOC', data = diff_ccl,
                          importance = 'permutation', keep.inbag = TRUE, seed = 1234)
soc_ranger

soc_features = diff_ccl[, -c('m_SOC')]
# SHAP COMPUTATION - N.B. TAKES A LONG TIME
print('Computing SHAP')
soc_SHAP       = explain(soc_ranger, X = soc_features, pred_wrapper = pred_fun, nsim = 10, adjust = TRUE,
                         shap_only = FALSE)
shv.soc_gl     = shapviz(soc_SHAP)
# initial visualization
soc_SHAP_gl_gg = sv_importance(shv.soc_gl)
soc_SHAP_be_gg = sv_importance(shv.soc_gl, kind = "beeswarm", show_numbers = FALSE, bee_width = 0.2)
# sv_dependence(shv.soc_gl, 'm_cr_NPP', color_var = NULL)
# sv_dependence(shv.soc_gl, 'm_sldcmp', color_var = NULL)
sv_dependence(shv.soc_gl, 'm_annet', color_var = NULL)
sv_dependence(shv.soc_gl, 'm_sldcmp', color_var = NULL)
