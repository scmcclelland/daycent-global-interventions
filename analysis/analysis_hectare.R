# filename:    analysis_aggregated.R
# created:     06 March 2023
# updated:     15 August 2024
# author:      S.C. McClelland
# description: This file analyzes data from DayCent simulations
#              for global cropland soil N2O, CO2 over time for different 
#              interventions and SSPs.
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
fig.path  = paste(base.path, 'manuscript-figures', sep = '/')
#-------------------------------------------------------------------------------
# LOAD RELATIVE GCM DATA
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
# CONVERT Biomass Units: g C m-2 to kg ha-1 yr-1 or Mg ha-1
#-------------------------------------------------------------------------------
Mg_ha = 100L
kg_ha = 10L
C_gr  = 0.42 # Ma et al. 2018 | value for 'reproductive organs'
# RESIDUE
residue_gcm_flux_dt[, m_cr_grain := (m_cr_grain*kg_ha)/C_gr]
residue_gcm_flux_dt[, s_cr_grain := (s_cr_grain/Mg_ha)/C_gr]
# NTILL
ntill_gcm_flux_dt[, m_cr_grain := (m_cr_grain*kg_ha)/C_gr]
ntill_gcm_flux_dt[, s_cr_grain := (s_cr_grain/Mg_ha)/C_gr]
# NTILL-RES
ntill_res_gcm_flux_dt[, m_cr_grain := (m_cr_grain*kg_ha)/C_gr]
ntill_res_gcm_flux_dt[, s_cr_grain := (s_cr_grain/Mg_ha)/C_gr]
# CCG
ccg_gcm_flux_dt[, m_cr_grain := (m_cr_grain*kg_ha)/C_gr]
ccg_gcm_flux_dt[, s_cr_grain := (s_cr_grain/Mg_ha)/C_gr]
# CCG-RES
ccg_res_gcm_flux_dt[, m_cr_grain := (m_cr_grain*kg_ha)/C_gr]
ccg_res_gcm_flux_dt[, s_cr_grain := (s_cr_grain/Mg_ha)/C_gr]
# CCL
ccl_gcm_flux_dt[, m_cr_grain := (m_cr_grain*kg_ha)/C_gr]
ccl_gcm_flux_dt[, s_cr_grain := (s_cr_grain/Mg_ha)/C_gr]
# CCL-RES
ccl_res_gcm_flux_dt[, m_cr_grain := (m_cr_grain*kg_ha)/C_gr]
ccl_res_gcm_flux_dt[, s_cr_grain := (s_cr_grain/Mg_ha)/C_gr]
# CCG-NTILL-RES
ccg_ntill_gcm_flux_dt[, m_cr_grain := (m_cr_grain*kg_ha)/C_gr]
ccg_ntill_gcm_flux_dt[, s_cr_grain := (s_cr_grain/Mg_ha)/C_gr]
# CCL-NTILL-RES
ccl_ntill_gcm_flux_dt[, m_cr_grain := (m_cr_grain*kg_ha)/C_gr]
ccl_ntill_gcm_flux_dt[, s_cr_grain := (s_cr_grain/Mg_ha)/C_gr]
#-------------------------------------------------------------------------------
# MEAN GRID CELL RELATIVE RESPONSE BY DECADE | GLOBAL
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
# LOAD ABSOLUTE DATA
#-------------------------------------------------------------------------------
# CONV
conv_abs_flux_dt      = readRDS(paste(data.path, 'ensemble-absolute-responses-weighted-mean-conv.rds', sep = '/'))
#-------------------------------------------------------------------------------
# UPDATE REGION NAMES
#-------------------------------------------------------------------------------
# CONV
conv_abs_flux_dt     = conv_abs_flux_dt[crop_area_dt[, .(WB_NAME, IPCC_NAME)], 
                                           on = .(WB_NAME = WB_NAME) ]
conv_abs_flux_dt     = conv_abs_flux_dt[!is.na(gcm)]
#-------------------------------------------------------------------------------
# CONVERT Biomass Units
#-------------------------------------------------------------------------------
Mg_ha = 100L
kg_ha = 10L
C_gr  = 0.42 # Ma et al. 2018 | value for 'reproductive organs'
# CONV
conv_abs_flux_dt[, m_cr_grain := (m_cr_grain*kg_ha)/C_gr]
conv_abs_flux_dt[, s_cr_grain := (s_cr_grain/Mg_ha)/C_gr]
#-------------------------------------------------------------------------------
# MEAN GRID CELL RELATIVE RESPONSE BY DECADE | GLOBAL
#-------------------------------------------------------------------------------
  # GLOBAL
# CONV
d_conv_gl_mn  = global_decadal_mean(conv_abs_flux_dt)
d_conv_gl_mCI = global_decadal_mCI(d_conv_gl_mn)
d_conv_gl_mn[,  scenario := 'conv']
d_conv_gl_mCI[, scenario := 'conv']

# global response mean, CI
d_conv_gl_mCI[, IPCC_NAME := 'GLOBE']
d_conv_gl_mCI[, m_GHG_units := 'Mg CO2-eq ha-1 yr-1']
d_conv_gl_mCI[, s_GHG_units := 'Mg CO2-eq ha-1']
d_conv_gl_mCI[, m_biomass_units := 'kg C ha-1 yr-1']
d_conv_gl_mCI[, s_biomass_units := 'Mg C ha-1']
setcolorder(d_conv_gl_mCI, c('IPCC_NAME','ssp', 'scenario', 'y_block'))
fwrite(d_conv_gl_mCI, paste(data.path, 'ensemble-decadal-global-crop-soil-responses-reference-case.csv', sep = '/'))

  # REGION CI
# CONV
d_conv_r_mn  = regional_decadal_mean(conv_abs_flux_dt)
d_conv_r_mCI = regional_decadal_mCI(d_conv_r_mn)
d_conv_r_mn[,  scenario := 'conv']
d_conv_r_mCI[, scenario := 'conv']

# regional response mean, CI
d_conv_r_mCI[, m_GHG_units := 'Mg CO2-eq ha-1 yr-1']
d_conv_r_mCI[, s_GHG_units := 'Mg CO2-eq ha-1']
d_conv_r_mCI[, m_biomass_units := 'kg C ha-1 yr-1']
d_conv_r_mCI[, s_biomass_units := 'Mg C ha-1']
setcolorder(d_conv_r_mCI, c('IPCC_NAME','ssp', 'scenario', 'y_block'))
fwrite(d_conv_r_mCI, paste(data.path, 'ensemble-decadal-IPCC-region-crop-soil-responses-reference-case.csv', sep = '/'))