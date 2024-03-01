# filename:    data_imputed.R
# created:     22 February 2024
# updated:     23 February 2024
# author:      S.C. McClelland
# description: This file imputes data from DayCent simulations
#              to all global cropland soil N2O, CO2 over time.
#-------------------------------------------------------------------------------
library(data.table)
library(rstudioapi)
library(terra)
#-------------------------------------------------------------------------------
source('results_functions.R')
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
# RESIDUE
#-------------------------------------------------------------------------------
# ENSEMBLE
res_relative_flux = readRDS(paste(data.path, 'ensemble-relative-responses-weighted-mean-res.rds', sep = '/'))
res_ens_imp_dt    = impute_missing_ens(res_relative_flux, 
                                       rast(paste(lu.path, raster, sep = '/')))
saveRDS(res_ens_imp_dt, file = paste(data.path, 'imputed-ensemble-relative-responses-weighted-mean-res.rds', sep = '/'))

# GCM
res_relative_gcm_flux = readRDS(paste(data.path, 'gcm-relative-responses-weighted-mean-res.rds', sep = '/'))
res_gcm_imp_hist      = impute_missing_gcm(res_relative_gcm_flux,
                                           rast(paste(lu.path, raster, sep = '/')),
                                           'historical', 'historical')
res_gcm_imp_126       = impute_missing_gcm(res_relative_gcm_flux,
                                           rast(paste(lu.path, raster, sep = '/')),
                                           unique(res_relative_gcm_flux[ssp %in% 'ssp126', gcm]), 
                                           'ssp126')
res_gcm_imp_370       = impute_missing_gcm(res_relative_gcm_flux,
                                           rast(paste(lu.path, raster, sep = '/')),
                                           unique(res_relative_gcm_flux[ssp %in% 'ssp370', gcm]), 
                                           'ssp370')
res_gcm_imp_dt        = rbind(res_gcm_imp_hist, res_gcm_imp_126, res_gcm_imp_370)
gc()
saveRDS(res_gcm_imp_dt, file = paste(data.path, 'imputed-gcm-relative-responses-weighted-mean-res.rds', sep = '/'))
#-------------------------------------------------------------------------------
# NTILL
#-------------------------------------------------------------------------------
# ENSEMBLE
ntill_relative_flux = readRDS(paste(data.path, 'ensemble-relative-responses-weighted-mean-ntill.rds', sep = '/'))
ntill_ens_imp_dt    = impute_missing_ens(ntill_relative_flux, 
                                       rast(paste(lu.path, raster, sep = '/')))
saveRDS(ntill_ens_imp_dt, file = paste(data.path, 'imputed-ensemble-relative-responses-weighted-mean-ntill.rds', sep = '/'))

# GCM
ntill_relative_gcm_flux = readRDS(paste(data.path, 'gcm-relative-responses-weighted-mean-ntill.rds', sep = '/'))
ntill_gcm_imp_hist      = impute_missing_gcm(ntill_relative_gcm_flux,
                                           rast(paste(lu.path, raster, sep = '/')),
                                           'historical', 'historical')
ntill_gcm_imp_126       = impute_missing_gcm(ntill_relative_gcm_flux,
                                           rast(paste(lu.path, raster, sep = '/')),
                                           unique(ntill_relative_gcm_flux[ssp %in% 'ssp126', gcm]), 
                                           'ssp126')
ntill_gcm_imp_370       = impute_missing_gcm(ntill_relative_gcm_flux,
                                           rast(paste(lu.path, raster, sep = '/')),
                                           unique(ntill_relative_gcm_flux[ssp %in% 'ssp370', gcm]), 
                                           'ssp370')
ntill_gcm_imp_dt        = rbind(ntill_gcm_imp_hist, ntill_gcm_imp_126, ntill_gcm_imp_370)
gc()
saveRDS(ntill_gcm_imp_dt, file = paste(data.path, 'imputed-gcm-relative-responses-weighted-mean-ntill.rds', sep = '/'))
#-------------------------------------------------------------------------------
# NTILL-RES
#-------------------------------------------------------------------------------
# ENSEMBLE
ntill_res_relative_flux = readRDS(paste(data.path, 'ensemble-relative-responses-weighted-mean-ntill-res.rds', sep = '/'))
ntill_res_ens_imp_dt    = impute_missing_ens(ntill_res_relative_flux, 
                                         rast(paste(lu.path, raster, sep = '/')))
saveRDS(ntill_res_ens_imp_dt, file = paste(data.path, 'imputed-ensemble-relative-responses-weighted-mean-ntill-res.rds', sep = '/'))

# GCM
ntill_res_relative_gcm_flux = readRDS(paste(data.path, 'gcm-relative-responses-weighted-mean-ntill-res.rds', sep = '/'))
ntill_res_gcm_imp_hist      = impute_missing_gcm(ntill_res_relative_gcm_flux,
                                             rast(paste(lu.path, raster, sep = '/')),
                                             'historical', 'historical')
ntill_res_gcm_imp_126       = impute_missing_gcm(ntill_res_relative_gcm_flux,
                                             rast(paste(lu.path, raster, sep = '/')),
                                             unique(ntill_res_relative_gcm_flux[ssp %in% 'ssp126', gcm]), 
                                             'ssp126')
ntill_res_gcm_imp_370       = impute_missing_gcm(ntill_res_relative_gcm_flux,
                                             rast(paste(lu.path, raster, sep = '/')),
                                             unique(ntill_res_relative_gcm_flux[ssp %in% 'ssp370', gcm]), 
                                             'ssp370')
ntill_res_gcm_imp_dt        = rbind(ntill_res_gcm_imp_hist, ntill_res_gcm_imp_126, ntill_res_gcm_imp_370)
gc()
saveRDS(ntill_res_gcm_imp_dt, file = paste(data.path, 'imputed-gcm-relative-responses-weighted-mean-ntill-res.rds', sep = '/'))
#-------------------------------------------------------------------------------
# CCG
#-------------------------------------------------------------------------------
# ENSEMBLE
ccg_relative_flux = readRDS(paste(data.path, 'ensemble-relative-responses-weighted-mean-ccg.rds', sep = '/'))
ccg_ens_imp_dt    = impute_missing_ens(ccg_relative_flux, 
                                         rast(paste(lu.path, raster_cc, sep = '/')))
saveRDS(ccg_ens_imp_dt, file = paste(data.path, 'imputed-ensemble-relative-responses-weighted-mean-ccg.rds', sep = '/'))

# GCM
ccg_relative_gcm_flux = readRDS(paste(data.path, 'gcm-relative-responses-weighted-mean-ccg.rds', sep = '/'))
ccg_gcm_imp_hist      = impute_missing_gcm(ccg_relative_gcm_flux,
                                             rast(paste(lu.path, raster_cc, sep = '/')),
                                             'historical', 'historical')
ccg_gcm_imp_126       = impute_missing_gcm(ccg_relative_gcm_flux,
                                             rast(paste(lu.path, raster_cc, sep = '/')),
                                             unique(ccg_relative_gcm_flux[ssp %in% 'ssp126', gcm]), 
                                             'ssp126')
ccg_gcm_imp_370       = impute_missing_gcm(ccg_relative_gcm_flux,
                                             rast(paste(lu.path, raster_cc, sep = '/')),
                                             unique(ccg_relative_gcm_flux[ssp %in% 'ssp370', gcm]), 
                                             'ssp370')
ccg_gcm_imp_dt        = rbind(ccg_gcm_imp_hist, ccg_gcm_imp_126, ccg_gcm_imp_370)
gc()
saveRDS(ccg_gcm_imp_dt, file = paste(data.path, 'imputed-gcm-relative-responses-weighted-mean-ccg.rds', sep = '/'))
#-------------------------------------------------------------------------------
# CCG-RES
#-------------------------------------------------------------------------------
# ENSEMBLE
ccg_res_relative_flux = readRDS(paste(data.path, 'ensemble-relative-responses-weighted-mean-ccg-res.rds', sep = '/'))
ccg_res_ens_imp_dt    = impute_missing_ens(ccg_res_relative_flux, 
                                       rast(paste(lu.path, raster_cc, sep = '/')))
saveRDS(ccg_res_ens_imp_dt, file = paste(data.path, 'imputed-ensemble-relative-responses-weighted-mean-ccg-res.rds', sep = '/'))

# GCM
ccg_res_relative_gcm_flux = readRDS(paste(data.path, 'gcm-relative-responses-weighted-mean-ccg-res.rds', sep = '/'))
ccg_res_gcm_imp_hist      = impute_missing_gcm(ccg_res_relative_gcm_flux,
                                           rast(paste(lu.path, raster_cc, sep = '/')),
                                           'historical', 'historical')
ccg_res_gcm_imp_126       = impute_missing_gcm(ccg_res_relative_gcm_flux,
                                           rast(paste(lu.path, raster_cc, sep = '/')),
                                           unique(ccg_res_relative_gcm_flux[ssp %in% 'ssp126', gcm]), 
                                           'ssp126')
ccg_res_gcm_imp_370       = impute_missing_gcm(ccg_res_relative_gcm_flux,
                                           rast(paste(lu.path, raster_cc, sep = '/')),
                                           unique(ccg_res_relative_gcm_flux[ssp %in% 'ssp370', gcm]), 
                                           'ssp370')
ccg_res_gcm_imp_dt        = rbind(ccg_res_gcm_imp_hist, ccg_res_gcm_imp_126, ccg_res_gcm_imp_370)
gc()
saveRDS(ccg_res_gcm_imp_dt, file = paste(data.path, 'imputed-gcm-relative-responses-weighted-mean-ccg-res.rds', sep = '/'))
#-------------------------------------------------------------------------------
# CCL
#-------------------------------------------------------------------------------
# ENSEMBLE
ccl_relative_flux = readRDS(paste(data.path, 'ensemble-relative-responses-weighted-mean-ccl.rds', sep = '/'))
ccl_ens_imp_dt    = impute_missing_ens(ccl_relative_flux, 
                                       rast(paste(lu.path, raster_cc, sep = '/')))
saveRDS(ccl_ens_imp_dt, file = paste(data.path, 'imputed-ensemble-relative-responses-weighted-mean-ccl.rds', sep = '/'))

# GCM
ccl_relative_gcm_flux = readRDS(paste(data.path, 'gcm-relative-responses-weighted-mean-ccl.rds', sep = '/'))
ccl_gcm_imp_hist      = impute_missing_gcm(ccl_relative_gcm_flux,
                                           rast(paste(lu.path, raster_cc, sep = '/')),
                                           'historical', 'historical')
ccl_gcm_imp_126       = impute_missing_gcm(ccl_relative_gcm_flux,
                                           rast(paste(lu.path, raster_cc, sep = '/')),
                                           unique(ccl_relative_gcm_flux[ssp %in% 'ssp126', gcm]), 
                                           'ssp126')
ccl_gcm_imp_370       = impute_missing_gcm(ccl_relative_gcm_flux,
                                           rast(paste(lu.path, raster_cc, sep = '/')),
                                           unique(ccl_relative_gcm_flux[ssp %in% 'ssp370', gcm]), 
                                           'ssp370')
ccl_gcm_imp_dt        = rbind(ccl_gcm_imp_hist, ccl_gcm_imp_126, ccl_gcm_imp_370)
gc()
saveRDS(ccl_gcm_imp_dt, file = paste(data.path, 'imputed-gcm-relative-responses-weighted-mean-ccl.rds', sep = '/'))
#-------------------------------------------------------------------------------
# CCL-RES
#-------------------------------------------------------------------------------
# ENSEMBLE
ccl_res_relative_flux = readRDS(paste(data.path, 'ensemble-relative-responses-weighted-mean-ccl-res.rds', sep = '/'))
ccl_res_ens_imp_dt    = impute_missing_ens(ccl_res_relative_flux, 
                                           rast(paste(lu.path, raster_cc, sep = '/')))
saveRDS(ccl_res_ens_imp_dt, file = paste(data.path, 'imputed-ensemble-relative-responses-weighted-mean-ccl-res.rds', sep = '/'))

# GCM
ccl_res_relative_gcm_flux = readRDS(paste(data.path, 'gcm-relative-responses-weighted-mean-ccl-res.rds', sep = '/'))
ccl_res_gcm_imp_hist      = impute_missing_gcm(ccl_res_relative_gcm_flux,
                                               rast(paste(lu.path, raster_cc, sep = '/')),
                                               'historical', 'historical')
ccl_res_gcm_imp_126       = impute_missing_gcm(ccl_res_relative_gcm_flux,
                                               rast(paste(lu.path, raster_cc, sep = '/')),
                                               unique(ccl_res_relative_gcm_flux[ssp %in% 'ssp126', gcm]), 
                                               'ssp126')
ccl_res_gcm_imp_370       = impute_missing_gcm(ccl_res_relative_gcm_flux,
                                               rast(paste(lu.path, raster_cc, sep = '/')),
                                               unique(ccl_res_relative_gcm_flux[ssp %in% 'ssp370', gcm]), 
                                               'ssp370')
ccl_res_gcm_imp_dt        = rbind(ccl_res_gcm_imp_hist, ccl_res_gcm_imp_126, ccl_res_gcm_imp_370)
gc()
saveRDS(ccl_res_gcm_imp_dt, file = paste(data.path, 'imputed-gcm-relative-responses-weighted-mean-ccl-res.rds', sep = '/'))
#-------------------------------------------------------------------------------
# CCG-NTILL-RES
#-------------------------------------------------------------------------------
# ENSEMBLE
ccg_ntill_relative_flux = readRDS(paste(data.path, 'ensemble-relative-responses-weighted-mean-ccg-ntill.rds', sep = '/'))
ccg_ntill_ens_imp_dt    = impute_missing_ens(ccg_ntill_relative_flux, 
                                           rast(paste(lu.path, raster_cc, sep = '/')))
saveRDS(ccg_ntill_ens_imp_dt, file = paste(data.path, 'imputed-ensemble-relative-responses-weighted-mean-ccg-ntill.rds', sep = '/'))

# GCM
ccg_ntill_relative_gcm_flux = readRDS(paste(data.path, 'gcm-relative-responses-weighted-mean-ccg-ntill.rds', sep = '/'))
ccg_ntill_gcm_imp_hist      = impute_missing_gcm(ccg_ntill_relative_gcm_flux,
                                               rast(paste(lu.path, raster_cc, sep = '/')),
                                               'historical', 'historical')
ccg_ntill_gcm_imp_126       = impute_missing_gcm(ccg_ntill_relative_gcm_flux,
                                               rast(paste(lu.path, raster_cc, sep = '/')),
                                               unique(ccg_ntill_relative_gcm_flux[ssp %in% 'ssp126', gcm]), 
                                               'ssp126')
ccg_ntill_gcm_imp_370       = impute_missing_gcm(ccg_ntill_relative_gcm_flux,
                                               rast(paste(lu.path, raster_cc, sep = '/')),
                                               unique(ccg_ntill_relative_gcm_flux[ssp %in% 'ssp370', gcm]), 
                                               'ssp370')
ccg_ntill_gcm_imp_dt        = rbind(ccg_ntill_gcm_imp_hist, ccg_ntill_gcm_imp_126, ccg_ntill_gcm_imp_370)
gc()
saveRDS(ccg_ntill_gcm_imp_dt, file = paste(data.path, 'imputed-gcm-relative-responses-weighted-mean-ccg-ntill.rds', sep = '/'))
#-------------------------------------------------------------------------------
# CCL-NTILL-RES
#-------------------------------------------------------------------------------
# ENSEMBLE
ccl_ntill_relative_flux = readRDS(paste(data.path, 'ensemble-relative-responses-weighted-mean-ccl-ntill.rds', sep = '/'))
ccl_ntill_ens_imp_dt    = impute_missing_ens(ccl_ntill_relative_flux, 
                                             rast(paste(lu.path, raster_cc, sep = '/')))
saveRDS(ccl_ntill_ens_imp_dt, file = paste(data.path, 'imputed-ensemble-relative-responses-weighted-mean-ccl-ntill.rds', sep = '/'))

# GCM
ccl_ntill_relative_gcm_flux = readRDS(paste(data.path, 'gcm-relative-responses-weighted-mean-ccl-ntill.rds', sep = '/'))
ccl_ntill_gcm_imp_hist      = impute_missing_gcm(ccl_ntill_relative_gcm_flux,
                                                 rast(paste(lu.path, raster_cc, sep = '/')),
                                                 'historical', 'historical')
ccl_ntill_gcm_imp_126       = impute_missing_gcm(ccl_ntill_relative_gcm_flux,
                                                 rast(paste(lu.path, raster_cc, sep = '/')),
                                                 unique(ccl_ntill_relative_gcm_flux[ssp %in% 'ssp126', gcm]), 
                                                 'ssp126')
ccl_ntill_gcm_imp_370       = impute_missing_gcm(ccl_ntill_relative_gcm_flux,
                                                 rast(paste(lu.path, raster_cc, sep = '/')),
                                                 unique(ccl_ntill_relative_gcm_flux[ssp %in% 'ssp370', gcm]), 
                                                 'ssp370')
ccl_ntill_gcm_imp_dt        = rbind(ccl_ntill_gcm_imp_hist, ccl_ntill_gcm_imp_126, ccl_ntill_gcm_imp_370)
gc()
saveRDS(ccl_ntill_gcm_imp_dt, file = paste(data.path, 'imputed-gcm-relative-responses-weighted-mean-ccl-ntill.rds', sep = '/'))
#-------------------------------------------------------------------------------