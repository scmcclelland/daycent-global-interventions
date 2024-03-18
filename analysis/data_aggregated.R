# filename:    data_aggregated.R
# created:     06 March 2023
# updated:     13 March 2024
# author:      S.C. McClelland
# description: This file generates data from DayCent simulations
#              for global cropland soil N2O, CO2 over time.
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
# ALL Cropland Area DT
#-------------------------------------------------------------------------------
crop_area_t_dt = country_area(lu.path, shp, raster) # 525,346,032

# All Cropland Area
country_names  = crop_area_t_dt[crop_sum_2015 > 0, WB_NAME]
country_names  = sort(country_names)

# IPCC Region Names (AR6 & Roe et al. 2021)
# Africa and Middle East
AME   = c('Congo, Democratic Republic of', 'Nigeria', 'Tanzania', 'South Africa', 'Congo, Rep. of', 'Zambia',
          'Angola', 'Cameroon', 'Ethiopia', 'Mozambique', 'Iran, Islamic Republic of', 'Uganda',
          'Central African Republic', 'Gabon', 'Sudan', "Côte d'Ivoire", 'Kenya', 'Egypt, Arab Republic of',
          'Ghana', 'Zimbabwe', 'Mali', 'Namibia', 'South Sudan', 'Chad', 'Morocco', 'Botswana', 'Burkina Faso',
          'Niger', 'Guinea', 'Algeria', 'Liberia', 'Malawi', 'Senegal', 'Somalia', 'Saudi Arabia', 'Benin', 
          'Sierra Leone', 'Iraq', 'Rwanda', 'Eritrea', 'eSwatini', 'Benin', 'Burundi', 'Djibouti', 'Equatorial Guinea',
          'Madagascar', 'Mauritania', 'Tunisia', 'Syrian Arab Republic', 'Lebanon', 'Jordan', 'Libya', 'Israel', 
          'West Bank and Gaza', 'Kuwait', 'Oman', 'Qatar', 'United Arab Emirates', 'Yemen, Republic of', 'Cabo Verde',
          'Guinea-Bissau', 'Togo', 'Comoros', 'Mauritius', 'Lesotho')
ADP   = c('China', 'Indonesia', 'India', 'Myanmar', 'Vietnam', 'Malaysia', 'Thailand', 'Pakistan', 'Papua New Guinea',
          'Philippines', 'Bangladesh', 'Cambodia', "Lao People's Democratic Republic", 'Mongolia', 'Korea, Republic of',
          'Afghanistan', 'Nepa', 'Sri Lanka', "Korea, Democratic People's Republic of", 'Solomon Islands', 'Bhutan',
          'Timor-Leste', 'Fiji', 'Nepal', 'Hong Kong (SAR, China)', 'Brunei Darussalam', 'Samoa', 'Vanuatu', 'Tonga')
DEV   = c('United States of America', 'Canada', 'Austria', 'Belgium', 'Bulgaria', 'Croatia', 'Czech Republic', 'Denmark',
          'Estonia', 'Finland', 'France', 'Germany', 'Greece', 'Hungary', 'Ireland', 'Italy', 'Latvia', 'Lithuania', 'Luxembourg',
          'Netherlands', 'Poland', 'Portugal', 'Romania', 'Slovak Republic', 'Slovenia','Spain', 'Sweden', 'United Kingdom', 'Australia', 'Ukraine',
          'Japan', 'Turkey', 'New Zealand', 'Norway', 'Iceland', 'Greenland (Den.)', 'Faroe Islands (Den.)', 'Switzerland', 'Saint-Pierre-et-Miquelon (Fr.)',
          'Cyprus', 'Puerto Rico (US)', 'American Samoa (US)', 'Saint Helena, Ascension and Tristan da Cunha (UK)', 'New Caledonia (Fr.)',
          'French Southern and Antarctic Lands (Fr.)', 'Falkland Islands (UK)/Islas Malvinas', 'South Georgia and South Sandwich Islands (UK)')
EEWCA = c('Russian Federation', 'Kazakhstan', 'Belarus', 'Uzbekistan', 'Turkmenistan', 'Kyrgyz Republic', 'Azerbaijan',
          'Moldova', 'Tajikistan', 'Armenia', 'Serbia', 'Bosnia and Herzegovina', 'Georgia', 'Montenegro', 'Kosovo', 'Albania',
          'North Macedonia')
LAC   = c('Brazil', 'Colombia', 'Mexico', 'Argentina', 'Bolivia', 'Peru', 'Venezuela', 'Paraguay', 'Ecuador', 'Chile', 'Guyana', 'Suriname',
          'Cuba', 'Uruguay', 'Honduras', 'Nicaragua', 'Guatemala', 'Guyana', 'Costa Rica', 'Panama', 'Dominican Republic', 'El Salvador', 'Belize',
          'Bahamas, The', 'Haiti', 'Turks and Caicos Islands (UK)', 'Jamaica', 'Venezuela, Republica Bolivariana de', 'Trinidad and Tobago')
crop_area_t_dt[WB_NAME %in% AME, IPCC_NAME := 'AME']
crop_area_t_dt[WB_NAME %in% ADP, IPCC_NAME := 'ADP']
crop_area_t_dt[WB_NAME %in% DEV, IPCC_NAME := 'DEV']
crop_area_t_dt[WB_NAME %in% EEWCA, IPCC_NAME := 'EEWCA']
crop_area_t_dt[WB_NAME %in% LAC, IPCC_NAME := 'LAC']
saveRDS(crop_area_t_dt, file = paste(data.path, 'crop-area-country-ipcc-region.rds', sep = '/'))
#-------------------------------------------------------------------------------
# COVER CROP ONLY Cropland Area DT
#-------------------------------------------------------------------------------
crop_area_tcc_dt = country_area(lu.path, shp, raster_cc) # 405,585,947
country_names  = crop_area_tcc_dt[crop_sum_2015 > 0, WB_NAME]
country_names  = sort(country_names)

# IPCC Region Names (AR6 & Roe et al. 2021)

# check for Na and assign, esp for ME
# Africa and Middle East
AME   = c('Congo, Democratic Republic of', 'Nigeria', 'Tanzania', 'South Africa', 'Congo, Rep. of', 'Zambia',
          'Angola', 'Cameroon', 'Ethiopia', 'Mozambique', 'Iran, Islamic Republic of', 'Uganda',
          'Central African Republic', 'Gabon', 'Sudan', "Côte d'Ivoire", 'Kenya', 'Egypt, Arab Republic of',
          'Ghana', 'Zimbabwe', 'Mali', 'Namibia', 'South Sudan', 'Chad', 'Morocco', 'Botswana', 'Burkina Faso',
          'Niger', 'Guinea', 'Algeria', 'Liberia', 'Malawi', 'Senegal', 'Somalia', 'Saudi Arabia', 'Benin', 
          'Sierra Leone', 'Iraq', 'Rwanda', 'Eritrea', 'eSwatini', 'Benin', 'Burundi', 'Djibouti', 'Equatorial Guinea',
          'Madagascar', 'Mauritania', 'Tunisia', 'Syrian Arab Republic', 'Lebanon', 'Jordan', 'Libya', 'Israel', 
          'West Bank and Gaza', 'Kuwait', 'Oman', 'Qatar', 'United Arab Emirates', 'Yemen, Republic of', 'Cabo Verde',
          'Guinea-Bissau', 'Togo', 'Comoros', 'Mauritius', 'Lesotho')
ADP   = c('China', 'Indonesia', 'India', 'Myanmar', 'Vietnam', 'Malaysia', 'Thailand', 'Pakistan', 'Papua New Guinea',
          'Philippines', 'Bangladesh', 'Cambodia', "Lao People's Democratic Republic", 'Mongolia', 'Korea, Republic of',
          'Afghanistan', 'Nepa', 'Sri Lanka', "Korea, Democratic People's Republic of", 'Solomon Islands', 'Bhutan',
          'Timor-Leste', 'Fiji', 'Nepal', 'Hong Kong (SAR, China)', 'Brunei Darussalam', 'Samoa', 'Vanuatu', 'Tonga', 'Kiribati')
DEV   = c('United States of America', 'Canada', 'Austria', 'Belgium', 'Bulgaria', 'Croatia', 'Czech Republic', 'Denmark',
          'Estonia', 'Finland', 'France', 'Germany', 'Greece', 'Hungary', 'Ireland', 'Italy', 'Latvia', 'Lithuania', 'Luxembourg',
          'Netherlands', 'Poland', 'Portugal', 'Romania', 'Slovak Republic', 'Slovenia','Spain', 'Sweden', 'United Kingdom', 'Australia', 'Ukraine',
          'Japan', 'Turkey', 'New Zealand', 'Norway', 'Iceland', 'Greenland (Den.)', 'Faroe Islands (Den.)', 'Switzerland', 'Saint-Pierre-et-Miquelon (Fr.)',
          'Cyprus', 'Puerto Rico (US)', 'American Samoa (US)', 'Saint Helena, Ascension and Tristan da Cunha (UK)', 'New Caledonia (Fr.)',
          'French Southern and Antarctic Lands (Fr.)', 'Falkland Islands (UK)/Islas Malvinas', 'South Georgia and South Sandwich Islands (UK)')
EEWCA = c('Russian Federation', 'Kazakhstan', 'Belarus', 'Uzbekistan', 'Turkmenistan', 'Kyrgyz Republic', 'Azerbaijan',
          'Moldova', 'Tajikistan', 'Armenia', 'Serbia', 'Bosnia and Herzegovina', 'Georgia', 'Montenegro', 'Kosovo', 'Albania',
          'North Macedonia')
LAC   = c('Brazil', 'Colombia', 'Mexico', 'Argentina', 'Bolivia', 'Peru', 'Venezuela', 'Paraguay', 'Ecuador', 'Chile', 'Guyana', 'Suriname',
          'Cuba', 'Uruguay', 'Honduras', 'Nicaragua', 'Guatemala', 'Guyana', 'Costa Rica', 'Panama', 'Dominican Republic', 'El Salvador', 'Belize',
          'Bahamas, The', 'Haiti', 'Turks and Caicos Islands (UK)', 'Jamaica', 'Venezuela, Republica Bolivariana de', 'Trinidad and Tobago')
crop_area_tcc_dt[WB_NAME %in% AME, IPCC_NAME := 'AME']
crop_area_tcc_dt[WB_NAME %in% ADP, IPCC_NAME := 'ADP']
crop_area_tcc_dt[WB_NAME %in% DEV, IPCC_NAME := 'DEV']
crop_area_tcc_dt[WB_NAME %in% EEWCA, IPCC_NAME := 'EEWCA']
crop_area_tcc_dt[WB_NAME %in% LAC, IPCC_NAME := 'LAC']
saveRDS(crop_area_tcc_dt, file = paste(data.path, 'cover-crop-crop-area-country-ipcc-region.rds', sep = '/'))

#-------------------------------------------------------------------------------
# BAD RUN YR BY CROP AND IRR AND SCENARIO
#-------------------------------------------------------------------------------
# REMOVE INCOMPLETE RUN YEARS AND SAVE TO FILTER GRIDID BY CROP & IRR & SCENARIO
# ---> CHECKS
# 1.   Simulations that did not run for 85 yrs, by crop + irr + scenario
# N.B. From this data it is not possible to determine reason for all failed runs
#      Likely due to: (a) bad weather files and (b) bad site files
store.path1 = '/mnt2/Documents/soilfutures-back-up/Soil-Futures-Data/'
store.path2 = '/mnt2/Documents/soilfutures-back-up/Soil-Futures-Data-2/analysis-files'
# CONV
bad_run_year('conv', base.path, paste0(store.path1, '/oct-23/analysis-files'), store.path2)
# RESIDUE
bad_run_year('res', base.path, paste0(store.path1, '/oct-23/analysis-files'), store.path2)
# NTILL
bad_run_year('ntill', base.path, paste0(store.path1, '/oct-23/analysis-files'), store.path2)
# CCG
bad_run_year('ccg', base.path, paste0(store.path1, '/oct-23/analysis-files'), store.path2)
# CCL
bad_run_year('ccl', base.path, paste0(store.path1, '/oct-23/analysis-files'), store.path2)
# NTILL-RES
bad_run_year('ntill-res', base.path, paste0(store.path1, '/feb-24/analysis-files'), store.path2)
# CCG-RES
bad_run_year('ccg-res', base.path, paste0(store.path1, '/feb-24/analysis-files'), store.path2)
# CCL-RES
bad_run_year('ccl-res', base.path, paste0(store.path1, '/feb-24/analysis-files'), store.path2)
# CCG-NTILL-RES
bad_run_year('ccg-ntill', base.path, paste0(store.path1, '/oct-23/analysis-files'), store.path2)
# CCL-NTILL-RES
bad_run_year('ccl-ntill', base.path, paste0(store.path1, '/oct-23/analysis-files'), store.path2)
#-------------------------------------------------------------------------------
# RESIDUE
#-------------------------------------------------------------------------------
conv_bad_gr   = fread(paste(base.path, 'data/bad-run-years-for-conv-gridid.csv', sep = '/'))
res_bad_gr    = fread(paste(base.path, 'data/bad-run-years-for-res-gridid.csv', sep = '/'))
res_bad_run_f = bad_run_match(res_bad_gr, conv_bad_gr)
# ENSEMBLE
res_relative_flux         = cell_crop_area_wmean(base.path, 'res', lu.path, raster,
                                         res_bad_run_f)
saveRDS(res_relative_flux, file = paste(data.path, 'ensemble-relative-responses-weighted-mean-res.rds', sep = '/'))
# GCM
res_relative_gcm_flux     = cell_crop_area_gcm_wmean(base.path, 'res', lu.path, raster,
                                                 res_bad_run_f)
saveRDS(res_relative_gcm_flux, file = paste(data.path, 'gcm-relative-responses-weighted-mean-res.rds', sep = '/'))
# GCM PRICE
res_relative_gcm_price    = cell_crop_area_gcm_price(base.path, 'res', lu.path, raster,
                                                         res_bad_run_f, 'OECD-Grain-Outlook-2023-2032.csv')
saveRDS(res_relative_gcm_price, file = paste(data.path, 'gcm-relative-grain-price-weighted-mean-res.rds', sep = '/'))
# ABSOLUTE
res_absolute_flux         = cell_crop_area_wmean_abs(base.path, 'res', lu.path, raster,
                                                 res_bad_run_f)
saveRDS(res_absolute_flux, file = paste(data.path, 'ensemble-absolute-responses-weighted-mean-res.rds', sep = '/'))
#-------------------------------------------------------------------------------
# NTILL
#-------------------------------------------------------------------------------
conv_bad_gr     = fread(paste(base.path, 'data/bad-run-years-for-conv-gridid.csv', sep = '/'))
ntill_bad_gr    = fread(paste(base.path, 'data/bad-run-years-for-ntill-gridid.csv', sep = '/'))
ntill_bad_run_f = bad_run_match(ntill_bad_gr, conv_bad_gr)
# ENSEMBLE
ntill_relative_flux     = cell_crop_area_wmean(base.path, 'ntill', lu.path, raster,
                                             ntill_bad_run_f)
saveRDS(ntill_relative_flux, file = paste(data.path, 'ensemble-relative-responses-weighted-mean-ntill.rds', sep = '/'))
# GCM
ntill_relative_gcm_flux = cell_crop_area_gcm_wmean(base.path, 'ntill', lu.path, raster,
                                                 ntill_bad_run_f)
saveRDS(ntill_relative_gcm_flux, file = paste(data.path, 'gcm-relative-responses-weighted-mean-ntill.rds', sep = '/'))
# GCM PRICE
ntill_relative_gcm_price    = cell_crop_area_gcm_price(base.path, 'ntill', lu.path, raster,
                                                     ntill_bad_run_f, 'OECD-Grain-Outlook-2023-2032.csv')
saveRDS(ntill_relative_gcm_price, file = paste(data.path, 'gcm-relative-grain-price-weighted-mean-ntill.rds', sep = '/'))
# ABSOLUTE
ntill_absolute_flux     = cell_crop_area_wmean_abs(base.path, 'ntill', lu.path, raster,
                                               ntill_bad_run_f)
saveRDS(ntill_absolute_flux, file = paste(data.path, 'ensemble-absolute-responses-weighted-mean-ntill.rds', sep = '/'))
#-------------------------------------------------------------------------------
# NTILL-RES
#-------------------------------------------------------------------------------
conv_bad_gr         = fread(paste(base.path, 'data/bad-run-years-for-conv-gridid.csv', sep = '/'))
ntill_res_bad_gr    = fread(paste(base.path, 'data/bad-run-years-for-ntill-res-gridid.csv', sep = '/'))
ntill_res_bad_run_f = bad_run_match(ntill_res_bad_gr, conv_bad_gr)
# ENSEMBLE
ntill_res_relative_flux         = cell_crop_area_wmean(base.path, 'ntill-res', lu.path, raster,
                                               ntill_res_bad_run_f)
saveRDS(ntill_res_relative_flux, file = paste(data.path, 'ensemble-relative-responses-weighted-mean-ntill-res.rds', sep = '/'))
# GCM
ntill_res_relative_gcm_flux     = cell_crop_area_gcm_wmean(base.path, 'ntill-res', lu.path, raster,
                                                   ntill_res_bad_run_f)
saveRDS(ntill_res_relative_gcm_flux, file = paste(data.path, 'gcm-relative-responses-weighted-mean-ntill-res.rds', sep = '/'))
# GCM PRICE
ntill_res_relative_gcm_price    = cell_crop_area_gcm_price(base.path, 'ntill-res', lu.path, raster,
                                                       ntill_res_bad_run_f, 'OECD-Grain-Outlook-2023-2032.csv')
saveRDS(ntill_res_relative_gcm_price, file = paste(data.path, 'gcm-relative-grain-price-weighted-mean-ntill-res.rds', sep = '/'))
# ABSOLUTE
ntill_res_absolute_flux         = cell_crop_area_wmean_abs(base.path, 'ntill-res', lu.path, raster,
                                                   ntill_res_bad_run_f)
saveRDS(ntill_res_absolute_flux, file = paste(data.path, 'ensemble-absolute-responses-weighted-mean-ntill-res.rds', sep = '/'))
#-------------------------------------------------------------------------------
# CCG
#-------------------------------------------------------------------------------
conv_bad_gr   = fread(paste(base.path, 'data/bad-run-years-for-conv-gridid.csv', sep = '/'))
ccg_bad_gr    = fread(paste(base.path, 'data/bad-run-years-for-ccg-gridid.csv', sep = '/'))
ccg_bad_run_f = bad_run_match(ccg_bad_gr, conv_bad_gr)
# ENSEMBLE
ccg_relative_flux         = cell_crop_area_wmean(base.path, 'ccg', lu.path, raster_cc,
                                               ccg_bad_run_f)
saveRDS(ccg_relative_flux, file = paste(data.path, 'ensemble-relative-responses-weighted-mean-ccg.rds', sep = '/'))
# GCM
ccg_relative_gcm_flux     = cell_crop_area_gcm_wmean(base.path, 'ccg', lu.path, raster_cc,
                                                   ccg_bad_run_f)
saveRDS(ccg_relative_gcm_flux, file = paste(data.path, 'gcm-relative-responses-weighted-mean-ccg.rds', sep = '/'))
# GCM PRICE
ccg_relative_gcm_price    = cell_crop_area_gcm_price(base.path, 'ccg', lu.path, raster_cc,
                                                 ccg_bad_run_f, 'OECD-Grain-Outlook-2023-2032.csv')
saveRDS(ccg_relative_gcm_price, file = paste(data.path, 'gcm-relative-grain-price-weighted-mean-ccg.rds', sep = '/'))
# ABSOLUTE
ccg_absolute_flux         = cell_crop_area_wmean_abs(base.path, 'ccg', lu.path, raster_cc,
                                                 ccg_bad_run_f)
saveRDS(ccg_absolute_flux, file = paste(data.path, 'ensemble-absolute-responses-weighted-mean-ccg.rds', sep = '/'))
#-------------------------------------------------------------------------------
# CCG-RES
#-------------------------------------------------------------------------------
conv_bad_gr       = fread(paste(base.path, 'data/bad-run-years-for-conv-gridid.csv', sep = '/'))
ccg_res_bad_gr    = fread(paste(base.path, 'data/bad-run-years-for-ccg-res-gridid.csv', sep = '/'))
ccg_res_bad_run_f = bad_run_match(ccg_res_bad_gr, conv_bad_gr)
# ENSEMBLE
ccg_res_relative_flux         = cell_crop_area_wmean(base.path, 'ccg-res', lu.path, raster_cc,
                                                   ccg_res_bad_run_f)
saveRDS(ccg_res_relative_flux, file = paste(data.path, 'ensemble-relative-responses-weighted-mean-ccg-res.rds', sep = '/'))
# GCM
ccg_res_relative_gcm_flux = cell_crop_area_gcm_wmean(base.path, 'ccg-res', lu.path, raster_cc,
                                                       ccg_res_bad_run_f)
saveRDS(ccg_res_relative_gcm_flux, file = paste(data.path, 'gcm-relative-responses-weighted-mean-ccg-res.rds', sep = '/'))
# GCM PRICE
ccg_res_relative_gcm_price    = cell_crop_area_gcm_price(base.path, 'ccg-res', lu.path, raster_cc,
                                                     ccg_res_bad_run_f, 'OECD-Grain-Outlook-2023-2032.csv')
saveRDS(ccg_res_relative_gcm_price, file = paste(data.path, 'gcm-relative-grain-price-weighted-mean-ccg-res.rds', sep = '/'))
# ABSOLUTE
ccg_res_absolute_flux         = cell_crop_area_wmean_abs(base.path, 'ccg-res', lu.path, raster_cc,
                                                     ccg_res_bad_run_f)
saveRDS(ccg_res_absolute_flux, file = paste(data.path, 'ensemble-absolute-responses-weighted-mean-ccg-res.rds', sep = '/'))
#-------------------------------------------------------------------------------
# CCL
#-------------------------------------------------------------------------------
conv_bad_gr   = fread(paste(base.path, 'data/bad-run-years-for-conv-gridid.csv', sep = '/'))
ccl_bad_gr    = fread(paste(base.path, 'data/bad-run-years-for-ccl-gridid.csv', sep = '/'))
ccl_bad_run_f = bad_run_match(ccl_bad_gr, conv_bad_gr)
# ENSEMBLE
ccl_relative_flux     = cell_crop_area_wmean(base.path, 'ccl', lu.path, raster_cc,
                                             ccl_bad_run_f)
saveRDS(ccl_relative_flux, file = paste(data.path, 'ensemble-relative-responses-weighted-mean-ccl.rds', sep = '/'))
# GCM
ccl_relative_gcm_flux = cell_crop_area_gcm_wmean(base.path, 'ccl', lu.path, raster_cc,
                                                 ccl_bad_run_f)
saveRDS(ccl_relative_gcm_flux, file = paste(data.path, 'gcm-relative-responses-weighted-mean-ccl.rds', sep = '/'))
# GCM PRICE
ccl_relative_gcm_price    = cell_crop_area_gcm_price(base.path, 'ccl', lu.path, raster_cc,
                                                     ccl_bad_run_f, 'OECD-Grain-Outlook-2023-2032.csv')
saveRDS(ccl_relative_gcm_price, file = paste(data.path, 'gcm-relative-grain-price-weighted-mean-ccl.rds', sep = '/'))
# ABSOLUTE
ccl_absolute_flux     = cell_crop_area_wmean_abs(base.path, 'ccl', lu.path, raster_cc,
                                             ccl_bad_run_f)
saveRDS(ccl_absolute_flux, file = paste(data.path, 'ensemble-absolute-responses-weighted-mean-ccl.rds', sep = '/'))
#-------------------------------------------------------------------------------
# CCL-RES
#-------------------------------------------------------------------------------
conv_bad_gr         = fread(paste(base.path, 'data/bad-run-years-for-conv-gridid.csv', sep = '/'))
ccl_res_bad_gr    = fread(paste(base.path, 'data/bad-run-years-for-ccl-res-gridid.csv', sep = '/'))
ccl_res_bad_run_f = bad_run_match(ccl_res_bad_gr, conv_bad_gr)
# ENSEMBLE
ccl_res_relative_flux     = cell_crop_area_wmean(base.path, 'ccl-res', lu.path, raster_cc,
                                                 ccl_res_bad_run_f)
saveRDS(ccl_res_relative_flux, file = paste(data.path, 'ensemble-relative-responses-weighted-mean-ccl-res.rds', sep = '/'))
# GCM
ccl_res_relative_gcm_flux = cell_crop_area_gcm_wmean(base.path, 'ccl-res', lu.path, raster_cc,
                                                     ccl_res_bad_run_f)
saveRDS(ccl_res_relative_gcm_flux, file = paste(data.path, 'gcm-relative-responses-weighted-mean-ccl-res.rds', sep = '/'))
# GCM PRICE
ccl_res_relative_gcm_price    = cell_crop_area_gcm_price(base.path, 'ccl-res', lu.path, raster_cc,
                                                         ccl_res_bad_run_f, 'OECD-Grain-Outlook-2023-2032.csv')
saveRDS(ccl_res_relative_gcm_price, file = paste(data.path, 'gcm-relative-grain-price-weighted-mean-ccl-res.rds', sep = '/'))
# ABSOLUTE
ccl_res_absolute_flux     = cell_crop_area_wmean_abs(base.path, 'ccl-res', lu.path, raster_cc,
                                                 ccl_res_bad_run_f)
saveRDS(ccl_res_absolute_flux, file = paste(data.path, 'ensemble-absolute-responses-weighted-mean-ccl-res.rds', sep = '/'))
#-------------------------------------------------------------------------------
# CCG-NTILL-RES
#-------------------------------------------------------------------------------
conv_bad_gr         = fread(paste(base.path, 'data/bad-run-years-for-conv-gridid.csv', sep = '/'))
ccg_ntill_bad_gr    = fread(paste(base.path, 'data/bad-run-years-for-ccg-ntill-gridid.csv', sep = '/'))
ccg_ntill_bad_run_f = bad_run_match(ccg_ntill_bad_gr, conv_bad_gr)
# ENSEMBLE
ccg_ntill_relative_flux     = cell_crop_area_wmean(base.path, 'ccg-ntill', lu.path, raster_cc,
                                             ccg_ntill_bad_run_f)
saveRDS(ccg_ntill_relative_flux, file = paste(data.path, 'ensemble-relative-responses-weighted-mean-ccg-ntill.rds', sep = '/'))
# GCM
ccg_ntill_relative_gcm_flux = cell_crop_area_gcm_wmean(base.path, 'ccg-ntill', lu.path, raster_cc,
                                                 ccg_ntill_bad_run_f)
saveRDS(ccg_ntill_relative_gcm_flux, file = paste(data.path, 'gcm-relative-responses-weighted-mean-ccg-ntill.rds', sep = '/'))
# GCM PRICE
ccg_ntill_relative_gcm_price    = cell_crop_area_gcm_price(base.path, 'ccg-ntill', lu.path, raster_cc,
                                                         ccg_ntill_bad_run_f, 'OECD-Grain-Outlook-2023-2032.csv')
saveRDS(ccg_ntill_relative_gcm_price, file = paste(data.path, 'gcm-relative-grain-price-weighted-mean-ccg-ntill.rds', sep = '/'))
# ABSOLUTE
ccg_ntill_absolute_flux     = cell_crop_area_wmean_abs(base.path, 'ccg-ntill', lu.path, raster_cc,
                                                   ccg_ntill_bad_run_f)
saveRDS(ccg_ntill_absolute_flux, file = paste(data.path, 'ensemble-absolute-responses-weighted-mean-ccg-ntill.rds', sep = '/'))
#-------------------------------------------------------------------------------
# CCL-NTILL-RES
#-------------------------------------------------------------------------------
conv_bad_gr         = fread(paste(base.path, 'data/bad-run-years-for-conv-gridid.csv', sep = '/'))
ccl_ntill_bad_gr    = fread(paste(base.path, 'data/bad-run-years-for-ccl-ntill-gridid.csv', sep = '/'))
ccl_ntill_bad_run_f = bad_run_match(ccl_ntill_bad_gr, conv_bad_gr)
# ENSEMBLE
ccl_ntill_relative_flux     = cell_crop_area_wmean(base.path, 'ccl-ntill', lu.path, raster_cc,
                                                   ccl_ntill_bad_run_f)
saveRDS(ccl_ntill_relative_flux, file = paste(data.path, 'ensemble-relative-responses-weighted-mean-ccl-ntill.rds', sep = '/'))
# GCM
ccl_ntill_relative_gcm_flux = cell_crop_area_gcm_wmean(base.path, 'ccl-ntill', lu.path, raster_cc,
                                                       ccl_ntill_bad_run_f)
saveRDS(ccl_ntill_relative_gcm_flux, file = paste(data.path, 'gcm-relative-responses-weighted-mean-ccl-ntill.rds', sep = '/'))
# GCM PRICE
ccl_ntill_relative_gcm_price    = cell_crop_area_gcm_price(base.path, 'ccl-ntill', lu.path, raster_cc,
                                                         ccl_ntill_bad_run_f, 'OECD-Grain-Outlook-2023-2032.csv')
saveRDS(ccl_ntill_relative_gcm_price, file = paste(data.path, 'gcm-relative-grain-price-weighted-mean-ccl-ntill.rds', sep = '/'))
# ABSOLUTE
ccl_ntill_absolute_flux     = cell_crop_area_wmean_abs(base.path, 'ccl-ntill', lu.path, raster_cc,
                                                   ccl_ntill_bad_run_f)
saveRDS(ccl_ntill_absolute_flux, file = paste(data.path, 'ensemble-absolute-responses-weighted-mean-ccl-ntill.rds', sep = '/'))
#-------------------------------------------------------------------------------