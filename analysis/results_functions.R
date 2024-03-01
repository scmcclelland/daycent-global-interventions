# DayCent Simulation Model Output Functions for Results Processing for Analysis #

# file name:    results_functions.R
# created:      24 July 2023
# last updated: 28 February 2024
#-----------------------------------------------------------------------------------------
library(data.table)
library(sf)
library(terra)
options(scipen = 999, digits = 4)
#-----------------------------------------------------------------------------------------
# DATA #
#-----------------------------------------------------------------------------------------
country_area             = function(.lu_path, .shp_f, .raster) {
  country.sf  = st_read(paste(.lu_path, .shp_f, sep = '/'))
  country.sf_dt = setDT(as.data.frame(country.sf))
  # CREATE raster
  shp_r       = rast(ext(country.sf), nrow = 360, ncol = 720)
  # CREATE shp as raster
  country_r   = terra::rasterize(country.sf, shp_r, fun = 'sum', "OBJECTID")
  # MATCH resolution of simulation data, dimensions the same
  target.r    = rast(nrow = 360, ncol = 720, resolution = 0.5)
  country_r   = resample(country_r, target.r, method = "near")
  country_r   = focal(country_r, w=9, fun = "modal", na.policy = "only", na.rm = TRUE) # needed to capture all gridid
  names(country_r) = "OBJECTID"
  # CREATE data.frame, merge
  country_r.dt    = as.data.frame(country_r, cells=TRUE, xy=TRUE)
  country_r.dt    = setDT(country_r.dt)
  country_n       = data.table(WB_NAME = country.sf_dt$WB_NAME, ID = country.sf_dt$OBJECTID,
                               WB_REGION = country.sf_dt$WB_REGION)
  # BIND to cell numbers
  country_r.dt = country_r.dt[country_n, on = .(OBJECTID = ID)]
  # JOIN with crop area table
  crop_area_r           = rast(paste(.lu_path, .raster, sep = '/'))
  crop_area_r$cell_area = cellSize(crop_area_r, mask=FALSE, lyrs=FALSE, unit="ha")
  crop_area_dt          = as.data.table(terra::as.data.frame(crop_area_r, xy = TRUE, cells = TRUE))
  
  crop_area_dt = crop_area_dt[country_r.dt[, .(cell, WB_NAME, WB_REGION)], on = .(cell = cell)]
  gc()
  setorder(crop_area_dt, cell)
  crop_area_dt = crop_area_dt[!is.na(x),]
  gc()
  
  # NA to 0
  crop_area_dt[is.na(maize_rainfed_2015), maize_rainfed_2015 := 0]
  crop_area_dt[is.na(maize_irrigated_2015), maize_irrigated_2015 := 0]
  crop_area_dt[is.na(soybean_rainfed_2015), soybean_rainfed_2015 := 0]
  crop_area_dt[is.na(soybean_irrigated_2015), soybean_irrigated_2015 := 0]
  crop_area_dt[is.na(wheat_rainfed_2015), wheat_rainfed_2015 := 0]
  crop_area_dt[is.na(wheat_irrigated_2015), wheat_irrigated_2015 := 0]
  
  # SUM AREA BY COUNTRY, REGION
  sum_cols     = c('maize_rainfed_2015', 'maize_irrigated_2015', 'soybean_rainfed_2015', 'soybean_irrigated_2015',
                   'wheat_rainfed_2015', 'wheat_irrigated_2015')
  crop_area_dt = crop_area_dt[, lapply(.SD, sum), .SDcols = sum_cols, by = .(WB_NAME, WB_REGION)]
  crop_area_dt[, crop_sum_2015 := maize_rainfed_2015 + maize_irrigated_2015 + soybean_rainfed_2015 +
                 soybean_irrigated_2015 + wheat_rainfed_2015 + wheat_irrigated_2015]
  return(crop_area_dt)
}
ipcc_name                = function(.lu_path, .shp_f, .raster) {
  country.sf  = st_read(paste(.lu_path, .shp_f, sep = '/'))
  country.sf_dt = setDT(as.data.frame(country.sf))
  # CREATE raster
  shp_r       = rast(ext(country.sf), nrow = 360, ncol = 720)
  # CREATE shp as raster
  country_r   = terra::rasterize(country.sf, shp_r, fun = 'sum', "OBJECTID")
  # MATCH resolution of simulation data, dimensions the same
  target.r    = rast(nrow = 360, ncol = 720, resolution = 0.5)
  country_r   = resample(country_r, target.r, method = "near")
  country_r   = focal(country_r, w=9, fun = "modal", na.policy = "only", na.rm = TRUE) # needed to capture all gridid
  names(country_r) = "OBJECTID"
  # CREATE data.frame, merge
  country_r.dt    = as.data.frame(country_r, cells=TRUE, xy=TRUE)
  country_r.dt    = setDT(country_r.dt)
  country_n       = data.table(WB_NAME = country.sf_dt$WB_NAME, ID = country.sf_dt$OBJECTID,
                               WB_REGION = country.sf_dt$WB_REGION)
  # BIND to cell numbers
  country_r.dt = country_r.dt[country_n, on = .(OBJECTID = ID)]
  # JOIN with crop area table
  crop_area_r           = rast(paste(.lu_path, .raster, sep = '/'))
  crop_area_r$cell_area = cellSize(crop_area_r, mask=FALSE, lyrs=FALSE, unit="ha")
  crop_area_dt          = as.data.table(terra::as.data.frame(crop_area_r, xy = TRUE, cells = TRUE))
  
  crop_area_dt = crop_area_dt[country_r.dt[, .(cell, WB_NAME, WB_REGION)], on = .(cell = cell)]
  gc()
  setorder(crop_area_dt, cell)
  crop_area_dt = crop_area_dt[!is.na(x),]
  gc()
  
  # NA to 0
  crop_area_dt[is.na(maize_rainfed_2015), maize_rainfed_2015 := 0]
  crop_area_dt[is.na(maize_irrigated_2015), maize_irrigated_2015 := 0]
  crop_area_dt[is.na(soybean_rainfed_2015), soybean_rainfed_2015 := 0]
  crop_area_dt[is.na(soybean_irrigated_2015), soybean_irrigated_2015 := 0]
  crop_area_dt[is.na(wheat_rainfed_2015), wheat_rainfed_2015 := 0]
  crop_area_dt[is.na(wheat_irrigated_2015), wheat_irrigated_2015 := 0]
  
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
  crop_area_dt[WB_NAME %in% AME, IPCC_NAME := 'AME']
  crop_area_dt[WB_NAME %in% ADP, IPCC_NAME := 'ADP']
  crop_area_dt[WB_NAME %in% DEV, IPCC_NAME := 'DEV']
  crop_area_dt[WB_NAME %in% EEWCA, IPCC_NAME := 'EEWCA']
  crop_area_dt[WB_NAME %in% LAC, IPCC_NAME := 'LAC']
  return(crop_area_dt)
}
bad_run_year             = function(.scenario, .base.path, .store.path1, .store.path2) {
  
  if(!.scenario %like% '-res') {
    # gcm
    print('This is not a -res scenario. Loading run 1 dt.')
    load(paste0(.store.path1, '/', .scenario,'-processed-output.Rdata'))
    rds_dt1     = rds_dt
    rm(rds_dt)
    gc()
    rds_dt1[crop %in% 'wwht', crop := 'wht']
    rds_dt1[crop %in% 'swht', crop := 'wht']
    # second run
    print('This is not a -res scenario. Loading run 2 dt.')
    load(paste0(.store.path2, '/', .scenario,'-processed-output-2.Rdata'))
    rds_dt[crop %in% 'wwht', crop := 'wht']
    rds_dt[crop %in% 'swht', crop := 'wht']
    gc()
    
    # COMBINE
    rds_dt     = rbind(rds_dt1, rds_dt)
    rm(rds_dt1)
    gc()
    setorder(rds_dt, gridid)
    gc()
  } else {
    print('This is a -res scenario. Loading dt.')
    load(paste0(.store.path1, '/', .scenario,'-processed-output.Rdata'))
    rds_dt[crop %in% 'wwht', crop := 'wht']
    rds_dt[crop %in% 'swht', crop := 'wht']
    setorder(rds_dt, gridid)
    gc()
  }
  
  # DISAGGREGATED INCOMPLETE
  bad_run_yrs_m0  = unique(rds_dt[run_yrs < 85 & crop %in% 'maiz' & irr == 0, gridid])
  bad_run_yrs_m1  = unique(rds_dt[run_yrs < 85 & crop %in% 'maiz' & irr == 1, gridid])
  bad_run_yrs_s0  = unique(rds_dt[run_yrs < 85 & crop %in% 'soyb' & irr == 0, gridid])
  bad_run_yrs_s1  = unique(rds_dt[run_yrs < 85 & crop %in% 'soyb' & irr == 1, gridid])
  bad_run_yrs_w0  = unique(rds_dt[run_yrs < 85 & crop %in% 'wht' & irr == 0, gridid])
  bad_run_yrs_w1  = unique(rds_dt[run_yrs < 85 & crop %in% 'wht' & irr == 1, gridid])
  
  # CREATE DT
  bad_run_yrs_m0  = data.table(gridid = bad_run_yrs_m0, crop = 'maiz', irr = 0)
  bad_run_yrs_m1  = data.table(gridid = bad_run_yrs_m1, crop = 'maiz', irr = 1)
  bad_run_yrs_s0  = data.table(gridid = bad_run_yrs_s0, crop = 'soyb', irr = 0)
  bad_run_yrs_s1  = data.table(gridid = bad_run_yrs_s1, crop = 'soyb', irr = 1)
  bad_run_yrs_w0  = data.table(gridid = bad_run_yrs_w0, crop = 'wht', irr = 0)
  bad_run_yrs_w1  = data.table(gridid = bad_run_yrs_w1, crop = 'wht', irr = 1)
  
  # RBIND
  bad_run_yrs     = rbind(bad_run_yrs_m0, bad_run_yrs_m1, bad_run_yrs_s0,
                          bad_run_yrs_s1, bad_run_yrs_w0, bad_run_yrs_w1)
  bad_run_yrs[, scenario := .scenario]
  
  # SAVE
  fwrite(bad_run_yrs, paste0(.base.path, '/data/bad-run-years-for-', .scenario,'-gridid.csv'))
  rm(rds_dt)
  gc()
}
bad_run_match            = function(.scen_bad_run_dt, .conv_bad_run_dt) {
  # MAIZE 0
    # SCENARIO
  scen_m0 = .scen_bad_run_dt[crop %in% 'maiz' & irr == 0, gridid]
    # CONV
  conv_m0 = .conv_bad_run_dt[crop %in% 'maiz' & irr == 0, gridid]
  if (identical(scen_m0, conv_m0)) {
    print('Gridid for m0 are the same. Defaulting to scenario gridid.')
  } else {
    print('Gridid for m0 are different. Checking for non-matches.')
    if(length(conv_m0[!conv_m0 %in% scen_m0]) > 0L) {
      print('Conv contains additional gridid. Adding conv gridid to scenario gridid.')
      add     = conv_m0[!conv_m0 %in% scen_m0]
      scen_m0 = c(scen_m0, add)
    } else if (length(scen_m0[!scen_m0 %in% conv_m0]) > 0L) {
      print('Scenario contains additional gridid. Adding scenario gridid to conv gridid.')
      add     = scen_m0[!scen_m0 %in% conv_m0]
      scen_m0 = c(conv_m0, add)
    } else {
      print('This case should not exist! Recheck function.')
    }
  } 
  # MAIZE 1
    # SCENARIO
  scen_m1 = .scen_bad_run_dt[crop %in% 'maiz' & irr == 1, gridid]
    # CONV
  conv_m1 = .conv_bad_run_dt[crop %in% 'maiz' & irr == 1, gridid]
  if (identical(scen_m1, conv_m1)) {
    print('Gridid for m1 are the same. Defaulting to scenario gridid.')
  } else {
    print('Gridid for m1 are different. Checking for non-matches.')
    if(length(conv_m1[!conv_m1 %in% scen_m1]) > 0L) {
      print('Conv contains additional gridid. Adding conv gridid to scenario gridid.')
      add     = conv_m1[!conv_m1 %in% scen_m1]
      scen_m1 = c(scen_m1, add)
    } else if (length(scen_m1[!scen_m1 %in% conv_m1]) > 0L) {
      print('Scenario contains additional gridid. Adding scenario gridid to conv gridid.')
      add     = scen_m1[!scen_m1 %in% conv_m1]
      scen_m1 = c(conv_m1, add)
    } else {
      print('This case should not exist! Recheck function.')
    }
  } 
  # SOYBEAN 0
    # SCENARIO
  scen_s0 = .scen_bad_run_dt[crop %in% 'soyb' & irr == 0, gridid]
    # CONV
  conv_s0 = .conv_bad_run_dt[crop %in% 'soyb' & irr == 0, gridid]
  if (identical(scen_s0, conv_s0)) {
    print('Gridid for s0 are the same. Defaulting to scenario gridid.')
  } else {
    print('Gridid for s0 are different. Checking for non-matches.')
    if(length(conv_s0[!conv_s0 %in% scen_s0]) > 0L) {
      print('Conv contains additional gridid. Adding conv gridid to scenario gridid.')
      add     = conv_s0[!conv_s0 %in% scen_s0]
      scen_s0 = c(scen_s0, add)
    } else if (length(scen_s0[!scen_s0 %in% conv_s0]) > 0L) {
      print('Scenario contains additional gridid. Adding scenario gridid to conv gridid.')
      add     = scen_s0[!scen_s0 %in% conv_s0]
      scen_s0 = c(conv_s0, add)
    } else {
      print('This case should not exist! Recheck function.')
    }
  } 
  # SOYBEAN 1
    # SCENARIO
  scen_s1 = .scen_bad_run_dt[crop %in% 'soyb' & irr == 1, gridid]
    # CONV
  conv_s1 = .conv_bad_run_dt[crop %in% 'soyb' & irr == 1, gridid]
  if (identical(scen_s1, conv_s1)) {
    print('Gridid for s1 are the same. Defaulting to scenario gridid.')
  } else {
    print('Gridid for s1 are different. Checking for non-matches.')
    if(length(conv_s1[!conv_s1 %in% scen_s1]) > 0L) {
      print('Conv contains additional gridid. Adding conv gridid to scenario gridid.')
      add     = conv_s1[!conv_s1 %in% scen_s1]
      scen_s1 = c(scen_s1, add)
    } else if (length(scen_s1[!scen_s1 %in% conv_s1]) > 0L) {
      print('Scenario contains additional gridid. Adding scenario gridid to conv gridid.')
      add     = scen_s1[!scen_s1 %in% conv_s1]
      scen_s1 = c(conv_s1, add)
    } else {
      print('This case should not exist! Recheck function.')
    }
  } 
  # WHT 0
    # SCENARIO
  scen_w0 = .scen_bad_run_dt[crop %in% 'wht' & irr == 0, gridid]
    # CONV
  conv_w0 = .conv_bad_run_dt[crop %in% 'wht' & irr == 0, gridid]
  if (identical(scen_w0, conv_w0)) {
    print('Gridid for w0 are the same. Defaulting to scenario gridid.')
  } else {
    print('Gridid for w0 are different. Checking for non-matches.')
    if(length(conv_w0[!conv_w0 %in% scen_w0]) > 0L) {
      print('Conv contains additional gridid. Adding conv gridid to scenario gridid.')
      add     = conv_w0[!conv_w0 %in% scen_w0]
      scen_w0 = c(scen_w0, add)
    } else if (length(scen_w0[!scen_w0 %in% conv_w0]) > 0L) {
      print('Scenario contains additional gridid. Adding scenario gridid to conv gridid.')
      add     = scen_w0[!scen_w0 %in% conv_w0]
      scen_w0 = c(conv_w0, add)
    } else {
      print('This case should not exist! Recheck function.')
    }
  } 
  # WHT 1
    # SCENARIO
  scen_w1 = .scen_bad_run_dt[crop %in% 'wht' & irr == 1, gridid]
    # CONV
  conv_w1 = .conv_bad_run_dt[crop %in% 'wht' & irr == 1, gridid]
  if (identical(scen_w1, conv_w1)) {
    print('Gridid for w1 are the same. Defaulting to scenario gridid.')
  } else {
    print('Gridid for w1 are different. Checking for non-matches.')
    if(length(conv_w1[!conv_w1 %in% scen_w1]) > 0L) {
      print('Conv contains additional gridid. Adding conv gridid to scenario gridid.')
      add     = conv_w1[!conv_w1 %in% scen_w1]
      scen_w1 = c(scen_w1, add)
    } else if (length(scen_w1[!scen_w1 %in% conv_w1]) > 0L) {
      print('Scenario contains additional gridid. Adding scenario gridid to conv gridid.')
      add     = scen_w1[!scen_w1 %in% conv_w1]
      scen_w1 = c(conv_w1, add)
    } else {
      print('This case should not exist! Recheck function.')
    }
  } 
  crop_irr_gr = list(m0 = scen_m0, m1 = scen_m1, s0 = scen_s0, s1 = scen_s1,
                     w0 = scen_w0, w1 = scen_w1)
  return(crop_irr_gr)
  
}
cell_crop_area_wmean     = function(.base.path, .scenario, .lu_path, .raster,
                                    .bad_run_match) {
  if(!.scenario %like% '-res') {
    # gcm
    #first run
    print('This is not a -res scenario. Loading run 1 dt.')
    load(paste0(.base.path, '/data/', .scenario,'-ensemble-relative-responses.Rdata'))
    rds_dt_r_ens1 = rds_dt_r_s_ens
    rm(rds_dt_r_s_ens)
    gc()
    # second run
    print('This is not a -res scenario. Loading run 2 dt.')
    load(paste0(.base.path, '/data/', .scenario,'-ensemble-relative-responses-2.Rdata'))
    gc()
    
    # COMBINE
    rds_dt_r_s_ens = rbind(rds_dt_r_ens1, rds_dt_r_s_ens)
    rm(rds_dt_r_ens1)
    gc()
    setorder(rds_dt_r_s_ens, gridid)
  } else {
    print('This is a -res scenario. Loading dt.')
    load(paste0(.base.path, '/data/', .scenario,'-ensemble-relative-responses.Rdata'))
    setorder(rds_dt_r_s_ens, gridid)
    gc()
  }
  
  # CREATE dt
  crop_area_r           = rast(paste(.lu_path, .raster, sep = '/'))
  crop_area_r$cell_area = cellSize(crop_area_r, mask=FALSE, lyrs=FALSE, unit="ha")
  crop_area_dt          = as.data.table(terra::as.data.frame(crop_area_r, xy = TRUE, cells = TRUE))
  
  # FILTER raster by relevant gridid
  gridid_all       = unique(rds_dt_r_s_ens[,gridid])
  print(paste0('Number of gridid in input table is ', 
               length(gridid_all),
               '.'))
  crop_area_dt_f   = crop_area_dt[cell %in% gridid_all,]
  
  # ADD relevant crop area by crop / irrigation / gridid
  wht_rn = rds_dt_r_s_ens[crop %in% 'wht' & irr == 0,]
  wht_rn = wht_rn[crop_area_dt_f[,c('cell','wheat_rainfed_2015', 'cell_area')], on = .(gridid = cell)]
  wht_rn = wht_rn[!is.na(x),] 
  names(wht_rn)[names(wht_rn) == 'wheat_rainfed_2015'] = 'crop_area_ha'
  
  maiz_rn = rds_dt_r_s_ens[crop %in% 'maiz' & irr == 0,]
  maiz_rn = maiz_rn[crop_area_dt_f[,c('cell','maize_rainfed_2015', 'cell_area')], on = .(gridid = cell)]
  maiz_rn = maiz_rn[!is.na(x),] 
  names(maiz_rn)[names(maiz_rn) == 'maize_rainfed_2015'] = 'crop_area_ha'
  
  soyb_rn = rds_dt_r_s_ens[crop %in% 'soyb' & irr == 0,]
  soyb_rn = soyb_rn[crop_area_dt_f[,c('cell','soybean_rainfed_2015', 'cell_area')], on = .(gridid = cell)]
  soyb_rn = soyb_rn[!is.na(x),]
  names(soyb_rn)[names(soyb_rn) == 'soybean_rainfed_2015'] = 'crop_area_ha'
  
  wht_ir = rds_dt_r_s_ens[crop %in% 'wht' & irr == 1,]
  wht_ir = wht_ir[crop_area_dt_f[,c('cell','wheat_irrigated_2015', 'cell_area')], on = .(gridid = cell)]
  wht_ir = wht_ir[!is.na(x),]
  names(wht_ir)[names(wht_ir) == 'wheat_irrigated_2015'] = 'crop_area_ha'
  
  maiz_ir = rds_dt_r_s_ens[crop %in% 'maiz' & irr == 1,]
  maiz_ir = maiz_ir[crop_area_dt_f[,c('cell','maize_irrigated_2015', 'cell_area')], on = .(gridid = cell)]
  maiz_ir = maiz_ir[!is.na(x),]
  names(maiz_ir)[names(maiz_ir) == 'maize_irrigated_2015'] = 'crop_area_ha'
  
  soyb_ir = rds_dt_r_s_ens[crop %in% 'soyb' & irr == 1,]
  soyb_ir = soyb_ir[crop_area_dt_f[,c('cell','soybean_irrigated_2015', 'cell_area')], on = .(gridid = cell)]
  soyb_ir = soyb_ir[!is.na(x),]
  names(soyb_ir)[names(soyb_ir) == 'soybean_irrigated_2015'] = 'crop_area_ha'
  
  .dt_crop_area = rbind(maiz_rn, soyb_rn, wht_rn, maiz_ir, soyb_ir, wht_ir)
  rm(rds_dt_r_s_ens)
  gc()
  
  # ROUND area
  .dt_crop_area[, crop_area_ha       := round(crop_area_ha, digits = 2)]
  .dt_crop_area[, cell_total_area_ha := round(cell_area, digits = 2)]
  setorder(.dt_crop_area, gridid)
  
  # FILTER gridid >= 1
  print(paste0('Number of gridid by crop and irr with equal or more than 1 ha cropland is ', 
               length(unique(.dt_crop_area[crop_area_ha >= 1, gridid])),
               '.'))
  .dt_crop_area = .dt_crop_area[crop_area_ha >= 1,]
  
  # FILTER BAD RUN YEARS
  .dt_crop_area = .dt_crop_area[!crop %in% 'maiz' | !irr == 0 |
                                  !gridid %in% .bad_run_match$m0]
  .dt_crop_area = .dt_crop_area[!crop %in% 'maiz' | !irr == 1 |
                                  !gridid %in% .bad_run_match$m1]
  .dt_crop_area = .dt_crop_area[!crop %in% 'soyb' | !irr == 0 |
                                  !gridid %in% .bad_run_match$s0]
  .dt_crop_area = .dt_crop_area[!crop %in% 'soyb' | !irr == 1 |
                                  !gridid %in% .bad_run_match$s1]
  .dt_crop_area = .dt_crop_area[!crop %in% 'wht'  | !irr == 0 |
                                  !gridid %in% .bad_run_match$w0]
  .dt_crop_area = .dt_crop_area[!crop %in% 'wht'  | !irr == 1 |
                                  !gridid %in% .bad_run_match$w1]
  print(paste0('Number of 1 ha filtered gridid by crop and irr and good run years is ', 
               length(unique(.dt_crop_area[, gridid])),
               '.'))
  # REMOVE KNOWN OUTLIERS
  .dt_crop_area = .dt_crop_area[!gridid %in% c(103186,164414),]
  print(paste0('Number of 1 ha filtered gridid by crop and irr, good run years, and known outliers is ', 
               length(unique(.dt_crop_area[, gridid])),
               '.'))
  
  # WEIGHTED mean by crop area in gridid
  .dt_crop_area[, total_crop_area_ha := sum(crop_area_ha), by = .(y_block, ssp, gcm, gridid, scenario)]
  .dt_crop_area[, weights            := (crop_area_ha/total_crop_area_ha)*100L]
  cols_for_mean = c('m_cr_rootC', 'm_cr_shootC', 'm_cr_shootN', 'm_cr_grain', 'm_cr_grainN','m_cr_NPP', 'm_cr_residC', 
                    'm_cr_residN', 's_cr_rootC', 's_cr_shootC', 's_cr_shootN', 's_cr_grain', 's_cr_grainN','s_cr_NPP',
                    's_cr_residC',  's_cr_residN', 'm_cc_rootC', 'm_cc_shootC', 'm_cc_shootN','s_cc_rootC', 's_cc_shootC', 
                    's_cc_shootN', 'm_SOC', 'm_N2O', 'm_iN2O', 'm_dN2O','m_GHG', 'm_annet','m_sfdcmp', 'm_sldcmp',  'm_cr_irr',
                    'm_sC.N','m_nfix','s_cr_irr', 's_annet','s_SOC', 's_N2O', 's_iN2O', 's_dN2O', 
                    's_GHG','s_gr_nit')
  print('Compute weighted mean by crop area in gridcell. Note: Calculation takes several minutes.')
  .dt_crop_area = .dt_crop_area[, lapply(.SD, weighted.mean, weights), .SDcols = cols_for_mean, 
                                by = .(gridid, x, y, scenario, y_block, gcm, ssp, WB_NAME, WB_REGION,total_crop_area_ha, cell_area)]
  .dt_crop_area = .dt_crop_area[, lapply(.SD, round, digits = 2), .SDcols = cols_for_mean,
                                by = .(gridid, x, y, scenario, y_block, gcm, ssp, WB_NAME, WB_REGION, total_crop_area_ha, cell_area)]
  setorder(.dt_crop_area, gridid)
  
  .dt_crop_area[, total_crop_area_ha := NULL]
  .dt_crop_area[, cell_area := NULL]
  
  print(paste0('Number of crop area weighted gridid is ', 
               length(unique(.dt_crop_area[, gridid])),
               '.'))
  return(.dt_crop_area)
}
cell_crop_area_gcm_wmean     = function(.base.path, .scenario, .lu_path, .raster,
                                    .bad_run_match) {
  if(!.scenario %like% '-res') {
    # gcm
    #first run
    print('This is not a -res scenario. Loading run 1 dt.')
    load(paste0(.base.path, '/data/', 'relative-estimates-scenario-counterfactual-', .scenario,
                '-scenario.Rdata'))
    rds_dt_r_s1 = rds_dt_r_s
    rm(rds_dt_r_s)
    gc()
    # second run
    print('This is not a -res scenario. Loading run 2 dt.')
    load(paste0(.base.path, '/data/', 'relative-estimates-scenario-counterfactual-', .scenario,
                '-scenario-2.Rdata'))
    gc()
    
    # COMBINE
    rds_dt_r_s = rbind(rds_dt_r_s1, rds_dt_r_s)
    rm(rds_dt_r_s1)
    gc()
    setorder(rds_dt_r_s, gridid)
    
    # REMOVE INCOMPLETE GCM + SSP RUNS
    rds_dt_r_s = rds_dt_r_s[!gcm %in% 'ACCESS-ESM1-5' | !ssp %in% 'ssp126',]
    gc()
    rds_dt_r_s = rds_dt_r_s[!gcm %in% 'CMCC-ESM2' | !ssp %in% 'ssp370',]
    gc()
    
    # REDUCE COLUMNS
    cols = c('gridid', 'x', 'y', 'scenario', 'y_block', 'gcm', 'ssp', 'WB_NAME', 
             'WB_REGION', 'crop', 'irr','m_cr_grain','m_SOC','m_N2O','m_iN2O',
             'm_dN2O','m_GHG','s_cr_grain','s_SOC','s_N2O','s_iN2O','s_dN2O','s_GHG')
    rds_dt_r_s = rds_dt_r_s[, ..cols]
    gc()
  } else {
    print('This is a -res scenario. Loading dt.')
    load(paste0(.base.path, '/data/', 'relative-estimates-scenario-counterfactual-', .scenario,
                '-scenario.Rdata'))
    setorder(rds_dt_r_s, gridid)
    gc()
    
    # REMOVE INCOMPLETE GCM + SSP RUNS
    rds_dt_r_s = rds_dt_r_s[!gcm %in% 'ACCESS-ESM1-5' | !ssp %in% 'ssp126',]
    gc()
    rds_dt_r_s = rds_dt_r_s[!gcm %in% 'CMCC-ESM2' | !ssp %in% 'ssp370',]
    gc()
    
    # REDUCE COLUMNS
    cols = c('gridid', 'x', 'y', 'scenario', 'y_block', 'gcm', 'ssp', 'WB_NAME', 
             'WB_REGION', 'crop', 'irr','m_cr_grain','m_SOC','m_N2O','m_iN2O',
             'm_dN2O','m_GHG','s_cr_grain','s_SOC','s_N2O','s_iN2O','s_dN2O','s_GHG')
    rds_dt_r_s = rds_dt_r_s[, ..cols]
    gc()
  }
  
  # CREATE dt
  crop_area_r           = rast(paste(.lu_path, .raster, sep = '/'))
  crop_area_r$cell_area = cellSize(crop_area_r, mask=FALSE, lyrs=FALSE, unit="ha")
  crop_area_dt          = as.data.table(terra::as.data.frame(crop_area_r, xy = TRUE, cells = TRUE))
  
  # FILTER raster by relevant gridid
  gridid_all       = unique(rds_dt_r_s[,gridid])
  print(paste0('Number of gridid in input table is ', 
               length(gridid_all),
               '.'))
  crop_area_dt_f   = crop_area_dt[cell %in% gridid_all,]
  
  # ADD relevant crop area by crop / irrigation / gridid
  wht_rn = rds_dt_r_s[crop %in% 'wht' & irr == 0,]
  wht_rn = wht_rn[crop_area_dt_f[,c('cell','wheat_rainfed_2015', 'cell_area')], on = .(gridid = cell)]
  wht_rn = wht_rn[!is.na(x),] 
  names(wht_rn)[names(wht_rn) == 'wheat_rainfed_2015'] = 'crop_area_ha'
  gc()
  
  maiz_rn = rds_dt_r_s[crop %in% 'maiz' & irr == 0,]
  maiz_rn = maiz_rn[crop_area_dt_f[,c('cell','maize_rainfed_2015', 'cell_area')], on = .(gridid = cell)]
  maiz_rn = maiz_rn[!is.na(x),] 
  names(maiz_rn)[names(maiz_rn) == 'maize_rainfed_2015'] = 'crop_area_ha'
  gc()
  
  soyb_rn = rds_dt_r_s[crop %in% 'soyb' & irr == 0,]
  soyb_rn = soyb_rn[crop_area_dt_f[,c('cell','soybean_rainfed_2015', 'cell_area')], on = .(gridid = cell)]
  soyb_rn = soyb_rn[!is.na(x),]
  names(soyb_rn)[names(soyb_rn) == 'soybean_rainfed_2015'] = 'crop_area_ha'
  gc()
  
  wht_ir = rds_dt_r_s[crop %in% 'wht' & irr == 1,]
  wht_ir = wht_ir[crop_area_dt_f[,c('cell','wheat_irrigated_2015', 'cell_area')], on = .(gridid = cell)]
  wht_ir = wht_ir[!is.na(x),]
  names(wht_ir)[names(wht_ir) == 'wheat_irrigated_2015'] = 'crop_area_ha'
  gc()
  
  maiz_ir = rds_dt_r_s[crop %in% 'maiz' & irr == 1,]
  maiz_ir = maiz_ir[crop_area_dt_f[,c('cell','maize_irrigated_2015', 'cell_area')], on = .(gridid = cell)]
  maiz_ir = maiz_ir[!is.na(x),]
  names(maiz_ir)[names(maiz_ir) == 'maize_irrigated_2015'] = 'crop_area_ha'
  gc()
  
  soyb_ir = rds_dt_r_s[crop %in% 'soyb' & irr == 1,]
  soyb_ir = soyb_ir[crop_area_dt_f[,c('cell','soybean_irrigated_2015', 'cell_area')], on = .(gridid = cell)]
  soyb_ir = soyb_ir[!is.na(x),]
  names(soyb_ir)[names(soyb_ir) == 'soybean_irrigated_2015'] = 'crop_area_ha'
  gc()
  
  .dt_crop_area = rbind(maiz_rn, soyb_rn, wht_rn, maiz_ir, soyb_ir, wht_ir)
  rm(rds_dt_r_s)
  gc()
  
  # ROUND area
  .dt_crop_area[, crop_area_ha       := round(crop_area_ha, digits = 2)]
  .dt_crop_area[, cell_total_area_ha := round(cell_area, digits = 2)]
  setorder(.dt_crop_area, gridid)
  
  # FILTER gridid >= 1
  print(paste0('Number of gridid by crop and irr with equal or more than 1 ha cropland is ', 
               length(unique(.dt_crop_area[crop_area_ha >= 1, gridid])),
               '.'))
  .dt_crop_area = .dt_crop_area[crop_area_ha >= 1,]
  
  # FILTER BAD RUN YEARS
  .dt_crop_area = .dt_crop_area[!crop %in% 'maiz' | !irr == 0 |
                                  !gridid %in% .bad_run_match$m0]
  .dt_crop_area = .dt_crop_area[!crop %in% 'maiz' | !irr == 1 |
                                  !gridid %in% .bad_run_match$m1]
  .dt_crop_area = .dt_crop_area[!crop %in% 'soyb' | !irr == 0 |
                                  !gridid %in% .bad_run_match$s0]
  .dt_crop_area = .dt_crop_area[!crop %in% 'soyb' | !irr == 1 |
                                  !gridid %in% .bad_run_match$s1]
  .dt_crop_area = .dt_crop_area[!crop %in% 'wht'  | !irr == 0 |
                                  !gridid %in% .bad_run_match$w0]
  .dt_crop_area = .dt_crop_area[!crop %in% 'wht'  | !irr == 1 |
                                  !gridid %in% .bad_run_match$w1]
  print(paste0('Number of 1 ha filtered gridid by crop and irr and good run years is ', 
               length(unique(.dt_crop_area[, gridid])),
               '.'))
  # REMOVE KNOWN OUTLIERS
  .dt_crop_area = .dt_crop_area[!gridid %in% c(103186,164414),]
  print(paste0('Number of 1 ha filtered gridid by crop and irr, good run years, and known outliers is ', 
               length(unique(.dt_crop_area[, gridid])),
               '.'))
  
  # WEIGHTED mean by crop area in gridid
  .dt_crop_area[, total_crop_area_ha := sum(crop_area_ha), by = .(y_block, ssp, gcm, gridid, scenario)]
  .dt_crop_area[, weights            := (crop_area_ha/total_crop_area_ha)*100L]
  cols_for_mean = c('m_cr_grain','m_SOC','m_N2O','m_iN2O','m_dN2O',
                    'm_GHG','s_cr_grain','s_SOC','s_N2O','s_iN2O','s_dN2O','s_GHG')
  print('Compute weighted mean by crop area in gridcell. Note: Calculation takes several minutes.')
  .dt_crop_area = .dt_crop_area[, lapply(.SD, weighted.mean, weights), .SDcols = cols_for_mean, 
                                by = .(gridid, x, y, scenario, y_block, gcm, ssp, WB_NAME, WB_REGION, total_crop_area_ha, cell_area)]
  gc()
  .dt_crop_area = .dt_crop_area[, lapply(.SD, round, digits = 2), .SDcols = cols_for_mean,
                                by = .(gridid, x, y, scenario, y_block, gcm, ssp, WB_NAME, WB_REGION, total_crop_area_ha, cell_area)]
  gc()
  setorder(.dt_crop_area, gridid)
  
  .dt_crop_area[, total_crop_area_ha := NULL]
  .dt_crop_area[, cell_area := NULL]
  
  print(paste0('Number of crop area weighted gridid is ', 
               length(unique(.dt_crop_area[, gridid])),
               '.'))
  return(.dt_crop_area)
}
cell_crop_area_gcm_var_wmean = function(.base.path, .scenario, .lu_path, .raster,
                                    .bad_run_match) {
  if(!.scenario %like% '-res') {
    # gcm
    #first run
    print('This is not a -res scenario. Loading run 1 dt.')
    load(paste0(.base.path, '/data/', 'relative-estimates-scenario-counterfactual-', .scenario,
                '-scenario.Rdata'))
    rds_dt_r_s1 = rds_dt_r_s
    rm(rds_dt_r_s)
    gc()
    # second run
    print('This is not a -res scenario. Loading run 2 dt.')
    load(paste0(.base.path, '/data/', 'relative-estimates-scenario-counterfactual-', .scenario,
                '-scenario-2.Rdata'))
    gc()
    
    # COMBINE
    rds_dt_r_s = rbind(rds_dt_r_s1, rds_dt_r_s)
    rm(rds_dt_r_s1)
    gc()
    setorder(rds_dt_r_s, gridid)
    
    # REMOVE INCOMPLETE GCM + SSP RUNS
    rds_dt_r_s = rds_dt_r_s[!gcm %in% 'ACCESS-ESM1-5' | !ssp %in% 'ssp126',]
    gc()
    rds_dt_r_s = rds_dt_r_s[!gcm %in% 'CMCC-ESM2' | !ssp %in% 'ssp370',]
    gc()
    
    # REDUCE COLUMNS
    cols = c('gridid', 'x', 'y', 'scenario', 'y_block', 'gcm', 'ssp', 'WB_NAME', 
             'WB_REGION', 'crop', 'irr','s_cr_grain', 's_cr_NPP', 's_cr_residC','s_cc_rootC', 's_cc_shootC', 
             'm_sfdcmp', 'm_sldcmp', 'm_sC.N','m_nfix','s_annet', 's_gr_nit')
    rds_dt_r_s = rds_dt_r_s[, ..cols]
    gc()
  } else {
    print('This is a -res scenario. Loading dt.')
    load(paste0(.base.path, '/data/', 'relative-estimates-scenario-counterfactual-', .scenario,
                '-scenario.Rdata'))
    setorder(rds_dt_r_s, gridid)
    gc()
    
    # REMOVE INCOMPLETE GCM + SSP RUNS
    rds_dt_r_s = rds_dt_r_s[!gcm %in% 'ACCESS-ESM1-5' | !ssp %in% 'ssp126',]
    gc()
    rds_dt_r_s = rds_dt_r_s[!gcm %in% 'CMCC-ESM2' | !ssp %in% 'ssp370',]
    gc()
    
    # REDUCE COLUMNS
    cols = c('gridid', 'x', 'y', 'scenario', 'y_block', 'gcm', 'ssp', 'WB_NAME', 
             'WB_REGION', 'crop', 'irr','s_cr_grain', 's_cr_NPP', 's_cr_residC','s_cc_rootC', 's_cc_shootC', 
             'm_sfdcmp', 'm_sldcmp', 'm_sC.N','m_nfix','s_annet', 's_gr_nit')
    rds_dt_r_s = rds_dt_r_s[, ..cols]
    gc()
  }
  
  # CREATE dt
  crop_area_r           = rast(paste(.lu_path, .raster, sep = '/'))
  crop_area_r$cell_area = cellSize(crop_area_r, mask=FALSE, lyrs=FALSE, unit="ha")
  crop_area_dt          = as.data.table(terra::as.data.frame(crop_area_r, xy = TRUE, cells = TRUE))
  
  # FILTER raster by relevant gridid
  gridid_all       = unique(rds_dt_r_s[,gridid])
  print(paste0('Number of gridid in input table is ', 
               length(gridid_all),
               '.'))
  crop_area_dt_f   = crop_area_dt[cell %in% gridid_all,]
  
  # ADD relevant crop area by crop / irrigation / gridid
  wht_rn = rds_dt_r_s[crop %in% 'wht' & irr == 0,]
  wht_rn = wht_rn[crop_area_dt_f[,c('cell','wheat_rainfed_2015', 'cell_area')], on = .(gridid = cell)]
  wht_rn = wht_rn[!is.na(x),] 
  names(wht_rn)[names(wht_rn) == 'wheat_rainfed_2015'] = 'crop_area_ha'
  gc()
  
  maiz_rn = rds_dt_r_s[crop %in% 'maiz' & irr == 0,]
  maiz_rn = maiz_rn[crop_area_dt_f[,c('cell','maize_rainfed_2015', 'cell_area')], on = .(gridid = cell)]
  maiz_rn = maiz_rn[!is.na(x),] 
  names(maiz_rn)[names(maiz_rn) == 'maize_rainfed_2015'] = 'crop_area_ha'
  gc()
  
  soyb_rn = rds_dt_r_s[crop %in% 'soyb' & irr == 0,]
  soyb_rn = soyb_rn[crop_area_dt_f[,c('cell','soybean_rainfed_2015', 'cell_area')], on = .(gridid = cell)]
  soyb_rn = soyb_rn[!is.na(x),]
  names(soyb_rn)[names(soyb_rn) == 'soybean_rainfed_2015'] = 'crop_area_ha'
  gc()
  
  wht_ir = rds_dt_r_s[crop %in% 'wht' & irr == 1,]
  wht_ir = wht_ir[crop_area_dt_f[,c('cell','wheat_irrigated_2015', 'cell_area')], on = .(gridid = cell)]
  wht_ir = wht_ir[!is.na(x),]
  names(wht_ir)[names(wht_ir) == 'wheat_irrigated_2015'] = 'crop_area_ha'
  gc()
  
  maiz_ir = rds_dt_r_s[crop %in% 'maiz' & irr == 1,]
  maiz_ir = maiz_ir[crop_area_dt_f[,c('cell','maize_irrigated_2015', 'cell_area')], on = .(gridid = cell)]
  maiz_ir = maiz_ir[!is.na(x),]
  names(maiz_ir)[names(maiz_ir) == 'maize_irrigated_2015'] = 'crop_area_ha'
  gc()
  
  soyb_ir = rds_dt_r_s[crop %in% 'soyb' & irr == 1,]
  soyb_ir = soyb_ir[crop_area_dt_f[,c('cell','soybean_irrigated_2015', 'cell_area')], on = .(gridid = cell)]
  soyb_ir = soyb_ir[!is.na(x),]
  names(soyb_ir)[names(soyb_ir) == 'soybean_irrigated_2015'] = 'crop_area_ha'
  gc()
  
  .dt_crop_area = rbind(maiz_rn, soyb_rn, wht_rn, maiz_ir, soyb_ir, wht_ir)
  rm(rds_dt_r_s)
  gc()
  
  # ROUND area
  .dt_crop_area[, crop_area_ha       := round(crop_area_ha, digits = 2)]
  .dt_crop_area[, cell_total_area_ha := round(cell_area, digits = 2)]
  setorder(.dt_crop_area, gridid)
  
  # FILTER gridid >= 1
  print(paste0('Number of gridid by crop and irr with equal or more than 1 ha cropland is ', 
               length(unique(.dt_crop_area[crop_area_ha >= 1, gridid])),
               '.'))
  .dt_crop_area = .dt_crop_area[crop_area_ha >= 1,]
  
  # FILTER BAD RUN YEARS
  .dt_crop_area = .dt_crop_area[!crop %in% 'maiz' | !irr == 0 |
                                  !gridid %in% .bad_run_match$m0]
  .dt_crop_area = .dt_crop_area[!crop %in% 'maiz' | !irr == 1 |
                                  !gridid %in% .bad_run_match$m1]
  .dt_crop_area = .dt_crop_area[!crop %in% 'soyb' | !irr == 0 |
                                  !gridid %in% .bad_run_match$s0]
  .dt_crop_area = .dt_crop_area[!crop %in% 'soyb' | !irr == 1 |
                                  !gridid %in% .bad_run_match$s1]
  .dt_crop_area = .dt_crop_area[!crop %in% 'wht'  | !irr == 0 |
                                  !gridid %in% .bad_run_match$w0]
  .dt_crop_area = .dt_crop_area[!crop %in% 'wht'  | !irr == 1 |
                                  !gridid %in% .bad_run_match$w1]
  print(paste0('Number of 1 ha filtered gridid by crop and irr and good run years is ', 
               length(unique(.dt_crop_area[, gridid])),
               '.'))
  # REMOVE KNOWN OUTLIERS
  .dt_crop_area = .dt_crop_area[!gridid %in% c(103186,164414),]
  print(paste0('Number of 1 ha filtered gridid by crop and irr, good run years, and known outliers is ', 
               length(unique(.dt_crop_area[, gridid])),
               '.'))
  
  # WEIGHTED mean by crop area in gridid
  .dt_crop_area[, total_crop_area_ha := sum(crop_area_ha), by = .(y_block, ssp, gcm, gridid, scenario)]
  .dt_crop_area[, weights            := (crop_area_ha/total_crop_area_ha)*100L]
  cols_for_mean = c('s_cr_grain', 's_cr_NPP', 's_cr_residC','s_cc_rootC', 's_cc_shootC', 
                    'm_sfdcmp', 'm_sldcmp', 'm_sC.N','m_nfix','s_annet', 's_gr_nit')
  print('Compute weighted mean by crop area in gridcell. Note: Calculation takes several minutes.')
  .dt_crop_area = .dt_crop_area[, lapply(.SD, weighted.mean, weights), .SDcols = cols_for_mean, 
                                by = .(gridid, x, y, scenario, y_block, gcm, ssp, WB_NAME, WB_REGION, total_crop_area_ha, cell_area)]
  gc()
  .dt_crop_area = .dt_crop_area[, lapply(.SD, round, digits = 2), .SDcols = cols_for_mean,
                                by = .(gridid, x, y, scenario, y_block, gcm, ssp, WB_NAME, WB_REGION, total_crop_area_ha, cell_area)]
  gc()
  setorder(.dt_crop_area, gridid)
  
  .dt_crop_area[, total_crop_area_ha := NULL]
  .dt_crop_area[, cell_area := NULL]
  
  print(paste0('Number of crop area weighted gridid is ', 
               length(unique(.dt_crop_area[, gridid])),
               '.'))
  return(.dt_crop_area)
}
impute_missing_ens       = function(.dt, .crop_area) {

  # empty dt
  imp_gcm_dt = data.table(cell = numeric(), x = numeric(), y = numeric(),
                          y_block = numeric(), ssp = character(),
                          variable = factor(), value = numeric(), 
                          total_crop_area_ha = numeric())
  
  ssps  = c('historical', 'ssp126', 'ssp370')
  times = c('2030', '2050', '2100')

  for (.ssp in ssps) {
    print(paste0('Creating imputed table for ', .ssp, '.'))
    for (.time in times) {
      print(paste0('Extracting data for year ', .time, '.'))
          # FILTER DT
          dt_r = .dt[y_block == .time & ssp %in% .ssp,]
          # ROUND values
          dt_r[, s_GHG        := round(s_GHG, digits = 2)]
          dt_r[, s_SOC        := round(s_SOC, digits = 2)]
          dt_r[, s_N2O        := round(s_N2O, digits = 2)]
          dt_r[, s_dN2O       := round(s_dN2O, digits = 2)]
          dt_r[, s_iN2O       := round(s_iN2O, digits = 2)]
          dt_r[, s_cr_grain   := round(s_cr_grain, digits = 2)]
          dt_r[, m_GHG        := round(m_GHG, digits = 2)]
          dt_r[, m_SOC        := round(m_SOC, digits = 2)]
          dt_r[, m_N2O        := round(m_N2O, digits = 2)]
          dt_r[, m_dN2O       := round(m_dN2O, digits = 2)]
          dt_r[, m_iN2O       := round(m_iN2O, digits = 2)]
          dt_r[, m_cr_grain   := round(m_cr_grain, digits = 2)]
          # TO DATAFRAME
          dt_r = as.data.frame(dt_r, xy = TRUE)
          setnames(dt_r, 'gridid', 'cell')
          setDT(dt_r)
          
          # CREATE raster
          dt_r   = as.data.frame(dt_r, xy = TRUE)
          r      = rast(nrow = 360, ncol = 720, nlyr = 12, xmin = -180, xmax = 180, ymin = -90, ymax = 90)
          crs(r) = "epsg:4326"
          # original projection
          r_ghg = rast(res = 0.5, nlyr = 12, extent = ext(r), crs = crs(r)) 
          r_ghg[[1]][dt_r$cell]  = dt_r$s_GHG
          r_ghg[[2]][dt_r$cell]  = dt_r$s_SOC
          r_ghg[[3]][dt_r$cell]  = dt_r$s_N2O
          r_ghg[[4]][dt_r$cell]  = dt_r$s_dN2O
          r_ghg[[5]][dt_r$cell]  = dt_r$s_iN2O
          r_ghg[[6]][dt_r$cell]  = dt_r$s_cr_grain
          r_ghg[[7]][dt_r$cell]  = dt_r$m_GHG
          r_ghg[[8]][dt_r$cell]  = dt_r$m_SOC
          r_ghg[[9]][dt_r$cell]  = dt_r$m_N2O
          r_ghg[[10]][dt_r$cell] = dt_r$m_dN2O
          r_ghg[[11]][dt_r$cell] = dt_r$m_iN2O
          r_ghg[[12]][dt_r$cell] = dt_r$m_cr_grain
          
          # UPDATE NAMES
          lyrs     = c('s_GHG_', 's_SOC_', 's_N2O_','s_dN2O_', 's_iN2O_','s_cr_grain_',
                       'm_GHG_', 'm_SOC_', 'm_N2O_','m_dN2O_', 'm_iN2O_','m_cr_grain_')
          lyrs     = outer(lyrs, .time, paste0)
          lyrs     = outer(lyrs, paste0('_', .ssp), paste0)
          names(r_ghg) = lyrs
          
          # ITERATIVE MOVING WINDOW
          w = c(3,5,7,9,11,13,15,17,19,21,23,25,27,29,31)
          imp_w_dt  = c()
          imp_cells = c()
          for (i in w) {
            # IMPUTE with focal | need large window to capture ranges of table 
            r_all_imp = focal(r_ghg, w=i, fun = "mean", na.policy = "only", na.rm = TRUE)
            
            # TRANSFORM back to dt
            imp_dt    = setDT(terra::as.data.frame(r_all_imp, xy = TRUE, 
                                                   cells = TRUE))
            imp_dt    = imp_dt[!cell %in% imp_cells] # remove old cells
            print(paste0('The difference in cells for window ', i, ' is ', length(unique(imp_dt[, cell]))))
            imp_w_dt  = rbind(imp_w_dt, imp_dt)      # combine tables
            setorder(imp_w_dt, cell)
            imp_cells = unique(imp_w_dt[, cell])     # update cells
            print(paste0('New length of cells in iterative focal step is ',length(imp_cells), ' for window ', i, '.'))
            gc()
          }
          
          # WIDE TO LONG format
          imp_w_dt = melt(imp_w_dt,
                        id.vars = c("cell", "x", "y"),
                        measure.vars = patterns('_'))
          # ADJUST columns
          imp_w_dt[, y_block  := .time]
          imp_w_dt[, ssp      := .ssp]
          imp_w_dt[variable %like% 's_GHG', variable       := 's_GHG']
          imp_w_dt[variable %like% 'm_GHG', variable       := 'm_GHG']
          imp_w_dt[variable %like% 's_SOC', variable       := 's_SOC']
          imp_w_dt[variable %like% 'm_SOC', variable       := 'm_SOC']
          imp_w_dt[variable %like% 's_N2O', variable       := 's_N2O']
          imp_w_dt[variable %like% 'm_N2O', variable       := 'm_N2O']
          imp_w_dt[variable %like% 's_dN2O', variable      := 's_dN2O']
          imp_w_dt[variable %like% 'm_dN2O', variable      := 'm_dN2O']
          imp_w_dt[variable %like% 's_iN2O', variable      := 's_iN2O']
          imp_w_dt[variable %like% 'm_iN2O', variable      := 'm_iN2O']
          imp_w_dt[variable %like% 's_cr_grain', variable  := 's_cr_grain']
          imp_w_dt[variable %like% 'm_cr_grain', variable  := 'm_cr_grain']
          imp_w_dt[, value := round(value, digits = 2)] # round
          setcolorder(imp_w_dt, c('cell','x', 'y', 'y_block', 'ssp', 'variable', 'value'))
          
          # SUM crop area into one layer
          crop_sum = sum(.crop_area, na.rm = TRUE) 
          names(crop_sum) = 'total_crop_area_ha'
          # TRANSFORM to dt
          crop_dt  = setDT(terra::as.data.frame(crop_sum, xy = TRUE, 
                                                cells = TRUE, na.rm =TRUE))
          # JOIN
          imp_w_dt = imp_w_dt[crop_dt, on = .(cell = cell,
                                          x    = x,
                                          y    = y)]
          imp_w_dt = imp_w_dt[!is.na(y_block)]
          imp_w_dt = imp_w_dt[total_crop_area_ha > 0,]
          
          # COMBINE
          imp_gcm_dt = rbind(imp_gcm_dt, imp_w_dt)
          gc()
        }
      }
  return(imp_gcm_dt)
}
impute_missing_gcm       = function(.dt, .crop_area, gcms, .ssp) {
  print('Note: this function takes several minutes to run.')
  # empty dt
  imp_gcm_dt = data.table(cell = numeric(), x = numeric(), y = numeric(),
                          y_block = numeric(), ssp = character(),
                          gcm = character(), variable = factor(), 
                          value = numeric(), total_crop_area_ha = numeric())
  
  # set variables
  times = c('2030', '2050', '2100')
  
  for (.gcm in gcms) {
    print(paste0('Creating imputed table for ', .gcm, '.'))
    for (.time in times) {
      print(paste0('Extracting data for year ', .time, '.'))
      # FILTER DT
      dt_r = .dt[y_block == .time & ssp %in% .ssp & gcm %in% .gcm,]
      # ROUND values
      dt_r[, s_GHG        := round(s_GHG, digits = 2)]
      dt_r[, s_SOC        := round(s_SOC, digits = 2)]
      dt_r[, s_N2O        := round(s_N2O, digits = 2)]
      dt_r[, s_dN2O       := round(s_dN2O, digits = 2)]
      dt_r[, s_iN2O       := round(s_iN2O, digits = 2)]
      dt_r[, s_cr_grain   := round(s_cr_grain, digits = 2)]
      dt_r[, m_GHG        := round(m_GHG, digits = 2)]
      dt_r[, m_SOC        := round(m_SOC, digits = 2)]
      dt_r[, m_N2O        := round(m_N2O, digits = 2)]
      dt_r[, m_dN2O       := round(m_dN2O, digits = 2)]
      dt_r[, m_iN2O       := round(m_iN2O, digits = 2)]
      dt_r[, m_cr_grain   := round(m_cr_grain, digits = 2)]
      # TO DATAFRAME
      dt_r = as.data.frame(dt_r, xy = TRUE)
      setnames(dt_r, 'gridid', 'cell')
      setDT(dt_r)
      
      # CREATE raster
      dt_r   = as.data.frame(dt_r, xy = TRUE)
      r      = rast(nrow = 360, ncol = 720, nlyr = 12, xmin = -180, xmax = 180, ymin = -90, ymax = 90)
      crs(r) = "epsg:4326"
      # original projection
      r_ghg = rast(res = 0.5, nlyr = 12, extent = ext(r), crs = crs(r)) 
      r_ghg[[1]][dt_r$cell]  = dt_r$s_GHG
      r_ghg[[2]][dt_r$cell]  = dt_r$s_SOC
      r_ghg[[3]][dt_r$cell]  = dt_r$s_N2O
      r_ghg[[4]][dt_r$cell]  = dt_r$s_dN2O
      r_ghg[[5]][dt_r$cell]  = dt_r$s_iN2O
      r_ghg[[6]][dt_r$cell]  = dt_r$s_cr_grain
      r_ghg[[7]][dt_r$cell]  = dt_r$m_GHG
      r_ghg[[8]][dt_r$cell]  = dt_r$m_SOC
      r_ghg[[9]][dt_r$cell]  = dt_r$m_N2O
      r_ghg[[10]][dt_r$cell] = dt_r$m_dN2O
      r_ghg[[11]][dt_r$cell] = dt_r$m_iN2O
      r_ghg[[12]][dt_r$cell] = dt_r$m_cr_grain
      
      # UPDATE NAMES
      lyrs     = c('s_GHG_', 's_SOC_', 's_N2O_','s_dN2O_', 's_iN2O_','s_cr_grain_',
                   'm_GHG_', 'm_SOC_', 'm_N2O_','m_dN2O_', 'm_iN2O_','m_cr_grain_')
      lyrs     = outer(lyrs, .time, paste0)
      lyrs     = outer(lyrs, paste0('_', .ssp), paste0)
      names(r_ghg) = lyrs
      gc()

      # ITERATIVE MOVING WINDOW
      w = c(3,5,7,9,11,13,15,17,19,21,23,25,27,29,31)
      imp_w_dt  = c()
      imp_cells = c()
      for (i in w) {
        # IMPUTE with focal | need large window to capture ranges of table 
        r_all_imp = focal(r_ghg, w=i, fun = "mean", na.policy = "only", na.rm = TRUE)
        
        # TRANSFORM back to dt
        imp_dt    = setDT(terra::as.data.frame(r_all_imp, xy = TRUE, 
                                               cells = TRUE))
        imp_dt    = imp_dt[!cell %in% imp_cells] # remove old cells
        print(paste0('The difference in cells for window ', i, ' is ', length(unique(imp_dt[, cell]))))
        imp_w_dt  = rbind(imp_w_dt, imp_dt)      # combine tables
        setorder(imp_w_dt, cell)
        imp_cells = unique(imp_w_dt[, cell])     # update cells
        print(paste0('New length of cells in iterative focal step is ',length(imp_cells), ' for window ', i, '.'))
        gc()
      }
      
      # WIDE TO LONG format
      imp_w_dt = melt(imp_w_dt,
                    id.vars = c("cell", "x", "y"),
                    measure.vars = patterns('_'))
      # ADJUST columns
      imp_w_dt[, ssp := .ssp]
      imp_w_dt[, gcm := .gcm]
      imp_w_dt[variable %like% .time, y_block   := .time]
      imp_w_dt[variable %like% 's_GHG', variable       := 's_GHG']
      imp_w_dt[variable %like% 'm_GHG', variable       := 'm_GHG']
      imp_w_dt[variable %like% 's_SOC', variable       := 's_SOC']
      imp_w_dt[variable %like% 'm_SOC', variable       := 'm_SOC']
      imp_w_dt[variable %like% 's_N2O', variable       := 's_N2O']
      imp_w_dt[variable %like% 'm_N2O', variable       := 'm_N2O']
      imp_w_dt[variable %like% 's_dN2O', variable      := 's_dN2O']
      imp_w_dt[variable %like% 'm_dN2O', variable      := 'm_dN2O']
      imp_w_dt[variable %like% 's_iN2O', variable      := 's_iN2O']
      imp_w_dt[variable %like% 'm_iN2O', variable      := 'm_iN2O']
      imp_w_dt[variable %like% 's_cr_grain', variable  := 's_cr_grain']
      imp_w_dt[variable %like% 'm_cr_grain', variable  := 'm_cr_grain']
      imp_w_dt[, value := round(value, digits = 2)] # round
      setcolorder(imp_w_dt, c('cell','x', 'y', 'y_block', 'ssp', 'gcm','variable', 'value'))
      
      # SUM crop area into one layer
      crop_sum = sum(.crop_area, na.rm = TRUE) 
      names(crop_sum) = 'total_crop_area_ha'
      # TRANSFORM to dt
      crop_dt  = setDT(terra::as.data.frame(crop_sum, xy = TRUE, 
                                            cells = TRUE, na.rm =TRUE))
      # JOIN
      imp_w_dt = imp_w_dt[crop_dt, on = .(cell = cell,
                                      x    = x,
                                      y    = y)]
      imp_w_dt = imp_w_dt[!is.na(y_block)]
      imp_w_dt = imp_w_dt[total_crop_area_ha > 0,]
      
      # COMBINE
      imp_gcm_dt = rbind(imp_gcm_dt, imp_w_dt)
      gc()
    }
  }
  return(imp_gcm_dt)
}
# convert_biomass_units FUNCTION
#-----------------------------------------------------------------------------------------
# ANALYSIS #
#-----------------------------------------------------------------------------------------
range_check           = function(.var, .dt) {
  var_r     = range(.dt[, get(.var)])
  .dt_r     = data.table(variable = .var, range = c('lower', 'upper'),value = var_r)
}
regional_mean_CI      = function(.dt, .year) { # Tg
  if (.year ==  2030) {
    period = 15L
  } else if (.year == 2050) {
    period = 35L
  } else {
    print('This is not a valid year for this calculation.')
    stop()
  }
  Mg_to_Tg   = 1000000L
  dt_mCI = .dt[, .(
    GHG_m     = (mean(GHG_area)/period)/Mg_to_Tg,    # GHG mitigation potential
    GHG_sd    = (sd(GHG_area)/period)/Mg_to_Tg,
    GHG_n     = length(GHG_area),
    grain_m   = (mean(GRAIN_area)/period)/Mg_to_Tg,  # Grain yield 
    grain_sd  = (sd(GRAIN_area)/period)/Mg_to_Tg,
    grain_n   = length(GRAIN_area),
    area_m    = mean(total_crop_area_ha), # cropland area
    area_sd   = sd(total_crop_area_ha),
    area_n    = length(total_crop_area_ha)),
    by = .(ssp, y_block, IPCC_NAME)]
  dt_mCI[, GHG_ciL       := (GHG_m - (0.95*(GHG_sd))/sqrt(GHG_n))]
  dt_mCI[, GHG_ciH       := (GHG_m + (0.95*(GHG_sd))/sqrt(GHG_n))]
  dt_mCI[, grain_ciL     := (grain_m - (0.95*(grain_sd))/sqrt(grain_n))]
  dt_mCI[, grain_ciH     := (grain_m + (0.95*(grain_sd))/sqrt(grain_n))]
  dt_mCI[, area_ciL      := (area_m - (0.95*(area_sd))/sqrt(area_n))]
  dt_mCI[, area_ciH      := (area_m + (0.95*(area_sd))/sqrt(area_n))]
  
} 
global_mean_CI        = function(.dt, .year) { # Pg
  if (.year ==  2030) {
    period = 15L
  } else if (.year == 2050) {
    period = 35L
  } else {
    print('This is not a valid year for this calculation.')
    stop()
  }
  Mg_to_Pg   = 1000000000L
  dt_mCI = .dt[, .(
    GHG_m     = (mean(GHG_area)/period)/Mg_to_Pg,    # GHG mitigation potential
    GHG_sd    = (sd(GHG_area)/period)/Mg_to_Pg,
    GHG_n     = length(GHG_area),
    grain_m   = (mean(GRAIN_area)/period)/Mg_to_Pg,  # Grain yield 
    grain_sd  = (sd(GRAIN_area)/period)/Mg_to_Pg,
    grain_n   = length(GRAIN_area),
    area_m    = mean(total_crop_area_ha), # cropland area
    area_sd   = sd(total_crop_area_ha),
    area_n    = length(total_crop_area_ha)),
    by = .(ssp, y_block)]
  dt_mCI[, GHG_ciL       := (GHG_m - (0.95*(GHG_sd))/sqrt(GHG_n))]
  dt_mCI[, GHG_ciH       := (GHG_m + (0.95*(GHG_sd))/sqrt(GHG_n))]
  dt_mCI[, grain_ciL     := (grain_m - (0.95*(grain_sd))/sqrt(grain_n))]
  dt_mCI[, grain_ciH     := (grain_m + (0.95*(grain_sd))/sqrt(grain_n))]
  dt_mCI[, area_ciL      := (area_m - (0.95*(area_sd))/sqrt(area_n))]
  dt_mCI[, area_ciH      := (area_m + (0.95*(area_sd))/sqrt(area_n))]
  
} 
global_decadal_mean   = function(.dt) {
  .dt_s = copy(.dt)
  
  # 2. Global emissions
  r_relative_flux_dt = .dt_s[, .(
    m_SOC       = round((mean(m_SOC)), digits = 2),
    m_N2O       = round((mean(m_N2O)), digits = 2),
    m_dN2O      = round((mean(m_dN2O)), digits = 2),
    m_iN2O      = round((mean(m_iN2O)), digits = 2),
    m_GHG       = round((mean(m_GHG)), digits = 2),
    m_cr_grain  = round((mean(m_cr_grain)), digits = 2), # biomass
    # m_cr_residC = round((mean(m_cr_residC)), digits = 2),
    # m_cc_shootC = round((mean(m_cc_shootC)), digits = 2),
    s_SOC       = round((mean(s_SOC)), digits = 2),
    s_N2O       = round((mean(s_N2O)), digits = 2),
    s_dN2O      = round((mean(s_dN2O)), digits = 2),
    s_iN2O      = round((mean(s_iN2O)), digits = 2),
    s_GHG       = round((mean(s_GHG)), digits = 2),
    s_cr_grain  = round((mean(s_cr_grain)), digits = 2) # biomass
    # s_cr_residC = round((mean(s_cr_residC)), digits = 2),
    # s_cc_shootC = round((mean(s_cc_shootC)), digits = 2)
    ),
    by = .(gcm, ssp, y_block)]
  return(r_relative_flux_dt)
}
global_decadal_mCI    = function(.dt) {
  dt_mCI = .dt[, .(
    m_SOC_m   = mean(m_SOC),
    m_SOC_sd  = sd(m_SOC),
    m_SOC_n   = length(m_SOC),
    m_N2O_m   = mean(m_N2O),
    m_N2O_sd  = sd(m_N2O),
    m_N2O_n   = length(m_N2O),
    m_dN2O_m  = mean(m_dN2O),
    m_dN2O_sd = sd(m_dN2O),
    m_dN2O_n  = length(m_dN2O),
    m_iN2O_m    = mean(m_iN2O),
    m_iN2O_sd   = sd(m_iN2O),
    m_iN2O_n    = length(m_iN2O),
    m_GHG_m     = mean(m_GHG),
    m_GHG_sd    = sd(m_GHG),
    m_GHG_n     = length(m_GHG),
    m_grain_m   = mean(m_cr_grain), # biomass
    m_grain_sd  = sd(m_cr_grain),
    m_grain_n   = length(m_cr_grain),
    # m_residC_m  = mean(m_cr_residC),
    # m_residC_sd = sd(m_cr_residC),
    # m_residC_n  = length(m_cr_residC),
    # m_cc_shootC_m   = mean(m_cc_shootC),
    # m_cc_shootC_sd  = sd(m_cc_shootC),
    # m_cc_shootC_n   = length(m_cc_shootC),
    s_SOC_m   = mean(s_SOC),
    s_SOC_sd  = sd(s_SOC),
    s_SOC_n   = length(s_SOC),
    s_N2O_m   = mean(s_N2O),
    s_N2O_sd  = sd(s_N2O),
    s_N2O_n   = length(s_N2O),
    s_dN2O_m  = mean(s_dN2O),
    s_dN2O_sd = sd(s_dN2O),
    s_dN2O_n  = length(s_dN2O),
    s_iN2O_m    = mean(s_iN2O),
    s_iN2O_sd   = sd(s_iN2O),
    s_iN2O_n    = length(s_iN2O),
    s_GHG_m     = mean(s_GHG),
    s_GHG_sd    = sd(s_GHG),
    s_GHG_n     = length(s_GHG),
    s_grain_m   = mean(s_cr_grain), # biomass
    s_grain_sd  = sd(s_cr_grain),
    s_grain_n   = length(s_cr_grain)
    # s_residC_m  = mean(s_cr_residC),
    # s_residC_sd = sd(s_cr_residC),
    # s_residC_n  = length(s_cr_residC),
    # s_cc_shootC_m   = mean(s_cc_shootC),
    # s_cc_shootC_sd  = sd(s_cc_shootC),
    # s_cc_shootC_n   = length(s_cc_shootC)
    ),
    by = .(ssp, y_block)]
  dt_mCI[, Ts_SOC_ciL       := (s_SOC_m - (0.95*(s_SOC_sd))/sqrt(s_SOC_n))]
  dt_mCI[, Ts_SOC_ciH       := (s_SOC_m + (0.95*(s_SOC_sd))/sqrt(s_SOC_n))]
  dt_mCI[, Ts_N2O_ciL       := (s_N2O_m - (0.95*(s_N2O_sd))/sqrt(s_N2O_n))]
  dt_mCI[, Ts_N2O_ciH       := (s_N2O_m + (0.95*(s_N2O_sd))/sqrt(s_N2O_n))]
  dt_mCI[, Ts_dN2O_ciL      := (s_dN2O_m - (0.95*(s_dN2O_sd))/sqrt(s_dN2O_n))]
  dt_mCI[, Ts_dN2O_ciH      := (s_dN2O_m + (0.95*(s_dN2O_sd))/sqrt(s_dN2O_n))]
  dt_mCI[, Ts_iN2O_ciL      := (s_iN2O_m - (0.95*(s_iN2O_sd))/sqrt(s_iN2O_n))]
  dt_mCI[, Ts_iN2O_ciH      := (s_iN2O_m + (0.95*(s_iN2O_sd))/sqrt(s_iN2O_n))]
  dt_mCI[, Ts_GHG_ciL       := (s_GHG_m - (0.95*(s_GHG_sd))/sqrt(s_GHG_n))]
  dt_mCI[, Ts_GHG_ciH       := (s_GHG_m + (0.95*(s_GHG_sd))/sqrt(s_GHG_n))]
  dt_mCI[, Ts_grain_ciL       := (s_grain_m - (0.95*(s_grain_sd))/sqrt(s_grain_n))]
  dt_mCI[, Ts_grain_ciH       := (s_grain_m + (0.95*(s_grain_sd))/sqrt(s_grain_n))]
  # dt_mCI[, Ts_residC_ciL       := (s_residC_m - (0.95*(s_residC_sd))/sqrt(s_residC_n))]
  # dt_mCI[, Ts_residC_ciH       := (s_residC_m + (0.95*(s_residC_sd))/sqrt(s_residC_n))]
  # dt_mCI[, Ts_cc_shootC_ciL       := (s_cc_shootC_m - (0.95*(s_cc_shootC_sd))/sqrt(s_cc_shootC_n))]
  # dt_mCI[, Ts_cc_shootC_ciH       := (s_cc_shootC_m + (0.95*(s_cc_shootC_sd))/sqrt(s_cc_shootC_n))]
  dt_mCI[, Tm_SOC_ciL       := (m_SOC_m - (0.95*(m_SOC_sd))/sqrt(m_SOC_n))]
  dt_mCI[, Tm_SOC_ciH       := (m_SOC_m + (0.95*(m_SOC_sd))/sqrt(m_SOC_n))]
  dt_mCI[, Tm_N2O_ciL       := (m_N2O_m - (0.95*(m_N2O_sd))/sqrt(m_N2O_n))]
  dt_mCI[, Tm_N2O_ciH       := (m_N2O_m + (0.95*(m_N2O_sd))/sqrt(m_N2O_n))]
  dt_mCI[, Tm_dN2O_ciL      := (m_dN2O_m - (0.95*(m_dN2O_sd))/sqrt(m_dN2O_n))]
  dt_mCI[, Tm_dN2O_ciH      := (m_dN2O_m + (0.95*(m_dN2O_sd))/sqrt(m_dN2O_n))]
  dt_mCI[, Tm_iN2O_ciL      := (m_iN2O_m - (0.95*(m_iN2O_sd))/sqrt(m_iN2O_n))]
  dt_mCI[, Tm_iN2O_ciH      := (m_iN2O_m + (0.95*(m_iN2O_sd))/sqrt(m_iN2O_n))]
  dt_mCI[, Tm_GHG_ciL       := (m_GHG_m - (0.95*(m_GHG_sd))/sqrt(m_GHG_n))]
  dt_mCI[, Tm_GHG_ciH       := (m_GHG_m + (0.95*(m_GHG_sd))/sqrt(m_GHG_n))]
  dt_mCI[, Tm_grain_ciL       := (m_grain_m - (0.95*(m_grain_sd))/sqrt(m_grain_n))]
  dt_mCI[, Tm_grain_ciH       := (m_grain_m + (0.95*(m_grain_sd))/sqrt(m_grain_n))]
  # dt_mCI[, Tm_residC_ciL       := (m_residC_m - (0.95*(m_residC_sd))/sqrt(m_residC_n))]
  # dt_mCI[, Tm_residC_ciH       := (m_residC_m + (0.95*(m_residC_sd))/sqrt(m_residC_n))]
  # dt_mCI[, Tm_cc_shootC_ciL       := (m_cc_shootC_m - (0.95*(m_cc_shootC_sd))/sqrt(m_cc_shootC_n))]
  # dt_mCI[, Tm_cc_shootC_ciH       := (m_cc_shootC_m + (0.95*(m_cc_shootC_sd))/sqrt(m_cc_shootC_n))]
}
regional_decadal_mean = function(.dt) {
  .dt_s = copy(.dt)
  
  # 2. Global emissions
  r_relative_flux_dt = .dt_s[, .(
    m_SOC       = round((mean(m_SOC)), digits = 2),
    m_N2O       = round((mean(m_N2O)), digits = 2),
    m_dN2O      = round((mean(m_dN2O)), digits = 2),
    m_iN2O      = round((mean(m_iN2O)), digits = 2),
    m_GHG       = round((mean(m_GHG)), digits = 2),
    m_cr_grain  = round((mean(m_cr_grain)), digits = 2), # biomass
    # m_cr_residC = round((mean(m_cr_residC)), digits = 2),
    # m_cc_shootC = round((mean(m_cc_shootC)), digits = 2),
    s_SOC       = round((mean(s_SOC)), digits = 2),
    s_N2O       = round((mean(s_N2O)), digits = 2),
    s_dN2O      = round((mean(s_dN2O)), digits = 2),
    s_iN2O      = round((mean(s_iN2O)), digits = 2),
    s_GHG       = round((mean(s_GHG)), digits = 2),
    s_cr_grain  = round((mean(s_cr_grain)), digits = 2) # biomass
    # s_cr_residC = round((mean(s_cr_residC)), digits = 2),
    # s_cc_shootC = round((mean(s_cc_shootC)), digits = 2)
    ),
    by = .(gcm, ssp, y_block, IPCC_NAME)]
  return(r_relative_flux_dt)
}
regional_decadal_mCI  = function(.dt) {
  dt_mCI = .dt[, .(
    m_SOC_m   = mean(m_SOC),
    m_SOC_sd  = sd(m_SOC),
    m_SOC_n   = length(m_SOC),
    m_N2O_m   = mean(m_N2O),
    m_N2O_sd  = sd(m_N2O),
    m_N2O_n   = length(m_N2O),
    m_dN2O_m  = mean(m_dN2O),
    m_dN2O_sd = sd(m_dN2O),
    m_dN2O_n  = length(m_dN2O),
    m_iN2O_m    = mean(m_iN2O),
    m_iN2O_sd   = sd(m_iN2O),
    m_iN2O_n    = length(m_iN2O),
    m_GHG_m     = mean(m_GHG),
    m_GHG_sd    = sd(m_GHG),
    m_GHG_n     = length(m_GHG),
    m_grain_m   = mean(m_cr_grain), # biomass
    m_grain_sd  = sd(m_cr_grain),
    m_grain_n   = length(m_cr_grain),
    # m_residC_m  = mean(m_cr_residC),
    # m_residC_sd = sd(m_cr_residC),
    # m_residC_n  = length(m_cr_residC),
    # m_cc_shootC_m   = mean(m_cc_shootC),
    # m_cc_shootC_sd  = sd(m_cc_shootC),
    # m_cc_shootC_n   = length(m_cc_shootC),
    s_SOC_m   = mean(s_SOC),
    s_SOC_sd  = sd(s_SOC),
    s_SOC_n   = length(s_SOC),
    s_N2O_m   = mean(s_N2O),
    s_N2O_sd  = sd(s_N2O),
    s_N2O_n   = length(s_N2O),
    s_dN2O_m  = mean(s_dN2O),
    s_dN2O_sd = sd(s_dN2O),
    s_dN2O_n  = length(s_dN2O),
    s_iN2O_m    = mean(s_iN2O),
    s_iN2O_sd   = sd(s_iN2O),
    s_iN2O_n    = length(s_iN2O),
    s_GHG_m     = mean(s_GHG),
    s_GHG_sd    = sd(s_GHG),
    s_GHG_n     = length(s_GHG),
    s_grain_m   = mean(s_cr_grain), # biomass
    s_grain_sd  = sd(s_cr_grain),
    s_grain_n   = length(s_cr_grain)
    # s_residC_m  = mean(s_cr_residC),
    # s_residC_sd = sd(s_cr_residC),
    # s_residC_n  = length(s_cr_residC),
    # s_cc_shootC_m   = mean(s_cc_shootC),
    # s_cc_shootC_sd  = sd(s_cc_shootC),
    # s_cc_shootC_n   = length(s_cc_shootC)
    ),
    by = .(ssp, y_block, IPCC_NAME)]
  dt_mCI[, Ts_SOC_ciL       := (s_SOC_m - (0.95*(s_SOC_sd))/sqrt(s_SOC_n))]
  dt_mCI[, Ts_SOC_ciH       := (s_SOC_m + (0.95*(s_SOC_sd))/sqrt(s_SOC_n))]
  dt_mCI[, Ts_N2O_ciL       := (s_N2O_m - (0.95*(s_N2O_sd))/sqrt(s_N2O_n))]
  dt_mCI[, Ts_N2O_ciH       := (s_N2O_m + (0.95*(s_N2O_sd))/sqrt(s_N2O_n))]
  dt_mCI[, Ts_dN2O_ciL      := (s_dN2O_m - (0.95*(s_dN2O_sd))/sqrt(s_dN2O_n))]
  dt_mCI[, Ts_dN2O_ciH      := (s_dN2O_m + (0.95*(s_dN2O_sd))/sqrt(s_dN2O_n))]
  dt_mCI[, Ts_iN2O_ciL      := (s_iN2O_m - (0.95*(s_iN2O_sd))/sqrt(s_iN2O_n))]
  dt_mCI[, Ts_iN2O_ciH      := (s_iN2O_m + (0.95*(s_iN2O_sd))/sqrt(s_iN2O_n))]
  dt_mCI[, Ts_GHG_ciL       := (s_GHG_m - (0.95*(s_GHG_sd))/sqrt(s_GHG_n))]
  dt_mCI[, Ts_GHG_ciH       := (s_GHG_m + (0.95*(s_GHG_sd))/sqrt(s_GHG_n))]
  dt_mCI[, Ts_grain_ciL       := (s_grain_m - (0.95*(s_grain_sd))/sqrt(s_grain_n))]
  dt_mCI[, Ts_grain_ciH       := (s_grain_m + (0.95*(s_grain_sd))/sqrt(s_grain_n))]
  # dt_mCI[, Ts_residC_ciL       := (s_residC_m - (0.95*(s_residC_sd))/sqrt(s_residC_n))]
  # dt_mCI[, Ts_residC_ciH       := (s_residC_m + (0.95*(s_residC_sd))/sqrt(s_residC_n))]
  # dt_mCI[, Ts_cc_shootC_ciL       := (s_cc_shootC_m - (0.95*(s_cc_shootC_sd))/sqrt(s_cc_shootC_n))]
  # dt_mCI[, Ts_cc_shootC_ciH       := (s_cc_shootC_m + (0.95*(s_cc_shootC_sd))/sqrt(s_cc_shootC_n))]
  dt_mCI[, Tm_SOC_ciL       := (m_SOC_m - (0.95*(m_SOC_sd))/sqrt(m_SOC_n))]
  dt_mCI[, Tm_SOC_ciH       := (m_SOC_m + (0.95*(m_SOC_sd))/sqrt(m_SOC_n))]
  dt_mCI[, Tm_N2O_ciL       := (m_N2O_m - (0.95*(m_N2O_sd))/sqrt(m_N2O_n))]
  dt_mCI[, Tm_N2O_ciH       := (m_N2O_m + (0.95*(m_N2O_sd))/sqrt(m_N2O_n))]
  dt_mCI[, Tm_dN2O_ciL      := (m_dN2O_m - (0.95*(m_dN2O_sd))/sqrt(m_dN2O_n))]
  dt_mCI[, Tm_dN2O_ciH      := (m_dN2O_m + (0.95*(m_dN2O_sd))/sqrt(m_dN2O_n))]
  dt_mCI[, Tm_iN2O_ciL      := (m_iN2O_m - (0.95*(m_iN2O_sd))/sqrt(m_iN2O_n))]
  dt_mCI[, Tm_iN2O_ciH      := (m_iN2O_m + (0.95*(m_iN2O_sd))/sqrt(m_iN2O_n))]
  dt_mCI[, Tm_GHG_ciL       := (m_GHG_m - (0.95*(m_GHG_sd))/sqrt(m_GHG_n))]
  dt_mCI[, Tm_GHG_ciH       := (m_GHG_m + (0.95*(m_GHG_sd))/sqrt(m_GHG_n))]
  dt_mCI[, Tm_grain_ciL       := (m_grain_m - (0.95*(m_grain_sd))/sqrt(m_grain_n))]
  dt_mCI[, Tm_grain_ciH       := (m_grain_m + (0.95*(m_grain_sd))/sqrt(m_grain_n))]
  # dt_mCI[, Tm_residC_ciL       := (m_residC_m - (0.95*(m_residC_sd))/sqrt(m_residC_n))]
  # dt_mCI[, Tm_residC_ciH       := (m_residC_m + (0.95*(m_residC_sd))/sqrt(m_residC_n))]
  # dt_mCI[, Tm_cc_shootC_ciL       := (m_cc_shootC_m - (0.95*(m_cc_shootC_sd))/sqrt(m_cc_shootC_n))]
  # dt_mCI[, Tm_cc_shootC_ciH       := (m_cc_shootC_m + (0.95*(m_cc_shootC_sd))/sqrt(m_cc_shootC_n))]
}

# check still used
regional_decadal_mSE  = function(.dt) {
  dt_mSE = .dt[, .(
    m_SOC_m   = mean(m_SOC),
    m_SOC_sd  = sd(m_SOC),
    m_SOC_n   = length(m_SOC),
    m_N2O_m   = mean(m_N2O),
    m_N2O_sd  = sd(m_N2O),
    m_N2O_n   = length(m_N2O),
    m_dN2O_m  = mean(m_dN2O),
    m_dN2O_sd = sd(m_dN2O),
    m_dN2O_n  = length(m_dN2O),
    m_iN2O_m    = mean(m_iN2O),
    m_iN2O_sd   = sd(m_iN2O),
    m_iN2O_n    = length(m_iN2O),
    m_GHG_m     = mean(m_GHG),
    m_GHG_sd    = sd(m_GHG),
    m_GHG_n     = length(m_GHG),
    m_grain_m   = mean(m_cr_grain), # biomass
    m_grain_sd  = sd(m_cr_grain),
    m_grain_n   = length(m_cr_grain),
    m_residC_m  = mean(m_cr_residC),
    m_residC_sd = sd(m_cr_residC),
    m_residC_n  = length(m_cr_residC),
    m_cc_shootC_m   = mean(m_cc_shootC),
    m_cc_shootC_sd  = sd(m_cc_shootC),
    m_cc_shootC_n   = length(m_cc_shootC),
    s_SOC_m   = mean(s_SOC),
    s_SOC_sd  = sd(s_SOC),
    s_SOC_n   = length(s_SOC),
    s_N2O_m   = mean(s_N2O),
    s_N2O_sd  = sd(s_N2O),
    s_N2O_n   = length(s_N2O),
    s_dN2O_m  = mean(s_dN2O),
    s_dN2O_sd = sd(s_dN2O),
    s_dN2O_n  = length(s_dN2O),
    s_iN2O_m    = mean(s_iN2O),
    s_iN2O_sd   = sd(s_iN2O),
    s_iN2O_n    = length(s_iN2O),
    s_GHG_m     = mean(s_GHG),
    s_GHG_sd    = sd(s_GHG),
    s_GHG_n     = length(s_GHG),
    s_grain_m   = mean(s_cr_grain), # biomass
    s_grain_sd  = sd(s_cr_grain),
    s_grain_n   = length(s_cr_grain),
    s_residC_m  = mean(s_cr_residC),
    s_residC_sd = sd(s_cr_residC),
    s_residC_n  = length(s_cr_residC),
    s_cc_shootC_m   = mean(s_cc_shootC),
    s_cc_shootC_sd  = sd(s_cc_shootC),
    s_cc_shootC_n   = length(s_cc_shootC)),
    by = .(ssp, y_block, IPCC_NAME)]
  dt_mSE[, Ts_SOC_se       := s_SOC_sd/sqrt(s_SOC_n)]
  dt_mSE[, Ts_N2O_se       := s_N2O_sd/sqrt(s_N2O_n)]
  dt_mSE[, Ts_dN2O_se      := s_dN2O_sd/sqrt(s_dN2O_n)]
  dt_mSE[, Ts_iN2O_se      := s_iN2O_sd/sqrt(s_iN2O_n)]
  dt_mSE[, Ts_GHG_se       := s_GHG_sd/sqrt(s_GHG_n)]
  dt_mSE[, Ts_grain_se     := s_grain_sd/sqrt(s_grain_n)]
  dt_mSE[, Ts_residC_se    := s_residC_sd/sqrt(s_residC_n)]
  dt_mSE[, Ts_cc_shootC_se := s_cc_shootC_sd/sqrt(s_cc_shootC_n)]
  dt_mSE[, Tm_SOC_se       := m_SOC_sd/sqrt(m_SOC_n)]
  dt_mSE[, Tm_N2O_se       := m_N2O_sd/sqrt(m_N2O_n)]
  dt_mSE[, Tm_dN2O_se      := m_dN2O_sd/sqrt(m_dN2O_n)]
  dt_mSE[, Tm_iN2O_se      := m_iN2O_sd/sqrt(m_iN2O_n)]
  dt_mSE[, Tm_GHG_se       := m_GHG_sd/sqrt(m_GHG_n)]
  dt_mSE[, Tm_grain_se     := m_grain_sd/sqrt(m_grain_n)]
  dt_mSE[, Tm_residC_se    := m_residC_sd/sqrt(m_residC_n)]
  dt_mSE[, Tm_cc_shootC_se := m_cc_shootC_sd/sqrt(m_cc_shootC_n)]
}
country_decadal_mean  = function(.dt) {
  .dt_s = copy(.dt)
  
  # 2. Global emissions
  r_relative_flux_dt = .dt_s[, .(
    m_SOC       = round((mean(m_SOC)), digits = 2),
    m_N2O       = round((mean(m_N2O)), digits = 2),
    m_dN2O      = round((mean(m_dN2O)), digits = 2),
    m_iN2O      = round((mean(m_iN2O)), digits = 2),
    m_GHG       = round((mean(m_GHG)), digits = 2),
    m_cr_grain  = round((mean(m_cr_grain)), digits = 2), # biomass
    m_cr_residC = round((mean(m_cr_residC)), digits = 2),
    m_cc_shootC = round((mean(m_cc_shootC)), digits = 2),
    s_SOC       = round((mean(s_SOC)), digits = 2),
    s_N2O       = round((mean(s_N2O)), digits = 2),
    s_dN2O      = round((mean(s_dN2O)), digits = 2),
    s_iN2O      = round((mean(s_iN2O)), digits = 2),
    s_GHG       = round((mean(s_GHG)), digits = 2),
    s_cr_grain  = round((mean(s_cr_grain)), digits = 2), # biomass
    s_cr_residC = round((mean(s_cr_residC)), digits = 2),
    s_cc_shootC = round((mean(s_cc_shootC)), digits = 2)),
    by = .(gcm, ssp, y_block, WB_NAME, IPCC_NAME, crop_sum_2015)]
  return(r_relative_flux_dt)
}
country_decadal_mCI   = function(.dt) {
  dt_mCI = .dt[, .(
    m_SOC_m   = mean(m_SOC),
    m_SOC_sd  = sd(m_SOC),
    m_SOC_n   = length(m_SOC),
    m_N2O_m   = mean(m_N2O),
    m_N2O_sd  = sd(m_N2O),
    m_N2O_n   = length(m_N2O),
    m_dN2O_m  = mean(m_dN2O),
    m_dN2O_sd = sd(m_dN2O),
    m_dN2O_n  = length(m_dN2O),
    m_iN2O_m    = mean(m_iN2O),
    m_iN2O_sd   = sd(m_iN2O),
    m_iN2O_n    = length(m_iN2O),
    m_GHG_m     = mean(m_GHG),
    m_GHG_sd    = sd(m_GHG),
    m_GHG_n     = length(m_GHG),
    m_grain_m   = mean(m_cr_grain), # biomass
    m_grain_sd  = sd(m_cr_grain),
    m_grain_n   = length(m_cr_grain),
    m_residC_m  = mean(m_cr_residC),
    m_residC_sd = sd(m_cr_residC),
    m_residC_n  = length(m_cr_residC),
    m_cc_shootC_m   = mean(m_cc_shootC),
    m_cc_shootC_sd  = sd(m_cc_shootC),
    m_cc_shootC_n   = length(m_cc_shootC),
    s_SOC_m   = mean(s_SOC),
    s_SOC_sd  = sd(s_SOC),
    s_SOC_n   = length(s_SOC),
    s_N2O_m   = mean(s_N2O),
    s_N2O_sd  = sd(s_N2O),
    s_N2O_n   = length(s_N2O),
    s_dN2O_m  = mean(s_dN2O),
    s_dN2O_sd = sd(s_dN2O),
    s_dN2O_n  = length(s_dN2O),
    s_iN2O_m    = mean(s_iN2O),
    s_iN2O_sd   = sd(s_iN2O),
    s_iN2O_n    = length(s_iN2O),
    s_GHG_m     = mean(s_GHG),
    s_GHG_sd    = sd(s_GHG),
    s_GHG_n     = length(s_GHG),
    s_grain_m   = mean(s_cr_grain), # biomass
    s_grain_sd  = sd(s_cr_grain),
    s_grain_n   = length(s_cr_grain),
    s_residC_m  = mean(s_cr_residC),
    s_residC_sd = sd(s_cr_residC),
    s_residC_n  = length(s_cr_residC),
    s_cc_shootC_m   = mean(s_cc_shootC),
    s_cc_shootC_sd  = sd(s_cc_shootC),
    s_cc_shootC_n   = length(s_cc_shootC)),
    by = .(ssp, y_block, WB_NAME)]
  dt_mCI[, Ts_SOC_ciL       := (s_SOC_m - (0.95*(s_SOC_sd))/sqrt(s_SOC_n))]
  dt_mCI[, Ts_SOC_ciH       := (s_SOC_m + (0.95*(s_SOC_sd))/sqrt(s_SOC_n))]
  dt_mCI[, Ts_N2O_ciL       := (s_N2O_m - (0.95*(s_N2O_sd))/sqrt(s_N2O_n))]
  dt_mCI[, Ts_N2O_ciH       := (s_N2O_m + (0.95*(s_N2O_sd))/sqrt(s_N2O_n))]
  dt_mCI[, Ts_dN2O_ciL      := (s_dN2O_m - (0.95*(s_dN2O_sd))/sqrt(s_dN2O_n))]
  dt_mCI[, Ts_dN2O_ciH      := (s_dN2O_m + (0.95*(s_dN2O_sd))/sqrt(s_dN2O_n))]
  dt_mCI[, Ts_iN2O_ciL      := (s_iN2O_m - (0.95*(s_iN2O_sd))/sqrt(s_iN2O_n))]
  dt_mCI[, Ts_iN2O_ciH      := (s_iN2O_m + (0.95*(s_iN2O_sd))/sqrt(s_iN2O_n))]
  dt_mCI[, Ts_GHG_ciL       := (s_GHG_m - (0.95*(s_GHG_sd))/sqrt(s_GHG_n))]
  dt_mCI[, Ts_GHG_ciH       := (s_GHG_m + (0.95*(s_GHG_sd))/sqrt(s_GHG_n))]
  dt_mCI[, Ts_grain_ciL       := (s_grain_m - (0.95*(s_grain_sd))/sqrt(s_grain_n))]
  dt_mCI[, Ts_grain_ciH       := (s_grain_m + (0.95*(s_grain_sd))/sqrt(s_grain_n))]
  dt_mCI[, Ts_residC_ciL       := (s_residC_m - (0.95*(s_residC_sd))/sqrt(s_residC_n))]
  dt_mCI[, Ts_residC_ciH       := (s_residC_m + (0.95*(s_residC_sd))/sqrt(s_residC_n))]
  dt_mCI[, Ts_cc_shootC_ciL       := (s_cc_shootC_m - (0.95*(s_cc_shootC_sd))/sqrt(s_cc_shootC_n))]
  dt_mCI[, Ts_cc_shootC_ciH       := (s_cc_shootC_m + (0.95*(s_cc_shootC_sd))/sqrt(s_cc_shootC_n))]
  dt_mCI[, Tm_SOC_ciL       := (m_SOC_m - (0.95*(m_SOC_sd))/sqrt(m_SOC_n))]
  dt_mCI[, Tm_SOC_ciH       := (m_SOC_m + (0.95*(m_SOC_sd))/sqrt(m_SOC_n))]
  dt_mCI[, Tm_N2O_ciL       := (m_N2O_m - (0.95*(m_N2O_sd))/sqrt(m_N2O_n))]
  dt_mCI[, Tm_N2O_ciH       := (m_N2O_m + (0.95*(m_N2O_sd))/sqrt(m_N2O_n))]
  dt_mCI[, Tm_dN2O_ciL      := (m_dN2O_m - (0.95*(m_dN2O_sd))/sqrt(m_dN2O_n))]
  dt_mCI[, Tm_dN2O_ciH      := (m_dN2O_m + (0.95*(m_dN2O_sd))/sqrt(m_dN2O_n))]
  dt_mCI[, Tm_iN2O_ciL      := (m_iN2O_m - (0.95*(m_iN2O_sd))/sqrt(m_iN2O_n))]
  dt_mCI[, Tm_iN2O_ciH      := (m_iN2O_m + (0.95*(m_iN2O_sd))/sqrt(m_iN2O_n))]
  dt_mCI[, Tm_GHG_ciL       := (m_GHG_m - (0.95*(m_GHG_sd))/sqrt(m_GHG_n))]
  dt_mCI[, Tm_GHG_ciH       := (m_GHG_m + (0.95*(m_GHG_sd))/sqrt(m_GHG_n))]
  dt_mCI[, Tm_grain_ciL       := (m_grain_m - (0.95*(m_grain_sd))/sqrt(m_grain_n))]
  dt_mCI[, Tm_grain_ciH       := (m_grain_m + (0.95*(m_grain_sd))/sqrt(m_grain_n))]
  dt_mCI[, Tm_residC_ciL       := (m_residC_m - (0.95*(m_residC_sd))/sqrt(m_residC_n))]
  dt_mCI[, Tm_residC_ciH       := (m_residC_m + (0.95*(m_residC_sd))/sqrt(m_residC_n))]
  dt_mCI[, Tm_cc_shootC_ciL       := (m_cc_shootC_m - (0.95*(m_cc_shootC_sd))/sqrt(m_cc_shootC_n))]
  dt_mCI[, Tm_cc_shootC_ciH       := (m_cc_shootC_m + (0.95*(m_cc_shootC_sd))/sqrt(m_cc_shootC_n))]
}
#-----------------------------------------------------------------------------------------
# FIGURES #
#-----------------------------------------------------------------------------------------
IPCC_map          = function(.lu_path, .shp_f, .raster) {
  require(ggthemes)
  require(maptools)
  require(RColorBrewer)
  require(sf)
  require(terra)
  
  country.sf    = st_read(paste(.lu_path, .shp_f, sep = '/'))
  country.sf_dt = setDT(as.data.frame(country.sf))
  # CREATE raster
  shp_r       = rast(ext(country.sf), nrow = 360, ncol = 720)
  # CREATE shp as raster
  country_r   = terra::rasterize(country.sf, shp_r, fun = 'sum', "OBJECTID")
  # MATCH resolution of simulation data, dimensions the same
  target.r    = rast(nrow = 360, ncol = 720, resolution = 0.5)
  country_r   = resample(country_r, target.r, method = "near")
  country_r   = focal(country_r, w=9, fun = "modal", na.policy = "only", na.rm = TRUE) # needed to capture all gridid
  names(country_r) = "OBJECTID"
  # CREATE data.frame, merge
  country_r.dt    = as.data.frame(country_r, cells=TRUE, xy=TRUE)
  country_r.dt    = setDT(country_r.dt)
  country_n       = data.table(WB_NAME = country.sf_dt$WB_NAME, ID = country.sf_dt$OBJECTID,
                               WB_REGION = country.sf_dt$WB_REGION)
  # BIND to cell numbers
  country_r.dt = country_r.dt[country_n, on = .(OBJECTID = ID)]
  # JOIN with crop area table
  crop_area_r           = rast(paste(.lu_path, .raster, sep = '/'))
  crop_area_r$cell_area = cellSize(crop_area_r, mask=FALSE, lyrs=FALSE, unit="ha")
  crop_area_dt          = as.data.table(terra::as.data.frame(crop_area_r, xy = TRUE, cells = TRUE))
  
  crop_area_dt = crop_area_dt[country_r.dt[, .(cell, WB_NAME, WB_REGION)], on = .(cell = cell)]
  gc()
  setorder(crop_area_dt, cell)
  crop_area_dt = crop_area_dt[!is.na(x),]
  gc()
  
  # NA to 0
  crop_area_dt[is.na(maize_rainfed_2015), maize_rainfed_2015 := 0]
  crop_area_dt[is.na(maize_irrigated_2015), maize_irrigated_2015 := 0]
  crop_area_dt[is.na(soybean_rainfed_2015), soybean_rainfed_2015 := 0]
  crop_area_dt[is.na(soybean_irrigated_2015), soybean_irrigated_2015 := 0]
  crop_area_dt[is.na(wheat_rainfed_2015), wheat_rainfed_2015 := 0]
  crop_area_dt[is.na(wheat_irrigated_2015), wheat_irrigated_2015 := 0]
  
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
  crop_area_dt[WB_NAME %in% AME, IPCC_NAME := 'AME']
  crop_area_dt[WB_NAME %in% ADP, IPCC_NAME := 'ADP']
  crop_area_dt[WB_NAME %in% DEV, IPCC_NAME := 'DEV']
  crop_area_dt[WB_NAME %in% EEWCA, IPCC_NAME := 'EEWCA']
  crop_area_dt[WB_NAME %in% LAC, IPCC_NAME := 'LAC']
  
  crop_area_dt = unique(crop_area_dt[, c('WB_NAME', 'IPCC_NAME')])
  
  # WB_countries_Admin0_10m.shp
  country.sf    = st_read(paste(.lu_path, .shp_f, sep = '/'))
  country.sf    = st_transform(country.sf, crs = 4326)
  country.sf    = st_wrap_dateline(country.sf)
  country.sf    = st_transform(country.sf, crs = "+proj=eck4 +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs")
  country.sf    = merge(country.sf, crop_area_dt[, c("WB_NAME", "IPCC_NAME")], 
                        by.x = 'WB_NAME', by.y = 'WB_NAME', all.x = TRUE, sort = FALSE)
  WB_colors     = c('DEV' = c("#1B9E77"), 'AME' = c("#D95F02"), 'EEWCA' = c("#7570B3"),
                    'ADP' = c("#E7298A"), 'LAC' = c("#66A61E"))
  map = ggplot() +
    geom_sf() +
    geom_sf(data = country.sf, aes(alpha = 0.9, fill = IPCC_NAME, colour = IPCC_NAME),
            colour = "grey95", size = 0.2) +
    theme_map() +
    scale_fill_manual(values = WB_colors, na.value = "grey10") +
    theme(legend.position='none')
  map
  
  gg_maps = list(IPCC = map)
  return(gg_maps)
}
gl_bar_potential  = function(.dt_biophys, .dt_yield, .ssp) {
  scenario_lbl = c('ntill-res' = c('No-tillage'), 'ccg-res' = c('Grass CC'),
                   'ccl-res' = c('Legume CC'), 'ccg-ntill' = c('Grass CC + No-tillage'),
                   'ccl-ntill' = c('Legume CC + No-tillage'))
  # biophysical table
  .dt_biophys   = .dt_biophys[scenario %in% c('ntill-res', 'ccg-res', 'ccg-ntill', 'ccl-res', 'ccl-ntill')]
  .dt_biophys$scenario = factor(.dt_biophys$scenario, levels = c('ntill-res', 'ccg-res', 'ccg-ntill', 'ccl-res', 'ccl-ntill'))
  .dt_biophys[, type := 'biophysical']
  
  # yield-constrained table
  .dt_yield   = .dt_yield[scenario %in% c('ntill-res', 'ccg-res', 'ccg-ntill', 'ccl-res', 'ccl-ntill')]
  .dt_yield$scenario = factor(.dt_yield$scenario, levels = c('ntill-res', 'ccg-res', 'ccg-ntill', 'ccl-res', 'ccl-ntill'))
  .dt_yield[, type := 'yield_cst']
  
  # combine
  .dt = rbind(.dt_biophys, .dt_yield)
  
  plot_ghg = ggplot(data = .dt[ssp %in% .ssp], aes(x = scenario, 
                                                   y = abs(GHG_m), group = type, fill = type)) +
    geom_bar(stat = 'identity', position=position_dodge(), width = 0.8) +
    coord_flip() +
    geom_hline(yintercept = 0) +
    scale_x_discrete(labels = scenario_lbl, limits = c('ccg-res','ccl-res', 'ntill-res',
                                                       'ccg-ntill', 'ccl-ntill')) +
    ylab((expression(atop(paste(Annual~GHG~Mitigation~Potential), '('*Pg~CO[2]*-eq*~yr^-1*')')))) +
    scale_y_continuous(limits = c(0,1.5), breaks = seq(0,1.5, 0.25)) +
    scale_fill_manual(name = "Scenario", labels  = c('Biophysical', 'Yield Constrained'),
                        values = c('black', 'grey')) +
    theme_bw() +
    theme(text = element_text(color = 'black', size = 20),
          axis.text    = element_text(size = 18, color = 'black'),
          strip.text   = element_text(size = 18, color = 'black'),
          axis.title.y = element_blank(),
          axis.title.x = element_text(size = 20),
          panel.background = element_rect(fill='transparent'), 
          plot.background  = element_rect(fill='transparent', color=NA),
          legend.position  = 'none')
  
  ggplots = list(GHG = plot_ghg)
  return(ggplots)
}
rg_bar_potential  = function(.dt_biophys, .dt_yield, .ssp, .region) {
  scenario_lbl = c('ntill-res' = c('No-tillage'), 'ccg-res' = c('Grass CC'),
                   'ccl-res' = c('Legume CC'), 'ccg-ntill' = c('Grass CC + No-tillage'),
                   'ccl-ntill' = c('Legume CC + No-tillage'))
  # biophysical table
  .dt_biophys   = .dt_biophys[scenario %in% c('ntill-res', 'ccg-res', 'ccg-ntill', 'ccl-res', 'ccl-ntill')]
  .dt_biophys$scenario = factor(.dt_biophys$scenario, levels = c('ntill-res', 'ccg-res', 'ccg-ntill', 'ccl-res', 'ccl-ntill'))
  .dt_biophys[, type := 'biophysical']
  
  # yield-constrained table
  .dt_yield   = .dt_yield[scenario %in% c('ntill-res', 'ccg-res', 'ccg-ntill', 'ccl-res', 'ccl-ntill')]
  .dt_yield$scenario = factor(.dt_yield$scenario, levels = c('ntill-res', 'ccg-res', 'ccg-ntill', 'ccl-res', 'ccl-ntill'))
  .dt_yield[, type := 'yield_cst']
  
  # combine
  .dt = rbind(.dt_biophys, .dt_yield)
  
  plot_ghg = ggplot(data = .dt[ssp %in% .ssp & IPCC_NAME %in% .region], aes(x = scenario, 
                                                   y = abs(GHG_m), group = type, fill = type)) +
    geom_bar(stat = 'identity', position=position_dodge(), width = 0.8) +
    coord_flip() +
    geom_hline(yintercept = 0) +
    scale_x_discrete(labels = scenario_lbl, limits = c('ccg-res','ccl-res', 'ntill-res',
                                                       'ccg-ntill', 'ccl-ntill')) +
    ylab((expression(atop(paste(Annual~GHG~Mitigation~Potential), '('*Tg~CO[2]*-eq*~yr^-1*')')))) +
    scale_y_continuous(limits = c(0,500), breaks = seq(0,500, 100)) +
    scale_fill_manual(name = "Scenario", labels  = c('Biophysical', 'Yield Constrained'),
                      values = c('black', 'grey')) +
    theme_bw() +
    theme(text = element_text(color = 'black', size = 20),
          axis.text    = element_text(size = 18, color = 'black'),
          strip.text   = element_text(size = 18, color = 'black'),
          axis.title.y = element_blank(),
          axis.title.x = element_text(size = 20),
          panel.background = element_rect(fill='transparent'), 
          plot.background  = element_rect(fill='transparent', color=NA),
          legend.position  = 'none')
  
  ggplots = list(GHG = plot_ghg)
  return(ggplots)
}
annual_map        = function(dt_r, .time, .ssp) {
  require(ggthemes)
  require(gridExtra)
  require(maptools)
  require(RColorBrewer)
  require(scales)
  require(sf)
  require(terra)
  
  # divide by years
  if (.time == 2050) {
    period = 35L
  } else if (.time == 2030) {
    period = 15L
  } else {
    period = 85L
  }
  
  # filter dt
  dt_r = dt_r[y_block %in% .time & ssp %in% .ssp,]
  dt_r = as.data.frame(dt_r, xy = TRUE)
  setDT(dt_r)
  
  # modify dt
  Mg_to_kg = 1000L
  dt_r[, s_GHG_ann   := s_GHG/period]
  dt_r[, s_N2O_ann   := s_N2O/period]
  dt_r[, s_SOC_ann   := s_SOC/period]
  dt_r[, s_grain_ann := (s_cr_grain/period)*Mg_to_kg]

  # create manual breaks
  quants.GHG  = dt_r[, quantile(s_GHG_ann,seq.int(0,1,length.out = 6), na.rm = TRUE)]
  quants.N2O  = dt_r[, quantile(s_N2O_ann,seq.int(0,1,length.out = 6), na.rm = TRUE)]
  quants.SOC  = dt_r[, quantile(s_SOC_ann,seq.int(0,1,length.out = 6), na.rm = TRUE)]
  quants.gr   = dt_r[, quantile(s_grain_ann,seq.int(0,1,length.out = 6), na.rm = TRUE)]
  # make manual quantiles
  ghg.quants = c(-Inf,-3,-2,-1,-0.5,0,0.5,1,2,3,Inf)
  dt_r[, quant.cuts.GHG := cut(s_GHG_ann, breaks = ghg.quants)]
  n2o.quants = c(-Inf,-3,-2,-1,-0.5,0,0.5,1,2,3,Inf)
  dt_r[, quant.cuts.N2O := cut(s_N2O_ann, breaks = n2o.quants)]
  soc.quants = c(-Inf,-3,-2,-1,-0.5,0,0.5,1,2,3,Inf)
  dt_r[, quant.cuts.SOC := cut(s_SOC_ann, breaks = soc.quants)]
  gr.quants = c(-Inf,-750,-500,-250,-100,0,100,250,500,750,Inf)
  dt_r[, quant.cuts.gr := cut(s_grain_ann, breaks = gr.quants)]
  # create raster
  dt_r = as.data.frame(dt_r, xy = TRUE)
  r                = rast(nrow = 360, ncol = 720, nlyr = 4, xmin = -180, xmax = 180, ymin = -90, ymax = 90)
  crs(r)           = "epsg:4326"
  # original projection
  r_lat_ghg = rast(res = 0.5, nlyr = 4, extent = ext(r), crs = crs(r)) 
  r_lat_ghg[[1]][dt_r$cell] = dt_r$quant.cuts.GHG
  r_lat_ghg[[2]][dt_r$cell] = dt_r$quant.cuts.N2O
  r_lat_ghg[[3]][dt_r$cell] = dt_r$quant.cuts.SOC
  r_lat_ghg[[4]][dt_r$cell] = dt_r$quant.cuts.gr
  names(r_lat_ghg) = c('s_GHG_2100', 's_N2O_2100','s_SOC_2100', 's_gr_2100')
  # equal area projection
  newcrs = "+proj=eck4 +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"
  # get new output dimensions
  project(r_lat_ghg, newcrs)
  # dimensions  : 569, 1138, 3  (nrow, ncol, nlyr)
  # resolution  : 29727.52, 29727.52  (x, y)
  # extent      : -16921203, 16908719, -8454359, 8460601  (xmin, xmax, ymin, ymax)
  r_newcrs   = rast(ncols = 1138, nrows = 569, nlyr = 4, xmin = -16921203, xmax = 16908719, ymin = -8454359, ymax = 8460601, crs = newcrs)
  r_eckiv    = project(r_lat_ghg, r_newcrs)
  
  r_eckiv_dt = as.data.frame(r_eckiv, cells = TRUE, xy = TRUE)
  r_eckiv_dt = setDT(r_eckiv_dt)
  r_eckiv_dt_GHG = r_eckiv_dt[,.(cell, x, y, s_GHG_2100, s_N2O_2100,s_SOC_2100, s_gr_2100)]
  r_eckiv_dt_GHG[, s_GHG_2100 := round(s_GHG_2100, digits = 0)]
  r_eckiv_dt_GHG[, s_GHG_2100 := as.character(s_GHG_2100)]
  r_eckiv_dt_GHG[, s_SOC_2100 := round(s_SOC_2100, digits = 0)]
  r_eckiv_dt_GHG[, s_SOC_2100 := as.character(s_SOC_2100)]
  r_eckiv_dt_GHG[, s_N2O_2100 := round(s_N2O_2100, digits = 0)]
  r_eckiv_dt_GHG[, s_N2O_2100 := as.character(s_N2O_2100)]
  r_eckiv_dt_GHG[, s_gr_2100 := round(s_gr_2100, digits = 0)]
  r_eckiv_dt_GHG[, s_gr_2100 := as.character(s_gr_2100)]
  
  # create sf object
  data(wrld_simpl)
  wrld_simpl_sf = sf::st_as_sf(wrld_simpl)
  wrld_simpl_sf_eckiv = st_transform(wrld_simpl_sf, crs = newcrs)
  wrld_simpl_sf_eckiv = wrld_simpl_sf_eckiv[wrld_simpl_sf_eckiv$NAME != 'Antarctica',]
  
  small_islands       = wrld_simpl_sf[wrld_simpl_sf$AREA < 10000,]
  remove              = c('Antigua and Barbuda', 'American Samoa', 'Barbados', 'Bermuda',
                          'Bahamas', 'Solomon Islands', 'Cayman Islands', 'Comoros','Cook Islands', 'Cape Verde',
                          'Dominica', 'Fiji','Falkland Islands (Malvinas)', 'Micronesia, Federated States of', 'Grenada',
                          'New Caledonia', 'Niue', 'Anguilla','French Polynesia', 'Guam', 'Kiribati', 'Martinique','Maldives', 'Aruba', 'Northern Mariana Islands',
                          'Faroe Islands', 'Mayotte', 'Mauritius','Aaland Islands', 'Norfolk Island', 'Cocos (Keeling) Islands',
                          'Bouvet Island', 'French Southern and Antarctic Lands', 'Heard Island and McDonald Islands',
                          'British Indian Ocean Territory', 'Christmas Island', 'Vanuatu','United States Minor Outlying Islands',
                          'Nauru', 'Reunion', 'Saint Kitts and Nevis', 'Seychelles', 'Saint Lucia', 'Tokelau', 'Tonga',
                          'Tuvalu','Saint Vincent and the Grenadines', 'British Virgin Islands', 'United States Virgin Islands',
                          'Wallis and Futuna Islands', 'Samoa', 'Guadeloupe', 'Netherlands Antilles', 'Pitcairn Islands','
                          Palau', 'Marshall Islands', 'Saint Pierre and Miquelon', 'Saint Helena', 'San Marino',
                          'Turks and Caicos Islands', 'Svalbard', 'Saint Martin', 'Saint Barthelemy', 'South Georgia South Sandwich Islands',
                          'Guernsey', 'Jersey')
  wrld_simpl_sf_eckiv = wrld_simpl_sf_eckiv[!wrld_simpl_sf_eckiv$NAME %in% remove,]
  colors1 = c('1' = "#276419", '2' = "#4D9221",'3' = "#7FBC41", '4'= "#B8E186", 
              '5' = "#E6F5D0", '6' = "#FDE0EF",'7' = "#F1B6DA", '8' = "#DE77AE", '9' = "#C51B7D", '10'= "#8E0152")
  colors2 = c('1' = "#8E0152", '2' = "#C51B7D",'3' = "#DE77AE", '4'= "#F1B6DA", 
              '5' = "#FDE0EF", '6' = "#E6F5D0",'7' = "#B8E186", '8' = "#7FBC41", '9' = "#4D9221", '10'= "#276419")
  
  # create maps
  gg_GHG = ggplot() + 
    geom_sf(data = wrld_simpl_sf_eckiv, fill = "grey75",
            colour = "grey45", size = 0.2) +
    theme_map() +
    geom_tile(data = r_eckiv_dt_GHG[,.(cell, x, y, s_GHG_2100)],
              aes(x = x, y = y, fill = s_GHG_2100)) +
    scale_fill_manual(values = colors1) +
    theme(legend.position='none',
          plot.margin = unit(c(0,0,-1,0), "cm"))
  gg_GHG
  
  gg_N2O = ggplot() + 
    geom_sf(data = wrld_simpl_sf_eckiv, fill = "grey75",
            colour = "grey45", size = 0.2) +
    theme_map() +
    geom_tile(data = r_eckiv_dt_GHG[,.(cell, x, y, s_N2O_2100)],
              aes(x = x, y = y, fill = s_N2O_2100)) +
    scale_fill_manual(values = colors1) +
    theme(legend.position='none',
          plot.margin = unit(c(0,0,-1,0), "cm"))
  gg_N2O
  
  gg_SOC = ggplot() + 
    geom_sf(data = wrld_simpl_sf_eckiv, fill = "grey75",
            colour = "grey45", size = 0.2) +
    theme_map() +
    geom_tile(data = r_eckiv_dt_GHG[,.(cell, x, y, s_SOC_2100)],
              aes(x = x, y = y, fill = s_SOC_2100)) +
    scale_fill_manual(values = colors1) +
    theme(legend.position='none',
          plot.margin = unit(c(0,0,-1,0), "cm"))
  gg_SOC
  
  gg_gr = ggplot() + 
    geom_sf(data = wrld_simpl_sf_eckiv, fill = "grey75",
            colour = "grey45", size = 0.2) +
    theme_map() +
    geom_tile(data = r_eckiv_dt_GHG[,.(cell, x, y, s_gr_2100)],
              aes(x = x, y = y, fill = s_gr_2100)) +
    scale_fill_manual(values = colors2) +
    theme(legend.position='none',
          plot.margin = unit(c(0,0,-1,0), "cm"))
  gg_gr
  
  colors_bar1 = c("#276419", "#4D9221","#7FBC41","#B8E186","#E6F5D0",
                           "#FDE0EF","#F1B6DA","#DE77AE", "#C51B7D", "#8E0152")
                           
  colors_bar2 =c("#8E0152", "#C51B7D","#DE77AE", "#F1B6DA", 
                          "#FDE0EF", "#E6F5D0","#B8E186", "#7FBC41", "#4D9221", "#276419")
                          
  # create legend
  gg_legend1 = plot_discrete_cbar(c(-Inf,-3,-2,-1,-0.5,0,0.5,1,2,3,Inf),
                                  colors = colors_bar1,
                                  legend_title = expression(Soil~GHG~Emissions~Difference~'('~Mg~CO[2]*-eq~ha^-1*~yr^-1*')'),
                                  spacing = 'constant',
                                  font_size = 7)
  gg_legend1
  gg_legend2 = plot_discrete_cbar(c(-Inf,-750,-500,-250,-100,0,100,250,500,750,Inf),
                                  colors = colors_bar2,
                                  legend_title = expression(Yield~Difference~'('~kg~ha^-1*~yr^-1*')'),
                                  spacing = 'constant',
                                  font_size = 7)
  gg_legend2
  
  gg_maps = list(N2O = gg_N2O, SOC = gg_SOC, GHG = gg_GHG, Yield = gg_gr, legend1 = gg_legend1, legend2 = gg_legend2)
  return(gg_maps)
}
soc_gl_ha         = function(dt_gcm, .ssp, .time) {
  
  # divide by years
  if (.time == 2050) {
    period = 35L
  } else if (.time == 2030) {
    period = 15L
  } else {
    period = 85L
  }
  #update scenario names
  scenario_lbl = c('ntill-res' = c('No-tillage'), 'ccg-res' = c('Grass CC'),
                   'ccl-res' = c('Legume CC'), 'ccg-ntill' = c('Grass CC + No-tillage'),
                   'ccl-ntill' = c('Legume CC + No-tillage'))
  
  dt_gcm[, y_block := as.factor(y_block)]
  dt_gcm = dt_gcm[scenario %in% c('ntill-res', 'ccg-res', 'ccl-res',
                         'ccg-ntill', 'ccl-ntill')]

  bplot_soc = ggplot(data = dt_gcm[ssp %in% .ssp & y_block %in% .time],
                     aes(x = (s_SOC/period), y = scenario, fill = scenario)) +
    geom_density_ridges(bandwidth = 0.015) +
    geom_vline(xintercept = 0) +
    ylab('Scenario') +
    scale_y_discrete(labels = scenario_lbl, limits = c('ccg-res','ntill-res','ccl-res','ccg-ntill', 'ccl-ntill')) +
    scale_fill_manual(values = c('ccl-ntill' = "#66C2A5", 'ccg-ntill' = "#FC8D62", 'ccl-res' = "#8DA0CB",
                                 'ntill-res' = "#E78AC3", 'ccg-res' = "#A6D854")) +
    xlab(expression(atop(paste(Annual~SOC~Sequestration), '('*Mg~CO[2]*-eq*~ha^-1*~yr^-1*')'))) +
    scale_x_continuous(limits = c(-3, 0), breaks = seq(-3, 0, 0.5)) +
    theme_bw() +
    theme(text = element_text(color = 'black', size = 20),
          axis.text    = element_text(size = 18, color = 'black'),
          strip.text   = element_text(size = 18, color = 'black'),
          axis.title.x = element_text(size = 20),
          axis.title.y = element_blank(),
          legend.position = 'none')
  
  return(bplot_soc)
}
soc_rg_ha         = function(dt_gcm, .ssp, .time) {
  
  # divide by years
  if (.time == 2050) {
    period = 35L
  } else if (.time == 2030) {
    period = 15L
  } else {
    period = 85L
  }
  #update scenario names
  scenario_lbl = c('ntill-res' = c('No-tillage'), 'ccg-res' = c('Grass CC'),
                   'ccl-res' = c('Legume CC'), 'ccg-ntill' = c('Grass CC + No-tillage'),
                   'ccl-ntill' = c('Legume CC + No-tillage'))
  dt_gcm[, y_block := as.factor(y_block)]
  dt_gcm = dt_gcm[scenario %in% c('ntill-res', 'ccg-res', 'ccl-res',
                                  'ccg-ntill', 'ccl-ntill')]
  
  bplot_soc = ggplot(data = dt_gcm[ssp %in% .ssp & y_block %in% .time],
                     aes(x = (s_SOC/period), y = scenario, fill = scenario)) +
    geom_density_ridges(bandwidth = 0.015) +
    geom_vline(xintercept = 0) +
    facet_wrap(.~IPCC_NAME, nrow = 3, ncol = 2) +
    ylab('Scenario') +
    scale_y_discrete(labels = scenario_lbl, limits = c('ccg-res','ntill-res','ccl-res','ccg-ntill', 'ccl-ntill')) +
    scale_fill_manual(values = c('ccl-ntill' = "#66C2A5", 'ccg-ntill' = "#FC8D62", 'ccl-res' = "#8DA0CB",
                                 'ntill-res' = "#E78AC3", 'ccg-res' = "#A6D854")) +
    xlab(expression(atop(paste(Annual~SOC~Sequestration), '('*Mg~CO[2]*-eq*~ha^-1*~yr^-1*')'))) +
    scale_x_continuous(limits = c(-3.5, 0), breaks = seq(-3.5, 0, 0.5)) +
    theme_bw() +
    theme(text = element_text(color = 'black', size = 20),
          axis.text    = element_text(size = 16, color = 'black'),
          strip.text   = element_text(size = 16, color = 'black'),
          axis.title.x = element_text(size = 20),
          axis.title.y = element_blank(),
          legend.position = 'none')
  
  return(bplot_soc)
}
bmp_map           = function(dt_r, .ssp) {
  require(ggthemes)
  require(gridExtra)
  require(maptools)
  require(RColorBrewer)
  require(scales)
  require(sf)
  require(terra)
  
  # filter dt
  dt_r = dt_r[ssp %in% .ssp,]
  dt_r = as.data.frame(dt_r, xy = TRUE)
  setDT(dt_r)
  
  # update scenarios
  dt_r[scenario %in% 'BAU',       scenario := '1']
  dt_r[scenario %in% 'res',       scenario := '2']
  dt_r[scenario %in% 'ntill',     scenario := '3']
  dt_r[scenario %in% 'ntill-res', scenario := '4']
  dt_r[scenario %in% 'ccg',       scenario := '5']
  dt_r[scenario %in% 'ccg-res',   scenario := '6']
  dt_r[scenario %in% 'ccl',       scenario := '7']
  dt_r[scenario %in% 'ccl-res',   scenario := '8']
  dt_r[scenario %in% 'ccg-ntill', scenario := '9']
  dt_r[scenario %in% 'ccl-ntill', scenario := '10']
  dt_r[, scenario := as.numeric(scenario)]
  
  # create raster
  dt_r = as.data.frame(dt_r, xy = TRUE)
  r                = rast(nrow = 360, ncol = 720, nlyr = 1, xmin = -180, xmax = 180, ymin = -90, ymax = 90)
  crs(r)           = "epsg:4326"
  # original projection
  r_lat_ghg = rast(res = 0.5, nlyr = 1, extent = ext(r), crs = crs(r)) 
  r_lat_ghg[[1]][dt_r$cell] = dt_r$scenario

  names(r_lat_ghg) = c('scenario')
  # equal area projection
  newcrs = "+proj=eck4 +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"
  # get new output dimensions
  project(r_lat_ghg, newcrs)
  # dimensions  : 569, 1138, 1  (nrow, ncol, nlyr)
  # resolution  : 29727.52, 29727.52  (x, y)
  # extent      : -16921203, 16908719, -8454359, 8460601  (xmin, xmax, ymin, ymax)
  r_newcrs   = rast(ncols = 1138, nrows = 569, nlyr = 1, xmin = -16921203, xmax = 16908719, ymin = -8454359, ymax = 8460601, crs = newcrs)
  r_eckiv    = project(r_lat_ghg, r_newcrs)
  
  r_eckiv_dt = as.data.frame(r_eckiv, cells = TRUE, xy = TRUE)
  r_eckiv_dt = setDT(r_eckiv_dt)
  r_eckiv_dt_GHG = r_eckiv_dt[,.(cell, x, y, scenario)]
  r_eckiv_dt_GHG[, scenario := round(scenario, digits = 0)]
  r_eckiv_dt_GHG[, scenario := as.character(scenario)]
  
  # create sf object
  data(wrld_simpl)
  wrld_simpl_sf = sf::st_as_sf(wrld_simpl)
  wrld_simpl_sf_eckiv = st_transform(wrld_simpl_sf, crs = newcrs)
  wrld_simpl_sf_eckiv = wrld_simpl_sf_eckiv[wrld_simpl_sf_eckiv$NAME != 'Antarctica',]
  
  small_islands       = wrld_simpl_sf[wrld_simpl_sf$AREA < 10000,]
  remove              = c('Antigua and Barbuda', 'American Samoa', 'Barbados', 'Bermuda',
                          'Bahamas', 'Solomon Islands', 'Cayman Islands', 'Comoros','Cook Islands', 'Cape Verde',
                          'Dominica', 'Fiji','Falkland Islands (Malvinas)', 'Micronesia, Federated States of', 'Grenada',
                          'New Caledonia', 'Niue', 'Anguilla','French Polynesia', 'Guam', 'Kiribati', 'Martinique','Maldives', 'Aruba', 'Northern Mariana Islands',
                          'Faroe Islands', 'Mayotte', 'Mauritius','Aaland Islands', 'Norfolk Island', 'Cocos (Keeling) Islands',
                          'Bouvet Island', 'French Southern and Antarctic Lands', 'Heard Island and McDonald Islands',
                          'British Indian Ocean Territory', 'Christmas Island', 'Vanuatu','United States Minor Outlying Islands',
                          'Nauru', 'Reunion', 'Saint Kitts and Nevis', 'Seychelles', 'Saint Lucia', 'Tokelau', 'Tonga',
                          'Tuvalu','Saint Vincent and the Grenadines', 'British Virgin Islands', 'United States Virgin Islands',
                          'Wallis and Futuna Islands', 'Samoa', 'Guadeloupe', 'Netherlands Antilles', 'Pitcairn Islands','
                          Palau', 'Marshall Islands', 'Saint Pierre and Miquelon', 'Saint Helena', 'San Marino',
                          'Turks and Caicos Islands', 'Svalbard', 'Saint Martin', 'Saint Barthelemy', 'South Georgia South Sandwich Islands',
                          'Guernsey', 'Jersey')
  wrld_simpl_sf_eckiv = wrld_simpl_sf_eckiv[!wrld_simpl_sf_eckiv$NAME %in% remove,]
  colors = c('1' = "#7A0403FF", '2' = "#C82803FF",'3' = "#F36215FF", '4'= "#FCB036FF", 
             '5' = "#D6E635FF", '6' = "#8EFF49FF",'7' = "#2CF09EFF", '8' = "#23C3E4FF", 
             '9' = "#4681F7FF", '10'= "#3F3994FF")
    # create maps
  gg_BMP = ggplot() + 
    geom_sf(data = wrld_simpl_sf_eckiv, fill = "grey75",
            colour = "grey45", size = 0.2) +
    theme_map() +
    geom_tile(data = r_eckiv_dt_GHG[,.(cell, x, y, scenario)],
              aes(x = x, y = y, fill = scenario)) +
    scale_fill_manual(values = colors) +
    theme(legend.position='none',
          plot.margin = unit(c(0,0,-1,0), "cm"))
  gg_BMP
  
  gg_maps = list(bmp = gg_BMP)
  return(gg_maps)
}


# check still needed
cumul_map                = function(dt_r, .time, .ssp) {
  require(ggthemes)
  require(gridExtra)
  require(maptools)
  require(RColorBrewer)
  require(scales)
  require(sf)
  require(terra)
  
  # filter dt
  dt_r = dt_r[y_block %in% .time & ssp %in% .ssp,]
  dt_r = as.data.frame(dt_r, xy = TRUE)
  setnames(dt_r, 'gridid', 'cell')
  setDT(dt_r)
  
  # create manual breaks
  quants.GHG  = dt_r[, quantile(s_GHG,seq.int(0,1,length.out = 6), na.rm = TRUE)]
  quants.dN2O = dt_r[, quantile(s_dN2O,seq.int(0,1,length.out = 6), na.rm = TRUE)]
  quants.iN2O = dt_r[, quantile(s_iN2O,seq.int(0,1,length.out = 6), na.rm = TRUE)]
  quants.SOC  = dt_r[, quantile(s_SOC,seq.int(0,1,length.out = 6), na.rm = TRUE)]
  quants.gr   = dt_r[, quantile(s_cr_grain,seq.int(0,1,length.out = 6), na.rm = TRUE)]
  # make manual quantiles
  ghg.quants = c(-Inf,-50,-25,-10,-5,0,5,10,25,50,Inf)
  dt_r[, quant.cuts.GHG := cut(s_GHG, breaks = ghg.quants)]
  dn2o.quants = c(-Inf,-50,-25,-10,-5,0,5,10,25,50,Inf)
  dt_r[, quant.cuts.dN2O := cut(s_dN2O, breaks = dn2o.quants)]
  in2o.quants = c(-Inf,-50,-25,-10,-5,0,5,10,25,50,Inf)
  dt_r[, quant.cuts.iN2O := cut(s_iN2O, breaks = in2o.quants)]
  soc.quants = c(-Inf,-50,-25,-10,-5,0,5,10,25,50,Inf)
  dt_r[, quant.cuts.SOC := cut(s_SOC, breaks = soc.quants)]
  gr.quants = c(-Inf,-75,-50,-25,-10,0,10,25,50,75,Inf)
  dt_r[, quant.cuts.gr := cut(s_cr_grain, breaks = gr.quants)]
  # create raster
  dt_r = as.data.frame(dt_r, xy = TRUE)
  r                = rast(nrow = 360, ncol = 720, nlyr = 5, xmin = -180, xmax = 180, ymin = -90, ymax = 90)
  crs(r)           = "epsg:4326"
  # original projection
  r_lat_ghg = rast(res = 0.5, nlyr = 5, extent = ext(r), crs = crs(r)) 
  r_lat_ghg[[1]][dt_r$cell] = dt_r$quant.cuts.GHG
  r_lat_ghg[[2]][dt_r$cell] = dt_r$quant.cuts.dN2O
  r_lat_ghg[[3]][dt_r$cell] = dt_r$quant.cuts.iN2O
  r_lat_ghg[[4]][dt_r$cell] = dt_r$quant.cuts.SOC
  r_lat_ghg[[5]][dt_r$cell] = dt_r$quant.cuts.gr
  names(r_lat_ghg) = c('s_GHG_2100', 's_dN2O_2100', 's_iN2O_2100','s_SOC_2100', 's_gr_2100')
  # equal area projection
  newcrs = "+proj=eck4 +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"
  # get new output dimensions
  project(r_lat_ghg, newcrs)
  # dimensions  : 569, 1138, 3  (nrow, ncol, nlyr)
  # resolution  : 29727.52, 29727.52  (x, y)
  # extent      : -16921203, 16908719, -8454359, 8460601  (xmin, xmax, ymin, ymax)
  r_newcrs   = rast(ncols = 1138, nrows = 569, nlyr = 5, xmin = -16921203, xmax = 16908719, ymin = -8454359, ymax = 8460601, crs = newcrs)
  r_eckiv    = project(r_lat_ghg, r_newcrs)
  
  r_eckiv_dt = as.data.frame(r_eckiv, cells = TRUE, xy = TRUE)
  r_eckiv_dt = setDT(r_eckiv_dt)
  r_eckiv_dt_GHG = r_eckiv_dt[,.(cell, x, y, s_GHG_2100, s_dN2O_2100, s_iN2O_2100, s_SOC_2100, s_gr_2100)]
  r_eckiv_dt_GHG[, s_GHG_2100 := round(s_GHG_2100, digits = 0)]
  r_eckiv_dt_GHG[, s_GHG_2100 := as.character(s_GHG_2100)]
  r_eckiv_dt_GHG[, s_SOC_2100 := round(s_SOC_2100, digits = 0)]
  r_eckiv_dt_GHG[, s_SOC_2100 := as.character(s_SOC_2100)]
  r_eckiv_dt_GHG[, s_dN2O_2100 := round(s_dN2O_2100, digits = 0)]
  r_eckiv_dt_GHG[, s_dN2O_2100 := as.character(s_dN2O_2100)]
  r_eckiv_dt_GHG[, s_iN2O_2100 := round(s_iN2O_2100, digits = 0)]
  r_eckiv_dt_GHG[, s_iN2O_2100 := as.character(s_iN2O_2100)]
  r_eckiv_dt_GHG[, s_gr_2100 := round(s_gr_2100, digits = 0)]
  r_eckiv_dt_GHG[, s_gr_2100 := as.character(s_gr_2100)]
  
  # create sf object
  data(wrld_simpl)
  wrld_simpl_sf = sf::st_as_sf(wrld_simpl)
  wrld_simpl_sf_eckiv = st_transform(wrld_simpl_sf, crs = newcrs)
  wrld_simpl_sf_eckiv = wrld_simpl_sf_eckiv[wrld_simpl_sf_eckiv$NAME != 'Antarctica',]
  
  small_islands       = wrld_simpl_sf[wrld_simpl_sf$AREA < 10000,]
  remove              = c('Antigua and Barbuda', 'American Samoa', 'Barbados', 'Bermuda',
                          'Bahamas', 'Solomon Islands', 'Cayman Islands', 'Comoros','Cook Islands', 'Cape Verde',
                          'Dominica', 'Fiji','Falkland Islands (Malvinas)', 'Micronesia, Federated States of', 'Grenada',
                          'New Caledonia', 'Niue', 'Anguilla','French Polynesia', 'Guam', 'Kiribati', 'Martinique','Maldives', 'Aruba', 'Northern Mariana Islands',
                          'Faroe Islands', 'Mayotte', 'Mauritius','Aaland Islands', 'Norfolk Island', 'Cocos (Keeling) Islands',
                          'Bouvet Island', 'French Southern and Antarctic Lands', 'Heard Island and McDonald Islands',
                          'British Indian Ocean Territory', 'Christmas Island', 'Vanuatu','United States Minor Outlying Islands',
                          'Nauru', 'Reunion', 'Saint Kitts and Nevis', 'Seychelles', 'Saint Lucia', 'Tokelau', 'Tonga',
                          'Tuvalu','Saint Vincent and the Grenadines', 'British Virgin Islands', 'United States Virgin Islands',
                          'Wallis and Futuna Islands', 'Samoa', 'Guadeloupe', 'Netherlands Antilles', 'Pitcairn Islands','
                          Palau', 'Marshall Islands', 'Saint Pierre and Miquelon', 'Saint Helena', 'San Marino',
                          'Turks and Caicos Islands', 'Svalbard', 'Saint Martin', 'Saint Barthelemy', 'South Georgia South Sandwich Islands',
                          'Guernsey', 'Jersey')
  wrld_simpl_sf_eckiv = wrld_simpl_sf_eckiv[!wrld_simpl_sf_eckiv$NAME %in% remove,]
  colors1 = c('1' = "#276419", '2' = "#4D9221",'3' = "#7FBC41", '4'= "#B8E186", 
              '5' = "#E6F5D0", '6' = "#FDE0EF",'7' = "#F1B6DA", '8' = "#DE77AE", '9' = "#C51B7D", '10'= "#8E0152")
  colors2 = c('1' = "#8E0152", '2' = "#C51B7D",'3' = "#DE77AE", '4'= "#F1B6DA", 
              '5' = "#FDE0EF", '6' = "#E6F5D0",'7' = "#B8E186", '8' = "#7FBC41", '9' = "#4D9221", '10'= "#276419")
  
  # create maps
  gg_GHG = ggplot() + 
    geom_sf(data = wrld_simpl_sf_eckiv, fill = "grey75",
            colour = "grey45", size = 0.2) +
    theme_map() +
    geom_tile(data = r_eckiv_dt_GHG[,.(cell, x, y, s_GHG_2100)],
              aes(x = x, y = y, fill = s_GHG_2100)) +
    scale_fill_manual(values = colors1) +
    theme(legend.position='none',
          plot.margin = unit(c(0,0,-1,0), "cm"))
  gg_GHG
  
  gg_dN2O = ggplot() + 
    geom_sf(data = wrld_simpl_sf_eckiv, fill = "grey75",
            colour = "grey45", size = 0.2) +
    theme_map() +
    geom_tile(data = r_eckiv_dt_GHG[,.(cell, x, y, s_dN2O_2100)],
              aes(x = x, y = y, fill = s_dN2O_2100)) +
    scale_fill_manual(values = colors1) +
    theme(legend.position='none',
          plot.margin = unit(c(0,0,-1,0), "cm"))
  gg_dN2O
  
  gg_iN2O = ggplot() + 
    geom_sf(data = wrld_simpl_sf_eckiv, fill = "grey75",
            colour = "grey45", size = 0.2) +
    theme_map() +
    geom_tile(data = r_eckiv_dt_GHG[,.(cell, x, y, s_iN2O_2100)],
              aes(x = x, y = y, fill = s_iN2O_2100)) +
    scale_fill_manual(values = colors1) +
    theme(legend.position='none',
          plot.margin = unit(c(0,0,-1,0), "cm"))
  gg_iN2O
  
  gg_SOC = ggplot() + 
    geom_sf(data = wrld_simpl_sf_eckiv, fill = "grey75",
            colour = "grey45", size = 0.2) +
    theme_map() +
    geom_tile(data = r_eckiv_dt_GHG[,.(cell, x, y, s_SOC_2100)],
              aes(x = x, y = y, fill = s_SOC_2100)) +
    scale_fill_manual(values = colors1) +
    theme(legend.position='none',
          plot.margin = unit(c(0,0,-1,0), "cm"))
  gg_SOC
  
  gg_gr = ggplot() + 
    geom_sf(data = wrld_simpl_sf_eckiv, fill = "grey75",
            colour = "grey45", size = 0.2) +
    theme_map() +
    geom_tile(data = r_eckiv_dt_GHG[,.(cell, x, y, s_gr_2100)],
              aes(x = x, y = y, fill = s_gr_2100)) +
    scale_fill_manual(values = colors2) +
    theme(legend.position='none',
          plot.margin = unit(c(0,0,-1,0), "cm"))
  gg_gr
  
  colors_bar1 = c("#276419", "#4D9221","#7FBC41","#B8E186","#E6F5D0",
                           "#FDE0EF","#F1B6DA","#DE77AE", "#C51B7D", "#8E0152")
                           
  colors_bar2 =c("#8E0152", "#C51B7D","#DE77AE", "#F1B6DA", 
                          "#FDE0EF", "#E6F5D0","#B8E186", "#7FBC41", "#4D9221", "#276419")
                          
  # create legend
  gg_legend1 = plot_discrete_cbar(c(-Inf,-50,-25,-10,-5,0,5,10,25,50,Inf),
                                  colors = colors_bar1,
                                  legend_title = expression(Soil~GHG~Emissions~Difference~'('~Mg~CO[2]*-eq~ha^-1*~')'),
                                  spacing = 'constant',
                                  font_size = 7)
  gg_legend1
  gg_legend2 = plot_discrete_cbar(c(-Inf,-75,-50,-25,-10,0,10,25,50,75,Inf),
                                  colors = colors_bar2,
                                  legend_title = expression(Yield~Difference~'('~Mg~C~ha^-1*~')'),
                                  spacing = 'constant',
                                  font_size = 7)
  gg_legend2
  
  # gg_GHG  = grid.arrange(gg_GHG, gg_legend1, nrow = 2)
  # gg_SOC  = grid.arrange(gg_SOC, gg_legend1, nrow = 2)
  # gg_dN2O = grid.arrange(gg_dN2O, gg_legend1, nrow = 2)
  # gg_iN2O = grid.arrange(gg_iN2O, gg_legend1, nrow = 2)
  # gg_grain = grid.arrange(gg_gr, gg_legend2, nrow = 2)
  
  # gg_maps = list(dN2O = gg_dN2O, iN2O = gg_iN2O,SOC = gg_SOC, GHG = gg_GHG, Yield = gg_grain, legend1 = gg_legend1, legend2 = gg_legend2)
  gg_maps = list(dN2O = gg_dN2O, iN2O = gg_iN2O,SOC = gg_SOC, GHG = gg_GHG, Yield = gg_gr, legend1 = gg_legend1, legend2 = gg_legend2)
  return(gg_maps)
}
cumul_map_st             = function(dt_r, .time, .ssp) {
  require(ggthemes)
  require(gridExtra)
  require(maptools)
  require(RColorBrewer)
  require(scales)
  require(sf)
  require(terra)
  
  # filter dt
  dt_r = dt_r[y_block %in% .time & ssp %in% .ssp,]
  dt_r = as.data.frame(dt_r, xy = TRUE)
  setnames(dt_r, 'gridid', 'cell')
  setDT(dt_r)
  
  # create manual breaks
  quants.GHG = dt_r[, quantile(s_GHG,seq.int(0,1,length.out = 6), na.rm = TRUE)]
  quants.dN2O = dt_r[, quantile(s_dN2O,seq.int(0,1,length.out = 6), na.rm = TRUE)]
  quants.iN2O = dt_r[, quantile(s_iN2O,seq.int(0,1,length.out = 6), na.rm = TRUE)]
  quants.SOC = dt_r[, quantile(s_SOC,seq.int(0,1,length.out = 6), na.rm = TRUE)]
  quants.gr   = dt_r[, quantile(s_cr_grain,seq.int(0,1,length.out = 6), na.rm = TRUE)]
  # make manual quantiles
  ghg.quants = c(-Inf,-150,-100,-50,-10,0,10,50,100,150,Inf)
  dt_r[, quant.cuts.GHG := cut(s_GHG, breaks = ghg.quants)]
  dn2o.quants = c(-Inf,-150,-100,-50,-10,0,10,50,100,150,Inf)
  dt_r[, quant.cuts.dN2O := cut(s_dN2O, breaks = dn2o.quants)]
  in2o.quants = c(-Inf,-150,-100,-50,-10,0,10,50,100,150,Inf)
  dt_r[, quant.cuts.iN2O := cut(s_iN2O, breaks = in2o.quants)]
  soc.quants = c(-Inf,-150,-100,-50,-10,0,10,50,100,150,Inf)
  dt_r[, quant.cuts.SOC := cut(s_SOC, breaks = soc.quants)]
  gr.quants = c(-Inf,-75,-50,-25,-10,0,10,25,50,75,Inf)
  dt_r[, quant.cuts.gr := cut(s_cr_grain, breaks = gr.quants)]
  # create raster
  dt_r = as.data.frame(dt_r, xy = TRUE)
  r                = rast(nrow = 360, ncol = 720, nlyr = 5, xmin = -180, xmax = 180, ymin = -90, ymax = 90)
  crs(r)           = "epsg:4326"
  # original projection
  r_lat_ghg = rast(res = 0.5, nlyr = 5, extent = ext(r), crs = crs(r)) 
  r_lat_ghg[[1]][dt_r$cell] = dt_r$quant.cuts.GHG
  r_lat_ghg[[2]][dt_r$cell] = dt_r$quant.cuts.dN2O
  r_lat_ghg[[3]][dt_r$cell] = dt_r$quant.cuts.iN2O
  r_lat_ghg[[4]][dt_r$cell] = dt_r$quant.cuts.SOC
  r_lat_ghg[[5]][dt_r$cell] = dt_r$quant.cuts.gr
  names(r_lat_ghg) = c('s_GHG_2100', 's_dN2O_2100', 's_iN2O_2100','s_SOC_2100', 's_gr_2100')
  # equal area projection
  newcrs = "+proj=eck4 +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"
  # get new output dimensions
  project(r_lat_ghg, newcrs)
  # dimensions  : 569, 1138, 3  (nrow, ncol, nlyr)
  # resolution  : 29727.52, 29727.52  (x, y)
  # extent      : -16921203, 16908719, -8454359, 8460601  (xmin, xmax, ymin, ymax)
  r_newcrs   = rast(ncols = 1138, nrows = 569, nlyr = 5, xmin = -16921203, xmax = 16908719, ymin = -8454359, ymax = 8460601, crs = newcrs)
  r_eckiv    = project(r_lat_ghg, r_newcrs)
  
  r_eckiv_dt = as.data.frame(r_eckiv, cells = TRUE, xy = TRUE)
  r_eckiv_dt = setDT(r_eckiv_dt)
  r_eckiv_dt_GHG = r_eckiv_dt[,.(cell, x, y, s_GHG_2100, s_dN2O_2100, s_iN2O_2100, s_SOC_2100, s_gr_2100)]
  r_eckiv_dt_GHG[, s_GHG_2100 := round(s_GHG_2100, digits = 0)]
  r_eckiv_dt_GHG[, s_GHG_2100 := as.character(s_GHG_2100)]
  r_eckiv_dt_GHG[, s_SOC_2100 := round(s_SOC_2100, digits = 0)]
  r_eckiv_dt_GHG[, s_SOC_2100 := as.character(s_SOC_2100)]
  r_eckiv_dt_GHG[, s_dN2O_2100 := round(s_dN2O_2100, digits = 0)]
  r_eckiv_dt_GHG[, s_dN2O_2100 := as.character(s_dN2O_2100)]
  r_eckiv_dt_GHG[, s_iN2O_2100 := round(s_iN2O_2100, digits = 0)]
  r_eckiv_dt_GHG[, s_iN2O_2100 := as.character(s_iN2O_2100)]
  r_eckiv_dt_GHG[, s_gr_2100 := round(s_gr_2100, digits = 0)]
  r_eckiv_dt_GHG[, s_gr_2100 := as.character(s_gr_2100)]
  
  
  # create sf object
  data(wrld_simpl)
  wrld_simpl_sf = sf::st_as_sf(wrld_simpl)
  wrld_simpl_sf_eckiv = st_transform(wrld_simpl_sf, crs = newcrs)
  wrld_simpl_sf_eckiv = wrld_simpl_sf_eckiv[wrld_simpl_sf_eckiv$NAME != 'Antarctica',]
  
  small_islands       = wrld_simpl_sf[wrld_simpl_sf$AREA < 10000,]
  remove              = c('Antigua and Barbuda', 'American Samoa', 'Barbados', 'Bermuda',
                          'Bahamas', 'Solomon Islands', 'Cayman Islands', 'Comoros','Cook Islands', 'Cape Verde',
                          'Dominica', 'Fiji','Falkland Islands (Malvinas)', 'Micronesia, Federated States of', 'Grenada',
                          'New Caledonia', 'Niue', 'Anguilla','French Polynesia', 'Guam', 'Kiribati', 'Martinique','Maldives', 'Aruba', 'Northern Mariana Islands',
                          'Faroe Islands', 'Mayotte', 'Mauritius','Aaland Islands', 'Norfolk Island', 'Cocos (Keeling) Islands',
                          'Bouvet Island', 'French Southern and Antarctic Lands', 'Heard Island and McDonald Islands',
                          'British Indian Ocean Territory', 'Christmas Island', 'Vanuatu','United States Minor Outlying Islands',
                          'Nauru', 'Reunion', 'Saint Kitts and Nevis', 'Seychelles', 'Saint Lucia', 'Tokelau', 'Tonga',
                          'Tuvalu','Saint Vincent and the Grenadines', 'British Virgin Islands', 'United States Virgin Islands',
                          'Wallis and Futuna Islands', 'Samoa', 'Guadeloupe', 'Netherlands Antilles', 'Pitcairn Islands','
                          Palau', 'Marshall Islands', 'Saint Pierre and Miquelon', 'Saint Helena', 'San Marino',
                          'Turks and Caicos Islands', 'Svalbard', 'Saint Martin', 'Saint Barthelemy', 'South Georgia South Sandwich Islands',
                          'Guernsey', 'Jersey')
  wrld_simpl_sf_eckiv = wrld_simpl_sf_eckiv[!wrld_simpl_sf_eckiv$NAME %in% remove,]
  colors1 = c('1' = "#276419", '2' = "#4D9221",'3' = "#7FBC41", '4'= "#B8E186", 
              '5' = "#E6F5D0", '6' = "#FDE0EF",'7' = "#F1B6DA", '8' = "#DE77AE", '9' = "#C51B7D", '10'= "#8E0152")
  colors2 = c('1' = "#8E0152", '2' = "#C51B7D",'3' = "#DE77AE", '4'= "#F1B6DA", 
              '5' = "#FDE0EF", '6' = "#E6F5D0",'7' = "#B8E186", '8' = "#7FBC41", '9' = "#4D9221", '10'= "#276419")
  
  # create maps
  gg_GHG = ggplot() + 
    geom_sf(data = wrld_simpl_sf_eckiv, fill = "grey75",
            colour = "grey45", size = 0.2) +
    theme_map() +
    geom_tile(data = r_eckiv_dt_GHG[,.(cell, x, y, s_GHG_2100)],
              aes(x = x, y = y, fill = s_GHG_2100)) +
    scale_fill_manual(values = colors1) +
    theme(legend.position='none',
          plot.margin = unit(c(0,0,-1,0), "cm"))
  gg_GHG
  
  gg_dN2O = ggplot() + 
    geom_sf(data = wrld_simpl_sf_eckiv, fill = "grey75",
            colour = "grey45", size = 0.2) +
    theme_map() +
    geom_tile(data = r_eckiv_dt_GHG[,.(cell, x, y, s_dN2O_2100)],
              aes(x = x, y = y, fill = s_dN2O_2100)) +
    scale_fill_manual(values = colors1) +
    theme(legend.position='none',
          plot.margin = unit(c(0,0,-1,0), "cm"))
  gg_dN2O
  
  gg_iN2O = ggplot() + 
    geom_sf(data = wrld_simpl_sf_eckiv, fill = "grey75",
            colour = "grey45", size = 0.2) +
    theme_map() +
    geom_tile(data = r_eckiv_dt_GHG[,.(cell, x, y, s_iN2O_2100)],
              aes(x = x, y = y, fill = s_iN2O_2100)) +
    scale_fill_manual(values = colors1) +
    theme(legend.position='none',
          plot.margin = unit(c(0,0,-1,0), "cm"))
  gg_iN2O
  
  gg_SOC = ggplot() + 
    geom_sf(data = wrld_simpl_sf_eckiv, fill = "grey75",
            colour = "grey45", size = 0.2) +
    theme_map() +
    geom_tile(data = r_eckiv_dt_GHG[,.(cell, x, y, s_SOC_2100)],
              aes(x = x, y = y, fill = s_SOC_2100)) +
    scale_fill_manual(values = colors1) +
    theme(legend.position='none',
          plot.margin = unit(c(0,0,-1,0), "cm"))
  gg_SOC
  
  gg_gr = ggplot() + 
    geom_sf(data = wrld_simpl_sf_eckiv, fill = "grey75",
            colour = "grey45", size = 0.2) +
    theme_map() +
    geom_tile(data = r_eckiv_dt_GHG[,.(cell, x, y, s_gr_2100)],
              aes(x = x, y = y, fill = s_gr_2100)) +
    scale_fill_manual(values = colors2) +
    theme(legend.position='none',
          plot.margin = unit(c(0,0,-1,0), "cm"))
  gg_gr
  
  colors_bar1 = c("#276419", "#4D9221","#7FBC41","#B8E186","#E6F5D0",
                           "#FDE0EF","#F1B6DA","#DE77AE", "#C51B7D", "#8E0152")
                           
  colors_bar2 =c("#8E0152", "#C51B7D","#DE77AE", "#F1B6DA", 
                          "#FDE0EF", "#E6F5D0","#B8E186", "#7FBC41", "#4D9221", "#276419")
                          # create legend
  gg_legend1 = plot_discrete_cbar(c(-Inf,-150,-100,-50,-10,0,10,50,100,150,Inf),
                                  colors = colors_bar1,
                                  legend_title = expression(Soil~GHG~Emissions~Difference~'('~Mg~CO[2]*-eq~ha^-1*~')'),
                                  spacing = 'constant',
                                  font_size = 7)
  gg_legend1
  gg_legend2 = plot_discrete_cbar(c(-Inf,-75,-50,-25,-10,0,10,25,50,75,Inf),
                                  colors = colors_bar2,
                                  legend_title = expression(Yield~Difference~'('~Mg~ha^-1*~')'),
                                  spacing = 'constant',
                                  font_size = 7)
  gg_legend2
  
  
  # gg_GHG  = grid.arrange(gg_GHG, gg_legend1, nrow = 2)
  # gg_SOC  = grid.arrange(gg_SOC, gg_legend1, nrow = 2)
  # gg_dN2O = grid.arrange(gg_dN2O, gg_legend1, nrow = 2)
  # gg_iN2O = grid.arrange(gg_iN2O, gg_legend1, nrow = 2)
  # gg_grain = grid.arrange(gg_gr, gg_legend2, nrow = 2)
  
  # gg_maps = list(dN2O = gg_dN2O, iN2O = gg_iN2O,SOC = gg_SOC, GHG = gg_GHG, Yield = gg_grain, legend1 = gg_legend1, legend2 = gg_legend2)
  gg_maps = list(dN2O = gg_dN2O, iN2O = gg_iN2O,SOC = gg_SOC, GHG = gg_GHG, Yield = gg_gr, legend1 = gg_legend1, legend2 = gg_legend2)
  return(gg_maps)
}
cumul_var_map            = function(dt_r, .time, .ssp) {
  require(ggthemes)
  require(gridExtra)
  require(maptools)
  require(RColorBrewer)
  require(scales)
  require(sf)
  require(terra)
  
  # filter dt
  dt_r = dt_r[ssp %in% .ssp & y_block %in% .time,]
  dt_r = as.data.frame(dt_r, xy = TRUE)
  setnames(dt_r, 'gridid', 'cell')
  setDT(dt_r)
  
  # create manual breaks
  quants.NPP   = dt_r[, quantile(s_cr_NPP,seq.int(0,1,length.out = 6), na.rm = TRUE)]
  quants.ccbio = dt_r[, quantile(s_cc_shootC,seq.int(0,1,length.out = 6), na.rm = TRUE)]
  quants.annet = dt_r[, quantile(s_annet,seq.int(0,1,length.out = 6), na.rm = TRUE)] # cm yr
  # make manual quantiles
  NPP.quants = c(-Inf,-1000,-750,-500,-100,0,100,500,750,1000,Inf)
  dt_r[, quant.cuts.NPP := cut(s_cr_NPP, breaks = NPP.quants)]
  ccbio.quants = c(-Inf,-1000,-750,-500,-100,0,100,500,750,1000,Inf)
  dt_r[, quant.cuts.ccbio := cut(s_cc_shootC, breaks = ccbio.quants)]
  annet.quants = c(-Inf,-1000,-750,-500,-100,0,100,500,750,1000,Inf)
  dt_r[, quant.cuts.annet := cut(s_annet, breaks = annet.quants)]
  # create raster
  dt_r = as.data.frame(dt_r, xy = TRUE)
  r                = rast(nrow = 360, ncol = 720, nlyr = 3, xmin = -180, xmax = 180, ymin = -90, ymax = 90)
  crs(r)           = "epsg:4326"
  # original projection
  r_lat_var = rast(res = 0.5, nlyr = 3, extent = ext(r), crs = crs(r)) 
  r_lat_var[[1]][dt_r$cell] = dt_r$quant.cuts.NPP
  r_lat_var[[2]][dt_r$cell] = dt_r$quant.cuts.ccbio
  r_lat_var[[3]][dt_r$cell] = dt_r$quant.cuts.annet
  names(r_lat_var) = c('s_NPP_2100', 's_ccbio_2100', 's_annet_2100')
  # equal area projection
  newcrs = "+proj=eck4 +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"
  # get new output dimensions
  project(r_lat_var, newcrs)
  # dimensions  : 569, 1138, 3  (nrow, ncol, nlyr)
  # resolution  : 29727.52, 29727.52  (x, y)
  # extent      : -16921203, 16908719, -8454359, 8460601  (xmin, xmax, ymin, ymax)
  r_newcrs   = rast(ncols = 1138, nrows = 569, nlyr = 3, xmin = -16921203, xmax = 16908719, ymin = -8454359, ymax = 8460601, crs = newcrs)
  r_eckiv    = project(r_lat_var, r_newcrs)
  
  r_eckiv_dt = as.data.frame(r_eckiv, cells = TRUE, xy = TRUE)
  r_eckiv_dt = setDT(r_eckiv_dt)
  r_eckiv_dt_var = r_eckiv_dt[,.(cell, x, y, s_NPP_2100, s_ccbio_2100, s_annet_2100)]
  r_eckiv_dt_var[, s_NPP_2100 := round(s_NPP_2100, digits = 0)]
  r_eckiv_dt_var[, s_NPP_2100 := as.character(s_NPP_2100)]
  r_eckiv_dt_var[, s_ccbio_2100 := round(s_ccbio_2100, digits = 0)]
  r_eckiv_dt_var[, s_ccbio_2100 := as.character(s_ccbio_2100)]
  r_eckiv_dt_var[, s_annet_2100 := round(s_annet_2100, digits = 0)]
  r_eckiv_dt_var[, s_annet_2100 := as.character(s_annet_2100)]
 
  # create sf object
  data(wrld_simpl)
  wrld_simpl_sf = sf::st_as_sf(wrld_simpl)
  wrld_simpl_sf_eckiv = st_transform(wrld_simpl_sf, crs = newcrs)
  wrld_simpl_sf_eckiv = wrld_simpl_sf_eckiv[wrld_simpl_sf_eckiv$NAME != 'Antarctica',]
  
  small_islands       = wrld_simpl_sf[wrld_simpl_sf$AREA < 10000,]
  remove              = c('Antigua and Barbuda', 'American Samoa', 'Barbados', 'Bermuda',
                          'Bahamas', 'Solomon Islands', 'Cayman Islands', 'Comoros','Cook Islands', 'Cape Verde',
                          'Dominica', 'Fiji','Falkland Islands (Malvinas)', 'Micronesia, Federated States of', 'Grenada',
                          'New Caledonia', 'Niue', 'Anguilla','French Polynesia', 'Guam', 'Kiribati', 'Martinique','Maldives', 'Aruba', 'Northern Mariana Islands',
                          'Faroe Islands', 'Mayotte', 'Mauritius','Aaland Islands', 'Norfolk Island', 'Cocos (Keeling) Islands',
                          'Bouvet Island', 'French Southern and Antarctic Lands', 'Heard Island and McDonald Islands',
                          'British Indian Ocean Territory', 'Christmas Island', 'Vanuatu','United States Minor Outlying Islands',
                          'Nauru', 'Reunion', 'Saint Kitts and Nevis', 'Seychelles', 'Saint Lucia', 'Tokelau', 'Tonga',
                          'Tuvalu','Saint Vincent and the Grenadines', 'British Virgin Islands', 'United States Virgin Islands',
                          'Wallis and Futuna Islands', 'Samoa', 'Guadeloupe', 'Netherlands Antilles', 'Pitcairn Islands','
                          Palau', 'Marshall Islands', 'Saint Pierre and Miquelon', 'Saint Helena', 'San Marino',
                          'Turks and Caicos Islands', 'Svalbard', 'Saint Martin', 'Saint Barthelemy', 'South Georgia South Sandwich Islands',
                          'Guernsey', 'Jersey')
  wrld_simpl_sf_eckiv = wrld_simpl_sf_eckiv[!wrld_simpl_sf_eckiv$NAME %in% remove,]
  colors1 = c('1' = "#276419", '2' = "#4D9221",'3' = "#7FBC41", '4'= "#B8E186", 
              '5' = "#E6F5D0", '6' = "#FDE0EF",'7' = "#F1B6DA", '8' = "#DE77AE", '9' = "#C51B7D", '10'= "#8E0152")
  colors2 = c('1' = "#8E0152", '2' = "#C51B7D",'3' = "#DE77AE", '4'= "#F1B6DA", 
              '5' = "#FDE0EF", '6' = "#E6F5D0",'7' = "#B8E186", '8' = "#7FBC41", '9' = "#4D9221", '10'= "#276419")
  
  # create maps
  gg_NPP = ggplot() + 
    geom_sf(data = wrld_simpl_sf_eckiv, fill = "grey75",
            colour = "grey45", size = 0.2) +
    theme_map() +
    geom_tile(data = r_eckiv_dt_var[,.(cell, x, y, s_NPP_2100)],
              aes(x = x, y = y, fill = s_NPP_2100)) +
    scale_fill_manual(values = colors2) +
    theme(legend.position='none',
          plot.margin = unit(c(0,0,-1,0), "cm"))
  gg_NPP
  
  gg_ccbio = ggplot() + 
    geom_sf(data = wrld_simpl_sf_eckiv, fill = "grey75",
            colour = "grey45", size = 0.2) +
    theme_map() +
    geom_tile(data = r_eckiv_dt_var[,.(cell, x, y, s_ccbio_2100)],
              aes(x = x, y = y, fill = s_ccbio_2100)) +
    scale_fill_manual(values = colors2) +
    theme(legend.position='none',
          plot.margin = unit(c(0,0,-1,0), "cm"))
  gg_ccbio
  
  gg_annet = ggplot() + 
    geom_sf(data = wrld_simpl_sf_eckiv, fill = "grey75",
            colour = "grey45", size = 0.2) +
    theme_map() +
    geom_tile(data = r_eckiv_dt_var[,.(cell, x, y, s_annet_2100)],
              aes(x = x, y = y, fill = s_annet_2100)) +
    scale_fill_manual(values = colors2) +
    theme(legend.position='none',
          plot.margin = unit(c(0,0,-1,0), "cm"))
  gg_annet
  
  # colors_bar1 = c("#276419", "#4D9221","#7FBC41","#B8E186","#E6F5D0",
  #                          "#FDE0EF","#F1B6DA","#DE77AE", "#C51B7D", "#8E0152")
                           
  colors_bar2 =c("#8E0152", "#C51B7D","#DE77AE", "#F1B6DA", 
                          "#FDE0EF", "#E6F5D0","#B8E186", "#7FBC41", "#4D9221", "#276419")
                          
  # create legend
  gg_legend1 = plot_discrete_cbar(c(-Inf,-50,-25,-10,-5,0,5,10,25,50,Inf),
                                  colors = colors_bar1,
                                  legend_title = expression(Soil~GHG~Emissions~Difference~'('~Mg~CO[2]*-eq~ha^-1*~')'),
                                  spacing = 'constant',
                                  font_size = 7)
  gg_legend1
  gg_legend2 = plot_discrete_cbar(c(-Inf,-75,-50,-25,-10,0,10,25,50,75,Inf),
                                  colors = colors_bar2,
                                  legend_title = expression(Yield~Difference~'('~Mg~C~ha^-1*~')'),
                                  spacing = 'constant',
                                  font_size = 7)
  gg_legend2
  
  # gg_GHG  = grid.arrange(gg_GHG, gg_legend1, nrow = 2)
  # gg_SOC  = grid.arrange(gg_SOC, gg_legend1, nrow = 2)
  # gg_dN2O = grid.arrange(gg_dN2O, gg_legend1, nrow = 2)
  # gg_iN2O = grid.arrange(gg_iN2O, gg_legend1, nrow = 2)
  # gg_grain = grid.arrange(gg_gr, gg_legend2, nrow = 2)
  
  # gg_maps = list(dN2O = gg_dN2O, iN2O = gg_iN2O,SOC = gg_SOC, GHG = gg_GHG, Yield = gg_grain, legend1 = gg_legend1, legend2 = gg_legend2)
  gg_maps = list(dN2O = gg_dN2O, iN2O = gg_iN2O,SOC = gg_SOC, GHG = gg_GHG, Yield = gg_gr, legend1 = gg_legend1, legend2 = gg_legend2)
  return(gg_maps)
}
ann_var_map              = function(dt_r, .time, .ssp, .vars) {
  require(ggthemes)
  require(gridExtra)
  require(maptools)
  require(RColorBrewer)
  require(scales)
  require(sf)
  require(terra)
  
  # filter dt
  dt_r   = dt_r[ssp %in% .ssp & y_block %in% .time,]
  k_vars = c('gridid', 'x', 'y', 'y_block', 'ssp', .vars)
  dt_r   = dt_r[, ..k_vars]
  # mean response
  dt_r_m = dt_r[y_block <= .time, lapply(.SD, mean), .SDcols = .vars, by = .(gridid, x, y, ssp)]
  
  # to data frame
  dt_r_m = as.data.frame(dt_r_m, xy = TRUE)
  setnames(dt_r_m, 'gridid', 'cell')
  setDT(dt_r_m)
  
  # create manual breaks
  quants.CN     = dt_r_m[, quantile(m_sC.N,seq.int(0,1,length.out = 6), na.rm = TRUE)]
  quants.sldcmp = dt_r_m[, quantile(m_sldcmp,seq.int(0,1,length.out = 6), na.rm = TRUE)]
  quants.sfdcmp = dt_r_m[, quantile(m_sfdcmp,seq.int(0,1,length.out = 6), na.rm = TRUE)] 
  quants.annet  = dt_r_m[, quantile(m_annet,seq.int(0,1,length.out = 6), na.rm = TRUE)] # cm yr
  quants.NPP  = dt_r_m[, quantile(m_cr_NPP,seq.int(0,1,length.out = 6), na.rm = TRUE)] # cm yr
  # make manual quantiles
  CN.quants = c(-Inf,-3.0,-2.0,-1.0,-0.5,0,0.5,1.0,2.0,3.0,Inf)
  dt_r_m[, quant.cuts.CN := cut(m_sC.N, breaks = CN.quants)]
  sl.quants = c(-Inf,-0.10,-0.05,-0.01,-0.005,0,0.005,0.01,0.05,0.10,Inf)
  dt_r_m[, quant.cuts.sldcmp := cut(m_sldcmp, breaks = sl.quants)]
  sf.quants = c(-Inf,-0.10,-0.05,-0.01,-0.005,0,0.005,0.01,0.05,0.10,Inf)
  dt_r_m[, quant.cuts.sfdcmp := cut(m_sfdcmp, breaks = sf.quants)]
  an.quants = c(-Inf,-20,-10,-5,-1,0,1,5,10,20,Inf)
  dt_r_m[, quant.cuts.annet := cut(m_annet, breaks = an.quants)]
  npp.quants = c(-Inf,-15000,-10000,-1000,-500,0,500,1000,10000,15000,Inf)
  dt_r_m[, quant.cuts.NPP := cut(m_cr_NPP, breaks = npp.quants)]
  # create raster
  dt_r_m = as.data.frame(dt_r_m, xy = TRUE)
  r                = rast(nrow = 360, ncol = 720, nlyr = 5, xmin = -180, xmax = 180, ymin = -90, ymax = 90)
  crs(r)           = "epsg:4326"
  # original projection
  r_lat_var = rast(res = 0.5, nlyr = 5, extent = ext(r), crs = crs(r)) 
  r_lat_var[[1]][dt_r_m$cell] = dt_r_m$quant.cuts.CN
  r_lat_var[[2]][dt_r_m$cell] = dt_r_m$quant.cuts.sldcmp
  r_lat_var[[3]][dt_r_m$cell] = dt_r_m$quant.cuts.sfdcmp
  r_lat_var[[4]][dt_r_m$cell] = dt_r_m$quant.cuts.annet
  r_lat_var[[5]][dt_r_m$cell] = dt_r_m$quant.cuts.NPP
  names(r_lat_var) = c('CN_2100', 'sldcmp_2100', 'sfdcmp_2100', 'annet_2100', 'NPP_2100')
  # equal area projection
  newcrs = "+proj=eck4 +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"
  # get new output dimensions
  project(r_lat_var, newcrs)
  # dimensions  : 569, 1138, 5  (nrow, ncol, nlyr)
  # resolution  : 29727.52, 29727.52  (x, y)
  # extent      : -16921203, 16908719, -8454359, 8460601  (xmin, xmax, ymin, ymax)
  r_newcrs   = rast(ncols = 1138, nrows = 569, nlyr = 5, xmin = -16921203, xmax = 16908719, ymin = -8454359, ymax = 8460601, crs = newcrs)
  r_eckiv    = project(r_lat_var, r_newcrs)
  
  r_eckiv_dt = as.data.frame(r_eckiv, cells = TRUE, xy = TRUE)
  r_eckiv_dt = setDT(r_eckiv_dt)
  r_eckiv_dt_var = r_eckiv_dt[,.(cell, x, y, CN_2100, sldcmp_2100, sfdcmp_2100, annet_2100, NPP_2100)]
  r_eckiv_dt_var[, CN_2100 := round(CN_2100, digits = 0)]
  r_eckiv_dt_var[, CN_2100 := as.character(CN_2100)]
  r_eckiv_dt_var[, sldcmp_2100 := round(sldcmp_2100, digits = 0)]
  r_eckiv_dt_var[, sldcmp_2100 := as.character(sldcmp_2100)]
  r_eckiv_dt_var[, sfdcmp_2100 := round(sfdcmp_2100, digits = 0)]
  r_eckiv_dt_var[, sfdcmp_2100 := as.character(sfdcmp_2100)]
  r_eckiv_dt_var[, annet_2100 := round(annet_2100, digits = 0)]
  r_eckiv_dt_var[, annet_2100 := as.character(annet_2100)]
  r_eckiv_dt_var[, NPP_2100 := round(NPP_2100, digits = 0)]
  r_eckiv_dt_var[, NPP_2100 := as.character(NPP_2100)]
  # create sf object
  data(wrld_simpl)
  wrld_simpl_sf = sf::st_as_sf(wrld_simpl)
  wrld_simpl_sf_eckiv = st_transform(wrld_simpl_sf, crs = newcrs)
  wrld_simpl_sf_eckiv = wrld_simpl_sf_eckiv[wrld_simpl_sf_eckiv$NAME != 'Antarctica',]
  
  small_islands       = wrld_simpl_sf[wrld_simpl_sf$AREA < 10000,]
  remove              = c('Antigua and Barbuda', 'American Samoa', 'Barbados', 'Bermuda',
                          'Bahamas', 'Solomon Islands', 'Cayman Islands', 'Comoros','Cook Islands', 'Cape Verde',
                          'Dominica', 'Fiji','Falkland Islands (Malvinas)', 'Micronesia, Federated States of', 'Grenada',
                          'New Caledonia', 'Niue', 'Anguilla','French Polynesia', 'Guam', 'Kiribati', 'Martinique','Maldives', 'Aruba', 'Northern Mariana Islands',
                          'Faroe Islands', 'Mayotte', 'Mauritius','Aaland Islands', 'Norfolk Island', 'Cocos (Keeling) Islands',
                          'Bouvet Island', 'French Southern and Antarctic Lands', 'Heard Island and McDonald Islands',
                          'British Indian Ocean Territory', 'Christmas Island', 'Vanuatu','United States Minor Outlying Islands',
                          'Nauru', 'Reunion', 'Saint Kitts and Nevis', 'Seychelles', 'Saint Lucia', 'Tokelau', 'Tonga',
                          'Tuvalu','Saint Vincent and the Grenadines', 'British Virgin Islands', 'United States Virgin Islands',
                          'Wallis and Futuna Islands', 'Samoa', 'Guadeloupe', 'Netherlands Antilles', 'Pitcairn Islands','
                          Palau', 'Marshall Islands', 'Saint Pierre and Miquelon', 'Saint Helena', 'San Marino',
                          'Turks and Caicos Islands', 'Svalbard', 'Saint Martin', 'Saint Barthelemy', 'South Georgia South Sandwich Islands',
                          'Guernsey', 'Jersey')
  wrld_simpl_sf_eckiv = wrld_simpl_sf_eckiv[!wrld_simpl_sf_eckiv$NAME %in% remove,]
  colors1 = c('1' = "#276419", '2' = "#4D9221",'3' = "#7FBC41", '4'= "#B8E186", 
              '5' = "#E6F5D0", '6' = "#FDE0EF",'7' = "#F1B6DA", '8' = "#DE77AE", '9' = "#C51B7D", '10'= "#8E0152")
  colors2 = c('1' = "#8E0152", '2' = "#C51B7D",'3' = "#DE77AE", '4'= "#F1B6DA", 
              '5' = "#FDE0EF", '6' = "#E6F5D0",'7' = "#B8E186", '8' = "#7FBC41", '9' = "#4D9221", '10'= "#276419")
  
  # create maps
  gg_CN = ggplot() + 
    geom_sf(data = wrld_simpl_sf_eckiv, fill = "grey75",
            colour = "grey45", size = 0.2) +
    theme_map() +
    geom_tile(data = r_eckiv_dt_var[,.(cell, x, y, CN_2100)],
              aes(x = x, y = y, fill = CN_2100)) +
    scale_fill_manual(values = colors2) +
    theme(legend.position='none',
          plot.margin = unit(c(0,0,-1,0), "cm"))
  gg_CN
  
  gg_sldcmp = ggplot() + 
    geom_sf(data = wrld_simpl_sf_eckiv, fill = "grey75",
            colour = "grey45", size = 0.2) +
    theme_map() +
    geom_tile(data = r_eckiv_dt_var[,.(cell, x, y, sldcmp_2100)],
              aes(x = x, y = y, fill = sldcmp_2100)) +
    scale_fill_manual(values = colors1) +
    theme(legend.position='none',
          plot.margin = unit(c(0,0,-1,0), "cm"))
  gg_sldcmp
  
  gg_sfdcmp = ggplot() + 
    geom_sf(data = wrld_simpl_sf_eckiv, fill = "grey75",
            colour = "grey45", size = 0.2) +
    theme_map() +
    geom_tile(data = r_eckiv_dt_var[,.(cell, x, y, sfdcmp_2100)],
              aes(x = x, y = y, fill = sfdcmp_2100)) +
    scale_fill_manual(values = colors1) +
    theme(legend.position='none',
          plot.margin = unit(c(0,0,-1,0), "cm"))
  gg_sfdcmp
  
  gg_annet= ggplot() + 
    geom_sf(data = wrld_simpl_sf_eckiv, fill = "grey75",
            colour = "grey45", size = 0.2) +
    theme_map() +
    geom_tile(data = r_eckiv_dt_var[,.(cell, x, y, annet_2100)],
              aes(x = x, y = y, fill = annet_2100)) +
    scale_fill_manual(values = colors2) +
    theme(legend.position='none',
          plot.margin = unit(c(0,0,-1,0), "cm"))
  gg_annet
  
  gg_NPP = ggplot() + 
    geom_sf(data = wrld_simpl_sf_eckiv, fill = "grey75",
            colour = "grey45", size = 0.2) +
    theme_map() +
    geom_tile(data = r_eckiv_dt_var[,.(cell, x, y, NPP_2100)],
              aes(x = x, y = y, fill = NPP_2100)) +
    scale_fill_manual(values = colors2) +
    theme(legend.position='none',
          plot.margin = unit(c(0,0,-1,0), "cm"))
  gg_NPP
  
  colors_bar1 = c("#276419", "#4D9221","#7FBC41","#B8E186","#E6F5D0",
                           "#FDE0EF","#F1B6DA","#DE77AE", "#C51B7D", "#8E0152")

  colors_bar2 =c("#8E0152", "#C51B7D","#DE77AE", "#F1B6DA", 
                          "#FDE0EF", "#E6F5D0","#B8E186", "#7FBC41", "#4D9221", "#276419")
                          
  # create legend
  gg_legend1 = plot_discrete_cbar(c(-Inf,-3.0,-2.0,-1.0,-0.5,0,0.5,1.0,2.0,3.0,Inf),
                                  colors = colors_bar2,
                                  legend_title = expression(Mean~Difference~Soil~C:N),
                                  spacing = 'constant',
                                  font_size = 7)
  gg_legend1
  gg_legend2 = plot_discrete_cbar(c(-Inf,-0.10,-0.05,-0.01,-0.005,0,0.005,0.01,0.05,0.10,Inf),
                                  colors = colors_bar1,
                                  legend_title = expression(Mean~Difference~Decomposition~Fraction),
                                  spacing = 'constant',
                                  font_size = 7)
  gg_legend2
  gg_legend3 = plot_discrete_cbar(c(-Inf,-20,-10,-5,-1,0,1,5,10,20,Inf),
                                  colors = colors_bar2,
                                  legend_title = expression(Mean~Difference~Annual~ET~'('*cm~yr^-1*')'),
                                  spacing = 'constant',
                                  font_size = 7)
  gg_legend3
  gg_legend4 = plot_discrete_cbar(c(-Inf,-15000,-10000,-1000,-500,0,500,1000,10000,15000,Inf),
                                  colors = colors_bar2,
                                  legend_title = expression(Mean~Difference~NPP~'('~kg~ha^-1*')'),
                                  spacing = 'constant',
                                  font_size = 7)
  gg_legend4
  
  gg_maps = list(CN = gg_CN, sldcmp = gg_sldcmp, sfdcmp = gg_sfdcmp, annet = gg_annet, NPP = gg_NPP,
                 CN_legend = gg_legend1, dcmp_legend = gg_legend2, annet_legend = gg_legend3, NPP_legend = gg_legend4)
  return(gg_maps)
}
cumul_line               = function(dt, dt_gcm) {
  require(ggplot2)
  require(data.table)
  
  # add 2016, start at 0
  dt_s = data.table(ssp    = c('ssp126', 'ssp370', 'historical'), y_block = 2016,
                    Ts_SOC_m = 0L, Ts_N2O_m = 0L, Ts_GHG_m = 0L)
  cols = c('ssp', 'y_block','Ts_SOC_m', 'Ts_N2O_m', 'Ts_GHG_m')
  dt   = dt[, ..cols]
  dt   = rbind(dt, dt_s)
  setorder(dt, ssp, y_block)
  
  dt_s_gcm = data.table(ssp     = c(rep('ssp126',28), rep('ssp370',24), rep('historical',1)), 
                        gcm     = c(unique(dt_gcm[ssp %in% 'ssp126', gcm]),
                                    unique(dt_gcm[ssp %in% 'ssp370', gcm]),
                                    unique(dt_gcm[ssp %in% 'historical', gcm])),
                        y_block = rep(2016,53),
                        Ts_SOC  = rep(0L,53), Ts_N2O = rep(0L,53), Ts_GHG = rep(0L,53))
  cols   = c('ssp','gcm','y_block','Ts_SOC', 'Ts_N2O', 'Ts_GHG')
  dt_gcm = dt_gcm[, ..cols]
  dt_gcm = rbind(dt_gcm, dt_s_gcm)
  setorder(dt_gcm, ssp, y_block)
  
  ssp_lbl = c('ssp126' = c('SSP1-2.6'), 'ssp370' = c('SSP3-7.0'), 'historical' = c('No Change'))
  
  gg_n2o = ggplot(data = dt) +
    geom_line(aes(x = y_block, y = (Ts_N2O_m/1000L)), color = 'aquamarine4', size = 1.25) + # convert to Pg
    geom_hline(yintercept = 0) +
    facet_grid(.~ssp, labeller = labeller(ssp = ssp_lbl)) +
    xlab('Time') +
    ylab(expression(atop(paste(GHG~Emissions~Cumulative~Difference), '('*Pg~CO[2]*-eq*')'))) +
    scale_y_continuous(limits = c(-25, 25), breaks = seq(-25, 25, 5)) +
    scale_x_continuous(limits = c(2016,2100), breaks = c(2020, 2050, 2080, 2100)) +
    theme_bw() +
    theme(text = element_text(color = 'black', size = 16),
          axis.text.x  = element_text(angle = 45, hjust = 1),
          axis.text    = element_text(size = 14, color = 'black'),
          strip.text   = element_text(size = 14, color = 'black'),
          axis.title.x = element_text(size = 16),
          axis.title.y = element_text(size = 16),
          legend.position = 'none')
  
  gg_n2o_all = gg_n2o +
    geom_line(data = dt_gcm, aes(x = y_block, y = (Ts_N2O/1000L), group = gcm), color = 'aquamarine4', alpha = 0.2) + # convert to Pg
    theme(legend.position = 'none')
  
  gg_SOC = gg_n2o_all +
    geom_line(data = dt, aes(x = y_block, y = (Ts_SOC_m/1000L)), color = 'deeppink4', size = 1.25)
  
  gg_SOC_all = gg_SOC +
    geom_line(data = dt_gcm, aes(x = y_block, y = (Ts_SOC/1000L), group = gcm), color = 'deeppink4', alpha = 0.2) + # convert to Pg
    theme(legend.position = 'none')
  
  gg_GHG = gg_SOC_all +
    geom_line(data = dt, aes(x = y_block, y = (Ts_GHG_m/1000L)), size = 1.25) 
  
  gg_GHG_all = gg_GHG +
    geom_line(data = dt_gcm, aes(x = y_block, y = (Ts_GHG/1000L), group = gcm), color = 'black', alpha = 0.2) + # convert to Pg
    theme(legend.position = 'none')
  
  ggplots = list(GHG = gg_GHG_all)
  return(ggplots)
}
cumul_line_st            = function(dt, dt_gcm) {
  require(ggplot2)
  require(data.table)
  
  # add 2016, start at 0
  dt_s = data.table(ssp    = c('ssp126', 'ssp370', 'historical'), y_block = 2016,
                    Ts_SOC_m = 0L, Ts_N2O_m = 0L, Ts_GHG_m = 0L)
  cols = c('ssp', 'y_block','Ts_SOC_m', 'Ts_N2O_m', 'Ts_GHG_m')
  dt   = dt[, ..cols]
  dt   = rbind(dt, dt_s)
  setorder(dt, ssp, y_block)
  
  dt_s_gcm = data.table(ssp     = c(rep('ssp126',28), rep('ssp370',24), rep('historical',1)), 
                        gcm     = c(unique(dt_gcm[ssp %in% 'ssp126', gcm]),
                                    unique(dt_gcm[ssp %in% 'ssp370', gcm]),
                                    unique(dt_gcm[ssp %in% 'historical', gcm])),
                        y_block = rep(2016,53),
                        Ts_SOC  = rep(0L,53), Ts_N2O = rep(0L,53), Ts_GHG = rep(0L,53))
  cols   = c('ssp','gcm','y_block','Ts_SOC', 'Ts_N2O', 'Ts_GHG')
  dt_gcm = dt_gcm[, ..cols]
  dt_gcm = rbind(dt_gcm, dt_s_gcm)
  setorder(dt_gcm, ssp, y_block)
  
  ssp_lbl = c('ssp126' = c('SSP1-2.6'), 'ssp370' = c('SSP3-7.0'), 'historical' = c('No Change'))
  
  gg_n2o = ggplot(data = dt) +
    geom_line(aes(x = y_block, y = (Ts_N2O_m/1000L)), color = 'aquamarine4', size = 1.25) + # convert to Pg
    geom_hline(yintercept = 0) +
    facet_grid(.~ssp, labeller = labeller(ssp = ssp_lbl)) +
    xlab('Time') +
    ylab(expression(atop(paste(GHG~Emissions~Cumulative~Difference), '('*Pg~CO[2]*-eq*')'))) +
    scale_y_continuous(limits = c(-65, 20), breaks = seq(-65, 20, 10)) +
    scale_x_continuous(limits = c(2016,2100), breaks = c(2020, 2050, 2080, 2100)) +
    theme_bw() +
    theme(text = element_text(color = 'black', size = 16),
          axis.text.x  = element_text(angle = 45, hjust = 1),
          axis.text    = element_text(size = 14, color = 'black'),
          strip.text   = element_text(size = 14, color = 'black'),
          axis.title.x = element_text(size = 16),
          axis.title.y = element_text(size = 16),
          legend.position = 'none')
  
  gg_n2o_all = gg_n2o +
    geom_line(data = dt_gcm, aes(x = y_block, y = (Ts_N2O/1000L), group = gcm), color = 'aquamarine4', alpha = 0.2) + # convert to Pg
    theme(legend.position = 'none')
  
  gg_SOC = gg_n2o_all +
    geom_line(data = dt, aes(x = y_block, y = (Ts_SOC_m/1000L)), color = 'deeppink4', size = 1.25)
  
  gg_SOC_all = gg_SOC +
    geom_line(data = dt_gcm, aes(x = y_block, y = (Ts_SOC/1000L), group = gcm), color = 'deeppink4', alpha = 0.2) + # convert to Pg
    theme(legend.position = 'none')
  
  gg_GHG = gg_SOC_all +
    geom_line(data = dt, aes(x = y_block, y = (Ts_GHG_m/1000L)), size = 1.25) 
  
  gg_GHG_all = gg_GHG +
    geom_line(data = dt_gcm, aes(x = y_block, y = (Ts_GHG/1000L), group = gcm), color = 'black', alpha = 0.2) + # convert to Pg
    theme(legend.position = 'none')
  
  # gg_SOC1 = gg_SOC +
  # geom_line(aes(x = y_block, y = (Ts_SOC_ciH)/1000L)) + 
  # geom_line(aes(x = y_block, y = (Ts_SOC_ciL)/1000L)) +
  # geom_ribbon(aes(x = y_block, ymin=(Ts_SOC_ciL/1000L),ymax=(Ts_SOC_ciH/1000L)), fill="blue", alpha=0.5)
  
  ggplots = list(GHG = gg_GHG_all)
  return(ggplots)
}
cumul_rg_line            = function(dt, dt_gcm) {
  require(ggplot2)
  require(data.table)
  
  # add 2016, start at 0
  dt_s = data.table(ssp    = c('ssp126', 'ssp370', 'historical'), y_block = 2016,
                    Ts_SOC_m = 0L, Ts_N2O_m = 0L, Ts_GHG_m = 0L)
  cols = c('ssp', 'y_block','Ts_SOC_m', 'Ts_N2O_m', 'Ts_GHG_m')
  dt   = dt[, ..cols]
  dt   = rbind(dt, dt_s)
  setorder(dt, ssp, y_block)
  
  dt_s_gcm = data.table(ssp     = c(rep('ssp126',28), rep('ssp370',24), rep('historical',1)), 
                        gcm     = c(unique(dt_gcm[ssp %in% 'ssp126', gcm]),
                                    unique(dt_gcm[ssp %in% 'ssp370', gcm]),
                                    unique(dt_gcm[ssp %in% 'historical', gcm])),
                        y_block = rep(2016,53),
                        Ts_SOC  = rep(0L,53), Ts_N2O = rep(0L,53), Ts_GHG = rep(0L,53))
  cols   = c('ssp','gcm','y_block','Ts_SOC', 'Ts_N2O', 'Ts_GHG')
  dt_gcm = dt_gcm[, ..cols]
  dt_gcm = rbind(dt_gcm, dt_s_gcm)
  setorder(dt_gcm, ssp, y_block)
  
  ssp_lbl = c('ssp126' = c('SSP1-2.6'), 'ssp370' = c('SSP3-7.0'), 'historical' = c('No Change'))
  
  gg_n2o = ggplot(data = dt) +
    geom_line(aes(x = y_block, y = (Ts_N2O_m/1000L)), color = 'aquamarine4', size = 1.25) + # convert to Pg
    geom_hline(yintercept = 0) +
    facet_grid(.~ssp, labeller = labeller(ssp = ssp_lbl)) +
    xlab('Time') +
    ylab(expression(atop(paste(GHG~Emissions~Cumulative~Difference), '('*Pg~CO[2]*-eq*')'))) +
    scale_y_continuous(limits = c(-20, 15), breaks = seq(-20, 15, 5)) +
    scale_x_continuous(limits = c(2016,2100), breaks = c(2020, 2050, 2080, 2100)) +
    theme_bw() +
    theme(text = element_text(color = 'black', size = 16),
          axis.text.x  = element_text(angle = 45, hjust = 1),
          axis.text    = element_text(size = 14, color = 'black'),
          strip.text   = element_text(size = 14, color = 'black'),
          axis.title.x = element_text(size = 16),
          axis.title.y = element_text(size = 16),
          legend.position = 'none')
  
  gg_n2o_all = gg_n2o +
    geom_line(data = dt_gcm, aes(x = y_block, y = (Ts_N2O/1000L), group = gcm), color = 'aquamarine4', alpha = 0.2) + # convert to Pg
    theme(legend.position = 'none')
  
  gg_SOC = gg_n2o_all +
    geom_line(data = dt, aes(x = y_block, y = (Ts_SOC_m/1000L)), color = 'deeppink4', size = 1.25)
  
  gg_SOC_all = gg_SOC +
    geom_line(data = dt_gcm, aes(x = y_block, y = (Ts_SOC/1000L), group = gcm), color = 'deeppink4', alpha = 0.2) + # convert to Pg
    theme(legend.position = 'none')
  
  gg_GHG = gg_SOC_all +
    geom_line(data = dt, aes(x = y_block, y = (Ts_GHG_m/1000L)), size = 1.25) 
  
  gg_GHG_all = gg_GHG +
    geom_line(data = dt_gcm, aes(x = y_block, y = (Ts_GHG/1000L), group = gcm), color = 'black', alpha = 0.2) + # convert to Pg
    theme(legend.position = 'none')
  
  ggplots = list(GHG = gg_GHG_all)
  return(ggplots)
}

annual_total_rg_bp       = function(dt_gcm, .ssp, .region) {

  scenario_lbl = c('res' = c('Residue'), 'ntill' = c('No-tillage'), 'ccg' = c('Grass CC'),
                   'ccl' = c('Legume CC'), 'ccg-ntill' = c('Grass CC Combined Practices'),
                   'ccl-ntill' = c('Legume CC Combined Practices'))
  # wide to long format
  dt = melt(dt_gcm,
       id.vars = c("gcm", "ssp", "y_block", "scenario","IPCC_NAME"),
       measure.vars = colnames(dt_gcm)[5:20])
  dt[, y_block := as.factor(y_block)]
  
  dt$scenario = factor(dt$scenario, levels = c('res', 'ntill', 'ccg', 'ccg-ntill', 'ccl', 'ccl-ntill'))

  bplot_n2o = ggplot() +
    geom_boxplot(data = dt[ssp %in% .ssp & IPCC_NAME %in% .region &
                             y_block %in% c('2030', '2050', '2100') &
                             variable %in% 'Tm_N2O'],
                 aes(x = y_block, y = value, group = y_block, fill = variable), 
                 position=position_dodge(0.1), width = 0.8, fill = "grey60", alpha = 0.9) +
    facet_wrap(.~scenario, nrow = 4, ncol = 2, labeller = labeller(scenario = scenario_lbl)) +
    geom_hline(yintercept = 0) +
    xlab('Time') +
    ylab(expression(atop(paste(Annual~N[2]*O~Emissions~Difference), '('*Tg~CO[2]*-eq*~yr^-1*')'))) +
    scale_y_continuous(limits = c(-200, 200), breaks = seq(-200, 200, 100)) +
    theme_bw() +
    theme(text = element_text(color = 'black', size = 16),
          axis.text.x  = element_text(angle = 45, hjust = 1),
          axis.text    = element_text(size = 14, color = 'black'),
          strip.text   = element_text(size = 14, color = 'black'),
          axis.title.x = element_text(size = 16),
          axis.title.y = element_text(size = 16),
          legend.position = 'none')
  
  bplot_soc = ggplot() +
    geom_boxplot(data = dt[ssp %in% .ssp & IPCC_NAME %in% .region &
                             y_block %in% c('2030', '2050', '2100') &
                             variable %in% 'Tm_SOC'],
                 aes(x = y_block, y = value, group = y_block, fill = variable), 
                 position=position_dodge(0.1), width = 0.8, fill = "grey60", alpha = 0.9) +
    facet_wrap(.~scenario, nrow = 4, ncol = 2, labeller = labeller(scenario = scenario_lbl)) +
    geom_hline(yintercept = 0) +
    xlab('Time') +
    ylab(expression(atop(paste(Annual~SOC~Difference), '('*Tg~CO[2]*-eq*~yr^-1*')'))) +
    scale_y_continuous(limits = c(-500, 100), breaks = seq(-500, 100, 100)) +
    theme_bw() +
    theme(text = element_text(color = 'black', size = 16),
          axis.text.x  = element_text(angle = 45, hjust = 1),
          axis.text    = element_text(size = 14, color = 'black'),
          strip.text   = element_text(size = 14, color = 'black'),
          axis.title.x = element_text(size = 16),
          axis.title.y = element_text(size = 16),
          legend.position = 'none')
  
  bplot_ghg = ggplot() +
    geom_boxplot(data = dt[ssp %in% .ssp & IPCC_NAME %in% .region & 
                             y_block %in% c('2030', '2050', '2100') &
                             variable %in% 'Tm_GHG'],
                 aes(x = y_block, y = value, group = y_block, fill = variable), 
                 position=position_dodge(0.1), width = 0.5, fill = "grey60", alpha = 0.9) +
    facet_wrap(.~scenario, nrow = 4, ncol = 2, labeller = labeller(scenario = scenario_lbl)) +
    geom_hline(yintercept = 0) +
    xlab('Time') +
    ylab(expression(atop(paste(Annual~GHG~Emissions~Difference), '('*Tg~CO[2]*-eq*~yr^-1*')'))) +
    scale_y_continuous(limits = c(-500, 100), breaks = seq(-500, 100, 100)) +
    theme_bw() +
    theme(text = element_text(color = 'black', size = 16),
          axis.text.x  = element_text(angle = 45, hjust = 1),
          axis.text    = element_text(size = 14, color = 'black'),
          strip.text   = element_text(size = 14, color = 'black'),
          axis.title.x = element_text(size = 16),
          axis.title.y = element_text(size = 16),
          legend.position = 'none')
  
  bplot_grain = ggplot() +
    geom_boxplot(data = dt[ssp %in% .ssp & IPCC_NAME %in% .region &
                             y_block %in% c('2030', '2050', '2100') &
                             variable %in% 'Tm_grain'],
                 aes(x = y_block, y = value, group = y_block, fill = variable), 
                 position=position_dodge(0.1), width = 0.8, fill = "grey60", alpha = 0.9) +
    facet_wrap(.~scenario, nrow = 4, ncol = 2, labeller = labeller(scenario = scenario_lbl)) +
    geom_hline(yintercept = 0) +
    xlab('Time') +
    ylab(expression(atop(paste(Annual~Grain~Yield~Difference), '('*Tg~yr^-1*')'))) +
    scale_y_continuous(limits = c(-200, 200), breaks = seq(-200, 200, 100)) +
    theme_bw() +
    theme(text = element_text(color = 'black', size = 16),
          axis.text.x  = element_text(angle = 45, hjust = 1),
          axis.text    = element_text(size = 14, color = 'black'),
          strip.text   = element_text(size = 14, color = 'black'),
          axis.title.x = element_text(size = 16),
          axis.title.y = element_text(size = 16),
          legend.position = 'none')
  
  bplot_resid = ggplot() +
    geom_boxplot(data = dt[ssp %in% .ssp & IPCC_NAME %in% .region &
                             y_block %in% c('2030', '2050', '2100') &
                             variable %in% 'Tm_resid'],
                 aes(x = y_block, y = value, group = y_block, fill = variable), 
                 position=position_dodge(0.1), width = 0.8, fill = "grey60", alpha = 0.9) +
    facet_wrap(.~scenario, nrow = 4, ncol = 2, labeller = labeller(scenario = scenario_lbl)) +
    geom_hline(yintercept = 0) +
    xlab('Time') +
    ylab(expression(atop(paste(Annual~Crop~Residue~Difference), '('*Tg~yr^-1*')'))) +
    scale_y_continuous(limits = c(-200, 200), breaks = seq(-200, 200, 100)) +
    theme_bw() +
    theme(text = element_text(color = 'black', size = 16),
          axis.text.x  = element_text(angle = 45, hjust = 1),
          axis.text    = element_text(size = 14, color = 'black'),
          strip.text   = element_text(size = 14, color = 'black'),
          axis.title.x = element_text(size = 16),
          axis.title.y = element_text(size = 16),
          legend.position = 'none')
  
  bplot_ccbio = ggplot() +
    geom_boxplot(data = dt[ssp %in% .ssp & IPCC_NAME %in% .region &
                             y_block %in% c('2030', '2050', '2100') &
                             scenario %in% c('ccg', 'ccg-ntill', 'ccl', 'ccl-ntill') & 
                             variable %in% 'Tm_cc_shootC'],
                 aes(x = y_block, y = value, group = y_block, fill = variable), 
                 position=position_dodge(0.1), width = 0.8, fill = "grey60", alpha = 0.9) +
    facet_wrap(.~scenario, nrow = 2, ncol = 2, labeller = labeller(scenario = scenario_lbl)) +
    geom_hline(yintercept = 0) +
    xlab('Time') +
    ylab(expression(atop(paste(Annual~Cover~Crop~Biomass), '('*Tg~yr^-1*')'))) +
    scale_y_continuous(limits = c(0, 1250), breaks = seq(0, 1250, 250)) +
    theme_bw() +
    theme(text = element_text(color = 'black', size = 16),
          axis.text.x  = element_text(angle = 45, hjust = 1),
          axis.text    = element_text(size = 14, color = 'black'),
          strip.text   = element_text(size = 14, color = 'black'),
          axis.title.x = element_text(size = 16),
          axis.title.y = element_text(size = 16),
          legend.position = 'none')
  
  ggplots = list(N2O = bplot_n2o, SOC = bplot_soc, GHG = bplot_ghg, 
                 Yield = bplot_grain, Resid = bplot_resid, CC_bio = bplot_ccbio)
  return(ggplots)
}
rg_annual_sgl_mSE = function(dt, .ssp, .scenarios, .regions, .region) {
  require(ggplot2)
  require(data.table)
  
  # add 2016, start at 0
  dt_s = data.table(IPCC_NAME = rep_len(c(.regions), 20),
                    scenario  = rep_len(c(.scenarios),20),
                    ssp       = rep_len(.ssp, 20), 
                    y_block   = rep_len(2016L, 20),
                    m_SOC_m   = rep_len(0L, 20), 
                    m_N2O_m   = rep_len(0L, 20), 
                    m_GHG_m   = rep_len(0L, 20),
                    m_grain_m = rep_len(0L, 20),
                    Tm_SOC_se = rep_len(0L, 20),
                    Tm_N2O_se = rep_len(0L, 20),
                    Tm_GHG_se = rep_len(0L, 20),
                    Tm_grain_se = rep_len(0L, 20))
  setorder(dt_s, IPCC_NAME)
  cols = c('IPCC_NAME','scenario','ssp', 'y_block','m_SOC_m', 'm_N2O_m', 'm_GHG_m', 'm_grain_m', 
           'Tm_SOC_se', 'Tm_N2O_se', 'Tm_GHG_se', 'Tm_grain_se')
  dt   = dt[, ..cols]
  dt   = rbind(dt_s, dt)
  setorder(dt, IPCC_NAME, ssp, y_block, scenario)
  
  scenario_lbl = c('res' = c('Residue'), 'ntill' = c('No-tillage'), 'ccg' = c('Grass CC'),
                   'ccl' = c('Legume CC'), 'ccg-ntill' = c('Grass CC Combined Practices'),
                   'ccl-ntill' = c('Legume CC Combined Practices'))
  ssp_lbl = c('ssp126' = 'SSP1-2.6', 'ssp370' = 'SSP3-7.0', 'historical' = 'No Change')
  
  gg_n2o = ggplot(data = dt[ssp %in% .ssp & IPCC_NAME %in% .region & scenario %in% .scenarios]) +
    geom_point(aes(x = y_block, y = (m_N2O_m)), color = 'aquamarine4', size = 1.25) +
    geom_line(aes(x = y_block, y = (m_N2O_m)), color = 'aquamarine4') +
    geom_errorbar(aes(x = y_block, y = m_N2O_m, ymin=m_N2O_m-Tm_N2O_se, ymax=m_N2O_m+Tm_N2O_se), width=.2,
                  position=position_dodge(0.05), color = 'aquamarine4') +
    facet_grid(.~scenario, labeller = labeller(scenario = scenario_lbl)) +
    geom_hline(yintercept = 0) +
    xlab('Time') +
    ylab(expression(atop(paste(Annual~GHG~Emissions~Difference), '('*Mg~CO[2]*-eq*~ha^-1*~yr^-1*')'))) +
    scale_y_continuous(limits = c(-3, 1), breaks = seq(-3, 1, 0.5)) +
    scale_x_continuous(limits = c(2016,2100), breaks = c(2020, 2050, 2080, 2100)) +
    theme_bw() +
    theme(text = element_text(color = 'black', size = 16),
          axis.text.x  = element_text(angle = 45, hjust = 1),
          axis.text    = element_text(size = 14, color = 'black'),
          strip.text   = element_text(size = 14, color = 'black'),
          axis.title.x = element_text(size = 16),
          axis.title.y = element_text(size = 16),
          legend.position = 'none',
          panel.background = element_rect(fill = "transparent"), 
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          plot.background = element_rect(fill = "transparent", 
                                         color = NA))
  gg_SOC = gg_n2o +
    geom_point(data = dt[ssp %in% .ssp & IPCC_NAME %in% .region & scenario %in% .scenarios], aes(x = y_block, y = (m_SOC_m)), color = 'deeppink4', size = 1.25) +
    geom_line(aes(x = y_block, y = (m_SOC_m)), color = 'deeppink4') +
    geom_errorbar(aes(x = y_block, y = m_SOC_m, ymin=m_SOC_m-Tm_SOC_se, ymax=m_SOC_m+Tm_SOC_se), width=.2,
                  position=position_dodge(0.05), color = 'deeppink4')
  
  gg_GHG = gg_SOC +
    geom_point(data = dt[ssp %in% .ssp & IPCC_NAME %in% .region & scenario %in% .scenarios], aes(x = y_block, y = (m_GHG_m)), size = 1.25) +
    geom_line(aes(x = y_block, y = (m_GHG_m))) +
    geom_errorbar(aes(x = y_block, y = m_GHG_m, ymin=m_GHG_m-Tm_GHG_se, ymax=m_GHG_m+Tm_GHG_se), width=.2,
                  position=position_dodge(0.05))
  
  gg_grain = ggplot(data = dt[ssp %in% .ssp & IPCC_NAME %in% .region & scenario %in% .scenarios]) +
    facet_grid(.~scenario, labeller = labeller(scenario = scenario_lbl)) +
    geom_hline(yintercept = 0) +
    xlab('Time') +
    ylab(expression(atop(paste(Annual~Grain~Yield~Difference), '('*kg~ha^-1*~yr^-1*')'))) +
    scale_y_continuous(limits = c(-500, 1500), breaks = seq(-500, 1500, 250)) +
    scale_x_continuous(limits = c(2016,2100), breaks = c(2020, 2050, 2080, 2100)) +
    theme_bw() +
    theme(text = element_text(color = 'black', size = 16),
          axis.text.x  = element_text(angle = 45, hjust = 1),
          axis.text    = element_text(size = 14, color = 'black'),
          strip.text   = element_text(size = 14, color = 'black'),
          axis.title.x = element_text(size = 16),
          axis.title.y = element_text(size = 16),
          legend.position = 'none',
          panel.background = element_rect(fill = "transparent"), 
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          plot.background = element_rect(fill = "transparent", 
                                         color = NA))
  # gg_grain_se = gg_grain +
  #   geom_line(aes(x = y_block, y = (m_grain_m-(Tm_grain_se*10L))/C_to_mass), color = 'grey40', linetype = 'dashed') +
  #   geom_line(aes(x = y_block, y = (m_grain_m+(Tm_grain_se*10L))/C_to_mass), color = 'grey40', linetype = 'dashed')
  gg_grain_se = gg_grain +
    geom_ribbon(aes(x = y_block, y = m_grain_m, ymin = (m_grain_m-(Tm_grain_se)),
                    ymax = (m_grain_m+(Tm_grain_se))), fill = 'grey70', alpha = 0.6) +
    geom_point(aes(x = y_block, y = (m_grain_m)), size = 1.25) +
    geom_line(aes(x = y_block, y = (m_grain_m)))
    
  ggplots = list(GHG = gg_GHG, Yield = gg_grain_se)
  return(ggplots)
}
rg_annual_cmb_mSE = function(dt, .ssp, .scenarios, .regions, .region) {
  require(ggplot2)
  require(data.table)
  
  # add 2016, start at 0
  dt_s = data.table(IPCC_NAME = rep_len(c(.regions), 10),
                    scenario  = rep_len(c(.scenarios),10),
                    ssp       = rep_len(.ssp, 10), 
                    y_block   = rep_len(2016L, 10),
                    m_SOC_m   = rep_len(0L, 10), 
                    m_N2O_m   = rep_len(0L, 10), 
                    m_GHG_m   = rep_len(0L, 10),
                    m_grain_m = rep_len(0L, 10),
                    Tm_SOC_se = rep_len(0L, 10),
                    Tm_N2O_se = rep_len(0L, 10),
                    Tm_GHG_se = rep_len(0L, 10),
                    Tm_grain_se = rep_len(0L, 10))
  setorder(dt_s, IPCC_NAME)
  cols = c('IPCC_NAME','scenario','ssp', 'y_block','m_SOC_m', 'm_N2O_m', 'm_GHG_m', 'm_grain_m', 
           'Tm_SOC_se', 'Tm_N2O_se', 'Tm_GHG_se', 'Tm_grain_se')
  dt   = dt[, ..cols]
  dt   = rbind(dt_s, dt)
  setorder(dt, IPCC_NAME, ssp, y_block, scenario)
  
  scenario_lbl = c('res' = c('Residue'), 'ntill' = c('No-tillage'), 'ccg' = c('Grass CC'),
                   'ccl' = c('Legume CC'), 'ccg-ntill' = c('Grass CC Combined Practices'),
                   'ccl-ntill' = c('Legume CC Combined Practices'))
  ssp_lbl = c('ssp126' = 'SSP1-2.6', 'ssp370' = 'SSP3-7.0', 'historical' = 'No Change')
  
  gg_n2o = ggplot(data = dt[ssp %in% .ssp & IPCC_NAME %in% .region & scenario %in% .scenarios]) +
    geom_point(aes(x = y_block, y = (m_N2O_m)), color = 'aquamarine4', size = 1.25) +
    geom_line(aes(x = y_block, y = (m_N2O_m)), color = 'aquamarine4') +
    geom_errorbar(aes(x = y_block, y = m_N2O_m, ymin=m_N2O_m-Tm_N2O_se, ymax=m_N2O_m+Tm_N2O_se), width=.2,
                  position=position_dodge(0.05), color = 'aquamarine4') +
    facet_grid(.~scenario, labeller = labeller(scenario = scenario_lbl)) +
    geom_hline(yintercept = 0) +
    xlab('Time') +
    ylab(expression(atop(paste(Annual~Difference~GHG~Emissions), '('*Mg~CO[2]*-eq*~ha^-1*~yr^-1*')'))) +
    scale_y_continuous(limits = c(-5, 1), breaks = seq(-5, 1, 0.5)) +
    scale_x_continuous(limits = c(2016,2100), breaks = c(2020, 2050, 2080, 2100)) +
    theme_bw() +
    theme(text = element_text(color = 'black', size = 16),
          axis.text.x  = element_text(angle = 45, hjust = 1),
          axis.text    = element_text(size = 14, color = 'black'),
          strip.text   = element_text(size = 14, color = 'black'),
          axis.title.x = element_text(size = 16),
          axis.title.y = element_text(size = 16),
          legend.position = 'none',
          panel.background = element_rect(fill = "transparent"), 
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          plot.background = element_rect(fill = "transparent", 
                                         color = NA))
  gg_SOC = gg_n2o +
    geom_point(data = dt[ssp %in% .ssp & IPCC_NAME %in% .region & scenario %in% .scenarios], aes(x = y_block, y = (m_SOC_m)), color = 'deeppink4', size = 1.25) +
    geom_line(aes(x = y_block, y = (m_SOC_m)), color = 'deeppink4') +
    geom_errorbar(aes(x = y_block, y = m_SOC_m, ymin=m_SOC_m-Tm_SOC_se, ymax=m_SOC_m+Tm_SOC_se), width=.2,
                  position=position_dodge(0.05), color = 'deeppink4')
  
  gg_GHG = gg_SOC +
    geom_point(data = dt[ssp %in% .ssp & IPCC_NAME %in% .region & scenario %in% .scenarios], aes(x = y_block, y = (m_GHG_m)), size = 1.25) +
    geom_line(aes(x = y_block, y = (m_GHG_m))) +
    geom_errorbar(aes(x = y_block, y = m_GHG_m, ymin=m_GHG_m-Tm_GHG_se, ymax=m_GHG_m+Tm_GHG_se), width=.2,
                  position=position_dodge(0.05))
  
  gg_grain = ggplot(data = dt[ssp %in% .ssp & IPCC_NAME %in% .region & scenario %in% .scenarios]) +
    facet_grid(.~scenario, labeller = labeller(scenario = scenario_lbl)) +
    geom_hline(yintercept = 0) +
    xlab('Time') +
    ylab(expression(atop(paste(Annual~Grain~Yield~Difference), '('*kg~ha^-1*~yr^-1*')'))) +
    scale_y_continuous(limits = c(-750, 1500), breaks = seq(-750, 1500, 250)) +
    scale_x_continuous(limits = c(2016,2100), breaks = c(2020, 2050, 2080, 2100)) +
    theme_bw() +
    theme(text = element_text(color = 'black', size = 16),
          axis.text.x  = element_text(angle = 45, hjust = 1),
          axis.text    = element_text(size = 14, color = 'black'),
          strip.text   = element_text(size = 14, color = 'black'),
          axis.title.x = element_text(size = 16),
          axis.title.y = element_text(size = 16),
          legend.position = 'none',
          panel.background = element_rect(fill = "transparent"), 
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          plot.background = element_rect(fill = "transparent", 
                                         color = NA))
  # gg_grain_se = gg_grain +
  #   geom_line(aes(x = y_block, y = (m_grain_m-(Tm_grain_se*10L))/C_to_mass), color = 'grey40', linetype = 'dashed') +
  #   geom_line(aes(x = y_block, y = (m_grain_m+(Tm_grain_se*10L))/C_to_mass), color = 'grey40', linetype = 'dashed')
  gg_grain_se = gg_grain +
    geom_ribbon(aes(x = y_block, y = m_grain_m, ymin = (m_grain_m-(Tm_grain_se)),
                    ymax = (m_grain_m+(Tm_grain_se))), fill = 'grey70', alpha = 0.6) +
    geom_point(aes(x = y_block, y = (m_grain_m)), size = 1.25) +
    geom_line(aes(x = y_block, y = (m_grain_m)))
  
  ggplots = list(GHG = gg_GHG, Yield = gg_grain_se)
  return(ggplots)
}

#-----------------------------------------------------------------------------------------
# GRAVEYARD KEEP FOR NOW
#-----------------------------------------------------------------------------------------
# cell_crop_area_weights   = function(.dt, .lu_path, .raster) {
#   # CREATE dt
#   crop_area_r  = rast(paste(.lu_path, .raster, sep = '/'))
#   crop_area_r$cell_area = cellSize(crop_area_r, mask=FALSE, lyrs=FALSE, unit="ha")
#   crop_area_dt = as.data.table(terra::as.data.frame(crop_area_r, xy = TRUE, cells = TRUE))
#   
#   # FILTER raster by relevant gridid
#   gridid_all       = unique(.dt[,gridid])
#   crop_area_dt_f   = crop_area_dt[cell %in% gridid_all,]
#   
#   # ADD relevant crop area by crop / irrigation / gridid
#   wht_rn = .dt[crop %in% 'wht' & irr == 0,]
#   wht_rn = wht_rn[crop_area_dt_f[,c('cell','wheat_rainfed_2015', 'cell_area')], on = .(gridid = cell)]
#   wht_rn = wht_rn[!is.na(x),]
#   names(wht_rn)[names(wht_rn) == 'wheat_rainfed_2015'] = 'crop_area_ha'
#   
#   maiz_rn = .dt[crop %in% 'maiz' & irr == 0,]
#   maiz_rn = maiz_rn[crop_area_dt_f[,c('cell','maize_rainfed_2015', 'cell_area')], on = .(gridid = cell)]
#   maiz_rn = maiz_rn[!is.na(x),]
#   names(maiz_rn)[names(maiz_rn) == 'maize_rainfed_2015'] = 'crop_area_ha'
#   
#   soyb_rn = .dt[crop %in% 'soyb' & irr == 0,]
#   soyb_rn = soyb_rn[crop_area_dt_f[,c('cell','soybean_rainfed_2015', 'cell_area')], on = .(gridid = cell)]
#   soyb_rn = soyb_rn[!is.na(x),]
#   names(soyb_rn)[names(soyb_rn) == 'soybean_rainfed_2015'] = 'crop_area_ha'
#   
#   wht_ir = .dt[crop %in% 'wht' & irr == 1,]
#   wht_ir = wht_ir[crop_area_dt_f[,c('cell','wheat_irrigated_2015', 'cell_area')], on = .(gridid = cell)]
#   wht_ir = wht_ir[!is.na(x),]
#   names(wht_ir)[names(wht_ir) == 'wheat_irrigated_2015'] = 'crop_area_ha'
#   
#   maiz_ir = .dt[crop %in% 'maiz' & irr == 1,]
#   maiz_ir = maiz_ir[crop_area_dt_f[,c('cell','maize_irrigated_2015', 'cell_area')], on = .(gridid = cell)]
#   maiz_ir = maiz_ir[!is.na(x),]
#   names(maiz_ir)[names(maiz_ir) == 'maize_irrigated_2015'] = 'crop_area_ha'
#   
#   soyb_ir = .dt[crop %in% 'soyb' & irr == 1,]
#   soyb_ir = soyb_ir[crop_area_dt_f[,c('cell','soybean_irrigated_2015', 'cell_area')], on = .(gridid = cell)]
#   soyb_ir = soyb_ir[!is.na(x),]
#   names(soyb_ir)[names(soyb_ir) == 'soybean_irrigated_2015'] = 'crop_area_ha'
#   
#   .dt_crop_area = rbind(maiz_rn, soyb_rn, wht_rn, maiz_ir, soyb_ir, wht_ir)
#   
#   # ROUND area
#   .dt_crop_area[, crop_area_ha := round(crop_area_ha, digits = 2)]
#   .dt_crop_area[, cell_total_area_ha := round(cell_area, digits = 2)]
#   setorder(.dt_crop_area, gridid)
#   
#   # CHECK gridid > 0
#   print(length(unique(.dt_crop_area[crop_area_ha > 0, gridid])))
#   
#   # Weighted mean by crop area in gridid, and by irrigated area in gridid
#   .dt_crop_area[, total_crop_area_ha := sum(crop_area_ha), by = .(y_block, ssp, gcm, gridid, scenario)]
#   .dt_crop_area[, weights            := (crop_area_ha/total_crop_area_ha)]
#   .dt_crop_area[, weights            := round(weights, digits = 2)]
# 
#   .dt_crop_area = .dt_crop_area[y_block == 2020,] # most inclusive subset
#   .dt_crop_area = .dt_crop_area[, c('gridid', 'crop', 'irr' , 'weights', 'crop_area_ha')]
#   .dt_crop_area = unique(.dt_crop_area)
#   .dt_crop_area = .dt_crop_area[weights > 0,]
#   print(length(unique(.dt_crop_area[, gridid])))
#   setorder(.dt_crop_area, gridid)
#   
#   .dt_crop_area_c = dcast(.dt_crop_area,
#                gridid ~ crop + irr,
#                value.var = c("weights", "crop_area_ha"))
#   .dt_crop_area_c[is.na(weights_maiz_0), weights_maiz_0 := 0L]
#   .dt_crop_area_c[is.na(weights_maiz_1), weights_maiz_1 := 0L]
#   .dt_crop_area_c[is.na(weights_soyb_0), weights_soyb_0 := 0L]
#   .dt_crop_area_c[is.na(weights_soyb_1), weights_soyb_1 := 0L]
#   .dt_crop_area_c[is.na(weights_wht_0) , weights_wht_0  := 0L]
#   .dt_crop_area_c[is.na(weights_wht_1) , weights_wht_1  := 0L]
# 
#   .dt_crop_area_c[is.na(crop_area_ha_maiz_0), crop_area_ha_maiz_0 := 0L]
#   .dt_crop_area_c[is.na(crop_area_ha_maiz_1), crop_area_ha_maiz_1 := 0L]
#   .dt_crop_area_c[is.na(crop_area_ha_soyb_0), crop_area_ha_soyb_0 := 0L]
#   .dt_crop_area_c[is.na(crop_area_ha_soyb_1), crop_area_ha_soyb_1 := 0L]
#   .dt_crop_area_c[is.na(crop_area_ha_wht_0) , crop_area_ha_wht_0  := 0L]
#   .dt_crop_area_c[is.na(crop_area_ha_wht_1) , crop_area_ha_wht_1  := 0L]
# 
#   .dt_crop_area_c[, irr_frac  := (weights_maiz_1 + weights_soyb_1 + weights_wht_1)]
#   .dt_crop_area_c[, rain_frac := (weights_maiz_0 + weights_soyb_0 + weights_wht_0)]
#   .dt_crop_area_c[, maiz_frac := (weights_maiz_1 + weights_maiz_0)]
#   .dt_crop_area_c[, soyb_frac := (weights_soyb_1 + weights_soyb_0)]
#   .dt_crop_area_c[, wht_frac  := (weights_wht_1 + weights_wht_0)]
#   
#   return(.dt_crop_area_c)
# }
# cell_crop_area_wmean     = function(.dt, .lu_path, .raster) {
#   # CREATE dt
#   crop_area_r           = rast(paste(.lu_path, .raster, sep = '/'))
#   crop_area_r$cell_area = cellSize(crop_area_r, mask=FALSE, lyrs=FALSE, unit="ha")
#   crop_area_dt          = as.data.table(terra::as.data.frame(crop_area_r, xy = TRUE, cells = TRUE))
#   
#   # FILTER raster by relevant gridid
#   gridid_all       = unique(.dt[,gridid])
#   crop_area_dt_f   = crop_area_dt[cell %in% gridid_all,]
#   # ADD relevant crop area by crop / irrigation / gridid
#   wht_rn = .dt[crop %in% 'wht' & irr == 0,]
#   wht_rn = wht_rn[crop_area_dt_f[,c('cell','wheat_rainfed_2015', 'cell_area')], on = .(gridid = cell)]
#   wht_rn = wht_rn[!is.na(x),]
#   names(wht_rn)[names(wht_rn) == 'wheat_rainfed_2015'] = 'crop_area_ha'
#   
#   maiz_rn = .dt[crop %in% 'maiz' & irr == 0,]
#   maiz_rn = maiz_rn[crop_area_dt_f[,c('cell','maize_rainfed_2015', 'cell_area')], on = .(gridid = cell)]
#   maiz_rn = maiz_rn[!is.na(x),]
#   names(maiz_rn)[names(maiz_rn) == 'maize_rainfed_2015'] = 'crop_area_ha'
#   
#   soyb_rn = .dt[crop %in% 'soyb' & irr == 0,]
#   soyb_rn = soyb_rn[crop_area_dt_f[,c('cell','soybean_rainfed_2015', 'cell_area')], on = .(gridid = cell)]
#   soyb_rn = soyb_rn[!is.na(x),]
#   names(soyb_rn)[names(soyb_rn) == 'soybean_rainfed_2015'] = 'crop_area_ha'
#   
#   wht_ir = .dt[crop %in% 'wht' & irr == 1,]
#   wht_ir = wht_ir[crop_area_dt_f[,c('cell','wheat_irrigated_2015', 'cell_area')], on = .(gridid = cell)]
#   wht_ir = wht_ir[!is.na(x),]
#   names(wht_ir)[names(wht_ir) == 'wheat_irrigated_2015'] = 'crop_area_ha'
#   
#   maiz_ir = .dt[crop %in% 'maiz' & irr == 1,]
#   maiz_ir = maiz_ir[crop_area_dt_f[,c('cell','maize_irrigated_2015', 'cell_area')], on = .(gridid = cell)]
#   maiz_ir = maiz_ir[!is.na(x),]
#   names(maiz_ir)[names(maiz_ir) == 'maize_irrigated_2015'] = 'crop_area_ha'
#   
#   soyb_ir = .dt[crop %in% 'soyb' & irr == 1,]
#   soyb_ir = soyb_ir[crop_area_dt_f[,c('cell','soybean_irrigated_2015', 'cell_area')], on = .(gridid = cell)]
#   soyb_ir = soyb_ir[!is.na(x),]
#   names(soyb_ir)[names(soyb_ir) == 'soybean_irrigated_2015'] = 'crop_area_ha'
#   
#   .dt_crop_area = rbind(maiz_rn, soyb_rn, wht_rn, maiz_ir, soyb_ir, wht_ir)
#   
#   # ROUND area
#   .dt_crop_area[, crop_area_ha := round(crop_area_ha, digits = 2)]
#   .dt_crop_area[, cell_area := round(cell_area, digits = 2)]
#   
#   # CHECK gridid > 0
#   print(length(unique(.dt_crop_area[crop_area_ha > 0, gridid])))
#   
#   # Weighted mean by crop area in gridid
#   .dt_crop_area[, total_crop_area_ha := sum(crop_area_ha), by = .(y_block, ssp, gcm, gridid, scenario)]
#   .dt_crop_area[, weights            := (crop_area_ha/total_crop_area_ha)*100L]
#   cols_for_mean = c('m_cc_rootC', 'm_cc_shootC', 'm_cc_shootN', 's_cc_rootC', 's_cc_shootC',  's_cc_shootN',
#                     'm_cr_rootC', 'm_cr_shootC', 'm_cr_shootN', 'm_cr_grain', 'm_cr_grainN','m_cr_NPP', 'm_cr_residC', 
#                     'm_cr_residN', 's_cr_rootC', 's_cr_shootC', 's_cr_shootN', 's_cr_grain', 's_cr_grainN','s_cr_NPP',
#                     's_cr_residC',  's_cr_residN', 'm_SOC', 'm_N2O', 'm_iN2O', 'm_iN2O_v', 'm_iN2O_l','m_dN2O', 'm_dN2O_nit', 'm_dN2O_dnit',
#                     'm_GHG', 'm_N2_to_N2O','m_gr_nit','m_ANNPPT','m_PET_r', 'm_annet','m_sfdcmp', 'm_sldcmp', 'm_fert.N', 'm_omad.C', 'm_omad.N', 'm_cr_irr', 
#                     'm_sC.N','m_nfix', 's_cr_irr', 's_annet','s_SOC', 's_N2O', 's_iN2O', 's_iN2O_v', 's_iN2O_l','s_dN2O', 's_dN2O_nit', 's_dN2O_dnit','s_GHG',
#                     's_N2_to_N2O', 's_gr_nit')
#   print('Compute weighted mean by crop area in gridcell. Note: Calculation takes several minutes.')
#   .dt_crop_area = .dt_crop_area[, lapply(.SD, weighted.mean, weights), .SDcols = cols_for_mean, 
#                                 by = .(gridid, x, y, scenario, y_block, gcm, ssp, WB_NAME, WB_REGION,total_crop_area_ha, cell_area)]
#   .dt_crop_area = .dt_crop_area[, lapply(.SD, round, digits = 2), .SDcols = cols_for_mean,
#                                 by = .(gridid, x, y, scenario, y_block, gcm, ssp, WB_NAME, WB_REGION, total_crop_area_ha, cell_area)]
#   setorder(.dt_crop_area, gridid)
#   .dt_crop_area = .dt_crop_area[total_crop_area_ha > 0,]
#   print(length(unique(.dt_crop_area[, gridid])))
#   .dt_crop_area[, crop_area_weight := (total_crop_area_ha/cell_area)*100L]
#   return(.dt_crop_area)
# }
# cell_crop_area_gcm_wmean = function(.dt, .lu_path, .raster) {
#   # CREATE dt
#   crop_area_r           = rast(paste(.lu_path, .raster, sep = '/'))
#   crop_area_r$cell_area = cellSize(crop_area_r, mask=FALSE, lyrs=FALSE, unit="ha")
#   crop_area_dt          = as.data.table(terra::as.data.frame(crop_area_r, xy = TRUE, cells = TRUE))
#   
#   # FILTER raster by relevant gridid
#   gridid_all       = unique(.dt[,gridid])
#   crop_area_dt_f   = crop_area_dt[cell %in% gridid_all,]
#   
#   # Simplify dt
#   keep_cols     = c('gridid', 'x', 'y', 'scenario', 'crop', 'irr','y_block', 'gcm', 'ssp', 'WB_NAME', 'WB_REGION','m_SOC', 'm_N2O', 'm_iN2O', 'm_dN2O', 'm_cr_grain', 'm_cr_residC',
#                     'm_cc_shootC', 's_cc_shootC', 's_cr_grain', 's_cr_residC','m_GHG','s_SOC', 's_N2O', 's_iN2O', 's_dN2O', 's_GHG')
#   .dt = .dt[, ..keep_cols]
#   
#   # ADD relevant crop area by crop / irrigation / gridid
#   wht_rn = .dt[crop %in% 'wht' & irr == 0,]
#   wht_rn = wht_rn[crop_area_dt_f[,c('cell','wheat_rainfed_2015', 'cell_area')], on = .(gridid = cell)]
#   wht_rn = wht_rn[!is.na(x),]
#   names(wht_rn)[names(wht_rn) == 'wheat_rainfed_2015'] = 'crop_area_ha'
#   
#   maiz_rn = .dt[crop %in% 'maiz' & irr == 0,]
#   maiz_rn = maiz_rn[crop_area_dt_f[,c('cell','maize_rainfed_2015', 'cell_area')], on = .(gridid = cell)]
#   maiz_rn = maiz_rn[!is.na(x),]
#   names(maiz_rn)[names(maiz_rn) == 'maize_rainfed_2015'] = 'crop_area_ha'
#   
#   soyb_rn = .dt[crop %in% 'soyb' & irr == 0,]
#   soyb_rn = soyb_rn[crop_area_dt_f[,c('cell','soybean_rainfed_2015', 'cell_area')], on = .(gridid = cell)]
#   soyb_rn = soyb_rn[!is.na(x),]
#   names(soyb_rn)[names(soyb_rn) == 'soybean_rainfed_2015'] = 'crop_area_ha'
#   
#   wht_ir = .dt[crop %in% 'wht' & irr == 1,]
#   wht_ir = wht_ir[crop_area_dt_f[,c('cell','wheat_irrigated_2015', 'cell_area')], on = .(gridid = cell)]
#   wht_ir = wht_ir[!is.na(x),]
#   names(wht_ir)[names(wht_ir) == 'wheat_irrigated_2015'] = 'crop_area_ha'
#   
#   maiz_ir = .dt[crop %in% 'maiz' & irr == 1,]
#   maiz_ir = maiz_ir[crop_area_dt_f[,c('cell','maize_irrigated_2015', 'cell_area')], on = .(gridid = cell)]
#   maiz_ir = maiz_ir[!is.na(x),]
#   names(maiz_ir)[names(maiz_ir) == 'maize_irrigated_2015'] = 'crop_area_ha'
#   
#   soyb_ir = .dt[crop %in% 'soyb' & irr == 1,]
#   soyb_ir = soyb_ir[crop_area_dt_f[,c('cell','soybean_irrigated_2015', 'cell_area')], on = .(gridid = cell)]
#   soyb_ir = soyb_ir[!is.na(x),]
#   names(soyb_ir)[names(soyb_ir) == 'soybean_irrigated_2015'] = 'crop_area_ha'
#   
#   .dt_crop_area = rbind(maiz_rn, soyb_rn, wht_rn, maiz_ir, soyb_ir, wht_ir)
#   
#   # ROUND area
#   .dt_crop_area[, crop_area_ha := round(crop_area_ha, digits = 2)]
#   .dt_crop_area[, cell_area := round(cell_area, digits = 2)]
#   
#   # CHECK gridid > 0
#   print(length(unique(.dt_crop_area[crop_area_ha > 0, gridid])))
#   
#   # Weighted mean by crop area in gridid
#   .dt_crop_area[, total_crop_area_ha := sum(crop_area_ha), by = .(y_block, ssp, gcm, gridid, scenario)]
#   .dt_crop_area[, weights            := (crop_area_ha/total_crop_area_ha)*100L]
#   cols_for_mean = c('m_SOC', 'm_N2O', 'm_iN2O', 'm_dN2O', 'm_cr_grain', 'm_cr_residC',
#                     'm_cc_shootC', 's_cc_shootC', 's_cr_grain', 's_cr_residC','m_GHG','s_SOC', 's_N2O', 's_iN2O', 's_dN2O', 's_GHG')
#   print('Compute weighted mean by crop area in gridcell. Note: Calculation takes several minutes.')
#   .dt_crop_area = .dt_crop_area[, lapply(.SD, weighted.mean, weights), .SDcols = cols_for_mean, 
#                                 by = .(gridid, x, y, scenario, y_block, gcm, ssp, WB_NAME, WB_REGION,total_crop_area_ha, cell_area)]
#   .dt_crop_area = .dt_crop_area[, lapply(.SD, round, digits = 2), .SDcols = cols_for_mean,
#                                 by = .(gridid, x, y, scenario, y_block, gcm, ssp, WB_NAME, WB_REGION, total_crop_area_ha, cell_area)]
#   setorder(.dt_crop_area, gridid)
#   .dt_crop_area = .dt_crop_area[total_crop_area_ha > 0,]
#   print(length(unique(.dt_crop_area[, gridid])))
#   .dt_crop_area[, crop_area_weight := (total_crop_area_ha/cell_area)*100L]
#   return(.dt_crop_area)
# }
