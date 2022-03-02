
set_data_table <- function(outcome){
  
  agedat = outcome$agedat;
  PIdat = as.double(outcome$PIdat);
  asstdat = outcome$asstdat;
  beqdat = outcome$beqdat;
  MStatdat = outcome$MStatdat;
  obsdat = outcome$obsdat;
  mxdat = outcome$mxdat;
  mxobsdat = outcome$mxobsdat;
  hsdat = outcome$hsdat;
  hsobsdat = outcome$hsobsdat;
  totobs = outcome$totobs;
  datawgts = outcome$datawgts;
  agedat96 = outcome$agedat96;
  avgage96 = outcome$avgage96;
  HHIDdat = outcome$HHIDdat;
  
  #--------------------- Set up DATA FRAME ----------------------#
  
  
  # PI groups (PI in quantiles already) / CHANGED BEFORE VARS BY YEAR
  PIq = quantcut(PIdat, q=5, labels = FALSE)
  
  
  # Organizing variables
  ID    = data.frame(HHIDdat)
  cohort= data.frame(agedat96)
  PI    = data.frame(PIdat)
  PIq   = data.frame(PIq)
  alive = data.frame(hsobsdat)
  
  # Renaming
  setnames(ID,     c("HHIDdat"), c("ID"))
  setnames(cohort, c("agedat96"), c("cohort"))
  setnames(PI,     c("PIdat"), c("PI"))
  
  setnames(asstdat, c("V21","V22","V23","V24","V25","V26"),
           c("1996","1998","2000","2002","2004","2006"))
  setnames(mxdat,   c("V35","V36","V37","V38","V39","V40"),
           c("1996","1998","2000","2002","2004","2006"))
  setnames(hsdat,   c("V63","V64","V65","V66","V67","V68"),
           c("1996","1998","2000","2002","2004","2006"))
  setnames(alive,   c("V70","V71","V72","V73","V74","V75"),
           c("1996","1998","2000","2002","2004","2006"))
  setnames(agedat,  c("V77","V78","V79","V80","V81","V82"),
           c("1996","1998","2000","2002","2004","2006"))
  
  setnames(MStatdat,  c("V7", "V8", "V9", "V10", "V11", "V12"),
           c("1996","1998","2000","2002","2004","2006"))
  

  
  
  # Generate Data Table
  DATA <- data.table(ID)
  #DATA[, cohort := cohort]
  DATA[, PI  := PI]
  DATA[, gender:=(MStatdat[,1]==2)*1]
  
  
  
  DATA[,PI2:= 0]
  DATA[,PI3:= 0]
  DATA[,PI4:= 0]
  DATA[,PI5:= 0]
  
  DATA[PI > .2 & PI <= .4, PI2:= 1]
  DATA[PI > .4 & PI <= .6, PI3:= 1]
  DATA[PI > .6 & PI <= .8, PI4:= 1]
  DATA[PI > .8, PI5:= 1]
  
  # Cohort / PI groups
  #DATA[cohort < 77, cohortq:= 1]
  #DATA[cohort > 76 & cohort <= 81, cohortq:= 2]
  #DATA[cohort > 81 & cohort <= 86, cohortq:= 3]
  #DATA[cohort > 86 & cohort <= 91, cohortq:= 4]
  #DATA[cohort > 91, cohortq:= 5]
  
  
  
  DATA[PI < .2, PIq:= 1]
  DATA[PI > .2 & PI <= .4, PIq:= 2]
  DATA[PI > .4 & PI <= .6, PIq:= 3]
  DATA[PI > .6 & PI <= .8, PIq:= 4]
  DATA[PI > .8, PIq:= 5]
  
  # Variables by year
  
  DATA[, age_1996 := log(agedat[,1])]
  DATA[, age_1998 := agedat[,2]]
  DATA[, age_2000 := agedat[,3]]
  DATA[, age_2002 := agedat[,4]]
  DATA[, age_2004 := agedat[,5]]
  DATA[, age_2006 := agedat[,6]]
  
  # Alive pattern from MStatdat
  DATA[, alive_1996 := (MStatdat[,1]>0)*1]
  DATA[, alive_1998 := (MStatdat[,2]>0)*1]
  DATA[, alive_2000 := (MStatdat[,3]>0)*1]
  DATA[, alive_2002 := (MStatdat[,4]>0)*1]
  DATA[, alive_2004 := (MStatdat[,5]>0)*1]
  DATA[, alive_2006 := (MStatdat[,6]>0)*1]
  
  
  DATA[, obs_1996 := obsdat[,1]]
  DATA[, obs_1998 := obsdat[,2]]
  DATA[, obs_2000 := obsdat[,3]]
  DATA[, obs_2002 := obsdat[,4]]
  DATA[, obs_2004 := obsdat[,5]]
  DATA[, obs_2006 := obsdat[,6]]
  
  DATA[, obsal_1996 := obs_1996*alive_1996]
  DATA[, obsal_1998 := obs_1998*alive_1998]
  DATA[, obsal_2000 := obs_2000*alive_2000]
  DATA[, obsal_2002 := obs_2002*alive_2002]
  DATA[, obsal_2004 := obs_2004*alive_2004]
  DATA[, obsal_2006 := obs_2006*alive_2006]
  
  DATA[, health_1996 := (1-hsdat[,1])*obsdat[,1]*alive_1996]
  DATA[, health_1998 := (1-hsdat[,2])*obsdat[,2]*alive_1998]
  DATA[, health_2000 := (1-hsdat[,3])*obsdat[,3]*alive_2000]
  DATA[, health_2002 := (1-hsdat[,4])*obsdat[,4]*alive_2002]
  DATA[, health_2004 := (1-hsdat[,5])*obsdat[,5]*alive_2004]
  DATA[, health_2006 := (1-hsdat[,6])*obsdat[,6]*alive_2006]
  
  
  
  
  num_obs = sum(colSums(DATA[,.(obsal_1996,obsal_1998,obsal_2000,obsal_2002,obsal_2004,obsal_2006)]))
  mean_asset = sum(asstdat)/num_obs
  stdev_asset = sqrt(sum(asstdat^2)/num_obs -mean_asset^2)
  
  DATA[,health_bar:=(health_1996*alive_1996+health_1998*alive_1998+health_2000*alive_2000+health_2002*alive_2002+health_2004*alive_2004+health_2006*alive_2006)/(alive_1996+alive_1998+alive_2000+alive_2002+alive_2004+alive_2006)]
  
  
   DATA[, asset_1996 := (asstdat[,1]-mean_asset)/stdev_asset]
   DATA[, asset_1998 := (asstdat[,2]-mean_asset)/stdev_asset]
   DATA[, asset_2000 := (asstdat[,3]-mean_asset)/stdev_asset]
   DATA[, asset_2002 := (asstdat[,4]-mean_asset)/stdev_asset]
   DATA[, asset_2004 := (asstdat[,5]-mean_asset)/stdev_asset]
   DATA[, asset_2006 := (asstdat[,6]-mean_asset)/stdev_asset]
  
  # DATA[, asset_1996 := asstdat[,1]]
  # DATA[, asset_1998 := asstdat[,2]]
  # DATA[, asset_2000 := asstdat[,3]]
  # DATA[, asset_2002 := asstdat[,4]]
  # DATA[, asset_2004 := asstdat[,5]]
  # DATA[, asset_2006 := asstdat[,6]]
  
  
  DATA[obsal_1996==0, asset_1996 := 0]
  DATA[obsal_1998==0, asset_1998 := 0]
  DATA[obsal_2000==0, asset_2000 := 0]
  DATA[obsal_2002==0, asset_2002 := 0]
  DATA[obsal_2004==0, asset_2004 := 0]
  DATA[obsal_2006==0, asset_2006 := 0]
  
  DATA[,obsal_1996:=NULL]
  DATA[,obsal_1998:=NULL]
  DATA[,obsal_2000:=NULL]
  DATA[,obsal_2002:=NULL]
  DATA[,obsal_2004:=NULL]
  DATA[,obsal_2006:=NULL]
  
  DATA[,obs_1996:=NULL]
  DATA[,obs_1998:=NULL]
  DATA[,obs_2000:=NULL]
  DATA[,obs_2002:=NULL]
  DATA[,obs_2004:=NULL]
  DATA[,obs_2006:=NULL]
  
  
  
  #DATA[, medEx_1996 := mxdat[,1]]
  #DATA[, medEx_1998 := mxdat[,2]]
  #DATA[, medEx_2000 := mxdat[,3]]
  #DATA[, medEx_2002 := mxdat[,4]]
  #DATA[, medEx_2004 := mxdat[,5]]
  #DATA[, medEx_2006 := mxdat[,6]]
  
  # More variables (from dataprep2.dta) 
  # ALIVE, gender, marital status, assets, PI, medEx, health status
  # CHECK IT'S OK !
  #stata_vars  <- read_dta("G:/My Drive/3. NYU/Research/GAN/Nacho/Stata/dataprep2_edited_wide.dta")
  #more_vars   <- data.table(stata_vars)
  #DATA_merged <- merge(DATA, more_vars, by.x= "ID", by.y= "HHID")
  
  # Prepare data to reshape data: LONG -> WIDE
  # colA= paste("age",    seq(1996, 2006, by=2), sep= "_") 
  # colB= paste("alive",  seq(1996, 2006, by=2), sep= "_")
  # colC= paste("asset",  seq(1996, 2006, by=2), sep= "_")
  # colD= paste("health", seq(1996, 2006, by=2), sep= "_")
  # colE= paste("medEx",  seq(1996, 2006, by=2), sep= "_")
  # 
  # DATAlong= melt(DATA, measure    = list(colA, colB, colC, colD, colE), 
  #                value.name = c("age", "alive", "asset", "health", "medEx"))
  # 
  # setnames(DATAlong, c("variable"), c("period"))
  # 
  # # Convert from factor --> numeric
  # DATAlong[, period:= as.numeric(levels(period))[period]] 
  # 
  # # Select relevant data & check groups to be disregarded
  # DATAcut= DATAlong[alive==1]
  # DATAcut[, asset_p50:= median(asset), by=.(period, cohortq, PIq)]
  # DATAcut[, medEx_p50:= median(medEx), by=.(period, cohortq, PIq)]
  # DATAcut[, asset_p50_v2:= median(asset), by=.(period, PIq)] # Same but NO cohort distinction
  # DATAcut[, medEx_p50_v2:= median(medEx), by=.(period, PIq)]
  # DATAcut[, count:= sum(alive), by=.(period, cohortq, PIq)] # Use only groups w/ >10 obs
  # 
  # # Collapse data to see statistics
  # library(doBy)
  # #collapse <- summaryBy(asset + + period ~ PIq + cohortq + period, FUN=c(median, sum), data= DATAcut)
  # collapse    <- summaryBy(asset_p50 + medEx_p50 + count ~ PIq + cohortq + period, data= DATAcut)
  # collapse_v2 <- summaryBy(asset_p50_v2 + medEx_p50_v2 ~ PIq + period, data= DATAcut)
  # 
  # # Additional check
  # library(haven)
  # more_vars  <- read_dta("C:/wealth_singles/more_vars_v3wide.dta")
  # alive_stata <- data.table(more_vars)
  # 
  # DATAv2 <- data.table(ID)
  # 
  # DATAv2[, obs_1996 := obsdat[,1]]
  # DATAv2[, obs_1998 := obsdat[,2]]
  # DATAv2[, obs_2000 := obsdat[,3]]
  # DATAv2[, obs_2002 := obsdat[,4]]
  # DATAv2[, obs_2004 := obsdat[,5]]
  # DATAv2[, obs_2006 := obsdat[,6]]
  # 
  # DATAv2[, health_1996 := hsobsdat[,1]]
  # DATAv2[, health_1998 := hsobsdat[,2]]
  # DATAv2[, health_2000 := hsobsdat[,3]]
  # DATAv2[, health_2002 := hsobsdat[,4]]
  # DATAv2[, health_2004 := hsobsdat[,5]]
  # DATAv2[, health_2006 := hsobsdat[,6]]
  # 
  # DATAv2[, medEx_1996 := mxobsdat[,1]]
  # DATAv2[, medEx_1998 := mxobsdat[,2]]
  # DATAv2[, medEx_2000 := mxobsdat[,3]]
  # DATAv2[, medEx_2002 := mxobsdat[,4]]
  # DATAv2[, medEx_2004 := mxobsdat[,5]]
  # DATAv2[, medEx_2006 := mxobsdat[,6]]
  # 
  # CHECK <- merge(DATAv2, alive_stata, by.x= "ID", by.y= "HHID")
  # 
  # # More checks
  # #nacho <- as.data.table(DATAlong)
  # #nacho <- nacho[medEx > 0]
  # #nacho <- nacho[, period:=  as.numeric(levels(period))[period]]
  # #nacho[, asset.p50:= median(asset), by=.(period, cohortq, PIq)]
  # #nacho[, count:= sum(alive), by=.(period, cohortq, PIq)]
  # 
  # #nacho_collapse <- summaryBy(asset.p50 + count ~ PIq + cohortq + period, data= nacho)
  # 
  list_ret = list(DATA,mean_asset,stdev_asset)
  return(list_ret)
}


set_data_table2 <- function(outcome){
  
  agedat = outcome$agedat;
  PIdat = as.double(outcome$PIdat);
  asstdat = outcome$asstdat;
  beqdat = outcome$beqdat;
  MStatdat = outcome$MStatdat;
  obsdat = outcome$obsdat;
  mxdat = outcome$mxdat;
  mxobsdat = outcome$mxobsdat;
  hsdat = outcome$hsdat;
  hsobsdat = outcome$hsobsdat;
  totobs = outcome$totobs;
  datawgts = outcome$datawgts;
  agedat96 = outcome$agedat96;
  avgage96 = outcome$avgage96;
  HHIDdat = outcome$HHIDdat;
  
  #--------------------- Set up DATA FRAME ----------------------#
  
  
  # PI groups (PI in quantiles already) / CHANGED BEFORE VARS BY YEAR
  PIq = quantcut(PIdat, q=5, labels = FALSE)
  
  
  # Organizing variables
  ID    = data.frame(HHIDdat)
  cohort= data.frame(agedat96)
  PI    = data.frame(PIdat)
  PIq   = data.frame(PIq)
  alive = data.frame(hsobsdat)
  
  # Renaming
  setnames(ID,     c("HHIDdat"), c("ID"))
  setnames(cohort, c("agedat96"), c("cohort"))
  setnames(PI,     c("PIdat"), c("PI"))
  
  setnames(asstdat, c("V21","V22","V23","V24","V25","V26"),
           c("1996","1998","2000","2002","2004","2006"))
  setnames(mxdat,   c("V35","V36","V37","V38","V39","V40"),
           c("1996","1998","2000","2002","2004","2006"))
  setnames(hsdat,   c("V63","V64","V65","V66","V67","V68"),
           c("1996","1998","2000","2002","2004","2006"))
  setnames(alive,   c("V70","V71","V72","V73","V74","V75"),
           c("1996","1998","2000","2002","2004","2006"))
  setnames(agedat,  c("V77","V78","V79","V80","V81","V82"),
           c("1996","1998","2000","2002","2004","2006"))
  
  # Generate Data Table
  DATA <- data.table(ID)
  DATA[, cohort := cohort]
  DATA[, PI  := PI]
  
  # Cohort / PI groups
  DATA[cohort < 77, cohortq:= 1]
  DATA[cohort > 76 & cohort <= 81, cohortq:= 2]
  DATA[cohort > 81 & cohort <= 86, cohortq:= 3]
  DATA[cohort > 86 & cohort <= 91, cohortq:= 4]
  DATA[cohort > 91, cohortq:= 5]
  
  DATA[PI < .2, PIq:= 1]
  DATA[PI > .2 & PI <= .4, PIq:= 2]
  DATA[PI > .4 & PI <= .6, PIq:= 3]
  DATA[PI > .6 & PI <= .8, PIq:= 4]
  DATA[PI > .8, PIq:= 5]
  
  # Variables by year
  DATA[, age_1996 := agedat[,1]]
  DATA[, age_1998 := agedat[,2]]
  DATA[, age_2000 := agedat[,3]]
  DATA[, age_2002 := agedat[,4]]
  DATA[, age_2004 := agedat[,5]]
  DATA[, age_2006 := agedat[,6]]
  
  DATA[, alive_1996 := alive[,1]]
  DATA[, alive_1998 := alive[,2]]
  DATA[, alive_2000 := alive[,3]]
  DATA[, alive_2002 := alive[,4]]
  DATA[, alive_2004 := alive[,5]]
  DATA[, alive_2006 := alive[,6]]
  
  DATA[, asset_1996 := asstdat[,1]]
  DATA[, asset_1998 := asstdat[,2]]
  DATA[, asset_2000 := asstdat[,3]]
  DATA[, asset_2002 := asstdat[,4]]
  DATA[, asset_2004 := asstdat[,5]]
  DATA[, asset_2006 := asstdat[,6]]
  
  
  DATA[, health_1996 := (1-hsdat[,1])*obsdat[,1]*alive_1996]
  DATA[, health_1998 := (1-hsdat[,2])*obsdat[,2]*alive_1998]
  DATA[, health_2000 := (1-hsdat[,3])*obsdat[,3]*alive_2000]
  DATA[, health_2002 := (1-hsdat[,4])*obsdat[,4]*alive_2002]
  DATA[, health_2004 := (1-hsdat[,5])*obsdat[,5]*alive_2004]
  DATA[, health_2006 := (1-hsdat[,6])*obsdat[,6]*alive_2006]
  
  
  DATA[, medEx_1996 := mxdat[,1]]
  DATA[, medEx_1998 := mxdat[,2]]
  DATA[, medEx_2000 := mxdat[,3]]
  DATA[, medEx_2002 := mxdat[,4]]
  DATA[, medEx_2004 := mxdat[,5]]
  DATA[, medEx_2006 := mxdat[,6]]
  
  # More variables (from dataprep2.dta) 
  # ALIVE, gender, marital status, assets, PI, medEx, health status
  # CHECK IT'S OK !
  #stata_vars  <- read_dta("G:/My Drive/3. NYU/Research/GAN/Nacho/Stata/dataprep2_edited_wide.dta")
  #more_vars   <- data.table(stata_vars)
  #DATA_merged <- merge(DATA, more_vars, by.x= "ID", by.y= "HHID")
  
  # Prepare data to reshape data: LONG -> WIDE
  colA= paste("age",    seq(1996, 2006, by=2), sep= "_") 
  colB= paste("alive",  seq(1996, 2006, by=2), sep= "_")
  colC= paste("asset",  seq(1996, 2006, by=2), sep= "_")
  colD= paste("health", seq(1996, 2006, by=2), sep= "_")
  colE= paste("medEx",  seq(1996, 2006, by=2), sep= "_")
  
  DATAlong= melt(DATA, measure    = list(colA, colB, colC, colD, colE), 
                 value.name = c("age", "alive", "asset", "health", "medEx"))
  
  setnames(DATAlong, c("variable"), c("period"))
  
  # Convert from factor --> numeric
  DATAlong[, period:= as.numeric(levels(period))[period]] 
  
  # Select relevant data & check groups to be disregarded
  DATAcut= DATAlong[alive==1]
  DATAcut[, asset_p50:= median(asset), by=.(period, cohortq, PIq)]
  DATAcut[, medEx_p50:= median(medEx), by=.(period, cohortq, PIq)]
  DATAcut[, asset_p50_v2:= median(asset), by=.(period, PIq)] # Same but NO cohort distinction
  DATAcut[, medEx_p50_v2:= median(medEx), by=.(period, PIq)]
  DATAcut[, count:= sum(alive), by=.(period, cohortq, PIq)] # Use only groups w/ >10 obs
  
  # Collapse data to see statistics
  library(doBy)
  #collapse <- summaryBy(asset + + period ~ PIq + cohortq + period, FUN=c(median, sum), data= DATAcut)
  collapse    <- summaryBy(asset_p50 + medEx_p50 + count ~ PIq + cohortq + period, data= DATAcut)
  collapse_v2 <- summaryBy(asset_p50_v2 + medEx_p50_v2 ~ PIq + period, data= DATAcut)
  
  # Additional check
  library(haven)
  more_vars  <- read_dta("C:/wealth_singles/more_vars_v3wide.dta")
  alive_stata <- data.table(more_vars)
  
  DATAv2 <- data.table(ID)
  
  DATAv2[, obs_1996 := obsdat[,1]]
  DATAv2[, obs_1998 := obsdat[,2]]
  DATAv2[, obs_2000 := obsdat[,3]]
  DATAv2[, obs_2002 := obsdat[,4]]
  DATAv2[, obs_2004 := obsdat[,5]]
  DATAv2[, obs_2006 := obsdat[,6]]
  
  DATAv2[, health_1996 := hsobsdat[,1]]
  DATAv2[, health_1998 := hsobsdat[,2]]
  DATAv2[, health_2000 := hsobsdat[,3]]
  DATAv2[, health_2002 := hsobsdat[,4]]
  DATAv2[, health_2004 := hsobsdat[,5]]
  DATAv2[, health_2006 := hsobsdat[,6]]
  
  DATAv2[, medEx_1996 := mxobsdat[,1]]
  DATAv2[, medEx_1998 := mxobsdat[,2]]
  DATAv2[, medEx_2000 := mxobsdat[,3]]
  DATAv2[, medEx_2002 := mxobsdat[,4]]
  DATAv2[, medEx_2004 := mxobsdat[,5]]
  DATAv2[, medEx_2006 := mxobsdat[,6]]
  
  CHECK <- merge(DATAv2, alive_stata, by.x= "ID", by.y= "HHID")
  
  # More checks
  #nacho <- as.data.table(DATAlong)
  #nacho <- nacho[medEx > 0]
  #nacho <- nacho[, period:=  as.numeric(levels(period))[period]]
  #nacho[, asset.p50:= median(asset), by=.(period, cohortq, PIq)]
  #nacho[, count:= sum(alive), by=.(period, cohortq, PIq)]
  
  #nacho_collapse <- summaryBy(asset.p50 + count ~ PIq + cohortq + period, data= nacho)
  
  return(DATA)
}