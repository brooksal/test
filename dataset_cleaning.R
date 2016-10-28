######Cleaning ALL Lab data and putting together databases#####
library(lubridate)
library(ggplot2)
library(grid)
library(scales)
library(plyr)
library(dplyr)
library(tidyr)
library(foreign)
library(gplots)
library(reshape2)
library(plotly)




setwd("/Users/alexbrooks/Documents/MP/Data/Rdata/Sep_run/")
#datasets
#nitrate
nitrate <- read.csv("/Users/alexbrooks/Documents/MP/Data/Rdata/24Jun16/postlab/labid_no3n.csv")
#samples list
samples <- read.csv("/Users/alexbrooks/Documents/MP/Data/Rdata/24Jun16/postlab/samples-6.24.csv")
samples[,3] <- mdy_hm(samples[,3])
samples[,3]<- force_tz(samples[,3], "EST5EDT")
samples[,3] <- with_tz(samples[,3], "EST")
colnames(samples)[3] <-"datetime"
samples$datetime <- round_minute(samples$datetime, 10)
#run1 nitrate
samples_sel<- select(samples, c(site, lab_id, datetime)) #simply dataset
samples_sel <- filter(samples_sel, lab_id != "NA")  #remove samples not run for no3.n
all_nitrate <- merge(samples_sel, nitrate, by="lab_id") #simplified table
#run1 nh4 & po4
lat_samp1 <- read.csv("/Users/alexbrooks/Documents/MP/Data/Rdata/24Jun16/postlab/7-20-16_latchatp1.csv")
lat_samp1 <- lat_samp1[1:44,]
lat_samp2<- read.csv("/Users/alexbrooks/Documents/MP/Data/Rdata/24Jun16/postlab/7-20-16_latchatp2.csv")
lat_samp<- rbind.fill(lat_samp1,lat_samp2)
lat_samp <- mutate(lat_samp, Peak.Concentration = ifelse(Peak.Area<.102, -.01, Peak.Concentration))
lat_samp <- mutate(lat_samp, Peak.Concentration.1 = ifelse(Peak.Area.1<.121, -.01, Peak.Concentration.1))
lat_samp<- select(lat_samp, c(Cup.Number, Analyte.Name, Peak.Concentration, Analyte.Name.1, Peak.Concentration.1))
colnames(lat_samp)<- c("lab_id"," ","Ammonia"," ","PO4")
lat_samp<- lat_samp[,c(1,3,5)]
all_nh4_po4<- merge(samples_sel, lat_samp, by="lab_id")
#### Process TN data
tn_df <- read.csv("/Users/alexbrooks/Documents/MP/Data/Rdata/24Jun16/postlab/8-4-16_tn.csv")
tn_df[,3] <- as.numeric(as.character(tn_df$tn))
tndigest <- tn_df[25:116,] #digest samples
#blank and  correction
tndig_blank <-filter(tndigest, grepl("b", lab_id) )  #filtering for blanks
tndig_blank <- tndig_blank[2:4,]#ignoring first b1 since it is needed to be rerun
(b_mean <- mean(tndig_blank$tn))
(b_sd <- sd(tndig_blank$tn))
tndig_bcomp<- mutate(tndigest, tn_b = tn - b_mean ) #subtract mean blank from all values
tndigest_fin <- mutate(tndig_bcomp, tn_fin = tn_b * 1.2) %>%
select(lab_id, tn_fin)#account for dilution of samples by DI
#TP values
phos <- read.csv("/Users/alexbrooks/Documents/MP/Data/Rdata/24Jun16/postlab/8-5-16_TP_adj.csv")
phos_sel<- select( phos, Sample.ID, tp_adj)
phos_sel
phos_blank <-filter(phos_sel, grepl("B", Sample.ID) )
phos_blank <- phos_blank[2:4,] #removing instrument blank
b_mean <- mean(phos_blank$tp_adj)
phos_cor <- mutate(phos_sel, tp_cor = ((tp_adj - b_mean)*1.2))
phos_cor_arrange<- arrange(phos_cor, Sample.ID)
phos_cor_sel <- phos_cor_arrange[9:97,]
phos_cor_sel<- filter(phos_cor_sel, (Sample.ID != "2 ppm") & (Sample.ID != "2ppm"))
colnames(phos_cor_sel)[1] <- "lab_id"
tp<- select(phos_cor_sel, lab_id, tp_adj)
#tdn/NPOC - 1run
tdn_npoc_run1 <- read.csv("/Users/alexbrooks/Documents/MP/Data/Rdata/24Jun16/postlab/tdn_NPOC_run1.csv")
tdn_npoc_run1 <- mutate(tdn_npoc_run1, tdn.run1.mgl  = ifelse(tn.run1.mgl < 0.05, -.05 , tn.run1.mgl))%>%
mutate(npoc.run1.mgl = ifelse(npoc.run1.mgl > 50, -50, ifelse(npoc.run1.mgl <.25, -.25, npoc.run1.mgl)))%>%
select(lab_id = sample, npoc.run1.mgl, tdn.run1.mgl)
#1st run no3, tdn, tn, nh4, po4
nh4_sel<- select(all_nh4_po4, lab_id, Ammonia, PO4)
data1 <- merge(select(all_nitrate,lab_id, no3.n), nh4_sel,  by="lab_id", all.x=TRUE, all.y=TRUE)
data2 <- merge(data1, tndigest_fin, by="lab_id", all.x=TRUE, all.y=TRUE)
run1 <- merge(data2, tdn_npoc_run1, by="lab_id", all.x=TRUE, all.y=TRUE)
#run 2
##TDN
sep10raw <- read.csv("10Sep16_Shimadzu.csv")
sep09raw <- read.csv("09Sep16_Shimadzu.csv")
sep10raw.sel <- select(sep10raw, Sample.Name, Inj..No., Analysis.Inj.., Mean.Conc.) %>%
filter(Analysis.Inj.. == "TN")
sep10_TDN_uni<- distinct(sep10raw.sel,Mean.Conc., .keep_all = TRUE)
sep09raw.sel <- select(sep09raw, Sample.Name, Inj..No., Analysis.Inj.., Mean.Conc.) %>%
filter(Analysis.Inj.. == "TN")
sep09_TDN_uni <- distinct(sep09raw.sel,Mean.Conc., .keep_all = TRUE)
#combine
sep_TN <- rbind.fill(sep09_TDN_uni, sep10_TDN_uni)
sep_TN$Sample.Name<- as.character(sep_TN$Sample.Name)
sep_TDN_sel<- separate(sep_TN, Sample.Name, c("lab_id", "dilution"), sep="_", extra = "merge", fill = "right") %>%
separate(dilution, c("dilution", "extra"), sep="x", extra = "merge", fill = "right") %>%
mutate(tdn = Mean.Conc. * as.numeric(dilution))%>%
filter(extra != "NA")%>%
select(lab_id, tdn_run2 = tdn)
sep_TDN_sel #see clean_sep
sep10rawNPOC.sel <- select(sep10raw, Sample.Name, Inj..No., Analysis.Inj.., Mean.Conc.) %>%
filter(Analysis.Inj.. == "NPOC")
sep10_NPOC_uni<- distinct(sep10rawNPOC.sel,Mean.Conc., .keep_all = TRUE)
sep09rawNPOC.sel <- select(sep09raw, Sample.Name, Inj..No., Analysis.Inj.., Mean.Conc.) %>%
filter(Analysis.Inj.. == "NPOC")
sep09_NPOC_uni <- distinct(sep09rawNPOC.sel,Mean.Conc., .keep_all = TRUE)
#combine
sep_NPOC <- rbind.fill(sep09_NPOC_uni, sep10_NPOC_uni)
sep_NPOC$Sample.Name<- as.character(sep_NPOC$Sample.Name)
sep_NPOC_sel<- separate(sep_NPOC, Sample.Name, c("lab_id", "dilution"), sep="_", extra = "merge", fill = "right") %>%
separate(dilution, c("dilution", "extra"), sep="x", extra = "merge", fill = "right") %>%
mutate(npoc_2ndrun = ifelse(Mean.Conc. > 10, -10 * as.numeric(dilution) , ifelse(Mean.Conc.
<0.25, -.25 * as.numeric(dilution),  Mean.Conc. * as.numeric(dilution))))%>%
filter(extra != "NA")%>%
select(lab_id, npoc_2ndrun)
sep_NPOC_sel
#run 2 latchat
sep_lat_raw <- read.csv("9Sep16_latchat.csv")
colnames(sep_lat_raw)[1] <- "lab_id"
sep_lat_samp<- arrange(sep_lat_raw, lab_id)
sep_lat_samp <- sep_lat_samp[-1:-7,]
sep_lat_samp <- sep_lat_samp[-42:-47,]
sep_lat_samp$NO3.N..mg.L. <- as.numeric(as.character(sep_lat_samp$NO3.N..mg.L.))
sep_lat_samp[20,3] <- sep_lat_samp[21,3]*10
sep_lat_samp<- sep_lat_samp[-21,]
colnames(sep_lat_samp)<- c("lab_id", "nh4_2ndrun", "no3_latchat_2ndrun")
sep_lat_samp$lab_id <- as.numeric(as.character(sep_lat_samp$lab_id))
sep_lat_samp$nh4_2ndrun<- as.numeric(as.character(sep_lat_samp$nh4_2ndrun))
sep_lat_samp <- mutate(sep_lat_samp, nh4_2ndrun = ifelse(is.na(nh4_2ndrun), -.01, nh4_2ndrun)) %>%
mutate(no3_latchat_2ndrun = ifelse(no3_latchat_2ndrun <0.05, -.05, no3_latchat_2ndrun)) %>%
mutate(no3_latchat_2ndrun = ifelse(is.na(no3_latchat_2ndrun), -.05, no3_latchat_2ndrun))
#changing latchat values to reflect linear model at low end of standard curve
#-.1 is below dections for no3, -.01 below detection for latchat
lat_no3_updated<- c(.143,.2004,.101, .157, 089, .087, .198, .108)
lab_id<- c(223,227,232,234,235,244,250,261)
lat_updf<- data.frame(lab_id, lat_no3_updated)
sep_lat_samp <- merge(sep_lat_samp, lat_updf, by="lab_id", all.x=TRUE)
run2_latchat_no3nh4<- mutate(sep_lat_samp, lat_no3_updated = ifelse(is.na(lat_no3_updated), no3_latchat_2ndrun, lat_no3_updated)) %>%
select(lab_id, nh4_2ndrun, lat_no3_2ndrun = lat_no3_updated)
#ic - no3, so4, cl, run
IC_raw <- read.csv("2016_9-9_IC.csv")
IC_raw <- IC_raw[-97:-103,]
IC_sel_no3 <- select(IC_raw, lab_id = Sample.Name, no3_ic_run2 = NO3.N_adj)
so4_jul <- read.csv("7jul16_so4.csv")
so4_jul$lab_id<- as.character(so4_jul$lab_id)
so4_sep <- read.csv("9sep16_so4.csv")
so4 <- merge(so4_sep, so4_jul, by = "lab_id", all.x=TRUE, all.y=TRUE)
#delete? so4_withsites<- merge(so4, IC_sel_no3, by="lab_id")
cl_jul <- read.csv("jul_CL.csv")
cl_sep <- read.csv("sep_CL.csv")
cl <- merge(cl_jul, cl_sep, by = "lab_id", all.x=TRUE, all.y=TRUE)
colnames(so4)<- c("lab_id", "so4_run2", "so4_run1")
colnames(cl) <- c("lab_id", "cl_run1", "cl_run2")
##cations
cations_run1 <- read.csv("2016_7-20-MTM_cations.csv")
cations_run1 <- select(cations_run1, lab_id:ca.run1.mgl)
cations_run2 <- read.csv("2016_9-9-MTM_cations.csv")
cations <- merge(cations_run1, cations_run2, by="lab_id", all.x=TRUE, all.y=TRUE)

#run3
run3Raw<- read.csv("10-5_run3_no3_up.csv")
arrange(run3Raw, lab.id)
run3_no3_oldsamples <- filter(run3Raw, lab.id < 262) %>%
        select(lab.id, no3.n.ic.r3)
run3_no3_newsamples <- filter(run3Raw, lab.id > 261)
run3_no3_newsamples$site <- as.character(run3_no3_newsamples$site)
run3_no3_newsamples$site[run3_no3_newsamples$site == "LF "] <- "LF"
run3_no3_newsamples$site[run3_no3_newsamples$site == "Remington"] <- "REM1"
run3_no3_newsamples$site[run3_no3_newsamples$site == "HC"] <- "HC2"
run3_no3_newsamples$site[run3_no3_newsamples$site == "Sally Fork 2"] <- "SF2"
run3_no3_newsamples$site[run3_no3_newsamples$site == "Sally Fork "] <- "SF1"
run3_no3_newsamples$site[run3_no3_newsamples$site == "Ballard"] <- "BF"
run3_no3_newsamples$site[run3_no3_newsamples$site == "LB"] <- "LB2"
run3_no3_newsamples$site[run3_no3_newsamples$site == "Sugartree"] <- "Sugar Tree"

run3_no3_newsamples$datetime <-  mdy_hm(as.character(run3_no3_newsamples$datetime))
run3_no3_newsamples$datetime <- force_tz(run3_no3_newsamples$datetime, "EST5EDT")
run3_no3_newsamples$datetime<- with_tz(run3_no3_newsamples$datetime, "EST")
run3_no3_newsamples_sel <- select(run3_no3_newsamples, lab.id, site, datetime, no3.n.ic.r3 )

#merge data
data4 <- merge(sep_TDN_sel,run2_latchat_no3nh4, by="lab_id", all.x=TRUE, all.y=TRUE )
data5 <- merge(data4, IC_sel_no3, by="lab_id", all.x=TRUE, all.y=TRUE )
data6 <- merge(data5, so4, by="lab_id", all.x=TRUE, all.y=TRUE )
data7 <- merge(data6, cl, by="lab_id", all.x=TRUE, all.y=TRUE )
data8 <- merge(data7, tp, by="lab_id", all.x=TRUE, all.y=TRUE )
data9 <-  merge(data8, sep_NPOC_sel, by="lab_id", all.x=TRUE, all.y=TRUE )
run2<- merge(data9, cations, by="lab_id", all.x=TRUE, all.y=TRUE)
run2$lab_id <- as.factor(data8$lab_id)
#run2<- select(data7, lab_id:no3_latchat_2ndrun, datetime:cl_run2)
#run2 <- merge(data7, run2_sites, by="lab_id", all.x=TRUE, all.y=TRUE )
#combine (still not all vf data)
fulldata1 <- merge(run1, run2, by=c("lab_id"), all.x=TRUE, all.y=TRUE )
fulldata2 <- merge(fulldata1,run3_no3_oldsamples, by.x="lab_id", by.y="lab.id", all.x=TRUE)
#add site names
run1_site <- select(all_nitrate, lab_id, site)
run2_site <- select(IC_raw, lab_id = Sample.Name, site = Site)
sitelist1 <- merge(run1_site, run2_site, by="lab_id", all.x=TRUE, all.y=TRUE)
sitelist1<- mutate(sitelist1, site = ifelse(is.na(site.x), as.character(site.y), as.character(site.x))) %>%
select(lab_id, site)
sitelist1$site <- as.character(sitelist1$site)
fulldata_withsites<- (merge(sitelist1, fulldata2, by ="lab_id"))
fulldata_withsites <- select(fulldata_withsites, lab.id = lab_id, site, no3.n.ic.r1.mgl =no3.n,
nh4.n.r1.mgl = Ammonia, tn.r1.mgl = tn_fin, tdn.r1.mgl = tdn.run1.mgl, tdn.r2.mgl = tdn_run2,
nh4.r2.mgl = nh4_2ndrun, no3.n.lat.r2.mgl = lat_no3_2ndrun, no3.n.ic.r2.mgl = no3_ic_run2,
so4.r2.mgl = so4_run2, so4.r1.mgl = so4_run1, cl.r1.mgl = cl_run1, cl.r2.mgl = cl_run2,
na.run1.mgl:ca.run2.mgl, po4 = PO4, tp.mgl = tp_adj, npoc.r1.mgl = npoc.run1.mgl, npoc.r2.mgl = npoc_2ndrun, no3.n.ic.r3 )
fulldata_withsites$site[fulldata_withsites$site == "Ballard Fork"] <- "BF"
fulldata_withsites$site[fulldata_withsites$site == "Ballard Top"] <- "BF"
fulldata_withsites$site[fulldata_withsites$site == "HC Top"] <- "HC2"
fulldata_withsites$site[fulldata_withsites$site == "HC"] <- "HC2"
fulldata_withsites$site[fulldata_withsites$site == "LB"] <- "LB2"
fulldata_withsites$site[fulldata_withsites$site == "LB2 Top"] <- "LB2"
fulldata_withsites$site[fulldata_withsites$site == "LB2 Top"] <- "LB2"
fulldata_withsites$site[fulldata_withsites$site == "MR14 "] <- "MR14"
fulldata_withsites$site[fulldata_withsites$site == "Muellers Branch"] <- "MB"
fulldata_withsites$site[fulldata_withsites$site == "Muellers Top"] <- "MB"
fulldata_withsites$site[fulldata_withsites$site == "Muellers Top"] <- "MB"
fulldata_withsites$site[fulldata_withsites$site == "Sally Fork 1"] <- "SF1"
fulldata_withsites$site[fulldata_withsites$site == "Sally Fork 2"] <- "SF2"
fulldata_withsites$site[fulldata_withsites$site == "Spring"] <-"Spring Branch"
fulldata_withsites$site[fulldata_withsites$site == " Stanley"] <-"Stanley Fork"
fulldata_withsites$site[fulldata_withsites$site == "Stanley"] <-"Stanley Fork"
fulldata_withsites$site[fulldata_withsites$site == "Stanley "] <-"Stanley Fork"
fulldata_withsites$site[fulldata_withsites$site == "Sugar"] <-"Sugar Tree"
arrange(fulldata_withsites, lab.id)
filter(fulldata_withsites, duplicated(fulldata_withsites$lab.id) == TRUE)
####datetime -- still need to make sure dates before March 13th are not DLS
run2_datetime<- select(IC_raw, lab.id = Sample.Name, datetime)
run1_datetime<- select(all_nitrate, lab.id= lab_id, datetime )
run2_datetime$datetime <- mdy_hm(as.character(run2_datetime$datetime))
run2_datetime$datetime <- force_tz(run2_datetime$datetime , "EST5EDT")
run2_datetime$datetime <- with_tz(run2_datetime$datetime , "EST")
datetime <- rbind(run1_datetime, run2_datetime)
datetime<- arrange(unique(datetime), lab.id)
fdata1 <- (merge(fulldata_withsites, datetime,by ="lab.id"))
str(fdata1)
#add new run3 sep baseflow data!
fdata <- rbind.fill(fdata1,run3_no3_newsamples_sel  )
fdata <- mutate(fdata, no3.n.mgl= ifelse(is.na(no3.n.ic.r2.mgl), ifelse(is.na(no3.n.ic.r3),no3.n.ic.r1.mgl,no3.n.ic.r3 ), no3.n.ic.r2.mgl))
fdata <- mutate(fdata, note= ifelse(no3.n.mgl == no3.n.ic.r1.mgl, "original ic value used" , ""))
fdata <- mutate(fdata, nh4.n.mgl = ifelse(is.na(nh4.r2.mgl), nh4.n.r1.mgl, nh4.r2.mgl))
fdata <- mutate(fdata, tdn.mgl = ifelse(is.na(tdn.r2.mgl), tdn.r1.mgl, tdn.r2.mgl))
fdata_nitrogen_sel <- select(fdata, lab.id, site,datetime, no3.n.mgl, nh4.n.mgl, tdn.mgl, tn.mgl = tn.r1.mgl,note )
############### need to include other ions in fdata ###########


write.csv(fdata, file="fulldata.csv")
save(fdata_nitrogen_sel, file= "alex_nitrogen.RData")
#melt
fdata_long <- melt(fdata, id=c("lab.id", "site", "datetime"))



#### Bring in ross, lindberg, and pond
setwd("/Users/alexbrooks/Documents/MP/Data/Rdata/analysis/")


##Pond Data
load("alex_nitrogen.RData")
pond_raw <- read.csv("pond_data.csv")
pond_raw_sel <- select(pond_raw, stationID = StationID, site= StreamName, lab.id = SAMPLE.., date = SAMPLE.DATE, no3.n.mgl = VALUE, BelowDet, REPORTING.LEVEL)
## NEED TO IDENTIFY IF SITES ARE CLOSE TO CURRENT SITES ###

#pick sites
studysites <- c("MT01", "MT13", "MT14", "MT15", "MT18", "MT23", "MT24", "MT60")
pond_raw_sel_fil <- filter(pond_raw_sel, stationID %in% studysites )
#remove duplicate 
pond_raw_sel_fil_nodup <- filter(pond_raw_sel_fil, !(is.na(no3.n.mgl) & REPORTING.LEVEL == 0.05))

#correct site naems (note this is PROVISIONAL until GIS analysis is done) -- 
#two stanfley fork sites, don't know which is closer to our current site (or lindberg site)
pond_raw_sel_fil_nodup$stationID<- as.character(pond_raw_sel_fil_nodup$stationID)
pond_raw_sel_fil_nodup$site<- as.character(pond_raw_sel_fil_nodup$site)
pond_raw_sel_fil_nodup$site[pond_raw_sel_fil_nodup$stationID == "MT01"] <- "MR2" ##probably wrong
pond_raw_sel_fil_nodup$site[pond_raw_sel_fil_nodup$stationID == "MT13"] <- "Spring Branch"
pond_raw_sel_fil_nodup$site[pond_raw_sel_fil_nodup$stationID == "MT14"] <- "BF" ##this one is definitely wrong!
pond_raw_sel_fil_nodup$site[pond_raw_sel_fil_nodup$stationID == "MT18"] <- "Sugar Tree" 
pond_raw_sel_fil_nodup$site[pond_raw_sel_fil_nodup$stationID == "MT60"] <- "LF" 
pond_raw_sel_fil_nodup$site[pond_raw_sel_fil_nodup$stationID == "MT23"] <- "MR14" ##probably wrong
pond_raw_sel_fil_nodup$date<- mdy(pond_raw_sel_fil_nodup$date)
#new lab.id
id <- paste("P",1:42,sep="")
pond_raw_sel_fil_nodup <- mutate(pond_raw_sel_fil_nodup, lab.id = id)
pond_df <- pond_raw_sel_fil_nodup

#LINDBERG DATA
lind_raw<- read.csv("Lindberg2011_data.csv")
#select nitrogen daata
lind_nit_sel<- select(lind_raw, site= Sample.Name,Month,Year, no3.n.mgl = NO3.N.mg.l, tdn.mgl = TN.mg.l, Conductivity = Conductivity.uS )
#select
sites_lind <- c("MR2", "MR7", "BF1", "LB1", "MB1", "SF1","ST1")
lind_nit_sel_fil<- filter(lind_nit_sel, site %in% sites_lind)
lind_nit_sel_fil$date<- myd(with(lind_nit_sel_fil, paste(Month,"-",Year, "-1",sep=""))) #made up day for datetime
lind_nitrogen <- lind_nit_sel_fil
#fix site names
lind_nitrogen$site <- as.character(lind_nitrogen$site)
lind_nitrogen$site[lind_nitrogen$site == "LB1"] <- "LB2" #probably wrong location
lind_nitrogen$site[lind_nitrogen$site == "MB1"] <- "MB" #listed as mullin branch not muellers branch but likely same site (google maps seems to confirm)
lind_nitrogen$site[lind_nitrogen$site == "BF1"] <- "BF"  #mainsteam site not at same location as new BF site
lind_nitrogen$site[lind_nitrogen$site == "SF1"] <- "Stanley Fork" #need to confim location 
lind_nitrogen$site[lind_nitrogen$site == "ST1"] <- "Sugar Tree"#need to confim location 


fdata_nitrogen_sel$date <- date(fdata_nitrogen_sel$datetime)

####IC data from October 14 - January 2016 -- adj csv to fill in sites/datetimes that were impromperly stored

allIC_raw <- read.csv("IC.May.2016_abadj.csv")
allIC_nitrate <- select(allIC_raw, lab.id = Sample, site = Site, datetime = DateTime, no3.n.mgl= NO3.N )
#convert datattime
allIC_nitrate$datetime <- mdy_hm(allIC_nitrate$datetime)
allIC_nitrate <- mutate(allIC_nitrate, date = date(allIC_nitrate$datetime))
#convert to chars
allIC_nitrate$site <- as.character(allIC_nitrate$site)
allIC_nitrate$lab.id <- as.character(allIC_nitrate$lab.id)
allIC_nitrate <- mutate(allIC_nitrate, site = ifelse(site == "", strsplit(lab.id, "_")[[1]],site))
#fix location names -- notes I 
allIC_nitrate$site[allIC_nitrate$site == "LF_2"] <- "LF"
allIC_nitrate$site[allIC_nitrate$site == "RB_2"] <- "RB2"
allIC_nitrate$site[allIC_nitrate$site == "RB_1"] <- "RB1"
allIC_nitrate$site[allIC_nitrate$site == "REM_2"] <- "REM1"
allIC_nitrate$site[allIC_nitrate$site == "LB_1"] <- "LB1"
allIC_nitrate$site[allIC_nitrate$site == "LB_2"] <- "LB2"
allIC_nitrate$site[allIC_nitrate$site == "MR_7"] <- "MR7"
allIC_nitrate$site[allIC_nitrate$site == "MR_14"] <- "MR14"
allIC_nitrate$site[allIC_nitrate$site == "MR_2"] <- "MR2"
allIC_nitrate$site[allIC_nitrate$site == "LB_-T"] <- "LB_T"
allIC_nitrate$site[allIC_nitrate$site == "Rem_2"] <- "REM2"
allIC_nitrate$site[allIC_nitrate$site == "HC"] <- "HC2"
allIC_nitrate$site[allIC_nitrate$site == "SpringBranch"] <- "Spring Branch"
allIC_nitrate$site[allIC_nitrate$site == "SallyFork_2"] <- "SF2"
allIC_nitrate$site[allIC_nitrate$site == "SallyFork1"] <- "SF1"
allIC_nitrate$site[allIC_nitrate$site == "StanleyFork_1"] <- "Stanley Fork"
allIC_nitrate$site[allIC_nitrate$site == "SugarTree"] <- "Sugar Tree"
allIC_nitrate$site[allIC_nitrate$site == "MuellersBranch"] <- "MB"
allIC_nitrate$site[allIC_nitrate$site == "BallardFork"] <- "BF"



#merge datasets
comp1<- rbind.fill(lind_nitrogen, fdata_nitrogen_sel)
comp2<- rbind.fill(comp1, pond_df)
comp_nitrogen<- rbind.fill(comp2, allIC_nitrate)
choosesites = c("LB2", "Sugar Tree", "BF", "Stanley Fork", "MB", "LF", "Spring Branch", "MR2", "MR14", "MR7")
comp_nitrate<- select(comp_nitrogen, site, date, no3.n.mgl, note)%>%
        mutate(year = year(date)) %>%
        mutate(year = ifelse(year == 2001, 2000 ,year)) %>%
        filter(site %in% choosesites)
comp_nitrate_mm<- melt(select(comp_nitrate, site,year, no3.n.mgl), id=c("year", "site"))

dt <- data.frame(comp_nitrate_mm)
grp <- group_by(dt, site,year)
sum <- summarise(grp, mean=mean(value, na.rm=TRUE), sd=sd(value, na.rm=TRUE), n=length(value))
print(sum, n=80)

gcomp <- ggplot(comp_nitrate_mm, aes(year, value, shape=site, col = site))+ geom_point()
gcomp <- gcomp + geom_line(data= sum, aes(year,mean, col = site))
ggplotly(gcomp)






load("solinst_2016_sc.Rdata")
load("hobo_2016_data_sc.Rdata")

unique(fdata$site)









##################fooling around
#start <- mdy_hm("6/23/16 0:00", tz="EST")
#end <- mdy_hm("6/25/16 14:00", tz ="EST")
#rough converstion to storm/baseflow
#filter(fdata, datetime >= start & datetime <= end)





#plotting
#+ scale_x_log10() + scale_y_log10()
gbox <- ggplot(fdata, aes(no3.n.ic.r2.mgl, no3.n.ic.r3, color=site))+ geom_point() +geom_abline()
ggplotly(gbox)
gtime <- ggplot(fdata2, aes(x=no3.n.ic.r1.mgl,y=, color = site)) + geom_point()
ggplotly(gtime)
str(fdata2)

dev.off()
gtime <- ggplot(fulldata2, aes(x=datetime,y=npoc_2ndrun, color = site)) + geom_point() + geom_abline()


arrange(select(fdata, site, no3.n.ic.r2.mgl, so4.r2.mgl, cl.r1.mgl,  mg.run1.mgl, ca.run1.mgl    ), no3.n.ic.r2.mgl )
