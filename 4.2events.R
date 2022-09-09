load("2. Baseline Characteristics/final_baseline.RData")
baseline_ENA$Index.Drug <- 0
baseline_SV$Index.Drug <- 1
baseline <- rbind(baseline_SV, baseline_ENA)

load("Data/Combined data/IP_entresto.Rdata")
IP_SV <- Entresto_IP
rm(Entresto_IP)
IP_SV$Discharge.Date..yyyy.mm.dd.. <- as.Date(IP_SV$Discharge.Date..yyyy.mm.dd..)
load("Data/Combined data/IP_enalapril.Rdata")
IP_ENA <- Enalapril_IP
rm(Enalapril_IP)
IP_ENA$Admission.Date..yyyy.mm.dd.. <- as.Date(IP_ENA$Admission.Date..yyyy.mm.dd..)

# Chris amendments
IP_2019_2021 <- readRDS("Data/Chris follow up/combined/2019_2021_ip_combined.rds")
IP_2019_2021$Discharge.Date..yyyy.mm.dd.. <- as.Date(IP_2019_2021$Discharge.Date..yyyy.mm.dd..)
library(plyr)
IP_all <- rbind.fill(IP_SV, IP_ENA, IP_2019_2021)
IP_all <- IP_all[IP_all$Reference.Key. %in% baseline$Reference.Key.,]

# Select only admissions after index date and sort them ascendingly
IP_all <- merge(IP_all, baseline[,c("Reference.Key.","Index.Date")], by="Reference.Key.")
IP_all <- IP_all[IP_all$Admission.Date..yyyy.mm.dd.. > IP_all$Index.Date,]
IP_all <- IP_all[order(IP_all$Reference.Key., IP_all$Admission.Date..yyyy.mm.dd..),]

# Chris amendments
# load("Data/Combined data/Demo_entresto.Rdata")
# load("Data/Combined data/Demo_enalapril.Rdata")
# DEMO <- unique(rbind(entresto_demo, enalapril_demo))
DEMO <- readRDS("Data/Chris follow up/combined/2021_11_demo.rds")
baseline <- baseline[,-(3:5)]
baseline <- merge(baseline, DEMO[, c("Reference.Key.", "Sex.", "Date.of.Birth..yyyy.mm.dd..", "Date.of.Registered.Death.")], by="Reference.Key.")

for(pt in DEMO[duplicated(DEMO$Reference.Key.),]$Reference.Key.) {
  DEMO <- DEMO[!(DEMO$Reference.Key.==pt & is.na(DEMO$Date.of.Registered.Death.)),]
}
rm(list=c("entresto_demo", "enalapril_demo", "pt"))
DEMO <- DEMO[DEMO$Reference.Key. %in% baseline$Reference.Key.,]

IP <- IP_all[!duplicated(IP_all[,c("Reference.Key.")]),]
IP <- cbind.data.frame(Reference.Key.=IP$Reference.Key., Admission.first.any=IP$Admission.Date..yyyy.mm.dd..)
baseline <- merge(baseline, IP, by="Reference.Key.", all.x = T)

ICDCM_CHF <- "398.91|402.01|402.11|402.91|404.01|404.03|404.11|404.13|404.91|404.93|428"
#ICDCM_CHF <- "428.0|428.1|428.9"
IP <- IP_all[grepl(ICDCM_CHF, paste(IP_all$Dx.Px.code.,IP_all$Principal.Diagnosis.Code.,IP_all$Diagnosis..rank.2..,IP_all$Diagnosis..rank.3..,IP_all$Diagnosis..rank.4..,IP_all$Diagnosis..rank.5..,IP_all$Diagnosis..rank.6..,IP_all$Diagnosis..rank.7..,IP_all$Diagnosis..rank.8..,IP_all$Diagnosis..rank.9..,IP_all$Diagnosis..rank.10..,IP_all$Diagnosis..rank.11..,IP_all$Diagnosis..rank.12..,IP_all$Diagnosis..rank.13..,IP_all$Diagnosis..rank.14..,IP_all$Diagnosis..rank.15..)),]
IP <- IP[!duplicated(IP[,c("Reference.Key.")]),]
IP <- cbind.data.frame(Reference.Key.=IP$Reference.Key., Admission.first.HF=IP$Admission.Date..yyyy.mm.dd..)
baseline <- merge(baseline, IP, by="Reference.Key.", all.x = T)

baseline$Death.any <- ifelse(baseline$Date.of.Registered.Death. >= baseline$Index.Date, baseline$Date.of.Registered.Death., NA)
class(baseline$Death.any) <- "Date"
ICDCM_CVDEATH <- "I1[0-5]|I2[0-8]|I[346789][0-9]|I5[0-2]|I0[1-9]"
cvdeath <- DEMO[grepl(ICDCM_CVDEATH, paste(DEMO$Death.Cause..Main.Cause.., DEMO$Death.Cause..Supplementary.Cause..)),]
baseline$Death.cv <- ifelse(baseline$Reference.Key. %in% cvdeath$Reference.Key., baseline$Death.any, NA)
class(baseline$Death.cv) <- "Date"
rm(cvdeath)

load("1. Cohort Identification/final_cohort.RData")
load("Data/Combined data/Rx2019.Rdata")
# Chris amendments
RX_2019_2021 <- readRDS("Data/Chris follow up/combined/2019_2021_rx_combined.rds")
RX_all <- unique(rbind(final_ENA[,1:5], final_SV[,1:5], Rx2019, RX_2019_2021))
rm(list=ls(pattern="^final_|^Rx2019$"))
RX_all <- merge(RX_all, baseline[,c("Reference.Key.","Index.Drug","Index.Date")], by="Reference.Key.")
RX_all <- RX_all[RX_all$Prescription.Start.Date. >= RX_all$Index.Date,]
RX_all <- RX_all[order(RX_all$Reference.Key., RX_all$Prescription.Start.Date., RX_all$Prescription.End.Date.),]

SV_grep <- function(df) { return(grepl("^ENTR", df$Drug.Item.Code.) | grepl("ENTRESTO|SACUBITRIL|LCZ", df$Drug.Name., ignore.case = TRUE)) }
ENA_grep <- function(df) { return(grepl("^ENAL", df$Drug.Item.Code.) | grepl("ENALAPRIL|ANALAPRIL|ENAP|LAPRIL|RENITEC", df$Drug.Name., ignore.case = TRUE)) }
RX_all$Drug <- NA
RX_all[SV_grep(RX_all),]$Drug <- 1
RX_all[ENA_grep(RX_all),]$Drug <- 0
RX_switched <- RX_all[RX_all$Drug != RX_all$Index.Drug,]
RX_switched <- RX_switched[!duplicated(RX_switched[,c("Reference.Key.")]),]
RX_switched <- cbind.data.frame(Reference.Key.=RX_switched$Reference.Key., Date.Switched=RX_switched$Prescription.Start.Date.-1)
baseline <- merge(baseline, RX_switched, by="Reference.Key.", all.x = T)

txdiscon <- function(outcol) {
  # Must have RX_all and baseline in workspace
  RX <- RX_all #[RX_all$Drug == RX_all$Index.Drug,]
  RX <- RX[order(RX$Reference.Key., RX$Prescription.Start.Date., RX$Prescription.End.Date.),]
  
  lastrefkey <- -1
  for(i in 1:nrow(RX)) {
    row <- RX[i,]
    if(row$Reference.Key. != lastrefkey) {
      # new pt
      if(lastrefkey!=-1) {
        baseline[baseline$Reference.Key.==lastrefkey, outcol] <<- discondate
      }
      discondate <- NA
      lastenddate <- NA
      lastrefkey <- row$Reference.Key.
    }
    if(!is.na(lastenddate) && (row$Prescription.Start.Date. - lastenddate > 30) && is.na(discondate)) discondate <- lastenddate
    if(is.na(lastenddate) || row$Prescription.End.Date. > lastenddate) lastenddate <- row$Prescription.End.Date.
  }
  baseline[baseline$Reference.Key.==lastrefkey, outcol] <<- discondate
}
baseline$Date.Discontd <- NA
txdiscon("Date.Discontd")
class(baseline$Date.Discontd) <- "Date"

# Chris amendments
# baseline$Date.StudyEnd <- as.Date("2019-07-31")
baseline$Date.StudyEnd <- as.Date("2021-11-01")


# For readmission analysis -------------------------------------------
IP_all <- rbind(IP_SV, IP_ENA)
IP_all <- IP_all[IP_all$Reference.Key. %in% baseline$Reference.Key.,]
IP_all <- merge(IP_all, baseline[,c("Reference.Key.","Index.Date")], by="Reference.Key.")
IP <- IP_all[grepl(ICDCM_CHF, paste(IP_all$Dx.Px.code.,IP_all$Principal.Diagnosis.Code.,IP_all$Diagnosis..rank.2..,IP_all$Diagnosis..rank.3..,IP_all$Diagnosis..rank.4..,IP_all$Diagnosis..rank.5..,IP_all$Diagnosis..rank.6..,IP_all$Diagnosis..rank.7..,IP_all$Diagnosis..rank.8..,IP_all$Diagnosis..rank.9..,IP_all$Diagnosis..rank.10..,IP_all$Diagnosis..rank.11..,IP_all$Diagnosis..rank.12..,IP_all$Diagnosis..rank.13..,IP_all$Diagnosis..rank.14..,IP_all$Diagnosis..rank.15..)),]
IP <- IP[IP$Index.Date >= IP$Admission.Date..yyyy.mm.dd.. & IP$Index.Date <= IP$Discharge.Date..yyyy.mm.dd..,]
IP <- unique(IP)
IP <- IP[order(IP$Reference.Key., IP$Admission.Date..yyyy.mm.dd.., -xtfrm(IP$Discharge.Date..yyyy.mm.dd..)),]
IP <- IP[!duplicated(IP[,c("Reference.Key.")]),]
IP <- cbind.data.frame(Reference.Key.=IP$Reference.Key., Readmit.Index.Date=IP$Discharge.Date..yyyy.mm.dd..)
baseline <- merge(baseline, IP, by="Reference.Key.", all.x = T)

IP <- merge(IP_all, IP, by="Reference.Key.")
IP <- IP[(IP$Admission.Date..yyyy.mm.dd.. - IP$Readmit.Index.Date) <= 30 & (IP$Admission.Date..yyyy.mm.dd.. - IP$Readmit.Index.Date) > 0,]
IP <- IP[order(IP$Reference.Key., IP$Admission.Date..yyyy.mm.dd.., IP$Discharge.Date..yyyy.mm.dd..),]
baseline$Readmit.any <- baseline$Reference.Key. %in% IP$Reference.Key.
IP_dates <- IP[!duplicated(IP[,c("Reference.Key.")]),]
IP_dates <- cbind.data.frame(Reference.Key.=IP_dates$Reference.Key., Readmit.any.Date=IP_dates$Admission.Date..yyyy.mm.dd..)
baseline <- merge(baseline, IP_dates, by="Reference.Key.", all.x = T)

IP <- IP[grepl(ICDCM_CHF, paste(IP$Dx.Px.code.,IP$Principal.Diagnosis.Code.,IP$Diagnosis..rank.2..,IP$Diagnosis..rank.3..,IP$Diagnosis..rank.4..,IP$Diagnosis..rank.5..,IP$Diagnosis..rank.6..,IP$Diagnosis..rank.7..,IP$Diagnosis..rank.8..,IP$Diagnosis..rank.9..,IP$Diagnosis..rank.10..,IP$Diagnosis..rank.11..,IP$Diagnosis..rank.12..,IP$Diagnosis..rank.13..,IP$Diagnosis..rank.14..,IP$Diagnosis..rank.15..)),]
baseline$Readmit.HF <- baseline$Reference.Key. %in% IP$Reference.Key.
IP_dates <- IP[!duplicated(IP[,c("Reference.Key.")]),]
IP_dates <- cbind.data.frame(Reference.Key.=IP_dates$Reference.Key., Readmit.HF.Date=IP_dates$Admission.Date..yyyy.mm.dd..)
baseline <- merge(baseline, IP_dates, by="Reference.Key.", all.x = T)

# --------------------------------------------------------------------


save(baseline, file="final_events_2021.RData")
