load("1. Cohort Identification/final_cohort.RData")

load("Data/Combined data/Ddx_enalapril.Rdata")
load("Data/Combined data/Ddx_entresto.Rdata")
DX_all <- unique(rbind(entresto_ddx, enalapril_ddx))

load("Data/Combined data/IP_entresto.Rdata")
IP_SV <- Entresto_IP
rm(Entresto_IP)
IP_SV$Discharge.Date..yyyy.mm.dd.. <- as.Date(IP_SV$Discharge.Date..yyyy.mm.dd..)
load("Data/Combined data/IP_enalapril.Rdata")
IP_ENA <- Enalapril_IP
rm(Enalapril_IP)
IP_ENA$Admission.Date..yyyy.mm.dd.. <- as.Date(IP_ENA$Admission.Date..yyyy.mm.dd..)
IP_all <- rbind(IP_SV, IP_ENA)

index.drug = "ENA" #### SET THIS

####### To be repeated for each index.drug #######################################

pts <- get(paste0("final_pts_", index.drug))
rx <- get(paste0("final_", index.drug))

pts$Age <- as.numeric((pts$Index.Date - pts$Date.of.Birth..yyyy.mm.dd..) / 365.25)
median(pts$Age)
pts$Sex. <- relevel(as.factor(pts$Sex.), ref="M")
summary(pts$Sex.)



DX <- DX_all[DX_all$Reference.Key. %in% pts$Reference.Key.,]
DX <- merge(DX, pts[,c("Reference.Key.","Index.Date")], by="Reference.Key.")
DX <- DX[DX$Reference.Date. <= DX$Index.Date,]

com.icm = "414.8"
com.htn = "401|402|403|404|405|437.2"
com.dm = "250"
com.af = "427.3"
com.ihd = "410|411|412|413|414"
com.mi = "410|412"
com.ist = "433.01|433.11|433.21|433.31|433.81|433.91|434|436|437.0|437.1"

com_greps <- c("com.icm", "com.htn", "com.dm", "com.af", "com.ihd", "com.mi", "com.ist")
for(comtyp in com_greps) {
  dx <- DX[grepl(get(comtyp), DX$All.Diagnosis.Code..ICD9..),]
  pts[,comtyp] = as.numeric(pts$Reference.Key. %in% dx$Reference.Key.)
}
rm(list=ls(pattern="^com|^dx$"))

##### CCI #####

cci.mi <- "410|412"
cci.chf <- "398.91|402.01|402.11|402.91|404.01|404.03|404.11|404.13|404.91|404.93|428"
cci.pvd <- "441|443.9|785.4|V43.4"
cci.cbd <- "43[0-8]"
cci.copd <- "49[0-6]|50[0-5]|506.4"
cci.dementia <- "290"
cci.paralysis <- "342|344.1"
cci.dm <- "250.[01237]"
cci.dm_ <- "250.[456]"
cci.crf <- "582|583.0|583.1|583.2|583.4|583.6|583.7|585|586|588"
cci.liver_mild <- "571.[2456]"
cci.liver_modsevere <- "456.0|456.1|456.2|572.2|572.3|572.4|572.8"
cci.ulcer <- "53[1-4]"
cci.ra <- "710.0|710.1|710.4|714.0|714.1|714.2|714.81|725"
cci.aids <- "042"
cci.cancer <- "1[458][0-9]|16[0-5]|17[0124569]|1[69][0-5]|20[0-8]"
cci.cancer_mets <- "19[6-9]"

cci_greps <- c("cci.mi", "cci.chf", "cci.pvd", "cci.cbd", "cci.copd", "cci.dementia", "cci.paralysis", "cci.dm", "cci.dm_", "cci.crf", "cci.liver_mild", "cci.liver_modsevere", "cci.ulcer", "cci.ra", "cci.aids", "cci.cancer", "cci.cancer_mets")
cci <- pts
for(ccityp in cci_greps) {
  dx <- DX[grepl(get(ccityp), DX$All.Diagnosis.Code..ICD9..),]
  cci[,ccityp] = as.numeric(cci$Reference.Key. %in% dx$Reference.Key.)
}
cci$cci <- cci$cci.mi + cci$cci.chf + cci$cci.pvd + cci$cci.cbd + cci$cci.copd + cci$cci.dementia + cci$cci.paralysis + pmin(cci$cci.dm + cci$cci.dm_*2, 2) + cci$cci.crf*2 + pmin(cci$cci.liver_mild + cci$cci.liver_modsevere*3, 3) + cci$cci.ulcer + cci$cci.ra + cci$cci.aids*6 + cci$cci.cancer*2 + cci$cci.cancer_mets*6
pts <- merge(pts, cci[,c("Reference.Key.","cci")], by="Reference.Key.")
#assign(paste0("debug_cci_", index.drug), cci)
rm(list=ls(pattern="^cci|^dx$"))

##### END CCI #####

IP <- IP_all[IP_all$Reference.Key. %in% pts$Reference.Key.,]
IP <- merge(IP, pts[,c("Reference.Key.","Index.Date")], by="Reference.Key.")
IP <- IP[order(IP$Reference.Key., IP$Admission.Date..yyyy.mm.dd.., -xtfrm(IP$Discharge.Date..yyyy.mm.dd..)),]
IP <- IP[!duplicated(IP[,c("Reference.Key.","Admission.Date..yyyy.mm.dd..")]),]
IP <- IP[((IP$Admission.Date..yyyy.mm.dd.. <= IP$Index.Date) & (IP$Admission.Date..yyyy.mm.dd.. >= IP$Index.Date - 365)) | ((IP$Discharge.Date..yyyy.mm.dd.. <= IP$Index.Date) & (IP$Discharge.Date..yyyy.mm.dd.. >= IP$Index.Date - 365)),]
IP <- IP[!is.na(IP$Reference.Key.),]

numadmissions <- function(outcol) {
  # Must have IP and pts in workspace
  IP <<- IP[order(IP$Reference.Key., IP$Admission.Date..yyyy.mm.dd.., -xtfrm(IP$Discharge.Date..yyyy.mm.dd..)),]
  lastrefkey <- -1
  for(i in 1:nrow(IP)) {
    row <- IP[i,]
    if(row$Reference.Key. != lastrefkey) {
      # new pt
      if(lastrefkey!=-1) pts[pts$Reference.Key.==lastrefkey, outcol] <<- admissions
      admissions <- 0
      lastdcdate <- NA
      lastrefkey <- row$Reference.Key.
    }
    if(is.na(lastdcdate) || row$Admission.Date..yyyy.mm.dd.. - lastdcdate > 1) admissions <- admissions + 1
    lastdcdate <- row$Discharge.Date..yyyy.mm.dd..
  }
  pts[pts$Reference.Key.==lastrefkey, outcol] <<- admissions
}

pts$Admissions.last.year <- 0
numadmissions("Admissions.last.year")


IP <- IP_all[IP_all$Reference.Key. %in% pts$Reference.Key.,]
IP <- merge(IP, pts[,c("Reference.Key.","Index.Date")], by="Reference.Key.")
IP <- IP[order(IP$Reference.Key., IP$Admission.Date..yyyy.mm.dd.., -xtfrm(IP$Discharge.Date..yyyy.mm.dd..)),]
IP <- IP[!duplicated(IP[,c("Reference.Key.","Admission.Date..yyyy.mm.dd..")]),]
IP <- IP[IP$Admission.Date..yyyy.mm.dd.. <= IP$Index.Date,]
ICDCM_CHF <- "398.91|402.01|402.11|402.91|404.01|404.03|404.11|404.13|404.91|404.93|428"
IP <- IP[grepl(ICDCM_CHF, paste(IP$Dx.Px.code.,IP$Principal.Diagnosis.Code.,IP$Diagnosis..rank.2..,IP$Diagnosis..rank.3..,IP$Diagnosis..rank.4..,IP$Diagnosis..rank.5..,IP$Diagnosis..rank.6..,IP$Diagnosis..rank.7..,IP$Diagnosis..rank.8..,IP$Diagnosis..rank.9..,IP$Diagnosis..rank.10..,IP$Diagnosis..rank.11..,IP$Diagnosis..rank.12..,IP$Diagnosis..rank.13..,IP$Diagnosis..rank.14..,IP$Diagnosis..rank.15..)),]
pts$Admissions.HF <- 0
numadmissions("Admissions.HF")

##### Export results #####

assign(paste0("baseline_", index.drug), pts)
rm(list=ls(pattern="^index[.]drug$|^pts$|^rx$|^DX$|^dx$|^IP$"))

#######################################################################

save_SV <- cbind.data.frame(Reference.Key.=baseline_SV$Reference.Key., CCI=baseline_SV$cci, Adm.last.yr=baseline_SV$Admissions.last.year, Adm.HF=baseline_SV$Admissions.HF)
save_ENA <- cbind.data.frame(Reference.Key.=baseline_ENA$Reference.Key., CCI=baseline_ENA$cci, Adm.last.yr=baseline_ENA$Admissions.last.year, Adm.HF=baseline_ENA$Admissions.HF)
save(save_SV, save_ENA, file="Baseline_cci_adm.RData")
save(baseline_SV, baseline_ENA, file="final_baseline.RData")
