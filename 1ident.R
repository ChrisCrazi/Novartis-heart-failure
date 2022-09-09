npts <- function (df) { return(length(unique(df$Reference.Key.))) }
printcount <- function(df) { print(sprintf("Pts: %i, Obs: %i", npts(df), nrow(df))) }
SV_grep <- function(df) { return(grepl("^ENTR", df$Drug.Item.Code.) | grepl("ENTRESTO|SACUBITRIL|LCZ", df$Drug.Name., ignore.case = TRUE)) }
ENA_grep <- function(df) { return(grepl("^ENAL", df$Drug.Item.Code.) | grepl("ENALAPRIL|ANALAPRIL|ENAP|LAPRIL|RENITEC", df$Drug.Name., ignore.case = TRUE)) }

load("Data/Combined data/Acohort_entresto.Rdata")
load("Data/Combined data/Rx_entresto.Rdata")
SV <- unique(rbind(Entresto[, 1:5], Entresto_Rx[SV_grep(Entresto_Rx),]))
SV <- SV[!grepl("ENTRECTINIB", SV$Drug.Name.),]

load("Data/Combined data/Acohort_enalapril.Rdata")
load("Data/Combined data/Rx_enalapril.Rdata")
#ENA <- unique(rbind(Enalapril[, 1:5], Enalapril_rx[ENA_grep(Enalapril_rx),]))
#ENA <- ENA[!grepl("ZENAPAX|ACETOMENAPHTHONE|CO-RENITEC", ENA$Drug.Name.),]
ENA <- unique(Enalapril[, 1:5])

cohort <- unique(rbind(SV, ENA))
cohort_ <- cohort[(cohort$Prescription.Start.Date.>=as.Date("2016-07-01")) & (cohort$Prescription.Start.Date.<=as.Date("2019-06-30")), ]
cohort_ <- cohort_[complete.cases(cohort_),]

load("Data/Combined data/Demo_entresto.Rdata")
load("Data/Combined data/Demo_enalapril.Rdata")
DEMO <- unique(rbind(entresto_demo, enalapril_demo))
for(pt in DEMO[duplicated(DEMO$Reference.Key.),]$Reference.Key.) {
  DEMO <- DEMO[!(DEMO$Reference.Key.==pt & is.na(DEMO$Date.of.Registered.Death.)),]
}
cohort_clean <- merge(cohort_, DEMO[, c("Reference.Key.", "Sex.", "Date.of.Birth..yyyy.mm.dd..", "Date.of.Registered.Death.")], by="Reference.Key.")
cohort_clean <- cohort_clean[!(cohort_clean$Prescription.Start.Date. > cohort_clean$Prescription.End.Date.),]
cohort_clean <- cohort_clean[!(cohort_clean$Prescription.Start.Date. < cohort_clean$Date.of.Birth..yyyy.mm.dd..),]
cohort_clean <- cohort_clean[!(!is.na(cohort_clean$Date.of.Registered.Death.) & (cohort_clean$Prescription.Start.Date. > cohort_clean$Date.of.Registered.Death.)),]
cohort_clean <- cohort_clean[!is.na(cohort_clean$Prescription.Start.Date.) & !is.na(cohort_clean$Prescription.End.Date.),]

cohort_firstRx <- cohort_clean[order(cohort_clean$Prescription.Start.Date.),]
cohort_firstRx <- cohort_firstRx[match(unique(cohort_firstRx$Reference.Key.), cohort_firstRx$Reference.Key.),]
pts_SV <- cohort_firstRx[SV_grep(cohort_firstRx),]
pts_ENA <- cohort_firstRx[ENA_grep(cohort_firstRx),]
cohort_SV <- cohort_clean[cohort_clean$Reference.Key. %in% pts_SV$Reference.Key.,]
cohort_ENA <- cohort_clean[cohort_clean$Reference.Key. %in% pts_ENA$Reference.Key.,]

# dup2 <- merge(cohort_firstRx, cohort_clean, by=c("Reference.Key.", "Prescription.Start.Date."))
# dup2 <- dup2[dup2$Drug.Name..x != dup2$Drug.Name..y,]
# => num_obs=7, irrelevant

ICDCM_CHF <- "398.91|402.01|402.11|402.91|404.01|404.03|404.11|404.13|404.91|404.93|428"
load("Data/Combined data/Ddx_enalapril.Rdata")
DX <- enalapril_ddx
DX <- DX[grepl(ICDCM_CHF, DX$All.Diagnosis.Code..ICD9..),]
load("Data/Combined data/IP_enalapril.Rdata")
IP <- Enalapril_IP
IP <- IP[grepl(ICDCM_CHF, paste(IP$Dx.Px.code.,IP$Principal.Diagnosis.Code.,IP$Diagnosis..rank.2..,IP$Diagnosis..rank.3..,IP$Diagnosis..rank.4..,IP$Diagnosis..rank.5..,IP$Diagnosis..rank.6..,IP$Diagnosis..rank.7..,IP$Diagnosis..rank.8..,IP$Diagnosis..rank.9..,IP$Diagnosis..rank.10..,IP$Diagnosis..rank.11..,IP$Diagnosis..rank.12..,IP$Diagnosis..rank.13..,IP$Diagnosis..rank.14..,IP$Diagnosis..rank.15..)),]

cohort_SV <- cohort_SV[!is.na(cohort_SV$Sex.) & !is.na(cohort_SV$Date.of.Birth..yyyy.mm.dd..),]
cohort_ENA <- cohort_ENA[!is.na(cohort_ENA$Sex.) & !is.na(cohort_ENA$Date.of.Birth..yyyy.mm.dd..),]

chf_dx <- merge(DX, pts_ENA, by="Reference.Key.")
chf_dx <- chf_dx[chf_dx$Reference.Date. <= chf_dx$Prescription.Start.Date.,]
chf_ip <- merge(IP, pts_ENA, by="Reference.Key.")
chf_ip <- chf_ip[chf_ip$Admission.Date..yyyy.mm.dd.. <= chf_ip$Prescription.Start.Date.,]
pts_chf <- unique(c(chf_dx$Reference.Key., chf_ip$Reference.Key.))
cohort_ENA_CHF <- cohort_ENA[cohort_ENA$Reference.Key. %in% pts_chf,]

pts_SV$Index.Date <- pts_SV$Prescription.Start.Date.
ind_SV <- cbind.data.frame(Reference.Key.=pts_SV$Reference.Key., Index.Date=pts_SV$Index.Date)
pts_ENA$Index.Date <- pts_ENA$Prescription.Start.Date.
ind_ENA <- cbind.data.frame(Reference.Key.=pts_ENA$Reference.Key., Index.Date=pts_ENA$Index.Date)

final_ENA <- merge(cohort_ENA_CHF, ind_ENA, by="Reference.Key.")
final_SV <- merge(cohort_SV, ind_SV, by="Reference.Key.")


ENA_ <- Enalapril_rx[ENA_grep(Enalapril_rx),]
ENA_ <- ENA_[!grepl("ZENAPAX|ACETOMENAPHTHONE|CO-RENITEC", ENA_$Drug.Name.),]
ENA_ <- ENA_[complete.cases(ENA_),]
RX_ENA <- merge(ENA_, ind_ENA[ind_ENA$Reference.Key. %in% final_ENA$Reference.Key.,], by="Reference.Key.")
RX_ENA <- RX_ENA[(RX_ENA$Prescription.Start.Date. < RX_ENA$Index.Date) & (RX_ENA$Prescription.Start.Date. > RX_ENA$Index.Date - 365),]

SV_ENA <- Entresto_Rx[ENA_grep(Entresto_Rx),]
SV_ENA <- SV_ENA[!grepl("ZENAPAX|ACETOMENAPHTHONE|CO-RENITEC", SV_ENA$Drug.Name.),]
SV_ENA <- SV_ENA[complete.cases(SV_ENA),]
RX_SV <- merge(SV_ENA, ind_SV[ind_SV$Reference.Key. %in% final_SV$Reference.Key.,], by="Reference.Key.")
RX_SV <- RX_SV[(RX_SV$Prescription.Start.Date. < RX_SV$Index.Date) & (RX_SV$Prescription.Start.Date. > RX_SV$Index.Date - 365),]

final_ENA <- final_ENA[!final_ENA$Reference.Key. %in% RX_ENA$Reference.Key.,]
final_SV <- final_SV[!final_SV$Reference.Key. %in% RX_SV$Reference.Key.,]
final_ENA$Age <- (final_ENA$Index.Date - final_ENA$Date.of.Birth..yyyy.mm.dd..) / 365.25
final_SV$Age <- (final_SV$Index.Date - final_SV$Date.of.Birth..yyyy.mm.dd..) / 365.25
final_ENA <- final_ENA[final_ENA$Age > 18,]
final_SV <- final_SV[final_SV$Age > 18,]

final_pts_ENA <- cbind.data.frame(Reference.Key.=final_ENA$Reference.Key., Index.Date=final_ENA$Index.Date)
final_pts_SV <- cbind.data.frame(Reference.Key.=final_SV$Reference.Key., Index.Date=final_SV$Index.Date)

final_pts_ENA <- unique(final_pts_ENA)
final_pts_SV <- unique(final_pts_SV)

final_pts_ENA <- merge(final_pts_ENA, DEMO[, c("Reference.Key.", "Sex.", "Date.of.Birth..yyyy.mm.dd..", "Date.of.Registered.Death.")], by="Reference.Key.")
final_pts_SV <- merge(final_pts_SV, DEMO[, c("Reference.Key.", "Sex.", "Date.of.Birth..yyyy.mm.dd..", "Date.of.Registered.Death.")], by="Reference.Key.")

printcount(final_pts_ENA)
printcount(final_pts_SV)
printcount(final_ENA)
printcount(final_SV)
save(final_pts_ENA, final_pts_SV, final_ENA, final_SV, file="final_cohort.RData")

