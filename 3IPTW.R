
load("2. Baseline Characteristics/final_baseline.RData")
baseline_ENA$Index.Drug <- 0
baseline_SV$Index.Drug <- 1
baseline <- rbind(baseline_SV, baseline_ENA)
rm(list=c("baseline_ENA", "baseline_SV"))

load("2. Baseline Characteristics/supp_meds.RData")
meds$Index.Date <- NULL
meds$Index.Drug <- NULL
baseline <- merge(baseline, meds, by="Reference.Key.")

baseline$Date.of.Birth..yyyy.mm.dd.. <- NULL
baseline$Date.of.Registered.Death. <- NULL
for(name in grep("com[.]", names(baseline), value = TRUE)) {
  baseline[,name] <- as.factor(baseline[,name])
}
rm(name)
for(name in grep("meds[.]", names(baseline), value = TRUE)) {
  baseline[,name] <- as.factor(baseline[,name])
}
rm(name)


# IPTW ----------------------------------------------------------------------------
library(WeightIt)
W <- weightit(Index.Drug ~ Age + Sex. + com.icm + com.htn + com.dm + com.af + com.ihd + com.mi + com.ist + cci + Admissions.last.year + Admissions.HF + meds.arrthy + meds.bbk + meds.ccb + meds.digoxin + meds.diuretic + meds.mra + meds.ras + meds.statins + meds.throm + meds.dm,
                  data = baseline, estimand = "ATE", method = "ps")
#W.trim <- trim(W, at=.99, lower=T)
baseline$wts <- W$weights


# GLM -----------------------------------------------------------------------------
PS <- glm(Index.Drug ~ Age + Sex. + com.icm + com.htn + com.dm + com.af + com.ihd + com.mi + com.ist + cci + Admissions.last.year + Admissions.HF + meds.arrthy + meds.bbk + meds.ccb + meds.digoxin + meds.diuretic + meds.mra + meds.ras + meds.statins + meds.throm + meds.dm,
          data = baseline, family = "binomial")
baseline$ps <- PS$fitted.values

# Table 1 --------------------------------------------------------------------------
library(survey)
baseline.design <- svydesign(ids = ~ 1, data = baseline, weights = ~ wts)

library(tableone)
printtable <- function(df) {
  tbl1 <- svyCreateTableOne(c("Age", "Sex.", "com.icm", "com.htn", "com.dm", "com.af", "com.ihd", "com.mi", "com.ist", "cci", "Admissions.last.year", "Admissions.HF", "meds.arrthy", "meds.bbk", "meds.ccb", "meds.digoxin", "meds.diuretic", "meds.mra", "meds.ras", "meds.statins", "meds.throm", "meds.dm"), 
                            c("Index.Drug"), df, test=FALSE)
  print(tbl1, smd=TRUE, test=FALSE)
}

printtable(baseline.design)
save(baseline, file="final_iptw.RData")

#fits.ipw.all.svy <- fits(baseline.ipw.design)

tbl2 <- CreateTableOne(c("Age", "Sex.", "com.icm", "com.htn", "com.dm", "com.af", "com.ihd", "com.mi", "com.ist", "cci", "Admissions.last.year", "Admissions.HF", "meds.arrthy", "meds.bbk", "meds.ccb", "meds.digoxin", "meds.diuretic", "meds.mra", "meds.ras", "meds.statins", "meds.throm", "meds.dm"), 
                          c("Index.Drug"), baseline, test=FALSE)
print(tbl2, smd=TRUE, test=FALSE)
