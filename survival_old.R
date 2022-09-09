
load("C:/Users/ProfWongRA6.ProfWongRA6-PC/Desktop/Vincent/Entresto Enalapril/4. Main Analysis/final_events.RData")

#load("C:/Users/ProfWongRA6.ProfWongRA6-PC/Desktop/Vincent/Entresto Enalapril/4. Main Analysis/final_events_CHFonly.RData")


# ONLY FOR AFTER MATCHING
load("Y:/Swathi/Entresto/Data/Cross-check/FINAL_base and matched.RData")
matched_pts <- matched_model$Reference.Key.
rm(list = c("base_model", "index", "matched_model"))
baseline <- baseline[baseline$Reference.Key. %in% matched_pts,]
#########################

baseline$Index.Drug.Name <- as.factor(ifelse(baseline$Index.Drug, "Entresto", "Enalapril"))

#baseline$Date.censored <- pmin(baseline$Admission.first.any, baseline$Admission.first.HF, baseline$Death.any, baseline$Death.cv, baseline$Date.Switched, baseline$Date.Discontd, baseline$Date.StudyEnd, na.rm = T)
baseline$Censored.Admission.first.any <- pmin(baseline$Admission.first.any, baseline$Death.any, baseline$Date.Switched, baseline$Date.Discontd, baseline$Date.StudyEnd, na.rm = T)
baseline$Censored.Admission.first.HF <- pmin(baseline$Admission.first.HF, baseline$Death.any, baseline$Date.Switched, baseline$Date.Discontd, baseline$Date.StudyEnd, na.rm = T)
baseline$Censored.Death.any <- pmin(baseline$Death.any, baseline$Date.Switched, baseline$Date.Discontd, baseline$Date.StudyEnd, na.rm = T)
baseline$Censored.Death.cv <- pmin(baseline$Death.cv, baseline$Death.any, baseline$Date.Switched, baseline$Date.Discontd, baseline$Date.StudyEnd, na.rm = T)
baseline$Censored.Composite <- pmin(baseline$Censored.Admission.first.HF, baseline$Censored.Death.cv)

baseline$Outcome.Admission.first.any <- !is.na(baseline$Admission.first.any) & baseline$Admission.first.any<=baseline$Date.StudyEnd
baseline$Outcome.Admission.first.HF <- !is.na(baseline$Admission.first.HF) & baseline$Admission.first.HF<=baseline$Date.StudyEnd
baseline$Outcome.Death.any <- !is.na(baseline$Death.any) & baseline$Death.any<=baseline$Date.StudyEnd
baseline$Outcome.Death.cv <- !is.na(baseline$Death.cv) & baseline$Death.cv<=baseline$Date.StudyEnd
baseline$Outcome.Composite <- baseline$Outcome.Admission.first.HF | baseline$Outcome.Death.cv

baseline$Time.to.Admission.first.any <- as.numeric(baseline$Censored.Admission.first.any) - as.numeric(baseline$Index.Date) + 1
baseline$Time.to.Admission.first.HF <- as.numeric(baseline$Censored.Admission.first.HF) - as.numeric(baseline$Index.Date) + 1
baseline$Time.to.Death.any <- as.numeric(baseline$Censored.Death.any) - as.numeric(baseline$Index.Date) + 1
baseline$Time.to.Death.cv <- as.numeric(baseline$Censored.Death.cv) - as.numeric(baseline$Index.Date) + 1
baseline$Time.to.Composite <- as.numeric(baseline$Censored.Composite) - as.numeric(baseline$Index.Date) + 1

#baseline$wts <- 1

library(survival)

printhr <- function(df, name) {
  num_evt_tx <- sum(df[, paste0("Outcome.", name)] & df$Index.Drug)
  num_evt_ctrl <- sum(df[, paste0("Outcome.", name)] & !df$Index.Drug)
  if(is.null(df$wts)) fit <- coxph(Surv(get(paste0("Time.to.", name)), get(paste0("Outcome.", name))) ~ Index.Drug.Name, df)
  else fit <- coxph(Surv(get(paste0("Time.to.", name)), get(paste0("Outcome.", name))) ~ Index.Drug.Name + cluster(Reference.Key.), df, df$wts)
  hr <- as.numeric(exp(fit$coefficients[1]))
  ci <- exp(confint(fit))
  print(sprintf("%-25s - #events ctrl=%i tx=%i | HR %f (%f - %f)", name, num_evt_ctrl, num_evt_tx, hr, ci[1,1], ci[1,2]))
  return(fit)
}

printhrsvy <- function(df, name) {
  design <- df
  df <- design$variables
  num_evt_tx <- sum(df[, paste0("Outcome.", name)] & df$Index.Drug)
  num_evt_ctrl <- sum(df[, paste0("Outcome.", name)] & !df$Index.Drug)
  fit <- svycoxph(as.formula(paste0("Surv(", paste0("Time.to.", name), ", ", paste0("Outcome.", name), ") ~ Index.Drug.Name")), design=design)
  hr <- as.numeric(exp(fit$coefficients[1]))
  ci <- exp(confint(fit))
  print(sprintf("%-25s - #events ctrl=%i tx=%i | HR %f (%f - %f)", name, num_evt_ctrl, num_evt_tx, hr, ci[1,1], ci[1,2]))
  return(fit)
}

printhrcov <- function(df, name) {
  num_evt_tx <- sum(df[, paste0("Outcome.", name)] & df$Index.Drug)
  num_evt_ctrl <- sum(df[, paste0("Outcome.", name)] & !df$Index.Drug)
  fit <- coxph(Surv(get(paste0("Time.to.", name)), get(paste0("Outcome.", name))) ~ Index.Drug.Name + Age + Sex. + com.icm + com.htn + com.dm + com.af + com.ihd + com.mi + com.ist + cci + Admissions.last.year + Admissions.HF, df)
  hr <- as.numeric(exp(fit$coefficients[1]))
  ci <- exp(confint(fit))
  print(sprintf("%-25s - #events ctrl=%i tx=%i | HR %f (%f - %f)", name, num_evt_ctrl, num_evt_tx, hr, ci[1,1], ci[1,2]))
  return(fit)
}

fits <- function(df) {
  fun <- printhr
  if(grepl("survey.design", class(df)[1])) fun <- printhrsvy
  fit <- list()
  fit[["Composite"]] <- fun(df, "Composite")
  fit[["Admission.first.any"]] <- fun(df, "Admission.first.any")
  fit[["Admission.first.HF"]] <- fun(df, "Admission.first.HF")
  fit[["Death.any"]] <- fun(df, "Death.any")
  fit[["Death.cv"]] <- fun(df, "Death.cv")
  return(fit)
}

fits.all <- fits(baseline)

baseline.young <- baseline[baseline$Age < 65, ]
baseline.old <- baseline[baseline$Age >= 65, ]
print(sprintf("Age<65: %i, Age>=65: %i", nrow(baseline.young), nrow(baseline.old)))
fits.old <- fits(baseline.old)
fits.young <- fits(baseline.young)

#kaplan meier-----------------------------------------------------------------------------------------------

library(survminer)

plotkm <- function(df, name) {
  fm <- as.formula(paste0("Surv(", paste0("Time.to.", name), ", ", paste0("Outcome.", name), ") ~ Index.Drug.Name"))
  fit <- survfit(fm, data=df)
  fit$call$formula <- fm
  summary(fit)
  ggsurvplot(fit, data=df, conf.int=T, pval=T, risk.table=T, legend.title="Drug", legend.labs=levels(df$Index.Drug.Name))
}

plots <- function(df) {
  fun <- plotkm
  plot <- list()
  plot[["Composite"]] <- fun(df, "Composite")
  plot[["Admission.first.any"]] <- fun(df, "Admission.first.any")
  plot[["Admission.first.HF"]] <- fun(df, "Admission.first.HF")
  plot[["Death.any"]] <- fun(df, "Death.any")
  plot[["Death.cv"]] <- fun(df, "Death.cv")
  return(plot)
}

plots.all <- plots(baseline)

# Sensitivity analysis (IPTW) ------------------------------------------------------------------------------

library(ipw)
#wts <- ipwpoint(Index.Drug, "binomial", "logit", numerator = ~ 1, denominator = ~ Age + Sex. + com.icm + com.htn + com.dm + com.af + com.ihd + com.mi + com.ist + cci + Admissions.last.year + Admissions.HF, baseline)
wts.trunc <- ipwpoint(Index.Drug, "binomial", "logit", numerator = ~ 1, denominator = ~ Age + Sex. + com.icm + com.htn + com.dm + com.af + com.ihd + com.mi + com.ist + cci + Admissions.last.year + Admissions.HF, baseline,
                      trunc=0.05)
baseline.ipw <- baseline
baseline.ipw$wts <- wts.trunc$weights.trunc

fits.ipw.all <- fits(baseline.ipw)
baseline.ipw.young <- baseline.ipw[baseline.ipw$Age < 65, ]
baseline.ipw.old <- baseline.ipw[baseline.ipw$Age >= 65, ]
fits.ipw.old <- fits(baseline.ipw.old)
fits.ipw.young <- fits(baseline.ipw.young)

library(survey)
baseline.ipw.design <- svydesign(ids = ~ 1, data = baseline.ipw, weights = ~ wts)
fits.ipw.all.svy <- fits(baseline.ipw.design)

library(tableone)
printtable <- function(df) {
  tbl1 <- svyCreateTableOne(c("Age", "Sex.", "com.icm", "com.htn", "com.dm", "com.af", "com.ihd", "com.mi", "com.ist", "cci", "Admissions.last.year", "Admissions.HF"), 
                         c("Index.Drug"), df, test=FALSE)
  print(tbl1, smd=TRUE, test=FALSE)
}

printtable(baseline.ipw.design)


# Readmission analysis -------------------------------------------------------------------------------------

library(pscl)
baseline[!is.na(baseline$Readmit.Index.Date) & !is.na(baseline$Date.Switched) & baseline$Readmit.Index.Date > baseline$Date.Switched,]$Readmit.Index.Date <- NA
baseline$Readmit.any.Censored <- pmin(baseline$Readmit.any.Date, baseline$Death.any, baseline$Date.Switched, baseline$Readmit.Index.Date+30, na.rm = T)
baseline$Readmit.HF.Censored <- pmin(baseline$Readmit.HF.Date, baseline$Death.any, baseline$Date.Switched, baseline$Readmit.Index.Date+30, na.rm = T)
baseline$Readmit.any.Incidence <-  as.integer(baseline$Readmit.any / as.numeric((baseline$Readmit.any.Censored-baseline$Readmit.Index.Date)/365.25) * 100)
baseline$Readmit.HF.Incidence <-  as.integer(baseline$Readmit.HF / as.numeric((baseline$Readmit.HF.Censored-baseline$Readmit.Index.Date)/365.25) * 100)
summary(zeroinfl(Readmit.any.Incidence ~ Index.Drug.Name, baseline[!is.na(baseline$Readmit.Index.Date),]))
summary(zeroinfl(Readmit.HF.Incidence ~ Index.Drug.Name, baseline[!is.na(baseline$Readmit.Index.Date),]))

