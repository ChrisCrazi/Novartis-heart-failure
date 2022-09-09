library(survminer)
library(tidyverse)
library(ggpubr)
library(patchwork)
library(readxl)

load("3. IPTW/final_iptw.RData")
baseline.iptw <- baseline
# rm(baseline)
# load("4. Main Analysis/final_events.RData")
load("4. Main Analysis/final_events_2021.RData") # Chris amendments

baseline$Index.Drug.Name <- as.factor(ifelse(baseline$Index.Drug, "Entresto", "Enalapril"))

# Time to event ---------------------------------------------------------------------------
# baseline$Date.censored <- pmin(baseline$Admission.first.any, baseline$Admission.first.HF, baseline$Death.any, baseline$Death.cv, baseline$Date.Switched, baseline$Date.Discontd, baseline$Date.StudyEnd, na.rm = T)
baseline$Censored.Admission.first.any <- pmin(baseline$Admission.first.any, baseline$Death.any, baseline$Date.Switched, baseline$Date.Discontd, baseline$Date.StudyEnd, na.rm = T)
baseline$Censored.Admission.first.HF <- pmin(baseline$Admission.first.HF, baseline$Death.any, baseline$Date.Switched, baseline$Date.Discontd, baseline$Date.StudyEnd, na.rm = T)
baseline$Censored.Death.any <- pmin(baseline$Death.any, baseline$Date.Switched, baseline$Date.Discontd, baseline$Date.StudyEnd, na.rm = T)
baseline$Censored.Death.cv <- pmin(baseline$Death.cv, baseline$Death.any, baseline$Date.Switched, baseline$Date.Discontd, baseline$Date.StudyEnd, na.rm = T)
baseline$Censored.Composite <- pmin(baseline$Censored.Admission.first.HF, baseline$Censored.Death.cv)

# # ITT
# baseline$Outcome.Admission.first.any <- !is.na(baseline$Admission.first.any) & baseline$Admission.first.any<=baseline$Date.StudyEnd
# baseline$Outcome.Admission.first.HF <- !is.na(baseline$Admission.first.HF) & baseline$Admission.first.HF<=baseline$Date.StudyEnd
# baseline$Outcome.Death.any <- !is.na(baseline$Death.any) & baseline$Death.any<=baseline$Date.StudyEnd
# baseline$Outcome.Death.cv <- !is.na(baseline$Death.cv) & baseline$Death.cv<=baseline$Date.StudyEnd
# baseline$Outcome.Composite <- baseline$Outcome.Admission.first.HF | baseline$Outcome.Death.cv

# per protocol
baseline$Outcome.Admission.first.any <- !is.na(baseline$Admission.first.any) & baseline$Admission.first.any<=baseline$Date.StudyEnd & (baseline$Admission.first.any<=baseline$Date.Switched|is.na(baseline$Date.Switched)) & (baseline$Admission.first.any<=baseline$Date.Discontd|is.na(baseline$Date.Discontd))
baseline$Outcome.Admission.first.HF <- !is.na(baseline$Admission.first.HF) & baseline$Admission.first.HF<=baseline$Date.StudyEnd & (baseline$Admission.first.HF<=baseline$Date.Switched|is.na(baseline$Date.Switched)) & (baseline$Admission.first.HF<=baseline$Date.Discontd|is.na(baseline$Date.Discontd))
baseline$Outcome.Death.any <- !is.na(baseline$Death.any) & baseline$Death.any<=baseline$Date.StudyEnd & (baseline$Death.any<=baseline$Date.Switched|is.na(baseline$Date.Switched)) & (baseline$Death.any<=baseline$Date.Discontd|is.na(baseline$Date.Discontd))
baseline$Outcome.Death.cv <- !is.na(baseline$Death.cv) & baseline$Death.cv<=baseline$Date.StudyEnd & (baseline$Death.cv<=baseline$Date.Switched|is.na(baseline$Date.Switched)) & (baseline$Death.cv<=baseline$Date.Discontd|is.na(baseline$Date.Discontd))
baseline$Outcome.Composite <- baseline$Outcome.Admission.first.HF | baseline$Outcome.Death.cv

# # per protocol same censor
# baseline$Outcome.Admission.first.any <- !is.na(baseline$Admission.first.any) & baseline$Admission.first.any<=baseline$Date.censored
# baseline$Outcome.Admission.first.HF <- !is.na(baseline$Admission.first.HF) & baseline$Admission.first.HF<=baseline$Date.censored
# baseline$Outcome.Death.any <- !is.na(baseline$Death.any) & baseline$Death.any<=baseline$Date.censored
# baseline$Outcome.Death.cv <- !is.na(baseline$Death.cv) & baseline$Death.cv<=baseline$Date.censored
# baseline$Outcome.Composite <- baseline$Outcome.Admission.first.HF | baseline$Outcome.Death.cv


baseline$Time.to.Admission.first.any <- as.numeric(baseline$Censored.Admission.first.any) - as.numeric(baseline$Index.Date) + 1
baseline$Time.to.Admission.first.HF <- as.numeric(baseline$Censored.Admission.first.HF) - as.numeric(baseline$Index.Date) + 1
baseline$Time.to.Death.any <- as.numeric(baseline$Censored.Death.any) - as.numeric(baseline$Index.Date) + 1
baseline$Time.to.Death.cv <- as.numeric(baseline$Censored.Death.cv) - as.numeric(baseline$Index.Date) + 1
baseline$Time.to.Composite <- as.numeric(baseline$Censored.Composite) - as.numeric(baseline$Index.Date) + 1

# IPTW weights ----------------------------------------------------------------------------
baseline <- merge(baseline, baseline.iptw[,c("Reference.Key.", "wts", "ps", "meds.arrthy", "meds.bbk", "meds.ccb", "meds.digoxin",
                                             "meds.diuretic", "meds.dm", "meds.mra", "meds.ras", "meds.statins", "meds.throm")], by="Reference.Key.")
rm(baseline.iptw)

# Survival analysis: Main -----------------------------------------------------------------
library(survival)

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

fits <- function(df, cov=F) {
  if(grepl("survey.design", class(df)[1])) fun <- printhrsvy
  else fun <- printhr
  if(cov) fun <- printhrcov
  fit <- list()
  fit[["Composite"]] <- fun(df, "Composite")
  fit[["Admission.first.any"]] <- fun(df, "Admission.first.any")
  fit[["Admission.first.HF"]] <- fun(df, "Admission.first.HF")
  fit[["Death.any"]] <- fun(df, "Death.any")
  fit[["Death.cv"]] <- fun(df, "Death.cv")
  return(fit)
}

library(survey)
baseline.design <- svydesign(ids = ~ 1, data = baseline) # crude result
fits.all <- fits(baseline.design)

baseline.design <- svydesign(ids = ~ 1, data = baseline, weights = ~ wts)
fits.all <- fits(baseline.design)

fit1<- survfit(Surv(Time.to.Composite, Outcome.Composite) ~ Index.Drug.Name, data = baseline)
fit2<- survfit(Surv(Time.to.Admission.first.any, Outcome.Admission.first.any) ~ Index.Drug.Name, data = baseline)
fit3<- survfit(Surv(Time.to.Admission.first.HF, Outcome.Admission.first.HF) ~ Index.Drug.Name, data = baseline)
fit4<- survfit(Surv(Time.to.Death.any, Outcome.Death.any) ~ Index.Drug.Name, data = baseline)
fit5<- survfit(Surv(Time.to.Death.cv, Outcome.Death.cv) ~ Index.Drug.Name, data = baseline)
plot1 <- ggsurvplot(fit1, data = baseline, xlab = "Time in days", legend.labs = c("Enalapril", "Sacubitril/valsartan"), title = "Composite of cardiovascular mortality\n and heart failure-related hospitalization")
plot2 <- ggsurvplot(fit2, data = baseline, xlab = "Time in days", legend.labs = c("Enalapril", "Sacubitril/valsartan"), title = "All-cause hospitalization")
plot3 <- ggsurvplot(fit3, data = baseline, xlab = "Time in days", legend.labs = c("Enalapril", "Sacubitril/valsartan"), title = "Heart failure-related hospitalization")
plot4 <- ggsurvplot(fit4, data = baseline, xlab = "Time in days", legend.labs = c("Enalapril", "Sacubitril/valsartan"), title = "All-cause mortality")
plot5 <- ggsurvplot(fit5, data = baseline, xlab = "Time in days", legend.labs = c("Enalapril", "Sacubitril/valsartan"), title = "Cardiovascular mortality")
layoutplot <- "
aabc
aade
"
plotlist <- list(a = plot1$plot, b = plot2$plot, c = plot3$plot, d = plot4$plot, e = plot5$plot)
tiff("result.tif", width = 2400, height = 1800, res = 160)
wrap_plots(plotlist, guides = 'collect', design = layoutplot)
dev.off()

# Subgroups: ejection fraction
HN_ENTRESTO <- read_excel("C:/Users/LabPC14CSMPR/Desktop/Chris/heart failure cohort (Novartis)/HN number/HN_ENTRESTO.xlsx")
HN_ENALAPRIL <- read_excel("C:/Users/LabPC14CSMPR/Desktop/Chris/heart failure cohort (Novartis)/HN number/HN_ENALAPRIL.xlsx")
HN_ENTRESTO <- HN_ENTRESTO %>% rename("ejection.fraction" = `ejection fraction`)
HN_ENTRESTO <- HN_ENTRESTO %>% mutate(ejection.fraction = as.numeric(ejection.fraction)) %>% select(reference.key, ejection.fraction)
HN_ENALAPRIL <- HN_ENALAPRIL %>% rename("ejection.fraction" = `ejection fraction`)
HN_ENALAPRIL <- HN_ENALAPRIL %>% mutate(ejection.fraction = as.numeric(ejection.fraction)) %>% select(reference.key, ejection.fraction)
baseline <- baseline %>% mutate(ejection.fraction.c = if_else(Reference.Key. %in% filter(rbind(HN_ENTRESTO, HN_ENALAPRIL), ejection.fraction <= 40)$reference.key, "r", "")) %>% 
  mutate(ejection.fraction.c = if_else(Reference.Key. %in% filter(rbind(HN_ENTRESTO, HN_ENALAPRIL), ejection.fraction > 40 | is.na(ejection.fraction))$reference.key, "p", ejection.fraction.c))

baseline.design <- svydesign(ids = ~ 1, data = filter(baseline, ejection.fraction.c == "r"|Index.Drug == 0), weights = ~ wts)
fits.all <- fits(baseline.design)

baseline.design <- svydesign(ids = ~ 1, data = filter(baseline, ejection.fraction.c == "p"|Index.Drug == 0), weights = ~ wts)
fits.all <- fits(baseline.design)

# Subgroups: cardio.specialty
library(tidyverse)
library(tidylog)
op_combined <- readRDS("Data/Chris follow up/combined/2016_2021_op_combined.rds")
op_patients <- op_combined %>% 
  filter(op.eis.sub.specialty == "CARDIO") %>%
  # filter(str_detect(drug.name, "ENTRESTO")) %>%
  .$reference.key %>% unique()
baseline <- baseline %>% 
  mutate(cardio.specialty = if_else(Reference.Key. %in% op_patients, "Y", "N"))

filter(baseline, cardio.specialty == "Y") %>% filter(Index.Drug.Name == "Entresto") %>% .$Outcome.Composite %>% sum
filter(baseline, cardio.specialty == "Y") %>% filter(Index.Drug.Name == "Entresto") %>% .$Time.to.Composite %>% sum
filter(baseline, cardio.specialty == "Y") %>% filter(Index.Drug.Name == "Enalapril") %>% .$Outcome.Composite %>% sum
filter(baseline, cardio.specialty == "Y") %>% filter(Index.Drug.Name == "Enalapril") %>% .$Time.to.Composite %>% sum

filter(baseline, cardio.specialty == "N") %>% filter(Index.Drug.Name == "Entresto") %>% .$Outcome.Composite %>% sum
filter(baseline, cardio.specialty == "N") %>% filter(Index.Drug.Name == "Entresto") %>% .$Time.to.Composite %>% sum
filter(baseline, cardio.specialty == "N") %>% filter(Index.Drug.Name == "Enalapril") %>% .$Outcome.Composite %>% sum
filter(baseline, cardio.specialty == "N") %>% filter(Index.Drug.Name == "Enalapril") %>% .$Time.to.Composite %>% sum

baseline.design <- svydesign(ids = ~ 1, data = filter(baseline, cardio.specialty == "Y"), weights = ~ wts)
fits.all <- fits(baseline.design)

baseline.design <- svydesign(ids = ~ 1, data = filter(baseline, cardio.specialty == "N"), weights = ~ wts)
fits.all <- fits(baseline.design)

fit1<- survfit(Surv(Time.to.Composite, Outcome.Composite) ~ Index.Drug.Name, data = filter(baseline, cardio.specialty == "Y"))
fit2<- survfit(Surv(Time.to.Admission.first.any, Outcome.Admission.first.any) ~ Index.Drug.Name, data = filter(baseline, cardio.specialty == "Y"))
fit3<- survfit(Surv(Time.to.Admission.first.HF, Outcome.Admission.first.HF) ~ Index.Drug.Name, data = filter(baseline, cardio.specialty == "Y"))
fit4<- survfit(Surv(Time.to.Death.any, Outcome.Death.any) ~ Index.Drug.Name, data = filter(baseline, cardio.specialty == "Y"))
fit5<- survfit(Surv(Time.to.Death.cv, Outcome.Death.cv) ~ Index.Drug.Name, data = filter(baseline, cardio.specialty == "Y"))
plot1 <- ggsurvplot(fit1, data = filter(baseline, cardio.specialty == "Y"), xlab = "Time in days", legend.labs = c("Enalapril", "Sacubitril/valsartan"), title = "Composite of cardiovascular mortality\n and heart failure-related hospitalization")
plot2 <- ggsurvplot(fit2, data = filter(baseline, cardio.specialty == "Y"), xlab = "Time in days", legend.labs = c("Enalapril", "Sacubitril/valsartan"), title = "All-cause hospitalization")
plot3 <- ggsurvplot(fit3, data = filter(baseline, cardio.specialty == "Y"), xlab = "Time in days", legend.labs = c("Enalapril", "Sacubitril/valsartan"), title = "Heart failure-related hospitalization")
plot4 <- ggsurvplot(fit4, data = filter(baseline, cardio.specialty == "Y"), xlab = "Time in days", legend.labs = c("Enalapril", "Sacubitril/valsartan"), title = "All-cause mortality")
plot5 <- ggsurvplot(fit5, data = filter(baseline, cardio.specialty == "Y"), xlab = "Time in days", legend.labs = c("Enalapril", "Sacubitril/valsartan"), title = "Cardiovascular mortality")
fit6<- survfit(Surv(Time.to.Composite, Outcome.Composite) ~ Index.Drug.Name, data = filter(baseline, cardio.specialty == "N"))
fit7<- survfit(Surv(Time.to.Admission.first.any, Outcome.Admission.first.any) ~ Index.Drug.Name, data = filter(baseline, cardio.specialty == "N"))
fit8<- survfit(Surv(Time.to.Admission.first.HF, Outcome.Admission.first.HF) ~ Index.Drug.Name, data = filter(baseline, cardio.specialty == "N"))
fit9<- survfit(Surv(Time.to.Death.any, Outcome.Death.any) ~ Index.Drug.Name, data = filter(baseline, cardio.specialty == "N"))
fit10<- survfit(Surv(Time.to.Death.cv, Outcome.Death.cv) ~ Index.Drug.Name, data = filter(baseline, cardio.specialty == "N"))
plot6 <- ggsurvplot(fit6, data = filter(baseline, cardio.specialty == "N"), xlab = "Time in days", legend.labs = c("Enalapril", "Sacubitril/valsartan"), title = "Composite of cardiovascular mortality\n and heart failure-related hospitalization")
plot7 <- ggsurvplot(fit7, data = filter(baseline, cardio.specialty == "N"), xlab = "Time in days", legend.labs = c("Enalapril", "Sacubitril/valsartan"), title = "All-cause hospitalization")
plot8 <- ggsurvplot(fit8, data = filter(baseline, cardio.specialty == "N"), xlab = "Time in days", legend.labs = c("Enalapril", "Sacubitril/valsartan"), title = "Heart failure-related hospitalization")
plot9 <- ggsurvplot(fit9, data = filter(baseline, cardio.specialty == "N"), xlab = "Time in days", legend.labs = c("Enalapril", "Sacubitril/valsartan"), title = "All-cause mortality")
plot10 <- ggsurvplot(fit10, data = filter(baseline, cardio.specialty == "N"), xlab = "Time in days", legend.labs = c("Enalapril", "Sacubitril/valsartan"), title = "Cardiovascular mortality")
col1 <- ggplot() + annotate(geom = 'text', size=10, x=1, y=1, label="Managed by cardiology specialist") + theme_void() 
col2 <- ggplot() + annotate(geom = 'text', size=10, x=1, y=1, label="Not managed by cardiology specialist") + theme_void() 
layoutplot <- "
xxxx
aabc
aade
yyyy
ffgh
ffij
"
plotlist <- list(a = plot1$plot, b = plot2$plot, c = plot3$plot, d = plot4$plot, e = plot5$plot,
                 f = plot6$plot, g = plot7$plot, h = plot8$plot, i = plot9$plot, j = plot10$plot,
                 x = col1, y = col2)
tiff("result1.tif", width = 2400, height = 2400, res = 160)
wrap_plots(plotlist, guides = 'collect', design = layoutplot)
dev.off()




# Subgroups: recurrent and new onset HF
ICDCM_CHF <- "398.91|402.01|402.11|402.91|404.01|404.03|404.11|404.13|404.91|404.93|428"
load("Data/Combined data/Ddx_enalapril.Rdata")
load("Data/Combined data/Ddx_entresto.Rdata")
dx_combined <- readRDS("Data/Chris follow up/combined/2019_2021_dx_combined.rds")
recurrent_patients <- unique(rbind(entresto_ddx, mutate(enalapril_ddx, Reference.Date. = as.Date(Reference.Date.)), dx_combined)) %>% filter(str_detect(All.Diagnosis.Code..ICD9.., ICDCM_CHF)) %>% 
  filter(Reference.Key. %in% baseline$Reference.Key.) %>% 
  mutate(Index.Date = baseline$Index.Date[match(Reference.Key., baseline$Reference.Key.)]) %>% 
  filter(Reference.Date. <= Index.Date) %>% 
  filter(Index.Date - Reference.Date. > 30) %>% 
  select(Reference.Date., Index.Date, everything()) %>% .$Reference.Key. %>% unique()
baseline <- baseline %>% 
  mutate(recurrent.patient = if_else(Reference.Key. %in% recurrent_patients, "Y", "N"))

baseline.design <- svydesign(ids = ~ 1, data = filter(baseline, recurrent.patient == "Y"), weights = ~ wts)
fits.all <- fits(baseline.design)

baseline.design <- svydesign(ids = ~ 1, data = filter(baseline, recurrent.patient == "N"), weights = ~ wts)
fits.all <- fits(baseline.design)

fit1<- survfit(Surv(Time.to.Composite, Outcome.Composite) ~ Index.Drug.Name, data = filter(baseline, recurrent.patient == "Y"))
fit2<- survfit(Surv(Time.to.Admission.first.any, Outcome.Admission.first.any) ~ Index.Drug.Name, data = filter(baseline, recurrent.patient == "Y"))
fit3<- survfit(Surv(Time.to.Admission.first.HF, Outcome.Admission.first.HF) ~ Index.Drug.Name, data = filter(baseline, recurrent.patient == "Y"))
fit4<- survfit(Surv(Time.to.Death.any, Outcome.Death.any) ~ Index.Drug.Name, data = filter(baseline, recurrent.patient == "Y"))
fit5<- survfit(Surv(Time.to.Death.cv, Outcome.Death.cv) ~ Index.Drug.Name, data = filter(baseline, recurrent.patient == "Y"))
plot1 <- ggsurvplot(fit1, data = filter(baseline, recurrent.patient == "Y"), xlab = "Time in days", legend.labs = c("Enalapril", "Sacubitril/valsartan"), title = "Composite of cardiovascular mortality\n and heart failure-related hospitalization")
plot2 <- ggsurvplot(fit2, data = filter(baseline, recurrent.patient == "Y"), xlab = "Time in days", legend.labs = c("Enalapril", "Sacubitril/valsartan"), title = "All-cause hospitalization")
plot3 <- ggsurvplot(fit3, data = filter(baseline, recurrent.patient == "Y"), xlab = "Time in days", legend.labs = c("Enalapril", "Sacubitril/valsartan"), title = "Heart failure-related hospitalization")
plot4 <- ggsurvplot(fit4, data = filter(baseline, recurrent.patient == "Y"), xlab = "Time in days", legend.labs = c("Enalapril", "Sacubitril/valsartan"), title = "All-cause mortality")
plot5 <- ggsurvplot(fit5, data = filter(baseline, recurrent.patient == "Y"), xlab = "Time in days", legend.labs = c("Enalapril", "Sacubitril/valsartan"), title = "Cardiovascular mortality")
fit6<- survfit(Surv(Time.to.Composite, Outcome.Composite) ~ Index.Drug.Name, data = filter(baseline, recurrent.patient == "N"))
fit7<- survfit(Surv(Time.to.Admission.first.any, Outcome.Admission.first.any) ~ Index.Drug.Name, data = filter(baseline, recurrent.patient == "N"))
fit8<- survfit(Surv(Time.to.Admission.first.HF, Outcome.Admission.first.HF) ~ Index.Drug.Name, data = filter(baseline, recurrent.patient == "N"))
fit9<- survfit(Surv(Time.to.Death.any, Outcome.Death.any) ~ Index.Drug.Name, data = filter(baseline, recurrent.patient == "N"))
fit10<- survfit(Surv(Time.to.Death.cv, Outcome.Death.cv) ~ Index.Drug.Name, data = filter(baseline, recurrent.patient == "N"))
plot6 <- ggsurvplot(fit6, data = filter(baseline, recurrent.patient == "N"), xlab = "Time in days", legend.labs = c("Enalapril", "Sacubitril/valsartan"), title = "Composite of cardiovascular mortality\n and heart failure-related hospitalization")
plot7 <- ggsurvplot(fit7, data = filter(baseline, recurrent.patient == "N"), xlab = "Time in days", legend.labs = c("Enalapril", "Sacubitril/valsartan"), title = "All-cause hospitalization")
plot8 <- ggsurvplot(fit8, data = filter(baseline, recurrent.patient == "N"), xlab = "Time in days", legend.labs = c("Enalapril", "Sacubitril/valsartan"), title = "Heart failure-related hospitalization")
plot9 <- ggsurvplot(fit9, data = filter(baseline, recurrent.patient == "N"), xlab = "Time in days", legend.labs = c("Enalapril", "Sacubitril/valsartan"), title = "All-cause mortality")
plot10 <- ggsurvplot(fit10, data = filter(baseline, recurrent.patient == "N"), xlab = "Time in days", legend.labs = c("Enalapril", "Sacubitril/valsartan"), title = "Cardiovascular mortality")
col1 <- ggplot() + annotate(geom = 'text', size=10, x=1, y=1, label="Recurrent heart failure") + theme_void() 
col2 <- ggplot() + annotate(geom = 'text', size=10, x=1, y=1, label="New onset heart failure") + theme_void() 
layoutplot <- "
xxxx
aabc
aade
yyyy
ffgh
ffij
"
plotlist <- list(a = plot1$plot, b = plot2$plot, c = plot3$plot, d = plot4$plot, e = plot5$plot,
                 f = plot6$plot, g = plot7$plot, h = plot8$plot, i = plot9$plot, j = plot10$plot,
                 x = col1, y = col2)
tiff("result2.tif", width = 2400, height = 2400, res = 160)
wrap_plots(plotlist, guides = 'collect', design = layoutplot)
dev.off()



c# Subgroups: Age
baseline.young <- baseline[baseline$Age < 65,]
baseline.design.young <- svydesign(ids = ~ 1, data = baseline.young, weights = ~ wts)
fits.young <- fits(baseline.design.young)

baseline.old <- baseline[baseline$Age >= 65,]
baseline.design.old <- svydesign(ids = ~ 1, data = baseline.old, weights = ~ wts)
fits.old <- fits(baseline.design.old)


# Sesitivity analysis: PS matching ----------------------------------------------------
library(MatchIt)
MI <- matchit(Index.Drug ~ Age + Sex. + com.icm + com.htn + com.dm + com.af + com.ihd + com.mi + com.ist + cci + Admissions.last.year + Admissions.HF + meds.arrthy + meds.bbk + meds.ccb + meds.digoxin + meds.diuretic + meds.mra + meds.ras + meds.statins + meds.throm + meds.dm,
        data = baseline)
matched_pts <- match.data(MI)$Reference.Key. 
rm(matched_model)
baseline.ps <- baseline[baseline$Reference.Key. %in% matched_pts,]
baseline.ps$wts <- NULL

printhr <- function(df, name) {
  num_evt_tx <- sum(df[, paste0("Outcome.", name)] & df$Index.Drug)
  num_evt_ctrl <- sum(df[, paste0("Outcome.", name)] & !df$Index.Drug)
  if(is.null(df$wts)) fit <- coxph(Surv(get(paste0("Time.to.", name)), get(paste0("Outcome.", name))) ~ Index.Drug.Name, df)
  #else fit <- coxph(Surv(get(paste0("Time.to.", name)), get(paste0("Outcome.", name))) ~ Index.Drug.Name + cluster(Reference.Key.), df, df$wts)
  hr <- as.numeric(exp(fit$coefficients[1]))
  ci <- exp(confint(fit))
  print(sprintf("%-25s - #events ctrl=%i tx=%i | HR %f (%f - %f)", name, num_evt_ctrl, num_evt_tx, hr, ci[1,1], ci[1,2]))
  return(fit)
}

fits.ps.all <- fits(baseline.ps)

# Sesitivity analysis: Covariate adjustment ---
printhrcov <- function(df, name) {
  num_evt_tx <- sum(df[, paste0("Outcome.", name)] & df$Index.Drug)
  num_evt_ctrl <- sum(df[, paste0("Outcome.", name)] & !df$Index.Drug)
  fit <- coxph(Surv(get(paste0("Time.to.", name)), 
                    get(paste0("Outcome.", name))) ~ Index.Drug.Name + Age + Sex. + com.icm + com.htn + com.dm + com.af + com.ihd + com.mi + com.ist + cci + Admissions.last.year + Admissions.HF + meds.arrthy + meds.bbk + meds.ccb + meds.digoxin + meds.diuretic + meds.mra + meds.ras + meds.statins + meds.throm + meds.dm, 
                df)
  hr <- as.numeric(exp(fit$coefficients[1]))
  ci <- exp(confint(fit))
  print(sprintf("%-25s - #events ctrl=%i tx=%i | HR %f (%f - %f)", name, num_evt_ctrl, num_evt_tx, hr, ci[1,1], ci[1,2]))
  return(fit)
}

load("C:/Users/LabPC14CSMPR/Desktop/Chris/heart failure cohort (Novartis)/Entresto Enalapril_vincent/2. Baseline Characteristics/supp_meds.RData")
meds$Index.Date <- NULL
meds$Index.Drug <- NULL
baseline.cov <- merge(baseline, meds, by="Reference.Key.")
rm(meds)

fits.cov.all <- fits(baseline.cov, cov=T)

# Save ------------------------------------------------------------------------
save(fits.all, fits.old, fits.young, fits.ps.all, fits.cov.all, file="final_results.RData")


# Chris add, propensity score and cox regression
baseline


weightit(
         data = baseline, estimand = "ATE", method = "ps")