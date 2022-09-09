library(tidyverse)
library(lubridate)

load("C:/Users/LabPC14CSMPR/Desktop/Chris/heart failure cohort (Novartis)/Entresto Enalapril_Vincent_Chris/Data/Combined data/Acohort_entresto.Rdata")
load("C:/Users/LabPC14CSMPR/Desktop/Chris/heart failure cohort (Novartis)/Entresto Enalapril_Vincent_Chris/Data/Combined data/Ddx_entresto.Rdata")
load("C:/Users/LabPC14CSMPR/Desktop/Chris/heart failure cohort (Novartis)/Entresto Enalapril_Vincent_Chris/Data/Combined data/IP_entresto.Rdata")
load("4. Main Analysis/final_events_2021.RData") # Chris amendments

entresto_user <- Entresto %>% as_tibble() %>% arrange(Reference.Key., Prescription.Start.Date.) %>% group_by(Reference.Key.) %>% slice(1) %>% 
  ungroup() %>% mutate(index.date = Prescription.Start.Date.) %>% arrange(index.date) %>% 
  filter(Reference.Key. %in% filter(entresto_ddx, str_detect(All.Diagnosis.Code..ICD9.., "428"))$Reference.Key.)

dementia_patients <- filter(entresto_ddx, str_detect(All.Diagnosis.Code..ICD9.., "331.0|290.4|290.41|290.42|290.43|294.2|294.8|290.0|290.1|290.2|290.3|290.8|290.9|294.1|291.2|292.82|331.1|331.82")) %>% as_tibble() %>% 
  filter(Reference.Key. %in% entresto_user$Reference.Key.) %>% 
  mutate(index.date = entresto_user$index.date[match(Reference.Key., entresto_user$Reference.Key.)])

dementia_patients %>% arrange(Reference.Key., Reference.Date.) %>% group_by(Reference.Key.) %>% slice(1) %>% ungroup() %>% filter(Reference.Date. > index.date) %>% select(index.date, everything())

psychosis_codes <- "290.11|290.3|290.41|292.11|292.12|292.2|292.81|292.84|292.85|292.89|292.9|293.0|293.1|293.81|293.82|293.83|293.89|293.9|294.0|295.13|295.14|295.33|295.34|295.43|295.44|295.63|295.64|295.73|295.74|295.83|295.84|295.93|295.94|296.01|296.02|296.03|296.04|296.20|296.21|296.22|296.23|296.24|296.40|296.41|296.42|296.43|296.44|296.50|296.51|296.52|296.53|296.54|296.60|296.61|296.62|296.63|296.64|296.90|296.99|297.8|297.9|298.0|298.1|298.2|298.3|298.8|298.9|300.00|300.09|308.0|308.1|308.2|309.82|309.83|309.89|309.9|311|327.00|327.01|327.02|327.09|327.10|327.11|327.12|327.14|327.19|780.01|780.02|780.09|780.1|780.50|780.52|780.54|780.55|780.56|780.58|780.59|780.93|780.97|781.1|E930.3|E930.8|E930.9|E941.9|E947.8|E947.9|V40.0|V40.1|V40.2|V40.31|V40.39|V40.9"
psychosis_patients <- filter(Entresto_IP, str_detect(Principal.Diagnosis.Code., psychosis_codes)) %>% as_tibble() %>% 
  filter(Reference.Key. %in% entresto_user$Reference.Key.) %>% 
  mutate(index.date = entresto_user$index.date[match(Reference.Key., entresto_user$Reference.Key.)])

psychosis_patients %>% arrange(Reference.Key., Admission.Date..yyyy.mm.dd..) %>% group_by(Reference.Key.) %>% slice(1) %>% ungroup() %>% 
  select(index.date, Admission.Date..yyyy.mm.dd.., everything()) %>% filter(Admission.Date..yyyy.mm.dd.. > index.date)



entresto_rx <- Entresto %>% as_tibble() %>% select(-Drug.Item.Code., -Type.of.Patient..Drug..) %>% group_by(Reference.Key.) %>%
  arrange(Prescription.Start.Date.) %>%
  filter(as.numeric(Prescription.End.Date.) >= cummax(as.numeric(Prescription.End.Date.))) %>%
  group_by(Reference.Key., Prescription.End.Date.) %>%
  filter(Prescription.Start.Date. == min(Prescription.Start.Date.)) %>%
  slice(1) %>%
  ungroup() %>% 
  group_by(Reference.Key.) %>%
  mutate(previous.end.date = lag(Prescription.End.Date.)) %>%
  mutate(Prescription.Start.Date. = if_else(!is.na(previous.end.date)&as.numeric(Prescription.Start.Date. - previous.end.date) <= 7, previous.end.date, Prescription.Start.Date.)) %>%
  mutate(discontinuation = if_else(is.na(previous.end.date)|as.numeric(Prescription.Start.Date. - previous.end.date) <= 1, 0, 1)) %>%
  mutate(dis.group = cumsum(discontinuation)) %>%
  ungroup() %>%
  mutate(diff = Prescription.End.Date. - Prescription.Start.Date.) %>%
  group_by(Reference.Key., dis.group) %>%
  mutate(Prescription.Start.Date. = Prescription.End.Date. - cumsum(as.numeric(diff))) %>%
  group_by(Reference.Key., Prescription.Start.Date.) %>%
  filter(Prescription.End.Date. == max(Prescription.End.Date.)) %>%
  ungroup() %>% 
  mutate(duration = as.numeric(Prescription.End.Date. - Prescription.Start.Date. + 1))
entresto_rx$duration %>% summary() # 140[50-280]
Entresto %>% filter(Reference.Key. %in% filter(baseline, Index.Drug == 1)$Reference.Key.) %>% 
  as_tibble() %>% arrange(Reference.Key., Prescription.Start.Date.) %>% group_by(Reference.Key.) %>% slice(1) %>% 
  ungroup() %>% mutate(index.date = Prescription.Start.Date.) %>% arrange(index.date) %>% select(Dosage., Drug.Strength., Drug.Frequency.) %>% 
  mutate(Drug.Frequency. = if_else(Drug.Frequency. %in% c("TWICE DAILY", "THREE TIMES DAILY"), "TWICE DAILY", "DAILY")) %>%
  mutate(Dosage. = if_else(is.na(Dosage.), "1", Dosage.)) %>%
  mutate(Dosage. = str_replace(Dosage., " TAB.+", "")) %>% 
  mutate(Dosage. = as.numeric(Dosage.)) %>% 
  mutate(Drug.Strength. = str_replace(Drug.Strength., "MG", "")) %>% 
  mutate(Drug.Strength. = as.numeric(Drug.Strength.)) %>% 
  mutate(dosage = Dosage. * Drug.Strength.) %>% 
  filter(Drug.Frequency. == "DAILY") %>%
  # filter(Drug.Frequency. == "TWICE DAILY") %>%
  .$dosage %>% table()
  
Entresto$Dosage. %>% table()
baseline$Index.Drug %>% table()
