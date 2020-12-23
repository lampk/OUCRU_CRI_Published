# setup -------------------------------------------------------------------

library(tidyverse)
library(readxl)
library(lubridate)
library(Hmisc)
#library(Lmisc)

# load functions
source("functions.R")

data_path <- file.path("V:", "BIOSTATISTICS", "ClinicalProjectsNonRCT", "CRI_Dengue", "Data", "Data")

load(file.path(data_path, "Rdata", "data_analysis.Rdata")) ## dss_clinical, maindts, cri_final, cri0, cri48, vital_final, vital0, vital48
load(file.path(data_path, "Rdata", "data_analysis_final.Rdata")) ## final, final12, final24, final36

# new analysis ------------------------------------------------------------

## split follow-up time into intervals by shock/re-shock time
new_clinical_DSS <- final %>%
  select(StudyNo, ts0, ts1, ts2, ts3, ts4) %>%
  gather(key = "episode", value = "tstart", -1) %>%
  mutate(episode = as.numeric(gsub(pattern = "ts", replacement = "", x = episode)) + 1,
         tstop0 = tstart + 6) %>%
  filter(!is.na(tstart)) %>%
  arrange(StudyNo, episode)

## get the stop of each shock/reshock period
new_clinical_DSS2 <- new_clinical_DSS %>%
  group_by(StudyNo, episode) %>%
  do(cbind(., tstop = get_tstop(.))) %>%
  ungroup()

## combine overlapped episodes
new_clinical_DSS3 <- new_clinical_DSS2 %>%
  group_by(StudyNo) %>%
  mutate(tstop_lag = lag(tstop),
         tstart2 = ifelse(is.na(tstop_lag),
                          ifelse(is.na(lead(tstart)), tstart,
                                 ifelse(tstop > lead(tstart), NA, tstart)),
                          ifelse(tstop_lag > tstart, lag(tstart), tstart))) %>%
  ungroup()

new_clinical_DSS4 <- new_clinical_DSS3 %>%
  filter(!is.na(tstart2)) %>%
  group_by(StudyNo) %>%
  arrange(tstart2) %>%
  transmute(episode = 1:n(),
            tstart = tstart2,
            tstop = tstop) %>%
  ungroup() %>%
  arrange(StudyNo, episode) %>%
  group_by(StudyNo, tstop) %>%
  arrange(tstart) %>%
  slice(1) %>%
  mutate(episode = 1:n()) %>%
  ungroup() %>%
  arrange(StudyNo, episode)

## get the last vital sign measurement (PP)
last_vital <- vital_final %>%
  filter(fu == 0) %>%
  group_by(StudyNo) %>%
  arrange(t) %>%
  summarise(tlast = last(t)) %>%
  ungroup() %>%
  arrange(StudyNo)

new_clinical_DSS5 <- merge(new_clinical_DSS4,
                           last_vital,
                           by = "StudyNo", all.x = TRUE) %>%
  group_by(StudyNo) %>%
  arrange(tstart) %>%
  mutate(episode = 1:n()) %>%
  ungroup() %>%
  arrange(StudyNo, episode)

new_clinical_DSS5w <- merge(
  select(new_clinical_DSS5, -tstop, -tlast) %>% mutate(episode = paste("tstart", episode, sep = "_")) %>% spread(episode, tstart),
  select(new_clinical_DSS5, -tstart, -tlast) %>% mutate(episode = paste("tstop", episode, sep = "_")) %>% spread(episode, tstop),
  by = "StudyNo"
)

### new follow-up time (at risk intervals)

new_fup <- new_clinical_DSS5 %>%
  group_by(StudyNo) %>%
  mutate(tstart_lead = lead(tstart)) %>%
  transmute(
    tstart = tstop,
    tstop  = ifelse(is.na(tstart_lead), tlast, tstart_lead),
    length = tstop - tstart,
    reshock = ifelse(is.na(tstart_lead), 0, 1)
  ) %>%
  ungroup() %>%
  arrange(StudyNo, tstart)

### all reshock happened within the first 48 hours --> just need to truncate to 48h in those without event
new_fup_48 <- new_fup %>%
  filter(tstart < 48) %>%
  mutate(tstop = ifelse(tstop >= 48, 48, tstop),
         length = tstop - tstart) %>%
  group_by(StudyNo) %>%
  arrange(tstart) %>%
  mutate(interval = 1:n()) %>%
  ungroup() %>%
  arrange(StudyNo, tstart)

# landmark analysis (1): first episodes ------------------------------------------

## set prediction horizon at 6 hours
hor <- 6

## create stacked data based on landmark points
new_fup_48_first_LM <- new_fup_48 %>%
  filter(interval == 1) %>%
  group_by(StudyNo, interval) %>%
  do(cbind(., LM = get_LM(., LMs = 0:42))) %>%
  mutate(LM_start = LM,
         LM_stop = ifelse((LM + hor) >= tstop, tstop, LM + hor),
         status  = ifelse(LM_stop < tstop, 0, reshock)) %>%
  ungroup()



## split data at each event time
library(survival)
tt <- sort(unique(new_fup_48_first_LM$LM_stop[new_fup_48_first_LM$status == 1]))
new_fup_48_first_LM2 <- survSplit(data = new_fup_48_first_LM,
                                  cut = tt,
                                  end = "LM_stop", start = "LM_start", event = "status") %>%
  filter(!is.na(LM_start)) %>%
  mutate(d = LM_stop - LM_start)

## get CRI
new_fup_48_first_LM2_CRI <- new_fup_48_first_LM2 %>%
  mutate(idx = 1:n()) %>%
  group_by(idx) %>%
  do(cbind(., get_CRI(id = .$StudyNo, tp = .$LM, cri_dat = cri48, tol = 15/60))) %>%
  ungroup()

new_fup_48_first_LM2_CRI1 <- new_fup_48_first_LM2_CRI %>%
  filter(!is.na(CRI)) %>%
  mutate(CRI2 = CRI * 10,
         CRIt1 = CRI2 * d)

new_fup_48_first_LM2_CRI2 <- new_fup_48_first_LM2_CRI %>%
  filter(!is.na(CRI)) %>%
  mutate(CRI2 = CRI * 10,
         CRIt1 = CRI2 * d,
         CRIt2 = CRI2 * d^2)

new_fup_48_first_LM2_CRI3 <- new_fup_48_first_LM2_CRI %>%
  filter(!is.na(CRI)) %>%
  mutate(CRI2 = CRI * 10,
         CRIt1 = CRI2 * pmin(d, 1.25),
         CRIt2 = CRI2 * pmax(d, 1.25))

### extended landmark model

LMfit1 <- coxph(Surv(LM_start, LM_stop, status) ~ CRI2 + CRIt1 + strata(LM) + cluster(StudyNo), data = new_fup_48_first_LM2_CRI1, method = "breslow")
LMfit1

LMfit2 <- coxph(Surv(LM_start, LM_stop, status) ~ CRI2 + CRIt1 + CRIt2 + strata(LM) + cluster(StudyNo), data = new_fup_48_first_LM2_CRI2, method = "breslow")
LMfit2

LMfit3 <- coxph(Surv(LM_start, LM_stop, status) ~ CRI2 + CRIt1 + CRIt2 + strata(LM) + cluster(StudyNo), data = new_fup_48_first_LM2_CRI3, method = "breslow")
LMfit3

### plot

d <- seq(0, 6, by = 0.1)
y <- coef(LMfit1)[["CRI2"]] +
  coef(LMfit1)[["CRIt1"]] * d
plot(d, y, type = "l", lwd=2, xlab = "t - s", ylab = "Log hazard ratio")

d <- seq(0, 6, by = 0.1)
y <- coef(LMfit2)[["CRI2"]] +
  coef(LMfit2)[["CRIt1"]] * d +
  coef(LMfit2)[["CRIt2"]] * d^2
plot(d, y, type = "l", lwd=2, xlab = "t - s", ylab = "Log hazard ratio")

d <- seq(0, 6, by = 0.1)
y <- coef(LMfit3)[["CRI2"]] +
  coef(LMfit3)[["CRIt1"]] * pmin(d, 1.25) +
  coef(LMfit3)[["CRIt2"]] * pmax(d, 1.25)
plot(d, y, type = "l", lwd=2, xlab = "t - s", ylab = "Log hazard ratio")

## calculate CI

mean_current <- coef(LMfit3)
names(mean_current) <- c("x1", "x2", "x3")

vcov_current <- vcov(LMfit3)
names(vcov_current) <- c("x1", "x2", "x3")

library(msm)

plotdat <- data.frame(d = seq(0, 6, by = 0.1)) %>%
  mutate(form = sprintf("~ x1 + %f * x2 + %f * x3", pmin(d, 1.25), pmax(d, 1.25)),
         lHR =  coef(LMfit3)[["CRI2"]] + coef(LMfit3)[["CRIt1"]] * pmin(d, 1.25) + coef(LMfit3)[["CRIt2"]] * pmax(d, 1.25),
         se = sapply(form, function(x) deltamethod(g = as.formula(x), mean = mean_current, cov = vcov_current)),
         lHR_lo = lHR - 1.96 * se,
         lHR_hi = lHR + 1.96 * se)

ggplot(data = plotdat, aes(x = d)) +
  geom_point(aes(y = lHR)) +
  geom_errorbar(aes(ymin = lHR_lo, ymax = lHR_hi)) +
  geom_hline(yintercept = 0, linetype = 2) +
  theme_bw()

ggplot(data = subset(plotdat, d <= 2), aes(x = d)) +
  geom_point(aes(y = lHR)) +
  geom_errorbar(aes(ymin = lHR_lo, ymax = lHR_hi)) +
  geom_hline(yintercept = 0, linetype = 2) +
  theme_bw()
