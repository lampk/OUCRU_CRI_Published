
load(file.path("Data", "Rdata", "data_analysis.Rdata")) ## dss_clinical, maindts, cri_final, cri0, cri48, vital_final, vital0, vital48
load(file.path("Data", "Rdata", "data_analysis_final.Rdata")) ## final, final12, final24, final36


# old analysis ------------------------------------------------------------

## from CRI to PP
crilength <- cri48 %>%
  group_by(StudyNo, group) %>%
  arrange(t) %>%
  summarise(start = first(t),
            end = last(t),
            duration = end - start) %>%
  ungroup() %>%
  mutate(id_group = paste(StudyNo, group, sep = "_"))

crilength_15 <- crilength$id_group[crilength$duration >= (15/60)]

cridata_15 <- cri48 %>%
  mutate(id_group = paste(StudyNo, group, sep = "_")) %>%
  filter(id_group %in% crilength_15) %>%
  arrange(StudyNo, t)

library(geepack)

mygee_fit <- function(formula, data) {
  fitdat <- data[!is.na(data[, all.vars(formula)[2]]),] %>%
    mutate(StudyNo = as.factor(StudyNo))
  fit <- geeglm(formula = formula, family = binomial, corstr = "independence",
                id = StudyNo, data = fitdat)
  return(fit)
}

mygee_fit2 <- function(formula, data) {
  fitdat <- data[!is.na(data[, all.vars(formula)[2]]),] %>%
    mutate(StudyNo = as.factor(StudyNo))
  fit <- geeglm(formula = formula, family = gaussian, corstr = "independence",
                id = StudyNo, data = fitdat)
  return(fit)
}

mysummary <- function(fit) {
  est <- coef(fit)[[2]]
  se  <- summary(fit)$coef[, 2][2]
  pval <- summary(fit)$coef[2, 4]
  or <- 1/exp(est)
  or_lo <- 1/exp(est + 1.96 * se)
  or_hi <- 1/exp(est - 1.96 * se)

  ## return
  output <- t(c(or = or, or_lo = or_lo, or_hi = or_hi, p = pval))
  rownames(output) <- all.vars(formula(fit))[2]
  return(output)
}

mysummary2 <- function(fit) {
  est <- coef(fit)[[2]]
  se  <- summary(fit)$coef[, 2][2]
  pval <- summary(fit)$coef[2, 4]
  est_lo <- est - 1.96 * se
  est_hi <- est + 1.96 * se

  ## return
  output <- t(c(est = est, est_lo = est_lo, est_hi = est_hi, p = pval))
  rownames(output) <- all.vars(formula(fit))[2]
  return(output)
}

get_closest_pp <- function(id, tp, pp_dat = vital, tol = 30/60) {
  tmp <- subset(pp_dat, StudyNo == id & (t >= tp) & (t <= (tp + tol))) %>%
    arrange(t)
  if (nrow(tmp) == 0) {
    output <- NA
  } else {
    output <- cbind(PP = first(tmp$PP),
                    PP_t = last(tmp$t))
  }
  return(output)
}

tmp <- subset(vital48, !is.na(PP) & StudyNo %in% final$StudyNo) %>% arrange(StudyNo, t)

PP.m025 <- do.call(rbind, lapply(1:nrow(cridata_15), function(i) {get_closest_pp(id = cridata_15$StudyNo[i], tp = (cridata_15$t[i] + 0.25), pp_dat = tmp, tol = 30/60)}))
PP.m050 <- do.call(rbind, lapply(1:nrow(cridata_15), function(i) {get_closest_pp(id = cridata_15$StudyNo[i], tp = (cridata_15$t[i] + 0.5), pp_dat = tmp, tol = 30/60)}))
PP.m075 <- do.call(rbind, lapply(1:nrow(cridata_15), function(i) {get_closest_pp(id = cridata_15$StudyNo[i], tp = (cridata_15$t[i] + 0.75), pp_dat = tmp, tol = 30/60)}))
PP.m100 <- do.call(rbind, lapply(1:nrow(cridata_15), function(i) {get_closest_pp(id = cridata_15$StudyNo[i], tp = (cridata_15$t[i] + 1), pp_dat = tmp, tol = 30/60)}))
PP.m200 <- do.call(rbind, lapply(1:nrow(cridata_15), function(i) {get_closest_pp(id = cridata_15$StudyNo[i], tp = (cridata_15$t[i] + 2), pp_dat = tmp, tol = 30/60)}))
PP.m300 <- do.call(rbind, lapply(1:nrow(cridata_15), function(i) {get_closest_pp(id = cridata_15$StudyNo[i], tp = (cridata_15$t[i] + 3), pp_dat = tmp, tol = 30/60)}))
PP.m400 <- do.call(rbind, lapply(1:nrow(cridata_15), function(i) {get_closest_pp(id = cridata_15$StudyNo[i], tp = (cridata_15$t[i] + 4), pp_dat = tmp, tol = 30/60)}))
PP.m500 <- do.call(rbind, lapply(1:nrow(cridata_15), function(i) {get_closest_pp(id = cridata_15$StudyNo[i], tp = (cridata_15$t[i] + 5), pp_dat = tmp, tol = 30/60)}))
PP.m600 <- do.call(rbind, lapply(1:nrow(cridata_15), function(i) {get_closest_pp(id = cridata_15$StudyNo[i], tp = (cridata_15$t[i] + 6), pp_dat = tmp, tol = 30/60)}))
PP.m700 <- do.call(rbind, lapply(1:nrow(cridata_15), function(i) {get_closest_pp(id = cridata_15$StudyNo[i], tp = (cridata_15$t[i] + 7), pp_dat = tmp, tol = 30/60)}))

dat <- data.frame(cbind(select(cridata_15, StudyNo, age, sex, wt, t, CRI),
                        PP.m025 = PP.m025[,"PP"],
                        PP.m050 = PP.m050[,"PP"],
                        PP.m075 = PP.m075[,"PP"],
                        PP.m100 = PP.m100[,"PP"],
                        PP.m200 = PP.m200[,"PP"],
                        PP.m300 = PP.m300[,"PP"],
                        PP.m400 = PP.m400[,"PP"],
                        PP.m500 = PP.m500[,"PP"],
                        PP.m600 = PP.m600[,"PP"],
                        PP.m700 = PP.m700[,"PP"]))

cripp_dat2 <- gather(dat, key = "tlag", value = "PP", -1, -2, -3, -4, -5, -6) %>%
  mutate(pp.low = as.numeric(I(PP <= 20))) %>%
  arrange(StudyNo, t) %>%
  filter(!is.na(PP))

fit.m025 <- mygee_fit(formula = pp.low ~ I(CRI * 10), data = subset(cripp_dat2, tlag == "PP.m025"))
fit.m050 <- mygee_fit(formula = pp.low ~ I(CRI * 10), data = subset(cripp_dat2, tlag == "PP.m050"))
fit.m075 <- mygee_fit(formula = pp.low ~ I(CRI * 10), data = subset(cripp_dat2, tlag == "PP.m075"))
fit.m100 <- mygee_fit(formula = pp.low ~ I(CRI * 10), data = subset(cripp_dat2, tlag == "PP.m100"))
fit.m200 <- mygee_fit(formula = pp.low ~ I(CRI * 10), data = subset(cripp_dat2, tlag == "PP.m200"))
fit.m300 <- mygee_fit(formula = pp.low ~ I(CRI * 10), data = subset(cripp_dat2, tlag == "PP.m300"))
fit.m400 <- mygee_fit(formula = pp.low ~ I(CRI * 10), data = subset(cripp_dat2, tlag == "PP.m400"))
fit.m500 <- mygee_fit(formula = pp.low ~ I(CRI * 10), data = subset(cripp_dat2, tlag == "PP.m500"))
fit.m600 <- mygee_fit(formula = pp.low ~ I(CRI * 10), data = subset(cripp_dat2, tlag == "PP.m600"))
fit.m700 <- mygee_fit(formula = pp.low ~ I(CRI * 10), data = subset(cripp_dat2, tlag == "PP.m700"))

cripp_tab <- do.call(rbind, lapply(list(fit.m025, fit.m050, fit.m075, fit.m100, fit.m200, fit.m300, fit.m400, fit.m500, fit.m600, fit.m700), function(x) mysummary(x)))
round(cripp_tab, 2)

#---- after removing

cridata_152 <- merge(cri48,
                     new_clinical_DSS5w,
                     by = "StudyNo", all.x = TRUE) %>%
  mutate(id_group = paste(StudyNo, group, sep = "_")) %>%
  filter(id_group %in% crilength_15) %>%
  arrange(StudyNo, t) %>%
  filter(!is.na(tstart_1) & (t > tstop_1) &
           (is.na(tstart_2) |
              (!is.na(tstart_2) & ((t < tstart_2) | ((t > tstop_2))) &
                 (is.na(tstart_3) |
                    (!is.na(tstart_3) & ((t < tstart_3) | (t > tstop_3)) &
                       (is.na(tstart_4) |
                          (!is.na(tstart_4) & ((t < tstart_4) | (t > tstop_4))))))))) %>%
  arrange(StudyNo, t)

tmp <- subset(vital48, !is.na(PP) & StudyNo %in% final$StudyNo) %>% arrange(StudyNo, t)

PP.m0252 <- do.call(rbind, lapply(1:nrow(cridata_152), function(i) {get_closest_pp(id = cridata_152$StudyNo[i], tp = (cridata_152$t[i] + 0.25), pp_dat = tmp, tol = 30/60)}))
PP.m0502 <- do.call(rbind, lapply(1:nrow(cridata_152), function(i) {get_closest_pp(id = cridata_152$StudyNo[i], tp = (cridata_152$t[i] + 0.5), pp_dat = tmp, tol = 30/60)}))
PP.m0752 <- do.call(rbind, lapply(1:nrow(cridata_152), function(i) {get_closest_pp(id = cridata_152$StudyNo[i], tp = (cridata_152$t[i] + 0.75), pp_dat = tmp, tol = 30/60)}))
PP.m1002 <- do.call(rbind, lapply(1:nrow(cridata_152), function(i) {get_closest_pp(id = cridata_152$StudyNo[i], tp = (cridata_152$t[i] + 1), pp_dat = tmp, tol = 30/60)}))
PP.m2002 <- do.call(rbind, lapply(1:nrow(cridata_152), function(i) {get_closest_pp(id = cridata_152$StudyNo[i], tp = (cridata_152$t[i] + 2), pp_dat = tmp, tol = 30/60)}))
PP.m3002 <- do.call(rbind, lapply(1:nrow(cridata_152), function(i) {get_closest_pp(id = cridata_152$StudyNo[i], tp = (cridata_152$t[i] + 3), pp_dat = tmp, tol = 30/60)}))
PP.m4002 <- do.call(rbind, lapply(1:nrow(cridata_152), function(i) {get_closest_pp(id = cridata_152$StudyNo[i], tp = (cridata_152$t[i] + 4), pp_dat = tmp, tol = 30/60)}))
PP.m5002 <- do.call(rbind, lapply(1:nrow(cridata_152), function(i) {get_closest_pp(id = cridata_152$StudyNo[i], tp = (cridata_152$t[i] + 5), pp_dat = tmp, tol = 30/60)}))
PP.m6002 <- do.call(rbind, lapply(1:nrow(cridata_152), function(i) {get_closest_pp(id = cridata_152$StudyNo[i], tp = (cridata_152$t[i] + 6), pp_dat = tmp, tol = 30/60)}))
PP.m7002 <- do.call(rbind, lapply(1:nrow(cridata_152), function(i) {get_closest_pp(id = cridata_152$StudyNo[i], tp = (cridata_152$t[i] + 7), pp_dat = tmp, tol = 30/60)}))

dat2 <- data.frame(cbind(select(cridata_152, StudyNo, age, sex, wt, t, CRI),
                         PP.m025 = PP.m0252[,"PP"],
                         PP.m050 = PP.m0502[,"PP"],
                         PP.m075 = PP.m0752[,"PP"],
                         PP.m100 = PP.m1002[,"PP"],
                         PP.m200 = PP.m2002[,"PP"],
                         PP.m300 = PP.m3002[,"PP"],
                         PP.m400 = PP.m4002[,"PP"],
                         PP.m500 = PP.m5002[,"PP"],
                         PP.m600 = PP.m6002[,"PP"],
                         PP.m700 = PP.m7002[,"PP"]))

cripp_dat22 <- gather(dat2, key = "tlag", value = "PP", -1, -2, -3, -4, -5, -6) %>%
  mutate(pp.low = as.numeric(I(PP <= 20))) %>%
  arrange(StudyNo, t) %>%
  filter(!is.na(PP))

fit.m025 <- mygee_fit(formula = pp.low ~ I(CRI * 10), data = subset(cripp_dat22, tlag == "PP.m025"))
fit.m050 <- mygee_fit(formula = pp.low ~ I(CRI * 10), data = subset(cripp_dat22, tlag == "PP.m050"))
fit.m075 <- mygee_fit(formula = pp.low ~ I(CRI * 10), data = subset(cripp_dat22, tlag == "PP.m075"))
fit.m100 <- mygee_fit(formula = pp.low ~ I(CRI * 10), data = subset(cripp_dat22, tlag == "PP.m100"))
fit.m200 <- mygee_fit(formula = pp.low ~ I(CRI * 10), data = subset(cripp_dat22, tlag == "PP.m200"))
fit.m300 <- mygee_fit(formula = pp.low ~ I(CRI * 10), data = subset(cripp_dat22, tlag == "PP.m300"))
fit.m400 <- mygee_fit(formula = pp.low ~ I(CRI * 10), data = subset(cripp_dat22, tlag == "PP.m400"))
fit.m500 <- mygee_fit(formula = pp.low ~ I(CRI * 10), data = subset(cripp_dat22, tlag == "PP.m500"))
fit.m600 <- mygee_fit(formula = pp.low ~ I(CRI * 10), data = subset(cripp_dat22, tlag == "PP.m600"))
fit.m700 <- mygee_fit(formula = pp.low ~ I(CRI * 10), data = subset(cripp_dat22, tlag == "PP.m700"))

cripp_tab2 <- do.call(rbind, lapply(list(fit.m025, fit.m050, fit.m075, fit.m100, fit.m200, fit.m300, fit.m400, fit.m500, fit.m600, fit.m700), function(x) mysummary(x)))
round(cripp_tab2, 2)

# new analysis ------------------------------------------------------------


## new type of data

new_clinical_DSS <- final %>%
  select(StudyNo, ts0, ts1, ts2, ts3, ts4) %>%
  gather(key = "episode", value = "tstart", -1) %>%
  mutate(episode = as.numeric(gsub(pattern = "ts", replacement = "", x = episode)) + 1,
         tstop0 = tstart + 6) %>%
  filter(!is.na(tstart)) %>%
  arrange(StudyNo, episode)

tmp_func <- function(x){
  tmp <- filter(vital_final, (StudyNo == x$StudyNo) & (t >= x$tstart) & (t <= x$tstop0))
  output <- sum(tmp$PP <= 20)
  return(output)
}

get_tstop <- function(x){
  tmp1 <- filter(vital_final, (StudyNo == x$StudyNo) & (t >= x$tstart) & !is.na(PP)) %>%
    arrange(t)

  tmp2 <- sapply(1:nrow(tmp1), function(i){
    tmp_20 <- sum(tmp1$PP[tmp1$t >= (tmp1$t[i] - 6) & (tmp1$t <= tmp1$t[i])] <= 20)
    return(ifelse(tmp_20 > 0, 1, 0))
  })

  if (any(tmp2 == 0)) {
    output <- min(tmp1$t[tmp2 == 0 & tmp1$t >= (x$tstart + 6)])
  } else {
    output <- NA
  }
  return(output)
}

get_tstop2 <- function(x){
  tmp1 <- filter(vital_final, (StudyNo == x$StudyNo) & (t >= x$tstart) & !is.na(PP)) %>%
    arrange(t)

  tmp2 <- sapply(1:nrow(tmp1), function(i){
    tmp_20 <- sum(tmp1$PP[tmp1$t >= (tmp1$t[i] - 6) & (tmp1$t <= tmp1$t[i])] <= 20)
    return(ifelse(tmp_20 > 0, 1, 0))
  })

  if (any(tmp2 == 0)) {
    output <- min(tmp1$t[tmp2 == 0 & tmp1$t >= (x$tstart + 6)])
  } else {
    output <- NA
  }
  return(output)
}

#### get tstop
new_clinical_DSS2 <- new_clinical_DSS %>%
  group_by(StudyNo, episode) %>%
  do(cbind(., tstop = get_tstop(.))) %>%
  ungroup()

#### combine overlapped episode
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

### new follow-up time

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

# split participants by reshock status --------------------------------------

new_clinical_DSS6 <- new_clinical_DSS5 %>%
  mutate(duration = tstop - tstart)

final$reshock_type <- with(final,
                           ifelse(is.na(ts1), "no reshock",
                                  ifelse()))

# landmark analysis -------------------------------------------------------

hor <- 6

get_LM <- function(x, LMs) {
  LM <- LMs[LMs >= x$tstart & LMs <= x$tstop]
  if (length(LM) == 0) {
    return(NA)
  } else {
    return(LM)
  }
}

get_CRI <- function(id, tp, cri_dat = cridata, tol = 15/60, type = c("current", "slope")) {
  if (length(type) > 1) type <- "current"
  tmp <- subset(cri_dat, StudyNo == id & (t <= tp) & (t >= (tp - tol))) %>%
    arrange(t)
  if (nrow(tmp) == 0) {
    output <- cbind(CRI = NA, CRI_t = NA)
  } else {
    if (type == "current") {
      output <- cbind(CRI = last(tmp$CRI),
                      CRI_t = last(tmp$t))
    } else {
      if (type == "slope") {
        if (nrow(tmp) < 3) {
          output <- cbind(CRI = NA, CRI_t = NA)
        } else {
          output <- cbind(CRI = coef(lm(tmp$CRI ~ tmp$t))[2],
                          CRI_t = last(tmp$t))
        }
      }
    }

  }
  return(output)
}

get_CRI2 <- function(id, tp, cri_dat = cridata, tol = 15/60, type = c("current", "slope")) {
  if (length(type) > 1) type <- "current"
  tmp <- subset(cri_dat, StudyNo == id & (t <= tp) & (t >= (tp - tol))) %>%
    arrange(t)
  if (nrow(tmp) == 0) {
    output <- NA
  } else {
    if (type == "current") {
      output <- last(tmp$CRI)
    } else {
      if (type == "slope") {
        if (nrow(tmp) < 3) {
          output <- NA
        } else {
          output <- coef(lm(tmp$CRI ~ tmp$t))[2]
        }
      }
    }

  }
  return(output)
}

## all episodes

new_fup_48_LM <- new_fup_48 %>%
  group_by(StudyNo, interval) %>%
  do(cbind(., LM = get_LM(., LMs = 0:42))) %>%
  mutate(LM_start = LM,
         LM_stop = ifelse((LM + hor) >= tstop, tstop, LM + hor),
         status  = ifelse(LM_stop < tstop, 0, reshock)) %>%
  ungroup()

library(survival)
tt <- sort(unique(new_fup_48_LM$LM_stop[new_fup_48_LM$status == 1]))
new_fup_48_LM2 <- survSplit(data = new_fup_48_LM,
                            cut = tt,
                            end = "LM_stop", start = "LM_start", event = "status") %>%
  filter(!is.na(LM_start))

new_fup_48_LM_CRI <- new_fup_48_LM %>%
  mutate(idx = 1:n()) %>%
  group_by(idx) %>%
  do(cbind(., get_CRI(id = .$StudyNo, tp = .$LM, cri_dat = cridata, tol = 15/60))) %>%
  ungroup()

new_fup_48_LM2_CRI <- new_fup_48_LM2 %>%
  mutate(idx = 1:n()) %>%
  group_by(idx) %>%
  do(cbind(., get_CRI(id = .$StudyNo, tp = .$LM, cri_dat = cridata, tol = 15/60))) %>%
  ungroup()

new_fup_48_LM2_CRI$CRIt <- new_fup_48_LM2_CRI$CRI * (new_fup_48_LM2_CRI$LM_stop - new_fup_48_LM2_CRI$LM_start)
new_fup_48_LM2_CRI$CRIt2 <- new_fup_48_LM2_CRI$CRI * (new_fup_48_LM2_CRI$LM_stop - new_fup_48_LM2_CRI$LM_start)^2

### simple landmark model

LMfit_01 <- coxph(Surv(LM_start, LM_stop, status) ~ CRI + strata(LM) + cluster(StudyNo), data = new_fup_48_LM_CRI, method = "breslow")
LMfit_01

### extended landmark model

LMfit_02 <- coxph(Surv(LM_start, LM_stop, status) ~ CRI + CRIt + CRIt2 + strata(LM) + cluster(StudyNo), data = new_fup_48_LM2_CRI, method = "breslow")
LMfit_02

### plot

tseq <- seq(0, 6, by = 0.1)
plot(tseq, coef(LMfit_02)[["CRI"]] + coef(LMfit_02)[["CRIt"]] * tseq + coef(LMfit_02)[["CRIt2"]] * tseq^2,
     type = "l", lwd=2, xlab = "t - s", ylab = "Log hazard ratio")
lines(c(0, 6),rep(coef(LMfit_01)[["CRI"]], 2), type = "l", lty = 3)

## first episodes

hor <- 6

new_fup_48_first_LM <- new_fup_48 %>%
  filter(interval == 1) %>%
  group_by(StudyNo, interval) %>%
  do(cbind(., LM = get_LM(., LMs = 0:42))) %>%
  mutate(LM_start = LM,
         LM_stop = ifelse((LM + hor) >= tstop, tstop, LM + hor),
         status  = ifelse(LM_stop < tstop, 0, reshock)) %>%
  ungroup()

tt <- sort(unique(new_fup_48_first_LM$LM_stop[new_fup_48_first_LM$status == 1]))
new_fup_48_first_LM2 <- survSplit(data = new_fup_48_first_LM,
                                  cut = tt,
                                  end = "LM_stop", start = "LM_start", event = "status") %>%
  filter(!is.na(LM_start))

new_fup_48_first_LM_CRI <- new_fup_48_first_LM %>%
  mutate(idx = 1:n()) %>%
  group_by(idx) %>%
  do(cbind(., get_CRI(id = .$StudyNo, tp = .$LM, cri_dat = cridata, tol = 15/60))) %>%
  ungroup()

new_fup_48_first_LM2_CRI <- new_fup_48_first_LM2 %>%
  mutate(idx = 1:n()) %>%
  group_by(idx) %>%
  do(cbind(., get_CRI(id = .$StudyNo, tp = .$LM, cri_dat = cridata, tol = 15/60))) %>%
  ungroup()

new_fup_48_first_LM_CRI$CRI2 <- new_fup_48_first_LM_CRI$CRI * 10
new_fup_48_first_LM2_CRI$CRI2 <- new_fup_48_first_LM2_CRI$CRI * 10

### simple landmark model

LMfit_first_01 <- coxph(Surv(LM_start, LM_stop, status) ~ CRI2 + strata(LM) + cluster(StudyNo), data = new_fup_48_first_LM_CRI, method = "breslow")
LMfit_first_01

### extended landmark model

new_fup_48_first_LM2_CRI$CRIt <- new_fup_48_first_LM2_CRI$CRI2 * (new_fup_48_first_LM2_CRI$LM_stop - new_fup_48_first_LM2_CRI$LM_start)
new_fup_48_first_LM2_CRI$CRIt2 <- new_fup_48_first_LM2_CRI$CRI2 * (new_fup_48_first_LM2_CRI$LM_stop - new_fup_48_first_LM2_CRI$LM_start)^2

LMfit_first_02 <- coxph(Surv(LM_start, LM_stop, status) ~ CRI2 + CRIt + CRIt2 + strata(LM) + cluster(StudyNo), data = new_fup_48_first_LM2_CRI, method = "breslow")
LMfit_first_02


### plot

library(msm)
meanx <- coef(LMfit_first_02); names(meanx) <- c("x1", "x2", "x3")
covx <- survival:::vcov.coxph(LMfit_first_02); names(covx) <- c("x1", "x2", "x3")
tseq <- seq(0, 6, by = 0.1)
form <- sprintf("~ x1 + %f * x2 + %f * x3", tseq, tseq^2)
se <- sapply(form, function(x) deltamethod(g = as.formula(x), mean = meanx, cov = covx))

plotdat <- data.frame(t = tseq) %>%
  mutate(form = form,
         lHR = coef(LMfit_first_02)[["CRI2"]] + coef(LMfit_first_02)[["CRIt"]] * t + coef(LMfit_first_02)[["CRIt2"]] * t^2,
         lHR_se = se,
         lHR_lo = lHR - 1.96 * lHR_se,
         lHR_hi = lHR + 1.96 * lHR_se)

ggplot(data = plotdat, aes(x = t)) +
  geom_point(aes(y = lHR)) +
  geom_errorbar(aes(ymin = lHR_lo, ymax = lHR_hi)) +
  theme_bw()

ggplot(data = subset(plotdat, t <= 2), aes(x = t)) +
  geom_point(aes(y = lHR)) +
  geom_errorbar(aes(ymin = lHR_lo, ymax = lHR_hi)) +
  theme_bw()

## logistic regression

new_fup_48_first_LM_CRI2 <- new_fup_48_first_LM %>%
  mutate(idx = 1:n()) %>%
  group_by(idx) %>%
  do(cbind(.,
           CRI = 10 * get_CRI2(id = .$StudyNo, tp = .$LM, cri_dat = cridata, tol = 15/60, type = "current"),
           sCRI_15 = 10 * get_CRI2(id = .$StudyNo, tp = .$LM, cri_dat = cridata, tol = 15/60, type = "slope"),
           sCRI_30 = 10 * get_CRI2(id = .$StudyNo, tp = .$LM, cri_dat = cridata, tol = 30/60, type = "slope"),
           sCRI_60 = 10 * get_CRI2(id = .$StudyNo, tp = .$LM, cri_dat = cridata, tol = 60/60, type = "slope"))) %>%
  ungroup() %>%
  mutate(StudyNo = as.factor(StudyNo),
         dur = tstop - LM_start)

new_fup_48_first_LM_CRI2$CRIt <- new_fup_48_first_LM_CRI2$CRI * (new_fup_48_first_LM_CRI2$dur)
new_fup_48_first_LM_CRI2$CRIt2 <- new_fup_48_first_LM_CRI2$CRI * (new_fup_48_first_LM_CRI2$dur)^2
new_fup_48_first_LM_CRI2$sCRI_15t <- new_fup_48_first_LM_CRI2$sCRI_15 * (new_fup_48_first_LM_CRI2$dur)
new_fup_48_first_LM_CRI2$sCRI_15t2 <- new_fup_48_first_LM_CRI2$sCRI_15 * (new_fup_48_first_LM_CRI2$dur)^2
new_fup_48_first_LM_CRI2$sCRI_30t <- new_fup_48_first_LM_CRI2$sCRI_30 * (new_fup_48_first_LM_CRI2$dur)
new_fup_48_first_LM_CRI2$sCRI_30t2 <- new_fup_48_first_LM_CRI2$sCRI_30 * (new_fup_48_first_LM_CRI2$dur)^2

### fit models
library(geepack)

LMfit_glm_first_current <- geeglm(formula = I(reshock == 1) ~ CRI + CRIt + CRIt2,
                                  family = binomial, corstr = "independence", id = StudyNo, data = subset(new_fup_48_first_LM_CRI2, !is.na(CRI)))
LMfit_glm_first_slope15 <- geeglm(formula = I(reshock == 1) ~ sCRI_15 + sCRI_15t + sCRI_15t2,
                                  family = binomial, corstr = "independence", id = StudyNo, data = subset(new_fup_48_first_LM_CRI2, !is.na(sCRI_15)))
LMfit_glm_first_slope30 <- geeglm(formula = I(reshock == 1) ~ sCRI_30 + sCRI_30t + sCRI_30t2,
                                  family = binomial, corstr = "independence", id = StudyNo, data = subset(new_fup_48_first_LM_CRI2, !is.na(sCRI_30)))

### plot
mean_current <- coef(LMfit_glm_first_current)[-1]
mean_slope15 <- coef(LMfit_glm_first_slope15)[-1]
mean_slope30 <- coef(LMfit_glm_first_slope30)[-1]
names(mean_current) <- names(mean_slope15) <- names(mean_slope30) <- c("x1", "x2", "x3")

vcov_current <- vcov(LMfit_glm_first_current)[-1, -1]
vcov_slope15 <- vcov(LMfit_glm_first_slope15)[-1, -1]
vcov_slope30 <- vcov(LMfit_glm_first_slope30)[-1, -1]
names(vcov_current) <- names(vcov_slope15) <- names(vcov_slope30) <- c("x1", "x2", "x3")

plotdat <- expand.grid(t = seq(0, 48, by = 1),
                       type = c("current", "slope 15", "slope 30")) %>%
  mutate(form = sprintf("~ x1 + %f * x2 + %f * x3", t, t^2),
         lOR = ifelse(type == "current", coef(LMfit_glm_first_current)[["CRI"]] + coef(LMfit_glm_first_current)[["CRIt"]] * t +
                        coef(LMfit_glm_first_current)[["CRIt2"]] * t^2,
                      ifelse(type == "slope 15", coef(LMfit_glm_first_slope15)[["sCRI_15"]] + coef(LMfit_glm_first_slope15)[["sCRI_15t"]] * t +
                               coef(LMfit_glm_first_slope15)[["sCRI_15t2"]] * t^2,
                             coef(LMfit_glm_first_slope30)[["sCRI_30"]] + coef(LMfit_glm_first_slope30)[["sCRI_30t"]] * t +
                               coef(LMfit_glm_first_slope30)[["sCRI_30t2"]] * t^2)),
         se = ifelse(type == "current", sapply(form, function(x) deltamethod(g = as.formula(x), mean = mean_current, cov = vcov_current)),
                     ifelse(type == "slope 15", sapply(form, function(x) deltamethod(g = as.formula(x), mean = mean_slope15, cov = vcov_slope15)),
                            sapply(form, function(x) deltamethod(g = as.formula(x), mean = mean_slope30, cov = vcov_slope30)))),
         lOR_lo = lOR - 1.96 * se,
         lOR_hi = lOR + 1.96 * se)

ggplot(data = subset(plotdat, t <= 12), aes(x = t)) +
  geom_point(aes(y = lOR)) +
  geom_errorbar(aes(ymin = lOR_lo, ymax = lOR_hi)) +
  facet_wrap(type ~ .) +
  theme_bw()

## LM

new_fup_48_first_LM2_CRI2 <- new_fup_48_first_LM2 %>%
  mutate(idx = 1:n()) %>%
  group_by(idx) %>%
  do(cbind(.,
           CRI = 10 * get_CRI2(id = .$StudyNo, tp = .$LM, cri_dat = cridata, tol = 15/60, type = "current"),
           sCRI_15 = 10 * get_CRI2(id = .$StudyNo, tp = .$LM, cri_dat = cridata, tol = 15/60, type = "slope"),
           sCRI_30 = 10 * get_CRI2(id = .$StudyNo, tp = .$LM, cri_dat = cridata, tol = 30/60, type = "slope"),
           sCRI_60 = 10 * get_CRI2(id = .$StudyNo, tp = .$LM, cri_dat = cridata, tol = 60/60, type = "slope"))) %>%
  ungroup() %>%
  mutate(StudyNo = as.factor(StudyNo),
         dur = LM_stop - LM_start)

new_fup_48_first_LM2_CRI2$CRIt <- new_fup_48_first_LM2_CRI2$CRI * (new_fup_48_first_LM2_CRI2$dur)
new_fup_48_first_LM2_CRI2$CRIt2 <- new_fup_48_first_LM2_CRI2$CRI * (new_fup_48_first_LM2_CRI2$dur)^2
new_fup_48_first_LM2_CRI2$sCRI_15t <- new_fup_48_first_LM2_CRI2$sCRI_15 * (new_fup_48_first_LM2_CRI2$dur)
new_fup_48_first_LM2_CRI2$sCRI_15t2 <- new_fup_48_first_LM2_CRI2$sCRI_15 * (new_fup_48_first_LM2_CRI2$dur)^2
new_fup_48_first_LM2_CRI2$sCRI_30t <- new_fup_48_first_LM2_CRI2$sCRI_30 * (new_fup_48_first_LM2_CRI2$dur)
new_fup_48_first_LM2_CRI2$sCRI_30t2 <- new_fup_48_first_LM2_CRI2$sCRI_30 * (new_fup_48_first_LM2_CRI2$dur)^2

LMfit_cox_first_current <- coxph(Surv(LM_start, LM_stop, status) ~ CRI + CRIt + CRIt2 + strata(LM) + cluster(StudyNo), data = new_fup_48_first_LM2_CRI2, method = "breslow")
LMfit_cox_first_current
LMfit_cox_first_slope15 <- coxph(Surv(LM_start, LM_stop, status) ~ sCRI_15 + sCRI_15t + sCRI_15t2 + strata(LM) + cluster(StudyNo), data = new_fup_48_first_LM2_CRI2, method = "breslow")
LMfit_cox_first_slope15
LMfit_cox_first_slope30 <- coxph(Surv(LM_start, LM_stop, status) ~ sCRI_30 + sCRI_30t + sCRI_30t2 + strata(LM) + cluster(StudyNo), data = new_fup_48_first_LM2_CRI2, method = "breslow")
LMfit_cox_first_slope30

mean_current <- coef(LMfit_cox_first_current)
mean_slope15 <- coef(LMfit_cox_first_slope15)
mean_slope30 <- coef(LMfit_cox_first_slope30)
names(mean_current) <- names(mean_slope15) <- names(mean_slope30) <- c("x1", "x2", "x3")

vcov_current <- vcov(LMfit_cox_first_current)
vcov_slope15 <- vcov(LMfit_cox_first_slope15)
vcov_slope30 <- vcov(LMfit_cox_first_slope30)
names(vcov_current) <- names(vcov_slope15) <- names(vcov_slope30) <- c("x1", "x2", "x3")

plotdat <- expand.grid(t = seq(0, 6, by = 0.1),
                       type = c("current", "slope 15", "slope 30")) %>%
  mutate(form = sprintf("~ x1 + %f * x2 + %f * x3", t, t^2),
         lHR = ifelse(type == "current", coef(LMfit_cox_first_current)[["CRI"]] + coef(LMfit_cox_first_current)[["CRIt"]] * t +
                        coef(LMfit_cox_first_current)[["CRIt2"]] * t^2,
                      ifelse(type == "slope 15", coef(LMfit_cox_first_slope15)[["sCRI_15"]] + coef(LMfit_cox_first_slope15)[["sCRI_15t"]] * t +
                               coef(LMfit_cox_first_slope15)[["sCRI_15t2"]] * t^2,
                             coef(LMfit_cox_first_slope30)[["sCRI_30"]] + coef(LMfit_cox_first_slope30)[["sCRI_30t"]] * t +
                               coef(LMfit_cox_first_slope30)[["sCRI_30t2"]] * t^2)),
         se = ifelse(type == "current", sapply(form, function(x) deltamethod(g = as.formula(x), mean = mean_current, cov = vcov_current)),
                     ifelse(type == "slope 15", sapply(form, function(x) deltamethod(g = as.formula(x), mean = mean_slope15, cov = vcov_slope15)),
                            sapply(form, function(x) deltamethod(g = as.formula(x), mean = mean_slope30, cov = vcov_slope30)))),
         lHR_lo = lHR - 1.96 * se,
         lHR_hi = lHR + 1.96 * se)

ggplot(data = subset(plotdat, t <= 3), aes(x = t)) +
  geom_point(aes(y = lHR)) +
  geom_errorbar(aes(ymin = lHR_lo, ymax = lHR_hi)) +
  facet_wrap(type ~ .) +
  theme_bw()

# simulation --------------------------------------------------------------

t1 <- seq(from = 0, to = 48, by = 1)
t2 <- seq(from = 0, to = 48, by = 0.25)

pp <- 50 * cos(2 * pi * (1/20) * t2 + 0.6 * pi) + 50 + rnorm(n = length(t2), mean = 0, sd = 10)
#cri <- pmax(pmin((pp * 0.3 + rnorm(length(t2), mean = 0, sd = 1))/40 + 0.2, 1), 0)
cri <- pmax(pmin(rnorm(n = length(t2), mean = 0.5, sd = 0.25), 1), 0)
plot(x = t2, y = pp, type = "l")
plot(x = t2, y = cri, type = "l", col = "red")

pp_00 <- pp[t2 %in% t1]
plot(x = t1, y = pp_00, type = "l")

dat_pp <- data.frame(t = t1, pp = pp_00)
dat_cri <- data.frame(t = t2, cri = cri)

get_closest_pp_sim <- function(tp, pp_dat = dat_pp, tol = 15/60) {
  tmp <- subset(pp_dat, (t >= tp) & (t <= (tp + tol))) %>%
    arrange(t)
  if (nrow(tmp) == 0) {
    output <- NA
  } else {
    output <- cbind(PP = first(tmp$pp),
                    PP_t = last(tmp$t))
  }
  return(output)
}

pp_15 <- do.call(rbind, lapply(1:nrow(dat_cri), function(i) {get_closest_pp_sim(tp = (dat_cri$t[i] + 0.25), pp_dat = dat_pp, tol = 15/60)}))
pp_30 <- do.call(rbind, lapply(1:nrow(dat_cri), function(i) {get_closest_pp_sim(tp = (dat_cri$t[i] + 0.30), pp_dat = dat_pp, tol = 15/60)}))
pp_01 <- do.call(rbind, lapply(1:nrow(dat_cri), function(i) {get_closest_pp_sim(tp = (dat_cri$t[i] + 1), pp_dat = dat_pp, tol = 15/60)}))
pp_02 <- do.call(rbind, lapply(1:nrow(dat_cri), function(i) {get_closest_pp_sim(tp = (dat_cri$t[i] + 2), pp_dat = dat_pp, tol = 15/60)}))
pp_03 <- do.call(rbind, lapply(1:nrow(dat_cri), function(i) {get_closest_pp_sim(tp = (dat_cri$t[i] + 3), pp_dat = dat_pp, tol = 15/60)}))
pp_04 <- do.call(rbind, lapply(1:nrow(dat_cri), function(i) {get_closest_pp_sim(tp = (dat_cri$t[i] + 4), pp_dat = dat_pp, tol = 15/60)}))
pp_05 <- do.call(rbind, lapply(1:nrow(dat_cri), function(i) {get_closest_pp_sim(tp = (dat_cri$t[i] + 5), pp_dat = dat_pp, tol = 15/60)}))
pp_06 <- do.call(rbind, lapply(1:nrow(dat_cri), function(i) {get_closest_pp_sim(tp = (dat_cri$t[i] + 6), pp_dat = dat_pp, tol = 15/60)}))
pp_07 <- do.call(rbind, lapply(1:nrow(dat_cri), function(i) {get_closest_pp_sim(tp = (dat_cri$t[i] + 7), pp_dat = dat_pp, tol = 15/60)}))
pp_10 <- do.call(rbind, lapply(1:nrow(dat_cri), function(i) {get_closest_pp_sim(tp = (dat_cri$t[i] + 10), pp_dat = dat_pp, tol = 15/60)}))

dat <- data.frame(cbind(dat_cri,
                        PP.m025 = pp_15[,"PP"],
                        PP.m050 = pp_30[,"PP"],
                        PP.m100 = pp_01[,"PP"],
                        PP.m200 = pp_02[,"PP"],
                        PP.m300 = pp_03[,"PP"],
                        PP.m400 = pp_04[,"PP"],
                        PP.m500 = pp_05[,"PP"],
                        PP.m600 = pp_06[,"PP"],
                        PP.m700 = pp_07[,"PP"],
                        PP.m1000 = pp_10[,"PP"]))

cripp_dat <- gather(dat, key = "tlag", value = "PP", -1, -2) %>%
  mutate(pp.low = as.numeric(I(PP <= 20))) %>%
  arrange(t) %>%
  filter(!is.na(PP))

myfit <- function(formula, data) {
  fitdat <- data[!is.na(data[, all.vars(formula)[2]]),]
  fit <- glm(formula = formula, family = binomial, data = fitdat)
  return(fit)
}

mysummary <- function(fit) {
  est <- coef(fit)[[2]]
  se  <- summary(fit)$coef[, 2][2]
  pval <- summary(fit)$coef[2, 4]
  or <- 1/exp(est)
  or_lo <- 1/exp(est + 1.96 * se)
  or_hi <- 1/exp(est - 1.96 * se)

  ## return
  output <- t(c(or = or, or_lo = or_lo, or_hi = or_hi, p = pval))
  rownames(output) <- all.vars(formula(fit))[2]
  return(output)
}

fit.m025 <- myfit(formula = pp.low ~ I(cri * 10), data = subset(cripp_dat, tlag == "PP.m025"))
fit.m050 <- myfit(formula = pp.low ~ I(cri * 10), data = subset(cripp_dat, tlag == "PP.m050"))
fit.m100 <- myfit(formula = pp.low ~ I(cri * 10), data = subset(cripp_dat, tlag == "PP.m100"))
fit.m200 <- myfit(formula = pp.low ~ I(cri * 10), data = subset(cripp_dat, tlag == "PP.m200"))
fit.m300 <- myfit(formula = pp.low ~ I(cri * 10), data = subset(cripp_dat, tlag == "PP.m300"))
fit.m400 <- myfit(formula = pp.low ~ I(cri * 10), data = subset(cripp_dat, tlag == "PP.m400"))
fit.m500 <- myfit(formula = pp.low ~ I(cri * 10), data = subset(cripp_dat, tlag == "PP.m500"))
fit.m600 <- myfit(formula = pp.low ~ I(cri * 10), data = subset(cripp_dat, tlag == "PP.m600"))
fit.m700 <- myfit(formula = pp.low ~ I(cri * 10), data = subset(cripp_dat, tlag == "PP.m700"))
fit.m1000 <- myfit(formula = pp.low ~ I(cri * 10), data = subset(cripp_dat, tlag == "PP.m1000"))

cripp_tab <- do.call(rbind, lapply(list(fit.m025, fit.m050, fit.m100, fit.m200, fit.m300, fit.m400, fit.m500, fit.m600, fit.m700, fit.m1000), function(x) mysummary(x)))
round(cripp_tab, 2)

with(dat, cor.test(x = cri, y = PP.m025))
with(dat, cor.test(x = cri, y = PP.m050))
with(dat, cor.test(x = cri, y = PP.m100))
with(dat, cor.test(x = cri, y = PP.m200))
with(dat, cor.test(x = cri, y = PP.m300))
with(dat, cor.test(x = cri, y = PP.m400))
with(dat, cor.test(x = cri, y = PP.m500))
with(dat, cor.test(x = cri, y = PP.m600))
with(dat, cor.test(x = cri, y = PP.m700))
with(dat, cor.test(x = cri, y = PP.m1000))
