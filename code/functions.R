
# calculate fluid ---------------------------------------------------------

get_fluid_t_each <- function(data, t0, t1) {
  tmp <- subset(data, (t_end > t0) & (t_start < t1))
  if (nrow(tmp) == 0) {
    total <- crystalloid <- colloid <- blood <- 0
  } else {
    t02 <- pmax(tmp$t_start, t0)
    t12 <- pmin(tmp$t_end, t1)
    vol <- tmp$volume * (t12 - t02)/tmp$duration
    total <- sum(vol)
    crystalloid <- sum(vol[tmp$type %in% c("RL", "NS", "RingerLactate")])
    colloid <- sum(vol[tmp$type %in% c("0.1", "0.06", "HES10", "HES6")])
    blood <- sum(vol[tmp$type %in% c("DICH LOC MAU", "PRBC", "Platelets", "FFP")])
  }
  return(data.frame(t0 = t0, t1 = t1, total = total, crystalloid = crystalloid, colloid = colloid, blood = blood))
}

get_fluid_t <- function(id, t, type = c("den", "cum"), data) {
  tmp <- subset(data, (StudyNo == id) & (t_end > 0))
  out <- switch(type,
                den = get_fluid_t_each(data = tmp, t0 = t, t1 = t + 1),
                cum = get_fluid_t_each(data = tmp, t0 = 0, t1 = t))
  return(out)
}


# plot CRI ----------------------------------------------------------------

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

plot.cri_bp <- function(id, tmax) {
  tmp_cri <- subset(cridata, (StudyNo == id) & (t >= 0) & (t <= tmax))
  tmp_vital <- subset(vital, (StudyNo == id) & (t >= 0) & (t <= tmax))
  tmp_fluid <- subset(fluid_32dx_0_48, (StudyNo == id) & (t1 > 0) & (t0 < tmax)) %>%
    mutate(total_cat = factor(ifelse(total_per_kg == 0, "0",
                                     ifelse(total_per_kg > 0 & total_per_kg < 5, "0-<5",
                                            ifelse(total_per_kg >= 5 & total_per_kg < 10, "5-<10", "10-20"))),
                              levels = c("0", "0-<5", "5-<10", "10-20")),
           crysd_cat = factor(ifelse(crystalloid_per_kg == 0, "0",
                                     ifelse(crystalloid_per_kg > 0 & crystalloid_per_kg < 5, "0-<5",
                                            ifelse(crystalloid_per_kg >= 5 & crystalloid_per_kg < 10, "5-<10", "10-20"))),
                              levels = c("0", "0-<5", "5-<10", "10-20")),
           colld_cat = factor(ifelse(colloid_per_kg == 0, "0",
                                     ifelse(colloid_per_kg > 0 & colloid_per_kg < 5, "0-<5",
                                            ifelse(colloid_per_kg >= 5 & colloid_per_kg < 10, "5-<10", "10-20"))),
                              levels = c("0", "0-<5", "5-<10", "10-20")))
  tmp_demo <- subset(general, StudyNo == id)
  tmp_text <- paste("Patient ", tmp_cri$StudyNo[1], ", ", tolower(tmp_cri$sex[1]), ", ", tmp_cri$age[1], " yrs old, ", ifelse(tmp_cri$outpatient[1] == 1, "outpatient", "inpatient"),
                    ", enrolled (time 0) on ", tmp_cri$datetime_enrol[1], " (DOI ", tmp_cri$doi_enrol[1], ")", sep = "")
  p <- ggplot(tmp_cri, aes(x = t)) +
    geom_point(data = tmp_vital, aes(y = BloodSBP/200), colour = "grey60") +
    geom_line(data = tmp_vital, aes(y = BloodSBP/200, group = group), colour = "grey60", size = 1) +
    geom_point(data = tmp_vital, aes(y = BloodDBP/200), colour = "grey80") +
    geom_line(data = tmp_vital, aes(y = BloodDBP/200, group = group), colour = "grey80", size = 1) +
    geom_point(data = tmp_vital, aes(y = PP/200), colour = "grey") +
    geom_line(data = tmp_vital, aes(y = PP/200, group = group), colour = "grey", size = 1, linetype = 2) +
    geom_hline(yintercept = 20/200, colour = "grey", linetype = 2) +
    geom_point(aes(y = CRI), colour = "#F8766D") +
    geom_line(aes(y = CRI, group = group), colour = "#F8766D", size = 1) +

    scale_y_continuous(name = "CRI", breaks = seq(from = 0, to = 1, by = 0.2), limits = c(-0.1, 1),
                       sec.axis = sec_axis(~ . * 200, name = "Blood pressure (mmHg)", breaks = seq(from = 0, to = 200, by = 20))) +
    scale_x_continuous(name = "Hours from enrolment", breaks = seq(from = 0, to = tmax, by = 6), limits = c(-3, tmax)) +
    annotate(geom = "text", x = -3, y = 1, label = tmp_text, hjust = 0, vjust = 0) +
    theme(legend.position = "top",
          axis.text.y.left = element_text(color = "#F8766D", face = "bold"),
          axis.title.y.left = element_text(color = "#F8766D", face = "bold"))

  if (nrow(tmp_fluid) > 0) {
    p <- p +
      geom_segment(data = tmp_fluid, aes(y = 0, yend = 0, x = t0, xend = t1, colour = total_cat), size = 2) +
      annotate(geom = "text", x = -1, y = 0, label = "Total", hjust = 1, size = 3) +

      geom_segment(data = tmp_fluid, aes(y = -0.1, yend = -0.1, x = t0, xend = t1, colour = colld_cat), size = 2) +
      annotate(geom = "text", x = -1, y = -0.1, label = "Colloid", hjust = 1, size = 3)
  }

  return(p + scale_colour_manual(name = "Fluid type",
                                 breaks = c("0", "0-<5", "5-<10", "10-20"),
                                 values = c("#FFFFFF", "#00CCFF", "#6699FF", "#0033FF"),
                                 drop = FALSE) +
           theme(legend.key = element_rect(colour = "black")) +
           guides(colour = guide_legend(nrow = 1, byrow = TRUE)))
}

plot.cri_bp2 <- function(id, tmax) {
  tmp_general <- subset(general_sum2, StudyNo == id)
  tmp_cri <- subset(cridata, (StudyNo == id) & (t >= 0) & (t <= tmax))
  tmp_vital <- subset(vital, (StudyNo == id) & (t >= 0) & (t <= tmax))
  tmp_fluid <- subset(fluid_32dx_0_48, (StudyNo == id) & (t1 > 0) & (t0 < tmax)) %>%
    mutate(total_cat = factor(ifelse(total_per_kg == 0, "0",
                                     ifelse(total_per_kg > 0 & total_per_kg < 5, "0-<5",
                                            ifelse(total_per_kg >= 5 & total_per_kg < 10, "5-<10", "10-20"))),
                              levels = c("0", "0-<5", "5-<10", "10-20")),
           crysd_cat = factor(ifelse(crystalloid_per_kg == 0, "0",
                                     ifelse(crystalloid_per_kg > 0 & crystalloid_per_kg < 5, "0-<5",
                                            ifelse(crystalloid_per_kg >= 5 & crystalloid_per_kg < 10, "5-<10", "10-20"))),
                              levels = c("0", "0-<5", "5-<10", "10-20")),
           colld_cat = factor(ifelse(colloid_per_kg == 0, "0",
                                     ifelse(colloid_per_kg > 0 & colloid_per_kg < 5, "0-<5",
                                            ifelse(colloid_per_kg >= 5 & colloid_per_kg < 10, "5-<10", "10-20"))),
                              levels = c("0", "0-<5", "5-<10", "10-20")))
  tmp_demo <- subset(general, StudyNo == id)
  tmp_text <- paste("Patient ", tmp_cri$StudyNo[1], ", ", tolower(tmp_cri$sex[1]), ", ", tmp_cri$age[1], " yrs old, ", ifelse(tmp_cri$outpatient[1] == 1, "outpatient", "inpatient"),
                    ", enrolled (time 0) on ", tmp_cri$datetime_enrol[1], " (DOI ", tmp_cri$doi_enrol[1], ")", sep = "")
  p <- ggplot(tmp_cri, aes(x = t)) +
    geom_point(data = tmp_vital, aes(y = BloodSBP/200), colour = "grey60") +
    geom_line(data = tmp_vital, aes(y = BloodSBP/200, group = group), colour = "grey60", size = 1) +
    geom_point(data = tmp_vital, aes(y = BloodDBP/200), colour = "grey80") +
    geom_line(data = tmp_vital, aes(y = BloodDBP/200, group = group), colour = "grey80", size = 1) +
    geom_point(data = tmp_vital, aes(y = PP/200), colour = "grey") +
    geom_line(data = tmp_vital, aes(y = PP/200, group = group), colour = "grey", size = 1, linetype = 2) +
    geom_hline(yintercept = 20/200, colour = "grey", linetype = 2) +
    geom_point(aes(y = CRI), colour = "#F8766D") +
    geom_line(aes(y = CRI, group = group), colour = "#F8766D", size = 1) +
    geom_vline(data = tmp_general, aes(xintercept = ts0), linetype = 2) +
    geom_text(label = "ts0", x = tmp_general$ts0, y = 0) +
    geom_vline(data = tmp_general, aes(xintercept = ts1), linetype = 2) +
    geom_text(label = "ts1", x = tmp_general$ts1, y = 0) +
    geom_vline(data = tmp_general, aes(xintercept = ts2), linetype = 2) +
    geom_text(label = "ts2", x = tmp_general$ts2, y = 0) +
    geom_vline(data = tmp_general, aes(xintercept = ts3), linetype = 2) +
    geom_text(label = "ts3", x = tmp_general$ts3, y = 0) +

    scale_y_continuous(name = "CRI", breaks = seq(from = 0, to = 1, by = 0.2), limits = c(-0.1, 1),
                       sec.axis = sec_axis(~ . * 200, name = "Blood pressure (mmHg)", breaks = seq(from = 0, to = 200, by = 20))) +
    scale_x_continuous(name = "Hours from enrolment", breaks = seq(from = 0, to = tmax, by = 6), limits = c(-3, tmax)) +
    annotate(geom = "text", x = -3, y = 1, label = tmp_text, hjust = 0, vjust = 0) +
    theme(legend.position = "top",
          axis.text.y.left = element_text(color = "#F8766D", face = "bold"),
          axis.title.y.left = element_text(color = "#F8766D", face = "bold"))

  if (nrow(tmp_fluid) > 0) {
    p <- p +
      geom_segment(data = tmp_fluid, aes(y = 0, yend = 0, x = t0, xend = t1, colour = total_cat), size = 2) +
      annotate(geom = "text", x = -1, y = 0, label = "Total", hjust = 1, size = 3) +

      geom_segment(data = tmp_fluid, aes(y = -0.1, yend = -0.1, x = t0, xend = t1, colour = colld_cat), size = 2) +
      annotate(geom = "text", x = -1, y = -0.1, label = "Colloid", hjust = 1, size = 3)
  }

  return(p + scale_colour_manual(name = "Fluid type",
                                 breaks = c("0", "0-<5", "5-<10", "10-20"),
                                 values = c("#FFFFFF", "#00CCFF", "#6699FF", "#0033FF"),
                                 drop = FALSE) +
           theme(legend.key = element_rect(colour = "black")) +
           guides(colour = guide_legend(nrow = 1, byrow = TRUE)))
}

# get closest CRI ---------------------------------------------------------

## for each vital record, corresponding CRI is the closest value within 10 minutes before that record
get_closest_cri <- function(id, tp, cri_dat = cridata, tol = 10/60) {
  tmp <- subset(cri_dat, StudyNo == id & (t <= tp) & (t >= (tp - tol))) %>%
    arrange(t)
  if (nrow(tmp) == 0) {
    output <- NA
  } else {
    output <- last(tmp$CRI)
  }
  return(output)
}

get_closest_cri2 <- function(id, tp, cri_dat = cridata, tol = 10/60) {
  tmp <- subset(cri_dat, StudyNo == id & (t <= tp) & (t >= (tp - tol))) %>%
    arrange(t)
  if (nrow(tmp) == 0) {
    output <- NA
  } else {
    output <- cbind(CRI = last(tmp$CRI),
                    CRI_t = last(tmp$t))
  }
  return(output)
}

get_slope_cri <- function(id, tp, cri_dat = cridata, dur = 30/60) {
  tmp <- subset(cri_dat, StudyNo == id & (t <= tp) & (t >= (tp - dur))) %>%
    arrange(t)
  if (nrow(tmp) == 0) {
    output <- NA
  } else {
    output <- coef(lm(CRI ~ t, data = tmp))[[2]]
  }
  return(output)
}


get_closest_hr <- function(id, tp, cri_dat = cridata, tol = 10/60) {
  tmp <- subset(cri_dat, StudyNo == id & (t <= tp) & (t >= (tp - tol))) %>%
    arrange(t)
  if (nrow(tmp) == 0) {
    output <- NA
  } else {
    output <- last(tmp$HR)
  }
  return(output)
}

get_closest_hr2 <- function(id, tp, cri_dat = cridata, tol = 10/60) {
  tmp <- subset(cri_dat, StudyNo == id & (t <= tp) & (t >= (tp - tol))) %>%
    arrange(t)
  if (nrow(tmp) == 0) {
    output <- NA
  } else {
    output <- cbind(HR = last(tmp$HR),
                    HR_t = last(tmp$t))
  }
  return(output)
}


# calculate correlation ---------------------------------------------------

get_cor <- function(var1 = c("CRI", "HR", "fluid_total_kg", "fluid_colloid_kg", "BloodSBP", "BloodDBP", "PP"),
                    var2 = c("CRI", "HR", "fluid_total_kg", "fluid_colloid_kg", "BloodSBP", "BloodDBP", "PP"),
                    tlag = 0:24,
                    id) {
  tmp <- subset(vital, StudyNo == id, t >= 0)
  #browser()
  out <- do.call(rbind, lapply(1:length(tlag), function(l) {
    x1 <- tmp[, var1][[1]]
    x2 <- switch(var2,
                 CRI = sapply(1:nrow(tmp), function(i) {get_closest_cri(id = tmp$StudyNo[i], tp = (tmp$t[i] + tlag[l]))}),
                 HR  = sapply(1:nrow(tmp), function(i) {get_closest_hr(id = tmp$StudyNo[i], tp = (tmp$t[i] + tlag[l]))}),
                 fluid_total_kg = sapply(1:nrow(tmp), function(i) {get_fluid_t(id = tmp$StudyNo[i], t = (tmp$t[i] + tlag[l]), type = "den", data = fluid_32dx)$total})/tmp$iwt[1],
                 fluid_colloid_kg = sapply(1:nrow(tmp), function(i) {get_fluid_t(id = tmp$StudyNo[i], t = (tmp$t[i] + tlag[l]), type = "den", data = fluid_32dx)$colloid})/tmp$iwt[1],
                 BloodSBP = sapply(1:nrow(tmp), function(i) {ifelse(length(tmp$BloodSBP[tmp$t == (tmp$t[i] + tlag[l])]) == 0, NA, tmp$BloodSBP[tmp$t == (tmp$t[i] + tlag[l])])}),
                 BloodDBP = sapply(1:nrow(tmp), function(i) {ifelse(length(tmp$BloodDBP[tmp$t == (tmp$t[i] + tlag[l])]) == 0, NA, tmp$BloodDBP[tmp$t == (tmp$t[i] + tlag[l])])}),
                 PP = sapply(1:nrow(tmp), function(i) {ifelse(length(tmp$PP[tmp$t == (tmp$t[i] + tlag[l])]) == 0, NA, tmp$PP[tmp$t == (tmp$t[i] + tlag[l])])}))
    tmp_test <- cor.test(x1, x2, method = "pearson")
    return(cbind(cor_mean = tmp_test$estimate,
                 cor_lo   = tmp_test$conf.int[1],
                 cor_up   = tmp_test$conf.int[2]))
  }))
  return(data.frame(id = id, rel = paste(var1, var2, sep = " vs. "), lag = tlag, out))
}

get.rho <- function(x, y, method = "pearson") {
  ## get pearson correlation
  rho <- cor(x, y, method = method, use = "complete.obs")
  ## get z-transformation (http://www2.sas.com/proceedings/sugi31/170-31.pdf)
  z <- 0.5 * log((1 + rho)/(1 - rho))
  ## output
  return(c(rho = rho, z = z))
}

get_cor2 <- function(x, y, id = "StudyNo", data, method = "pearson", nBoot = 1000, ci = 0.95){
  #browser()

  ## get estimate of rho
  rho <- get.rho(x = data[, x], y = data[, y], method = method)

  ## get standard error of z
  cluster <- data[, id]
  if (is.factor(cluster)) cluster <- drop.levels(cluster)
  z_boot <- matrix(NA, nrow = nBoot, ncol = 2)

  #- set-up for cluster resampling
  clusters <- unique(cluster)
  nclust <- length(clusters)
  obsno <- split(1:length(cluster), cluster) # n = # total records, not subjects

  #- get the actual bootstrap samples and determine coefficients
  pb <- txtProgressBar(min = 0, max = nBoot, style = 3)
  for (i in 1:nBoot){
    setTxtProgressBar(pb, i)
    j <- sample(1:nclust, nclust, replace = TRUE)
    obs <- unlist(obsno[j])
    data.bsamp <- data[obs,]
    z_boot[i,] <- get.rho(x = data.bsamp[, x], y = data.bsamp[, y], method = method)
  }
  close(pb)
  se_z <- sd(z_boot[, 2])

  # prepare final summary
  rho["se_z"] <- se_z
  rho["z_lo"] <- rho["z"] - qnorm(0.5 * (ci + 1)) * se_z
  rho["z_hi"] <- rho["z"] + qnorm(0.5 * (ci + 1)) * se_z
  rho["rho_lo"] <- (exp(2 * rho["z_lo"]) - 1)/(exp(2 * rho["z_lo"]) + 1)
  rho["rho_hi"] <- (exp(2 * rho["z_hi"]) - 1)/(exp(2 * rho["z_hi"]) + 1)
  rho["p"] <- 2 * pnorm(-abs(rho["z"]/rho["se_z"]))

  return(rho)
}
