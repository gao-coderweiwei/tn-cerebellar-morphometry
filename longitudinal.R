# Inputs: baseline and followup, prepared analysis data.frames.
# Follow-up IDs identify the 49 patients with paired MRI.
# Output: longitudinal_results, a named list of statistical result tables.
# Additional dependency: boot. No MRI interval adjustment is performed.
source("statistical_functions.R")
config <- list(seed = 20260928L, bootstrap_reps = 5000L, alpha = 0.05)
set.seed(config$seed)
results <- list()
baseline <- as.data.frame(baseline)
followup <- as.data.frame(followup)
paired_id <- as.character(followup$participant_id)
if (anyDuplicated(paired_id) || anyDuplicated(baseline$participant_id))
  stop("Patient IDs must be unique for paired analysis")
pair_index <- match(paired_id, baseline$participant_id)
if (anyNA(pair_index)) stop("Follow-up patient has no matching baseline record")
pre <- baseline[pair_index, , drop = FALSE]
post <- followup

stems <- c("Lobules I-II", "Lobule III", "Lobule IV", "Lobule V", "Lobule VI",
           "Crus I", "Crus II", "Lobule VIIB", "Lobule VIIIA", "Lobule VIIIB", "Lobule IX", "Lobule X")
outcomes <- dplyr::bind_rows(lapply(c("volume", "thickness"), function(m) data.frame(
  roi_order = 1:12, roi = c("Lobules I-II", "Lobule III", "Lobule IV", "Lobule V", "Lobule VI",
    "Crus I", "Crus II", "Lobule VIIb", "Lobule VIIIa", "Lobule VIIIb", "Lobule IX", "Lobule X"), metric = m, unit = ifelse(m == "volume", "cm3", "mm"),
  outcome = paste(stems, if (m == "volume") "total volume cm3" else "mean thickness mm"))))
global_outcomes <- data.frame(roi_order = 1:3, roi = c("Whole cerebellum", "Gray matter", "White matter"),
  metric = "global_volume", unit = "cm3", outcome = paste(c("Cerebellum", "Grey Matter", "White Matter"), "total volume cm3"))
all_outcomes <- dplyr::bind_rows(outcomes, global_outcomes)
changes <- dplyr::bind_rows(lapply(seq_len(nrow(all_outcomes)), function(i) {
  info <- all_outcomes[i, ]
  a <- pre[[info$outcome]]; b <- post[[info$outcome]]
  data.frame(participant_id = paired_id, roi_order = info$roi_order, roi = info$roi,
    metric = info$metric, unit = info$unit, baseline = a, followup = b,
    change = b - a, percent_change = ifelse(a > 0, 100 * (b - a) / a, NA_real_),
    baseline_QC = pre$QC, followup_QC = post$QC,
    change_QC = post$QC - pre$QC,
    baseline_ICV = pre[["ICV total volume cm3"]], followup_ICV = post[["ICV total volume cm3"]],
    change_ICV = post[["ICV total volume cm3"]] - pre[["ICV total volume cm3"]])
}))
store_result(changes, "individual_changes")

# Primary paired t-tests and sensitivity analyses use all 12 lobules.
primary_rows <- wilcox_rows <- boot_rows <- loo_rows <- qc_rows <- list()
boot_stat <- function(d, indices) {
  x <- d[indices, , drop = FALSE]
  c(mean_change = mean(x$change), median_change = median(x$change),
    mean_percent_change = mean(x$percent_change), median_percent_change = median(x$percent_change),
    paired_dz = mean(x$change) / sd(x$change))
}
for (i in seq_len(nrow(all_outcomes))) {
  info <- all_outcomes[i, ]
  d <- changes[changes$roi == info$roi & changes$metric == info$metric, ]
  d <- d[complete.cases(d[, c("baseline", "followup", "change", "percent_change")]), ]
  if (nrow(d) < 3 || sd(d$change) == 0) stop("Insufficient variable paired observations: ", info$roi)
  fit <- t.test(d$followup, d$baseline, paired = TRUE, conf.level = .95)
  primary_rows[[i]] <- data.frame(info, n = nrow(d), n_missing_pairs = length(paired_id) - nrow(d),
    baseline_mean = mean(d$baseline), baseline_sd = sd(d$baseline), followup_mean = mean(d$followup),
    followup_sd = sd(d$followup), mean_change = mean(d$change), change_sd = sd(d$change),
    conf_low = fit$conf.int[1], conf_high = fit$conf.int[2], median_change = median(d$change),
    change_q25 = unname(quantile(d$change, .25)), change_q75 = unname(quantile(d$change, .75)),
    mean_percent_change = mean(d$percent_change), median_percent_change = median(d$percent_change),
    percent_change_q25 = unname(quantile(d$percent_change, .25)),
    percent_change_q75 = unname(quantile(d$percent_change, .75)),
    percent_change_of_group_mean = 100 * mean(d$change) / mean(d$baseline),
    paired_dz = mean(d$change) / sd(d$change), statistic = unname(fit$statistic),
    df = unname(fit$parameter), p_value = fit$p.value,
    difference_shapiro_p = shapiro.test(d$change)$p.value,
    test = "Two-sided paired t-test", change_definition = "followup minus baseline",
    ci_note = "Unadjusted 95% confidence interval", row.names = NULL)
  nonzero <- d$change[d$change != 0]
  exact_ok <- length(nonzero) < 50L && !anyDuplicated(abs(nonzero)) && length(nonzero) == nrow(d)
  w <- wilcox.test(d$followup, d$baseline, paired = TRUE, exact = exact_ok, correct = !exact_ok)
  wilcox_rows[[i]] <- data.frame(info, n = nrow(d), median_change = median(d$change),
    statistic = unname(w$statistic), p_value = w$p.value, exact = exact_ok,
    test = "Two-sided paired Wilcoxon signed-rank", row.names = NULL)
  set.seed(config$seed + i)
  b <- boot::boot(data = d, statistic = boot_stat, R = config$bootstrap_reps, stype = "i")
  for (j in seq_along(b$t0)) {
    if (any(!is.finite(b$t[, j]))) stop("Nonfinite bootstrap estimate")
    ci <- boot::boot.ci(b, type = "perc", conf = .95, index = j)$percent[1, 4:5]
    boot_rows[[paste(i, j)]] <- data.frame(info, estimand = names(b$t0)[j],
      estimate = b$t0[j], conf_low = ci[1], conf_high = ci[2],
      bootstrap_se = sd(b$t[, j]), n = nrow(d), R = config$bootstrap_reps,
      seed = config$seed + i, ci_method = "Patient-level percentile bootstrap; unadjusted 95%", row.names = NULL)
  }
  if (info$metric != "global_volume") {
    for (k in seq_len(nrow(d))) {
      dk <- d[-k, ]; tk <- t.test(dk$change)
      loo_rows[[paste(i, k)]] <- data.frame(info, omitted_participant = d$participant_id[k],
        n = nrow(dk), full_mean_change = mean(d$change), mean_change = mean(dk$change),
        conf_low = tk$conf.int[1], conf_high = tk$conf.int[2], p_value = tk$p.value)
    }
    f <- lm(change ~ change_ICV + change_QC, data = d, na.action = na.fail)
    ct <- hc3_table(f)
    qc_rows[[i]] <- data.frame(info[rep(1L, nrow(ct)), , drop = FALSE], n = nobs(f), ct,
      model = "change ~ change_ICV + change_QC; HC3",
      intercept_interpretation = "Expected imaging change at zero measured ICV and QC change",
      row.names = NULL)
  }
}
raw_primary <- dplyr::bind_rows(primary_rows)
primary <- bh(raw_primary[raw_primary$metric != "global_volume", ], "metric", 12)
global <- bh(raw_primary[raw_primary$metric == "global_volume", ], "metric", 3)
store_result(primary, "primary_paired_results")
store_result(global, "supportive_global_volumes")
w_all <- dplyr::bind_rows(wilcox_rows)
store_result(bh(w_all[w_all$metric != "global_volume", ], "metric", 12), "wilcoxon_sensitivity")
store_result(dplyr::bind_rows(boot_rows), "bootstrap_summary")
loo <- bh(dplyr::bind_rows(loo_rows), c("metric", "omitted_participant"), 12)
store_result(loo, "leave_one_out_subjects")
loo_summary <- loo |> dplyr::group_by(roi_order, roi, metric, unit) |>
  dplyr::summarise(full_mean_change = dplyr::first(full_mean_change),
    min_change = min(mean_change), max_change = max(mean_change),
    largest_absolute_shift = max(abs(mean_change - full_mean_change)),
    sign_reversal = any(sign(mean_change) != sign(full_mean_change)),
    min_q = min(q_value), max_q = max(q_value),
    significant_omissions = sum(significant), total_omissions = dplyr::n(), .groups = "drop")
store_result(loo_summary, "leave_one_out_summary")
store_result(bh(dplyr::bind_rows(qc_rows), c("metric", "term"), 12), "icv_qc_change_sensitivity")
qc_summary <- lapply(c("QC", "ICV total volume cm3", "Scale factor"), function(v) {
  a <- pre[[v]]; b <- post[[v]]; t <- t.test(b, a, paired = TRUE)
  data.frame(variable = v, n = length(a), baseline_mean = mean(a), baseline_sd = sd(a),
    followup_mean = mean(b), followup_sd = sd(b), mean_change = mean(b - a),
    conf_low = t$conf.int[1], conf_high = t$conf.int[2], p_value = t$p.value,
    inference = "Processing descriptives; unadjusted p; not test-retest reliability")
})
store_result(dplyr::bind_rows(qc_summary), "processing_metrics_by_timepoint")

# Associations with processing changes use imaging data only.
processing_correlations <- list()
for (i in seq_len(nrow(outcomes))) {
  info <- outcomes[i, ]
  d <- changes[changes$roi == info$roi & changes$metric == info$metric, ]
  for (exposure in c("change_QC", "change_ICV")) {
    dq <- d[complete.cases(d[, c("change", exposure)]), , drop = FALSE]
    dq <- dq[is.finite(dq$change) & is.finite(dq[[exposure]]), , drop = FALSE]
    if (nrow(dq) < 3L) stop("Insufficient processing-correlation observations")
    cq <- stats::cor.test(dq$change, dq[[exposure]], method = "spearman", exact = FALSE)
    processing_correlations[[paste(i, exposure)]] <- data.frame(info, exposure = exposure,
      n = nrow(dq), rho = unname(cq$estimate), p_value = cq$p.value,
      interpretation = "Exploratory processing association; cannot estimate scan-rescan error")
  }
}
store_result(bh(dplyr::bind_rows(processing_correlations), c("metric", "exposure"), 12),
             "processing_change_correlations")

longitudinal_results <- results
