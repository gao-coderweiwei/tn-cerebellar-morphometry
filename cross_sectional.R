# Input: baseline, the prepared baseline analysis data.frame (86 TN, 40 HC).
# Output: cross_sectional_results, a named list of statistical result tables.
# Additional dependency: boot. No input data are included in this package.
source("statistical_functions.R")
config <- list(seed = 20260927L, bootstrap_reps = 5000L, alpha = 0.05)
options(contrasts = c("contr.treatment", "contr.poly"))
set.seed(config$seed)
results <- list()

roi <- data.frame(roi_order = 1:12,
  roi = c("Lobules I-II", "Lobule III", "Lobule IV", "Lobule V", "Lobule VI", "Crus I",
          "Crus II", "Lobule VIIb", "Lobule VIIIa", "Lobule VIIIb", "Lobule IX", "Lobule X"),
  stem = c("Lobules I-II", "Lobule III", "Lobule IV", "Lobule V", "Lobule VI", "Crus I",
           "Crus II", "Lobule VIIB", "Lobule VIIIA", "Lobule VIIIB", "Lobule IX", "Lobule X"))
outcomes <- do.call(rbind, lapply(c("volume", "thickness"), function(m) {
  data.frame(roi, metric = m, outcome = paste0(roi$stem,
     if (m == "volume") " total volume cm3" else " mean thickness mm"))
}))
dat <- as.data.frame(baseline)
dat$group <- factor(dat$group, levels = c("HC", "TN"))
dat$sex <- factor(dat$sex, levels = c("Female", "Male"))
dat$icv <- dat[["ICV total volume cm3"]]
dat$duration <- dat$Disease_duration
dat$pain_side <- factor(dat$Pain_side, levels = c("R", "L"), labels = c("Right", "Left"))
dat$sindou <- factor(ifelse(dat$group == "HC", NA,
                            ifelse(dat$Sindou_grade <= 1, "0-1", as.character(dat$Sindou_grade))),
                      levels = c("0-1", "2", "3"))
dat$age_c <- dat$age - mean(dat$age)
dat$age_c2 <- dat$age_c^2
dat$icv_c <- dat$icv - mean(dat$icv)
knots <- unname(quantile(dat$age, c(.1, .5, .9), type = 7))
basis <- splines::ns(dat$age, knots = knots[2], Boundary.knots = knots[c(1,3)], intercept = FALSE)
# Re-express the natural-spline space with explicit linear age and one nonlinear term.
# This allows a one-df HC3 test of nonlinearity without changing fitted values.
candidate <- apply(basis, 2, function(b) residuals(lm(b ~ dat$age_c)))
j <- which.max(colSums(candidate^2))
dat$age_rcs_nonlinear <- candidate[,j] / sd(candidate[,j])
store_result(continuous_summary(dat,c("age","icv","QC")),"baseline_continuous")
store_result(as.data.frame(table(group=dat$group,sex=dat$sex)),"baseline_sex")
baseline_tests <- do.call(rbind,lapply(c("age","icv","QC"), function(v) {
  t <- t.test(dat[[v]][dat$group=="TN"],dat[[v]][dat$group=="HC"])
  data.frame(variable=v, test="Welch t; TN-HC; descriptive",estimate=unname(t$estimate[1]-t$estimate[2]),
             conf_low=t$conf.int[1],conf_high=t$conf.int[2],p_value=t$p.value)
}))
store_result(baseline_tests,"baseline_descriptive_tests")

age_results <- list(); model_store <- list(); diag_list <- list(); point_list <- list(); vif_list <- list()
nonlinear_tests <- list()
for (i in seq_len(nrow(outcomes))) {
  info <- outcomes[i,]
  for (age_model in c("linear","quadratic","spline")) {
    obj <- fit_group(dat,info$outcome,info$metric,age_model)
    key <- paste(info$metric,info$roi,age_model,sep="::")
    model_store[[key]] <- obj
    age_results[[key]] <- result_group(obj,info,age_model)
    if (age_model != "linear") {
      term <- if (age_model == "quadratic") "age_c2" else "age_rcs_nonlinear"
      tt <- hc3_table(obj$fit)
      nonlinear_tests[[key]] <- data.frame(roi=info$roi,metric=info$metric,age_model=age_model,
        tt[tt$term == term,c("estimate","std_error","conf_low","conf_high","statistic","df","p_value")])
    }
    if (age_model == "quadratic") {
      dg <- diagnostic_tables(obj$fit,obj$data$participant_id)
      diag_list[[key]] <- data.frame(roi=info$roi,metric=info$metric,dg$summary)
      point_list[[key]] <- data.frame(roi=info$roi,metric=info$metric,dg$details)
      vif_list[[key]] <- data.frame(roi=info$roi,metric=info$metric,dg$vif)
    }
  }
}
age_results <- bh(dplyr::bind_rows(age_results),c("age_model","metric"),12)
primary_results <- age_results[age_results$age_model=="quadratic",]
store_result(primary_results,"primary_results")
store_result(age_results,"age_model_results")
store_result(dplyr::bind_rows(nonlinear_tests),"age_nonlinearity_descriptive_tests")
store_result(dplyr::bind_rows(diag_list),"primary_model_diagnostics")
store_result(dplyr::bind_rows(point_list),"primary_influence_points")
store_result(dplyr::bind_rows(vif_list),"primary_vif")
supportive_info <- data.frame(roi_order=1:3,roi=c("Whole cerebellum","Gray matter","White matter"),
 metric="volume",outcome=paste0(c("Cerebellum","Grey Matter","White Matter")," total volume cm3"))
supportive <- dplyr::bind_rows(lapply(1:3,function(i) {
 info<-supportive_info[i,];result_group(fit_group(dat,info$outcome,"volume"),info,"quadratic")
}))
supportive <- bh(supportive,"metric",3)
store_result(supportive,"supportive_global_volumes")

icv_rows <- list(); qc_rows <- list(); correlation_rows <- list()
for (i in which(outcomes$metric == "thickness")) {
  info <- outcomes[i,]
  icv_rows[[info$roi]] <- result_group(fit_group(dat,info$outcome,"thickness",icv=TRUE),info,"quadratic")
  qc_rows[[info$roi]] <- result_group(fit_group(dat,info$outcome,"thickness",qc=TRUE),info,"quadratic")
  for (predictor in c("icv","QC")) for (population in c("All","HC","TN")) {
    d <- if (population == "All") dat else dat[dat$group == population,]
    test <- cor.test(d[[info$outcome]],d[[predictor]],method="spearman",exact=FALSE)
    correlation_rows[[paste(info$roi,predictor,population)]] <- data.frame(
      roi=info$roi,predictor=predictor,population=population,n=nrow(d),
      spearman_rho=unname(test$estimate),p_value=test$p.value)
  }
}
icv_sensitivity <- bh(dplyr::bind_rows(icv_rows),"metric",12)
qc_sensitivity <- bh(dplyr::bind_rows(qc_rows),"metric",12)
store_result(icv_sensitivity,"thickness_icv_sensitivity")
store_result(qc_sensitivity,"qc_adjusted_thickness")
store_result(bh(dplyr::bind_rows(correlation_rows),c("predictor","population"),12),
             "qc_icv_correlations")

selected <- primary_results[primary_results$significant,]
loo_rows <- list(); loo_details <- list(); boot_rows <- list()
# The statistic delegates OLS to stats::lm.fit; strata preserve each group's sample size.
make_boot_statistic <- function(X,y) {
  force(X); force(y)
  function(data,indices) {
    fit <- stats::lm.fit(X[indices,,drop=FALSE],y[indices])
    if (fit$rank < ncol(X)) return(NA_real_)
    unname(fit$coefficients["groupTN"])
  }
}
for (i in seq_len(nrow(selected))) {
  row <- selected[i,]
  key <- paste(row$metric,row$roi,"quadratic",sep="::")
  obj <- model_store[[key]]; d <- obj$data
  X <- model.matrix(obj$fit); y <- model.response(model.frame(obj$fit))
  loo <- vapply(seq_len(nrow(d)),function(j) {
    f <- lm.fit(X[-j,,drop=FALSE],y[-j]);stopifnot(f$rank==ncol(X))
    unname(f$coefficients["groupTN"])
  },numeric(1))
  loo_details[[key]] <- data.frame(roi=row$roi,metric=row$metric,
    omitted_participant=d$participant_id,estimate=loo,full_estimate=row$estimate,
    change=loo-row$estimate)
  loo_rows[[key]] <- data.frame(roi=row$roi,metric=row$metric,full_estimate=row$estimate,
    loo_min=min(loo),loo_max=max(loo),max_absolute_change=max(abs(loo-row$estimate)),
    sign_reversal=any(sign(loo)!=sign(row$estimate)),
    most_influential_omission=d$participant_id[which.max(abs(loo-row$estimate))])
  set.seed(config$seed+i)
  boot_fit <- boot::boot(d,statistic=make_boot_statistic(X,y),R=config$bootstrap_reps,
                         strata=d$group,stype="i")
  if (any(!is.finite(boot_fit$t))) stop("Invalid bootstrap replicate; investigate rank")
  ci <- boot::boot.ci(boot_fit,conf=.95,type="perc")$percent[1,4:5]
  boot_rows[[key]] <- data.frame(roi=row$roi,metric=row$metric,estimate=row$estimate,
     R=config$bootstrap_reps,seed=config$seed+i,bootstrap_mean=mean(boot_fit$t),
     bootstrap_se=sd(boot_fit$t[,1]),conf_low=ci[1],conf_high=ci[2],
     ci_method="Group-stratified pairs bootstrap percentile; unadjusted; selected outcomes")
}
loo_summary <- if(length(loo_rows)) dplyr::bind_rows(loo_rows) else
 data.frame(roi=character(),metric=character(),full_estimate=numeric(),loo_min=numeric(),
 loo_max=numeric(),max_absolute_change=numeric(),sign_reversal=logical(),most_influential_omission=character())
boot_summary <- if(length(boot_rows)) dplyr::bind_rows(boot_rows) else
 data.frame(roi=character(),metric=character(),estimate=numeric(),R=integer(),seed=integer(),
 bootstrap_mean=numeric(),bootstrap_se=numeric(),conf_low=numeric(),conf_high=numeric(),ci_method=character())
store_result(loo_summary,"leave_one_out_summary")
store_result(if(length(loo_details)) dplyr::bind_rows(loo_details) else
 data.frame(roi=character(),metric=character(),omitted_participant=character(),estimate=numeric(),
 full_estimate=numeric(),change=numeric()),"leave_one_out_subjects")
store_result(boot_summary,"bootstrap_summary")


# TN clinical associations and anatomical laterality.

secondary_required <- c(
  "participant_id", "group", "sex", "age", "icv", "duration",
  "pain_side", "sindou"
)
stopifnot(all(secondary_required %in% names(dat)), nrow(roi) == 12L)
secondary_tn <- as.data.frame(dat[!is.na(dat$group) & dat$group == "TN", ])
secondary_tn$sex <- factor(secondary_tn$sex, levels = c("Female", "Male"))
secondary_tn$pain_side <- factor(secondary_tn$pain_side,
                                 levels = c("Right", "Left"))
secondary_tn$sindou <- factor(secondary_tn$sindou, levels = c("0-1", "2", "3"))

# Fixed transformations use the full TN cohort, before outcome-specific omission.
secondary_age_center <- mean(secondary_tn$age, na.rm = TRUE)
secondary_icv_center <- mean(secondary_tn$icv, na.rm = TRUE)
secondary_duration_center <- mean(secondary_tn$duration, na.rm = TRUE)
secondary_duration_sd <- stats::sd(secondary_tn$duration, na.rm = TRUE)
if (!is.finite(secondary_duration_sd) || secondary_duration_sd <= 0) {
  stop("TN disease duration has no finite, positive standard deviation.")
}
secondary_tn$age_c <- secondary_tn$age - secondary_age_center
secondary_tn$age_c2 <- secondary_tn$age_c^2
secondary_tn$icv_c <- secondary_tn$icv - secondary_icv_center
secondary_tn$duration_z <- (secondary_tn$duration - secondary_duration_center) /
  secondary_duration_sd

secondary_metrics <- data.frame(
  metric = c("volume", "thickness"),
  unit = c("cm3", "mm"),
  bilateral_suffix = c(" total volume cm3", " mean thickness mm"),
  left_suffix = c(" left volume cm3", " left thickness mm"),
  right_suffix = c(" right volume cm3", " right thickness mm"),
  stringsAsFactors = FALSE
)

secondary_fit <- function(y, metric, roi_label, analysis) {
  x <- secondary_tn
  x$y <- y
  predictors <- c("duration_z", "pain_side", "sindou", "age_c", "age_c2", "sex")
  if (metric == "volume") predictors <- c(predictors, "icv_c")
  model_variables <- c("y", predictors)
  valid <- stats::complete.cases(x[, model_variables, drop = FALSE])
  for (variable in model_variables) {
    if (is.numeric(x[[variable]])) valid <- valid & is.finite(x[[variable]])
  }
  x_fit <- x[valid, , drop = FALSE]
  fit <- stats::lm(stats::reformulate(predictors, response = "y"),
                   data = x_fit, na.action = stats::na.fail,
                   model = TRUE, x = TRUE, y = TRUE)
  require_full_rank(fit)
  if (stats::df.residual(fit) < 1L) stop("Insufficient secondary-model residual df.")
  list(fit = fit, coefficients = hc3_table(fit), n = nrow(x_fit))
}

# An explicit contrast matrix avoids parsing factor-level names as expressions.
secondary_sindou_test <- function(fit) {
  terms <- names(stats::coef(fit))
  indices <- match(c("sindou2", "sindou3"), terms)
  if (anyNA(indices)) stop("Both Sindou contrasts are required for the 2-df test.")
  hypothesis <- matrix(0, nrow = 2L, ncol = length(terms),
                       dimnames = list(c("grade2_vs_0_1", "grade3_vs_0_1"), terms))
  hypothesis[cbind(seq_len(2L), indices)] <- 1
  test <- car::linearHypothesis(
    fit, hypothesis.matrix = hypothesis, rhs = c(0, 0),
    vcov. = sandwich::vcovHC(fit, type = "HC3"), test = "F"
  )
  c(statistic = unname(test[2L, "F"]),
    p_value = unname(test[2L, "Pr(>F)"]))
}

secondary_factor_rows <- list()
secondary_coefficient_rows <- list()
secondary_laterality_rows <- list()
secondary_models <- list()

for (metric_index in seq_len(nrow(secondary_metrics))) {
  specification <- secondary_metrics[metric_index, ]
  metric <- specification$metric
  for (roi_index in seq_len(nrow(roi))) {
    roi_label <- as.character(roi$roi[roi_index])
    stem <- as.character(roi$stem[roi_index])
    columns <- paste0(stem, c(specification$bilateral_suffix,
                             specification$left_suffix, specification$right_suffix))
    if (!all(columns %in% names(secondary_tn))) {
      stop("Missing secondary outcome columns: ",
           paste(setdiff(columns, names(secondary_tn)), collapse = ", "))
    }
    clinical <- secondary_fit(secondary_tn[[columns[1L]]], metric,
                              roi_label, "clinical_associations")
    model_key <- paste("clinical", metric, roi_index, sep = "_")
    secondary_models[[model_key]] <- clinical$fit
    clinical_coef <- clinical$coefficients
    for (factor_name in c("duration", "pain_side")) {
      term <- if (factor_name == "duration") "duration_z" else "pain_sideLeft"
      coefficient <- clinical_coef[clinical_coef$term == term, , drop = FALSE]
      if (nrow(coefficient) != 1L) stop("Required secondary coefficient not found: ", term)
      secondary_factor_rows[[length(secondary_factor_rows) + 1L]] <- data.frame(
        roi_order = roi$roi_order[roi_index], roi = roi_label, metric = metric,
        unit = specification$unit, factor = factor_name, n = clinical$n,
        contrast = if (factor_name == "duration") "per_TN_SD_of_duration_years" else "Left_vs_Right_pain",
        estimate = coefficient$estimate, std_error = coefficient$std_error,
        conf_low = coefficient$conf_low, conf_high = coefficient$conf_high,
        statistic = coefficient$statistic^2, statistic_type = "HC3_Wald_F",
        df_num = 1L, df_den = stats::df.residual(clinical$fit),
        p_value = coefficient$p_value, stringsAsFactors = FALSE
      )
    }
    sindou_test <- secondary_sindou_test(clinical$fit)
    secondary_factor_rows[[length(secondary_factor_rows) + 1L]] <- data.frame(
      roi_order = roi$roi_order[roi_index], roi = roi_label, metric = metric,
      unit = specification$unit, factor = "sindou", n = clinical$n,
      contrast = "overall_0_1_vs_2_vs_3", estimate = NA_real_, std_error = NA_real_,
      conf_low = NA_real_, conf_high = NA_real_,
      statistic = sindou_test[["statistic"]], statistic_type = "HC3_Wald_F",
      df_num = 2L, df_den = stats::df.residual(clinical$fit),
      p_value = sindou_test[["p_value"]], stringsAsFactors = FALSE
    )
    clinical_coef$roi_order <- roi$roi_order[roi_index]
    clinical_coef$roi <- roi_label
    clinical_coef$metric <- metric
    clinical_coef$unit <- specification$unit
    clinical_coef$n <- clinical$n
    clinical_coef$factor <- ifelse(
      clinical_coef$term == "duration_z", "duration",
      ifelse(clinical_coef$term == "pain_sideLeft", "pain_side",
             ifelse(grepl("^sindou", clinical_coef$term), "sindou", "covariate"))
    )
    secondary_coefficient_rows[[length(secondary_coefficient_rows) + 1L]] <- clinical_coef

    # L-R is anatomical asymmetry. beta(pain_sideLeft)/2 equals the
    # equal-pain-group-weighted adjusted ipsilateral-minus-contralateral contrast.
    # It is not the observed mean of individual ipsilateral-minus-contralateral values.
    lateral_y <- secondary_tn[[columns[2L]]] - secondary_tn[[columns[3L]]]
    lateral <- secondary_fit(lateral_y, metric, roi_label, "pain_laterality")
    secondary_models[[paste("laterality", metric, roi_index, sep = "_")]] <- lateral$fit
    lateral_coef <- lateral$coefficients[
      lateral$coefficients$term == "pain_sideLeft", , drop = FALSE
    ]
    if (nrow(lateral_coef) != 1L) stop("Pain-side coefficient not found in laterality model.")
    secondary_laterality_rows[[length(secondary_laterality_rows) + 1L]] <- data.frame(
      roi_order = roi$roi_order[roi_index], roi = roi_label, metric = metric,
      unit = specification$unit, n = lateral$n,
      n_right_pain = sum(stats::model.frame(lateral$fit)$pain_side == "Right"),
      n_left_pain = sum(stats::model.frame(lateral$fit)$pain_side == "Left"),
      anatomical_LminusR_pain_Left_vs_Right = lateral_coef$estimate,
      anatomical_contrast_std_error = lateral_coef$std_error,
      anatomical_contrast_conf_low = lateral_coef$conf_low,
      anatomical_contrast_conf_high = lateral_coef$conf_high,
      adjusted_ipsi_minus_contra_equal_weights = lateral_coef$estimate / 2,
      ipsi_minus_contra_std_error = lateral_coef$std_error / 2,
      ipsi_minus_contra_conf_low = lateral_coef$conf_low / 2,
      ipsi_minus_contra_conf_high = lateral_coef$conf_high / 2,
      statistic = lateral_coef$statistic, statistic_type = "HC3_t",
      df = stats::df.residual(lateral$fit), p_value = lateral_coef$p_value,
      stringsAsFactors = FALSE
    )
  }
}

secondary_clinical <- do.call(rbind, secondary_factor_rows)
secondary_laterality <- do.call(rbind, secondary_laterality_rows)
secondary_coefficients <- do.call(rbind, secondary_coefficient_rows)
for (metric in secondary_metrics$metric) {
  clinical_index <- secondary_clinical$metric == metric
  laterality_index <- secondary_laterality$metric == metric
  stopifnot(sum(clinical_index) == 36L, sum(laterality_index) == 12L)
  secondary_clinical$q_value[clinical_index] <- stats::p.adjust(
    secondary_clinical$p_value[clinical_index], method = "BH", n = 36L
  )
  secondary_laterality$q_value[laterality_index] <- stats::p.adjust(
    secondary_laterality$p_value[laterality_index], method = "BH", n = 12L
  )
}
secondary_clinical$fdr_family <- paste0("clinical_", secondary_clinical$metric, "_36")
secondary_laterality$fdr_family <- paste0("laterality_", secondary_laterality$metric, "_12")
secondary_clinical$significant_fdr_0_05 <- secondary_clinical$q_value < 0.05
secondary_laterality$significant_fdr_0_05 <- secondary_laterality$q_value < 0.05

# Coefficient p-values and CIs are unadjusted. For Sindou, the inferential
# decision belongs to the joint 2-df test, not a nominal coefficient p-value.
factor_keys <- paste(secondary_clinical$metric, secondary_clinical$roi,
                     secondary_clinical$factor, sep = "|")
coefficient_keys <- paste(secondary_coefficients$metric, secondary_coefficients$roi,
                          secondary_coefficients$factor, sep = "|")
secondary_coefficients$factor_q_value <- secondary_clinical$q_value[
  match(coefficient_keys, factor_keys)
]
secondary_coefficients$interpretation <- ifelse(
  secondary_coefficients$factor == "sindou",
  "Model coefficient only; interpret joint Sindou test and gated Holm contrasts.",
  ifelse(secondary_coefficients$factor == "covariate",
         "Adjustment coefficient; no discovery inference planned.",
         "Use factor_q_value for the clinical-association FDR decision.")
)

# The two original reference-group Sindou contrasts are interpreted only after
# the FDR-adjusted omnibus test passes. Holm controls these two comparisons.
secondary_sindou_posthoc <- data.frame(
  metric = character(), roi_order = integer(), roi = character(),
  contrast = character(), estimate = numeric(), std_error = numeric(),
  df = numeric(), conf_low = numeric(), conf_high = numeric(),
  statistic = numeric(), p_value = numeric(), p_holm = numeric(),
  omnibus_q_value = numeric(), stringsAsFactors = FALSE
)
eligible <- secondary_clinical[
  secondary_clinical$factor == "sindou" & secondary_clinical$q_value < 0.05, , drop = FALSE
]
if (nrow(eligible) > 0L) {
  for (index in seq_len(nrow(eligible))) {
    result <- eligible[index, ]
    roi_index <- match(result$roi_order, roi$roi_order)
    fit <- secondary_models[[paste("clinical", result$metric, roi_index, sep = "_")]]
    comparisons <- hc3_table(fit)
    comparisons <- comparisons[match(c("sindou2", "sindou3"), comparisons$term),
                               , drop = FALSE]
    new_rows <- data.frame(
      metric = result$metric, roi_order = result$roi_order, roi = result$roi,
      contrast = c("grade2_vs_0_1", "grade3_vs_0_1"), estimate = comparisons$estimate,
      std_error = comparisons$std_error, df = comparisons$df,
      conf_low = comparisons$conf_low, conf_high = comparisons$conf_high,
      statistic = comparisons$statistic, p_value = comparisons$p_value,
      p_holm = stats::p.adjust(comparisons$p_value, method = "holm", n = 2L),
      omnibus_q_value = result$q_value, stringsAsFactors = FALSE
    )
    secondary_sindou_posthoc <- rbind(secondary_sindou_posthoc, new_rows)
  }
}

store_result(secondary_clinical, "secondary_clinical_factor_tests")
store_result(secondary_coefficients, "secondary_clinical_coefficients")
store_result(secondary_laterality, "secondary_laterality")
store_result(secondary_sindou_posthoc, "secondary_sindou_posthoc")

design_cells <- as.data.frame(table(pain_side = secondary_tn$pain_side,
                                   sindou = secondary_tn$sindou))
store_result(design_cells, "secondary_design_cells")

cross_sectional_results <- results
