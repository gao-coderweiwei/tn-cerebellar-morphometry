# Inputs: baseline and followup, prepared analysis data.frames.
# Output: selection_results, baseline comparisons by paired MRI availability.
# Welch t tests, Fisher exact tests and continuous standardized differences.
source("statistical_functions.R")
results <- list()
paired_id <- as.character(followup$participant_id)

# Baseline comparisons by paired MRI availability.

selection_data <- as.data.frame(baseline[as.character(baseline$group) == "TN", , drop = FALSE])
selection_data$followup_available <- as.character(selection_data$participant_id) %in% as.character(paired_id)
selection_labels <- c("No paired follow-up dataset", "Available paired follow-up")
selection_data$followup_status <- selection_labels[as.integer(selection_data$followup_available) + 1L]
selection_stems <- c("Lobules I-II", "Lobule III", "Lobule IV", "Lobule V",
                     "Lobule VI", "Crus I", "Crus II", "Lobule VIIB",
                     "Lobule VIIIA", "Lobule VIIIB", "Lobule IX", "Lobule X")
selection_cont_vars <- c("age", "Disease_duration", "QC", "ICV total volume cm3",
                         paste0(selection_stems, " total volume cm3"),
                         paste0(selection_stems, " mean thickness mm"))

selection_finite <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  x[!is.finite(x)] <- NA_real_
  x
}
selection_smd <- function(a, b) {
  denom <- sqrt((stats::var(a, na.rm = TRUE) + stats::var(b, na.rm = TRUE))/2)
  if (!is.finite(denom) || denom == 0) return(NA_real_)
  (mean(a, na.rm = TRUE) - mean(b, na.rm = TRUE))/denom
}
selection_describe <- function(x) {
  ok <- x[!is.na(x)]
  qs <- if (length(ok)) stats::quantile(ok, c(.25, .5, .75), names = FALSE) else rep(NA_real_, 3)
  data.frame(n_total = length(x), n_observed = length(ok), n_missing = sum(is.na(x)),
             mean = if (length(ok)) mean(ok) else NA_real_,
             sd = if (length(ok) > 1L) stats::sd(ok) else NA_real_,
             median = qs[2], q25 = qs[1], q75 = qs[3])
}

selection_cont_rows <- list()
selection_test_rows <- list()
for (selection_var in selection_cont_vars) {
  selection_x <- if (selection_var %in% names(selection_data)) {
    selection_finite(selection_data[[selection_var]])
  } else rep(NA_real_, nrow(selection_data))
  selection_a <- selection_x[selection_data$followup_available]
  selection_b <- selection_x[!selection_data$followup_available]
  selection_smd_value <- selection_smd(selection_a, selection_b)
  for (selection_flag in c(FALSE, TRUE)) {
    selection_cont_rows[[length(selection_cont_rows) + 1L]] <- cbind(
      data.frame(variable = selection_var, followup_status = selection_labels[as.integer(selection_flag) + 1L]),
      selection_describe(selection_x[selection_data$followup_available == selection_flag]),
      standardized_mean_difference = selection_smd_value,
      smd_direction = "Available minus no paired follow-up")
  }
  selection_t <- tryCatch(stats::t.test(selection_a, selection_b, alternative = "two.sided"),
                          error = function(e) NULL)
  selection_test_rows[[length(selection_test_rows) + 1L]] <- data.frame(
    variable = selection_var, variable_type = "continuous", test = "Welch two-sample t test",
    n_available = sum(!is.na(selection_a)), n_no_paired_followup = sum(!is.na(selection_b)),
    mean_difference = if (all(is.na(selection_a)) || all(is.na(selection_b))) NA_real_ else
      mean(selection_a, na.rm = TRUE) - mean(selection_b, na.rm = TRUE),
    ci_low = if (is.null(selection_t)) NA_real_ else selection_t$conf.int[1],
    ci_high = if (is.null(selection_t)) NA_real_ else selection_t$conf.int[2],
    standardized_mean_difference = selection_smd_value,
    p_value_unadjusted = if (is.null(selection_t)) NA_real_ else selection_t$p.value,
    interpretation = "Descriptive only; available minus no paired follow-up; unadjusted 95% CI and P value",
    note = if (is.null(selection_t)) "Test not estimable from available observations" else "")
}
store_result(do.call(rbind, selection_cont_rows), "selection_continuous")

selection_cat_vars <- c("sex", "Pain_side", "Sindou_grade", "Sindou_group", "Pain_type", "Pharm_treatment")
selection_cat_rows <- list()
for (selection_var in selection_cat_vars) {
  selection_x <- if (selection_var %in% names(selection_data)) {
    trimws(as.character(selection_data[[selection_var]]))
  } else rep(NA_character_, nrow(selection_data))
  selection_x[is.na(selection_x) | selection_x %in% c("", "NA", "N/A", "n/a")] <- NA_character_
  selection_note <- "Source categories; missing values excluded from Fisher test"
  selection_original_missing <- is.na(selection_x)
  # Free-text medication descriptions are summarized by recording availability.
  if (selection_var == "Pharm_treatment" && length(unique(stats::na.omit(selection_x))) > 12L) {
    selection_x <- ifelse(selection_original_missing, "Not recorded", "Recorded")
    selection_note <- "Complex baseline medication text summarized as recorded/not recorded; not postoperative medication use"
  }
  selection_levels <- sort(unique(stats::na.omit(selection_x)))
  if (any(is.na(selection_x))) selection_levels <- c(selection_levels, "(Missing)")
  for (selection_flag in c(FALSE, TRUE)) {
    selection_g <- selection_x[selection_data$followup_available == selection_flag]
    selection_orig_missing_g <- selection_original_missing[selection_data$followup_available == selection_flag]
    for (selection_level in selection_levels) {
      selection_n <- if (selection_level == "(Missing)") sum(is.na(selection_g)) else sum(selection_g == selection_level, na.rm = TRUE)
      selection_cat_rows[[length(selection_cat_rows) + 1L]] <- data.frame(
        variable = selection_var, category = selection_level,
        followup_status = selection_labels[as.integer(selection_flag) + 1L],
        n = selection_n, n_group = length(selection_g), n_missing_source = sum(selection_orig_missing_g),
        percent_group = 100 * selection_n/length(selection_g),
        percent_observed = if (selection_level == "(Missing)" || !any(!is.na(selection_g))) NA_real_ else
          100 * selection_n/sum(!is.na(selection_g)),
        note = selection_note)
    }
  }
  selection_table <- table(selection_data$followup_available, selection_x)
  selection_fisher <- if (all(dim(selection_table) >= 2L)) {
    tryCatch(stats::fisher.test(selection_table, alternative = "two.sided", workspace = 20000000L),
             error = function(e) NULL)
  } else NULL
  selection_test_rows[[length(selection_test_rows) + 1L]] <- data.frame(
    variable = selection_var, variable_type = "categorical", test = "Fisher exact test (two-sided)",
    n_available = sum(!is.na(selection_x[selection_data$followup_available])),
    n_no_paired_followup = sum(!is.na(selection_x[!selection_data$followup_available])),
    mean_difference = NA_real_, ci_low = NA_real_, ci_high = NA_real_, standardized_mean_difference = NA_real_,
    p_value_unadjusted = if (is.null(selection_fisher)) NA_real_ else selection_fisher$p.value,
    interpretation = "Descriptive only; unadjusted P value; no covariate selection",
    note = paste(selection_note, if (is.null(selection_fisher)) "Test not estimable" else "", sep = "; "))
}
store_result(do.call(rbind, selection_cat_rows), "selection_categorical")
store_result(do.call(rbind, selection_test_rows), "selection_comparisons")
# Actual baseline exposure is distinct from availability of medication text.
medication <- trimws(as.character(selection_data$Pharm_treatment))
medication_recorded <- !is.na(medication) &
  !medication %in% c("", "NA", "N/A", "n/a")
medication_exposures <- list(
  "Any pharmacological treatment" = medication_recorded & tolower(medication) != "no",
  "Carbamazepine-containing regimen" = medication_recorded &
    grepl("(^|[^A-Za-z])CBZ([^A-Za-z]|$)", medication)
)
medication_rows <- lapply(names(medication_exposures), function(label) {
  exposed <- medication_exposures[[label]]
  paired_n <- sum(exposed & selection_data$followup_available)
  remaining_n <- sum(exposed & !selection_data$followup_available)
  paired_denominator <- sum(medication_recorded & selection_data$followup_available)
  remaining_denominator <- sum(medication_recorded & !selection_data$followup_available)
  table <- matrix(c(paired_n, paired_denominator - paired_n,
                    remaining_n, remaining_denominator - remaining_n), nrow = 2L)
  test <- if (all(rowSums(table) > 0) && all(colSums(table) > 0))
    stats::fisher.test(table) else NULL
  data.frame(characteristic = label, paired_n = paired_n,
    paired_denominator = paired_denominator, remaining_n = remaining_n,
    remaining_denominator = remaining_denominator,
    p_value = if (is.null(test)) NA_real_ else test$p.value,
    paired_missing = sum(!medication_recorded & selection_data$followup_available),
    remaining_missing = sum(!medication_recorded & !selection_data$followup_available),
    interpretation = "Baseline exposure; observed denominators; descriptive unadjusted Fisher P value")
})
store_result(dplyr::bind_rows(medication_rows), "selection_medication_exposure")

selection_results <- results
