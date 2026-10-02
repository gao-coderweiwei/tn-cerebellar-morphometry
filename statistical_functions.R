# Shared statistical functions. Dependencies: dplyr, sandwich, lmtest, car.
# Results are stored in memory; no files are read or written.
store_result <- function(x, name) {
  results[[name]] <<- as.data.frame(x)
  invisible(x)
}

require_full_rank <- function(fit) {
  if (fit$rank != ncol(model.matrix(fit))) stop("Rank-deficient regression")
  if (df.residual(fit) <= 0) stop("No residual degrees of freedom")
  invisible(fit)
}
hc3_table <- function(fit) {
  require_full_rank(fit)
  cov <- sandwich::vcovHC(fit, type = "HC3")
  tab <- lmtest::coeftest(fit, vcov. = cov, df = df.residual(fit))
  critical <- qt(0.975, df.residual(fit))
  data.frame(term = rownames(tab), estimate = tab[, 1], std_error = tab[, 2],
             conf_low = tab[, 1] - critical * tab[, 2],
             conf_high = tab[, 1] + critical * tab[, 2], statistic = tab[, 3],
             df = df.residual(fit), p_value = tab[, 4], row.names = NULL)
}
fit_group <- function(data, outcome, metric, age_model = "quadratic", icv = FALSE, qc = FALSE) {
  d <- data
  d$y <- d[[outcome]]
  age_terms <- switch(age_model, linear = "age_c", quadratic = c("age_c", "age_c2"),
                      spline = c("age_c", "age_rcs_nonlinear"), stop("Unknown age model"))
  predictors <- c("group", age_terms, "sex",
                  if (metric == "volume" || icv) "icv_c", if (qc) "QC")
  complete <- complete.cases(d[, c("y", predictors)])
  d <- d[complete, , drop = FALSE]
  fit <- lm(reformulate(predictors, response = "y"), data = d, na.action = na.fail)
  require_full_rank(fit)
  list(fit = fit, data = d, n_missing = sum(!complete))
}
result_group <- function(object, info, age_model) {
  fit <- object$fit; d <- object$data
  estimate <- hc3_table(fit)
  estimate <- estimate[estimate$term == "groupTN", setdiff(names(estimate), "term"), drop = FALSE]
  stopifnot(nrow(estimate) == 1)
  data.frame(roi_order = info$roi_order, roi = info$roi, metric = info$metric,
             unit = ifelse(info$metric == "volume", "cm3", "mm"), age_model = age_model,
             source_variable = info$outcome, n = nrow(d), n_missing = object$n_missing,
             n_hc = sum(d$group == "HC"), n_tn = sum(d$group == "TN"),
             hc_mean = mean(d$y[d$group == "HC"]), hc_sd = sd(d$y[d$group == "HC"]),
             tn_mean = mean(d$y[d$group == "TN"]), tn_sd = sd(d$y[d$group == "TN"]),
             estimate, aic = AIC(fit), row.names = NULL)
}
bh <- function(x, groups, expected_n) {
  x |>
    dplyr::group_by(dplyr::across(dplyr::all_of(groups))) |>
    dplyr::group_modify(function(.x, .y) {
      if (nrow(.x) != expected_n || any(!is.finite(.x$p_value)) || any(.x$p_value < 0 | .x$p_value > 1)) stop("Incomplete FDR family")
      .x$q_value <- p.adjust(.x$p_value, method = "BH")
      .x$significant <- .x$q_value < 0.05
      .x$fdr_family_size <- expected_n
      .x
    }) |> dplyr::ungroup()
}
continuous_summary <- function(d, variables) {
  do.call(rbind, lapply(variables, function(v) do.call(rbind, lapply(levels(d$group), function(g) {
    x <- d[[v]][d$group == g]; good <- x[!is.na(x)]
    data.frame(variable = v, group = g, n = length(good), missing = sum(is.na(x)),
               mean = mean(good), sd = sd(good), median = median(good),
               q25 = unname(quantile(good, .25)), q75 = unname(quantile(good, .75)),
               min = min(good), max = max(good))
  })) ))
}
diagnostic_tables <- function(fit, ids) {
  n <- nobs(fit); p <- fit$rank
  details <- data.frame(participant_id = ids, fitted = fitted(fit), residual = residuals(fit),
                        studentized = rstudent(fit), cooks_distance = cooks.distance(fit),
                        leverage = hatvalues(fit))
  details$flag_studentized <- abs(details$studentized) > 3
  details$flag_cook <- details$cooks_distance > 4/n
  details$flag_leverage <- details$leverage > 2*p/n
  summary <- data.frame(n = n, parameters = p, cook_threshold = 4/n,
                         leverage_threshold = 2*p/n, n_cook = sum(details$flag_cook),
                         n_studentized = sum(details$flag_studentized),
                         n_leverage = sum(details$flag_leverage),
                         shapiro_p = shapiro.test(residuals(fit))$p.value,
                         breusch_pagan_p = lmtest::bptest(fit)$p.value)
  vif <- car::vif(fit)
  if (is.matrix(vif)) {
    vif <- data.frame(term = rownames(vif), gvif = vif[,1], term_df = vif[,2],
                      gvif_adjusted = vif[,3], row.names = NULL)
  } else vif <- data.frame(term = names(vif), gvif = as.numeric(vif), term_df = 1,
                            gvif_adjusted = sqrt(as.numeric(vif)))
  list(summary = summary, details = details, vif = vif)
}
