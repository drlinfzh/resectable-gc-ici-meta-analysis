# Fig. 4. Survival by chemotherapy backbone.
library(survival)
library(grid)

input_dir <- "inputs"

output_dir <- "outputs/fig4"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

ink <- "#2C2E35"
blue <- "#0871B9"
orange <- "#FCAF17"
ci_grey <- "#959699"
row_grey <- "#B5B5B5"
font_family <- "Arial"

label <- function(text, x, y, size = 10, just = "centre", bold = FALSE,
                  rot = 0, colour = ink) {
  grid.text(text, unit(x, "mm"), unit(y, "mm"), just = just, rot = rot,
            gp = gpar(fontfamily = font_family, fontsize = size, col = colour,
                      fontface = if (bold) "bold" else "plain", lineheight = 1))
}

line <- function(x0, y0, x1, y1, col = ink, lwd = .7, lty = "solid") {
  grid.segments(unit(x0, "mm"), unit(y0, "mm"), unit(x1, "mm"), unit(y1, "mm"),
                gp = gpar(col = col, lwd = lwd, lty = lty, lineend = "butt"))
}

p_label <- function(p) {
  if (is.na(p)) return("P not reported")
  if (p < .001) "P < 0.001" else sprintf("P = %.3f", p)
}

ci_label <- function(lo, hi, digits = 2L) {
  if (is.na(lo) || is.na(hi)) return("95% CI not reported")
  sprintf("95%% CI %.*f–%.*f", digits, lo, digits, hi)
}

endpoint_name <- function(endpoint) {
  if (endpoint == "EFS") "Event-free survival" else "Overall survival"
}

pdf_start <- function(filename, width, height) {
  grDevices::cairo_pdf(file.path(output_dir, filename), width = width / 25.4,
                       height = height / 25.4, family = font_family)
  grid.newpage()
}

cox_estimate <- function(data, stratified = FALSE) {
  model <- if (stratified) {
    coxph(Surv(time, status) ~ exposure + strata(study), data = data, ties = "efron")
  } else {
    coxph(Surv(time, status) ~ exposure, data = data, ties = "efron")
  }
  s <- summary(model)
  data.frame(logHR = s$coefficients[1, "coef"], SE = s$coefficients[1, "se(coef)"],
             HR = s$conf.int[1, "exp(coef)"], lower95 = s$conf.int[1, "lower .95"],
             upper95 = s$conf.int[1, "upper .95"], P = s$coefficients[1, "Pr(>|z|)"])
}

pool_estimates <- function(studies) {
  valid <- is.finite(studies$logHR) & is.finite(studies$SE) & studies$SE > 0
  stopifnot(any(valid))
  fit <- metafor::rma.uni(yi = logHR, sei = SE, data = studies[valid, ],
                          method = if (sum(valid) > 1L) "REML" else "FE", test = "z")
  studies$weight_percent <- NA_real_
  studies$weight_percent[valid] <- as.numeric(weights(fit))
  pooled <- data.frame(HR = exp(fit$b[1]), lower95 = exp(fit$ci.lb), upper95 = exp(fit$ci.ub),
                       P = fit$pval, I2_percent = fit$I2, tau2 = fit$tau2, heterogeneity_P = fit$QEp)
  list(studies = studies, pooled = pooled)
}

trial_estimates <- function(data, trial_names) {
  ids <- names(trial_names)[names(trial_names) %in% data$study]
  do.call(rbind, lapply(ids, function(study_id) {
    trial <- data[data$study == study_id, ]
    e1 <- sum(trial$status[trial$exposure == 1L])
    e0 <- sum(trial$status[trial$exposure == 0L])
    result <- if (e1 > 0 && e0 > 0) cox_estimate(trial) else
      data.frame(logHR = NA_real_, SE = NA_real_, HR = NA_real_,
                 lower95 = NA_real_, upper95 = NA_real_, P = NA_real_)
    data.frame(study = trial_names[[study_id]], e1 = e1, n1 = sum(trial$exposure == 1L),
               e0 = e0, n0 = sum(trial$exposure == 0L), result)
  }))
}

trial_input <- function(endpoint) {
  d <- as.data.frame(readxl::read_excel(
    file.path(input_dir, "ipd", paste0(tolower(endpoint), "_ipd.xlsx"))))
  d$exposure <- as.integer(d$group == "ici_chemo")
  stopifnot(all(d$group %in% c("ici_chemo", "chemo")),
            !anyNA(d[c("study", "time", "status", "exposure")]))
  d[c("study", "time", "status", "exposure")]
}

draw_km <- function(data, estimate, filename, endpoint, title = "",
                    arms = c("ICI plus chemo", "Chemo alone"), reported = NULL,
                    citation = NULL, ci_digits = 2L, xmax = NULL) {
  stopifnot(!anyNA(data[c("time", "status", "exposure")]), all(data$time >= 0),
            all(data$status %in% 0:1), setequal(data$exposure, 0:1))
  data$arm <- factor(data$exposure, levels = c(1L, 0L))
  fit <- survfit(Surv(time, status) ~ arm, data = data)
  curve <- summary(fit, censored = TRUE, data.frame = TRUE)
  if (is.null(xmax)) xmax <- ceiling(max(data$time) / 6) * 6
  ticks <- seq(0, xmax, by = 6)
  risk <- summary(fit, times = ticks, extend = TRUE, data.frame = TRUE)
  stopifnot(sum(curve$n.event) == sum(data$status), xmax >= max(data$time))

  pdf_start(filename, width = 215, height = 130)
  left <- 42
  right <- 208
  bottom <- 41
  x <- function(months) left + (right - left) * months / xmax
  y <- function(probability) bottom + 72 * probability
  panel <- substr(filename, 1, 1)
  label(panel, 3, 124, size = 14, just = "left", bold = TRUE)
  label(title, left, 124, size = 11, just = "left", bold = TRUE)
  if (!is.null(citation)) label(citation, left, 119, size = 9, just = "left")

  for (probability in seq(.2, 1, .2))
    line(left, y(probability), right, y(probability), col = row_grey, lwd = .5)
  line(left, bottom, left, y(1.04), col = "#000000", lwd = 1)
  line(left, bottom, right, bottom, col = "#000000", lwd = 1)
  probabilities <- seq(0, 1, .2)
  line(left - 1.2, y(probabilities), left, y(probabilities), col = "#000000")
  label(round(100 * probabilities), left - 2.5, y(probabilities), just = "right")
  line(x(ticks), bottom, x(ticks), bottom - 1.4, col = "#000000")
  label(ticks, x(ticks), bottom - 4.5, size = 9)
  label(paste0(endpoint_name(endpoint), " (%)"), 15, y(.5), size = 11, rot = 90, bold = TRUE)
  label("Time (months)", (left + right) / 2, bottom - 10, size = 11, bold = TRUE)

  colours <- c(blue, orange)
  strata_names <- c("arm=1", "arm=0")
  at_time <- function(curve, t) c(1, curve$surv)[findInterval(t, curve$time) + 1L]
  curves <- lapply(strata_names, function(key) curve[curve$strata == key, ])
  shared_end <- min(vapply(curves, function(z) max(z$time), numeric(1)))
  anchor <- min(.62 * xmax, .8 * shared_end)
  survival_over <- function(from, to) {
    start <- xmax * (from - left) / (right - left)
    end <- xmax * (to - left) / (right - left)
    unlist(lapply(curves, function(z)
      c(at_time(z, start), at_time(z, end), z$surv[z$time >= start & z$time <= end])))
  }
  at_anchor <- vapply(curves, at_time, numeric(1), t = anchor)
  offsets <- if (at_anchor[1] >= at_anchor[2]) c(.052, -.052) else c(-.052, .052)

  label("Number at risk", 3, 22, just = "left", bold = TRUE)
  risk_labels <- sub(" complete response", "\ncomplete response", arms, fixed = TRUE)
  for (i in seq_along(arms)) {
    z <- curves[[i]]
    times <- c(0, z$time)
    survival <- c(1, z$surv)
    step_time <- rep(times, each = 2)[-1]
    step_survival <- head(rep(survival, each = 2), -1)
    grid.lines(unit(x(step_time), "mm"), unit(y(step_survival), "mm"),
               gp = gpar(col = colours[i], lwd = 1.2, lineend = "butt"))
    censored <- z$n.censor > 0
    line(x(z$time[censored]), y(z$surv[censored]), x(z$time[censored]),
         y(z$surv[censored]) + .9, col = colours[i], lwd = .55)

    text_width <- convertWidth(grobWidth(textGrob(arms[i],
      gp = gpar(fontfamily = font_family, fontsize = 9.5))), "mm", valueOnly = TRUE)
    text_left <- max(left + 3, min(x(anchor), x(shared_end) - text_width - 2))
    under_text <- survival_over(text_left, text_left + text_width)
    text_y <- if (offsets[i] > 0) max(under_text) + offsets[i] else min(under_text) + offsets[i]
    label(arms[i], text_left, y(min(1.035, max(.035, text_y))),
          size = 9.5, colour = colours[i], just = "left")

    risk_y <- c(15, 6.5)[i]
    label(risk_labels[i], 3, risk_y, size = 8.7, just = "left")
    n_at_risk <- risk$n.risk[risk$strata == strata_names[i]]
    label(n_at_risk, x(ticks), risk_y, size = 8.7)
  }

  values <- list(estimate)
  stat_x <- left + 4
  top <- y(.38)
  if (!is.null(reported)) {
    values <- list(reported, estimate)
    stat_x <- left + c(4, 64)
    above_stats <- min(survival_over(stat_x[1], stat_x[2] + 36))
    top <- min(top, y(above_stats) - 9)
    label(c("Original study", "Reconstructed data"), stat_x, top + 5.5,
          just = "left", bold = TRUE)
  }
  for (i in seq_along(values)) {
    s <- values[[i]]
    hr <- if (is.na(s$HR)) "HR not reported" else sprintf("HR = %.2f", s$HR)
    p <- if (identical(s$P_relation, "<") && !is.na(s$P)) sprintf("P < %.3f", s$P) else p_label(s$P)
    lines <- c(hr, ci_label(s$lower95, s$upper95, ci_digits), p)
    label(lines, stat_x[i], top - c(0, 4.5, 9), size = 9.5, just = "left")
  }
  invisible(dev.off())
}

draw_forest <- function(blocks, filename, title, ticks, arms = c("ICI plus chemo", "Chemo alone"),
                        metric = "HR", combined_counts = FALSE, favors_above = FALSE,
                        favors = arms, note = NULL, sparse_ci = NULL) {
  group_count <- length(blocks)
  group_rows <- vapply(blocks, function(block) nrow(block$studies), integer(1))
  group_headings <- nzchar(names(blocks))
  row_height <- 6
  axis_y <- if (favors_above) 16 else 27
  height <- axis_y + 30 + sum((group_rows + 4 + group_headings) * row_height) +
    if (is.null(note)) 0 else row_height
  width <- 208
  table_end <- 138
  graph_left <- 146
  graph_right <- 187
  weight_x <- 199
  table_y <- height - 16
  count_x <- if (combined_counts) c(56, 72) else c(53, 68, 83, 98)
  effect_x <- if (combined_counts) 86 else 106
  ci <- unlist(lapply(blocks, function(block)
    c(block$studies$lower95, block$studies$upper95, block$pooled$lower95, block$pooled$upper95)))
  stopifnot(all(ci[is.finite(ci)] > 0), all(ticks > 0))
  limits <- range(c(ticks, ci), na.rm = TRUE)
  x <- function(hr) graph_left + (graph_right - graph_left) *
    (log(hr) - log(limits[1])) / diff(log(limits))

  pdf_start(filename, width, height)
  panel <- substr(filename, 1, 1)
  label(panel, 3, height - 5, size = 14, just = "left", bold = TRUE)
  label(title, 9, height - 5, size = 11, just = "left", bold = TRUE)
  label("Trial", 3, table_y, just = "left", bold = TRUE)
  label(paste0(metric, " (95% CI)"), effect_x, table_y, just = "left", bold = TRUE)
  label("Weight, %", weight_x, table_y, bold = TRUE)
  if (combined_counts) {
    label(c("Events", "Total"), count_x, table_y, bold = TRUE)
  } else {
    group_label <- sub(" complete response", "\ncomplete response", arms, fixed = TRUE)
    label(group_label, c(mean(count_x[1:2]), mean(count_x[3:4])), table_y + 4.5,
          size = 8.5, bold = TRUE)
    line(count_x[c(1, 3)] - 6, table_y, count_x[c(2, 4)] + 6, table_y, lwd = .95)
    count_name <- if (metric == "RR") "Responses" else "Events"
    label(rep(c(count_name, "Total"), 2), count_x, table_y - 4.5, size = 9, bold = TRUE)
  }

  guide_x <- c(graph_left + 5, graph_right - 2)
  if (favors_above) {
    label("Favors", x(1), table_y + 5, bold = TRUE)
    label(c("ICI plus\nchemo", "Chemo\nalone"), x(1) + c(-15, 13), table_y - 2,
          size = 8.5, bold = TRUE)
  }
  line(3, table_y - 8, table_end, table_y - 8, lwd = 1.2)
  line(x(1), axis_y, x(1), table_y - 8, lwd = .85, lty = "24")

  y <- table_y - 11
  for (i in seq_len(group_count)) {
    studies <- blocks[[i]]$studies
    pooled <- blocks[[i]]$pooled
    if (group_headings[i]) {
      label(names(blocks)[i], 3, y, just = "left", bold = TRUE)
      y <- y - row_height
    }
    rows <- y - row_height * (seq_len(nrow(studies)) - 1L)
    valid <- is.finite(studies$HR)
    label(studies$study, 3, rows, just = "left", size = 9.2)
    if (combined_counts) {
      label(ifelse(is.na(studies$events), "NR", studies$events), count_x[1], rows)
      label(studies$n, count_x[2], rows)
    } else {
      counts <- if (metric == "RR") c("e0", "n0", "e1", "n1") else c("e1", "n1", "e0", "n0")
      for (j in seq_along(counts)) label(studies[[counts[j]]], count_x[j], rows)
    }
    digits <- if (is.null(sparse_ci)) 2L else
      ifelse(studies$study == sparse_ci[1] & names(blocks)[i] == sparse_ci[2], 3L, 2L)
    effects <- ifelse(valid, sprintf("%.2f (%.*f–%.*f)", studies$HR,
                     digits, studies$lower95, digits, studies$upper95), "Not estimable")
    label(effects, effect_x, rows, size = 9.2, just = "left")
    label(ifelse(valid, sprintf("%.1f", studies$weight_percent), "—"), weight_x, rows)
    line(x(studies$lower95[valid]), rows[valid], x(studies$upper95[valid]), rows[valid],
         col = ci_grey, lwd = .9)
    squares <- 4 * sqrt(studies$weight_percent[valid] / 100)
    grid.rect(unit(x(studies$HR[valid]), "mm"), unit(rows[valid], "mm"),
              width = unit(squares, "mm"), height = unit(squares, "mm"),
              gp = gpar(fill = "#398DC7", col = ink, lwd = .6))
    for (rule_y in rows - row_height / 2) line(3, rule_y, table_end, rule_y, col = row_grey)

    y <- min(rows) - row_height
    summary_label <- if (sum(valid) == 1L) "Single-trial estimate" else "Total (random effects)"
    label(summary_label, 3, y, just = "left", bold = TRUE, size = 9.2)
    if (combined_counts) {
      total_events <- if (anyNA(studies$events[valid])) "NA" else sum(studies$events[valid])
      label(c(total_events, sum(studies$n[valid])), count_x, y, bold = TRUE)
    } else {
      totals <- colSums(studies[valid, counts, drop = FALSE])
      label(totals, count_x, y, bold = TRUE)
    }
    estimate <- sprintf("%.2f (%.2f–%.2f)", pooled$HR, pooled$lower95, pooled$upper95)
    label(estimate, effect_x, y, size = 9.2, just = "left", bold = TRUE)
    label("100.0", weight_x, y, bold = TRUE)
    if (sum(valid) > 1L) {
      grid.polygon(unit(x(c(pooled$lower95, pooled$HR, pooled$upper95, pooled$HR)), "mm"),
                   unit(y + c(0, 1.6, 0, -1.6), "mm"), gp = gpar(fill = blue, col = ink))
      label(sprintf("Heterogeneity: I² = %.1f%%; τ² = %.4f; %s", pooled$I2_percent,
                    pooled$tau2, p_label(pooled$heterogeneity_P)),
            3, y - row_height, size = 9.2, just = "left")
    }
    label(paste("Overall effect:", p_label(pooled$P)), 3, y - 2 * row_height,
          size = 9.2, just = "left")
    y <- y - 3 * row_height
  }

  if (!is.null(note)) label(note, 3, y, just = "left", size = 9.2)
  line(graph_left, axis_y, graph_right, axis_y, lwd = 1)
  line(x(ticks), axis_y, x(ticks), axis_y - 1.4, lwd = 1)
  tick_text <- vapply(ticks, format, character(1), trim = TRUE, scientific = FALSE)
  label(tick_text, x(ticks), axis_y - 4.5, size = 9)
  label(paste0(metric, " (95% CI)"), mean(guide_x), axis_y - 10, bold = TRUE)
  if (!favors_above) {
    group_label <- sub(" complete response", "\ncomplete response", favors, fixed = TRUE)
    label(group_label, guide_x + c(-2, 5), axis_y - 16, size = 8.5, bold = TRUE)
    label("Favors", mean(guide_x), axis_y - 22, size = 9, bold = TRUE)
  }
  invisible(dev.off())
}

trial_names <- c(matterhorn = "MATTERHORN", kn585_flot = "KEYNOTE-585 FLOT",
                 kn585_main = "KEYNOTE-585 main", ns01 = "NEOSUMMIT-01", astrum006 = "ASTRUM-006")

for (endpoint in c("EFS", "OS")) {
  overall <- trial_input(endpoint)
  data <- overall[overall$study != "kn585", ]
  main <- read.csv(file.path(input_dir, "ipd",
                  paste0("kn585_", tolower(endpoint), "_main_ipd.csv")))
  flot <- read.csv(file.path(input_dir, "kmsubtraction",
                  paste0("kn585_flot_", tolower(endpoint), "_kmsubtraction_ipd.csv")))
  main$study <- "kn585_main"
  flot$study <- "kn585_flot"
  main$exposure <- main$treat
  flot$exposure <- flot$treat
  data <- rbind(data, main[names(data)], flot[names(data)])
  stopifnot(nrow(data) == nrow(overall))
  data$nonflot <- as.integer(!data$study %in% c("matterhorn", "kn585_flot"))

  for (regimen in 0:1) {
    subgroup <- data[data$nonflot == regimen, ]
    title <- c("FLOT backbone", "Non-FLOT backbone")[regimen + 1L]
    panel <- if (endpoint == "EFS") c("a", "c", "e", "g") else c("b", "d", "f", "h")
    one_stage <- cox_estimate(subgroup, stratified = TRUE)
    two_stage <- pool_estimates(trial_estimates(subgroup, trial_names))
    draw_km(subgroup, one_stage, paste0(panel[2L * regimen + 1L], "_km.pdf"), endpoint, title)
    draw_forest(setNames(list(two_stage), ""), paste0(panel[2L * regimen + 2L], "_forest.pdf"),
                paste(title, tolower(endpoint_name(endpoint)), sep = ", "), ticks = c(.2, .5, 1, 2))
  }
}
