# Supplementary Fig. 9. Event-free survival across clinical subgroups.
library(survival)
library(grid)

input_dir <- "inputs"

output_dir <- "outputs/sfig9"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

ink <- "#2C2E35"
blue <- "#0871B9"
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

endpoint <- "EFS"
axis_ticks <- list(c(.01, .1, 1, 2), c(.1, 1, 2, 5), c(.1, .5, 1, 2), c(.1, 1, 2, 5))
extended_precision <- c("NEOSUMMIT-01", "Diffuse")

published <- read.csv(text = 'feature,subgroup,study,n,HR,lower95,upper95
Age,Age <65 years,NEOSUMMIT-01,70,0.76,0.36,1.57
Age,Age ≥65 years,NEOSUMMIT-01,38,0.11,0.01,0.89
Sex,Male,NEOSUMMIT-01,70,0.62,0.27,1.4
Sex,Female,NEOSUMMIT-01,38,0.45,0.16,1.33
Tumor location,Gastric,NEOSUMMIT-01,71,0.55,0.25,1.21
Tumor location,Gastroesophageal junction,NEOSUMMIT-01,37,0.4,0.12,1.33
Age,Age <65 years,KEYNOTE-585,550,0.7,0.55,0.89
Age,Age ≥65 years,KEYNOTE-585,457,0.93,0.72,1.2
Sex,Female,KEYNOTE-585,283,0.73,0.52,1.01
Sex,Male,KEYNOTE-585,724,0.83,0.67,1.01
Tumor location,Gastric,KEYNOTE-585,762,0.76,0.62,0.93
Tumor location,Gastroesophageal junction,KEYNOTE-585,244,0.91,0.65,1.28
Sex,Male,MATTERHORN,682,0.65,0.51,0.83
Sex,Female,MATTERHORN,266,0.84,0.58,1.24
Age,Age <65 years,MATTERHORN,556,0.71,0.54,0.92
Age,Age ≥65 years,MATTERHORN,392,0.7,0.51,0.97
Tumor location,Gastric,MATTERHORN,640,0.76,0.59,0.97
Tumor location,Gastroesophageal junction,MATTERHORN,308,0.61,0.43,0.86
Age,Age <65 years,ASTRUM-006,400,0.72,0.53,0.99
Age,Age ≥65 years,ASTRUM-006,188,0.8,0.51,1.27
Sex,Female,ASTRUM-006,124,1.49,0.8,2.76
Sex,Male,ASTRUM-006,464,0.64,0.48,0.85
Lauren,Intestinal,KEYNOTE-585,468,0.82,0.63,1.08
Lauren,Diffuse,KEYNOTE-585,424,0.76,0.59,0.99
Lauren,Mixed,KEYNOTE-585,90,0.88,0.5,1.54
Lauren,Intestinal,MATTERHORN,483,0.66,0.48,0.89
Lauren,Diffuse,MATTERHORN,249,0.93,0.66,1.32
Lauren,Intestinal,ASTRUM-006,48,0.59,0.22,1.57
Lauren,Diffuse,ASTRUM-006,55,0.6,0.27,1.35')

event_counts <- read.csv(file.path(input_dir, "clinical", "subgroup_events.csv"),
                         encoding = "UTF-8", na.strings = c("NA", "NR"))
counts <- event_counts[event_counts$endpoint == endpoint & event_counts$trial != "Total", ]
counts$key <- with(counts, paste(feature, subgroup, trial))
published$key <- with(published, paste(feature, subgroup, study))
indices <- match(published$key, counts$key)
stopifnot(!anyNA(indices), all(published$n == counts$total[indices]))
published$events <- counts$events[indices]
published$logHR <- log(published$HR)
published$SE <- (log(published$upper95) - log(published$lower95)) / (2 * qnorm(.975))
published$key <- NULL

trial <- as.data.frame(readxl::read_excel(file.path(input_dir, "ipd", "clinical_cohort.xlsx"), sheet = "Sheet1"))
histology <- as.data.frame(readxl::read_excel(file.path(input_dir, "ipd", "histology_source.xlsx"),
                                              sheet = "NEOSUMMIT-01"))
record_id <- trimws(as.character(trial[["病历号"]]))
source_id <- trimws(as.character(histology[["病历号"]]))
rows <- match(record_id, source_id)
stopifnot(nrow(trial) == 108L, !anyDuplicated(record_id), !anyDuplicated(source_id),
          !anyNA(rows), all(trial$Group == histology$Group[rows]))
for (field in c("EFS_m", "EFS_status", "OS_m", "OS_status"))
  stopifnot(all(abs(trial[[field]] - histology[[field]][rows]) < 1e-9))
trial$subgroup <- ifelse(histology$Lauren_type[rows] == "Diffuse", "Diffuse",
                         trial$Lauren_type_cat3)
stopifnot(!anyNA(trial$subgroup),
          identical(as.integer(table(trial$subgroup)[c("Intestinal", "Diffuse", "Mixed")]),
                    c(44L, 38L, 26L)))
trial$exposure <- as.integer(trial$Group == "Toripalimab + chemotherapy")
trial$time <- trial[[paste0(endpoint, "_m")]]
trial$status <- trial[[paste0(endpoint, "_status")]]

lauren <- do.call(rbind, lapply(c("Intestinal", "Diffuse", "Mixed"), function(group) {
  data <- trial[trial$subgroup == group, ]
  estimate <- cox_estimate(data)
  data.frame(feature = "Lauren", subgroup = group, study = "NEOSUMMIT-01",
             n = nrow(data), HR = estimate$HR, lower95 = estimate$lower95,
             upper95 = estimate$upper95, events = sum(data$status),
             logHR = estimate$logHR, SE = estimate$SE)
}))
estimates <- rbind(published, lauren)
specs <- list(Age = c("Age <65 years", "Age ≥65 years"), Sex = c("Male", "Female"),
              "Tumor location" = c("Gastric", "Gastroesophageal junction"),
              Lauren = c("Intestinal", "Diffuse", "Mixed"))

for (i in seq_along(specs)) {
  feature <- names(specs)[i]
  subgroups <- specs[[i]]
  studies <- if (feature == "Lauren") {
    c("KEYNOTE-585", "MATTERHORN", "ASTRUM-006", "NEOSUMMIT-01")
  } else {
    c("NEOSUMMIT-01", "KEYNOTE-585", "MATTERHORN", "ASTRUM-006")
  }
  blocks <- lapply(setNames(subgroups, subgroups), function(group) {
    data <- estimates[estimates$feature == feature & estimates$subgroup == group, ]
    data <- data[order(match(data$study, studies)), ]
    pool_estimates(data)
  })
  first <- blocks[[1]]$studies
  second <- blocks[[2]]$studies
  stopifnot(identical(first$study, second$study))
  interaction <- metafor::rma.uni(yi = second$logHR - first$logHR,
                                   sei = sqrt(first$SE^2 + second$SE^2),
                                   method = "REML", test = "z")
  note <- paste("Interaction effect:", p_label(interaction$pval))
  if (feature == "Lauren") note <- paste(note, "(diffuse vs intestinal)")
  chart_title <- if (feature == "Lauren") "Lauren classification" else feature
  draw_forest(blocks, paste0(letters[i], "_forest.pdf"), chart_title, ticks = axis_ticks[[i]],
              combined_counts = TRUE, note = note, sparse_ci = extended_precision)
}
