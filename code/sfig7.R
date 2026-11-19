# Supplementary Fig. 7. Meta-analysis of pathological complete response.
library(grid)

output_dir <- "outputs/sfig7"
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

pool_estimates <- function(studies) {
  valid <- is.finite(studies$logRR) & is.finite(studies$SE) & studies$SE > 0
  stopifnot(any(valid))
  fit <- metafor::rma.uni(yi = logRR, sei = SE, data = studies[valid, ],
                          method = if (sum(valid) > 1L) "REML" else "FE", test = "z")
  studies$weight_percent <- NA_real_
  studies$weight_percent[valid] <- as.numeric(weights(fit))
  pooled <- data.frame(RR = exp(fit$b[1]), lower95 = exp(fit$ci.lb), upper95 = exp(fit$ci.ub),
                       P = fit$pval, I2_percent = fit$I2, tau2 = fit$tau2, heterogeneity_P = fit$QEp)
  list(studies = studies, pooled = pooled)
}

draw_forest <- function(blocks, filename, title, ticks, arms = c("ICI plus chemo", "Chemo alone"),
                        metric = "RR", combined_counts = FALSE, favors_above = FALSE,
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
    valid <- is.finite(studies$RR)
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
    effects <- ifelse(valid, sprintf("%.2f (%.*f–%.*f)", studies$RR,
                     digits, studies$lower95, digits, studies$upper95), "Not estimable")
    label(effects, effect_x, rows, size = 9.2, just = "left")
    label(ifelse(valid, sprintf("%.1f", studies$weight_percent), "—"), weight_x, rows)
    line(x(studies$lower95[valid]), rows[valid], x(studies$upper95[valid]), rows[valid],
         col = ci_grey, lwd = .9)
    squares <- 4 * sqrt(studies$weight_percent[valid] / 100)
    grid.rect(unit(x(studies$RR[valid]), "mm"), unit(rows[valid], "mm"),
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
    estimate <- sprintf("%.2f (%.2f–%.2f)", pooled$RR, pooled$lower95, pooled$upper95)
    label(estimate, effect_x, y, size = 9.2, just = "left", bold = TRUE)
    label("100.0", weight_x, y, bold = TRUE)
    if (sum(valid) > 1L) {
      grid.polygon(unit(x(c(pooled$lower95, pooled$RR, pooled$upper95, pooled$RR)), "mm"),
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

# Published counts. The KEYNOTE-585 FLOT cohort is the combined cohort minus the main cohort.
trials <- data.frame(
  study = c("MATTERHORN", "KEYNOTE-585", "KEYNOTE-585 FLOT", "KEYNOTE-585 main",
            "NEOSUMMIT-01", "ASTRUM-006"),
  e1 = c(91, 70, 70 - 54, 54, 12, 63), n1 = c(474, 492, 492 - 402, 402, 54, 292),
  e0 = c(34, 14, 14 - 8, 8, 4, 19), n0 = c(474, 495, 495 - 402, 402, 54, 296)
)
groups <- list(Overall = c(1, 2, 5, 6), "FLOT-based" = c(1, 3), "Non-FLOT-based" = c(4, 5, 6))
stopifnot(all(trials$e1 > 0 & trials$e1 <= trials$n1),
          all(trials$e0 > 0 & trials$e0 <= trials$n0))

estimates <- metafor::escalc(measure = "RR", ai = e1, n1i = n1, ci = e0, n2i = n0,
                             data = trials, add = 0)
estimates$logRR <- estimates$yi
estimates$SE <- sqrt(estimates$vi)
estimates$RR <- exp(estimates$yi)
estimates$lower95 <- exp(estimates$yi - qnorm(.975) * estimates$SE)
estimates$upper95 <- exp(estimates$yi + qnorm(.975) * estimates$SE)
blocks <- lapply(groups, function(index) pool_estimates(estimates[index, ]))

arms <- c("Chemo alone", "ICI plus chemo")
ticks <- c(.25, 1, 4, 16)
draw_forest(setNames(blocks["Overall"], ""), "a_forest.pdf", "Pathological complete response",
            ticks, arms = arms, metric = "RR")
draw_forest(blocks[c("FLOT-based", "Non-FLOT-based")], "b_forest.pdf",
            "Pathological complete response by chemotherapy backbone", ticks, arms = arms, metric = "RR")
