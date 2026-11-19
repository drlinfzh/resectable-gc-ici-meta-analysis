# Supplementary Fig. 2. Reconstructed trial-level survival curves.
library(survival)
library(grid)

input_dir <- "inputs"

output_dir <- "outputs/sfig2"
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

logrank_p <- function(data) {
  test <- survdiff(Surv(time, status) ~ exposure, data = data, rho = 0)
  pchisq(test$chisq, df = 1L, lower.tail = FALSE)
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

panels <- data.frame(
  stem = c("MATTERHORN_EFS", "MATTERHORN_OS",
           "KEYNOTE-585_Overall_EFS", "KEYNOTE-585_Overall_OS",
           "KEYNOTE-585_Main_EFS", "KEYNOTE-585_Main_OS", "ASTRUM-006_EFS",
           "KEYNOTE-585_FLOT_EFS", "KEYNOTE-585_FLOT_OS"),
  dataset = c("ipd/mth_efs_ipd.csv", "ipd/mth_os_ipd.csv",
              "ipd/kn585_efs_ipd.csv", "ipd/kn585_os_ipd.csv",
              "ipd/kn585_efs_main_ipd.csv", "ipd/kn585_os_main_ipd.csv", "ipd/ast006_efs_ipd.csv",
              "kmsubtraction/kn585_flot_efs_kmsubtraction_ipd.csv",
              "kmsubtraction/kn585_flot_os_kmsubtraction_ipd.csv"),
  endpoint = c("EFS", "OS", "EFS", "OS", "EFS", "OS", "EFS", "EFS", "OS"),
  title = c(rep("MATTERHORN", 2), rep("KEYNOTE-585 main plus FLOT cohort", 2),
            rep("KEYNOTE-585 main cohort", 2), "ASTRUM-006", rep("KEYNOTE-585 FLOT cohort", 2)),
  active = c(rep("Durva + FLOT", 2), rep("Pembro + XP/FP/FLOT", 2),
             rep("Pembro + XP/FP", 2), "Serplu + SOX", rep("Pembro + FLOT", 2)),
  comparator = c(rep("Placebo + FLOT", 2), rep("Placebo + XP/FP/FLOT", 2),
                 rep("Placebo + XP/FP", 2), "Placebo + SOX", rep("Placebo + FLOT", 2)),
  citation = c("Janjigian et al. 2025 N Engl J Med, Fig. 1",
               "Janjigian et al. 2026 Lancet, Fig. 2A",
               "Shitara et al. 2025 J Clin Oncol, Fig. 1B",
               "Shitara et al. 2025 J Clin Oncol, Fig. 1D",
               "Shitara et al. 2025 J Clin Oncol, Fig. 1A",
               "Shitara et al. 2025 J Clin Oncol, Fig. 1C",
               "Shen et al. 2026 Lancet, Fig. 2B",
               "Shitara et al. 2025 J Clin Oncol, Figs. 1A–B, A5",
               "Shitara et al. 2025 J Clin Oncol, Figs. 1C–D, 2B")
)
published <- read.csv(file.path(input_dir, "reference", "published_survival.csv"))
flot_published <- read.csv(file.path(input_dir, "reference", "published_flot_survival.csv"))
published <- rbind(published, flot_published)
for (i in seq_len(nrow(panels))) {
  panel <- panels[i, ]
  data <- read.csv(file.path(input_dir, panel$dataset))
  data$exposure <- data$treat
  estimate <- cox_estimate(data)
  # The FLOT curves use the two-sided Cox Wald P value.
  if (i <= 7L) estimate$P <- logrank_p(data)
  reported <- published[published$stem == panel$stem, ]
  stopifnot(nrow(reported) == 1L)
  draw_km(data, estimate, paste0(letters[i], "_km.pdf"), panel$endpoint,
          title = panel$title, arms = c(panel$active, panel$comparator),
          reported = reported, citation = panel$citation,
          xmax = if (panel$stem == "KEYNOTE-585_FLOT_OS") 78 else NULL)
}
