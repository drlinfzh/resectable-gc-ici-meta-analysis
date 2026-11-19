# Supplementary Fig. 3. KEYNOTE-585 FLOT matching diagnostics.

library(survival)
suppressPackageStartupMessages(library(KMSubtraction))

input_dir <- "inputs"
output_dir <- "outputs/sfig3"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
cohorts <- c("Main plus FLOT", "Main cohort")
colours <- c("#2C2E35", "#0093A1")

for (endpoint in c("EFS", "OS")) {
  endpoint_title <- if (endpoint == "EFS") "Event-free survival" else "Overall survival"
  overall <- read.csv(file.path(input_dir, "ipd", paste0("kn585_", tolower(endpoint), "_ipd.csv")))
  main <- read.csv(file.path(input_dir, "ipd", paste0("kn585_", tolower(endpoint), "_main_ipd.csv")))
  for (arm in c(1L, 0L)) {
    treatment <- if (arm == 1L) "ICI plus chemo" else "Chemo alone"
    title <- paste(endpoint_title, if (arm == 1L) treatment else tolower(treatment), sep = ", ")
    prefix <- paste(tolower(endpoint), if (arm == 1L) "ici_chemo" else "chemo", sep = "_")
    overall_arm <- overall[overall$treat == arm, c("time", "status")]
    main_arm <- main[main$treat == arm, c("time", "status")]
    matched <- KMSubtractionMatch(overall_arm, main_arm, matching = "bipartite")$data
    matched <- matched[!is.na(matched$matchid), ]
    matched$curve <- factor(matched$curve, levels = c("overall", "subgroup"))
    stopifnot(all(table(matched$curve) == nrow(main_arm)), nrow(matched) == 2L * nrow(main_arm))
    xmax <- ceiling(max(matched$time) / 6) * 6
    diagnostic_ticks <- seq(0, xmax, by = 24)

    for (status in c(1L, 0L)) {
      kind <- if (status == 1L) "event" else "censoring"
      kind_title <- if (status == 1L) "Event times" else "Censoring times"
      a <- matched[matched$status == status & matched$curve == "overall", ]
      b <- matched[matched$status == status & matched$curve == "subgroup", ]
      a <- a[order(a$matchid), ]
      b <- b[order(b$matchid), ]
      stopifnot(nrow(a) == nrow(b), nrow(a) > 0)
      agreement <- blandr::blandr.statistics(a$time, b$time)
      ks_p <- suppressWarnings(ks.test(a$time, b$time)$p.value)

      file <- file.path(output_dir, paste(prefix, kind, "cumulative.pdf", sep = "_"))
      grDevices::cairo_pdf(file, family = "Arial")
      plot(ecdf(a$time), xlim = c(0, xmax), verticals = TRUE, do.points = FALSE,
           xlab = "Follow-up (months)", ylab = "Cumulative probability",
           main = paste(title, kind_title, sep = "\n"), col = colours[1], xaxt = "n")
      axis(1, at = diagnostic_ticks)
      lines(ecdf(b$time), verticals = TRUE, do.points = FALSE, col = colours[2])
      legend("bottomright", cohorts, col = colours, lty = 1, bty = "n")
      invisible(dev.off())

      file <- file.path(output_dir, paste(prefix, kind, "bland_altman.pdf", sep = "_"))
      grDevices::cairo_pdf(file, family = "Arial")
      par(mar = c(6.5, 4.5, 4.5, 1.5))
      limits <- range(c(agreement$differences, agreement$lowerLOA, agreement$upperLOA))
      plot(agreement$means, agreement$differences, pch = 16, xlim = c(0, xmax), ylim = limits,
           main = paste(title, kind_title, "Bland–Altman", sep = "\n"),
           xlab = "Mean follow-up (months)", ylab = "Overall − main (months)", xaxt = "n")
      axis(1, at = diagnostic_ticks)
      abline(h = agreement$bias, col = colours[1])
      abline(h = c(agreement$lowerLOA, agreement$upperLOA), lty = 2, col = colours[1])
      mtext(sprintf("ΣO = %.3f; ΣM = %.3f; MAD = %.5f; KS P = %.3f",
                    sum(a$time), sum(b$time), mean(abs(agreement$differences)), ks_p),
            side = 1, line = 4.5, cex = .8)
      invisible(dev.off())
    }

    fit <- survfit(Surv(time, status) ~ curve, data = matched)
    grDevices::cairo_pdf(file.path(output_dir, paste0(prefix, "_km.pdf")), family = "Arial")
    par(mar = c(8, 8, 3, 1))
    times <- seq(0, xmax, by = 6)
    plot(fit, col = colours, xlim = c(0, xmax), ylim = c(0, 1.035), conf.int = FALSE,
         mark.time = TRUE, main = title, xlab = "Time (months)",
         ylab = paste0(endpoint_title, " (%)"), xaxt = "n", yaxt = "n")
    axis(1, at = times)
    axis(2, at = seq(0, 1, .2), labels = seq(0, 100, 20), las = 1)
    legend("bottomleft", cohorts, col = colours, lty = 1, bty = "n")
    risk_label_x <- -.05 * xmax
    mtext("Number at risk", side = 1, line = 4, at = risk_label_x, adj = 1, font = 2)
    for (i in seq_along(cohorts)) {
      followup <- matched$time[matched$curve == levels(matched$curve)[i]]
      n_risk <- vapply(times, function(t) sum(followup >= t), integer(1))
      mtext(cohorts[i], side = 1, line = 4 + i, at = risk_label_x, adj = 1, cex = .8)
      mtext(n_risk, side = 1, line = 4 + i, at = times, col = colours[i], cex = .8)
    }
    invisible(dev.off())
  }
}
