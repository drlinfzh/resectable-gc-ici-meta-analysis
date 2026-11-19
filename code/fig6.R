# Fig. 6. Adverse event percentages from the original aggregate counts.
output_dir <- "outputs/fig6"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

cohorts <- data.frame(
  cohort = c("MTH", "KN_main", "KN_FLOT", "NS01", "AST006"),
  n_ici = c(475, 399, 99, 54, 292),
  n_ctl = c(469, 400, 103, 54, 294)
)
# MATTERHORN: Lancet 2026, Table 2, and Ann Oncol 2026, abstract 515P.
# KEYNOTE-585: Lancet Oncol 2024, Table 2 and Table S6, and the final clinical study report.
# KEYNOTE-585 main-cohort all-grade anemia: J Clin Oncol 2025, Table A2.
# NEOSUMMIT-01: Nat Med 2024, Table 4, and the trial adverse-event records.
# ASTRUM-006: Lancet 2026, Table 3.

raw <- read.table(text = "
panel event cohort any_ici any_ctl high_ici high_ctl
A Any MTH 453 444 283 277
A Any KN_main 379 385 259 252
A Any KN_FLOT 98 101 75 65
A Any NS01 52 52 19 16
A Any AST006 284 288 136 172
A Nausea KN_main 218 234 24 22
A Nausea KN_FLOT 52 54 6 6
A Nausea NS01 4 9 0 1
A Nausea AST006 144 173 NA NA
A Neutropenia KN_main 112 109 70 68
A Neutropenia KN_FLOT 27 23 17 14
A Neutropenia NS01 31 25 9 7
A Diarrhea KN_main 113 104 15 16
A Diarrhea KN_FLOT 64 66 10 8
A Diarrhea NS01 8 10 0 1
A Anemia KN_main 101 108 18 30
A Anemia KN_FLOT 14 16 2 3
A Anemia AST006 147 152 18 26
A AST_increased KN_main NA NA 2 2
A AST_increased KN_FLOT 13 7 4 0
A AST_increased AST006 135 109 NA NA
B Any MTH 110 34 34 17
B Any KN_main 115 40 43 11
B Any KN_FLOT 45 29 21 9
B Any NS01 7 1 1 0
B Any AST006 117 19 18 3
B Hypothyroidism MTH 39 7 1 0
B Pneumonitis MTH 10 2 1 2
B Adrenal_insufficiency MTH 7 1 2 0
B Hyperthyroidism MTH 7 0 0 0
B Hypothyroidism KN_main 35 10 1 0
B Hypothyroidism KN_FLOT 6 5 0 0
B Hypothyroidism NS01 2 0 0 0
B Hyperthyroidism KN_main 15 5 0 0
B Hyperthyroidism KN_FLOT 6 4 0 0
B Adrenal_insufficiency KN_main 12 3 7 1
B Adrenal_insufficiency KN_FLOT 1 0 0 0
B Pneumonitis KN_main 16 3 2 0
B Pneumonitis KN_FLOT 3 3 1 0
B Hepatitis KN_main 7 1 7 1
B Hepatitis KN_FLOT 6 1 4 1
B Hepatitis NS01 1 0 0 0
B Hyperthyroidism NS01 NA NA NA NA
B Adrenal_insufficiency NS01 NA NA NA NA
B Pneumonitis NS01 NA NA NA NA
", header = TRUE)

d <- merge(raw, cohorts, by = "cohort", sort = FALSE)
stopifnot(!anyDuplicated(d[c("panel", "event", "cohort")]))
for (arm in c("ici", "ctl")) {
  a <- d[[paste0("any_", arm)]]
  h <- d[[paste0("high_", arm)]]
  n <- d[[paste0("n_", arm)]]
  stopifnot(all(a >= 0 & a <= n, na.rm = TRUE),
            all(h >= 0 & h <= n, na.rm = TRUE), all(h <= a, na.rm = TRUE))
}

# Sum events and denominators across reporting cohorts for each grade category.
all_grades <- aggregate(cbind(any_ici, any_ctl, n_ici, n_ctl) ~
                          panel + event, data = d, FUN = sum)
high_grades <- aggregate(cbind(high_ici, high_ctl, n_ici, n_ctl) ~
                           panel + event, data = d, FUN = sum)
pooled <- merge(all_grades, high_grades, by = c("panel", "event"),
                all.x = TRUE, suffixes = c("_any", "_high"))
events_by_panel <- list(
  A = c("Any", "Nausea", "Anemia", "AST_increased", "Diarrhea", "Neutropenia"),
  B = c("Any", "Hypothyroidism", "Hyperthyroidism", "Pneumonitis", "Hepatitis", "Adrenal_insufficiency")
)

base_colors <- c("#0871B9", "#FCAF17")
colors <- c(base_colors, adjustcolor(base_colors, alpha.f = 0.35))
panel_titles <- c(A = "Treatment-related adverse events",
                  B = "Immune-related adverse events")

for (panel in c("A", "B")) {
  grDevices::cairo_pdf(file.path(output_dir, paste0("adverse_events_", tolower(panel), ".pdf")),
                      width = 257 / 25.4, height = 98 / 25.4, family = "Arial", pointsize = 11)
  x <- pooled[pooled$panel == panel, ]
  x <- x[match(events_by_panel[[panel]], x$event), ]
  events <- rbind(x$any_ici, x$any_ctl, x$high_ici, x$high_ctl)
  n <- rbind(x$n_ici_any, x$n_ctl_any, x$n_ici_high, x$n_ctl_high)
  percentage <- 100 * events / n
  stopifnot(all(n > 0), all(percentage >= 0 & percentage <= 100))
  labels <- gsub("_", " ", x$event, fixed = TRUE)
  ymax <- if (panel == "A") 100 else 35
  par(mar = c(3.7, 4.2, 2.6, 0.5), mgp = c(2.7, 0.65, 0),
      family = "Arial", las = 1, cex = 1, xaxs = "i", yaxs = "i",
      col = "#2C2E35", fg = "#2C2E35", col.axis = "#2C2E35",
      col.lab = "#2C2E35", font.lab = 2)
  pos <- barplot(percentage, beside = TRUE, space = c(0.35, 0.7),
    xlim = c(0, ncol(percentage) * (4 + 0.7 + 3 * 0.35) + 0.5),
    col = colors, border = "#959699", lwd = 0.45,
    axisnames = FALSE, axes = FALSE,
    ylim = c(0, ymax * 1.11), ylab = "Percentage (%)", cex.lab = 1)
  axis(2, at = seq(0, ymax, if (panel == "A") 20 else 5), cex.axis = 0.86, tck = -0.015)
  axis(1, at = c(0, par("usr")[2]), labels = FALSE, tck = 0)
  text(colMeans(pos), -ymax * 0.065, labels, adj = c(0.5, 1), cex = 0.86, xpd = NA)
  label_y <- percentage + ymax * 0.012
  text(t(pos), t(label_y), sprintf("%.1f%%", t(percentage)),
       adj = c(0.5, 0), cex = 0.65)
  mtext(tolower(panel), side = 3, line = 1.1, adj = 0, cex = 1.2, font = 2)
  mtext(panel_titles[panel], side = 3, line = 1.1, adj = 0.065, cex = 1.05, font = 2)
  width <- diff(par("usr")[1:2])
  lx <- rep(c(0.78, 0.93) * width, each = 2)
  ly <- rep(c(0.80, 0.68) * ymax, 2)
  rect(lx - width * 0.026, ly - ymax * 0.035,
       lx + width * 0.026, ly + ymax * 0.035, col = colors,
       border = "#959699", lwd = 0.45)
  text(c(0.78, 0.93) * width, ymax * 0.94, c("All grades", "Grade ≥ 3"), cex = 0.86)
  text(width * 0.72, c(0.80, 0.68) * ymax, c("ICI plus chemo", "Chemo alone"), adj = 1, cex = 0.86)
  invisible(dev.off())
}
