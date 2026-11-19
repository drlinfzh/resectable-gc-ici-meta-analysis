# Fig. 3. Pathological complete response rates.
library(grid)

output_dir <- "outputs/fig3"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

ink <- "#2C2E35"
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

pdf_start <- function(filename, width, height) {
  grDevices::cairo_pdf(file.path(output_dir, filename), width = width / 25.4,
                       height = height / 25.4, family = font_family)
  grid.newpage()
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

count_columns <- c("e1", "n1", "e0", "n0")
add_total <- function(data) {
  total <- colSums(data[count_columns])
  rbind(data, data.frame(study = "Total", as.list(total)))
}

plot_rates <- function(data, labels, panel, backbone = FALSE) {
  width <- if (backbone) 256 else 218
  left <- 21
  right <- width - 5
  bottom <- 20
  y <- function(percent) bottom + 1.55 * percent
  centers <- left + (right - left) * (seq_len(nrow(data)) - .5) / nrow(data)
  colors <- c("#116FB0", "#F1AB22")

  pdf_start(paste0(panel, "_pcr_rates.pdf"), width, 105)
  label(panel, 4, 100, bold = TRUE, size = 14)
  for (tick in seq(0, 40, 10)) {
    line(left, y(tick), right, y(tick), col = row_grey, lwd = .5)
    line(left - 1.2, y(tick), left, y(tick), col = "#000000")
    label(tick, left - 2.4, y(tick), just = "right")
  }
  line(left, bottom, right, bottom, col = "#000000", lwd = 1)
  line(left, bottom, left, y(42), col = "#000000", lwd = 1)
  label("Pathological complete response (%)", 6, y(21), rot = 90, bold = TRUE, size = 11)
  label(labels, centers, bottom - 8, size = 9)
  if (backbone) {
    label(c("FLOT backbone", "Non-FLOT backbone"),
          c(mean(centers[1:3]), mean(centers[4:7])), 95, bold = TRUE)
    divider <- mean(centers[3:4])
    line(divider, bottom, divider, y(42), col = row_grey, lty = "24")
  } else {
    legend_x <- right - c(68, 29)
    grid.rect(unit(legend_x, "mm"), unit(c(94, 94), "mm"), width = unit(4, "mm"),
              height = unit(2.5, "mm"), gp = gpar(fill = colors, col = NA))
    label(c("ICI plus chemo", "Chemo alone"), legend_x + 4, 94, just = "left")
  }
  for (i in seq_len(nrow(data))) {
    events <- c(data$e1[i], data$e0[i])
    patients <- c(data$n1[i], data$n0[i])
    total <- data$study[i] == "Total"
    for (arm in 1:2) {
      rate <- 100 * events[arm] / patients[arm]
      x <- centers[i] + c(-6.3, 6.3)[arm]
      grid.rect(unit(x, "mm"), unit(bottom, "mm"), width = unit(7, "mm"),
                height = unit(y(rate) - bottom, "mm"), just = c("centre", "bottom"),
                gp = gpar(fill = colors[arm], col = NA))
      top <- rate
      if (!total) {
        ci <- 100 * binom.test(events[arm], patients[arm])$conf.int
        line(x, y(ci[1]), x, y(ci[2]), col = ci_grey)
        line(x - 1.2, y(ci), x + 1.2, y(ci), col = ci_grey)
        top <- ci[2]
      }
      label(sprintf("%.1f%%", rate), x, y(top) + 7.5, bold = TRUE, size = 9)
      label(sprintf("%d/%d", events[arm], patients[arm]), x, y(top) + 3.6, size = 8.5)
    }
  }
  invisible(dev.off())
}

overall <- add_total(trials[groups$Overall, ])
flot <- add_total(trials[groups[["FLOT-based"]], ])
nonflot <- add_total(trials[groups[["Non-FLOT-based"]], ])
stopifnot(all(unlist(tail(overall[count_columns], 1)) ==
                unlist(tail(flot[count_columns], 1) + tail(nonflot[count_columns], 1))))
plot_rates(overall, c("MATTERHORN", "KEYNOTE-585\nmain plus FLOT cohort",
                      "NEOSUMMIT-01", "ASTRUM-006", "Total"), "a")
plot_rates(rbind(flot, nonflot), c("MATTERHORN", "KEYNOTE-585\nFLOT cohort", "Total",
            "KEYNOTE-585\nmain cohort", "NEOSUMMIT-01", "ASTRUM-006", "Total"), "b", backbone = TRUE)
