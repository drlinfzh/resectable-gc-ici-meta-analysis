# Supplementary Fig. 4. Monte Carlo evaluation of reconstruction error.

library(survival)
suppressPackageStartupMessages(library(KMSubtraction))

input_dir <- "inputs"
output_dir <- "outputs/sfig4"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

methods <- c("bipartite", "maha", "logit")
method_labels <- c("Bipartite", "Mahalanobis", "Logistic")
colours <- c("#00897B", "#E69F00", "#4477AA")
repetitions <- as.integer(Sys.getenv("MC_REPETITIONS", "5000"))
stopifnot(repetitions >= 2L)

# KMSubtractionError simulation, adapted from KMSubtraction commit 82044ecd.
digitize_curve <- function(data) {
  max_time <- max(data$time)
  fit <- survfit(Surv(time, status) ~ 1, data = data)
  svg <- svglite::svgstring(width = 5, height = 5)
  plot(fit, col = "#f042ed", conf.int = FALSE, lwd = 1, xlim = c(0, max_time))
  points(c(0, max_time, 0, max_time), c(0, 0, 1, 1), col = "#f04242", pch = 16, cex = .3)
  invisible(dev.off())
  image <- magick::image_read_svg(svg(), width = 1000)
  raster <- magick::image_raster(image, tidy = FALSE)
  magick::image_destroy(image)
  colour <- as.vector(raster)
  k <- which(colour %in% c("#f04242ff", "#f042edff"))
  x <- scales::rescale((k - 1) %% ncol(raster) + 1, to = c(0, max_time))
  y <- scales::rescale(-((k - 1) %/% ncol(raster) + 1), to = c(0, 1))
  data.frame(x, y)[colour[k] == "#f042edff", ]
}

reconstruct_curve <- function(data, intervals) {
  fit <- survfit(Surv(time, status) ~ 1, data = data)
  times <- seq(0, max(fit$time), by = max(data$time) / intervals)
  risk <- survminer:::.get_timepoints_survsummary(fit, data, times)
  points <- digitize_curve(data)
  prepared <- IPDfromKM::preprocess(dat = points, trisk = risk$time, nrisk = risk$n.risk, maxy = 1)
  IPDfromKM::getIPD(prep = prepared, armID = 1, tot.events = NULL)$IPD
}

simulate_once <- function(iteration, parameters) {
  set.seed(parameters$seed + iteration)
  n <- parameters$n
  n_main <- round(parameters$subgroup.p * n)
  n_flot <- n - n_main
  censored_main <- round(parameters$censor_subgroup.p * n_main)
  censored_flot <- round(parameters$censor_overall.p * n) - censored_main
  status <- c(rep(0, censored_main), rep(1, n_main - censored_main),
              rep(0, censored_flot), rep(1, n_flot - censored_flot))
  overall <- data.frame(time = rweibull(n, shape = 1, scale = 5), status = status,
                        subgroup = rep(1:0, c(n_main, n_flot)))
  overall$subgroup[overall$subgroup == 0] <-
    wakefield::r_sample_binary(n_flot, x = c(NA, 0), prob = c(0, 1), name = "Binary")
  main <- overall[overall$subgroup == 1, ]
  full_ipd <- reconstruct_curve(overall, parameters$interval)
  main_ipd <- reconstruct_curve(main, parameters$interval)
  overall$strata <- "original"
  vapply(methods, function(method) {
    matched <- KMSubtractionMatch(full_ipd, main_ipd, matching = method)$data
    matched$strata <- "matched"
    compare <- dplyr::bind_rows(matched, overall)
    fit <- coxph(Surv(time, status) ~ strata, data = compare[which(compare$subgroup == 0), ])
    unname(coef(fit)[1])
  }, numeric(1))
}

parameters <- do.call(rbind, lapply(c("EFS", "OS"), function(endpoint) {
  overall <- read.csv(file.path(input_dir, "ipd", paste0("kn585_", tolower(endpoint), "_ipd.csv")))
  main <- read.csv(file.path(input_dir, "ipd", paste0("kn585_", tolower(endpoint), "_main_ipd.csv")))
  do.call(rbind, lapply(0:1, function(arm) {
    a <- overall[overall$treat == arm, ]
    b <- main[main$treat == arm, ]
    data.frame(endpoint, arm, n = nrow(a), subgroup.p = nrow(b) / nrow(a),
               censor_overall.p = mean(a$status == 0), censor_subgroup.p = mean(b$status == 0),
               interval = if (endpoint == "EFS") 12 else 13, missing.p = 0,
               seed = 585000 + if (endpoint == "EFS") arm * 10000 else 20000 + arm * 10000)
  }))
}))

reuse <- "--reuse-simulations" %in% commandArgs(trailingOnly = TRUE)
if (reuse) {
  saved <- readRDS(file.path(input_dir, "diagnostics", "monte_carlo_results.rds"))
  stopifnot(identical(saved$parameters, parameters), identical(saved$mc, repetitions))
}

for (j in seq_len(nrow(parameters))) {
  p <- parameters[j, ]
  if (reuse) {
    rows <- do.call(rbind, lapply(saved$results[[j]], "[[", "data"))
    errors <- lapply(methods, function(method) rows$TE[rows$matching == method])
  } else {
    draws <- vapply(seq_len(repetitions), simulate_once, numeric(length(methods)), parameters = p)
    errors <- lapply(seq_along(methods), function(i) draws[i, ])
  }
  stopifnot(all(lengths(errors) == repetitions), all(is.finite(unlist(errors))))
  endpoint <- if (p$endpoint == "EFS") "Event-free survival" else "Overall survival"
  arm <- if (p$arm == 0) "Chemo alone" else "ICI plus chemo"
  title <- paste(endpoint, arm, sep = "\n")
  prefix <- paste(tolower(p$endpoint), if (p$arm == 0) "chemo" else "ici_chemo", sep = "_")

  for (absolute in c(FALSE, TRUE)) {
    values <- if (absolute) lapply(errors, abs) else errors
    bandwidth <- bw.nrd0(unlist(values))
    limits <- if (absolute) c(0, max(unlist(values)) + 3 * bandwidth) else
      range(unlist(values)) + c(-3, 3) * bandwidth
    densities <- lapply(values, function(x) {
      d <- density(if (absolute) c(-x, x) else x, bw = bandwidth,
                   from = limits[1], to = limits[2], n = 1024)
      if (absolute) d$y <- 2 * d$y
      d
    })
    kind <- if (absolute) "absolute_error" else "signed_error"
    panel <- if (absolute) "Absolute error" else "Signed error"
    grDevices::cairo_pdf(file.path(output_dir, paste0(prefix, "_", kind, ".pdf")), family = "Arial")
    ymax <- max(vapply(densities, function(d) max(d$y), numeric(1)))
    plot(densities[[1]], type = "n", xlim = limits, ylim = c(0, 1.1 * ymax),
         main = paste(title, panel, sep = "\n"), xlab = if (absolute) "|ln(HR)|" else "ln(HR)",
         ylab = "Density")
    for (i in seq_along(methods)) {
      d <- densities[[i]]
      polygon(c(d$x[1], d$x, tail(d$x, 1)), c(0, d$y, 0),
              col = adjustcolor(colours[i], alpha.f = .15), border = NA)
    }
    for (i in seq_along(methods)) {
      lines(densities[[i]], col = colours[i], lty = i, lwd = 1.5)
      abline(v = mean(values[[i]]), col = "#2C2E35", lty = 2)
    }
    abline(v = 0, col = "#D64D4D", lty = 2)
    legend("topright", method_labels, col = colours, lty = seq_along(methods), bty = "n")
    invisible(dev.off())
  }

  convergence <- lapply(errors, function(x) {
    error <- abs(x)
    k <- seq_along(error)
    mean <- cumsum(error) / k
    sd <- sqrt(pmax(0, (cumsum(error^2) - k * mean^2) / pmax(1, k - 1L)))
    se <- sd / sqrt(k)
    data.frame(k, mean, lower = mean - 1.96 * se, upper = mean + 1.96 * se)
  })
  grDevices::cairo_pdf(file.path(output_dir, paste0(prefix, "_convergence.pdf")), family = "Arial")
  limits <- range(unlist(lapply(convergence, function(d) c(d$lower, d$upper))))
  plot(convergence[[1]]$k, convergence[[1]]$mean, type = "n", ylim = limits,
       main = paste(title, "Convergence", sep = "\n"),
       xlab = "Completed iteration", ylab = "Running mean |ln(HR)|")
  for (i in seq_along(methods)) {
    d <- convergence[[i]]
    polygon(c(d$k, rev(d$k)), c(d$lower, rev(d$upper)), border = NA,
            col = adjustcolor(colours[i], alpha.f = .12))
  }
  for (i in seq_along(methods))
    lines(convergence[[i]]$k, convergence[[i]]$mean, col = colours[i], lty = i, lwd = 1.5)
  legend("topright", method_labels, col = colours, lty = seq_along(methods), bty = "n")
  invisible(dev.off())
}
