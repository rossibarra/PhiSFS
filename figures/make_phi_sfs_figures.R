#!/usr/bin/env Rscript

# Deterministic vector versions of the two Phi-SFS explanatory figures.
# Run from the repository root with:
#   module load R/4.4.2
#   Rscript figures/make_phi_sfs_figures.R

orange <- "#E66101"
orange_fill <- "#F3C6AA"
teal <- "#257D7D"
teal_fill <- "#B9D7D7"
neutral <- "#9E9E9E"
neutral_fill <- "#D0D0D0"
ink <- "#202020"

normalize <- function(x) x / sum(x)

step_coordinates <- function(x, y) {
    list(
        x = c(x[1], rep(x[-1], each = 2)),
        y = c(y[1], as.vector(rbind(y[-length(y)], y[-1])))
    )
}

draw_sfs <- function(x, neutral_sfs, te_sfs, color, title, show_y = FALSE) {
    ymax <- max(c(neutral_sfs, te_sfs)) * 1.08
    plot(NA, xlim = c(0.02, 0.98), ylim = c(0, ymax), axes = FALSE,
         xlab = "", ylab = "", xaxs = "i", yaxs = "i")

    # The wider neutral bars remain a single uniform gray. Narrower opaque TE
    # bars are drawn in front, avoiding transparency-induced color changes.
    rect(x - 0.022, 0, x + 0.022, neutral_sfs,
         col = neutral_fill, border = neutral, lwd = 1.5)
    rect(x - 0.015, 0, x + 0.015, te_sfs,
         col = adjustcolor(color, alpha.f = 0.58), border = color, lwd = 2.1)

    axis(1, at = c(0.05, 0.95), labels = c("0.05", "0.95"),
         lwd = 2.2, lwd.ticks = 2.2, cex.axis = 1.2)
    if (show_y) {
        axis(2, at = c(0, ymax), labels = c("0", ""),
             lwd = 2.2, lwd.ticks = 2.2, cex.axis = 1.1)
    }
    box(bty = "l", lwd = 2.2)
    mtext("DAF", side = 1, line = 2.2, cex = 1.3, font = 2)
    title(main = title, col.main = color, cex.main = 1.35, font.main = 2,
          line = 0.45)
    legend("topright", legend = c("SNPs", "TEs"),
           fill = c(neutral_fill, adjustcolor(color, alpha.f = 0.58)),
           border = c(neutral, color), bty = "n", cex = 1.15,
           text.font = 2, x.intersp = 0.6, y.intersp = 0.9)
}

draw_cdf_arrow <- function(color) {
    plot.new()
    plot.window(xlim = c(0, 1), ylim = c(0, 1))
    arrows(0.47, 0.88, 0.47, 0.16, length = 0.13, angle = 28,
           lwd = 4.0, col = color)
    text(0.62, 0.53, "CDF", cex = 1.25, font = 2, col = ink)
}

draw_cdf <- function(x, neutral_sfs, te_sfs, color, fill,
                     show_y = FALSE, annotation_x = 0.48,
                     annotation_y = 0.34, arrow_x = 0.34) {
    neutral_cdf <- cumsum(neutral_sfs)
    te_cdf <- cumsum(te_sfs)
    neutral_step <- step_coordinates(x, neutral_cdf)
    te_step <- step_coordinates(x, te_cdf)

    plot(NA, xlim = c(0.02, 0.98), ylim = c(0, 1.02), axes = FALSE,
         xlab = "", ylab = "", xaxs = "i", yaxs = "i")
    polygon(c(neutral_step$x, rev(te_step$x)),
            c(neutral_step$y, rev(te_step$y)),
            col = fill, border = NA)
    lines(neutral_step$x, neutral_step$y, col = neutral, lwd = 3.0)
    lines(te_step$x, te_step$y, col = color, lwd = 3.2)

    axis(1, at = c(0.05, 0.95), labels = c("0.05", "0.95"),
         lwd = 2.2, lwd.ticks = 2.2, cex.axis = 1.2)
    if (show_y) {
        axis(2, at = c(0, 1), labels = c("0", "1"),
             lwd = 2.2, lwd.ticks = 2.2, cex.axis = 1.2)
    }
    box(bty = "l", lwd = 2.2)
    mtext("DAF", side = 1, line = 2.2, cex = 1.3, font = 2)

    target_bin <- max(which(x <= arrow_x))
    target_y <- (te_cdf[target_bin] + neutral_cdf[target_bin]) / 2
    text(annotation_x, annotation_y, expression(Phi[SFS]),
         cex = 1.45, font = 2, col = ink)
    arrows(annotation_x - 0.03, annotation_y + 0.08,
           arrow_x, target_y, length = 0.10, lwd = 2.4, col = ink)
}

draw_vertical_label <- function(label, cex = 1.25) {
    plot.new()
    plot.window(xlim = c(0, 1), ylim = c(0, 1))
    text(0.55, 0.5, label, srt = 90, cex = cex, font = 2, xpd = NA)
}

make_definition_figure <- function(path) {
    bins <- 1:19
    x <- bins / 20
    neutral_sfs <- normalize(1 / bins)
    rare_sfs <- normalize(1 / bins^1.55)
    high_sfs <- normalize(0.55 / bins + 0.055 * exp(0.42 * (bins - 12)))

    pdf(path, width = 11.2, height = 8.6, family = "Helvetica",
        useDingbats = FALSE)
    layout(matrix(1:9, nrow = 3, byrow = TRUE),
           heights = c(1.0, 0.24, 1.0), widths = c(0.13, 1, 1))

    par(mar = rep(0, 4))
    draw_vertical_label("Proportion")

    par(mar = c(3.8, 3.2, 2.5, 1.0), mgp = c(2.1, 0.65, 0),
        tcl = -0.35, las = 1)

    draw_sfs(x, neutral_sfs, rare_sfs, orange,
             "Excess rare variants", show_y = TRUE)
    draw_sfs(x, neutral_sfs, high_sfs, teal,
             "Excess high-frequency derived variants", show_y = FALSE)

    par(mar = rep(0, 4))
    plot.new()
    draw_cdf_arrow(orange)
    draw_cdf_arrow(teal)

    par(mar = rep(0, 4))
    draw_vertical_label("Proportion")

    par(mar = c(3.8, 3.2, 0.8, 1.0), mgp = c(2.1, 0.65, 0),
        tcl = -0.35, las = 1)
    draw_cdf(x, neutral_sfs, rare_sfs, orange, orange_fill,
             show_y = TRUE, annotation_x = 0.49,
             annotation_y = 0.31, arrow_x = 0.30)
    draw_cdf(x, neutral_sfs, high_sfs, teal, teal_fill,
             show_y = FALSE, annotation_x = 0.69,
             annotation_y = 0.31, arrow_x = 0.57)
    dev.off()
}

standardize <- function(x) as.numeric(scale(x))

draw_violin <- function(values, center, width = 0.30) {
    d <- density(values, from = -2.5, to = 4.0, n = 512, cut = 0,
                 bw = "nrd0")
    half_width <- width * d$y / max(d$y)
    polygon(c(center - half_width, rev(center + half_width)),
            c(d$x, rev(d$x)), col = "#D0D0D0", border = "#B5B5B5",
            lwd = 1.2)
}

make_null_figure <- function(path) {
    n <- 50000
    u <- ((1:n) - 0.5) / n
    nulls <- list(
        standardize(qgamma(u, shape = 3.0)),
        standardize(qt(u, df = 6)),
        standardize(c(qnorm(u[seq_len(n / 2)], -0.75, 0.30),
                      qnorm(u[(n / 2 + 1):n], 0.55, 0.23))),
        standardize(qgamma(u, shape = 1.25))
    )
    observed_z <- c(3.4, 2.5, 1.6, 0.9)
    minus_log10_p <- c(3.0, 2.15, 1.15, 0.45)
    categories <- c("In gene", "0-2 kb", "2-5 kb", ">5 kb")
    blue <- colorRampPalette(c("#DEEBF7", "#6BAED6", "#08519C"))(301)
    point_colors <- blue[1 + round(minus_log10_p / 3 * 300)]

    pdf(path, width = 10.2, height = 7.2, family = "Helvetica",
        useDingbats = FALSE)
    layout(matrix(c(1, 2, 3), nrow = 1), widths = c(0.38, 4.7, 1.35))
    par(oma = c(0, 0, 2.4, 0))

    par(mar = rep(0, 4))
    draw_vertical_label(
        expression(paste("Null-standardized ", Phi[SFS], " (Z)")),
        cex = 1.1
    )

    par(mar = c(4.7, 4.3, 1.0, 0.6), mgp = c(2.8, 0.75, 0),
        tcl = -0.35, las = 1)
    plot(NA, xlim = c(0.45, 4.55), ylim = c(-2.5, 4.0), axes = FALSE,
         xlab = "", ylab = "", xaxs = "i", yaxs = "i")
    abline(h = 0, col = "#8F8F8F", lty = 2, lwd = 1.5)
    for (i in seq_along(nulls)) draw_violin(nulls[[i]], i)
    lines(1:4, observed_z, col = "#79A6D2", lwd = 2.1)
    points(1:4, observed_z, pch = 21, bg = point_colors,
           col = "white", lwd = 1.4, cex = 1.75)
    axis(1, at = 1:4, labels = categories, lwd = 2.0,
         lwd.ticks = 2.0, cex.axis = 1.05)
    axis(2, at = seq(-2, 4, by = 0.5), lwd = 2.0,
         lwd.ticks = 2.0, cex.axis = 1.0)
    box(bty = "l", lwd = 2.0)
    mtext("Distance to nearest gene", side = 1, line = 3.1,
          cex = 1.2)

    par(mar = c(4.7, 0.4, 1.0, 0.5))
    plot.new()
    plot.window(xlim = c(0, 1), ylim = c(0, 1))
    text(0.5, 0.91, expression(-log[10](italic(P)*"-value")),
         cex = 1.05)
    edges <- seq(0.10, 0.90, length.out = length(blue) + 1)
    rect(edges[-length(edges)], 0.81, edges[-1], 0.865,
         col = blue, border = NA)
    rect(0.10, 0.81, 0.90, 0.865, border = ink, lwd = 1.2)
    tick_x <- 0.10 + (0:3) / 3 * 0.80
    segments(tick_x, 0.79, tick_x, 0.81, lwd = 1.1)
    text(tick_x, 0.755, labels = 0:3, cex = 0.95)
    rect(0.12, 0.61, 0.24, 0.67, col = "#D0D0D0",
         border = "#B5B5B5", lwd = 1.2)
    text(0.31, 0.64, "Null distribution", adj = 0, cex = 1.0)

    mtext("Illustrative example", side = 3, line = 0.5,
          outer = TRUE, cex = 1.55, font = 2)
    dev.off()
}

make_definition_figure("figures/phi_sfs_definition_schematic.pdf")
make_null_figure("figures/phi_sfs_null_standardization_example.pdf")
