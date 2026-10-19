# Shared layout and export helpers for main figures.
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(cowplot)
  library(yaml)
  library(here)
})

FS_TAG          <- 8
FS_TITLE        <- 6.5
FS_AXIS_TITLE   <- 5.5
FS_AXIS_TEXT    <- 5
FS_LEGEND       <- 5
FS_LEGEND_TITLE <- 5.5
FIG_W_MAX <- 180
FIG_H_MAX <- 170

source(here("scripts/panels/svg_export.R"))
cfg <- yaml::read_yaml(here("config.yaml"))
MODEL_COLORS <- unlist(cfg$MODEL_COLORS, use.names = TRUE)
ARCHS4_MODEL_COLORS   <- MODEL_COLORS[c("CLAMPfull", "CLAMPbase")]
COMPENDIUM_COLORS <- unlist(cfg$COMPENDIUM_COLORS, use.names = TRUE)
ARCHS4_DATASET_COLORS <- COMPENDIUM_COLORS[c("ARCHS4", "GTEx", "recount2")]
FRACTION_COLORS <- unlist(cfg$ARCHS4_FRACTION_COLORS, use.names = TRUE)
DATA_DIR <- here("output/03_model_biology/02_archs4")

theme_nm <- function(base_size = FS_AXIS_TEXT, base_line_size = 0.25) {
  theme_classic(base_size = base_size, base_family = "Helvetica") %+replace%
    theme(
      axis.line = element_line(linewidth = base_line_size, colour = "black"),
      axis.ticks = element_line(linewidth = base_line_size, colour = "black"),
      axis.ticks.length = unit(1.1, "mm"),
      axis.text = element_text(size = FS_AXIS_TEXT, colour = "black"),
      axis.title = element_text(size = FS_AXIS_TITLE, face = "plain", colour = "black"),
      legend.text = element_text(size = FS_LEGEND),
      legend.title = element_text(size = FS_LEGEND_TITLE, face = "plain"),
      legend.key = element_blank(),
      legend.background = element_blank(),
      strip.background = element_blank(),
      strip.text = element_text(size = FS_TITLE, face = "bold", margin = margin(b = 0.7)),
      panel.grid = element_blank(),
      plot.margin = margin(1.3, 1.3, 1.3, 1.3)
    )
}


add_tag <- function(p, label, x = 0.002, y = 0.998, top = 0.94) {
  ggdraw() +
    draw_plot(p, x = 0.015, y = 0, width = 0.985, height = top) +
    draw_label(tolower(label), x = x, y = y, hjust = 0, vjust = 1,
               size = FS_TAG, fontface = "bold", fontfamily = "Helvetica")
}

nm_p <- function(p) ifelse(p < .001, format(p, scientific = TRUE, digits = 2), sprintf("%.3f", p))

main_figure_dirs <- function(name) {
  out <- here("output/99_panels", name)
  src <- file.path(out, "source_data")
  dir.create(src, recursive = TRUE, showWarnings = FALSE)
  list(out = out, source = src)
}


clamp_fonts <- function(g) {
  if (inherits(g, "text")) {
    fs <- g$gp$fontsize
    if (is.null(fs)) fs <- 5
    g$gp$fontsize <- ifelse(fs == FS_TAG, FS_TAG, pmax(5, pmin(7, fs)))
    g$gp$fontfamily <- "Helvetica"
  }
  if (!is.null(g$grobs)) g$grobs <- lapply(g$grobs, clamp_fonts)
  if (!is.null(g$children)) for (i in seq_along(g$children)) g$children[[i]] <- clamp_fonts(g$children[[i]])
  g
}

save_main_figure <- function(fig, name, width, height) {
  stopifnot(width <= FIG_W_MAX, height <= FIG_H_MAX)
  out_dir <- main_figure_dirs(name)$out
  g <- clamp_fonts(grid::forceGrob(ggplotGrob(fig)))
  raw_pdf <- tempfile(fileext = ".pdf")
  ggsave(raw_pdf, g, width = width, height = height, units = "mm",
         device = cairo_pdf, bg = "white", limitsize = FALSE)
  dest <- file.path(out_dir, paste0(name, ".pdf"))
  if (nzchar(Sys.which("gs"))) {
    status <- system2("gs", c(
      "-q", "-dNOPAUSE", "-dBATCH", "-sDEVICE=pdfwrite",
      sprintf("-dDEVICEWIDTHPOINTS=%.6f", width / 25.4 * 72),
      sprintf("-dDEVICEHEIGHTPOINTS=%.6f", height / 25.4 * 72),
      "-dFIXEDMEDIA", "-dPDFFitPage", paste0("-sOutputFile=", dest), raw_pdf))
    stopifnot(status == 0L)
  } else {
    file.copy(raw_pdf, dest, overwrite = TRUE)
  }
  unlink(raw_pdf)
  ggsave(file.path(out_dir, paste0(name, ".svg")), g, width = width, height = height,
         units = "mm", device = svg_device, bg = "white", limitsize = FALSE)
  add_svg_font_stack(file.path(out_dir, paste0(name, ".svg")))
  ggsave(file.path(out_dir, paste0(name, ".png")), g, width = width, height = height,
         units = "mm", dpi = 600, device = ragg::agg_png, bg = "white", limitsize = FALSE)
  cat(sprintf("%s: %.0f x %.0f mm -> %s\n", name, width, height, out_dir))
  invisible(g)
}
