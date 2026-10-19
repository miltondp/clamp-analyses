# Export SVG figures while preserving bold text and spacing.
svg_device <- function(...) {
  svglite::svglite(..., system_fonts = list(Helvetica = "Liberation Sans"))
}

add_svg_font_stack <- function(path) {
  x <- readLines(path, warn = FALSE)
  x <- gsub("font-family: \"Liberation Sans\";",
            "font-family: \"Liberation Sans\", Helvetica, Arial, sans-serif;", x, fixed = TRUE)
  writeLines(x, path)
  invisible(path)
}
