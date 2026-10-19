# Format gene-module identifiers for publication panels.
gm_ids <- function(x) gsub("\\bLV([0-9]+)\\b", "GM\\1", x)
fread_gm <- function(...) {
  d <- data.table::fread(...)
  for (j in names(d)[vapply(d, is.character, logical(1))]) data.table::set(d, j = j, value = gm_ids(d[[j]]))
  d
}
