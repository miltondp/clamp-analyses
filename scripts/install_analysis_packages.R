# Install non-Conda comparison methods at the revisions used by this repository.
# Existing Conda environments are never removed or recreated by this script.
pins <- c(
  "wgmao/PLIER" = "fe4e9b23c47ee199afc5c984692df79ac5aabe80",
  "stephenslab/ebnm" = "825edf1e308af9f3e594a5482300b9fb99793129",
  "willwerscheid/flashier" = "45363f64c9ff01f0b63c49d5f3854fedef6ac1aa"
)
check_only <- "--check" %in% commandArgs(trailingOnly = TRUE)
options(repos = c(CRAN = "https://cloud.r-project.org"))
# CRAN dependencies absent from the Conda lock, pinned to the reference installation.
# Dependency order matters on a clean installation: fastTopics needs the pinned
# RcppParallel and ashr chain before its own installation.
cran_pins <- c(
  "RcppParallel" = "5.1.11-1",
  "mixsqp" = "0.3-54",
  "truncnorm" = "1.0-9",
  "trust" = "0.1-8",
  "deconvolveR" = "1.2-1",
  "quadprog" = "1.5-8",
  "etrunct" = "0.1",
  "invgamma" = "1.2",
  "ashr" = "2.2-63",
  "softImpute" = "1.4-3",
  "fastTopics" = "0.7-30",
  "Polychrome" = "1.5.4"
)
for (package in names(cran_pins)) {
  installed <- requireNamespace(package, quietly = TRUE)
  correct <- installed && utils::packageVersion(package) == package_version(cran_pins[[package]])
  if (!correct) {
    if (check_only) stop(package, " is missing or has a different version; run setup.sh")
    remotes::install_version(package, version = cran_pins[[package]], upgrade = "never", dependencies = FALSE)
    if (!requireNamespace(package, quietly = TRUE) || utils::packageVersion(package) != package_version(cran_pins[[package]])) {
      stop("Installation did not provide pinned version of ", package, ": ", cran_pins[[package]])
    }
  }
}
for (repo in names(pins)) {
  package <- basename(repo)
  description <- if (requireNamespace(package, quietly = TRUE)) utils::packageDescription(package) else NULL
  pinned <- !is.null(description) && identical(description$RemoteSha, unname(pins[[repo]]))
  if (!pinned) {
    if (check_only) stop(package, " is missing or has an unverified revision; run setup.sh")
    remotes::install_github(repo, ref = pins[[repo]], upgrade = "never", dependencies = FALSE, force = TRUE)
  }
  description <- utils::packageDescription(package)
  if (!identical(description$RemoteSha, unname(pins[[repo]]))) stop("Revision mismatch: ", package)
  message(package, " revision verified: ", pins[[repo]])
}

# Also catch dependency installers replacing an earlier pinned package.
for (package in names(cran_pins)) {
  if (!requireNamespace(package, quietly = TRUE) || utils::packageVersion(package) != package_version(cran_pins[[package]])) {
    stop("Final revision mismatch: ", package)
  }
}
