# Install the R packages used by the documented workflow.
pacotes <- c("readxl", "dplyr", "BMS", "spdep", "sf", "spatialreg")
faltantes <- pacotes[!vapply(pacotes, requireNamespace, logical(1), quietly = TRUE)]
repos <- getOption("repos")
if (!length(repos) || any(repos == "@CRAN@")) repos <- "https://cloud.r-project.org"
if (length(faltantes)) install.packages(faltantes, repos = repos)
for (p in pacotes) {
  tryCatch(loadNamespace(p), error = function(e)
    stop("Failed to load ", p, ": ", conditionMessage(e), call. = FALSE))
}
cat("Packages available. Record versions with: Rscript reproduzir.R ambiente\n")
