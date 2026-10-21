# Project directories, relative to this file's location.
local({
  arquivos <- Filter(Negate(is.null), lapply(sys.frames(), function(f) f$ofile))
  if (!length(arquivos)) stop("Load caminhos.R with source().")
  este <- normalizePath(tail(arquivos, 1)[[1]], winslash = "/", mustWork = TRUE)
  assign("RAIZ", dirname(este), envir = globalenv())
})
DADOS <- file.path(RAIZ, "dados")
ESTIMACOES <- file.path(RAIZ, "estimacoes")
ESPACIAL <- Sys.getenv("BMA_ESPACIAL", file.path(RAIZ, "espacial"))
MANUSCRITO <- file.path(RAIZ, "manuscrito")
FIGURAS <- file.path(RAIZ, "figuras")
DOCS <- file.path(RAIZ, "docs")
local({
  faltando <- c(ESTIMACOES, ESPACIAL)[!dir.exists(c(ESTIMACOES, ESPACIAL))]
  if (length(faltando)) stop("Missing directories: ", paste(faltando, collapse = ", "))
})
