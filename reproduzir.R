# See docs/EXECUTION.md for the order and cost of each stage.
args_full <- commandArgs(FALSE)
if (!isTRUE(l10n_info()[["UTF-8"]])) {
  for (loc in c(".UTF-8", "C.UTF-8", "en_US.UTF-8", "English_United States.utf8")) {
    suppressWarnings(try(Sys.setlocale("LC_CTYPE", loc), silent = TRUE))
    if (isTRUE(l10n_info()[["UTF-8"]])) break
  }
  if (!isTRUE(l10n_info()[["UTF-8"]])) stop("Set a UTF-8 locale before running the scripts.")
}
script <- sub("^--file=", "", args_full[grepl("^--file=", args_full)])
if (length(script) != 1L) stop("Run with Rscript reproduzir.R <stage>.")
root <- dirname(normalizePath(script, winslash = "/", mustWork = TRUE))
args <- commandArgs(TRUE)
etapa <- if (length(args)) args[1] else "verificar"
if (!etapa %in% c("verificar", "codigo", "ambiente", "preparar", "bma", "geometria", "auxiliares", "interpolacao", "validar", "teste", "estimar", "exportar", "dados-figuras"))
  stop("Unknown stage: ", etapa)
setwd(root)
run <- function(path) {
  before <- getwd(); on.exit(setwd(before))
  path <- normalizePath(file.path(root, path), winslash = "/", mustWork = TRUE)
  setwd(dirname(path))
  source(basename(path), local = new.env(parent = globalenv()), encoding = "UTF-8")
}
need <- function(packages) {
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("Install packages: ", paste(missing, collapse = ", "))
}
if (etapa == "verificar") run("scripts/verificar.R")
if (etapa == "codigo") run("scripts/sintaxe.R")
if (etapa == "ambiente") run("scripts/ambiente.R")
if (etapa == "bma") {
  need(c("readxl", "dplyr", "BMS"))
  dest <- file.path(root, "execucoes_locais", paste0("bma_", format(Sys.time(), "%Y%m%d_%H%M%S")))
  if (dir.exists(dest)) stop("Directory already exists: ", dest)
  Sys.setenv(COVID_BMA_SOMENTE_PREPARAR = "0", COVID_BMA_MODO_TESTE = "0", COVID_BMA_SAIDA = dest)
  run("estimacoes/covid_bma.R")
}
if (etapa == "preparar") {
  need(c("readxl", "dplyr", "BMS"))
  dest <- file.path(root, "execucoes_locais", "preparacao")
  Sys.setenv(COVID_BMA_SOMENTE_PREPARAR = "1", COVID_BMA_SAIDA = dest,
             COVID_BMA_MODO_TESTE = "0")
  run("estimacoes/covid_bma.R")
  grupos <- list.dirs(file.path(root, "estimacoes/resultados/estimacoes/matrizes_estimacao"), recursive = FALSE, full.names = FALSE)
  for (grupo in grupos) for (arquivo in list.files(file.path(root, "estimacoes/resultados/estimacoes/matrizes_estimacao", grupo), pattern = "[.]csv$")) {
    ref <- read.csv(file.path(root, "estimacoes/resultados/estimacoes/matrizes_estimacao", grupo, arquivo))
    novo <- read.csv(file.path(dest, "matrizes_estimacao", grupo, arquivo))
    cmp <- all.equal(ref, novo, tolerance = 1e-12)
    if (!isTRUE(cmp)) stop("Matrix mismatch: ", grupo, "/", arquivo, ": ", paste(cmp, collapse = "; "))
    cat("Rebuilt matrix matches: ", grupo, "/", arquivo, "\n", sep = "")
  }
}
if (etapa == "geometria") run("scripts/geometria.R")
if (etapa == "auxiliares") {
  need(c("BMS", "spdep", "spatialreg"))
  dest <- file.path(root, "execucoes_locais", "auxiliares")
  if (dir.exists(dest)) stop("Auxiliary directory already exists: ", dest)
  dir.create(file.path(dest, "saida"), recursive = TRUE)
  Sys.setenv(BMA_ESPACIAL = dest)
  for (s in c("03_logmarginal_sar.R", "04_otimizar_logmarginal.R", "11_hyperg_sar.R"))
    run(file.path("espacial", s))
  for (f in c("base_sar.rds", "grade_sar.rds", "hyperg_aux.rds")) {
    cmp <- all.equal(readRDS(file.path(root, "espacial/saida", f)),
                     readRDS(file.path(dest, "saida", f)), tolerance = 1e-12)
    if (!isTRUE(cmp)) stop("Auxiliary object mismatch: ", f, ": ", paste(cmp, collapse = "; "))
    cat("Rebuilt auxiliary object matches: ", f, "\n", sep = "")
  }
}
if (etapa == "interpolacao") {
  target <- file.path(root, "espacial/saida/tabela_f21_estendida.rds")
  if (file.exists(target)) stop("Table already exists. Use a clean copy to regenerate it.")
  Sys.setenv(N_TAB = "300001", Z_MAX_TAB = "0.9", ARQ_TAB = "tabela_f21_estendida.rds")
  run("espacial/15_otimizar_hyperg.R")
}
if (etapa == "validar") {
  need("spdep")
  run("espacial/72_validar_escala_tres.R")
  run("espacial/73_validar_troca_tres_vias.R")
}
if (etapa %in% c("teste", "estimar")) {
  need("spdep")
  if (etapa == "teste") {
    Sys.setenv(MC3_SORTEIOS = "30", MC3_BURN = "10")
    cat("SMOKE TEST: do not use these chains for inference.\n")
  } else {
    Sys.setenv(MC3_SORTEIOS = "1000000", MC3_BURN = "200000")
  }
  output <- file.path(root, "execucoes_locais", paste0(etapa, "_", format(Sys.time(), "%Y%m%d_%H%M%S")))
  if (dir.exists(output)) stop("Directory already exists: ", output)
  dir.create(output, recursive = TRUE)
  Sys.setenv(SAIDA_TRES = output)
  for (d in c("casos", "obitos", "letalidade")) for (w in c("queen", "rook")) {
    Sys.setenv(DESFECHO = d, VIZINHOS = w)
    run("espacial/74_probabilidade_tres_vias.R")
  }
  writeLines(capture.output(sessionInfo()), file.path(output, "sessionInfo.txt"))
}
if (etapa == "exportar") run("scripts/exportar.R")
if (etapa == "dados-figuras") {
  need("sf")
  Sys.setenv(SAIDA_DADOS_FIGURAS = file.path(root, "execucoes_locais", "dados_figuras"))
  run("espacial/50_exportar_dados_figuras.R")
}
