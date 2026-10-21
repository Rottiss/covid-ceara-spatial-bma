source("../caminhos.R")
for (p in c("sf", "spdep")) if (!requireNamespace(p, quietly = TRUE)) stop("Install ", p)
malha <- readRDS(file.path(ESPACIAL, "saida/malha_ce.rds"))
dados <- read.csv(file.path(ESTIMACOES, "resultados/estimacoes/matrizes_estimacao/principal/casos.csv"),
                  colClasses = c(codigo_ibge = "character"))
ix <- match(dados$codigo_ibge, as.character(malha$code_muni))
stopifnot(!anyNA(ix), !anyDuplicated(ix))
malha <- malha[ix, ]
dest <- file.path(RAIZ, "execucoes_locais", "geometria")
dir.create(dest, recursive = TRUE, showWarnings = FALSE)
for (w in c("queen", "rook")) {
  lw <- spdep::nb2listw(spdep::poly2nb(malha, queen = w == "queen"), style = "W")
  ref <- readRDS(file.path(ESTIMACOES, "resultados/diagnostico_espacial", paste0("pesos_", w, ".rds")))
  delta <- max(abs(spdep::listw2mat(lw) - spdep::listw2mat(ref)))
  if (delta > 1e-12) stop("Rebuilt weights differ: ", w, " difference ", delta)
  saveRDS(lw, file.path(dest, paste0("pesos_", w, ".rds")))
  cat("Rebuilt weights match: ", w, " | maximum difference ", delta, "\n", sep = "")
}
