source("../caminhos.R")
if (!requireNamespace("spdep", quietly = TRUE)) stop("Install spdep.")
cols <- c(casos = "z_log_casos100k", obitos = "z_log_obitos100k", letalidade = "z_logit_letal")
ids <- NULL
covs <- NULL
for (d in names(cols)) {
  x <- read.csv(file.path(ESTIMACOES, "resultados/estimacoes/matrizes_estimacao/principal", paste0(d, ".csv")),
                colClasses = c(codigo_ibge = "character"))
  stopifnot(nrow(x) == 184, ncol(x) == 33, !anyNA(x), !anyDuplicated(x$codigo_ibge))
  cc <- setdiff(names(x), c("codigo_ibge", "municipio", cols[[d]]))
  stopifnot(length(cc) == 30, all(is.finite(as.matrix(x[, c(cols[[d]], cc)]))))
  if (is.null(ids)) { ids <- x$codigo_ibge; covs <- cc }
  stopifnot(identical(ids, x$codigo_ibge), identical(cc, covs))
  for (dummy in c("metrop", "semiarido", "aliadogov")) stopifnot(all(x[[dummy]] %in% 0:1))
  cat("OK matrix: ", d, "\n", sep = "")
  for (w in c("queen", "rook")) {
    res <- readRDS(file.path(ESPACIAL, "saida/tres_vias", paste0(d, "_", w, ".rds")))
    stopifnot(res$desfecho == d, res$vizinhos == w,
              res$n_mc == 1000000, res$burn == 200000,
              identical(res$tabela$covariavel, cc),
              all(is.finite(res$p_esp)), all(res$p_esp >= 0), all(res$p_esp <= 1),
              abs(sum(res$p_esp) - 1) < 1e-12,
              max(abs(res$p_esp - (res$p_esp_a + res$p_esp_b) / 2)) < 1e-12)
    vals <- as.matrix(res$tabela[, -1])
    stopifnot(all(is.na(vals) | (vals >= 0 & vals <= 1)))
  }
}
for (w in c("queen", "rook")) {
  lw <- readRDS(file.path(ESTIMACOES, "resultados/diagnostico_espacial", paste0("pesos_", w, ".rds")))
  W <- spdep::listw2mat(lw)
  stopifnot(identical(dim(W), c(184L, 184L)), all(diag(W) == 0),
            max(abs(rowSums(W) - 1)) < 1e-12, all(W >= 0))
  cat("OK weights: ", w, "\n", sep = "")
}
cat("OK matrix, weight, and result contracts. This does not establish convergence.\n")
