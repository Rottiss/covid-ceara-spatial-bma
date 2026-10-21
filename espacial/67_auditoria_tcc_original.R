# Auditoria: as PIP das Tabelas 5 a 7 do TCC batem com as saidas do COVID_BMA?
#
# Este e o unico bloco do TCC que o projeto MC3 nunca conferiu: os numeros dos
# capitulos originais, produzidos pelo pacote BMS. A metodologia do TCC diz que
# as duas cadeias sao combinadas ponderando pelo numero de iteracoes retidas.
#
# Le apenas. Nao grava nada em COVID_BMA.

source(if (file.exists("../caminhos.R")) "../caminhos.R" else "caminhos.R")

suppressPackageStartupMessages(library(BMS))
RAIZ <- file.path(ESTIMACOES, "resultados/estimacoes/principal")
r <- readRDS(file.path(RAIZ, "resultados.rds"))

pip_combinada <- function(par) {
  saidas <- lapply(par, function(m) {
    co <- coef(m, order.by.pip = FALSE, include.constant = FALSE)
    setNames(co[, "PIP"], rownames(co))
  })
  n <- vapply(par, function(m) as.numeric(m$info$iter), numeric(1))
  nomes <- names(saidas[[1]])
  stopifnot(identical(nomes, names(saidas[[2]])))
  (saidas[[1]][nomes] * n[1] + saidas[[2]][nomes] * n[2]) / sum(n)
}

# Valores transcritos das Tabelas 5, 6 e 7 do TCC (Overleaf-TCC-MC3).
esperado <- list(
  casos = c(z_idosos = 0.9442, z_leitos1k = 0.8278, metrop = 0.6894,
            semiarido = 0.6293, z_ideb.fund = 0.6234, z_int.circ = 0.5790),
  obitos = c(z_idosos = 0.9737, z_ex.pobr = 0.5018),
  letalidade = c(metrop = 0.8089, z_eleitas.fem = 0.7473, z_imuni = 0.6551,
                 z_int.resp = 0.5992, aliadogov = 0.5991, z_bolsaf = 0.5916,
                 z_ideb.fund = 0.5559, z_sus1k = 0.5442, z_leitos1k = 0.5401,
                 z_ex.pobr = 0.5079, z_area = 0.5045, z_emprego = 0.5027,
                 z_ivs.infra = 0.5003)
)

cat("=== PIP do TCC contra a saida do BMS, duas cadeias ponderadas ===\n")
div <- 0L
for (d in names(esperado)) {
  p <- pip_combinada(r[[d]])
  cat("\n ", d, "| iteracoes retidas:",
      paste(vapply(r[[d]], function(m) as.numeric(m$info$iter), numeric(1)), collapse = " e "), "\n")
  for (v in names(esperado[[d]])) {
    if (!v %in% names(p)) { cat("    AUSENTE na saida:", v, "\n"); div <- div + 1L; next }
    dif <- abs(p[[v]] - esperado[[d]][[v]])
    marca <- if (dif > 5e-5) "  <-- DIVERGE" else ""
    cat(sprintf("    %-16s TCC %.4f | BMS %.6f | dif %.2e%s\n",
                v, esperado[[d]][[v]], p[[v]], dif, marca))
    if (dif > 5e-5) div <- div + 1L
  }
  acima <- sort(p[p >= 0.5], decreasing = TRUE)
  cat("    acima de 0,50 na saida:", length(acima),
      "| na tabela do TCC:", length(esperado[[d]]), "\n")
  if (length(acima) != length(esperado[[d]]))
    cat("    <-- CONTAGEM DIVERGE:", paste(names(acima), collapse = ", "), "\n")
}
cat("\n=== divergencias:", div, "===\n")
