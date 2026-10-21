# Auditoria: a coluna "Especificacoes >= 0,50" das Tabelas 5 a 7 do TCC.
#
# O TCC afirma, por covariavel, em quantas das sete especificacoes (oito em
# obitos) a PIP fica em 0,50 ou acima. Esta contagem nunca foi conferida contra
# as saidas. Le apenas os CSVs de comparacao do COVID_BMA.

source(if (file.exists("../caminhos.R")) "../caminhos.R" else "caminhos.R")

RAIZ <- file.path(ESTIMACOES, "resultados/estimacoes")

# A familia de escalas usa nomes proprios (casos_nivel, obitos_log1p, ...),
# por isso entra como par (pasta, sufixo do arquivo).
familias <- list(
  c("principal", "%s"),
  c("ampliado", "%s"),
  c("sensibilidade_escalas", "%s_nivel"),
  c("sensibilidade_sem_fortaleza", "%s"),
  c("sensibilidade_pop_10mil", "%s"),
  c("sensibilidade_pop_20mil", "%s"),
  c("sensibilidade_tempo_exposicao", "%s")
)
extra_obitos <- list(c("sensibilidade_escalas", "obitos_log1p"))

pip_familia <- function(par, d) {
  arq <- file.path(RAIZ, par[1], "tabelas_tecnicas",
                   paste0("comparacao_", sprintf(par[2], d), ".csv"))
  if (!file.exists(arq)) return(NULL)
  x <- read.csv(arq, stringsAsFactors = FALSE)
  setNames(x$PIP_media, x$variavel)
}

esperado <- list(
  casos = c(z_idosos = 7, z_leitos1k = 6, metrop = 7, semiarido = 6,
            z_ideb.fund = 7, z_int.circ = 5),
  obitos = c(z_idosos = 8, z_ex.pobr = 1),
  letalidade = c(metrop = 7, z_eleitas.fem = 6, z_imuni = 7, z_int.resp = 5,
                 aliadogov = 6, z_bolsaf = 5, z_ideb.fund = 4, z_sus1k = 4,
                 z_leitos1k = 2, z_ex.pobr = 2, z_area = 2, z_emprego = 3,
                 z_ivs.infra = 2)
)

cat("pastas extras encontradas (candidatas a log1p de obitos):",
    paste(extra_obitos, collapse = ", "), "\n\n")

div <- 0L
for (d in names(esperado)) {
  fams <- familias
  if (d == "obitos") fams <- c(familias, extra_obitos)
  tab <- list()
  for (f in fams) {
    p <- pip_familia(f, d)
    if (!is.null(p)) tab[[paste(f, collapse = "/")]] <- p
  }
  cat("===", d, "| familias com saida:", length(tab), "\n")
  cat("    ", paste(names(tab), collapse = ", "), "\n")
  for (v in names(esperado[[d]])) {
    cont <- sum(vapply(tab, function(p) isTRUE(p[[v]] >= 0.5), logical(1)), na.rm = TRUE)
    marca <- if (cont != esperado[[d]][[v]]) "  <-- DIVERGE" else ""
    cat(sprintf("    %-16s TCC %d de %d | recontado %d de %d%s\n",
                v, esperado[[d]][[v]], length(tab), cont, length(tab), marca))
    if (cont != esperado[[d]][[v]]) div <- div + 1L
  }
  cat("\n")
}
cat("=== divergencias de contagem:", div, "===\n")
