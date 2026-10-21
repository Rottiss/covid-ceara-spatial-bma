# Compara Queen e Rook na escolha da especificacao (33) e no fator de Bayes
# agregado (43 para Queen, 60 para Rook).
#
# A Secao 4.4 do artigo hoje reporta so a Queen. Este script mede quanto
# P(SAR | y), o rho a posteriori e o fator de Bayes agregado dependem da matriz.
#
# Saida: saida/comparacao_vizinhanca_especificacao.csv

source(if (file.exists("../caminhos.R")) "../caminhos.R" else "caminhos.R")

AQUI <- ESPACIAL
DESF <- c("casos", "obitos", "letalidade")

resumo <- function(d, viz) {
  arq <- file.path(AQUI, "saida", sprintf("especificacao_%s_%s.rds", d, viz))
  if (!file.exists(arq)) return(NULL)
  e <- readRDS(arq)
  w <- (e$a$w_rho + e$b$w_rho) / 2; w <- w / sum(w); ac <- cumsum(w)
  csv <- file.path(AQUI, "saida", if (viz == "queen")
    sprintf("sensibilidade_prior_rho_%s.csv", d) else
    sprintf("sensibilidade_prior_rho_%s_rook.csv", d))
  pr <- if (file.exists(csv)) read.csv(csv) else NULL
  data.frame(
    desfecho = d, vizinhos = viz,
    p_sar_cadeia = (e$a$p_sar + e$b$p_sar) / 2,
    visitas_mqo = e$a$n_ols + e$b$n_ols,
    rho_media = sum(w * e$rho),
    rho_inf = e$rho[which.max(ac >= 0.025)],
    rho_sup = e$rho[which.max(ac >= 0.975)],
    log_bf_agregado = if (is.null(pr)) NA_real_ else pr$log_bf_agregado_mc[pr$a0 == 1.01],
    p_sar_agregada = if (is.null(pr)) NA_real_ else pr$p_sar_agregada_mc[pr$a0 == 1.01],
    log_bf_cond_min = if (is.null(pr)) NA_real_ else min(pr$log_bf_cond_min),
    ess_mqo = if (is.null(pr)) NA_real_ else min(pr$ess_mqo),
    stringsAsFactors = FALSE)
}

linhas <- list()
for (d in DESF) for (v in c("queen", "rook")) linhas[[length(linhas) + 1]] <- resumo(d, v)
tab <- do.call(rbind, linhas)
arq <- file.path(AQUI, "saida", "comparacao_vizinhanca_especificacao.csv")
write.csv(tab, arq, row.names = FALSE)
op <- options(width = 200)
print(tab, digits = 4, row.names = FALSE)
options(op)
cat("\nGravado em", arq, "\n")
for (d in DESF) {
  q <- tab[tab$desfecho == d & tab$vizinhos == "queen", ]
  r <- tab[tab$desfecho == d & tab$vizinhos == "rook", ]
  if (nrow(q) == 1 && nrow(r) == 1)
    cat(sprintf("%-11s diferenca Queen-Rook: rho %+.4f | P(SAR) cadeia %+.4f | log BF %+.3f\n",
                d, q$rho_media - r$rho_media, q$p_sar_cadeia - r$p_sar_cadeia,
                q$log_bf_agregado - r$log_bf_agregado))
}
