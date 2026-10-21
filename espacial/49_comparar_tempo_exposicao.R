# Compara a sensibilidade de tempo de exposição rodada no 48 com duas referências.
#
#   1. A mesma especificação no TCC: BMS, uma cadeia, 100.000 de burn-in e
#      500.000 iterações. Mede quanto a configuração leve mudava as PIP.
#   2. As estimações principais (30 covariáveis) do 20, na mesma configuração de
#      cadeia. Mede o que os dias desde o primeiro registro mudam, sob MQO e SAR.
#
# Saída: saida/tempo_exposicao/comparacao.csv e um resumo no console.

source(if (file.exists("../caminhos.R")) "../caminhos.R" else "caminhos.R")

suppressPackageStartupMessages(library(BMS))

AQUI <- ESPACIAL
RAIZ <- ESTIMACOES
DEST <- file.path(AQUI, "saida", "tempo_exposicao")
DIAS <- "z_dias.primeiro.registro"

D <- function(p) { m <- mean(p); K <- length(p)
  sd(p) * sqrt((K - 1) / K) / sqrt(m * (1 - m)) }

tcc <- readRDS(file.path(RAIZ, "resultados/estimacoes/sensibilidade_tempo_exposicao/resultados.rds"))
linhas <- list()
for (d in c("casos", "obitos", "letalidade")) for (w in c("queen", "rook")) {
  arq <- file.path(DEST, sprintf("%s_%s.rds", d, w))
  if (!file.exists(arq)) { cat("falta:", arq, "\n"); next }
  te <- readRDS(arq)
  pr <- readRDS(file.path(AQUI, "saida", "grade", sprintf("%s_%s.rds", d, w)))
  cov <- te$tabela$covariavel
  mqo <- setNames(te$tabela$mqo, cov); sar <- setNames(te$tabela$sar, cov)

  cf <- coef(tcc[[d]]$cadeia_1, exact = FALSE)
  pip_tcc <- setNames(cf[, "PIP"], rownames(cf))[cov]
  stopifnot(!anyNA(pip_tcc))

  comuns <- setdiff(cov, DIAS)
  p_mqo <- setNames(pr$tabela$mqo, pr$tabela$covariavel)[comuns]
  p_sar <- setNames(pr$tabela$sar, pr$tabela$covariavel)[comuns]

  linhas[[length(linhas) + 1]] <- data.frame(
    desfecho = d, vizinhos = w,
    conv_mqo = te$res$MQO$conv, conv_sar = te$res$SAR$conv,
    dif_max_mqo_tcc = max(abs(mqo - pip_tcc)),
    acima_tcc = sum(pip_tcc >= 0.5), acima_mqo = sum(mqo >= 0.5), acima_sar = sum(sar >= 0.5),
    pip_dias_tcc = unname(pip_tcc[DIAS]), pip_dias_mqo = unname(mqo[DIAS]),
    pip_dias_sar = unname(sar[DIAS]),
    k_mqo = te$res$MQO$tam, k_sar = te$res$SAR$tam,
    D_mqo = D(mqo), D_sar = D(sar), razao_D = D(sar) / D(mqo),
    razao_D_principal = D(pr$res$SAR$pip) / D(pr$res$MQO$pip),
    acima_mqo_principal = sum(pr$res$MQO$pip >= 0.5),
    acima_sar_principal = sum(pr$res$SAR$pip >= 0.5),
    dif_max_mqo_vs_principal = max(abs(mqo[comuns] - p_mqo)),
    dif_max_sar_vs_principal = max(abs(sar[comuns] - p_sar)),
    stringsAsFactors = FALSE)

  cat(sprintf("\n=== %s x %s ===\n", d, w))
  cat(sprintf("  discordância entre cadeias: MQO %.4f | SAR %.4f\n", te$res$MQO$conv, te$res$SAR$conv))
  cat(sprintf("  MQO contra TCC (1 cadeia, 500.000): maior diferença %.4f\n", max(abs(mqo - pip_tcc))))
  cat(sprintf("  PIP dos dias: TCC %.4f | MQO %.4f | SAR %.4f\n", pip_tcc[DIAS], mqo[DIAS], sar[DIAS]))
  cat("  acima de 0,50, TCC:", paste(cov[pip_tcc >= 0.5], collapse = ", "), "\n")
  cat("  acima de 0,50, MQO:", paste(cov[mqo >= 0.5], collapse = ", "), "\n")
  cat("  acima de 0,50, SAR:", if (any(sar >= 0.5)) paste(cov[sar >= 0.5], collapse = ", ") else "nenhuma", "\n")
  cat(sprintf("  razão de D: %.3f (principal %.3f)\n", D(sar) / D(mqo),
              D(pr$res$SAR$pip) / D(pr$res$MQO$pip)))
}
tab <- do.call(rbind, linhas)
write.csv(tab, file.path(DEST, "comparacao.csv"), row.names = FALSE)
cat("\ngravado:", file.path(DEST, "comparacao.csv"), "\n")
