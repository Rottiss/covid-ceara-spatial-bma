# Summarize Queen-weight results by outcome, including PIP rankings,
# normalized dispersion, mean model size, and outcome Moran's I.
# The comparisons are descriptive, not tests of a causal explanation.

source(if (file.exists("../caminhos.R")) "../caminhos.R" else "caminhos.R")

AQUI <- ESPACIAL
D <- function(p) { m <- mean(p); K <- length(p)
  sd(p) * sqrt((K - 1) / K) / sqrt(m * (1 - m)) }

arqs <- list.files(file.path(AQUI, "saida", "grade"), pattern = "queen\\.rds$",
                   full.names = TRUE)

for (a in arqs) {
  g <- readRDS(a)
  t <- g$tabela[order(-g$tabela$mqo), ]
  cat(sprintf("\n=== %s, Queen ===  (I de Moran do desfecho: %.4f)\n",
              g$desfecho, g$mi_y))
  cat(sprintf("  tamanho medio %.2f -> %.2f | D %.4f -> %.4f | linha de base %.4f -> %.4f\n",
              sum(g$res$MQO$pip), sum(g$res$SAR$pip),
              D(g$res$MQO$pip), D(g$res$SAR$pip),
              mean(g$res$MQO$pip), mean(g$res$SAR$pip)))
  cat(sprintf("  %-16s %8s %8s %8s\n", "covariavel", "MQO", "SAR", "dif"))
  for (i in seq_len(min(10, nrow(t)))) {
    marca <- ""
    if (t$mqo[i] >= 0.5 && t$sar[i] < 0.5) marca <- "  <- cai abaixo de 0,50"
    if (t$mqo[i] < 0.5 && t$sar[i] >= 0.5) marca <- "  <- sobe acima de 0,50"
    cat(sprintf("  %-16s %8.4f %8.4f %8.4f%s\n",
                t$covariavel[i], t$mqo[i], t$sar[i], t$dif[i], marca))
  }
  entram <- t$covariavel[t$mqo < 0.5 & t$sar >= 0.5]
  if (length(entram)) cat("  entram no SAR:", paste(entram, collapse = ", "), "\n")
}

# ---------------------------- teste da explicacao
todos <- list.files(file.path(AQUI, "saida", "grade"), pattern = "\\.rds$",
                    full.names = TRUE)
dd <- do.call(rbind, lapply(todos, function(a) {
  g <- readRDS(a)
  data.frame(desfecho = g$desfecho, vizinhos = g$vizinhos, mi_y = g$mi_y,
             k_mqo = sum(g$res$MQO$pip), D_mqo = D(g$res$MQO$pip),
             razao_D = D(g$res$SAR$pip) / D(g$res$MQO$pip),
             n_marginais = sum(g$res$MQO$pip >= 0.5 & g$res$MQO$pip < 0.7),
             stringsAsFactors = FALSE)
}))

cat("\n\n=== A severidade acompanha o que? ===\n")
print(transform(dd, mi_y = round(mi_y, 3), k_mqo = round(k_mqo, 2),
                D_mqo = round(D_mqo, 4), razao_D = round(razao_D, 3)),
      row.names = FALSE)

cat("\n  correlacoes com a razao de D (6 pontos, 3 independentes):\n")
cat(sprintf("    I de Moran do desfecho : %+.4f\n", cor(dd$mi_y, dd$razao_D)))
cat(sprintf("    tamanho medio sob MQO  : %+.4f\n", cor(dd$k_mqo, dd$razao_D)))
cat(sprintf("    D sob MQO              : %+.4f\n", cor(dd$D_mqo, dd$razao_D)))
cat(sprintf("    n de PIP entre 0,5 e 0,7: %+.4f\n", cor(dd$n_marginais, dd$razao_D)))

cat("\n  Com tres desfechos apenas, isto e sugestivo e nao conclusivo.\n")
cat("  Serve para orientar a extensao, nao para sustentar afirmacao.\n")

write.csv(dd, file.path(AQUI, "saida", "severidade.csv"), row.names = FALSE)
