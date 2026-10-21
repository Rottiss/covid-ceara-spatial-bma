# Compare PIP dispersion, model size, and rankings across the empirical grid.

source(if (file.exists("../caminhos.R")) "../caminhos.R" else "caminhos.R")

AQUI <- ESPACIAL

D <- function(p) {
  m <- mean(p); K <- length(p)
  sd(p) * sqrt((K - 1) / K) / sqrt(m * (1 - m))
}

arqs <- list.files(file.path(AQUI, "saida", "grade"), pattern = "\\.rds$",
                   full.names = TRUE)
if (length(arqs) == 0L) stop("Nenhum resultado da grade ainda.")

linhas <- list()
for (a in arqs) {
  g <- readRDS(a)
  t <- g$tabela
  topo_mqo <- head(t[order(-t$mqo), "covariavel"], 5)
  topo_sar <- head(t[order(-t$sar), "covariavel"], 5)
  acima_mqo <- t$covariavel[t$mqo >= 0.5]
  acima_sar <- t$covariavel[t$sar >= 0.5]

  cat(sprintf("=== %s x %s ===\n", g$desfecho, g$vizinhos))
  cat(sprintf("  I de Moran do desfecho: %.4f\n", g$mi_y))
  cat(sprintf("  tamanho medio: %.3f -> %.3f\n",
              sum(g$res$MQO$pip), sum(g$res$SAR$pip)))
  cat(sprintf("  D: %.4f -> %.4f  (%.0f%% do original)\n",
              D(g$res$MQO$pip), D(g$res$SAR$pip),
              100 * D(g$res$SAR$pip) / D(g$res$MQO$pip)))
  cat(sprintf("  acima de 0,50 no MQO (%d): %s\n", length(acima_mqo),
              ifelse(length(acima_mqo), paste(acima_mqo, collapse = ", "), "nenhuma")))
  cat(sprintf("  acima de 0,50 no SAR (%d): %s\n", length(acima_sar),
              ifelse(length(acima_sar), paste(acima_sar, collapse = ", "), "nenhuma")))
  cat(sprintf("  interseccao dos cinco maiores: %d de 5\n",
              length(intersect(topo_mqo, topo_sar))))
  cat(sprintf("  correlacao de Spearman entre as PIP: %.4f\n\n",
              cor(t$mqo, t$sar, method = "spearman")))

  linhas[[length(linhas) + 1]] <- data.frame(
    desfecho = g$desfecho, vizinhos = g$vizinhos, mi_y = g$mi_y,
    k_mqo = sum(g$res$MQO$pip), k_sar = sum(g$res$SAR$pip),
    D_mqo = D(g$res$MQO$pip), D_sar = D(g$res$SAR$pip),
    razao_D = D(g$res$SAR$pip) / D(g$res$MQO$pip),
    n_acima_mqo = length(acima_mqo), n_acima_sar = length(acima_sar),
    interseccao_topo5 = length(intersect(topo_mqo, topo_sar)),
    spearman = cor(t$mqo, t$sar, method = "spearman"),
    conv_mqo = g$res$MQO$conv, conv_sar = g$res$SAR$conv,
    stringsAsFactors = FALSE)
}

tab <- do.call(rbind, linhas)
tab <- tab[order(tab$desfecho, tab$vizinhos), ]

cat("=== Quadro comparativo ===\n")
print(transform(tab[, c("desfecho", "vizinhos", "mi_y", "k_mqo", "k_sar",
                        "D_mqo", "D_sar", "razao_D", "n_acima_mqo",
                        "n_acima_sar", "interseccao_topo5", "spearman")],
                mi_y = round(mi_y, 3), k_mqo = round(k_mqo, 2),
                k_sar = round(k_sar, 2), D_mqo = round(D_mqo, 4),
                D_sar = round(D_sar, 4), razao_D = round(razao_D, 3),
                spearman = round(spearman, 3)), row.names = FALSE)

cat("\n=== Descriptive comparison rules ===\n")
cat("  Report a D ratio below 0.50 first; otherwise report top-five overlap\n")
cat("  of at most two variables. All remaining rows are labelled Other.\n")
cat("  These rules are not formal tests of model performance.\n")
cat("\n")
for (i in seq_len(nrow(tab))) {
  padrao <- if (tab$razao_D[i] < 0.5) "D ratio <0.50" else
    if (tab$interseccao_topo5[i] <= 2) "Top-five overlap <=2" else "Other"
  cat(sprintf("  %-11s %-6s razao de D %.2f | interseccao %d/5 -> %s\n",
              tab$desfecho[i], tab$vizinhos[i], tab$razao_D[i],
              tab$interseccao_topo5[i], padrao))
}

cat("\n=== Convergencia ===\n")
cat("  maior discordancia entre cadeias: MQO", round(max(tab$conv_mqo), 4),
    "| SAR", round(max(tab$conv_sar), 4), "\n")

write.csv(tab, file.path(AQUI, "saida", "resumo_grade.csv"), row.names = FALSE)
cat("\nGravado em saida/resumo_grade.csv\n")
