# Estatistica D em todos os cenarios simulados, para comparacao direta com os
# dados reais.

source(if (file.exists("../caminhos.R")) "../caminhos.R" else "caminhos.R")

AQUI <- ESPACIAL
D <- function(p) { m <- mean(p); K <- length(p)
  sd(p) * sqrt((K - 1) / K) / sqrt(m * (1 - m)) }

cat("=== Cenarios com rho = 0 e o controle ===\n")
sim <- readRDS(file.path(AQUI, "saida", "experimento_simulacao.rds"))
for (nm in names(sim)) {
  r <- sim[[nm]]
  do <- mean(apply(r$ols, 1, D)); ds <- mean(apply(r$sar, 1, D))
  cat(sprintf("  %-10s rho %.2f | D %.4f -> %.4f (%3.0f%%) | verdadeiras %.4f -> %.4f | falsas %.4f -> %.4f\n",
              nm, r$rho, do, ds, 100 * ds / do,
              mean(r$ols[, r$idx]), mean(r$sar[, r$idx]),
              mean(r$ols[, -r$idx]), mean(r$sar[, -r$idx])))
}

cat("\n=== Cenario realista ===\n")
cr <- readRDS(file.path(AQUI, "saida", "cenario_realista.rds"))
cat(sprintf("  rho do DGP: %.4f\n", cr$rho))
for (nm in names(cr$resultados)) {
  r <- cr$resultados[[nm]]
  do <- mean(apply(r$ols, 1, D)); ds <- mean(apply(r$sar, 1, D))
  cat(sprintf("  %-8s (%d verd.) I de Moran de y %.3f | D %.4f -> %.4f (%3.0f%%)\n",
              nm, length(r$idx), r$mi, do, ds, 100 * ds / do))
  cat(sprintf("           verdadeiras %.4f -> %.4f | falsas %.4f -> %.4f\n",
              mean(r$ols[, r$idx]), mean(r$sar[, r$idx]),
              mean(r$ols[, -r$idx]), mean(r$sar[, -r$idx])))
  cat(sprintf("           linha de nao informacao: MQO %.4f | SAR %.4f\n",
              mean(r$tam[, 1]) / 30, mean(r$tam[, 2]) / 30))
  cat(sprintf("           razao verdadeiras/linha de base: MQO %.2f | SAR %.2f\n",
              mean(r$ols[, r$idx]) / (mean(r$tam[, 1]) / 30),
              mean(r$sar[, r$idx]) / (mean(r$tam[, 2]) / 30)))
}

cat("\n=== Comparacao com os dados reais ===\n")
real <- read.csv(file.path(AQUI, "saida", "comparacao_1milhao.csv"),
                 stringsAsFactors = FALSE)
cat(sprintf("  casos: D %.4f -> %.4f (%.0f%%)\n",
            D(real$mqo), D(real$sar), 100 * D(real$sar) / D(real$mqo)))
cat(sprintf("  maior PIP sob SAR nos dados reais: %.4f (%s)\n",
            max(real$sar), real$covariavel[which.max(real$sar)]))
cat(sprintf("  linha de nao informacao nos dados reais sob SAR: %.4f\n",
            mean(real$sar)))
cat(sprintf("  razao maior PIP / linha de base: %.2f\n",
            max(real$sar) / mean(real$sar)))

cat("\n  Nas simulacoes com processo gerador conhecido e dependencia espacial realista, as\n")
cat("  verdadeiras ficam entre 1,5 e 2,9 vezes a linha de base sob SAR. Nos dados\n")
cat("  reais, nem a maior chega a 1,6 vezes.\n")
