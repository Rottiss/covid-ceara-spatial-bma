# Gera a tabela da Secao 5.3 direto dos arquivos, com uma definicao so.
#
# HAVIA DUAS. O D pode ser calculado por replicacao e depois promediado, ou
# calculado uma vez sobre a PIP media das replicacoes. As duas respondem
# perguntas diferentes e dao numeros diferentes: no cenario realista de quatro
# verdadeiras, 0,863 contra 0,825. A primeira e a que descreve o que acontece
# numa amostra, que e a situacao de quem estima uma vez, e e a que usamos em
# todo o artigo. Esta rotina existe para que a tabela do texto nunca mais
# discorde do arquivo.

source(if (file.exists("../caminhos.R")) "../caminhos.R" else "caminhos.R")

AQUI <- ESPACIAL
D <- function(p) { m <- mean(p); K <- length(p)
  sd(p) * sqrt((K - 1) / K) / sqrt(m * (1 - m)) }

linha <- function(rotulo, r, rho) {
  # Os cenarios da primeira rodada nao guardaram o I de Moran no .rds; ele fica
  # so no log. Quem tem, usa; quem nao tem, entra como NA e o texto cita o log.
  mi <- if (is.null(r$mi)) NA_real_ else r$mi
  data.frame(cenario = rotulo, rho = rho, mi = mi,
             verd_ols = mean(r$ols[, r$idx]), verd_sar = mean(r$sar[, r$idx]),
             D_ols = mean(apply(r$ols, 1, D)), D_sar = mean(apply(r$sar, 1, D)),
             reps = nrow(r$ols), stringsAsFactors = FALSE)
}

L <- list()
sim <- readRDS(file.path(AQUI, "saida", "experimento_simulacao.rds"))
for (nm in names(sim)) L[[length(L) + 1]] <- linha(nm, sim[[nm]], sim[[nm]]$rho)

arqs <- c("cenario_realista.rds", "cenario_realista_oito.rds",
          "cenario_rho070_quatro.rds", "cenario_rho070_oito.rds")
for (a in arqs) {
  cam <- file.path(AQUI, "saida", a)
  if (!file.exists(cam)) next
  cr <- readRDS(cam)
  for (nm in names(cr$resultados))
    L[[length(L) + 1]] <- linha(sprintf("realista, %s verd.", nm),
                                cr$resultados[[nm]], cr$rho)
}
tab <- do.call(rbind, L)
tab$razao <- tab$D_sar / tab$D_ols
tab <- tab[order(-tab$razao), ]

cat("=== Secao 5.3, tabela para o texto ===\n\n")
cat("| cenário | $\\rho$ do DGP | I de Moran de $y$ | verdadeiras MQO -> SAR | $D$ MQO -> SAR | razão de $D$ |\n")
cat("|---|---|---|---|---|---|\n")
for (i in seq_len(nrow(tab))) {
  z <- tab[i, ]
  cat(sprintf("| %s | %s | %.3f | %.4f -> %.4f | %.3f -> %.3f | %.2f |\n",
              z$cenario, format(round(z$rho, 3), nsmall = 3), z$mi,
              z$verd_ols, z$verd_sar, z$D_ols, z$D_sar, z$razao))
}

real <- readRDS(file.path(AQUI, "saida", "grade", "casos_queen.rds"))
cat(sprintf("| **dados reais, casos** | 0,708 (posteriori) | %.3f | não se aplica | %.4f -> %.4f | %.2f |\n",
            real$mi_y, D(real$res$MQO$pip), D(real$res$SAR$pip),
            D(real$res$SAR$pip) / D(real$res$MQO$pip)))

cat(sprintf("\n  cenarios com processo gerador conhecido: %d | faixa da razao: %.2f a %.2f\n",
            nrow(tab), min(tab$razao), max(tab$razao)))
write.csv(tab, file.path(AQUI, "saida", "tabela_simulacoes.csv"), row.names = FALSE)
cat("\nGravado em saida/tabela_simulacoes.csv\n")
