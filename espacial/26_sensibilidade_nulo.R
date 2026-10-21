# Sensibilidade do piso de D ao parametro de persistencia.
#
# Registro: a primeira versao de 25_nulo_de_D.R escalava o alvo de discordancia
# no sentido inverso, nao convergia na calibracao e caia em troca = 0,01. A
# auditoria de 10/09/2026 corrigiu o 25. Este script segue util como varredura:
# mostra quanto o piso depende da persistencia. Aqui a discordancia simulada da
# cadeia curta e levada a cadeia longa por sqrt(T_curto/T_real), que e o sentido
# correto, e comparada com a discordancia observada nas cadeias longas.
#
# Varre a persistencia em ordens de grandeza e reporta o piso em cada uma.

source(if (file.exists("../caminhos.R")) "../caminhos.R" else "caminhos.R")

AQUI <- ESPACIAL

D <- function(pip) {
  p <- mean(pip); K <- length(pip)
  if (p <= 0 || p >= 1) return(NA_real_)
  sd(pip) * sqrt((K - 1) / K) / sqrt(p * (1 - p))
}

cadeia_nula <- function(K, T_guardados, k_medio, troca) {
  inclusao <- numeric(K)
  atual <- sample(K, min(K, max(0, rpois(1, k_medio))))
  for (t in seq_len(T_guardados)) {
    if (runif(1) < troca) atual <- sample(K, min(K, max(0, rpois(1, k_medio))))
    if (length(atual) > 0L) inclusao[atual] <- inclusao[atual] + 1
  }
  inclusao / T_guardados
}

K <- 30L; T_CURTO <- 20000L; T_REAL <- 800000L
fator <- sqrt(T_CURTO / T_REAL)
REPS <- 200L

cat("=== Piso de D contra persistencia da cadeia ===\n")
cat("  Um valor de 'troca' menor significa cadeia mais persistente, menos\n")
cat("  informacao efetiva e piso mais alto.\n\n")
cat(sprintf("  %-10s %14s %14s %14s %14s\n",
            "troca", "p95 (MQO)", "discord. MQO", "p95 (SAR)", "discord. SAR"))

trocas <- c(0.001, 0.005, 0.01, 0.05, 0.2, 0.5, 1.0)
saida <- data.frame()
for (tr in trocas) {
  linha <- c(tr)
  for (km in c(13.0954, 7.1880)) {
    set.seed(20200731L)
    Ds <- replicate(REPS, D(cadeia_nula(K, T_CURTO, km, tr)))
    set.seed(4242L)
    disc <- mean(replicate(30L, {
      a <- cadeia_nula(K, T_CURTO, km, tr); b <- cadeia_nula(K, T_CURTO, km, tr)
      max(abs(a - b))
    })) * fator
    linha <- c(linha, quantile(Ds, 0.95) * fator, disc)
  }
  cat(sprintf("  %-10.4f %14.5f %14.5f %14.5f %14.5f\n",
              linha[1], linha[2], linha[3], linha[4], linha[5]))
  saida <- rbind(saida, data.frame(troca = linha[1], p95_mqo = linha[2],
                                   disc_mqo = linha[3], p95_sar = linha[4],
                                   disc_sar = linha[5]))
}

cat("\n=== Discordancia observada nas cadeias reais ===\n")
cat("  MQO: 0,0224 | SAR: 0,0125\n")
cat("  Escolhe-se a linha cuja discordancia simulada mais se aproxima.\n\n")

i_mqo <- which.min(abs(saida$disc_mqo - 0.0224))
i_sar <- which.min(abs(saida$disc_sar - 0.0125))
cat(sprintf("  MQO: troca %.4f reproduz discordancia %.5f, piso p95 = %.5f\n",
            saida$troca[i_mqo], saida$disc_mqo[i_mqo], saida$p95_mqo[i_mqo]))
cat(sprintf("  SAR: troca %.4f reproduz discordancia %.5f, piso p95 = %.5f\n",
            saida$troca[i_sar], saida$disc_sar[i_sar], saida$p95_sar[i_sar]))

cat(sprintf("\n  D observado no MQO: 0,3235 | razao com o piso: %.1f vezes\n",
            0.3235 / saida$p95_mqo[i_mqo]))
cat(sprintf("  D observado no SAR: 0,0973 | razao com o piso: %.1f vezes\n",
            0.0973 / saida$p95_sar[i_sar]))

cat("\n  Faixa do piso ao longo de toda a varredura:\n")
cat(sprintf("    MQO: %.5f a %.5f | SAR: %.5f a %.5f\n",
            min(saida$p95_mqo), max(saida$p95_mqo),
            min(saida$p95_sar), max(saida$p95_sar)))
cat(sprintf("  Mesmo no pior caso da varredura, o D do SAR fica %.1f vezes acima.\n",
            0.0973 / max(saida$p95_sar)))

write.csv(saida, file.path(AQUI, "saida", "sensibilidade_nulo.csv"), row.names = FALSE)
cat("\nGravado.\n")
