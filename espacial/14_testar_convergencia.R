# Compara o tamanho medio posterior e as PIPs entre cadeias com diferentes
# comprimentos. O teste usa a configuracao de referencia de 1.000.000 de
# iteracoes e 200.000 de descarte para investigar diferencas de Monte Carlo.

source(if (file.exists("../caminhos.R")) "../caminhos.R" else "caminhos.R")

AQUI <- ESPACIAL
base <- readRDS(file.path(AQUI, "saida", "base_sar.rds"))
aux <- readRDS(file.path(AQUI, "saida", "grade_sar.rds"))
hyp <- readRDS(file.path(AQUI, "saida", "hyperg_aux.rds"))

X_total <- base$X; n <- base$n; covariaveis <- base$covariaveis
yc <- aux$yc; ta <- aux$ta
K <- length(covariaveis)
f21a <- hyp$f21a; a_par <- hyp$a_par; L <- hyp$L; J_MAX <- hyp$J_MAX

log_f21 <- function(z, k) {
  cc <- (f21a + k) / 2
  j_pico <- max(0, ceiling((z * a_par - cc) / max(1 - z, 1e-12)))
  J <- min(J_MAX, max(300L, as.integer(j_pico * 4)))
  lt <- c(0, z_seq <- log(z) * seq_len(J) + L[[k + 1]][seq_len(J)])
  m <- max(lt)
  m + log(sum(exp(lt - m)))
}

loglik <- function(ymy, k) {
  nl <- (1 - n) / 2 * log(ta)
  if (k == 0L) return(nl)
  nl + log(f21a - 2) - log(k + f21a - 2) + log_f21(1 - ymy / ta, k)
}

ymy_de <- function(indices) {
  if (length(indices) == 0L) return(ta)
  Xc <- X_total[, indices, drop = FALSE]
  Xc <- sweep(Xc, 2, colMeans(Xc), "-")
  Q <- qr.Q(qr(Xc))
  ta - sum((crossprod(Q, yc))^2)
}

lpm <- vapply(0:K, function(k) lgamma(1 + k) + lgamma(1 + K - k), numeric(1))
tipos_disp <- function(k) {
  t <- character(0)
  if (k < K) t <- c(t, "birth")
  if (k > 0L) t <- c(t, "death")
  if (k > 0L && k < K) t <- c(t, "move")
  t
}
log_q <- function(k, tipo) {
  -log(length(tipos_disp(k))) -
    log(switch(tipo, birth = K - k, death = k, move = k * (K - k)))
}

# Roda a cadeia registrando o tamanho medio acumulado em pontos de corte, para
# ver se ele ainda esta subindo.
mc3 <- function(n_sorteios, burn, semente, cortes) {
  set.seed(semente)
  atual <- sort(sample(seq_len(K), sample(0:K, 1)))
  post <- loglik(ymy_de(atual), length(atual)) + lpm[length(atual) + 1]
  inclusao <- numeric(K); guardados <- 0L; soma_k <- 0
  historico <- list()
  for (it in seq_len(n_sorteios)) {
    k <- length(atual)
    tipos <- tipos_disp(k)
    tipo <- tipos[sample.int(length(tipos), 1L)]
    fora <- setdiff(seq_len(K), atual)
    prop <- switch(tipo,
                   birth = sort(c(atual, fora[sample.int(length(fora), 1L)])),
                   death = atual[-sample.int(k, 1L)],
                   move  = sort(c(atual[-sample.int(k, 1L)],
                                  fora[sample.int(length(fora), 1L)])))
    post_p <- loglik(ymy_de(prop), length(prop)) + lpm[length(prop) + 1]
    inv <- switch(tipo, birth = "death", death = "birth", move = "move")
    la <- post_p - post + log_q(length(prop), inv) - log_q(k, tipo)
    if (is.finite(la) && log(runif(1)) < la) { atual <- prop; post <- post_p }
    if (it > burn) {
      guardados <- guardados + 1L
      if (length(atual) > 0L) inclusao[atual] <- inclusao[atual] + 1
      soma_k <- soma_k + length(atual)
      if (it %in% cortes) {
        historico[[as.character(it)]] <- list(
          iteracao = it, guardados = guardados,
          tamanho_medio = soma_k / guardados,
          pip = inclusao / guardados)
      }
    }
  }
  list(pip = inclusao / guardados, tamanho_medio = soma_k / guardados,
       historico = historico)
}

N <- 1000000L
BURN <- 200000L
cortes <- c(300000L, 400000L, 600000L, 800000L, 1000000L)

cat("=== Configuracao identica a do TCC ===\n")
cat("  2 cadeias, ", format(N, big.mark = ""), " iteracoes, ", format(BURN, big.mark = ""),
    " de descarte\n", sep = "")
cat("  (inventario_estimacoes.csv: cadeias = 2, iteracoes = 1e+06, burn_in = 2e+05)\n\n")

t0 <- proc.time()[3]
c1 <- mc3(N, BURN, 20200731L, cortes)
cat("  cadeia 1 em", round(proc.time()[3] - t0, 1), "s\n")
t0 <- proc.time()[3]
c2 <- mc3(N, BURN, 987654L, cortes)
cat("  cadeia 2 em", round(proc.time()[3] - t0, 1), "s\n\n")

pub <- read.csv(file.path(ESTIMACOES,
                          "resultados/estimacoes/principal/tabelas_tecnicas/resultado_casos.csv"),
                stringsAsFactors = FALSE)
pip_pub <- setNames(pub$PIP, pub$variavel)[covariaveis]
soma_pub <- sum(pip_pub)

cat("=== Tamanho medio do modelo ao longo da cadeia ===\n")
cat("  (a soma das PIP e o tamanho medio a posteriori)\n\n")
cat(sprintf("  %-12s %14s %14s %14s\n", "ate a iter.", "cadeia 1", "cadeia 2", "media"))
for (ct in as.character(cortes)) {
  t1 <- c1$historico[[ct]]$tamanho_medio
  t2 <- c2$historico[[ct]]$tamanho_medio
  cat(sprintf("  %-12s %14.4f %14.4f %14.4f\n", ct, t1, t2, (t1 + t2) / 2))
}
cat(sprintf("  %-12s %14s %14s %14.4f\n", "TCC", "", "", soma_pub))
cat(sprintf("  %-12s %14s %14s %14.4f\n", "priori E[k]", "", "", K / 2))

pip_media <- (c1$pip + c2$pip) / 2
comp <- data.frame(covariavel = covariaveis,
                   publicada = as.numeric(pip_pub),
                   minha = pip_media)
comp$dif <- comp$minha - comp$publicada
comp <- comp[order(-comp$publicada), ]

cat("\n=== PIP com 1.000.000 de iteracoes ===\n")
print(transform(comp, publicada = round(publicada, 4), minha = round(minha, 4),
                dif = round(dif, 4)), row.names = FALSE)

cat("\n  media da diferenca:", round(mean(comp$dif), 4), "\n")
cat("  negativas:", sum(comp$dif < 0), "de", K, "\n")
cat("  maior |diferenca|:", round(max(abs(comp$dif)), 4), "\n")
cat("  discordancia entre cadeias:", round(max(abs(c1$pip - c2$pip)), 4), "\n")
cat("\n  para comparar, com 400.000 iteracoes: media -0,0158, todas negativas\n")

write.csv(comp, file.path(AQUI, "saida", "convergencia_1m.csv"), row.names = FALSE)
cat("\nGravado.\n")
