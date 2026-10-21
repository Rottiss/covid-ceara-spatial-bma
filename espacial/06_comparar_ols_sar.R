# Compare OLS/SAR PIPs for the 30 covariates using MC3.
#
# Scope:
#   - priori uniforme sobre modelos e g fixo (UIP), NAO a configuracao do TCC
#     (hyper-g UIP com Binomial-Beta). A comparacao aqui e MQO contra SAR sob a
#     MESMA priori, que e o que isola o efeito do termo espacial;
#   - por isso as PIP de MQO abaixo nao sao, e nao devem ser lidas como, as PIP
#     publicadas no capitulo 4;
#   - desfecho: casos. W: Queen.

source(if (file.exists("../caminhos.R")) "../caminhos.R" else "caminhos.R")

AQUI <- ESPACIAL
base <- readRDS(file.path(AQUI, "saida", "base_sar.rds"))
aux <- readRDS(file.path(AQUI, "saida", "grade_sar.rds"))

y <- base$y; X_total <- base$X; n <- base$n; covariaveis <- base$covariaveis
grade <- aux$grade; ta <- aux$ta; tb <- aux$tb; tc <- aux$tc
yc <- aux$yc; Wyc <- aux$Wyc
K <- length(covariaveis)
g_uip <- n

estatisticas <- function(indices) {
  if (length(indices) == 0L) {
    eo <- yc; ed <- Wyc
  } else {
    Xc <- X_total[, indices, drop = FALSE]
    Xc <- sweep(Xc, 2, colMeans(Xc), "-")
    Q <- qr.Q(qr(Xc))
    eo <- yc - Q %*% crossprod(Q, yc)
    ed <- Wyc - Q %*% crossprod(Q, Wyc)
  }
  eo <- as.numeric(eo); ed <- as.numeric(ed)
  c(k = length(indices), a = sum(eo^2), b = sum(eo * ed), c = sum(ed^2))
}

lm_sar <- function(est, g = g_uip) {
  rss <- est["a"] - 2 * grade$rho * est["b"] + grade$rho2 * est["c"]
  tss <- ta - 2 * grade$rho * tb + grade$rho2 * tc
  lf <- grade$ld_mais_prior -
    ((n - 1) / 2) * log((g / (1 + g)) * rss + (1 / (1 + g)) * tss)
  m <- max(lf)
  unname(-(est["k"] / 2) * log(1 + g) + m +
           log(grade$h / 3 * sum(grade$pesos * exp(lf - m))))
}

# Com rho = 0 nao ha integracao: e a log-marginal de minimos quadrados,
# validada na Etapa 1.
lm_ols <- function(est, g = g_uip) {
  unname(-(est["k"] / 2) * log(1 + g) -
           ((n - 1) / 2) * log(ta - (g / (1 + g)) * (ta - est["a"])))
}

tipos_disponiveis <- function(k, K) {
  tipos <- character(0)
  if (k < K) tipos <- c(tipos, "birth")
  if (k > 0L) tipos <- c(tipos, "death")
  if (k > 0L && k < K) tipos <- c(tipos, "move")
  tipos
}
log_q <- function(k, K, tipo) {
  -log(length(tipos_disponiveis(k, K))) -
    log(switch(tipo, birth = K - k, death = k, move = k * (K - k)))
}

mc3 <- function(log_marginal, n_sorteios, burn, semente) {
  set.seed(semente)
  atual <- sort(sample(seq_len(K), sample(0:K, 1)))
  lm_atual <- log_marginal(estatisticas(atual))
  aceitos <- 0L; inclusao <- numeric(K); guardados <- 0L
  tamanhos <- integer(0); unicos <- new.env(hash = TRUE, parent = emptyenv())

  for (it in seq_len(n_sorteios)) {
    k <- length(atual)
    tipos <- tipos_disponiveis(k, K)
    tipo <- tipos[sample.int(length(tipos), 1L)]
    fora <- setdiff(seq_len(K), atual)
    proposto <- switch(
      tipo,
      birth = sort(c(atual, fora[sample.int(length(fora), 1L)])),
      death = atual[-sample.int(k, 1L)],
      move  = sort(c(atual[-sample.int(k, 1L)], fora[sample.int(length(fora), 1L)]))
    )
    lm_prop <- log_marginal(estatisticas(proposto))
    inverso <- switch(tipo, birth = "death", death = "birth", move = "move")
    log_alfa <- lm_prop - lm_atual +
      log_q(length(proposto), K, inverso) - log_q(k, K, tipo)
    if (log(runif(1)) < log_alfa) {
      atual <- proposto; lm_atual <- lm_prop; aceitos <- aceitos + 1L
    }
    if (it > burn) {
      guardados <- guardados + 1L
      if (length(atual) > 0L) inclusao[atual] <- inclusao[atual] + 1
      tamanhos <- c(tamanhos, length(atual))
      chave <- if (length(atual) == 0L) "vazio" else paste(atual, collapse = ",")
      assign(chave, TRUE, envir = unicos)
    }
  }
  list(pip = inclusao / guardados, taxa = aceitos / n_sorteios,
       tamanho_medio = mean(tamanhos), unicos = length(ls(unicos)))
}

N_SORTEIOS <- 300000L
BURN <- 60000L

cat("=== Configuracao ===\n")
cat("  desfecho: casos | W: Queen | K =", K, "| espaco de modelos: 2^", K, "\n", sep = "")
cat("  sorteios:", N_SORTEIOS, "| descarte:", BURN, "\n")
cat("  g = UIP fixo | priori de modelo: uniforme\n\n")

resultados <- list()
for (rotulo in c("MQO", "SAR")) {
  f <- if (rotulo == "MQO") lm_ols else lm_sar
  cat("Rodando", rotulo, "... ")
  t0 <- proc.time()[3]
  a <- mc3(f, N_SORTEIOS, BURN, 20200731L)
  b <- mc3(f, N_SORTEIOS, BURN, 987654L)
  cat(round(proc.time()[3] - t0, 1), "s\n")
  cat("  aceitacao:", round(a$taxa, 3), "e", round(b$taxa, 3),
      "| tamanho medio:", round(a$tamanho_medio, 2), "e", round(b$tamanho_medio, 2),
      "| modelos unicos:", a$unicos, "e", b$unicos, "\n")
  cat("  maior diferenca de PIP entre as duas cadeias:",
      signif(max(abs(a$pip - b$pip)), 3), "\n")
  resultados[[rotulo]] <- list(a = a, b = b, pip = (a$pip + b$pip) / 2,
                               conv = max(abs(a$pip - b$pip)))
}

comparacao <- data.frame(
  covariavel = covariaveis,
  pip_mqo = resultados$MQO$pip,
  pip_sar = resultados$SAR$pip,
  stringsAsFactors = FALSE
)
comparacao$diferenca <- comparacao$pip_sar - comparacao$pip_mqo
comparacao <- comparacao[order(-comparacao$pip_sar), ]

cat("\n=== PIP: MQO vs. SAR, 30 covariaveis, desfecho casos ===\n")
print(transform(comparacao,
                pip_mqo = round(pip_mqo, 4),
                pip_sar = round(pip_sar, 4),
                diferenca = round(diferenca, 4)), row.names = FALSE)

acima_mqo <- comparacao$covariavel[comparacao$pip_mqo >= 0.5]
acima_sar <- comparacao$covariavel[comparacao$pip_sar >= 0.5]

cat("\n=== Quem passa de 0,50 ===\n")
cat("  MQO:", if (length(acima_mqo)) paste(acima_mqo, collapse = ", ") else "nenhuma", "\n")
cat("  SAR:", if (length(acima_sar)) paste(acima_sar, collapse = ", ") else "nenhuma", "\n")
cat("  so no MQO:", paste(setdiff(acima_mqo, acima_sar), collapse = ", "), "\n")
cat("  so no SAR:", paste(setdiff(acima_sar, acima_mqo), collapse = ", "), "\n")

cat("\n=== Magnitude das diferencas ===\n")
cat("  maior diferenca absoluta:", round(max(abs(comparacao$diferenca)), 4),
    "em", comparacao$covariavel[which.max(abs(comparacao$diferenca))], "\n")
cat("  diferencas acima de 0,10:",
    sum(abs(comparacao$diferenca) > 0.10), "de", K, "\n")
cat("  diferencas acima de 0,05:",
    sum(abs(comparacao$diferenca) > 0.05), "de", K, "\n")
cat("  ruido de Monte Carlo (maior discordancia entre cadeias): MQO",
    signif(resultados$MQO$conv, 3), "| SAR", signif(resultados$SAR$conv, 3), "\n")

write.csv(comparacao, file.path(AQUI, "saida", "comparacao_ols_sar.csv"), row.names = FALSE)
saveRDS(resultados, file.path(AQUI, "saida", "mc3_resultados.rds"))
cat("\nGravado em", file.path(AQUI, "saida"), "\n")
