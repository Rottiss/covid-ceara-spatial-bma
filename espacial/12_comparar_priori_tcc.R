# Etapa 4 refeita sob a priori do TCC: hyper-g UIP com priori de modelo
# Binomial-Beta (mprior.size = K/2).
#
# A versao anterior (06_comparar_ols_sar.R) usou g fixo e priori uniforme, o que
# isolava o efeito espacial mas nao permitia dizer se o achado sobrevive na
# especificacao publicada. Este script responde isso.
#
# Componentes validados antes de rodar:
#   10_validar_hyperg.R  - a hyper-g + Binomial-Beta reproduz o BMS a 2e-15
#   11_hyperg_sar.R      - a hipergeometrica em log bate com a serie direta a
#                          7e-15, e o SAR em rho = 0 devolve o MQO
#
# Desfecho: casos. W: Queen. Igual a Etapa 4 original, para ser comparavel.

source(if (file.exists("../caminhos.R")) "../caminhos.R" else "caminhos.R")

AQUI <- ESPACIAL
base <- readRDS(file.path(AQUI, "saida", "base_sar.rds"))
aux <- readRDS(file.path(AQUI, "saida", "grade_sar.rds"))
hyp <- readRDS(file.path(AQUI, "saida", "hyperg_aux.rds"))

y <- base$y; X_total <- base$X; n <- base$n; covariaveis <- base$covariaveis
grade <- aux$grade; ta <- aux$ta; tb <- aux$tb; tc <- aux$tc
yc <- aux$yc; Wyc <- aux$Wyc
K <- length(covariaveis)
f21a <- hyp$f21a; a_par <- hyp$a_par; L <- hyp$L; J_MAX <- hyp$J_MAX

log_f21_vec <- function(z, k) {
  cc <- (f21a + k) / 2
  j_pico <- max(0, ceiling(max((z * a_par - cc) / pmax(1 - z, 1e-12))))
  J <- min(J_MAX, max(300L, as.integer(j_pico * 4)))
  lt <- cbind(0, outer(log(z), seq_len(J)) +
                rep(L[[k + 1]][seq_len(J)], each = length(z)))
  m <- apply(lt, 1, max)
  m + log(rowSums(exp(lt - m)))
}

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

lm_ols <- function(est) {
  k <- unname(est["k"])
  nl <- (1 - n) / 2 * log(ta)
  if (k == 0L) return(nl)
  unname(nl + log(f21a - 2) - log(k + f21a - 2) +
           log_f21_vec(1 - est["a"] / ta, k))
}

lm_sar <- function(est) {
  k <- unname(est["k"])
  rss <- est["a"] - 2 * grade$rho * est["b"] + grade$rho2 * est["c"]
  tss <- ta - 2 * grade$rho * tb + grade$rho2 * tc
  if (any(rss <= 0) || any(tss <= 0)) return(-Inf)
  R2 <- pmax(0, 1 - rss / tss)
  lf <- grade$ld_mais_prior + (1 - n) / 2 * log(tss)
  if (k > 0L) lf <- lf + log(f21a - 2) - log(k + f21a - 2) + log_f21_vec(R2, k)
  m <- max(lf)
  if (!is.finite(m)) return(-Inf)
  unname(m + log(grade$h / 3 * sum(grade$pesos * exp(lf - m))))
}

# Priori de modelo Binomial-Beta, a menos de constante.
lpm <- vapply(0:K, function(k) lgamma(1 + k) + lgamma(1 + K - k), numeric(1))

tipos_disponiveis <- function(k) {
  t <- character(0)
  if (k < K) t <- c(t, "birth")
  if (k > 0L) t <- c(t, "death")
  if (k > 0L && k < K) t <- c(t, "move")
  t
}
log_q <- function(k, tipo) {
  -log(length(tipos_disponiveis(k))) -
    log(switch(tipo, birth = K - k, death = k, move = k * (K - k)))
}

mc3 <- function(log_marginal, n_sorteios, burn, semente) {
  set.seed(semente)
  atual <- sort(sample(seq_len(K), sample(0:K, 1)))
  post_atual <- log_marginal(estatisticas(atual)) + lpm[length(atual) + 1]
  aceitos <- 0L; inclusao <- numeric(K); guardados <- 0L; soma_k <- 0
  unicos <- new.env(hash = TRUE, parent = emptyenv())
  for (it in seq_len(n_sorteios)) {
    k <- length(atual)
    tipos <- tipos_disponiveis(k)
    tipo <- tipos[sample.int(length(tipos), 1L)]
    fora <- setdiff(seq_len(K), atual)
    proposto <- switch(
      tipo,
      birth = sort(c(atual, fora[sample.int(length(fora), 1L)])),
      death = atual[-sample.int(k, 1L)],
      move  = sort(c(atual[-sample.int(k, 1L)], fora[sample.int(length(fora), 1L)]))
    )
    post_prop <- log_marginal(estatisticas(proposto)) + lpm[length(proposto) + 1]
    inverso <- switch(tipo, birth = "death", death = "birth", move = "move")
    log_alfa <- post_prop - post_atual + log_q(length(proposto), inverso) - log_q(k, tipo)
    if (is.finite(log_alfa) && log(runif(1)) < log_alfa) {
      atual <- proposto; post_atual <- post_prop; aceitos <- aceitos + 1L
    }
    if (it > burn) {
      guardados <- guardados + 1L
      if (length(atual) > 0L) inclusao[atual] <- inclusao[atual] + 1
      soma_k <- soma_k + length(atual)
      assign(if (length(atual) == 0L) "vazio" else paste(atual, collapse = ","),
             TRUE, envir = unicos)
    }
  }
  list(pip = inclusao / guardados, taxa = aceitos / n_sorteios,
       tamanho_medio = soma_k / guardados, unicos = length(ls(unicos)))
}

N <- as.integer(Sys.getenv("MC3_SORTEIOS", "100000"))
BURN <- as.integer(N * 0.2)

cat("=== Configuracao ===\n")
cat("  priori: hyper-g UIP (a =", format(f21a, digits = 10),
    ") + Binomial-Beta com mprior.size = K/2\n")
cat("  desfecho: casos | W: Queen | K =", K, "\n")
cat("  sorteios:", N, "| descarte:", BURN, "\n\n")

res <- list()
for (rotulo in c("MQO", "SAR")) {
  f <- if (rotulo == "MQO") lm_ols else lm_sar
  cat("Rodando", rotulo, "... "); flush.console()
  t0 <- proc.time()[3]
  a <- mc3(f, N, BURN, 20200731L)
  b <- mc3(f, N, BURN, 987654L)
  cat(round(proc.time()[3] - t0, 1), "s\n")
  cat("  aceitacao:", round(a$taxa, 3), "e", round(b$taxa, 3),
      "| tamanho medio:", round(a$tamanho_medio, 2), "e", round(b$tamanho_medio, 2),
      "| unicos:", a$unicos, "e", b$unicos, "\n")
  cat("  maior diferenca de PIP entre cadeias:", signif(max(abs(a$pip - b$pip)), 3), "\n")
  res[[rotulo]] <- list(pip = (a$pip + b$pip) / 2, conv = max(abs(a$pip - b$pip)))
}

tab <- data.frame(covariavel = covariaveis,
                  pip_mqo_tcc = res$MQO$pip,
                  pip_sar_tcc = res$SAR$pip,
                  stringsAsFactors = FALSE)
tab$diferenca <- tab$pip_sar_tcc - tab$pip_mqo_tcc

# Junta com a versao de priori uniforme, para mostrar o efeito da priori.
arq <- file.path(AQUI, "saida", "comparacao_ols_sar.csv")
if (file.exists(arq)) {
  ant <- read.csv(arq, stringsAsFactors = FALSE)
  names(ant)[names(ant) == "pip_mqo"] <- "pip_mqo_unif"
  names(ant)[names(ant) == "pip_sar"] <- "pip_sar_unif"
  ant$diferenca <- NULL
  tab <- merge(tab, ant, by = "covariavel")
}
tab <- tab[order(-tab$pip_mqo_tcc), ]

cat("\n=== PIP sob a priori do TCC (hyper-g + Binomial-Beta) ===\n")
print(transform(tab[, intersect(c("covariavel", "pip_mqo_tcc", "pip_sar_tcc",
                                  "diferenca", "pip_mqo_unif", "pip_sar_unif"),
                                names(tab))],
                pip_mqo_tcc = round(pip_mqo_tcc, 4),
                pip_sar_tcc = round(pip_sar_tcc, 4),
                diferenca = round(diferenca, 4)), row.names = FALSE)

cat("\n=== Quem passa de 0,50 ===\n")
cat("  MQO:", paste(tab$covariavel[tab$pip_mqo_tcc >= 0.5], collapse = ", "), "\n")
cat("  SAR:", paste(tab$covariavel[tab$pip_sar_tcc >= 0.5], collapse = ", "), "\n")
cat("\n  maior queda:", round(min(tab$diferenca), 4), "em",
    tab$covariavel[which.min(tab$diferenca)], "\n")
cat("  maior alta :", round(max(tab$diferenca), 4), "em",
    tab$covariavel[which.max(tab$diferenca)], "\n")
cat("  mudam mais de 0,10:", sum(abs(tab$diferenca) > 0.10), "de", K, "\n")
cat("  ruido de Monte Carlo: MQO", signif(res$MQO$conv, 3),
    "| SAR", signif(res$SAR$conv, 3), "\n")

write.csv(tab, file.path(AQUI, "saida", "comparacao_priori_tcc.csv"), row.names = FALSE)
cat("\nGravado em", file.path(AQUI, "saida"), "\n")
