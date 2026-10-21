# Estimate posterior specification probabilities jointly with covariate inclusion.
#
# Aqui a especificacao entra no espaco amostral. A cadeia percorre pares
# (modelo de covariaveis, especificacao), com um quarto tipo de passo que troca
# MQO por SAR mantendo o conjunto de covariaveis. A frequencia posterior da
# especificacao espacial e a razao de massas
#
#   P(SAR | y) / P(MQO | y) = soma_M p(y|M,SAR) P(M) / soma_M p(y|M,MQO) P(M),
#
# com a mesma priori de modelo dos dois lados e priori 1/2 para a especificacao.
# As duas log-marginais ja estao na mesma escala: mesmo y, mesma priori sobre
# beta e sigma. No SAR ha ainda a integracao em rho com priori Beta(1,01; 1,01)
# em (-1, 1), entao a razao ja penaliza o parametro a mais.
#
# A cadeia tambem acumula a posteriori de rho ponderada por modelo, que e a
# posteriori de rho sob media bayesiana de modelos.
#
#   DESFECHO = casos | obitos | letalidade      VIZINHOS = queen | rook

source(if (file.exists("../caminhos.R")) "../caminhos.R" else "caminhos.R")

suppressPackageStartupMessages(library(spdep))

AQUI <- ESPACIAL
RAIZ <- ESTIMACOES
hyp <- readRDS(file.path(AQUI, "saida", "hyperg_aux.rds"))
tb  <- readRDS(file.path(AQUI, "saida", "tabela_f21_estendida.rds"))

DESFECHO <- Sys.getenv("DESFECHO", "casos")
VIZINHOS <- Sys.getenv("VIZINHOS", "queen")
N_MC <- as.integer(Sys.getenv("MC3_SORTEIOS", "1000000"))
BURN <- as.integer(Sys.getenv("MC3_BURN", "200000"))

COLUNA <- c(casos = "z_log_casos100k", obitos = "z_log_obitos100k",
            letalidade = "z_logit_letal")[[DESFECHO]]

dados <- read.csv(file.path(RAIZ, "resultados/estimacoes/matrizes_estimacao/principal",
                            paste0(DESFECHO, ".csv")),
                  stringsAsFactors = FALSE, colClasses = c(codigo_ibge = "character"))
covariaveis <- setdiff(names(dados), c("codigo_ibge", "municipio", COLUNA))
y <- as.numeric(dados[[COLUNA]])
X <- as.matrix(dados[, covariaveis])
n <- length(y); K <- length(covariaveis)
Xc <- sweep(X, 2, colMeans(X), "-")

listw <- readRDS(file.path(RAIZ, "resultados/diagnostico_espacial",
                           paste0("pesos_", VIZINHOS, ".rds")))
W <- spdep::listw2mat(listw); dimnames(W) <- NULL
stopifnot(nrow(W) == n, max(abs(rowSums(W) - 1)) < 1e-9)

f21a <- hyp$f21a; a_par <- hyp$a_par; L <- hyp$L; J_MAX <- hyp$J_MAX
TAB <- tb$TAB; passo <- tb$passo; N_TAB <- tb$N_TAB; Z_MAX <- tb$Z_MAX

fora_tab <- 0L
log_f21_serie1 <- function(zi, k) {
  cc <- (f21a + k) / 2
  j_pico <- max(0, ceiling((zi * a_par - cc) / max(1 - zi, 1e-12)))
  J <- min(J_MAX, max(300L, as.integer(j_pico * 4)))
  lt <- c(0, log(zi) * seq_len(J) + L[[k + 1]][seq_len(J)])
  m <- max(lt); m + log(sum(exp(lt - m)))
}
log_f21 <- function(z, k) {
  acima <- z > Z_MAX
  if (any(acima)) {
    fora_tab <<- fora_tab + sum(acima)
    out <- numeric(length(z))
    if (any(!acima)) out[!acima] <- log_f21(z[!acima], k)
    out[acima] <- vapply(z[acima], log_f21_serie1, numeric(1), k = k)
    return(out)
  }
  pos <- z / passo; i <- as.integer(pos) + 1L
  i[i >= N_TAB] <- N_TAB - 1L
  frac <- pos - (i - 1L); off <- k * N_TAB
  TAB[i + off] * (1 - frac) + TAB[i + 1L + off] * frac
}

mu <- Re(eigen(W, only.values = TRUE)$values)
ng <- 201L
rho <- seq(-0.999, 0.999, length.out = ng)
ld <- as.numeric(outer(rho, mu, function(r, m) log1p(-r * m)) %*% rep(1, length(mu)))
a0 <- 1.01
lpr <- (a0 - 1) * log1p(rho) + (a0 - 1) * log1p(-rho) -
  (2 * a0 - 1) * log(2) - lbeta(a0, a0)
ld_mais_prior <- ld + lpr
pesos <- c(1, rep(c(4, 2), length.out = ng - 2), 1); h <- rho[2] - rho[1]

yc <- y - mean(y); Wy <- as.numeric(W %*% y); Wyc <- Wy - mean(Wy)
ta <- sum(yc^2); tbb <- sum(yc * Wyc); tcc <- sum(Wyc^2)

estat <- function(idx) {
  if (length(idx) == 0L) { eo <- yc; ed <- Wyc } else {
    Q <- qr.Q(qr(Xc[, idx, drop = FALSE]))
    eo <- yc - Q %*% crossprod(Q, yc); ed <- Wyc - Q %*% crossprod(Q, Wyc)
  }
  eo <- as.numeric(eo); ed <- as.numeric(ed)
  c(k = length(idx), a = sum(eo^2), b = sum(eo * ed), c = sum(ed^2))
}
lm_ols <- function(e) {
  k <- unname(e["k"]); nl <- (1 - n) / 2 * log(ta)
  if (k == 0L) return(nl)
  unname(nl + log(f21a - 2) - log(k + f21a - 2) + log_f21(1 - e["a"] / ta, k))
}
# Devolve a log-marginal do SAR e, opcionalmente, o vetor de pesos em rho.
lm_sar_pesos <- function(e, quer_rho = FALSE) {
  k <- unname(e["k"])
  rss <- e["a"] - 2 * rho * e["b"] + rho^2 * e["c"]
  tss <- ta - 2 * rho * tbb + rho^2 * tcc
  if (any(rss <= 0) || any(tss <= 0)) return(list(lm = -Inf, w = NULL))
  lf <- ld_mais_prior + (1 - n) / 2 * log(tss)
  if (k > 0L) lf <- lf + log(f21a - 2) - log(k + f21a - 2) +
      log_f21(pmax(0, 1 - rss / tss), k)
  m <- max(lf); if (!is.finite(m)) return(list(lm = -Inf, w = NULL))
  ex <- exp(lf - m)
  val <- m + log(h / 3 * sum(pesos * ex))
  list(lm = unname(val), w = if (quer_rho) ex * pesos / sum(ex * pesos) else NULL)
}
lm_sar <- function(e) lm_sar_pesos(e)$lm

lpm <- vapply(0:K, function(k) lgamma(1 + k) + lgamma(1 + K - k), numeric(1))

# Tipos de passo. "troca" muda so a especificacao e esta sempre disponivel, de
# modo que o numero de tipos continua dependendo apenas de k e a razao de
# Hastings mantem a forma usada no resto do projeto.
tipos_disp <- function(k) {
  t <- "troca"
  if (k < K) t <- c(t, "birth"); if (k > 0L) t <- c(t, "death")
  if (k > 0L && k < K) t <- c(t, "move"); t
}
log_q <- function(k, tipo) {
  -log(length(tipos_disp(k))) -
    log(switch(tipo, troca = 1, birth = K - k, death = k, move = k * (K - k)))
}

marg <- function(idx, esp) if (esp == "SAR") lm_sar(estat(idx)) else lm_ols(estat(idx))

mc3 <- function(semente) {
  set.seed(semente)
  atual <- sort(sample(seq_len(K), sample(0:K, 1)))
  esp <- "MQO"
  post <- marg(atual, esp) + lpm[length(atual) + 1]
  inc_sar <- numeric(K); inc_ols <- numeric(K)
  n_sar <- 0L; n_ols <- 0L; g <- 0L; ace <- 0L
  w_rho <- numeric(ng)
  for (it in seq_len(N_MC)) {
    k <- length(atual); disp <- tipos_disp(k)
    tp <- disp[sample.int(length(disp), 1L)]
    if (tp == "troca") {
      prop <- atual; esp_p <- if (esp == "SAR") "MQO" else "SAR"
    } else {
      esp_p <- esp
      fora <- setdiff(seq_len(K), atual)
      prop <- switch(tp,
                     birth = sort(c(atual, fora[sample.int(length(fora), 1L)])),
                     death = atual[-sample.int(k, 1L)],
                     move  = sort(c(atual[-sample.int(k, 1L)],
                                    fora[sample.int(length(fora), 1L)])))
    }
    pp <- marg(prop, esp_p) + lpm[length(prop) + 1]
    inv <- switch(tp, troca = "troca", birth = "death", death = "birth", move = "move")
    la <- pp - post + log_q(length(prop), inv) - log_q(k, tp)
    if (is.finite(la) && log(runif(1)) < la) {
      atual <- prop; esp <- esp_p; post <- pp; ace <- ace + 1L
    }
    if (it > BURN) {
      g <- g + 1L
      if (esp == "SAR") {
        n_sar <- n_sar + 1L
        if (length(atual) > 0L) inc_sar[atual] <- inc_sar[atual] + 1
        w <- lm_sar_pesos(estat(atual), quer_rho = TRUE)$w
        if (!is.null(w)) w_rho <- w_rho + w
      } else {
        n_ols <- n_ols + 1L
        if (length(atual) > 0L) inc_ols[atual] <- inc_ols[atual] + 1
      }
    }
  }
  list(p_sar = n_sar / g, n_sar = n_sar, n_ols = n_ols, taxa = ace / N_MC,
       pip_sar = if (n_sar > 0) inc_sar / n_sar else rep(NA_real_, K),
       pip_ols = if (n_ols > 0) inc_ols / n_ols else rep(NA_real_, K),
       w_rho = if (n_sar > 0) w_rho / n_sar else NULL)
}

cat("=== Especificacao dentro do espaco amostral ===\n")
cat("  desfecho:", DESFECHO, "| vizinhanca:", VIZINHOS, "\n")
cat("  sorteios:", N_MC, "| descarte:", BURN, "\n\n")
t0 <- proc.time()[3]
a <- mc3(20200731L); b <- mc3(987654L)
cat("  duas cadeias em", round(proc.time()[3] - t0, 1), "s\n")
cat("  aceitacao:", round(a$taxa, 3), "e", round(b$taxa, 3), "\n\n")

cat("  visitas ao SAR: cadeia 1", a$n_sar, "de", a$n_sar + a$n_ols,
    "| cadeia 2", b$n_sar, "de", b$n_sar + b$n_ols, "\n")
cat("  P(SAR | y):", signif(a$p_sar, 6), "e", signif(b$p_sar, 6), "\n")

# Se a cadeia nunca visita o MQO, a frequencia observada e 1, mas nao fornece um
# limite deterministico para a probabilidade posterior. A comparacao abaixo usa
# fatores de Bayes condicionais a conjuntos de covariaveis sorteados da priori.
# A mediana desses fatores nao e o fator de Bayes agregado entre especificacoes.
if (a$n_ols == 0L || b$n_ols == 0L) {
  cat("\n  A cadeia nao visitou o MQO depois do descarte. Isso indica separacao\n")
  cat("  forte, mas nao autoriza o limite P(MQO | y) < 1/N.\n")
  cat("  Diferenca de log-marginal em conjuntos sorteados da priori de modelos:\n")
  set.seed(1L)
  amostra <- replicate(200, sort(sample(K, sample(0:K, 1))), simplify = FALSE)
  dif <- vapply(amostra, function(idx) lm_sar(estat(idx)) - lm_ols(estat(idx)),
                numeric(1))
  cat(sprintf("    mediana %.1f | minimo %.1f | maximo %.1f (em log)\n",
              median(dif), min(dif), max(dif)))
  cat("    fator de Bayes condicional mediano:",
      format(exp(median(dif)), digits = 3), "\n")
}

if (!is.null(a$w_rho)) {
  wr <- (a$w_rho + b$w_rho) / 2; wr <- wr / sum(wr)
  m_rho <- sum(wr * rho)
  ac <- cumsum(wr)
  q <- rho[c(which.max(ac >= 0.025), which.max(ac >= 0.5), which.max(ac >= 0.975))]
  cat("\n  Posteriori de rho sob media de modelos:\n")
  cat(sprintf("    media %.4f | mediana %.4f | intervalo de 95%%: [%.4f, %.4f]\n",
              m_rho, q[2], q[1], q[3]))
  cat(sprintf("    massa em rho > 0: %.4f\n", sum(wr[rho > 0])))
}

saveRDS(list(desfecho = DESFECHO, vizinhos = VIZINHOS, a = a, b = b,
             rho = rho, covariaveis = covariaveis),
        file.path(AQUI, "saida", sprintf("especificacao_%s_%s.rds", DESFECHO, VIZINHOS)))
cat("\nGravado em saida/\n")
