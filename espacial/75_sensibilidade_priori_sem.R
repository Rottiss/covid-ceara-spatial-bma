# Sensibilidade a priori sob SEM.
#
# O 32 variou as prioris sob MQO e SAR. O SEM so foi estimado com a combinacao
# adotada no TCC (31_sem_hyperg.R). Este script aplica ao SEM o mesmo desenho em
# estrela do 32:
#
#   PRIORI_M = binbeta | binbeta5 | binomial5 | uniforme   (com PRIORI_G = hyperg)
#   PRIORI_G = hyperg | uip                                (com PRIORI_M = binbeta)
#
# As definicoes das prioris sao as do 32, linha por linha.
#
# MARGINAL. A do SEM e a do 31, com a centragem depois da filtragem. Para g fixo
# em n, a forma e a do SAR com g fixo no 32, trocando rho por lambda:
#   log p(y|M,lambda) = -(k/2) log(1+g) - ((n-1)/2) log(g/(1+g) RSS + 1/(1+g) TSS)
# mais o log-determinante e a priori de lambda, integrada por Simpson.
#
# DIFERENCA EM RELACAO AO 31, so de implementacao: cache das marginais, como no
# 74. O cache nao consome numeros aleatorios, entao com binbeta + hyperg e as
# mesmas sementes a cadeia e a do 31 passo a passo. Por isso a combinacao de base
# tambem e rodada aqui: ela deve reproduzir saida/sem_hyperg_<desfecho>_<viz>.csv,
# e a comparacao e feita no fim do script e registrada no log.
#
# VERIFICACOES ANTES DA CADEIA (abortam se falharem):
#   - RSS pre-computado contra ajuste explicito com intercepto, em varios lambda;
#   - em lambda = 0 a marginal do SEM, sem a priori de lambda, iguala a do MQO
#     da mesma priori de g (hyper-g ou g fixo).
#
#   DESFECHO = casos | obitos | letalidade     VIZINHOS = queen | rook
#   MC3_SORTEIOS padrao 1.000.000    MC3_BURN padrao 200.000
#   SAIDA_PRIORI_SEM padrao saida/priori_sem
#
# GRAVA <SAIDA_PRIORI_SEM>/<desfecho>_<viz>_<m>_<g>.rds e aborta se ja existir.

source(if (file.exists("../caminhos.R")) "../caminhos.R" else "caminhos.R")

suppressPackageStartupMessages(library(spdep))

AQUI <- ESPACIAL
RAIZ <- ESTIMACOES
DESFECHO <- Sys.getenv("DESFECHO", "casos")
VIZINHOS <- Sys.getenv("VIZINHOS", "queen")
PRIORI_M <- Sys.getenv("PRIORI_M", "binbeta")
PRIORI_G <- Sys.getenv("PRIORI_G", "hyperg")
N_MC <- as.integer(Sys.getenv("MC3_SORTEIOS", "1000000"))
BURN <- as.integer(Sys.getenv("MC3_BURN", "200000"))
DEST <- Sys.getenv("SAIDA_PRIORI_SEM", file.path(AQUI, "saida", "priori_sem"))
ROTULO <- paste(DESFECHO, VIZINHOS, PRIORI_M, PRIORI_G, sep = "_")
stopifnot(BURN < N_MC, PRIORI_G %in% c("hyperg", "uip"))
arq <- file.path(DEST, paste0(ROTULO, ".rds"))
if (file.exists(arq)) stop("ja existe: ", arq)

hyp <- readRDS(file.path(AQUI, "saida", "hyperg_aux.rds"))
tb  <- readRDS(file.path(AQUI, "saida", "tabela_f21_estendida.rds"))
f21a <- hyp$f21a; a_par <- hyp$a_par; L <- hyp$L; J_MAX <- hyp$J_MAX
TAB <- tb$TAB; passo <- tb$passo; N_TAB <- tb$N_TAB; Z_MAX <- tb$Z_MAX

COLUNA <- c(casos = "z_log_casos100k", obitos = "z_log_obitos100k",
            letalidade = "z_logit_letal")[[DESFECHO]]
dados <- read.csv(file.path(RAIZ, "resultados/estimacoes/matrizes_estimacao/principal",
                            paste0(DESFECHO, ".csv")),
                  stringsAsFactors = FALSE, colClasses = c(codigo_ibge = "character"))
covariaveis <- setdiff(names(dados), c("codigo_ibge", "municipio", COLUNA))
y <- as.numeric(dados[[COLUNA]])
X <- as.matrix(dados[, covariaveis])
n <- length(y); K <- length(covariaveis)
X <- sweep(X, 2, colMeans(X), "-")

listw <- readRDS(file.path(RAIZ, "resultados/diagnostico_espacial",
                           paste0("pesos_", VIZINHOS, ".rds")))
W <- spdep::listw2mat(listw); dimnames(W) <- NULL
stopifnot(nrow(W) == n, max(abs(rowSums(W) - 1)) < 1e-9, max(abs(diag(W))) < 1e-12)

# --------------------------------------------------------- hipergeometrica (31)
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

# --------------------------------------------------------- grade de lambda (31)
mu <- Re(eigen(W, only.values = TRUE)$values)
ng <- 101L
lambda <- seq(-0.999, 0.999, length.out = ng)
ld <- as.numeric(outer(lambda, mu, function(r, m) log1p(-r * m)) %*% rep(1, length(mu)))
a0 <- 1.01
lpr <- (a0 - 1) * log1p(lambda) + (a0 - 1) * log1p(-lambda) -
  (2 * a0 - 1) * log(2) - lbeta(a0, a0)
ld_mais_prior <- ld + lpr
pesos <- c(1, rep(c(4, 2), length.out = ng - 2), 1); h <- lambda[2] - lambda[1]
j0 <- which(abs(lambda) < 1e-12); stopifnot(length(j0) == 1L)

Wy <- as.numeric(W %*% y); WX <- W %*% X
A00 <- crossprod(X); A01 <- crossprod(X, WX); A11 <- crossprod(WX)
b0 <- as.numeric(crossprod(X, y));  b1 <- as.numeric(crossprod(X, Wy))
c0 <- as.numeric(crossprod(WX, y)); c1 <- as.numeric(crossprod(WX, Wy))
s0 <- colSums(X); s1 <- colSums(WX)
yy <- sum(y^2); yWy <- sum(y * Wy); WyWy <- sum(Wy^2); sy <- sum(y); sWy <- sum(Wy)
ys_ys <- yy - 2 * lambda * yWy + lambda^2 * WyWy
soma_ys <- sy - lambda * sWy
tss_s <- ys_ys - n * (soma_ys / n)^2

rss_lambda <- function(idx) {
  k <- length(idx)
  if (k == 0L) return(tss_s)
  A00i <- A00[idx, idx, drop = FALSE]; A01i <- A01[idx, idx, drop = FALSE]
  A11i <- A11[idx, idx, drop = FALSE]
  b0i <- b0[idx]; b1i <- b1[idx]; c0i <- c0[idx]; c1i <- c1[idx]
  s0i <- s0[idx]; s1i <- s1[idx]
  sim <- A01i + t(A01i); soma_bc <- b1i + c0i
  out <- numeric(ng)
  for (j in seq_len(ng)) {
    Lj <- lambda[j]
    Xti <- s0i - Lj * s1i
    XtX <- A00i - Lj * sim + Lj^2 * A11i - tcrossprod(Xti) / n
    Xty <- b0i - Lj * soma_bc + Lj^2 * c1i - Xti * (soma_ys[j] / n)
    beta <- tryCatch(solve(XtX, Xty), error = function(e) NULL)
    if (is.null(beta)) { out[j] <- NA_real_; next }
    out[j] <- tss_s[j] - sum(beta * Xty)
  }
  out
}

g_uip <- n
# log p(y | M, lambda) sem a priori de lambda nem o log-determinante
lf_sem_hyperg <- function(k, rss) {
  lf <- (1 - n) / 2 * log(tss_s)
  if (k > 0L) lf <- lf + log(f21a - 2) - log(k + f21a - 2) +
      log_f21(pmax(0, pmin(1 - rss / tss_s, 1 - 1e-12)), k)
  lf
}
lf_sem_uip <- function(k, rss) {
  sq <- (g_uip / (1 + g_uip)) * rss + (1 / (1 + g_uip)) * tss_s
  if (any(sq <= 0)) return(rep(-Inf, ng))
  -(k / 2) * log(1 + g_uip) - ((n - 1) / 2) * log(sq)
}
lf_sem <- if (PRIORI_G == "uip") lf_sem_uip else lf_sem_hyperg

lm_sem_calc <- function(idx) {
  rss <- rss_lambda(idx)
  if (anyNA(rss) || any(rss <= 0) || any(tss_s <= 0)) return(-Inf)
  lf <- ld_mais_prior + lf_sem(length(idx), rss)
  m <- max(lf); if (!is.finite(m)) return(-Inf)
  m + log(h / 3 * sum(pesos * exp(lf - m)))
}
cache <- new.env(hash = TRUE)
lm_sem <- function(idx) {
  chave <- if (length(idx)) paste(idx, collapse = ",") else "vazio"
  v <- cache[[chave]]
  if (!is.null(v)) return(v)
  v <- lm_sem_calc(idx)
  assign(chave, v, envir = cache)
  v
}

# --------------------------------------------------------- verificacoes
yc <- y - mean(y); ta <- sum(yc^2)
lm_ols <- function(idx) {
  k <- length(idx)
  a <- if (k == 0L) ta else {
    e <- lm.fit(cbind(1, X[, idx, drop = FALSE]), y)$residuals; sum(e^2) }
  if (PRIORI_G == "uip")
    return(-(k / 2) * log(1 + g_uip) - ((n - 1) / 2) * log(ta - (g_uip / (1 + g_uip)) * (ta - a)))
  nl <- (1 - n) / 2 * log(ta)
  if (k == 0L) return(nl)
  nl + log(f21a - 2) - log(k + f21a - 2) + log_f21(1 - a / ta, k)
}
cat("=== Verificacoes ===\n")
set.seed(20200731L)
am <- c(list(integer(0)), replicate(12, sort(sample(K, sample(1:12, 1))), simplify = FALSE))
js <- unique(c(1L, round(ng * c(0.1, 0.3)), j0, round(ng * c(0.7, 0.9)), ng))
pior_rss <- 0; pior_escala <- 0
for (idx in am) {
  r <- rss_lambda(idx)
  for (j in js) {
    ys <- y - lambda[j] * Wy
    Xs <- X - lambda[j] * WX
    e <- if (length(idx)) lm.fit(cbind(1, Xs[, idx, drop = FALSE]), ys)$residuals else ys - mean(ys)
    pior_rss <- max(pior_rss, abs(r[j] - sum(e^2)))
  }
  pior_escala <- max(pior_escala, abs(lf_sem(length(idx), r)[j0] - lm_ols(idx)))
}
cat("  maior diferenca de RSS:", format(pior_rss, scientific = TRUE), "\n")
cat("  SEM em lambda = 0 contra MQO,", PRIORI_G, ":", format(pior_escala, scientific = TRUE), "\n")
if (pior_rss > 1e-6 || pior_escala > 1e-8) stop("verificacao reprovada")
cat("  veredito: APROVADO\n\n")

# --------------------------------------------------------- priori de modelo (32)
lpm_binbeta <- function(m) {
  b <- (K - m) / m
  vapply(0:K, function(k) lgamma(1 + k) + lgamma(b + K - k), numeric(1))
}
lpm <- switch(PRIORI_M,
  binbeta   = lpm_binbeta(K / 2),
  binbeta5  = lpm_binbeta(5),
  uniforme  = numeric(K + 1),
  binomial5 = { th <- 5 / K
                vapply(0:K, function(k) k * log(th) + (K - k) * log(1 - th),
                       numeric(1)) },
  stop("PRIORI_M desconhecida"))
stopifnot(all(is.finite(lpm)))

# --------------------------------------------------------- amostrador (31)
tipos_disp <- function(k) { t <- character(0)
  if (k < K) t <- c(t, "birth"); if (k > 0L) t <- c(t, "death")
  if (k > 0L && k < K) t <- c(t, "move"); t }
log_q <- function(k, tipo) -log(length(tipos_disp(k))) -
  log(switch(tipo, birth = K - k, death = k, move = k * (K - k)))

mc3 <- function(semente) {
  set.seed(semente)
  atual <- sort(sample(seq_len(K), sample(0:K, 1)))
  post <- lm_sem(atual) + lpm[length(atual) + 1]
  inc <- numeric(K); g <- 0L; sk <- 0; ace <- 0L
  for (it in seq_len(N_MC)) {
    k <- length(atual); disp <- tipos_disp(k)
    tp <- disp[sample.int(length(disp), 1L)]
    fora <- setdiff(seq_len(K), atual)
    prop <- switch(tp,
                   birth = sort(c(atual, fora[sample.int(length(fora), 1L)])),
                   death = atual[-sample.int(k, 1L)],
                   move  = sort(c(atual[-sample.int(k, 1L)],
                                  fora[sample.int(length(fora), 1L)])))
    pp <- lm_sem(prop) + lpm[length(prop) + 1]
    inv <- switch(tp, birth = "death", death = "birth", move = "move")
    la <- pp - post + log_q(length(prop), inv) - log_q(k, tp)
    if (is.finite(la) && log(runif(1)) < la) {
      atual <- prop; post <- pp; ace <- ace + 1L
    }
    if (it > BURN) { g <- g + 1L
      if (length(atual) > 0L) inc[atual] <- inc[atual] + 1
      sk <- sk + length(atual) }
    if (it %% 50000L == 0L) {
      cat(sprintf("    cadeia %d: %d de %d | k = %d | cache %d | %.0f s\n", semente, it,
                  N_MC, length(atual), length(ls(cache)), proc.time()[3] - T0))
      flush.console()
    }
  }
  list(pip = inc / g, tam = sk / g, taxa = ace / N_MC)
}

cat("=== SEM:", ROTULO, "===\n")
cat("  sorteios:", N_MC, "| descarte:", BURN, "\n\n")
T0 <- proc.time()[3]
a <- mc3(20200731L); b <- mc3(987654L)
pip <- (a$pip + b$pip) / 2
names(pip) <- names(a$pip) <- names(b$pip) <- covariaveis

Dst <- function(p) { m <- mean(p); Kk <- length(p)
  sd(p) * sqrt((Kk - 1) / Kk) / sqrt(m * (1 - m)) }
cat("\n  tempo total:", round(proc.time()[3] - T0, 1), "s\n")
cat("  aceitacao:", round(a$taxa, 3), "e", round(b$taxa, 3), "\n")
cat("  tamanho medio:", round(a$tam, 3), "e", round(b$tam, 3), "\n")
cat("  maior diferenca de PIP entre cadeias:", signif(max(abs(a$pip - b$pip)), 3), "\n")
cat("  avaliacoes fora da tabela de interpolacao:", fora_tab, "\n")
cat("  D:", round(Dst(pip), 4), "| acima de 0,50:", sum(pip >= 0.5),
    if (any(pip >= 0.5)) paste0("(", paste(covariaveis[pip >= 0.5], collapse = ", "), ")") else "", "\n")

# a combinacao de base deve reproduzir o 31
reproducao <- NA_real_
if (PRIORI_M == "binbeta" && PRIORI_G == "hyperg" && N_MC == 1000000L && BURN == 200000L) {
  ref <- read.csv(file.path(AQUI, "saida", sprintf("sem_hyperg_%s_%s.csv", DESFECHO, VIZINHOS)),
                  stringsAsFactors = FALSE)
  stopifnot(setequal(ref$covariavel, covariaveis))
  reproducao <- max(abs(pip[ref$covariavel] - ref$sem))
  cat("  reproducao do 31, maior diferenca de PIP:", format(reproducao, scientific = TRUE), "\n")
}

dir.create(DEST, showWarnings = FALSE, recursive = TRUE)
stopifnot(!file.exists(arq))
saveRDS(list(desfecho = DESFECHO, vizinhos = VIZINHOS, priori_m = PRIORI_M,
             priori_g = PRIORI_G, n_mc = N_MC, burn = BURN, covariaveis = covariaveis,
             pip = pip, a = a, b = b, tam = (a$tam + b$tam) / 2, D = Dst(pip),
             conv = max(abs(a$pip - b$pip)), fora_tab = fora_tab,
             reproducao_31 = reproducao), arq)
cat("Gravado em", arq, "\n")
