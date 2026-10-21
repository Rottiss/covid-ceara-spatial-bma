# Compare D ratios in municipal subsamples of different sizes.
# This does not establish how results would change in larger samples.
#
# LIMITACAO ASSUMIDA. Subamostrar nao e o mesmo que acrescentar unidades:
# remover municipios rompe contiguidades e reduz a densidade da vizinhanca, o
# que por si so pode alterar o resultado. O exercicio mostra tendencia, nao
# magnitude, e a extrapolacao para n maior nao esta demonstrada.
#
# A matriz de vizinhanca e reconstruida em cada subamostra a partir da
# contiguidade binaria original restrita ao subconjunto, renormalizada por linha.
# Subamostras que produzam ilhas sao descartadas e resorteadas.

source(if (file.exists("../caminhos.R")) "../caminhos.R" else "caminhos.R")

suppressPackageStartupMessages(library(spdep))

AQUI <- ESPACIAL
base <- readRDS(file.path(AQUI, "saida", "base_sar.rds"))
hyp  <- readRDS(file.path(AQUI, "saida", "hyperg_aux.rds"))

X <- base$X; W_full <- base$W; y_full <- base$y
covariaveis <- base$covariaveis; K <- length(covariaveis)
f21a <- hyp$f21a; a_par <- hyp$a_par; L <- hyp$L; J_MAX <- hyp$J_MAX

# Contiguidade binaria original, antes da normalizacao por linha.
B_full <- (W_full > 0) * 1

D <- function(p) { m <- mean(p); Kk <- length(p)
  sd(p) * sqrt((Kk - 1) / Kk) / sqrt(m * (1 - m)) }

log_f21_serie <- function(zi, k, n) {
  cc <- (f21a + k) / 2; ap <- (n - 1) / 2
  j_pico <- max(0, ceiling((zi * ap - cc) / max(1 - zi, 1e-12)))
  J <- min(J_MAX, max(300L, as.integer(j_pico * 4)))
  lt <- c(0, cumsum(log(ap + 0:(J - 1)) - log(cc + 0:(J - 1))) + log(zi) * seq_len(J))
  m <- max(lt); m + log(sum(exp(lt - m)))
}

# O parametro a = (n-1)/2 da hipergeometrica muda com n, entao a tabela de
# interpolacao precisa ser reconstruida por tamanho de amostra. Sem tabela, cada
# avaliacao percorre a serie 201 vezes e o experimento fica inviavel.
Z_MAX <- 0.90; N_TAB <- 60001L
z_tab <- seq(0, Z_MAX, length.out = N_TAB); passo_tab <- z_tab[2] - z_tab[1]
TABS <- new.env(parent = emptyenv())
tabela_para <- function(n) {
  chave <- as.character(n)
  if (!is.null(TABS[[chave]])) return(TABS[[chave]])
  M <- vapply(0:K, function(k)
    vapply(z_tab, log_f21_serie, numeric(1), k = k, n = n), numeric(N_TAB))
  assign(chave, M, envir = TABS)
  M
}
log_f21_tab <- function(z, k, TAB) {
  pos <- z / passo_tab; i <- as.integer(pos) + 1L
  i[i >= N_TAB] <- N_TAB - 1L
  frac <- pos - (i - 1L); off <- k * N_TAB
  TAB[i + off] * (1 - frac) + TAB[i + 1L + off] * frac
}

rodar <- function(idx_mun, n_sorteios, semente) {
  n <- length(idx_mun)
  TAB <- tabela_para(n)
  Bs <- B_full[idx_mun, idx_mun, drop = FALSE]
  if (any(rowSums(Bs) == 0)) return(NULL)
  Ws <- Bs / rowSums(Bs)
  Xs <- X[idx_mun, , drop = FALSE]
  Xc <- sweep(Xs, 2, colMeans(Xs), "-")
  y <- y_full[idx_mun]
  yc <- y - mean(y); Wy <- as.numeric(Ws %*% y); Wyc <- Wy - mean(Wy)
  ta <- sum(yc^2); tbb <- sum(yc * Wyc); tcc <- sum(Wyc^2)

  mu <- Re(eigen(Ws, only.values = TRUE)$values)
  ng <- 201L
  rho <- seq(-0.999, 0.999, length.out = ng)
  ld <- as.numeric(outer(rho, mu, function(r, m) log1p(-r * m)) %*% rep(1, length(mu)))
  a0 <- 1.01
  lpr <- (a0 - 1) * log1p(rho) + (a0 - 1) * log1p(-rho) -
    (2 * a0 - 1) * log(2) - lbeta(a0, a0)
  pesos <- c(1, rep(c(4, 2), length.out = ng - 2), 1); h <- rho[2] - rho[1]

  est <- function(ix) {
    if (length(ix) == 0L) { eo <- yc; ed <- Wyc } else {
      Q <- qr.Q(qr(Xc[, ix, drop = FALSE]))
      eo <- yc - Q %*% crossprod(Q, yc); ed <- Wyc - Q %*% crossprod(Q, Wyc)
    }
    eo <- as.numeric(eo); ed <- as.numeric(ed)
    c(k = length(ix), a = sum(eo^2), b = sum(eo * ed), c = sum(ed^2))
  }
  ols <- function(e) {
    k <- unname(e["k"]); nl <- (1 - n) / 2 * log(ta)
    if (k == 0L) return(nl)
    nl + log(f21a - 2) - log(k + f21a - 2) + log_f21_tab(1 - e["a"] / ta, k, TAB)
  }
  sar <- function(e) {
    k <- unname(e["k"])
    rss <- e["a"] - 2 * rho * e["b"] + rho^2 * e["c"]
    tss <- ta - 2 * rho * tbb + rho^2 * tcc
    if (any(rss <= 0) || any(tss <= 0)) return(-Inf)
    R2 <- pmax(0, 1 - rss / tss)
    lf <- ld + lpr + (1 - n) / 2 * log(tss)
    if (k > 0L) lf <- lf + log(f21a - 2) - log(k + f21a - 2) +
        log_f21_tab(pmin(R2, Z_MAX), k, TAB)
    m <- max(lf); if (!is.finite(m)) return(-Inf)
    m + log(h / 3 * sum(pesos * exp(lf - m)))
  }

  lpm <- vapply(0:K, function(k) lgamma(1 + k) + lgamma(1 + K - k), numeric(1))
  tipos <- function(k) { t <- character(0)
    if (k < K) t <- c(t, "birth"); if (k > 0L) t <- c(t, "death")
    if (k > 0L && k < K) t <- c(t, "move"); t }
  lq <- function(k, tp) -log(length(tipos(k))) -
    log(switch(tp, birth = K - k, death = k, move = k * (K - k)))

  mc3 <- function(marg, sem) {
    set.seed(sem)
    atual <- sort(sample(K, sample(0:K, 1)))
    post <- marg(est(atual)) + lpm[length(atual) + 1]
    inc <- numeric(K); g <- 0L; sk <- 0
    burn <- as.integer(n_sorteios * 0.2)
    for (it in seq_len(n_sorteios)) {
      k <- length(atual); tp <- tipos(k)[sample.int(length(tipos(k)), 1L)]
      fora <- setdiff(seq_len(K), atual)
      prop <- switch(tp,
                     birth = sort(c(atual, fora[sample.int(length(fora), 1L)])),
                     death = atual[-sample.int(k, 1L)],
                     move  = sort(c(atual[-sample.int(k, 1L)],
                                    fora[sample.int(length(fora), 1L)])))
      pp <- marg(est(prop)) + lpm[length(prop) + 1]
      inv <- switch(tp, birth = "death", death = "birth", move = "move")
      la <- pp - post + lq(length(prop), inv) - lq(k, tp)
      if (is.finite(la) && log(runif(1)) < la) { atual <- prop; post <- pp }
      if (it > burn) { g <- g + 1L
        if (length(atual) > 0L) inc[atual] <- inc[atual] + 1
        sk <- sk + length(atual) }
    }
    list(pip = inc / g, tam = sk / g)
  }

  a <- mc3(ols, 20200731L); b <- mc3(sar, 20200731L)
  list(D_ols = D(a$pip), D_sar = D(b$pip),
       k_ols = a$tam, k_sar = b$tam, n = n,
       ligacoes = sum(Bs) / n)
}

N_MC <- as.integer(Sys.getenv("MC3_SORTEIOS", "150000"))
REPS <- as.integer(Sys.getenv("REPLICACOES", "6"))
tamanhos <- c(80L, 120L, 160L, 184L)

cat("=== Discriminacao contra tamanho da amostra ===\n")
cat("  desfecho: casos | W reconstruida em cada subamostra | sorteios:", N_MC, "\n")
cat("  replicacoes por tamanho:", REPS, "\n\n")
cat(sprintf("  %-6s %8s %10s %10s %10s %10s %10s %8s %s\n",
            "n", "ligacoes", "k MQO", "k SAR", "D MQO", "D SAR", "razao",
            "dp", "razoes por replicacao"))

linhas <- list()
for (nn in tamanhos) {
  res <- list(); tentativas <- 0L
  set.seed(1000L + nn)
  while (length(res) < ifelse(nn == 184L, 1L, REPS) && tentativas < 60L) {
    tentativas <- tentativas + 1L
    idx <- if (nn == 184L) seq_len(184L) else sort(sample(184L, nn))
    r <- rodar(idx, N_MC, 20200731L + tentativas)
    if (!is.null(r)) res[[length(res) + 1]] <- r
  }
  if (length(res) == 0L) {
    cat(sprintf("  %-6d nenhuma das %d tentativas deu subamostra sem municipio isolado\n",
                nn, tentativas)); next }
  m <- function(f) mean(vapply(res, function(z) z[[f]], numeric(1)))
  # A razao por replicacao, e nao a razao das medias, e o que permite dizer se
  # a diferenca entre tamanhos supera a dispersao entre subamostras.
  razoes <- vapply(res, function(z) z$D_sar / z$D_ols, numeric(1))
  dp <- if (length(razoes) > 1L) sd(razoes) else NA_real_
  cat(sprintf("  %-6d %8.2f %10.2f %10.2f %10.4f %10.4f %10.3f %8s %s\n",
              nn, m("ligacoes"), m("k_ols"), m("k_sar"),
              m("D_ols"), m("D_sar"), mean(razoes),
              if (is.na(dp)) "-" else sprintf("%.3f", dp),
              paste(sprintf("%.2f", razoes), collapse = " ")))
  flush.console()
  linhas[[length(linhas) + 1]] <- data.frame(
    n = nn, reps = length(res), ligacoes = m("ligacoes"),
    k_ols = m("k_ols"), k_sar = m("k_sar"),
    D_ols = m("D_ols"), D_sar = m("D_sar"),
    razao = mean(razoes), dp_razao = dp,
    razoes = paste(sprintf("%.4f", razoes), collapse = ";"),
    stringsAsFactors = FALSE)
}

tab <- do.call(rbind, linhas)
cat("\n  correlacao entre n e a razao de D:", round(cor(tab$n, tab$razao), 4), "\n")
cat("\n  A associacao entre n e a razao de D nao isola o efeito do tamanho\n")
cat("  da amostra. A subamostragem tambem altera a estrutura espacial.\n")
cat("  Subamostrar reduz tambem a densidade da vizinhanca, entao os dois efeitos\n")
cat("  estao confundidos; a coluna de ligacoes por municipio deixa isso visivel.\n")

write.csv(tab, file.path(AQUI, "saida", "escala_n.csv"), row.names = FALSE)
cat("\nGravado em saida/escala_n.csv\n")
