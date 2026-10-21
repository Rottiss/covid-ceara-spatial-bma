# Cenario realista: DGP com estrutura espacial igual a do desfecho observado.
#
# Os cenarios anteriores tinham rho = 0 e y com I de Moran de 0,13, contra 0,4834
# do desfecho real. O rho que reproduz a estrutura real e 0,5285 (calibrado em
# 21_calibrar_rho_e_colinearidade.R).
#
# Compare PIPs under a calibrated DGP with known nonzero coefficients.
# These scenarios do not establish why PIPs differ in the observed data.
#
# Roda tambem com 8 verdadeiras em vez de 4, porque o modelo a posteriori dos
# dados reais tem tamanho medio de 13 e o cenario de 4 pode ser esparso demais.

source(if (file.exists("../caminhos.R")) "../caminhos.R" else "caminhos.R")

suppressPackageStartupMessages(library(spdep))

AQUI <- ESPACIAL
base <- readRDS(file.path(AQUI, "saida", "base_sar.rds"))
aux  <- readRDS(file.path(AQUI, "saida", "grade_sar.rds"))
hyp  <- readRDS(file.path(AQUI, "saida", "hyperg_aux.rds"))
tb   <- readRDS(file.path(AQUI, "saida", "tabela_f21_estendida.rds"))
des  <- readRDS(file.path(AQUI, "saida", "desenho_simulacao.rds"))
crho <- readRDS(file.path(AQUI, "saida", "calibracao_rho.rds"))

X <- base$X; W <- base$W; n <- base$n; covariaveis <- base$covariaveis
listw <- base$listw; grade <- aux$grade; K <- length(covariaveis)
f21a <- hyp$f21a; a_par <- hyp$a_par; L <- hyp$L; J_MAX <- hyp$J_MAX
TAB <- tb$TAB; passo <- tb$passo; N_TAB <- tb$N_TAB; Z_MAX <- tb$Z_MAX
Xc <- sweep(X, 2, colMeans(X), "-")

# O rho do DGP. O padrao e o calibrado para reproduzir o I de Moran observado
# (0,5285). A posteriori de rho sob media de modelos, estimada em
# 33_probabilidade_espacial.R, tem media 0,70, mais alta que a calibrada: casar
# o I de Moran de y nao e o mesmo que casar o rho posterior, porque parte da
# suavidade espacial de y vem das proprias covariaveis. RHO_DGP permite rodar o
# cenario no rho que os dados indicam, que e o caso mais adverso.
RHO <- as.numeric(Sys.getenv("RHO_DGP", crho$rho_alvo))
S_inv <- solve(diag(n) - RHO * W)

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

lpm <- vapply(0:K, function(k) lgamma(1 + k) + lgamma(1 + K - k), numeric(1))
tipos_disp <- function(k) {
  t <- character(0)
  if (k < K) t <- c(t, "birth"); if (k > 0L) t <- c(t, "death")
  if (k > 0L && k < K) t <- c(t, "move"); t
}
log_q <- function(k, tipo) {
  -log(length(tipos_disp(k))) -
    log(switch(tipo, birth = K - k, death = k, move = k * (K - k)))
}

montar <- function(y) {
  yc <- y - mean(y); Wy <- as.numeric(W %*% y); Wyc <- Wy - mean(Wy)
  ta <- sum(yc^2); tbb <- sum(yc * Wyc); tcc <- sum(Wyc^2)
  est <- function(idx) {
    if (length(idx) == 0L) { eo <- yc; ed <- Wyc } else {
      Q <- qr.Q(qr(Xc[, idx, drop = FALSE]))
      eo <- yc - Q %*% crossprod(Q, yc); ed <- Wyc - Q %*% crossprod(Q, Wyc)
    }
    eo <- as.numeric(eo); ed <- as.numeric(ed)
    c(k = length(idx), a = sum(eo^2), b = sum(eo * ed), c = sum(ed^2))
  }
  ols <- function(e) {
    k <- unname(e["k"]); nl <- (1 - n) / 2 * log(ta)
    if (k == 0L) return(nl)
    unname(nl + log(f21a - 2) - log(k + f21a - 2) + log_f21(1 - e["a"] / ta, k))
  }
  sar <- function(e) {
    k <- unname(e["k"])
    rss <- e["a"] - 2 * grade$rho * e["b"] + grade$rho2 * e["c"]
    tss <- ta - 2 * grade$rho * tbb + grade$rho2 * tcc
    if (any(rss <= 0) || any(tss <= 0)) return(-Inf)
    lf <- grade$ld_mais_prior + (1 - n) / 2 * log(tss)
    if (k > 0L) lf <- lf + log(f21a - 2) - log(k + f21a - 2) +
        log_f21(pmax(0, 1 - rss / tss), k)
    m <- max(lf); if (!is.finite(m)) return(-Inf)
    unname(m + log(grade$h / 3 * sum(grade$pesos * exp(lf - m))))
  }
  list(est = est, ols = ols, sar = sar)
}

mc3 <- function(fs, marginal, n_sorteios, burn, semente) {
  set.seed(semente)
  atual <- sort(sample(seq_len(K), sample(0:K, 1)))
  post <- marginal(fs$est(atual)) + lpm[length(atual) + 1]
  inclusao <- numeric(K); guardados <- 0L; soma_k <- 0
  for (it in seq_len(n_sorteios)) {
    k <- length(atual); tipos <- tipos_disp(k)
    tipo <- tipos[sample.int(length(tipos), 1L)]
    fora <- setdiff(seq_len(K), atual)
    prop <- switch(tipo,
                   birth = sort(c(atual, fora[sample.int(length(fora), 1L)])),
                   death = atual[-sample.int(k, 1L)],
                   move  = sort(c(atual[-sample.int(k, 1L)],
                                  fora[sample.int(length(fora), 1L)])))
    post_p <- marginal(fs$est(prop)) + lpm[length(prop) + 1]
    inv <- switch(tipo, birth = "death", death = "birth", move = "move")
    la <- post_p - post + log_q(length(prop), inv) - log_q(k, tipo)
    if (is.finite(la) && log(runif(1)) < la) { atual <- prop; post <- post_p }
    if (it > burn) {
      guardados <- guardados + 1L
      if (length(atual) > 0L) inclusao[atual] <- inclusao[atual] + 1
      soma_k <- soma_k + length(atual)
    }
  }
  list(pip = inclusao / guardados, tam = soma_k / guardados)
}

# Conjuntos verdadeiros: as 4 e as 8 continuas de maior I de Moran.
cont <- covariaveis[grepl("^z_", covariaveis)]
mc_ord <- sort(des$moran[cont], decreasing = TRUE)
conjuntos <- list(quatro = names(mc_ord)[1:4], oito = names(mc_ord)[1:8])
# Permite rodar um conjunto por vez, para nao perder o outro se a execucao cair.
so <- Sys.getenv("CONJUNTO", "")
if (nzchar(so)) conjuntos <- conjuntos[so]
ARQ <- Sys.getenv("ARQ_SAIDA", "cenario_realista.rds")

# Calibra beta para R2 do modelo verdadeiro proximo de 0,35 sob este rho.
calibrar <- function(idx, alvo = 0.35, reps = 40L) {
  f <- function(b) {
    set.seed(555L)
    mean(replicate(reps, {
      y <- as.numeric(S_inv %*% (as.numeric(Xc[, idx, drop = FALSE] %*% rep(b, length(idx))) + rnorm(n)))
      yc <- y - mean(y); Q <- qr.Q(qr(Xc[, idx, drop = FALSE]))
      sum((crossprod(Q, yc))^2) / sum(yc^2)
    })) - alvo
  }
  uniroot(f, c(0.02, 1.2))$root
}

N_MC <- as.integer(Sys.getenv("MC3_SORTEIOS", "120000"))
BURN <- as.integer(N_MC * 0.2)
R <- as.integer(Sys.getenv("REPLICACOES", "10"))

cat("=== Cenario realista ===\n")
cat("  rho do DGP:", round(RHO, 4), "(reproduz I de Moran de",
    round(crho$mi_real, 4), ")\n")
cat("  replicacoes:", R, "| sorteios por cadeia:", N_MC, "\n\n")

resultados <- list()
for (nome in names(conjuntos)) {
  vs <- conjuntos[[nome]]; idx <- match(vs, covariaveis)
  b <- calibrar(idx)
  cat("### ", nome, ": ", length(idx), " verdadeiras, beta = ", round(b, 4), "\n", sep = "")
  cat("  ", paste(vs, collapse = ", "), "\n", sep = "")
  set.seed(999L)
  mi <- mean(replicate(20, {
    y <- as.numeric(S_inv %*% (as.numeric(Xc[, idx, drop = FALSE] %*% rep(b, length(idx))) + rnorm(n)))
    spdep::moran.test(y, listw)$estimate[["Moran I statistic"]]
  }))
  cat("  I de Moran do y simulado:", round(mi, 4), "\n"); flush.console()

  M_ols <- matrix(NA_real_, R, K); M_sar <- matrix(NA_real_, R, K)
  tam <- matrix(NA_real_, R, 2)
  t0 <- proc.time()[3]
  for (r in seq_len(R)) {
    set.seed(30000L + 100L * which(names(conjuntos) == nome) + r)
    y <- as.numeric(S_inv %*% (as.numeric(Xc[, idx, drop = FALSE] %*% rep(b, length(idx))) + rnorm(n)))
    fs <- montar(y)
    a <- mc3(fs, fs$ols, N_MC, BURN, 20200731L + r)
    bb <- mc3(fs, fs$sar, N_MC, BURN, 20200731L + r)
    M_ols[r, ] <- a$pip; M_sar[r, ] <- bb$pip
    tam[r, ] <- c(a$tam, bb$tam)
  }
  cat("  tempo:", round(proc.time()[3] - t0, 1), "s\n")

  fals <- setdiff(seq_len(K), idx)
  for (met in list(c("MQO", "M_ols"), c("SAR", "M_sar"))) {
    M <- get(met[2])
    cat(sprintf("  %-4s verdadeiras: %.4f | falsas: %.4f | razao: %.2f | acima de 0,50: %.0f%% vs %.0f%%\n",
                met[1], mean(M[, idx]), mean(M[, fals]),
                mean(M[, idx]) / mean(M[, fals]),
                100 * mean(M[, idx] >= 0.5), 100 * mean(M[, fals] >= 0.5)))
  }
  cat(sprintf("       tamanho medio: MQO %.2f | SAR %.2f\n", mean(tam[, 1]), mean(tam[, 2])))
  cat(sprintf("       desvio das PIP: MQO %.4f | SAR %.4f\n",
              mean(apply(M_ols, 1, sd)), mean(apply(M_sar, 1, sd))))
  cat("\n"); flush.console()
  resultados[[nome]] <- list(ols = M_ols, sar = M_sar, idx = idx, vars = vs,
                             beta = b, mi = mi, tam = tam)
  # Salva a cada conjunto concluido: uma interrupcao no meio do segundo nao
  # apaga o primeiro.
  saveRDS(list(resultados = resultados, rho = RHO),
          file.path(AQUI, "saida", ARQ))
}

cat("  R2 fora da tabela:", fora_tab, "\n")
saveRDS(list(resultados = resultados, rho = RHO),
        file.path(AQUI, "saida", ARQ))
cat("\nGravado em saida/", ARQ, "\n", sep = "")
