# Passo 1: experimento com processo gerador conhecido.
#
# Tres cenarios, todos com X real e as 30 covariaveis como candidatas no BMA:
#
#   suaves   y = X_suaves  * beta + eps        rho = 0
#   rugosas  y = X_rugosas * beta + eps        rho = 0
#   controle y = (I - 0,6 W)^-1 (X_suaves * beta + eps)   rho = 0,6
#
# Os dois primeiros perguntam: o SAR destroi PIP de covariavel genuinamente
# causal? O terceiro e controle na direcao oposta: quando ha de fato defasagem
# espacial, o SAR tem que recuperar as verdadeiras e o MQO tem que falhar. Sem
# esse controle nao da para distinguir "o SAR e cego" de "o SAR esta certo".
#
# Reporta, por cenario e por metodo, a PIP media das 4 verdadeiras e das 26
# falsas, mais a taxa de acerto no limiar de 0,50.

source(if (file.exists("../caminhos.R")) "../caminhos.R" else "caminhos.R")

AQUI <- ESPACIAL
base <- readRDS(file.path(AQUI, "saida", "base_sar.rds"))
aux  <- readRDS(file.path(AQUI, "saida", "grade_sar.rds"))
hyp  <- readRDS(file.path(AQUI, "saida", "hyperg_aux.rds"))
tb   <- readRDS(file.path(AQUI, "saida", "tabela_f21_estendida.rds"))
cal  <- readRDS(file.path(AQUI, "saida", "calibracao_sinal.rds"))
des  <- readRDS(file.path(AQUI, "saida", "desenho_simulacao.rds"))

X <- base$X; W <- base$W; n <- base$n; covariaveis <- base$covariaveis
grade <- aux$grade; K <- length(covariaveis)
f21a <- hyp$f21a; a_par <- hyp$a_par; L <- hyp$L; J_MAX <- hyp$J_MAX
TAB <- tb$TAB; passo <- tb$passo; N_TAB <- tb$N_TAB; Z_MAX <- tb$Z_MAX
Xc <- sweep(X, 2, colMeans(X), "-")

RHO_CTRL <- 0.6
S_inv <- solve(diag(n) - RHO_CTRL * W)

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
  # Indexacao direta na matriz. Extrair a coluna inteira com TAB[, k+1] copia
  # 300.001 doubles (2,4 MB) a cada avaliacao: medido em 3,52 ms contra 0,02 ms
  # da indexacao direta, 176 vezes mais lento. Era o gargalo do experimento.
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

# Fabrica de funcoes de log-marginal para um y especifico.
montar <- function(y) {
  yc <- y - mean(y)
  Wy <- as.numeric(W %*% y); Wyc <- Wy - mean(Wy)
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

N_MC <- as.integer(Sys.getenv("MC3_SORTEIOS", "200000"))
BURN <- as.integer(N_MC * 0.2)
R <- as.integer(Sys.getenv("REPLICACOES", "10"))

# Betas calibrados em 19_calibrar_sinal.R para igualar a detectabilidade sob
# MQO: 0,81 nas suaves e 0,76 nas rugosas, ambas proximas da media de 0,72 das
# seis covariaveis selecionadas nos dados reais. Sem igualar isso, a comparacao
# entre cenarios confundiria o efeito do SAR com diferenca de forca do sinal.
b_suaves  <- cal$suaves$beta_escolhido
b_rugosas <- cal$rugosas$beta_escolhido
cenarios <- list(
  suaves   = list(vars = des$config$suaves$vars,  idx = des$config$suaves$idx,
                  beta = rep(b_suaves, 4),  rho = 0),
  rugosas  = list(vars = des$config$rugosas$vars, idx = des$config$rugosas$idx,
                  beta = rep(b_rugosas, 4), rho = 0),
  controle = list(vars = des$config$suaves$vars,  idx = des$config$suaves$idx,
                  beta = rep(b_suaves, 4),  rho = RHO_CTRL)
)
cat("  betas calibrados: suaves", b_suaves, "| rugosas", b_rugosas, "
")

cat("=== Experimento ===\n")
cat("  replicacoes por cenario:", R, "| sorteios por cadeia:", N_MC, "\n")
cat("  n =", n, "| candidatas no BMA:", K, "| verdadeiras: 4\n\n")

resultados <- list()
for (nome in names(cenarios)) {
  cn <- cenarios[[nome]]
  cat("### cenario", nome, "| rho do DGP =", cn$rho, "\n")
  cat("  verdadeiras:", paste(cn$vars, collapse = ", "), "\n"); flush.console()
  acumulado <- list(ols = matrix(NA_real_, R, K), sar = matrix(NA_real_, R, K),
                    tam_ols = numeric(R), tam_sar = numeric(R))
  t0 <- proc.time()[3]
  for (r in seq_len(R)) {
    set.seed(1000L * which(names(cenarios) == nome) + r)
    sinal <- as.numeric(Xc[, cn$idx, drop = FALSE] %*% cn$beta)
    y <- if (cn$rho == 0) sinal + rnorm(n) else as.numeric(S_inv %*% (sinal + rnorm(n)))
    fs <- montar(y)
    a <- mc3(fs, fs$ols, N_MC, BURN, 20200731L + r)
    b <- mc3(fs, fs$sar, N_MC, BURN, 20200731L + r)
    acumulado$ols[r, ] <- a$pip; acumulado$sar[r, ] <- b$pip
    acumulado$tam_ols[r] <- a$tam; acumulado$tam_sar[r] <- b$tam
    cat("    replicacao", r, "de", R, "\n"); flush.console()
  }
  cat("  tempo:", round(proc.time()[3] - t0, 1), "s\n")
  resultados[[nome]] <- c(acumulado, list(idx = cn$idx, vars = cn$vars, rho = cn$rho))

  verd <- cn$idx; fals <- setdiff(seq_len(K), verd)
  for (met in c("ols", "sar")) {
    M <- acumulado[[met]]
    cat(sprintf("  %-4s PIP media verdadeiras: %.4f | falsas: %.4f | razao: %.2f\n",
                toupper(met), mean(M[, verd]), mean(M[, fals]),
                mean(M[, verd]) / mean(M[, fals])))
    cat(sprintf("       verdadeiras acima de 0,50: %.1f%% | falsas acima de 0,50: %.1f%%\n",
                100 * mean(M[, verd] >= 0.5), 100 * mean(M[, fals] >= 0.5)))
    cat(sprintf("       tamanho medio do modelo: %.2f\n",
                mean(acumulado[[paste0("tam_", met)]])))
  }
  cat("\n"); flush.console()
}

cat("  avaliacoes de R2 fora da tabela:", fora_tab, "\n\n")

cat("=== PIP media por covariavel verdadeira, cenario a cenario ===\n")
for (nome in names(resultados)) {
  rr <- resultados[[nome]]
  cat(" ", nome, "(rho do DGP =", rr$rho, "):\n")
  for (j in seq_along(rr$idx)) {
    v <- rr$idx[j]
    cat(sprintf("    %-16s MQO %.4f | SAR %.4f | I de Moran da covariavel %.4f\n",
                covariaveis[v], mean(rr$ols[, v]), mean(rr$sar[, v]),
                des$moran[covariaveis[v]]))
  }
}

saveRDS(resultados, file.path(AQUI, "saida", "experimento_simulacao.rds"))
cat("\nGravado em saida/experimento_simulacao.rds\n")
