# Grade empirica completa: 3 desfechos x 2 matrizes de vizinhanca.
#
# Uma combinacao por execucao, escolhida por variavel de ambiente, para permitir
# rodar as seis em paralelo. Configuracao identica a do TCC: 2 cadeias,
# 1.000.000 de iteracoes, 200.000 de descarte, hyper-g UIP com Binomial-Beta.
#
#   DESFECHO = casos | obitos | letalidade
#   VIZINHOS = queen | rook
#
# Saida: saida/grade/<desfecho>_<vizinhos>.rds e um log por combinacao.

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
stopifnot(nrow(W) == n, max(abs(rowSums(W) - 1)) < 1e-9, max(abs(diag(W))) < 1e-12)

f21a <- hyp$f21a; a_par <- hyp$a_par; L <- hyp$L; J_MAX <- hyp$J_MAX
TAB <- tb$TAB; passo <- tb$passo; N_TAB <- tb$N_TAB; Z_MAX <- tb$Z_MAX

# Grade de rho e log-determinante, especificos desta matriz W.
mu <- Re(eigen(W, only.values = TRUE)$values)
n_grade <- 201L
rho <- seq(-0.999, 0.999, length.out = n_grade)
ld <- as.numeric(outer(rho, mu, function(r, m) log1p(-r * m)) %*% rep(1, length(mu)))
a0 <- 1.01
lprior_rho <- (a0 - 1) * log1p(rho) + (a0 - 1) * log1p(-rho) -
  (2 * a0 - 1) * log(2) - lbeta(a0, a0)
grade <- list(rho = rho, rho2 = rho^2, h = rho[2] - rho[1],
              ld_mais_prior = ld + lprior_rho,
              pesos = c(1, rep(c(4, 2), length.out = n_grade - 2), 1))

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
lm_sar <- function(e) {
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
mc3 <- function(marg, semente) {
  set.seed(semente)
  atual <- sort(sample(seq_len(K), sample(0:K, 1)))
  post <- marg(estat(atual)) + lpm[length(atual) + 1]
  aceitos <- 0L; inclusao <- numeric(K); guardados <- 0L; soma_k <- 0
  for (it in seq_len(N_MC)) {
    k <- length(atual); tipos <- tipos_disp(k)
    tipo <- tipos[sample.int(length(tipos), 1L)]
    fora <- setdiff(seq_len(K), atual)
    prop <- switch(tipo,
                   birth = sort(c(atual, fora[sample.int(length(fora), 1L)])),
                   death = atual[-sample.int(k, 1L)],
                   move  = sort(c(atual[-sample.int(k, 1L)],
                                  fora[sample.int(length(fora), 1L)])))
    post_p <- marg(estat(prop)) + lpm[length(prop) + 1]
    inv <- switch(tipo, birth = "death", death = "birth", move = "move")
    la <- post_p - post + log_q(length(prop), inv) - log_q(k, tipo)
    if (is.finite(la) && log(runif(1)) < la) {
      atual <- prop; post <- post_p; aceitos <- aceitos + 1L
    }
    if (it > BURN) {
      guardados <- guardados + 1L
      if (length(atual) > 0L) inclusao[atual] <- inclusao[atual] + 1
      soma_k <- soma_k + length(atual)
    }
  }
  list(pip = inclusao / guardados, taxa = aceitos / N_MC,
       tam = soma_k / guardados)
}

cat("=== ", DESFECHO, " x ", VIZINHOS, " ===\n", sep = "")
cat("  n =", n, "| K =", K, "| sorteios:", N_MC, "| descarte:", BURN, "\n")
mi_y <- spdep::moran.test(y, listw)$estimate[["Moran I statistic"]]
cat("  I de Moran do desfecho:", round(mi_y, 4), "\n")
cat("  autovalores de W: min", round(min(mu), 4), "max", round(max(mu), 4), "\n\n")

res <- list()
for (rot in c("MQO", "SAR")) {
  f <- if (rot == "MQO") lm_ols else lm_sar
  t0 <- proc.time()[3]
  a <- mc3(f, 20200731L); b <- mc3(f, 987654L)
  cat(sprintf("  %-4s %6.1f s | aceitacao %.3f/%.3f | tamanho %.3f/%.3f | discordancia %.4f\n",
              rot, proc.time()[3] - t0, a$taxa, b$taxa, a$tam, b$tam,
              max(abs(a$pip - b$pip))))
  res[[rot]] <- list(pip = (a$pip + b$pip) / 2, conv = max(abs(a$pip - b$pip)),
                     tam = (a$tam + b$tam) / 2)
}

tabela <- data.frame(covariavel = covariaveis,
                     mqo = res$MQO$pip, sar = res$SAR$pip,
                     stringsAsFactors = FALSE)
tabela$dif <- tabela$sar - tabela$mqo
tabela <- tabela[order(-tabela$mqo), ]

cat("\n  acima de 0,50 no MQO:",
    paste(tabela$covariavel[tabela$mqo >= 0.5], collapse = ", "), "\n")
cat("  acima de 0,50 no SAR:",
    ifelse(any(tabela$sar >= 0.5),
           paste(tabela$covariavel[tabela$sar >= 0.5], collapse = ", "), "nenhuma"), "\n")
cat("  tamanho medio: MQO", round(res$MQO$tam, 3), "| SAR", round(res$SAR$tam, 3), "\n")
cat("  linha de base (tamanho/K): MQO", round(res$MQO$tam / K, 4),
    "| SAR", round(res$SAR$tam / K, 4), "\n")
cat("  desvio-padrao das PIP: MQO", round(sd(res$MQO$pip), 4),
    "| SAR", round(sd(res$SAR$pip), 4), "\n")
cat("  R2 fora da tabela:", fora_tab, "\n")

dir.create(file.path(AQUI, "saida", "grade"), showWarnings = FALSE, recursive = TRUE)
saveRDS(list(desfecho = DESFECHO, vizinhos = VIZINHOS, tabela = tabela,
             mi_y = mi_y, res = res, K = K, n = n, mu = mu),
        file.path(AQUI, "saida", "grade", paste0(DESFECHO, "_", VIZINHOS, ".rds")))
cat("\nGravado.\n")
