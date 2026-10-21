# Otimizacao da log-marginal do SAR antes do MC3.
#
# A versao do 03_logmarginal_sar.R recalcula ln|I - rho W| para cada valor da
# grade a cada modelo avaliado. Esse termo NAO depende do modelo: depende so de
# W e da grade. Pre-computa-lo uma vez elimina a maior parte do custo.
#
# Tambem reduz a grade: a verificacao C mostrou que 101 pontos ja convergem ate
# 1e-13, entao os 1.001 pontos usados antes eram folga desnecessaria.
#
# Confere que a versao rapida devolve exatamente o mesmo numero da lenta.

source(if (file.exists("../caminhos.R")) "../caminhos.R" else "caminhos.R")

suppressPackageStartupMessages(library(spdep))

AQUI <- ESPACIAL
base <- readRDS(file.path(AQUI, "saida", "base_sar.rds"))

y <- base$y; X_total <- base$X; W <- base$W; mu <- base$mu
n <- base$n; covariaveis <- base$covariaveis

# ------------------------------------------------- 1. grade pre-computada
montar_grade <- function(n_grade = 201L, limites = c(-0.999, 0.999), a0 = 1.01) {
  if (n_grade %% 2L == 0L) n_grade <- n_grade + 1L
  rho <- seq(limites[1], limites[2], length.out = n_grade)
  # ln|I - rho W| = sum_i ln(1 - rho mu_i), calculado de uma vez para toda a grade
  ld <- as.numeric(outer(rho, mu, function(r, m) log1p(-r * m)) %*% rep(1, length(mu)))
  log_prior <- (a0 - 1) * log1p(rho) + (a0 - 1) * log1p(-rho) -
    (2 * a0 - 1) * log(2) - lbeta(a0, a0)
  pesos <- c(1, rep(c(4, 2), length.out = n_grade - 2), 1)
  list(rho = rho, rho2 = rho^2, h = rho[2] - rho[1],
       ld_mais_prior = ld + log_prior, pesos = pesos, n_grade = n_grade)
}

# ------------------------------------------------ 2. estatisticas do modelo
# Tudo o que depende do modelo se resume a seis escalares.
Wy <- as.numeric(W %*% y)
yc <- y - mean(y)
Wyc <- Wy - mean(Wy)
ta <- sum(yc^2); tb <- sum(yc * Wyc); tc <- sum(Wyc^2)

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

# ------------------------------------------------------ 3. versao rapida
log_marginal_rapida <- function(est, g, grade) {
  rss <- est["a"] - 2 * grade$rho * est["b"] + grade$rho2 * est["c"]
  tss <- ta - 2 * grade$rho * tb + grade$rho2 * tc
  lf <- grade$ld_mais_prior -
    ((n - 1) / 2) * log((g / (1 + g)) * rss + (1 / (1 + g)) * tss)
  m <- max(lf)
  unname(-(est["k"] / 2) * log(1 + g) + m +
           log(grade$h / 3 * sum(grade$pesos * exp(lf - m))))
}

# --------------------------------------- 4. versao lenta, para conferencia
log_marginal_lenta <- function(indices, g, n_grade = 201L, a0 = 1.01) {
  if (n_grade %% 2L == 0L) n_grade <- n_grade + 1L
  rho <- seq(-0.999, 0.999, length.out = n_grade)
  est <- estatisticas(indices)
  rss <- est["a"] - 2 * rho * est["b"] + rho^2 * est["c"]
  tss <- ta - 2 * rho * tb + rho^2 * tc
  ld <- vapply(rho, function(r) sum(log1p(-r * mu)), numeric(1))
  lp <- (a0 - 1) * log1p(rho) + (a0 - 1) * log1p(-rho) -
    (2 * a0 - 1) * log(2) - lbeta(a0, a0)
  lf <- -(est["k"] / 2) * log(1 + g) + ld -
    ((n - 1) / 2) * log((g / (1 + g)) * rss + (1 / (1 + g)) * tss) + lp
  m <- max(lf)
  h <- rho[2] - rho[1]
  pesos <- c(1, rep(c(4, 2), length.out = n_grade - 2), 1)
  unname(m + log(h / 3 * sum(pesos * exp(lf - m))))
}

# ------------------------------------------------------- 5. equivalencia
g_uip <- n
grade <- montar_grade(201L)

set.seed(20200731L)
amostras <- replicate(50, sort(sample(seq_along(covariaveis), sample(0:12, 1))),
                      simplify = FALSE)

difs <- vapply(amostras, function(idx) {
  log_marginal_rapida(estatisticas(idx), g_uip, grade) - log_marginal_lenta(idx, g_uip)
}, numeric(1))

cat("=== Equivalencia entre a versao rapida e a lenta ===\n")
cat("  modelos testados:", length(amostras), "\n")
cat("  maior diferenca absoluta:", format(max(abs(difs)), scientific = TRUE), "\n")
cat("  veredito:", if (max(abs(difs)) < 1e-12) "APROVADO" else "REPROVADO", "\n\n")

# ------------------------------------------------------------ 6. velocidade
cat("=== Velocidade ===\n")
medir <- function(f, r = 300L) {
  set.seed(1L)
  idxs <- replicate(r, sample(seq_along(covariaveis), sample(1:15, 1)), simplify = FALSE)
  t0 <- proc.time()[3]
  for (idx in idxs) invisible(f(idx))
  (proc.time()[3] - t0) / r
}

t_lenta <- medir(function(idx) log_marginal_lenta(idx, g_uip))
t_rapida <- medir(function(idx) log_marginal_rapida(estatisticas(idx), g_uip, grade))

velocidades <- data.frame(
  versao = c("lenta (log-det por modelo)", "rapida (log-det pre-computado)"),
  ms_por_modelo = c(t_lenta, t_rapida) * 1000,
  min_para_250k = c(t_lenta, t_rapida) * 250000 / 60
)
print(transform(velocidades, ms_por_modelo = round(ms_por_modelo, 3),
                min_para_250k = round(min_para_250k, 1)), row.names = FALSE)
cat("\n  ganho:", round(t_lenta / t_rapida, 1), "vezes\n\n")

# Onde sobra tempo agora: a decomposicao QR de X_j, que e inevitavel por modelo.
cat("=== Decomposicao do custo da versao rapida ===\n")
set.seed(1L)
idxs <- replicate(300L, sample(seq_along(covariaveis), sample(1:15, 1)), simplify = FALSE)
t0 <- proc.time()[3]; for (idx in idxs) invisible(estatisticas(idx))
t_est <- (proc.time()[3] - t0) / 300
ests <- lapply(idxs, estatisticas)
t0 <- proc.time()[3]; for (e in ests) invisible(log_marginal_rapida(e, g_uip, grade))
t_int <- (proc.time()[3] - t0) / 300
cat("  estatisticas do modelo (QR):", round(t_est * 1000, 3), "ms\n")
cat("  integracao em rho          :", round(t_int * 1000, 3), "ms\n\n")

saveRDS(list(grade = grade, ta = ta, tb = tb, tc = tc, Wy = Wy, yc = yc, Wyc = Wyc),
        file.path(AQUI, "saida", "grade_sar.rds"))

cat("Grade e constantes gravadas em saida/grade_sar.rds\n")
