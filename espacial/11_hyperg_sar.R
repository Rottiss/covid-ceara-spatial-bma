# Log-marginal do SAR sob a priori do TCC (hyper-g UIP + Binomial-Beta).
#
# CONSTRUCAO. Condicional a rho, o SAR e um modelo linear nos dados
# transformados y*(rho) = (I - rho W) y. Entao a marginal condicional e
# exatamente a hyper-g validada em 10_validar_hyperg.R, aplicada a y*(rho), mais
# o jacobiano da transformacao:
#
#   log p(y | M, rho) = log|I - rho W|
#                       + (1-N)/2 * log TSS(rho)
#                       + log(a-2) - log(k+a-2)
#                       + log 2F1( (N-1)/2, 1; (a+k)/2; R2(rho) )
#
# e integra-se em rho por Simpson, com priori Beta(1,01; 1,01) em [-1,1].
#
# HIPERGEOMETRICA EM ESPACO LOGARITMICO. Os valores de 2F1 chegam a 1e35 nesta
# amostra e podem estourar em R2 alto. A serie e acumulada em log, com
# log-sum-exp, e vetorizada sobre a grade de rho: os log-quocientes
# log(a+j) - log(c+j) nao dependem de rho, so de k, entao sao pre-computados
# uma vez por k.

source(if (file.exists("../caminhos.R")) "../caminhos.R" else "caminhos.R")

AQUI <- ESPACIAL
base <- readRDS(file.path(AQUI, "saida", "base_sar.rds"))
aux <- readRDS(file.path(AQUI, "saida", "grade_sar.rds"))

y <- base$y; X_total <- base$X; n <- base$n; covariaveis <- base$covariaveis
grade <- aux$grade; ta <- aux$ta; tb <- aux$tb; tc <- aux$tc
yc <- aux$yc; Wyc <- aux$Wyc
K <- length(covariaveis)

f21a <- 2 + 2 / n
a_par <- (n - 1) / 2
J_MAX <- 4000L

# L[[k+1]] = cumsum_{i<j} [log(a+i) - log(c_k+i)], j = 1..J_MAX
L <- lapply(0:K, function(k) {
  cc <- (f21a + k) / 2
  cumsum(log(a_par + 0:(J_MAX - 1)) - log(cc + 0:(J_MAX - 1)))
})

# log 2F1(a,1;c;z) para um vetor z e um k, por log-sum-exp sobre os termos.
log_f21_vec <- function(z, k) {
  if (any(z >= 1) || any(z < 0)) return(rep(NA_real_, length(z)))
  lz <- log(z)
  # Onde os termos param de crescer: a razao entre termos consecutivos e
  # z*(a+j)/(c+j), que cai abaixo de 1 quando j > (z*a - c)/(1 - z). Truncar
  # antes desse pico subestima a soma. A primeira versao deste script tinha o
  # sinal invertido nessa conta e cortava a serie no meio da subida, errando o
  # log da soma em ate 12 unidades; a conferencia contra a serie direta pegou.
  cc <- (f21a + k) / 2
  j_pico <- max(0, ceiling(max((z * a_par - cc) / pmax(1 - z, 1e-12))))
  J <- min(J_MAX, max(300L, as.integer(j_pico * 4)))
  lt <- outer(lz, seq_len(J)) + rep(L[[k + 1]][seq_len(J)], each = length(z))
  lt <- cbind(0, lt)                      # termo j = 0 vale 1
  m <- apply(lt, 1, max)
  list(valor = m + log(rowSums(exp(lt - m))),
       ultimo_relativo = max(exp(lt[, J + 1] - m)))
}

# Verificacao contra a versao direta em espaco natural, onde ela nao estoura.
f21_direto <- function(z, k) {
  cc <- (f21a + k) / 2
  termo <- 1; soma <- 1; j <- 0L
  repeat {
    termo <- termo * z * (a_par + j) / (cc + j)
    soma <- soma + termo; j <- j + 1L
    if (abs(termo) < 1e-15 * abs(soma) || j >= 100000L) break
  }
  soma
}

cat("=== Conferencia: log-sum-exp vs. serie direta ===\n")
zs <- c(0.05, 0.2, 0.4, 0.6, 0.8, 0.9)
pior <- 0
for (k in c(1L, 5L, 15L, 30L)) {
  lv <- log_f21_vec(zs, k)$valor
  dv <- log(vapply(zs, f21_direto, numeric(1), k = k))
  d <- max(abs(lv - dv))
  pior <- max(pior, d)
  cat(sprintf("  k = %2d  maior |diferenca| em log: %.3e\n", k, d))
}
cat("  pior caso:", format(pior, scientific = TRUE), "\n")
if (pior > 1e-9) {
  stop("A hipergeometrica em log nao reproduz a serie direta. Nao prosseguir.")
}
cat("  veredito: APROVADO\n\n")

# --------------------------------------------------- estatisticas do modelo
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

log_prior_modelo <- function(k) lgamma(1 + k) + lgamma(1 + K - k)

# ----------------------------------------------------- MQO com hyper-g
lm_ols_hyperg <- function(est) {
  k <- unname(est["k"])
  null_lik <- (1 - n) / 2 * log(ta)
  if (k == 0L) return(null_lik)
  R2 <- 1 - est["a"] / ta
  unname(null_lik + log(f21a - 2) - log(k + f21a - 2) +
           log_f21_vec(R2, k)$valor)
}

# ----------------------------------------------------- SAR com hyper-g
lm_sar_hyperg <- function(est) {
  k <- unname(est["k"])
  rss <- est["a"] - 2 * grade$rho * est["b"] + grade$rho2 * est["c"]
  tss <- ta - 2 * grade$rho * tb + grade$rho2 * tc
  if (any(rss <= 0) || any(tss <= 0)) return(-Inf)
  R2 <- 1 - rss / tss
  R2[R2 < 0] <- 0
  lf <- grade$ld_mais_prior + (1 - n) / 2 * log(tss)
  if (k > 0L) {
    lf <- lf + log(f21a - 2) - log(k + f21a - 2) + log_f21_vec(R2, k)$valor
  }
  m <- max(lf)
  if (!is.finite(m)) return(-Inf)
  unname(m + log(grade$h / 3 * sum(grade$pesos * exp(lf - m))))
}

# ------------------------------------- verificacao: SAR em rho = 0 vs. MQO
cat("=== Verificacao: a densidade do SAR em rho = 0 bate com o MQO ===\n")
i_zero <- which.min(abs(grade$rho))
cat("  rho da grade mais proximo de zero:", grade$rho[i_zero], "\n")
set.seed(20200731L)
amostra <- replicate(20, sort(sample(K, sample(0:8, 1))), simplify = FALSE)
difs <- vapply(amostra, function(idx) {
  est <- estatisticas(idx)
  k <- unname(est["k"])
  rss <- est["a"] - 2 * grade$rho * est["b"] + grade$rho2 * est["c"]
  tss <- ta - 2 * grade$rho * tb + grade$rho2 * tc
  lf <- (1 - n) / 2 * log(tss[i_zero])
  if (k > 0L) {
    R2 <- 1 - rss[i_zero] / tss[i_zero]
    lf <- lf + log(f21a - 2) - log(k + f21a - 2) + log_f21_vec(R2, k)$valor
  }
  lf - lm_ols_hyperg(est)
}, numeric(1))
cat("  modelos testados:", length(amostra), "\n")
cat("  maior diferenca:", format(max(abs(difs)), scientific = TRUE), "\n")
cat("  veredito:", if (max(abs(difs)) < 1e-9) "APROVADO" else "REPROVADO", "\n\n")

# --------------------------------------------------------------- custo
cat("=== Custo por modelo ===\n")
set.seed(1)
idxs <- replicate(200L, sample(K, sample(1:15, 1)), simplify = FALSE)
ests <- lapply(idxs, estatisticas)
t0 <- proc.time()[3]; for (e in ests) invisible(lm_ols_hyperg(e))
t_ols <- (proc.time()[3] - t0) / 200
t0 <- proc.time()[3]; for (e in ests) invisible(lm_sar_hyperg(e))
t_sar <- (proc.time()[3] - t0) / 200
t0 <- proc.time()[3]; for (idx in idxs) invisible(estatisticas(idx))
t_est <- (proc.time()[3] - t0) / 200
cat("  estatisticas (QR):", round(t_est * 1000, 3), "ms\n")
cat("  MQO hyper-g      :", round(t_ols * 1000, 3), "ms\n")
cat("  SAR hyper-g      :", round(t_sar * 1000, 3), "ms\n")
cat("  projecao SAR para 100.000 sorteios:",
    round((t_est + t_sar) * 100000 / 60, 1), "min por cadeia\n\n")

saveRDS(list(f21a = f21a, a_par = a_par, L = L, J_MAX = J_MAX),
        file.path(AQUI, "saida", "hyperg_aux.rds"))
cat("Auxiliares gravados.\n")
