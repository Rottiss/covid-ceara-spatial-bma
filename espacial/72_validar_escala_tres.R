# Etapa 1 da cadeia de tres vias: as tres verossimilhancas marginais estao na
# mesma escala?
#
# Um fator de Bayes entre MQO, SAR e SEM so faz sentido se as tres marginais
# omitirem as mesmas constantes de normalizacao. As formas sao identicas nos
# scripts 33 e 31, mas o SEM passa por uma correcao de centragem que o SAR nao
# tem, e um deslocamento de escala ali sairia como P(SEM | y) viesado sem que
# nada no resultado denunciasse.
#
# O teste e exato e nao usa MCMC: em rho = 0 o integrando do SAR tem de devolver
# a marginal do MQO, e em lambda = 0 o do SEM tambem, a menos do termo de priori
# do parametro espacial, que e conhecido. Enumeracao completa de 2^10 = 1024
# modelos sobre um subconjunto de 10 covariaveis.
#
# Configuracao do TCC: hiper-g com calibracao UIP. Zero é ponto exato das duas
# grades, porque ambas sao simetricas com numero impar de pontos.
#
# So le. Nao grava nada.

source(if (file.exists("../caminhos.R")) "../caminhos.R" else "caminhos.R")

suppressPackageStartupMessages(library(spdep))
AQUI <- ESPACIAL
RAIZ <- ESTIMACOES
hyp <- readRDS(file.path(AQUI, "saida", "hyperg_aux.rds"))
tb  <- readRDS(file.path(AQUI, "saida", "tabela_f21_estendida.rds"))
f21a <- hyp$f21a; a_par <- hyp$a_par; L <- hyp$L; J_MAX <- hyp$J_MAX
TAB <- tb$TAB; passo <- tb$passo; N_TAB <- tb$N_TAB; Z_MAX <- tb$Z_MAX

DESFECHO <- "casos"; COLUNA <- "z_log_casos100k"
dados <- read.csv(file.path(RAIZ, "resultados/estimacoes/matrizes_estimacao/principal",
                            paste0(DESFECHO, ".csv")),
                  stringsAsFactors = FALSE, colClasses = c(codigo_ibge = "character"))
todas <- setdiff(names(dados), c("codigo_ibge", "municipio", COLUNA))
set.seed(20200731L)
sel <- sort(sample(seq_along(todas), 10))
covs <- todas[sel]
y <- as.numeric(dados[[COLUNA]])
X <- as.matrix(dados[, covs]); n <- length(y); K <- ncol(X)
Xc <- sweep(X, 2, colMeans(X), "-")
listw <- readRDS(file.path(RAIZ, "resultados/diagnostico_espacial", "pesos_queen.rds"))
W <- spdep::listw2mat(listw); dimnames(W) <- NULL
stopifnot(nrow(W) == n, max(abs(rowSums(W) - 1)) < 1e-9)

cat("covariaveis do subconjunto:", paste(covs, collapse = ", "), "\n")
cat("modelos enumerados:", 2^K, "\n\n")

# --------------------------------------------------- hipergeometrica
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
a0 <- 1.01
grade <- function(ngp) {
  v <- seq(-0.999, 0.999, length.out = ngp)
  ld <- as.numeric(outer(v, mu, function(r, m) log1p(-r * m)) %*% rep(1, length(mu)))
  lpr <- (a0 - 1) * log1p(v) + (a0 - 1) * log1p(-v) -
    (2 * a0 - 1) * log(2) - lbeta(a0, a0)
  j0 <- which(abs(v) < 1e-12)
  stopifnot(length(j0) == 1L)
  list(v = v, ld = ld, lpr = lpr, j0 = j0)
}
gr <- grade(201L)   # rho, como no 33
gl <- grade(101L)   # lambda, como no 31
cat(sprintf("zero e ponto exato das grades: rho[%d]=%g, lambda[%d]=%g\n",
            gr$j0, gr$v[gr$j0], gl$j0, gl$v[gl$j0]))
cat(sprintf("log|I - 0*W| = %g (rho) e %g (lambda), tem de ser zero\n\n",
            gr$ld[gr$j0], gl$ld[gl$j0]))

# --------------------------------------------------- MQO e SAR, forma do 33
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
# integrando do SAR, ponto a ponto da grade
lf_sar <- function(e) {
  k <- unname(e["k"])
  rss <- e["a"] - 2 * gr$v * e["b"] + gr$v^2 * e["c"]
  tss <- ta - 2 * gr$v * tbb + gr$v^2 * tcc
  lf <- gr$ld + gr$lpr + (1 - n) / 2 * log(tss)
  if (k > 0L) lf <- lf + log(f21a - 2) - log(k + f21a - 2) +
      log_f21(pmax(0, 1 - rss / tss), k)
  lf
}

# --------------------------------------------------- SEM, forma do 31
WX <- W %*% X
A00 <- crossprod(X); A01 <- crossprod(X, WX); A11 <- crossprod(WX)
b0 <- as.numeric(crossprod(X, y));  b1 <- as.numeric(crossprod(X, Wy))
c0 <- as.numeric(crossprod(WX, y)); c1 <- as.numeric(crossprod(WX, Wy))
s0 <- colSums(X); s1 <- colSums(WX)
yy <- sum(y^2); yWy <- sum(y * Wy); WyWy <- sum(Wy^2)
sy <- sum(y); sWy <- sum(Wy)
ys_ys <- yy - 2 * gl$v * yWy + gl$v^2 * WyWy
soma_ys <- sy - gl$v * sWy
tss_s <- ys_ys - n * (soma_ys / n)^2
rss_lambda <- function(idx) {
  if (length(idx) == 0L) return(tss_s)
  A00i <- A00[idx, idx, drop = FALSE]; A01i <- A01[idx, idx, drop = FALSE]
  A11i <- A11[idx, idx, drop = FALSE]
  b0i <- b0[idx]; b1i <- b1[idx]; c0i <- c0[idx]; c1i <- c1[idx]
  s0i <- s0[idx]; s1i <- s1[idx]
  sim <- A01i + t(A01i); soma_bc <- b1i + c0i
  out <- numeric(length(gl$v))
  for (j in seq_along(gl$v)) {
    Lj <- gl$v[j]
    Xti <- s0i - Lj * s1i
    XtX <- A00i - Lj * sim + Lj^2 * A11i - tcrossprod(Xti) / n
    Xty <- b0i - Lj * soma_bc + Lj^2 * c1i - Xti * (soma_ys[j] / n)
    beta <- tryCatch(solve(XtX, Xty), error = function(e) NULL)
    out[j] <- if (is.null(beta)) NA_real_ else tss_s[j] - sum(beta * Xty)
  }
  out
}
lf_sem <- function(idx) {
  k <- length(idx)
  rss <- rss_lambda(idx)
  lf <- gl$ld + gl$lpr + (1 - n) / 2 * log(tss_s)
  if (k > 0L) lf <- lf + log(f21a - 2) - log(k + f21a - 2) +
      log_f21(pmax(0, pmin(1 - rss / tss_s, 1 - 1e-12)), k)
  lf
}

# --------------------------------------------------- enumeracao
modelos <- lapply(0:(2^K - 1), function(m) which(bitwAnd(m, 2^(0:(K - 1))) > 0))
d_sar <- d_sem <- numeric(length(modelos))
t0 <- proc.time()[3]
for (i in seq_along(modelos)) {
  idx <- modelos[[i]]
  ref <- lm_ols(estat(idx))
  d_sar[i] <- lf_sar(estat(idx))[gr$j0] - gr$lpr[gr$j0] - ref
  d_sem[i] <- lf_sem(idx)[gl$j0]        - gl$lpr[gl$j0] - ref
}
cat(sprintf("enumeracao: %.1f s\n\n", proc.time()[3] - t0))

cat("=== SAR em rho = 0 contra a marginal do MQO ===\n")
cat(sprintf("  diferenca maxima absoluta: %.3e | mediana: %.3e\n",
            max(abs(d_sar)), median(abs(d_sar))))
cat("=== SEM em lambda = 0 contra a marginal do MQO ===\n")
cat(sprintf("  diferenca maxima absoluta: %.3e | mediana: %.3e\n",
            max(abs(d_sem)), median(abs(d_sem))))
cat(sprintf("\navaliacoes fora da tabela de f21: %d\n", fora_tab))

TOL <- 1e-8
ok <- max(abs(d_sar)) < TOL && max(abs(d_sem)) < TOL
cat(sprintf("\nescala comum (tolerancia %.0e): %s\n", TOL, if (ok) "SIM" else "NAO"))
if (!ok) {
  pior <- which.max(pmax(abs(d_sar), abs(d_sem)))
  cat(sprintf("  pior modelo: k = %d, d_sar = %.3e, d_sem = %.3e\n",
              length(modelos[[pior]]), d_sar[pior], d_sem[pior]))
  stop("As marginais nao estao na mesma escala.")
}
