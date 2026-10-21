# Prepara e valida uma tabela de interpolacao para log 2F1, antes de rodar o SAR
# com 1.000.000 de iteracoes.
#
# MOTIVACAO. No SAR, log 2F1 e avaliada em cada um dos pontos da grade de rho, a
# cada modelo proposto. A serie tem centenas de termos, e montar a matriz
# (pontos de rho) x (termos) domina o custo: 4,3 ms por modelo, o que dava
# 7.478 s para 400.000 sorteios em duas cadeias.
#
# Mas log 2F1((N-1)/2, 1; (a+k)/2; z) depende apenas de z e de k. Como k vai de 0
# a 30 e z fica num intervalo limitado, cabe uma tabela: calcula-se a serie uma
# vez numa grade fina de z, para cada k, e depois so se interpola.
#
# Este script mede o intervalo de z que o SAR realmente visita, monta a tabela e
# confere o erro de interpolacao contra a serie exata. So depois disso vale rodar.

source(if (file.exists("../caminhos.R")) "../caminhos.R" else "caminhos.R")

AQUI <- ESPACIAL
base <- readRDS(file.path(AQUI, "saida", "base_sar.rds"))
aux <- readRDS(file.path(AQUI, "saida", "grade_sar.rds"))
hyp <- readRDS(file.path(AQUI, "saida", "hyperg_aux.rds"))

X_total <- base$X; n <- base$n; covariaveis <- base$covariaveis
grade <- aux$grade; ta <- aux$ta; tb <- aux$tb; tc <- aux$tc
yc <- aux$yc; Wyc <- aux$Wyc
K <- length(covariaveis)
f21a <- hyp$f21a; a_par <- hyp$a_par; L <- hyp$L; J_MAX <- hyp$J_MAX

cat("=== Grade de rho em uso ===\n")
cat("  pontos:", grade$n_grade, "| de", round(min(grade$rho), 4),
    "a", round(max(grade$rho), 4), "\n\n")

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

R2_de <- function(est) {
  rss <- est["a"] - 2 * grade$rho * est["b"] + grade$rho2 * est["c"]
  tss <- ta - 2 * grade$rho * tb + grade$rho2 * tc
  pmax(0, 1 - rss / tss)
}

# ------------------------------------ 1. intervalo de R2 realmente visitado
cat("=== Intervalo de R2 visitado pelo SAR ===\n")
set.seed(20200731L)
amostras <- c(lapply(1:K, function(k) sample(K, k)),
              replicate(300, sort(sample(K, sample(0:K, 1))), simplify = FALSE))
r2s <- unlist(lapply(amostras, function(idx) R2_de(estatisticas(idx))))
cat("  modelos amostrados:", length(amostras), "\n")
cat("  R2 minimo:", round(min(r2s), 6), "| maximo:", round(max(r2s), 6), "\n")
cat("  quantis 0,99 e 0,999:", round(quantile(r2s, c(0.99, 0.999)), 6), "\n\n")

Z_MAX <- as.numeric(Sys.getenv("Z_MAX_TAB", min(0.995, max(r2s) + 0.10)))
cat("  limite superior da tabela:", round(Z_MAX, 4), "\n\n")

# ------------------------------------------------- 2. serie exata (referencia)
log_f21_serie <- function(z, k) {
  cc <- (f21a + k) / 2
  vapply(z, function(zi) {
    j_pico <- max(0, ceiling((zi * a_par - cc) / max(1 - zi, 1e-12)))
    J <- min(J_MAX, max(300L, as.integer(j_pico * 4)))
    lt <- c(0, log(zi) * seq_len(J) + L[[k + 1]][seq_len(J)])
    m <- max(lt)
    m + log(sum(exp(lt - m)))
  }, numeric(1))
}

# ---------------------------------------------------------- 3. tabela
N_TAB <- as.integer(Sys.getenv("N_TAB", "80001"))
z_tab <- seq(0, Z_MAX, length.out = N_TAB)
cat("=== Construindo a tabela ===\n")
cat("  ", N_TAB, "pontos de z x", K + 1, "valores de k\n")
t0 <- proc.time()[3]
TAB <- vapply(0:K, function(k) log_f21_serie(z_tab, k), numeric(N_TAB))
cat("  construida em", round(proc.time()[3] - t0, 1), "s\n\n")

passo <- z_tab[2] - z_tab[1]
log_f21_tab <- function(z, k) {
  pos <- z / passo
  i <- as.integer(pos) + 1L
  i[i >= N_TAB] <- N_TAB - 1L
  frac <- pos - (i - 1L)
  col <- TAB[, k + 1L]
  col[i] * (1 - frac) + col[i + 1L] * frac
}

# ---------------------------------------------- 4. erro de interpolacao
cat("=== Erro da interpolacao contra a serie exata ===\n")
set.seed(1)
pior <- 0
for (k in c(0L, 1L, 3L, 7L, 15L, 22L, 30L)) {
  z_teste <- runif(400, 0, Z_MAX)
  e <- max(abs(log_f21_tab(z_teste, k) - log_f21_serie(z_teste, k)))
  pior <- max(pior, e)
  cat(sprintf("  k = %2d  maior |erro| em log: %.3e\n", k, e))
}
cat("  pior caso:", format(pior, scientific = TRUE), "\n")
if (pior > 1e-7) stop("Erro de interpolacao alto demais. Aumentar N_TAB.")
cat("  veredito: APROVADO\n\n")

# Efeito do erro na log-marginal integrada, que e o que importa de fato.
cat("=== Efeito na log-marginal integrada ===\n")
lm_com <- function(est, f21fun) {
  k <- unname(est["k"])
  rss <- est["a"] - 2 * grade$rho * est["b"] + grade$rho2 * est["c"]
  tss <- ta - 2 * grade$rho * tb + grade$rho2 * tc
  R2 <- pmax(0, 1 - rss / tss)
  lf <- grade$ld_mais_prior + (1 - n) / 2 * log(tss)
  if (k > 0L) lf <- lf + log(f21a - 2) - log(k + f21a - 2) + f21fun(R2, k)
  m <- max(lf)
  m + log(grade$h / 3 * sum(grade$pesos * exp(lf - m)))
}
set.seed(2)
alvo <- replicate(60, sort(sample(K, sample(0:K, 1))), simplify = FALSE)
d <- vapply(alvo, function(idx) {
  e <- estatisticas(idx)
  unname(lm_com(e, log_f21_tab) - lm_com(e, function(z, k) log_f21_serie(z, k)))
}, numeric(1))
cat("  modelos testados:", length(alvo), "\n")
cat("  maior |diferenca| na log-marginal:", format(max(abs(d)), scientific = TRUE), "\n")
if (max(abs(d)) > 1e-6) stop("Diferenca na log-marginal alta demais.")
cat("  veredito: APROVADO\n\n")

# ------------------------------------------------------------ 5. ganho
cat("=== Velocidade ===\n")
set.seed(3)
idxs <- replicate(150L, sort(sample(K, sample(1:20, 1))), simplify = FALSE)
ests <- lapply(idxs, estatisticas)
t0 <- proc.time()[3]; for (e in ests) invisible(lm_com(e, log_f21_tab))
t_tab <- (proc.time()[3] - t0) / 150
t0 <- proc.time()[3]
for (e in ests) invisible(lm_com(e, function(z, k) log_f21_serie(z, k)))
t_ser <- (proc.time()[3] - t0) / 150
t0 <- proc.time()[3]; for (idx in idxs) invisible(estatisticas(idx))
t_est <- (proc.time()[3] - t0) / 150

cat("  serie   :", round(t_ser * 1000, 3), "ms por modelo\n")
cat("  tabela  :", round(t_tab * 1000, 3), "ms por modelo\n")
cat("  QR      :", round(t_est * 1000, 3), "ms por modelo\n")
cat("  ganho na parte da hipergeometrica:", round(t_ser / t_tab, 1), "vezes\n")
cat("  projecao para 1.000.000 x 2 cadeias:",
    round((t_est + t_tab) * 2e6 / 60, 1), "min\n\n")

saveRDS(list(TAB = TAB, z_tab = z_tab, passo = passo, N_TAB = N_TAB, Z_MAX = Z_MAX),
        file.path(AQUI, "saida", Sys.getenv("ARQ_TAB", "tabela_f21.rds")))
cat("Tabela gravada em saida/tabela_f21.rds\n")
