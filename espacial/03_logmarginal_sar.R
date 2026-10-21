# Etapa 2 do roteiro LeSage e Parent (2007): log-marginal do SAR com
# integracao numerica em rho.
#
# Implementa a Equacao (11) do artigo, com a vetorizacao do Apendice A. Faz tres
# verificacoes independentes:
#
#   A. com rho = 0, a expressao colapsa exatamente na log-marginal validada na
#      Etapa 1 (diferenca esperada: precisao de maquina);
#   B. com g -> infinito, o modo do posterior de rho coincide com a estimativa
#      de maxima verossimilhanca do SAR obtida por spatialreg::lagsarlm;
#   C. a integracao por Simpson converge: refinar a grade nao muda o resultado.
#
# CONVENCAO DE g. O artigo escreve a priori como N[b0, sigma^2 (g X'X)^-1], de
# modo que o g dele multiplica X'X e vale 1/max{n, k^2}. Aqui uso a convencao
# usual, em que a priori e N[b0, sigma^2 g (X'X)^-1] e g = n. As duas se
# correspondem por g_artigo = 1/g_usual, o que troca os papeis dos fatores
# 1/(g+1) e g/(g+1) na Equacao (11). Todas as formulas abaixo estao na
# convencao usual.

source(if (file.exists("../caminhos.R")) "../caminhos.R" else "caminhos.R")

suppressPackageStartupMessages({
  library(BMS)
  library(spdep)
  library(spatialreg)
})

AQUI <- ESPACIAL
RAIZ <- ESTIMACOES
MATRIZ <- file.path(RAIZ, "resultados/estimacoes/matrizes_estimacao/principal/casos.csv")
PESOS <- file.path(RAIZ, "resultados/diagnostico_espacial/pesos_queen.rds")
SAIDA <- file.path(AQUI, "saida")
dir.create(SAIDA, recursive = TRUE, showWarnings = FALSE)

# ---------------------------------------------------------------- 1. dados
dados <- read.csv(MATRIZ, stringsAsFactors = FALSE, colClasses = c(codigo_ibge = "character"))
covariaveis <- setdiff(names(dados), c("codigo_ibge", "municipio", "z_log_casos100k"))
y <- as.numeric(dados$z_log_casos100k)
X_total <- as.matrix(dados[, covariaveis])
n <- length(y)

peso_queen <- readRDS(PESOS)
W <- spdep::listw2mat(peso_queen)
dimnames(W) <- NULL

stopifnot(nrow(W) == n, ncol(W) == n,
          max(abs(rowSums(W) - 1)) < 1e-9,
          max(abs(diag(W))) < 1e-12)

cat("=== Dados ===\n")
cat("  n =", n, "| covariaveis disponiveis =", length(covariaveis), "\n")
cat("  W: Queen, normalizada por linha, de diagnostico_espacial.R\n\n")

# ------------------------------------------------- 2. autovalores de W
# O termo caro da Equacao (11) e ln|I - rho W|. Com n = 184 ele sai exato a
# partir dos autovalores, calculados uma unica vez:
#     ln|I - rho W| = sum_i ln(1 - rho * mu_i)
# Isso dispensa os algoritmos esparsos do Apendice A, pensados para n na casa
# dos milhares.
t0 <- proc.time()[3]
mu <- Re(eigen(W, only.values = TRUE)$values)
tempo_autovalores <- proc.time()[3] - t0

mu_min <- min(mu); mu_max <- max(mu)
cat("=== Autovalores de W ===\n")
cat("  calculados em", round(tempo_autovalores, 3), "s\n")
cat("  minimo =", round(mu_min, 6), "| maximo =", round(mu_max, 6), "\n")
cat("  suporte teorico de rho: [", round(1 / mu_min, 4), ",", round(1 / mu_max, 4), "]\n\n")

log_det <- function(rho) sum(log1p(-rho * mu))

# ------------------------------------------- 3. log-marginal do SAR
# Vetorizacao do Apendice A, equacoes (A2) e (A3):
#   e(rho) = e_o - rho * e_d, com
#   e_o = y - X b_o - alpha_o,  b_o = (X'X)^-1 X'y,      alpha_o = media(y)
#   e_d = Wy - X b_d - alpha_d, b_d = (X'X)^-1 X'Wy,     alpha_d = media(Wy)
# Assim RSS(rho) = ||e(rho)||^2 sai de tres produtos internos pre-computados,
# sem refazer minimos quadrados a cada valor de rho.

preparar_modelo <- function(indices, g) {
  Wy <- as.numeric(W %*% y)
  yc <- y - mean(y)
  Wyc <- Wy - mean(Wy)

  if (length(indices) == 0L) {
    eo <- yc; ed <- Wyc
  } else {
    Xj <- X_total[, indices, drop = FALSE]
    Xc <- sweep(Xj, 2, colMeans(Xj), "-")
    Q <- qr.Q(qr(Xc))
    eo <- yc - Q %*% crossprod(Q, yc)
    ed <- Wyc - Q %*% crossprod(Q, Wyc)
  }
  eo <- as.numeric(eo); ed <- as.numeric(ed)

  list(
    k = length(indices),
    # RSS(rho) = a - 2*rho*b + rho^2*c
    a = sum(eo^2), b = sum(eo * ed), c = sum(ed^2),
    # TSS(rho) = ||(I - rho W)y - media||^2 = ta - 2*rho*tb + rho^2*tc
    ta = sum(yc^2), tb = sum(yc * Wyc), tc = sum(Wyc^2),
    g = g
  )
}

# log da densidade nao normalizada p(rho | y, X, W, M), a menos de constante
# comum a todos os modelos.
log_dens_rho <- function(mod, rho, a0 = 1.01) {
  g <- mod$g
  rss <- mod$a - 2 * rho * mod$b + rho^2 * mod$c
  tss <- mod$ta - 2 * rho * mod$tb + rho^2 * mod$tc
  sq <- (g / (1 + g)) * rss + (1 / (1 + g)) * tss
  ld <- vapply(rho, log_det, numeric(1))
  # priori Beta(a0, a0) reescalada para [-1, 1], como na Equacao (9)
  log_prior <- (a0 - 1) * log1p(rho) + (a0 - 1) * log1p(-rho) -
    (2 * a0 - 1) * log(2) - lbeta(a0, a0)
  -(mod$k / 2) * log(1 + g) + ld - ((n - 1) / 2) * log(sq) + log_prior
}

# Integracao por regra de Simpson sobre uma grade de rho, com escalonamento pelo
# maximo antes de exponenciar (o artigo recomenda isso na pagina 266).
log_marginal_sar <- function(mod, n_grade = 1001L, limites = c(-0.999, 0.999), a0 = 1.01) {
  if (n_grade %% 2L == 0L) n_grade <- n_grade + 1L
  rho <- seq(limites[1], limites[2], length.out = n_grade)
  h <- rho[2] - rho[1]
  lf <- log_dens_rho(mod, rho, a0)
  m <- max(lf)
  f <- exp(lf - m)
  pesos <- c(1, rep(c(4, 2), length.out = n_grade - 2), 1)
  m + log(h / 3 * sum(pesos * f))
}

# --------------------------------- 4. verificacao A: rho = 0 reproduz a Etapa 1
cat("=== Verificacao A: com rho = 0, colapsa na log-marginal da Etapa 1 ===\n")

g_uip <- n
log_marginal_ols <- function(indices, g) {
  yc <- y - mean(y)
  if (length(indices) == 0L) {
    quad <- 0
  } else {
    Xj <- X_total[, indices, drop = FALSE]
    Xc <- sweep(Xj, 2, colMeans(Xj), "-")
    Q <- qr.Q(qr(Xc))
    quad <- sum((crossprod(Q, yc))^2)
  }
  -(length(indices) / 2) * log(1 + g) -
    ((n - 1) / 2) * log(sum(yc^2) - (g / (1 + g)) * quad)
}

set.seed(20200731L)
amostra <- replicate(20, sort(sample(seq_along(covariaveis), sample(0:6, 1))), simplify = FALSE)
difs <- vapply(amostra, function(idx) {
  mod <- preparar_modelo(idx, g_uip)
  # a densidade em rho = 0, sem o termo de priori, deve ser a log-marginal OLS
  parte_sar <- log_dens_rho(mod, 0) -
    ((1.01 - 1) * 0 + (1.01 - 1) * 0 - (2 * 1.01 - 1) * log(2) - lbeta(1.01, 1.01))
  parte_sar - log_marginal_ols(idx, g_uip)
}, numeric(1))

cat("  modelos testados:", length(amostra), "\n")
cat("  maior diferenca absoluta:", format(max(abs(difs)), scientific = TRUE), "\n")
cat("  veredito:", if (max(abs(difs)) < 1e-10) "APROVADO" else "REPROVADO", "\n\n")

# --------------------- 5. verificacao B: modo em rho contra maxima verossimilhanca
cat("=== Verificacao B: modo do posterior de rho vs. MV do SAR ===\n")
cat("  Com g grande o termo de encolhimento some e a expressao fica proporcional\n")
cat("  a log|I - rho W| - ((n-1)/2) log RSS(rho), que difere da log-verossimilhanca\n")
cat("  concentrada do SAR apenas por (n-1)/2 contra n/2. Os modos devem quase coincidir.\n\n")

indices_todos <- seq_along(covariaveis)
modo_de <- function(g) {
  mod <- preparar_modelo(indices_todos, g)
  opt <- optimize(function(r) log_dens_rho(mod, r, a0 = 1), c(-0.99, 0.99), maximum = TRUE)
  opt$maximum
}

formula_sar <- as.formula(paste("y ~", paste(covariaveis, collapse = " + ")))
quadro <- data.frame(y = y, X_total, check.names = FALSE)
ajuste_mv <- spatialreg::lagsarlm(formula_sar, data = quadro, listw = peso_queen)
rho_mv <- unname(ajuste_mv$rho)

modos <- data.frame(
  g = c(g_uip, 1e4, 1e8),
  modo_rho = vapply(c(g_uip, 1e4, 1e8), modo_de, numeric(1))
)
modos$dif_para_mv <- modos$modo_rho - rho_mv

cat("  rho por maxima verossimilhanca (lagsarlm, K = 30):", round(rho_mv, 6), "\n\n")
print(transform(modos, modo_rho = round(modo_rho, 6), dif_para_mv = signif(dif_para_mv, 3)),
      row.names = FALSE)

dif_limite <- abs(modos$dif_para_mv[nrow(modos)])
cat("\n  diferenca com g = 1e8:", signif(dif_limite, 4), "\n")
cat("  veredito:", if (dif_limite < 5e-3) "APROVADO" else "INVESTIGAR", "\n\n")

# ------------------------------- 6. verificacao C: convergencia da integracao
cat("=== Verificacao C: convergencia da regra de Simpson ===\n")
mod_teste <- preparar_modelo(indices_todos, g_uip)
grades <- c(101L, 201L, 401L, 801L, 1601L, 3201L)
valores <- vapply(grades, function(ng) log_marginal_sar(mod_teste, ng), numeric(1))
conv <- data.frame(pontos = grades, log_marginal = valores)
conv$dif_para_mais_fina <- valores - valores[length(valores)]
print(transform(conv, log_marginal = round(log_marginal, 10),
                dif_para_mais_fina = signif(dif_para_mais_fina, 3)), row.names = FALSE)

dif_1001 <- abs(log_marginal_sar(mod_teste, 1001L) - valores[length(valores)])
cat("\n  diferenca entre 1.001 pontos e 3.201 pontos:", signif(dif_1001, 3), "\n")
cat("  veredito:", if (dif_1001 < 1e-6) "APROVADO" else "INVESTIGAR", "\n\n")

# ------------------------------------------------------- 7. custo por modelo
cat("=== Custo computacional ===\n")
t0 <- proc.time()[3]
for (i in 1:200) {
  idx <- sample(indices_todos, sample(1:10, 1))
  invisible(log_marginal_sar(preparar_modelo(idx, g_uip), 1001L))
}
por_modelo <- (proc.time()[3] - t0) / 200
cat("  tempo medio por modelo (grade de 1.001 pontos):",
    round(por_modelo * 1000, 2), "ms\n")
cat("  projecao para 250.000 sorteios de MC3:",
    round(por_modelo * 250000 / 60, 1), "min\n\n")

saveRDS(list(mu = mu, W = W, y = y, X = X_total, covariaveis = covariaveis,
             listw = peso_queen, n = n),
        file.path(SAIDA, "base_sar.rds"))
write.csv(conv, file.path(SAIDA, "convergencia_simpson.csv"), row.names = FALSE)

cat("Gravado em", SAIDA, "\n")
cat("R:", R.version.string, "| spatialreg:",
    as.character(packageVersion("spatialreg")), "\n")
