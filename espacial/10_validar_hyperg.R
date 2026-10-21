# Valida a implementacao da priori do TCC: hyper-g UIP com priori de modelo
# Binomial-Beta. Sem esta validacao, comparar MQO e SAR sob "a priori do TCC"
# seria afirmacao sem lastro.
#
# As definicoes abaixo foram extraidas do proprio codigo do BMS 0.3.5
# (.gprior.hyperg.init, .f21_4hyperg, .mprior.randomt.init), nao de memoria.
#
# HYPER-G UIP. O parametro e a = 2 + 2/N. A verossimilhanca marginal integra g
# analiticamente, resultando numa hipergeometrica de Gauss:
#
#   log p(y|M_k) = (1-N)/2 * log(yty) + log(a-2) - log(k+a-2)
#                  + log 2F1( (N-1)/2, 1; (a+k)/2; R2 )
#
#   com R2 = 1 - RSS/TSS. Para k = 0, so o primeiro termo.
#
# Note que nao ha fator (1+g)^(-k/2): ele foi absorvido na integral em g.
#
# BINOMIAL-BETA. Com mprior.size = K/2, o parametro b da Beta vale
# (K - K/2)/(K/2) = 1, ou seja, w ~ Beta(1,1). A priori de modelo fica
#
#   log P(M) = lgamma(1 + k) + lgamma(1 + K - k)   (a menos de constante)

source(if (file.exists("../caminhos.R")) "../caminhos.R" else "caminhos.R")

suppressPackageStartupMessages(library(BMS))

AQUI <- ESPACIAL
MATRIZ <- file.path(ESTIMACOES,
                    "resultados/estimacoes/matrizes_estimacao/principal/casos.csv")

dados <- read.csv(MATRIZ, stringsAsFactors = FALSE, colClasses = c(codigo_ibge = "character"))
covariaveis <- setdiff(names(dados), c("codigo_ibge", "municipio", "z_log_casos100k"))
y <- as.numeric(dados$z_log_casos100k)
X_total <- as.matrix(dados[, covariaveis])
n <- length(y)

K_TESTE <- 10L
X <- X_total[, seq_len(K_TESTE), drop = FALSE]
K <- ncol(X)

yc <- y - mean(y)
Xc <- sweep(X, 2, colMeans(X), "-")
yty <- sum(yc^2)

f21a <- 2 + 2 / n

# 2F1(a, 1; c; z) por serie, com truncamento adaptativo. Com b = 1 os termos
# viram produto cumulativo de (a+j)/(c+j)*z.
f21 <- function(z, k, a_par = (n - 1) / 2, tol = 1e-15, max_termos = 100000L) {
  cc <- (f21a + k) / 2
  termo <- 1; soma <- 1; j <- 0L
  repeat {
    termo <- termo * z * (a_par + j) / (cc + j)
    soma <- soma + termo
    j <- j + 1L
    if (abs(termo) < tol * abs(soma) || j >= max_termos) break
  }
  list(valor = soma, termos = j)
}

log_marginal_hyperg <- function(indices) {
  k <- length(indices)
  null_lik <- (1 - n) / 2 * log(yty)
  if (k == 0L) return(null_lik)
  Q <- qr.Q(qr(Xc[, indices, drop = FALSE]))
  ymy <- yty - sum((crossprod(Q, yc))^2)
  R2 <- 1 - ymy / yty
  f <- f21(R2, k)
  null_lik + log(f21a - 2) - log(k + f21a - 2) + log(f$valor)
}

log_prior_modelo <- function(k, K) lgamma(1 + k) + lgamma(1 + K - k)

cat("=== Configuracao ===\n")
cat("  n =", n, "| K =", K, "| a = 2 + 2/N =", format(f21a, digits = 10), "\n")
cat("  priori de modelo: Binomial-Beta com mprior.size = K/2, isto e, Beta(1,1)\n\n")

# ----------------------------------------------- conferir a hipergeometrica
cat("=== Conferencia da hipergeometrica contra BMS::f21hyper ===\n")
set.seed(1)
casos <- expand.grid(z = c(0.05, 0.2, 0.4, 0.6), k = c(1L, 5L, 10L))
casos$propria <- mapply(function(z, k) f21(z, k)$valor, casos$z, casos$k)
casos$bms <- mapply(function(z, k) BMS::f21hyper((n - 1) / 2, 1, (f21a + k) / 2, z),
                    casos$z, casos$k)
casos$dif_relativa <- abs(casos$propria / casos$bms - 1)
print(transform(casos, propria = signif(propria, 10), bms = signif(bms, 10),
                dif_relativa = signif(dif_relativa, 3)), row.names = FALSE)
cat("  maior diferenca relativa:", format(max(casos$dif_relativa), scientific = TRUE), "\n\n")

# --------------------------------------------- enumeracao completa
cat("=== Enumeracao completa dos", 2^K, "modelos ===\n")
combinacoes <- lapply(0:(2^K - 1), function(i) which(bitwAnd(i, 2^(0:(K - 1))) > 0))
logml <- vapply(combinacoes, log_marginal_hyperg, numeric(1))
lprior <- vapply(combinacoes, function(s) log_prior_modelo(length(s), K), numeric(1))
lpost <- logml + lprior
m <- max(lpost)
pmp <- exp(lpost - m) / sum(exp(lpost - m))
pip_propria <- vapply(seq_len(K), function(j)
  sum(pmp[vapply(combinacoes, function(s) j %in% s, logical(1))]), numeric(1))

# ------------------------------------------------------------- BMS
cat("Rodando BMS com a configuracao do TCC...\n")
ajuste <- BMS::bms(cbind(y, X), mcmc = "enumerate", g = "hyper=UIP",
                   mprior = "random", mprior.size = K / 2,
                   nmodel = 2^K, user.int = FALSE)
cat("  a usado pelo BMS:", format(ajuste$gprior.info$hyper.parameter, digits = 10), "\n")
cat("  mprior:", ajuste$mprior.info$mp.mode, "| tamanho:", ajuste$mprior.info$mp.msize, "\n\n")

pip_bms <- coef(ajuste, order.by.pip = FALSE, exact = TRUE)[colnames(X), "PIP"]

comparacao <- data.frame(
  covariavel = colnames(X),
  pip_propria = pip_propria,
  pip_bms = as.numeric(pip_bms)
)
comparacao$diferenca <- comparacao$pip_propria - comparacao$pip_bms

cat("=== PIP: implementacao propria vs. BMS, priori do TCC ===\n")
print(transform(comparacao, pip_propria = round(pip_propria, 8),
                pip_bms = round(pip_bms, 8), diferenca = signif(diferenca, 3)),
      row.names = FALSE)

max_dif <- max(abs(comparacao$diferenca))
cat("\n  maior diferenca absoluta em PIP:", format(max_dif, scientific = TRUE), "\n")
cat("  veredito:", if (max_dif < 1e-8) "APROVADO" else "REPROVADO", "\n\n")

# Tamanho medio do modelo, mais sensivel a erro na priori de modelo do que as PIP.
tam_propria <- sum(pmp * vapply(combinacoes, length, numeric(1)))
tam_bms <- ajuste$info$msize
cat("=== Tamanho medio a posteriori (sensivel a priori de modelo) ===\n")
cat("  proprio:", round(tam_propria, 8), "| BMS:", round(tam_bms, 8),
    "| diferenca:", signif(tam_propria - tam_bms, 3), "\n")

write.csv(comparacao, file.path(AQUI, "saida", "validacao_hyperg.csv"), row.names = FALSE)
cat("\nR:", R.version.string, "| BMS:", as.character(packageVersion("BMS")), "\n")
