# Etapa 1 do roteiro LeSage e Parent (2007): validar a log-marginal.
#
# Antes de acrescentar o termo espacial, e preciso ter certeza de que a
# log-marginal implementada aqui reproduz a do BMS no caso nao espacial. Com
# rho = 0 a Equacao (11) do artigo colapsa no caso de minimos quadrados, entao
# esta e a verificacao que sustenta todo o resto.
#
# Estrategia: subconjunto de 10 covariaveis, 2^10 = 1024 modelos, enumeracao
# completa dos dois lados. Sem MCMC, sem aproximacao, sem semente.
#
# Priori usada na comparacao: g-priori de Zellner com g fixo (UIP, g = N) e
# priori uniforme sobre o espaco de modelos. NAO e a configuracao do TCC
# (hyper-g UIP com Binomial-Beta): hyper-g acrescenta uma integral em g e
# mprior="random" muda o peso dos modelos por tamanho. Ambas seriam ruido nesta
# etapa, cujo objetivo e conferir a algebra da verossimilhanca marginal.
#
# Saida: mc3_espacial/saida/validacao_logmarginal.csv e um resumo no terminal.

source(if (file.exists("../caminhos.R")) "../caminhos.R" else "caminhos.R")

suppressPackageStartupMessages(library(BMS))

AQUI <- ESPACIAL
MATRIZ <- file.path(
  ESTIMACOES,
  "resultados/estimacoes/matrizes_estimacao/principal/casos.csv"
)
SAIDA <- file.path(AQUI, "saida")
dir.create(SAIDA, recursive = TRUE, showWarnings = FALSE)

K_TESTE <- 10L

# ---------------------------------------------------------------- 1. dados
dados <- read.csv(MATRIZ, stringsAsFactors = FALSE, colClasses = c(codigo_ibge = "character"))
ids <- c("codigo_ibge", "municipio")
desfecho <- "z_log_casos100k"
covariaveis <- setdiff(names(dados), c(ids, desfecho))

stopifnot(nrow(dados) == 184L, length(covariaveis) == 30L)

y <- as.numeric(dados[[desfecho]])
X_total <- as.matrix(dados[, covariaveis])
X <- X_total[, seq_len(K_TESTE), drop = FALSE]

n <- length(y)
K <- ncol(X)
g <- n  # UIP

cat("=== Configuracao ===\n")
cat("  n =", n, "| K de teste =", K, "| modelos =", 2^K, "\n")
cat("  g =", g, "(UIP, fixo) | priori de modelo: uniforme\n")
cat("  covariaveis:", paste(colnames(X), collapse = ", "), "\n\n")

# ------------------------------------------------- 2. log-marginal propria
# Com y e X centrados, priori nao informativa no intercepto e em sigma^2, e
# g-priori de Zellner, a verossimilhanca marginal do modelo M_j e proporcional a
#
#   (1 + g)^(-k_j/2) * [ y'y - (g/(1+g)) y'P_j y ]^(-(n-1)/2)
#
# com P_j o projetor sobre as colunas de X_j centradas. O modelo nulo tem
# k_j = 0 e P_j = 0.

centrar <- function(m) sweep(m, 2, colMeans(m), "-")

yc <- y - mean(y)
Xc <- centrar(X)
yty <- sum(yc^2)

log_marginal <- function(indices) {
  k <- length(indices)
  if (k == 0L) {
    quad <- 0
  } else {
    Xj <- Xc[, indices, drop = FALSE]
    # Via decomposicao QR: y'P_j y = ||Q'y||^2, estavel e sem inverter X'X.
    Q <- qr.Q(qr(Xj))
    quad <- sum((crossprod(Q, yc))^2)
  }
  sr <- g / (1 + g)
  -(k / 2) * log(1 + g) - ((n - 1) / 2) * log(yty - sr * quad)
}

# Enumeracao completa dos 2^K modelos.
combinacoes <- lapply(0:(2^K - 1), function(i) which(bitwAnd(i, 2^(0:(K - 1))) > 0))
logml <- vapply(combinacoes, log_marginal, numeric(1))

# Probabilidades posteriores com priori uniforme sobre modelos.
m <- max(logml)
pmp_propria <- exp(logml - m) / sum(exp(logml - m))

# PIP de cada covariavel: soma das PMP dos modelos que a contem.
pip_propria <- vapply(seq_len(K), function(j) {
  sum(pmp_propria[vapply(combinacoes, function(s) j %in% s, logical(1))])
}, numeric(1))

# --------------------------------------------------------------- 3. BMS
cat("=== Rodando BMS com enumeracao completa ===\n")
dados_bms <- cbind(y, X)
ajuste <- BMS::bms(
  dados_bms,
  mcmc = "enumerate",
  g = "UIP",
  mprior = "uniform",
  nmodel = 2^K_TESTE,
  user.int = FALSE
)

coefs <- coef(ajuste, order.by.pip = FALSE, exact = TRUE)
pip_bms <- coefs[colnames(X), "PIP"]

# ------------------------------------------------------------ 4. comparacao
comparacao <- data.frame(
  covariavel = colnames(X),
  pip_propria = pip_propria,
  pip_bms = as.numeric(pip_bms),
  stringsAsFactors = FALSE
)
comparacao$diferenca <- comparacao$pip_propria - comparacao$pip_bms

cat("\n=== PIP: implementacao propria vs. BMS ===\n")
print(
  transform(
    comparacao,
    pip_propria = round(pip_propria, 8),
    pip_bms = round(pip_bms, 8),
    diferenca = signif(diferenca, 3)
  ),
  row.names = FALSE
)

max_dif <- max(abs(comparacao$diferenca))
cat("\n  maior diferenca absoluta em PIP:", format(max_dif, scientific = TRUE), "\n")

# Comparacao tambem no nivel das probabilidades de modelo, que e mais exigente.
tabela_bms <- BMS::pmp.bma(ajuste, oldstyle = TRUE)
pmp_bms_top <- head(tabela_bms[order(-tabela_bms[, 1]), 1, drop = FALSE], 10)

ordem <- order(-pmp_propria)
top_propria <- data.frame(
  posicao = 1:10,
  variaveis = vapply(combinacoes[ordem][1:10], function(s) {
    if (length(s) == 0L) "(modelo nulo)" else paste(colnames(X)[s], collapse = " + ")
  }, character(1)),
  pmp_propria = pmp_propria[ordem][1:10],
  stringsAsFactors = FALSE
)

cat("\n=== Dez modelos de maior probabilidade posterior (implementacao propria) ===\n")
print(transform(top_propria, pmp_propria = round(pmp_propria, 6)), row.names = FALSE)

cat("\n=== Dez maiores PMP segundo o BMS ===\n")
print(round(pmp_bms_top, 6))

write.csv(comparacao, file.path(SAIDA, "validacao_logmarginal.csv"), row.names = FALSE)

cat("\n=== Veredito ===\n")
if (max_dif < 1e-8) {
  cat("  APROVADO: a log-marginal reproduz o BMS ate 1e-8 em todas as PIP.\n")
  cat("  A Etapa 2 (integracao em rho) pode partir desta base.\n")
} else if (max_dif < 1e-4) {
  cat("  ATENCAO: diferenca de", format(max_dif, scientific = TRUE),
      "- pequena, mas maior que erro de maquina. Investigar antes de seguir.\n")
} else {
  cat("  REPROVADO: diferenca de", format(max_dif, scientific = TRUE),
      "- a implementacao nao reproduz o BMS. Nao seguir para a Etapa 2.\n")
}

cat("\nR:", R.version.string, "| BMS:", as.character(packageVersion("BMS")), "\n")
