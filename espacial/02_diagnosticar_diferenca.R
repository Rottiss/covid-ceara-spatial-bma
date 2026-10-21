# Diagnostica a diferenca sistematica entre a log-marginal propria e a do BMS.
#
# A comparacao anterior deu diferenca de ~9e-4 nas PIP, com o mesmo ordenamento
# de modelos. Isso indica divergencia de convencao, nao erro numerico. Este
# script isola onde ela esta: extrai as PMP exatas do BMS modelo a modelo,
# converte em log-marginal implicita e verifica se a diferenca depende apenas do
# numero de covariaveis k.

source(if (file.exists("../caminhos.R")) "../caminhos.R" else "caminhos.R")

suppressPackageStartupMessages(library(BMS))

MATRIZ <- file.path(
  ESTIMACOES,
  "resultados/estimacoes/matrizes_estimacao/principal/casos.csv"
)

dados <- read.csv(MATRIZ, stringsAsFactors = FALSE, colClasses = c(codigo_ibge = "character"))
covariaveis <- setdiff(names(dados), c("codigo_ibge", "municipio", "z_log_casos100k"))
y <- as.numeric(dados$z_log_casos100k)
X <- as.matrix(dados[, covariaveis])[, 1:10, drop = FALSE]
n <- length(y); K <- ncol(X); g <- n

ajuste <- BMS::bms(cbind(y, X), mcmc = "enumerate", g = "UIP",
                   mprior = "uniform", user.int = FALSE)

# topmodels.bma devolve a matriz de inclusao (K linhas) mais PMP nas ultimas
# linhas. Pega todos os 2^K modelos.
tm <- BMS::topmodels.bma(ajuste)
inclusao <- tm[seq_len(K), , drop = FALSE]
pmp_bms <- tm[K + 1, ]  # PMP exata

cat("modelos devolvidos pelo BMS:", ncol(tm), "de", 2^K, "\n")
cat("soma das PMP do BMS:", sum(pmp_bms), "\n\n")

yc <- y - mean(y)
Xc <- sweep(X, 2, colMeans(X), "-")
yty <- sum(yc^2)

quad_de <- function(idx) {
  if (length(idx) == 0L) return(0)
  Q <- qr.Q(qr(Xc[, idx, drop = FALSE]))
  sum((crossprod(Q, yc))^2)
}

k_vec <- integer(ncol(tm))
logml_propria <- numeric(ncol(tm))
for (i in seq_len(ncol(tm))) {
  idx <- which(inclusao[, i] == 1)
  k_vec[i] <- length(idx)
  logml_propria[i] <- -(k_vec[i] / 2) * log(1 + g) -
    ((n - 1) / 2) * log(yty - (g / (1 + g)) * quad_de(idx))
}

# Log-marginal implicita do BMS, a menos de constante: log(PMP) com priori
# uniforme. Ancorada no modelo nulo (ou no primeiro modelo, se o nulo nao
# aparecer) para eliminar a constante.
logml_bms <- log(pmp_bms)
ancora <- which.max(pmp_bms)
d_propria <- logml_propria - logml_propria[ancora]
d_bms <- logml_bms - logml_bms[ancora]
delta <- d_propria - d_bms

cat("=== Diferenca nas log-marginais relativas, por numero de covariaveis ===\n")
tab <- aggregate(delta, by = list(k = k_vec), function(z) {
  c(n = length(z), media = mean(z), dp = sd(z), min = min(z), max = max(z))
})
resumo <- data.frame(k = tab$k, n_modelos = tab$x[, "n"],
                     media = tab$x[, "media"], desvio = tab$x[, "dp"],
                     minimo = tab$x[, "min"], maximo = tab$x[, "max"])
print(transform(resumo, media = signif(media, 6), desvio = signif(desvio, 3),
                minimo = signif(minimo, 6), maximo = signif(maximo, 6)),
      row.names = FALSE)

cat("\nSe 'desvio' for ~0 dentro de cada k, a divergencia e so no termo que\n")
cat("depende de k, isto e, na penalizacao (1+g)^(-k/2).\n\n")

# Ajusta delta ~ k para descobrir a constante por covariavel.
k_ancora <- k_vec[ancora]
fit <- lm(delta ~ I(k_vec - k_ancora) + 0)
inclinacao <- unname(coef(fit)[1])
cat("=== Inclinacao de delta em relacao a k ===\n")
cat("  estimada          :", signif(inclinacao, 10), "\n")
cat("  R^2               :", signif(summary(fit)$r.squared, 8), "\n")
cat("  -0,5*log(1+g)     :", signif(-0.5 * log(1 + g), 10), "  (o que uso)\n")
cat("  -0,5*log(g)       :", signif(-0.5 * log(g), 10), "\n")
cat("  diferenca entre os dois:", signif(-0.5 * log(1 + g) + 0.5 * log(g), 10), "\n\n")

g_implicito <- exp(-2 * inclinacao) - 1
cat("  g implicito no BMS, se a forma for (1+g)^(-k/2):", signif(g_implicito, 10), "\n")
cat("  n =", n, "\n")

# Teste direto: refaz a comparacao trocando (1+g)^(-k/2) por g^(-k/2)... e
# tambem variando o expoente dos graus de liberdade.
cat("\n=== Teste de variantes ===\n")
variantes <- list(
  "(1+g)^(-k/2), df = (n-1)/2" = function(k, q) -(k/2)*log(1+g) - ((n-1)/2)*log(yty - (g/(1+g))*q),
  "(1+g)^(-k/2), df = (n-1-k)/2" = function(k, q) -(k/2)*log(1+g) - ((n-1-k)/2)*log(yty - (g/(1+g))*q),
  "g^(-k/2),     df = (n-1)/2" = function(k, q) -(k/2)*log(g) - ((n-1)/2)*log(yty - (g/(1+g))*q)
)
quads <- vapply(seq_len(ncol(tm)), function(i) quad_de(which(inclusao[, i] == 1)), numeric(1))

for (nome in names(variantes)) {
  lm_v <- variantes[[nome]](k_vec, quads)
  p <- exp(lm_v - max(lm_v)); p <- p / sum(p)
  pip_v <- vapply(seq_len(K), function(j) sum(p[inclusao[j, ] == 1]), numeric(1))
  pip_b <- vapply(seq_len(K), function(j) sum(pmp_bms[inclusao[j, ] == 1]), numeric(1))
  cat(sprintf("  %-32s maior |dif| em PIP: %.3e\n", nome, max(abs(pip_v - pip_b))))
}
