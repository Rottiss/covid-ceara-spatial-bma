# Passo 1, fase de desenho: escolher os cenarios e calibrar o sinal antes de
# rodar o experimento.
#
# PERGUNTA. Quando as PIP colapsam sob o SAR, isso significa que as covariaveis
# nao tinham conteudo proprio, ou que o termo espacial destroi PIP mesmo de
# covariaveis genuinamente causais? Os dados reais nao distinguem.
#
# DESENHO. Mantem X real, com toda a estrutura espacial que ele tem, e gera y:
#     y = X_verdadeiro * beta + epsilon,   com rho = 0 por construcao.
# Nao ha dependencia espacial no desfecho alem da que as proprias covariaveis
# carregam. Se o SAR ainda assim derrubar as PIP das verdadeiras, o metodo nao
# distingue.
#
# DOIS CENARIOS, em vez de um. A hipotese e que a sobrevivencia da PIP sob o SAR
# dependa de quao espacialmente suave e a covariavel verdadeira. Testar isso
# exige contraste:
#   - "suaves": as quatro continuas de maior I de Moran
#   - "rugosas": as quatro de menor I de Moran
# Se as suaves perderem PIP e as rugosas nao, o mecanismo esta identificado.
#
# CALIBRACAO. O R2 alvo e o do MODELO VERDADEIRO, nao o do modelo cheio. A
# primeira versao deste script media o R2 da regressao sobre as 30 covariaveis,
# que infla por sobreajuste (30 regressores para 184 observacoes) e deu 0,54
# quando o alvo era 0,35.

source(if (file.exists("../caminhos.R")) "../caminhos.R" else "caminhos.R")

suppressPackageStartupMessages(library(spdep))

AQUI <- ESPACIAL
base <- readRDS(file.path(AQUI, "saida", "base_sar.rds"))
X <- base$X; n <- base$n; covariaveis <- base$covariaveis
listw <- base$listw; K <- length(covariaveis)
y_real <- base$y
Xc <- sweep(X, 2, colMeans(X), "-")

R2_de <- function(y, indices) {
  yc <- y - mean(y)
  if (length(indices) == 0L) return(0)
  Q <- qr.Q(qr(Xc[, indices, drop = FALSE]))
  sum((crossprod(Q, yc))^2) / sum(yc^2)
}

# --------------------------------- 1. autocorrelacao de cada covariavel
moran <- vapply(seq_len(K), function(j) {
  spdep::moran.test(X[, j], listw, randomisation = TRUE)$estimate[["Moran I statistic"]]
}, numeric(1))
names(moran) <- covariaveis
mi_y <- spdep::moran.test(y_real, listw)$estimate[["Moran I statistic"]]

cat("=== Referencias dos dados reais ===\n")
cat("  I de Moran do desfecho (casos):", round(mi_y, 4), "\n")
cat("  I de Moran das covariaveis: de", round(min(moran), 4), "a",
    round(max(moran), 4), ", mediana", round(median(moran), 4), "\n")
cat("  R2 do modelo verdadeiro-equivalente (as 6 com PIP > 0,5):",
    round(R2_de(y_real, match(c("z_idosos", "z_leitos1k", "metrop", "semiarido",
                                "z_ideb.fund", "z_int.circ"), covariaveis)), 4), "\n\n")

# ----------------------------------------- 2. os dois cenarios
continuas <- covariaveis[grepl("^z_", covariaveis)]
mc <- sort(moran[continuas], decreasing = TRUE)

cenarios <- list(
  suaves  = names(mc)[1:4],
  rugosas = names(mc)[(length(mc) - 3):length(mc)]
)

cat("=== Cenarios ===\n")
for (nome in names(cenarios)) {
  vs <- cenarios[[nome]]
  cat(" ", nome, ":\n", sep = "")
  for (v in vs) cat(sprintf("    %-16s I de Moran = %7.4f\n", v, moran[v]))
  cat(sprintf("    I de Moran medio: %.4f\n", mean(moran[vs])))
}
cat("\n  As outras", K - 4, "covariaveis entram no BMA como candidatas com beta zero.\n\n")

# --------------------------------------------- 3. calibrar o sinal
R2_ALVO <- 0.35
cat("=== Calibracao, R2 alvo do modelo verdadeiro =", R2_ALVO, "===\n")

calibrar <- function(vs) {
  idx <- match(vs, covariaveis)
  Xa <- Xc[, idx, drop = FALSE]
  v <- var(rowSums(Xa))
  b <- sqrt(R2_ALVO / (1 - R2_ALVO) / v)
  list(idx = idx, beta = rep(b, length(idx)), b = b)
}

set.seed(20200731L)
config <- list()
for (nome in names(cenarios)) {
  cal <- calibrar(cenarios[[nome]])
  # Confere numa amostra de realizacoes, sobre o modelo VERDADEIRO.
  r2s <- replicate(200, {
    yy <- as.numeric(Xc[, cal$idx, drop = FALSE] %*% cal$beta) + rnorm(n)
    R2_de(yy, cal$idx)
  })
  mis <- replicate(50, {
    yy <- as.numeric(Xc[, cal$idx, drop = FALSE] %*% cal$beta) + rnorm(n)
    spdep::moran.test(yy, listw)$estimate[["Moran I statistic"]]
  })
  cat(sprintf("  %-8s beta = %.4f | R2 medio = %.4f (dp %.4f) | I de Moran de y = %.4f\n",
              nome, cal$b, mean(r2s), sd(r2s), mean(mis)))
  config[[nome]] <- c(cal, list(vars = cenarios[[nome]],
                                r2_medio = mean(r2s), mi_medio = mean(mis)))
}

cat("\n  Para comparar, o desfecho real tem I de Moran de", round(mi_y, 4), ".\n")
cat("  O cenario 'suaves' se aproxima dele sem qualquer termo espacial no DGP;\n")
cat("  o 'rugosas' fica proximo de zero. E esse o contraste que interessa.\n\n")

saveRDS(list(config = config, moran = moran, mi_y_real = mi_y, R2_alvo = R2_ALVO),
        file.path(AQUI, "saida", "desenho_simulacao.rds"))
cat("Desenho gravado em saida/desenho_simulacao.rds\n")
