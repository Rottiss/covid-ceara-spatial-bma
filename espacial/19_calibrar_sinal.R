# Recalibra a forca do sinal do experimento.
#
# PROBLEMA DETECTADO NO TESTE DE FUMACA. Com R2 alvo de 0,35 repartido entre 4
# covariaveis, cada uma carrega cerca de 0,088 de R2 e e detectada com PIP de
# 0,95 a 0,99 tanto pelo MQO quanto pelo SAR. Nesse regime o teste nao informa:
# o sinal e forte demais para que qualquer coisa o derrube.
#
# Nos dados reais, as seis covariaveis selecionadas explicam 0,2259 de R2
# JUNTAS, cerca de 0,038 cada, e tem PIP entre 0,57 e 0,94. O experimento
# precisa reproduzir essa dificuldade, nao uma mais facil.
#
# ESTRATEGIA. Buscar o beta que faz a PIP media das verdadeiras sob MQO cair
# perto de 0,75, que e o meio da faixa observada. So MQO nesta etapa, que e
# barato; o SAR entra depois.

source(if (file.exists("../caminhos.R")) "../caminhos.R" else "caminhos.R")

AQUI <- ESPACIAL
base <- readRDS(file.path(AQUI, "saida", "base_sar.rds"))
hyp  <- readRDS(file.path(AQUI, "saida", "hyperg_aux.rds"))
des  <- readRDS(file.path(AQUI, "saida", "desenho_simulacao.rds"))

X <- base$X; n <- base$n; covariaveis <- base$covariaveis
K <- length(covariaveis)
f21a <- hyp$f21a; a_par <- hyp$a_par; L <- hyp$L; J_MAX <- hyp$J_MAX
Xc <- sweep(X, 2, colMeans(X), "-")

log_f21 <- function(zi, k) {
  cc <- (f21a + k) / 2
  j_pico <- max(0, ceiling((zi * a_par - cc) / max(1 - zi, 1e-12)))
  J <- min(J_MAX, max(300L, as.integer(j_pico * 4)))
  lt <- c(0, log(zi) * seq_len(J) + L[[k + 1]][seq_len(J)])
  m <- max(lt); m + log(sum(exp(lt - m)))
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

mc3_ols <- function(y, n_sorteios, burn, semente) {
  yc <- y - mean(y); ta <- sum(yc^2)
  rss <- function(idx) {
    if (length(idx) == 0L) return(ta)
    Q <- qr.Q(qr(Xc[, idx, drop = FALSE]))
    ta - sum((crossprod(Q, yc))^2)
  }
  marg <- function(idx) {
    k <- length(idx); nl <- (1 - n) / 2 * log(ta)
    if (k == 0L) return(nl)
    nl + log(f21a - 2) - log(k + f21a - 2) + log_f21(1 - rss(idx) / ta, k)
  }
  set.seed(semente)
  atual <- sort(sample(seq_len(K), sample(0:K, 1)))
  post <- marg(atual) + lpm[length(atual) + 1]
  inclusao <- numeric(K); guardados <- 0L
  for (it in seq_len(n_sorteios)) {
    k <- length(atual); tipos <- tipos_disp(k)
    tipo <- tipos[sample.int(length(tipos), 1L)]
    fora <- setdiff(seq_len(K), atual)
    prop <- switch(tipo,
                   birth = sort(c(atual, fora[sample.int(length(fora), 1L)])),
                   death = atual[-sample.int(k, 1L)],
                   move  = sort(c(atual[-sample.int(k, 1L)],
                                  fora[sample.int(length(fora), 1L)])))
    post_p <- marg(prop) + lpm[length(prop) + 1]
    inv <- switch(tipo, birth = "death", death = "birth", move = "move")
    la <- post_p - post + log_q(length(prop), inv) - log_q(k, tipo)
    if (is.finite(la) && log(runif(1)) < la) { atual <- prop; post <- post_p }
    if (it > burn) {
      guardados <- guardados + 1L
      if (length(atual) > 0L) inclusao[atual] <- inclusao[atual] + 1
    }
  }
  inclusao / guardados
}

cat("=== Referencia dos dados reais ===\n")
sel <- c("z_idosos", "z_leitos1k", "metrop", "semiarido", "z_ideb.fund", "z_int.circ")
pub <- read.csv(file.path(ESTIMACOES,
                          "resultados/estimacoes/principal/tabelas_tecnicas/resultado_casos.csv"),
                stringsAsFactors = FALSE)
pips_sel <- setNames(pub$PIP, pub$variavel)[sel]
cat("  PIP das seis selecionadas:", paste(round(pips_sel, 3), collapse = ", "), "\n")
cat("  media:", round(mean(pips_sel), 4), "| faixa:", round(min(pips_sel), 3),
    "a", round(max(pips_sel), 3), "\n")
cat("  R2 conjunto delas: 0,2259, ou cerca de 0,038 por covariavel\n\n")

ALVO_PIP <- 0.75
N_MC <- 40000L
R_CAL <- 4L

cat("=== Busca do beta, alvo de PIP media", ALVO_PIP, "sob MQO ===\n")
cat("  ", R_CAL, "replicacoes por valor, ", N_MC, " sorteios cada\n\n", sep = "")

avaliar <- function(vars, idx, b) {
  beta <- rep(b, length(idx))
  pips <- vapply(seq_len(R_CAL), function(r) {
    set.seed(7000L + r)
    y <- as.numeric(Xc[, idx, drop = FALSE] %*% beta) + rnorm(n)
    mean(mc3_ols(y, N_MC, as.integer(N_MC * 0.2), 900L + r)[idx])
  }, numeric(1))
  mean(pips)
}

resultado <- list()
for (nome in c("suaves", "rugosas")) {
  cfg <- des$config[[nome]]
  cat("### cenario", nome, "\n"); flush.console()
  candidatos <- if (nome == "suaves") c(0.32, 0.38, 0.44, 0.50) else c(0.19, 0.22, 0.25)
  valores <- numeric(length(candidatos))
  for (i in seq_along(candidatos)) {
    valores[i] <- avaliar(cfg$vars, cfg$idx, candidatos[i])
    Xa <- Xc[, cfg$idx, drop = FALSE]
    r2 <- candidatos[i]^2 * var(rowSums(Xa)) / (candidatos[i]^2 * var(rowSums(Xa)) + 1)
    cat(sprintf("  beta = %.3f | R2 do modelo verdadeiro = %.3f | PIP media = %.4f\n",
                candidatos[i], r2, valores[i])); flush.console()
  }
  melhor <- candidatos[which.min(abs(valores - ALVO_PIP))]
  cat("  escolhido:", melhor, "\n\n")
  resultado[[nome]] <- list(beta_escolhido = melhor,
                            candidatos = candidatos, pips = valores)
}

saveRDS(resultado, file.path(AQUI, "saida", "calibracao_sinal.rds"))
cat("Gravado em saida/calibracao_sinal.rds\n")
