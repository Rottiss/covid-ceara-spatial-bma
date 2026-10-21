# MQO e SAR com 1.000.000 de iteracoes e 200.000 de descarte, a configuracao
# exata registrada no inventario_estimacoes.csv do TCC.
#
# Substitui a execucao de 400.000 do ETAPA4B, que carregava o vies de
# convergencia diagnosticado em DIAGNOSTICO_VIES.md.
#
# Usa a tabela de interpolacao de log 2F1 construida e validada em
# 15_otimizar_hyperg.R (erro de 4,3e-9 no log; 2,5e-9 na log-marginal
# integrada). Valores de R2 acima do limite da tabela caem de volta na serie
# exata, e sao contados.

source(if (file.exists("../caminhos.R")) "../caminhos.R" else "caminhos.R")

AQUI <- ESPACIAL
base <- readRDS(file.path(AQUI, "saida", "base_sar.rds"))
aux <- readRDS(file.path(AQUI, "saida", "grade_sar.rds"))
hyp <- readRDS(file.path(AQUI, "saida", "hyperg_aux.rds"))
tb <- readRDS(file.path(AQUI, "saida", "tabela_f21.rds"))

X_total <- base$X; n <- base$n; covariaveis <- base$covariaveis
grade <- aux$grade; ta <- aux$ta; tbb <- aux$tb; tcc <- aux$tc
yc <- aux$yc; Wyc <- aux$Wyc
K <- length(covariaveis)
f21a <- hyp$f21a; a_par <- hyp$a_par; L <- hyp$L; J_MAX <- hyp$J_MAX
TAB <- tb$TAB; passo <- tb$passo; N_TAB <- tb$N_TAB; Z_MAX <- tb$Z_MAX

fora_da_tabela <- 0L

log_f21_serie1 <- function(zi, k) {
  cc <- (f21a + k) / 2
  j_pico <- max(0, ceiling((zi * a_par - cc) / max(1 - zi, 1e-12)))
  J <- min(J_MAX, max(300L, as.integer(j_pico * 4)))
  lt <- c(0, log(zi) * seq_len(J) + L[[k + 1]][seq_len(J)])
  m <- max(lt)
  m + log(sum(exp(lt - m)))
}

log_f21 <- function(z, k) {
  acima <- z > Z_MAX
  if (any(acima)) {
    fora_da_tabela <<- fora_da_tabela + sum(acima)
    out <- numeric(length(z))
    out[!acima] <- log_f21(z[!acima], k)
    out[acima] <- vapply(z[acima], log_f21_serie1, numeric(1), k = k)
    return(out)
  }
  pos <- z / passo
  i <- as.integer(pos) + 1L
  i[i >= N_TAB] <- N_TAB - 1L
  frac <- pos - (i - 1L)
  col <- TAB[, k + 1L]
  col[i] * (1 - frac) + col[i + 1L] * frac
}

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

lm_ols <- function(est) {
  k <- unname(est["k"])
  nl <- (1 - n) / 2 * log(ta)
  if (k == 0L) return(nl)
  unname(nl + log(f21a - 2) - log(k + f21a - 2) + log_f21(1 - est["a"] / ta, k))
}

lm_sar <- function(est) {
  k <- unname(est["k"])
  rss <- est["a"] - 2 * grade$rho * est["b"] + grade$rho2 * est["c"]
  tss <- ta - 2 * grade$rho * tbb + grade$rho2 * tcc
  if (any(rss <= 0) || any(tss <= 0)) return(-Inf)
  R2 <- pmax(0, 1 - rss / tss)
  lf <- grade$ld_mais_prior + (1 - n) / 2 * log(tss)
  if (k > 0L) lf <- lf + log(f21a - 2) - log(k + f21a - 2) + log_f21(R2, k)
  m <- max(lf)
  if (!is.finite(m)) return(-Inf)
  unname(m + log(grade$h / 3 * sum(grade$pesos * exp(lf - m))))
}

lpm <- vapply(0:K, function(k) lgamma(1 + k) + lgamma(1 + K - k), numeric(1))
tipos_disp <- function(k) {
  t <- character(0)
  if (k < K) t <- c(t, "birth")
  if (k > 0L) t <- c(t, "death")
  if (k > 0L && k < K) t <- c(t, "move")
  t
}
log_q <- function(k, tipo) {
  -log(length(tipos_disp(k))) -
    log(switch(tipo, birth = K - k, death = k, move = k * (K - k)))
}

mc3 <- function(log_marginal, n_sorteios, burn, semente, cortes) {
  set.seed(semente)
  atual <- sort(sample(seq_len(K), sample(0:K, 1)))
  post <- log_marginal(estatisticas(atual)) + lpm[length(atual) + 1]
  aceitos <- 0L; inclusao <- numeric(K); guardados <- 0L; soma_k <- 0
  historico <- list()
  for (it in seq_len(n_sorteios)) {
    k <- length(atual)
    tipos <- tipos_disp(k)
    tipo <- tipos[sample.int(length(tipos), 1L)]
    fora <- setdiff(seq_len(K), atual)
    prop <- switch(tipo,
                   birth = sort(c(atual, fora[sample.int(length(fora), 1L)])),
                   death = atual[-sample.int(k, 1L)],
                   move  = sort(c(atual[-sample.int(k, 1L)],
                                  fora[sample.int(length(fora), 1L)])))
    post_p <- log_marginal(estatisticas(prop)) + lpm[length(prop) + 1]
    inv <- switch(tipo, birth = "death", death = "birth", move = "move")
    la <- post_p - post + log_q(length(prop), inv) - log_q(k, tipo)
    if (is.finite(la) && log(runif(1)) < la) {
      atual <- prop; post <- post_p; aceitos <- aceitos + 1L
    }
    if (it > burn) {
      guardados <- guardados + 1L
      if (length(atual) > 0L) inclusao[atual] <- inclusao[atual] + 1
      soma_k <- soma_k + length(atual)
      if (it %in% cortes) {
        historico[[as.character(it)]] <- soma_k / guardados
      }
    }
  }
  list(pip = inclusao / guardados, taxa = aceitos / n_sorteios,
       tamanho_medio = soma_k / guardados, historico = historico)
}

N <- 1000000L; BURN <- 200000L
cortes <- c(400000L, 600000L, 800000L, 1000000L)

cat("=== Configuracao ===\n")
cat("  priori: hyper-g UIP + Binomial-Beta, mprior.size = K/2\n")
cat("  desfecho: casos | W: Queen | K =", K, "\n")
cat("  2 cadeias, 1.000.000 de iteracoes, 200.000 de descarte\n")
cat("  (identica ao inventario_estimacoes.csv do TCC)\n\n")

res <- list()
for (rotulo in c("MQO", "SAR")) {
  f <- if (rotulo == "MQO") lm_ols else lm_sar
  cat("Rodando", rotulo, "...\n"); flush.console()
  t0 <- proc.time()[3]
  a <- mc3(f, N, BURN, 20200731L, cortes)
  cat("  cadeia 1:", round(proc.time()[3] - t0, 1), "s\n"); flush.console()
  t0 <- proc.time()[3]
  b <- mc3(f, N, BURN, 987654L, cortes)
  cat("  cadeia 2:", round(proc.time()[3] - t0, 1), "s\n")
  cat("  aceitacao:", round(a$taxa, 3), "e", round(b$taxa, 3), "\n")
  cat("  tamanho medio:", round(a$tamanho_medio, 4), "e", round(b$tamanho_medio, 4), "\n")
  cat("  tamanho medio acumulado (cadeia 1 | cadeia 2):\n")
  for (ct in as.character(cortes)) {
    cat(sprintf("    ate %9s: %8.4f | %8.4f\n", ct,
                a$historico[[ct]], b$historico[[ct]]))
  }
  cat("  discordancia entre cadeias:", round(max(abs(a$pip - b$pip)), 4), "\n\n")
  res[[rotulo]] <- list(pip = (a$pip + b$pip) / 2,
                        conv = max(abs(a$pip - b$pip)),
                        tam = (a$tamanho_medio + b$tamanho_medio) / 2)
}

cat("  avaliacoes de R2 fora da tabela (caidas na serie exata):",
    fora_da_tabela, "\n\n")

pub <- read.csv(file.path(ESTIMACOES,
                          "resultados/estimacoes/principal/tabelas_tecnicas/resultado_casos.csv"),
                stringsAsFactors = FALSE)
pip_pub <- setNames(pub$PIP, pub$variavel)[covariaveis]

tab <- data.frame(covariavel = covariaveis,
                  publicada = as.numeric(pip_pub),
                  mqo = res$MQO$pip,
                  sar = res$SAR$pip,
                  stringsAsFactors = FALSE)
tab$mqo_menos_pub <- tab$mqo - tab$publicada
tab$sar_menos_mqo <- tab$sar - tab$mqo
tab <- tab[order(-tab$publicada), ]

cat("=== PIP: publicada, MQO reimplementado, SAR ===\n")
print(transform(tab, publicada = round(publicada, 4), mqo = round(mqo, 4),
                sar = round(sar, 4), mqo_menos_pub = round(mqo_menos_pub, 4),
                sar_menos_mqo = round(sar_menos_mqo, 4)), row.names = FALSE)

cat("\n=== Replicacao do capitulo 4 ===\n")
cat("  media da diferenca MQO - publicada:", round(mean(tab$mqo_menos_pub), 4), "\n")
cat("  maior |diferenca|:", round(max(abs(tab$mqo_menos_pub)), 4), "\n")
cat("  negativas:", sum(tab$mqo_menos_pub < 0), "de", K, "\n")
cat("  tamanho medio do modelo: publicado", round(sum(pip_pub), 4),
    "| MQO", round(res$MQO$tam, 4), "| SAR", round(res$SAR$tam, 4), "\n")

cat("\n=== Efeito espacial ===\n")
cat("  acima de 0,50 no MQO:",
    paste(tab$covariavel[tab$mqo >= 0.5], collapse = ", "), "\n")
cat("  acima de 0,50 no SAR:",
    ifelse(any(tab$sar >= 0.5),
           paste(tab$covariavel[tab$sar >= 0.5], collapse = ", "), "nenhuma"), "\n")
cat("  maior queda:", round(min(tab$sar_menos_mqo), 4), "em",
    tab$covariavel[which.min(tab$sar_menos_mqo)], "\n")
cat("  mudam mais de 0,10:", sum(abs(tab$sar_menos_mqo) > 0.10), "de", K, "\n")
cat("  ruido de Monte Carlo: MQO", round(res$MQO$conv, 4),
    "| SAR", round(res$SAR$conv, 4), "\n")

write.csv(tab, file.path(AQUI, "saida", "comparacao_1milhao.csv"), row.names = FALSE)
cat("\nGravado em saida/comparacao_1milhao.csv\n")
