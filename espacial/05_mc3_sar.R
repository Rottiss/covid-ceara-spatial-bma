# Etapa 3 do roteiro LeSage e Parent (2007): amostrador MC3 para o SAR.
#
# Implementa o MC3 de Madigan e York (1995) com os tres tipos de proposta que os
# autores usam: birth, death e move, cada um com probabilidade 1/3.
#
# RAZAO DE HASTINGS. O artigo afirma que a aceitacao min[1, p(M'|y)/p(M|y)] e
# valida "desde que birth, death e move tenham probabilidade igual de 1/3". Isso
# vale para a formulacao original de Madigan e York, em que se sorteia
# uniformemente na vizinhanca nbd(M) -- e como |nbd(M)| = K + 1 para todo M, a
# proposta e simetrica e a razao se cancela.
#
# Nao vale quando se sorteia primeiro o TIPO de passo e depois a variavel: nesse
# caso q(M -> M') = (1/n_tipos) * (1/n_escolhas) depende de k, e a proposta e
# assimetrica. Aqui a razao de Hastings entra explicitamente. A verificacao
# contra enumeracao exata abaixo mostra que sem ela o amostrador nao converge
# para o posterior correto.
#
# Verificacao: subconjunto de 10 covariaveis, em que os 2^10 = 1024 modelos podem
# ser enumerados com a propria log-marginal do SAR. As PIP da cadeia tem que
# convergir para as exatas.

source(if (file.exists("../caminhos.R")) "../caminhos.R" else "caminhos.R")

AQUI <- ESPACIAL
base <- readRDS(file.path(AQUI, "saida", "base_sar.rds"))
aux <- readRDS(file.path(AQUI, "saida", "grade_sar.rds"))

y <- base$y; X_total <- base$X; W <- base$W; mu <- base$mu
n <- base$n; covariaveis <- base$covariaveis
grade <- aux$grade; ta <- aux$ta; tb <- aux$tb; tc <- aux$tc
yc <- aux$yc; Wyc <- aux$Wyc

g_uip <- n

estatisticas <- function(indices, Xref) {
  if (length(indices) == 0L) {
    eo <- yc; ed <- Wyc
  } else {
    Xc <- Xref[, indices, drop = FALSE]
    Xc <- sweep(Xc, 2, colMeans(Xc), "-")
    Q <- qr.Q(qr(Xc))
    eo <- yc - Q %*% crossprod(Q, yc)
    ed <- Wyc - Q %*% crossprod(Q, Wyc)
  }
  eo <- as.numeric(eo); ed <- as.numeric(ed)
  c(k = length(indices), a = sum(eo^2), b = sum(eo * ed), c = sum(ed^2))
}

log_marginal <- function(est, g = g_uip) {
  rss <- est["a"] - 2 * grade$rho * est["b"] + grade$rho2 * est["c"]
  tss <- ta - 2 * grade$rho * tb + grade$rho2 * tc
  lf <- grade$ld_mais_prior -
    ((n - 1) / 2) * log((g / (1 + g)) * rss + (1 / (1 + g)) * tss)
  m <- max(lf)
  unname(-(est["k"] / 2) * log(1 + g) + m +
           log(grade$h / 3 * sum(grade$pesos * exp(lf - m))))
}

# ------------------------------------------------------------ amostrador
# Tipos disponiveis em cada estado: nao existe death com k = 0, nem birth com
# k = K, nem move em nenhum dos dois extremos.
tipos_disponiveis <- function(k, K) {
  tipos <- character(0)
  if (k < K) tipos <- c(tipos, "birth")
  if (k > 0L) tipos <- c(tipos, "death")
  if (k > 0L && k < K) tipos <- c(tipos, "move")
  tipos
}

# log q(M -> M'), para a razao de Hastings.
log_q <- function(k, K, tipo) {
  nt <- length(tipos_disponiveis(k, K))
  escolhas <- switch(tipo,
                     birth = K - k,
                     death = k,
                     move = k * (K - k))
  -log(nt) - log(escolhas)
}

mc3 <- function(Xref, K, n_sorteios, burn, semente, inicio = NULL, hastings = TRUE) {
  set.seed(semente)
  atual <- if (is.null(inicio)) sort(sample(seq_len(K), sample(0:K, 1))) else inicio
  est_atual <- estatisticas(atual, Xref)
  lm_atual <- log_marginal(est_atual)

  contagem <- new.env(hash = TRUE, parent = emptyenv())
  aceitos <- 0L
  inclusao <- numeric(K)
  guardados <- 0L

  for (it in seq_len(n_sorteios)) {
    k <- length(atual)
    tipos <- tipos_disponiveis(k, K)
    tipo <- tipos[sample.int(length(tipos), 1L)]
    fora <- setdiff(seq_len(K), atual)

    proposto <- switch(
      tipo,
      birth = sort(c(atual, fora[sample.int(length(fora), 1L)])),
      death = atual[-sample.int(k, 1L)],
      move  = sort(c(atual[-sample.int(k, 1L)], fora[sample.int(length(fora), 1L)]))
    )

    est_prop <- estatisticas(proposto, Xref)
    lm_prop <- log_marginal(est_prop)

    # Priori uniforme sobre modelos: a razao de prioris vale 1.
    log_alfa <- lm_prop - lm_atual
    if (hastings) {
      tipo_inverso <- switch(tipo, birth = "death", death = "birth", move = "move")
      log_alfa <- log_alfa +
        log_q(length(proposto), K, tipo_inverso) - log_q(k, K, tipo)
    }

    if (log(runif(1)) < log_alfa) {
      atual <- proposto; lm_atual <- lm_prop
      aceitos <- aceitos + 1L
    }

    if (it > burn) {
      guardados <- guardados + 1L
      if (length(atual) > 0L) inclusao[atual] <- inclusao[atual] + 1
      chave <- if (length(atual) == 0L) "vazio" else paste(atual, collapse = ",")
      contagem[[chave]] <- (contagem[[chave]] %||% 0L) + 1L
    }
  }

  list(pip = inclusao / guardados,
       taxa_aceitacao = aceitos / n_sorteios,
       modelos_unicos = length(ls(contagem)),
       contagem = contagem,
       guardados = guardados)
}

`%||%` <- function(a, b) if (is.null(a)) b else a

# ---------------------------------------- verificacao contra enumeracao
K_TESTE <- 10L
X_teste <- X_total[, seq_len(K_TESTE), drop = FALSE]

cat("=== Referencia: enumeracao completa dos 2^10 = 1024 modelos ===\n")
t0 <- proc.time()[3]
combinacoes <- lapply(0:(2^K_TESTE - 1),
                      function(i) which(bitwAnd(i, 2^(0:(K_TESTE - 1))) > 0))
logml <- vapply(combinacoes, function(idx) log_marginal(estatisticas(idx, X_teste)),
                numeric(1))
m <- max(logml)
pmp <- exp(logml - m) / sum(exp(logml - m))
pip_exata <- vapply(seq_len(K_TESTE), function(j)
  sum(pmp[vapply(combinacoes, function(s) j %in% s, logical(1))]), numeric(1))
cat("  tempo:", round(proc.time()[3] - t0, 1), "s\n")
cat("  modelo de maior PMP:",
    paste(colnames(X_teste)[combinacoes[[which.max(pmp)]]], collapse = " + "),
    "com PMP", round(max(pmp), 4), "\n\n")

cat("=== Cadeias MC3 ===\n")
n_sorteios <- 250000L; burn <- 50000L

t0 <- proc.time()[3]
c1 <- mc3(X_teste, K_TESTE, n_sorteios, burn, semente = 20200731L)
c2 <- mc3(X_teste, K_TESTE, n_sorteios, burn, semente = 987654L)
tempo <- proc.time()[3] - t0
cat("  duas cadeias de", format(n_sorteios, big.mark = "."), "sorteios em",
    round(tempo, 1), "s\n")
cat("  taxa de aceitacao:", round(c1$taxa_aceitacao, 3), "e",
    round(c2$taxa_aceitacao, 3), "\n")
cat("  modelos unicos visitados:", c1$modelos_unicos, "e", c2$modelos_unicos,
    "de 1024\n\n")

comparacao <- data.frame(
  covariavel = colnames(X_teste),
  pip_exata = pip_exata,
  pip_cadeia1 = c1$pip,
  pip_cadeia2 = c2$pip
)
comparacao$erro1 <- comparacao$pip_cadeia1 - comparacao$pip_exata
comparacao$erro2 <- comparacao$pip_cadeia2 - comparacao$pip_exata

cat("=== PIP: exata vs. cadeias ===\n")
print(transform(comparacao,
                pip_exata = round(pip_exata, 5),
                pip_cadeia1 = round(pip_cadeia1, 5),
                pip_cadeia2 = round(pip_cadeia2, 5),
                erro1 = signif(erro1, 3), erro2 = signif(erro2, 3)),
      row.names = FALSE)

erro_max <- max(abs(c(comparacao$erro1, comparacao$erro2)))
entre_cadeias <- max(abs(comparacao$pip_cadeia1 - comparacao$pip_cadeia2))
cat("\n  maior erro contra a enumeracao:", signif(erro_max, 3), "\n")
cat("  maior diferenca entre as duas cadeias:", signif(entre_cadeias, 3), "\n")
cat("  veredito:", if (erro_max < 0.01) "APROVADO" else "INVESTIGAR", "\n\n")

# ------------------------------- o que acontece sem a razao de Hastings
cat("=== Controle: a mesma cadeia sem a razao de Hastings ===\n")
c_sem <- mc3(X_teste, K_TESTE, n_sorteios, burn, semente = 20200731L, hastings = FALSE)
erro_sem <- max(abs(c_sem$pip - pip_exata))
cat("  maior erro contra a enumeracao:", signif(erro_sem, 3), "\n")
cat("  (com a razao de Hastings:", signif(erro_max, 3), ")\n")
cat("  razao entre os erros:", round(erro_sem / erro_max, 1), "vezes\n\n")

controle <- data.frame(
  covariavel = colnames(X_teste),
  pip_exata = round(pip_exata, 5),
  com_hastings = round(c1$pip, 5),
  sem_hastings = round(c_sem$pip, 5)
)
print(controle, row.names = FALSE)

write.csv(comparacao, file.path(AQUI, "saida", "mc3_validacao.csv"), row.names = FALSE)
write.csv(controle, file.path(AQUI, "saida", "mc3_controle_hastings.csv"), row.names = FALSE)

cat("\nGravado em", file.path(AQUI, "saida"), "\n")
