# Diagnostica a diferenca sistematica de -0,0158 entre as PIP da implementacao
# propria e as publicadas no capitulo 4.
#
# HIPOTESE 1 (estimador da PIP) JA ESTA DESCARTADA. O covid_bma.R chama
# coef(..., exact = FALSE), que no BMS significa PIP por frequencia de visita --
# o mesmo estimador usado aqui. Nao ha diferenca de estimador.
#
# HIPOTESE 2: truncamento da hipergeometrica no BMS. O .f21_4hyperg corta a
# serie em 200, 600, 1.400 ou 3.000 termos por heuristica. Com K = 30 os R^2 sao
# maiores do que no teste de validacao com K = 10.
#
# Este script compara as duas log-marginais DIRETAMENTE, modelo a modelo, sem
# MCMC. Se elas baterem, a hipotese 2 tambem cai e o vies vem de outro lugar.
#
# HIPOTESE 3, testada em seguida: mistura da cadeia. Recalcula as PIP proprias
# ponderando os modelos unicos visitados pela marginal analitica, em vez de por
# frequencia. Se as duas versoes divergirem, a cadeia nao esta bem misturada.

source(if (file.exists("../caminhos.R")) "../caminhos.R" else "caminhos.R")

suppressPackageStartupMessages(library(BMS))

AQUI <- ESPACIAL
base <- readRDS(file.path(AQUI, "saida", "base_sar.rds"))
aux <- readRDS(file.path(AQUI, "saida", "grade_sar.rds"))
hyp <- readRDS(file.path(AQUI, "saida", "hyperg_aux.rds"))

y <- base$y; X_total <- base$X; n <- base$n; covariaveis <- base$covariaveis
yc <- aux$yc; ta <- aux$ta
K <- length(covariaveis)
f21a <- hyp$f21a; a_par <- hyp$a_par; L <- hyp$L; J_MAX <- hyp$J_MAX

log_f21_vec <- function(z, k) {
  cc <- (f21a + k) / 2
  j_pico <- max(0, ceiling(max((z * a_par - cc) / pmax(1 - z, 1e-12))))
  J <- min(J_MAX, max(300L, as.integer(j_pico * 4)))
  lt <- cbind(0, outer(log(z), seq_len(J)) +
                rep(L[[k + 1]][seq_len(J)], each = length(z)))
  m <- apply(lt, 1, max)
  m + log(rowSums(exp(lt - m)))
}

minha_loglik <- function(ymy, k) {
  nl <- (1 - n) / 2 * log(ta)
  if (k == 0L) return(nl)
  unname(nl + log(f21a - 2) - log(k + f21a - 2) + log_f21_vec(1 - ymy / ta, k))
}

# ------------------- log-marginal do proprio BMS, sem passar por bms()
gpo <- BMS:::.gprior.hyperg.init(g = "hyper=UIP", N = n, K = K, yty = ta,
                                 null.lik = (1 - n) / 2 * log(ta))
bms_loglik <- gpo$lprobcalc$just.loglik
cat("=== Objeto de priori g do BMS ===\n")
cat("  a =", format(gpo$hyper.parameter, digits = 12), "\n")
cat("  a da implementacao =", format(f21a, digits = 12), "\n\n")

# --------------------- comparacao modelo a modelo, sobre modelos reais
cat("=== Log-marginais, modelo a modelo, K variando de 1 a 30 ===\n")
set.seed(20200731L)
amostras <- c(
  lapply(1:30, function(k) seq_len(k)),                       # aninhados
  replicate(120, sort(sample(K, sample(1:30, 1))), simplify = FALSE)
)

linhas <- lapply(amostras, function(idx) {
  k <- length(idx)
  Xc <- X_total[, idx, drop = FALSE]
  Xc <- sweep(Xc, 2, colMeans(Xc), "-")
  Q <- qr.Q(qr(Xc))
  ymy <- ta - sum((crossprod(Q, yc))^2)
  data.frame(k = k, R2 = 1 - ymy / ta,
             minha = minha_loglik(ymy, k),
             bms = bms_loglik(ymy, k))
})
tab <- do.call(rbind, linhas)
tab$dif <- tab$minha - tab$bms

cat("  modelos comparados:", nrow(tab), "\n")
cat("  R2 varia de", round(min(tab$R2), 4), "a", round(max(tab$R2), 4), "\n")
cat("  maior |diferenca| em log-marginal:", format(max(abs(tab$dif)), scientific = TRUE), "\n")
cat("  diferenca media:", format(mean(tab$dif), scientific = TRUE), "\n\n")

por_k <- aggregate(dif ~ k, data = tab, FUN = function(z) max(abs(z)))
names(por_k)[2] <- "max_abs_dif"
por_k$R2_max <- aggregate(R2 ~ k, data = tab, FUN = max)$R2
piores <- por_k[order(-por_k$max_abs_dif), ][1:8, ]
cat("  oito piores k:\n")
print(transform(piores, max_abs_dif = signif(max_abs_dif, 3),
                R2_max = round(R2_max, 4)), row.names = FALSE)

if (max(abs(tab$dif)) < 1e-8) {
  cat("\n  HIPOTESE 2 DESCARTADA: as duas log-marginais sao a mesma funcao.\n\n")
} else {
  cat("\n  HIPOTESE 2 SUSTENTADA: as log-marginais divergem.\n\n")
}

# ---------------------------------- hipotese 3: mistura da cadeia
cat("=== Hipotese 3: PIP por frequencia vs. por peso analitico ===\n")
cat("  Mesma cadeia, dois estimadores. Se divergirem, a cadeia nao misturou.\n\n")

lpm <- vapply(0:K, function(k) lgamma(1 + k) + lgamma(1 + K - k), numeric(1))

estat_ymy <- function(indices) {
  if (length(indices) == 0L) return(ta)
  Xc <- X_total[, indices, drop = FALSE]
  Xc <- sweep(Xc, 2, colMeans(Xc), "-")
  Q <- qr.Q(qr(Xc))
  ta - sum((crossprod(Q, yc))^2)
}

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

mc3_com_pesos <- function(n_sorteios, burn, semente) {
  set.seed(semente)
  atual <- sort(sample(seq_len(K), sample(0:K, 1)))
  post_atual <- minha_loglik(estat_ymy(atual), length(atual)) + lpm[length(atual) + 1]
  inclusao <- numeric(K); guardados <- 0L
  reg <- new.env(hash = TRUE, parent = emptyenv())
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
    post_prop <- minha_loglik(estat_ymy(prop), length(prop)) + lpm[length(prop) + 1]
    inv <- switch(tipo, birth = "death", death = "birth", move = "move")
    la <- post_prop - post_atual + log_q(length(prop), inv) - log_q(k, tipo)
    if (is.finite(la) && log(runif(1)) < la) {
      atual <- prop; post_atual <- post_prop
    }
    if (it > burn) {
      guardados <- guardados + 1L
      if (length(atual) > 0L) inclusao[atual] <- inclusao[atual] + 1
      chave <- if (length(atual) == 0L) "vazio" else paste(atual, collapse = ",")
      if (is.null(reg[[chave]])) assign(chave, post_atual, envir = reg)
    }
  }
  chaves <- ls(reg)
  lpost <- vapply(chaves, function(ch) get(ch, envir = reg), numeric(1))
  m <- max(lpost)
  w <- exp(lpost - m); w <- w / sum(w)
  pip_analitica <- numeric(K)
  for (i in seq_along(chaves)) {
    if (chaves[i] == "vazio") next
    idx <- as.integer(strsplit(chaves[i], ",")[[1]])
    pip_analitica[idx] <- pip_analitica[idx] + w[i]
  }
  list(freq = inclusao / guardados, analitica = pip_analitica,
       unicos = length(chaves))
}

N <- as.integer(Sys.getenv("MC3_SORTEIOS", "200000"))
t0 <- proc.time()[3]
a <- mc3_com_pesos(N, as.integer(N * 0.2), 20200731L)
b <- mc3_com_pesos(N, as.integer(N * 0.2), 987654L)
cat("  duas cadeias de", N, "sorteios em", round(proc.time()[3] - t0, 1), "s\n")
cat("  modelos unicos:", a$unicos, "e", b$unicos, "\n\n")

pub <- read.csv(file.path(ESTIMACOES,
                          "resultados/estimacoes/principal/tabelas_tecnicas/resultado_casos.csv"),
                stringsAsFactors = FALSE)
pip_pub <- setNames(pub$PIP, pub$variavel)[covariaveis]

comp <- data.frame(
  covariavel = covariaveis,
  publicada = as.numeric(pip_pub),
  freq = (a$freq + b$freq) / 2,
  analitica = (a$analitica + b$analitica) / 2
)
comp$dif_freq <- comp$freq - comp$publicada
comp$dif_analitica <- comp$analitica - comp$publicada
comp <- comp[order(-comp$publicada), ]

print(transform(comp, publicada = round(publicada, 4), freq = round(freq, 4),
                analitica = round(analitica, 4),
                dif_freq = round(dif_freq, 4),
                dif_analitica = round(dif_analitica, 4)), row.names = FALSE)

cat("\n  media da diferenca, por frequencia:", round(mean(comp$dif_freq), 4),
    "| negativas:", sum(comp$dif_freq < 0), "de", K, "\n")
cat("  media da diferenca, analitica    :", round(mean(comp$dif_analitica), 4),
    "| negativas:", sum(comp$dif_analitica < 0), "de", K, "\n")
cat("  maior |diferenca| freq:", round(max(abs(comp$dif_freq)), 4),
    "| analitica:", round(max(abs(comp$dif_analitica)), 4), "\n")

write.csv(comp, file.path(AQUI, "saida", "diagnostico_vies.csv"), row.names = FALSE)
write.csv(tab, file.path(AQUI, "saida", "logmarginais_bms_vs_propria.csv"), row.names = FALSE)
cat("\nGravado em", file.path(AQUI, "saida"), "\n")
