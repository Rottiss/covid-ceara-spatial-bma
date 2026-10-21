# Etapa 2 da cadeia de tres vias: o passo de troca entre MQO, SAR e SEM esta
# correto?
#
# Com duas especificacoes o passo de troca tem um unico destino, a proposta e
# simetrica e a razao de Hastings e 1. Com tres, o passo escolhe entre dois
# destinos, e a pergunta e se a razao continua 1. Ela continua, desde que o
# destino seja sorteado uniformemente: q depende so de k, que a troca nao muda.
# Mas esse e exatamente o tipo de raciocinio que ja falhou neste projeto uma vez
# (a razao de Hastings dos passos de covariavel, erro de ate 0,185 na PIP), e a
# unica prova que vale e contra enumeracao exata.
#
# Desenho: K = 10, logo 1024 modelos x 3 especificacoes = 3072 estados. Todas as
# marginais sao pre-computadas, de modo que a cadeia consulta uma tabela. Isso
# isola o que se quer testar, que e o amostrador, do que a etapa 1 ja validou,
# que e a marginal.
#
# Tres cadeias sao comparadas com o exato:
#   duas vias  (MQO, SAR)        controle: reproduz o 33, que ja e usado no TCC
#   tres vias  (MQO, SAR, SEM)   o objeto do teste
#   tres vias por ciclo          erro deliberado, para mostrar que o teste pega
#
# Nota sobre o termo log(ne - 1) em log_q: ele cancela. A troca nao muda k, logo
# log_q(k, "troca") aparece identico nos dois lados da razao de Hastings e some,
# qualquer que seja o numero de especificacoes. A primeira versao deste script
# usava justamente esse termo como erro deliberado e as tres cadeias saiam
# identicas, o que mostrou que o suposto erro era equivalencia algebrica. O erro
# que de fato quebra o balanco detalhado e assimetria na escolha do destino, e e
# esse que o terceiro caso reproduz.
#
# So le. Nao grava nada.

source(if (file.exists("../caminhos.R")) "../caminhos.R" else "caminhos.R")

suppressPackageStartupMessages(library(spdep))
AQUI <- ESPACIAL
RAIZ <- ESTIMACOES
hyp <- readRDS(file.path(AQUI, "saida", "hyperg_aux.rds"))
tb  <- readRDS(file.path(AQUI, "saida", "tabela_f21_estendida.rds"))
f21a <- hyp$f21a; a_par <- hyp$a_par; L <- hyp$L; J_MAX <- hyp$J_MAX
TAB <- tb$TAB; passo <- tb$passo; N_TAB <- tb$N_TAB; Z_MAX <- tb$Z_MAX

DESFECHO <- "casos"; COLUNA <- "z_log_casos100k"
dados <- read.csv(file.path(RAIZ, "resultados/estimacoes/matrizes_estimacao/principal",
                            paste0(DESFECHO, ".csv")),
                  stringsAsFactors = FALSE, colClasses = c(codigo_ibge = "character"))
todas <- setdiff(names(dados), c("codigo_ibge", "municipio", COLUNA))
set.seed(20200731L)
covs <- todas[sort(sample(seq_along(todas), 10))]
y <- as.numeric(dados[[COLUNA]])
X <- as.matrix(dados[, covs]); n <- length(y); K <- ncol(X)
Xc <- sweep(X, 2, colMeans(X), "-")
listw <- readRDS(file.path(RAIZ, "resultados/diagnostico_espacial", "pesos_queen.rds"))
W <- spdep::listw2mat(listw); dimnames(W) <- NULL

fora_tab <- 0L
log_f21_serie1 <- function(zi, k) {
  cc <- (f21a + k) / 2
  j_pico <- max(0, ceiling((zi * a_par - cc) / max(1 - zi, 1e-12)))
  J <- min(J_MAX, max(300L, as.integer(j_pico * 4)))
  lt <- c(0, log(zi) * seq_len(J) + L[[k + 1]][seq_len(J)])
  m <- max(lt); m + log(sum(exp(lt - m)))
}
log_f21 <- function(z, k) {
  acima <- z > Z_MAX
  if (any(acima)) {
    fora_tab <<- fora_tab + sum(acima)
    out <- numeric(length(z))
    if (any(!acima)) out[!acima] <- log_f21(z[!acima], k)
    out[acima] <- vapply(z[acima], log_f21_serie1, numeric(1), k = k)
    return(out)
  }
  pos <- z / passo; i <- as.integer(pos) + 1L
  i[i >= N_TAB] <- N_TAB - 1L
  frac <- pos - (i - 1L); off <- k * N_TAB
  TAB[i + off] * (1 - frac) + TAB[i + 1L + off] * frac
}
mu <- Re(eigen(W, only.values = TRUE)$values); a0 <- 1.01
grade <- function(ngp) {
  v <- seq(-0.999, 0.999, length.out = ngp)
  ld <- as.numeric(outer(v, mu, function(r, m) log1p(-r * m)) %*% rep(1, length(mu)))
  lpr <- (a0 - 1) * log1p(v) + (a0 - 1) * log1p(-v) -
    (2 * a0 - 1) * log(2) - lbeta(a0, a0)
  list(v = v, ldp = ld + lpr, pes = c(1, rep(c(4, 2), length.out = ngp - 2), 1),
       h = v[2] - v[1])
}
gr <- grade(201L); gl <- grade(101L)

yc <- y - mean(y); Wy <- as.numeric(W %*% y); Wyc <- Wy - mean(Wy)
ta <- sum(yc^2); tbb <- sum(yc * Wyc); tcc <- sum(Wyc^2)
estat <- function(idx) {
  if (length(idx) == 0L) { eo <- yc; ed <- Wyc } else {
    Q <- qr.Q(qr(Xc[, idx, drop = FALSE]))
    eo <- yc - Q %*% crossprod(Q, yc); ed <- Wyc - Q %*% crossprod(Q, Wyc)
  }
  eo <- as.numeric(eo); ed <- as.numeric(ed)
  c(k = length(idx), a = sum(eo^2), b = sum(eo * ed), c = sum(ed^2))
}
lm_ols <- function(e) {
  k <- unname(e["k"]); nl <- (1 - n) / 2 * log(ta)
  if (k == 0L) return(nl)
  unname(nl + log(f21a - 2) - log(k + f21a - 2) + log_f21(1 - e["a"] / ta, k))
}
lm_sar <- function(e) {
  k <- unname(e["k"])
  rss <- e["a"] - 2 * gr$v * e["b"] + gr$v^2 * e["c"]
  tss <- ta - 2 * gr$v * tbb + gr$v^2 * tcc
  lf <- gr$ldp + (1 - n) / 2 * log(tss)
  if (k > 0L) lf <- lf + log(f21a - 2) - log(k + f21a - 2) +
      log_f21(pmax(0, 1 - rss / tss), k)
  m <- max(lf); m + log(gr$h / 3 * sum(gr$pes * exp(lf - m)))
}
WX <- W %*% X
A00 <- crossprod(X); A01 <- crossprod(X, WX); A11 <- crossprod(WX)
b0 <- as.numeric(crossprod(X, y));  b1 <- as.numeric(crossprod(X, Wy))
c0 <- as.numeric(crossprod(WX, y)); c1 <- as.numeric(crossprod(WX, Wy))
s0 <- colSums(X); s1 <- colSums(WX)
yy <- sum(y^2); yWy <- sum(y * Wy); WyWy <- sum(Wy^2)
sy <- sum(y); sWy <- sum(Wy)
ys_ys <- yy - 2 * gl$v * yWy + gl$v^2 * WyWy
soma_ys <- sy - gl$v * sWy
tss_s <- ys_ys - n * (soma_ys / n)^2
lm_sem <- function(idx) {
  k <- length(idx)
  if (k == 0L) rss <- tss_s else {
    A00i <- A00[idx, idx, drop = FALSE]; A01i <- A01[idx, idx, drop = FALSE]
    A11i <- A11[idx, idx, drop = FALSE]
    b0i <- b0[idx]; b1i <- b1[idx]; c0i <- c0[idx]; c1i <- c1[idx]
    s0i <- s0[idx]; s1i <- s1[idx]
    sim <- A01i + t(A01i); soma_bc <- b1i + c0i
    rss <- numeric(length(gl$v))
    for (j in seq_along(gl$v)) {
      Lj <- gl$v[j]
      Xti <- s0i - Lj * s1i
      XtX <- A00i - Lj * sim + Lj^2 * A11i - tcrossprod(Xti) / n
      Xty <- b0i - Lj * soma_bc + Lj^2 * c1i - Xti * (soma_ys[j] / n)
      bt <- tryCatch(solve(XtX, Xty), error = function(e) NULL)
      rss[j] <- if (is.null(bt)) NA_real_ else tss_s[j] - sum(bt * Xty)
    }
  }
  if (anyNA(rss)) return(-Inf)
  lf <- gl$ldp + (1 - n) / 2 * log(tss_s)
  if (k > 0L) lf <- lf + log(f21a - 2) - log(k + f21a - 2) +
      log_f21(pmax(0, pmin(1 - rss / tss_s, 1 - 1e-12)), k)
  m <- max(lf); m + log(gl$h / 3 * sum(gl$pes * exp(lf - m)))
}

# ------------------------------------------------- tabela completa de estados
ESP <- c("MQO", "SAR", "SEM")
modelos <- lapply(0:(2^K - 1), function(m) which(bitwAnd(m, 2^(0:(K - 1))) > 0))
lpm <- vapply(0:K, function(k) lgamma(1 + k) + lgamma(1 + K - k), numeric(1))
cat("pre-computando", length(modelos) * 3, "estados...\n")
t0 <- proc.time()[3]
M <- matrix(0, nrow = length(modelos), ncol = 3, dimnames = list(NULL, ESP))
for (i in seq_along(modelos)) {
  idx <- modelos[[i]]; e <- estat(idx)
  M[i, "MQO"] <- lm_ols(e); M[i, "SAR"] <- lm_sar(e); M[i, "SEM"] <- lm_sem(idx)
}
tam <- vapply(modelos, length, integer(1))
POST <- M + lpm[tam + 1]          # priori de especificacao uniforme: constante
cat(sprintf("  %.1f s | avaliacoes fora da tabela de f21: %d\n\n",
            proc.time()[3] - t0, fora_tab))

exato <- function(cols) {
  p <- POST[, cols, drop = FALSE]
  w <- exp(p - max(p)); w <- w / sum(w)
  pip <- vapply(seq_len(K), function(j)
    sum(w[vapply(modelos, function(s) j %in% s, logical(1)), , drop = FALSE]), numeric(1))
  list(pip = pip, esp = colSums(w))
}

# ------------------------------------------------- amostrador
id_de <- function(idx) if (!length(idx)) 1L else sum(2^(idx - 1)) + 1L
tipos_disp <- function(k) {
  t <- "troca"
  if (k < K) t <- c(t, "birth"); if (k > 0L) t <- c(t, "death")
  if (k > 0L && k < K) t <- c(t, "move"); t
}
# bug = "ciclo": em vez de sortear o destino entre os outros dois, percorre
# MQO -> SAR -> SEM -> MQO. E um erro plausivel de implementacao e viola o
# balanco detalhado, porque o movimento inverso tem probabilidade zero.
mc3 <- function(cols, semente, n_mc, burn, bug = FALSE) {
  set.seed(semente)
  ne <- length(cols)
  log_q <- function(k, tipo)
    -log(length(tipos_disp(k))) -
    log(switch(tipo, troca = ne - 1, birth = K - k, death = k, move = k * (K - k)))
  atual <- sort(sample(seq_len(K), sample(0:K, 1))); esp <- cols[1]
  post <- POST[id_de(atual), esp]
  inc <- numeric(K); ce <- setNames(numeric(ne), cols); g <- 0L
  for (it in seq_len(n_mc)) {
    k <- length(atual); disp <- tipos_disp(k)
    tp <- disp[sample.int(length(disp), 1L)]
    if (tp == "troca") {
      prop <- atual
      alt <- setdiff(cols, esp)
      esp_p <- if (bug) cols[(match(esp, cols) %% ne) + 1L]
               else alt[sample.int(length(alt), 1L)]
    } else {
      esp_p <- esp; fora <- setdiff(seq_len(K), atual)
      prop <- switch(tp,
                     birth = sort(c(atual, fora[sample.int(length(fora), 1L)])),
                     death = atual[-sample.int(k, 1L)],
                     move  = sort(c(atual[-sample.int(k, 1L)],
                                    fora[sample.int(length(fora), 1L)])))
    }
    pp <- POST[id_de(prop), esp_p]
    inv <- switch(tp, troca = "troca", birth = "death", death = "birth", move = "move")
    la <- pp - post + log_q(length(prop), inv) - log_q(k, tp)
    if (is.finite(la) && log(runif(1)) < la) { atual <- prop; esp <- esp_p; post <- pp }
    if (it > burn) {
      g <- g + 1L
      if (length(atual)) inc[atual] <- inc[atual] + 1
      ce[esp] <- ce[esp] + 1
    }
  }
  list(pip = inc / g, esp = ce / g)
}

N_MC <- 400000L; BURN <- 80000L
casos <- list(
  list(rot = "duas vias (controle)", cols = c("MQO", "SAR"), bug = FALSE),
  list(rot = "tres vias",            cols = ESP,             bug = FALSE),
  list(rot = "tres vias por ciclo (erro)", cols = ESP,      bug = TRUE))

for (cs in casos) {
  ex <- exato(cs$cols)
  a <- mc3(cs$cols, 20200731L, N_MC, BURN, cs$bug)
  b <- mc3(cs$cols, 987654L,  N_MC, BURN, cs$bug)
  pip <- (a$pip + b$pip) / 2; esp <- (a$esp + b$esp) / 2
  # Margem numerica para este teste de Monte Carlo com sementes fixas.
  # O caso com erro deliberado e apenas um controle negativo.
  if (!cs$bug && (max(abs(pip - ex$pip)) > 0.03 || max(abs(esp - ex$esp)) > 0.03))
    stop("A cadeia de validacao divergiu da enumeracao: ", cs$rot)
  cat(sprintf("=== %s ===\n", cs$rot))
  cat(sprintf("  erro maximo de PIP contra o exato: %.5f\n", max(abs(pip - ex$pip))))
  cat(sprintf("  discordancia entre as duas cadeias: %.5f\n", max(abs(a$pip - b$pip))))
  for (e in cs$cols)
    cat(sprintf("  P(%s|y)  exato %.5f | cadeia %.5f | erro %+.5f\n",
                e, ex$esp[[e]], esp[[e]], esp[[e]] - ex$esp[[e]]))
  cat("\n")
}
