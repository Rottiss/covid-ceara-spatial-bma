# SEM sob a priori do TCC (hyper-g UIP + Binomial-Beta), na mesma configuracao
# da grade empirica.
#
# POR QUE ISTO IMPORTA. Toda a evidencia do artigo vem do modelo de defasagem
# espacial. Um leitor pode objetar que o resultado e propriedade do SAR e nao da
# dependencia espacial: a defasagem Wy compete diretamente com as covariaveis
# pela mesma variacao suave no espaco. O SEM poe a dependencia no erro, sem
# competidor no lado direito. Se as PIP tambem colapsarem, o achado e sobre
# dependencia espacial; se sobreviverem, e sobre a defasagem.
#
# O 07_mc3_sem.R ja fazia isso com priori uniforme e g fixo, em 120 mil
# sorteios. Aqui a comparacao fica direta com a grade: mesma priori, mesmo
# numero de sorteios, mesmas duas cadeias.
#
#   DESFECHO = casos | obitos | letalidade      VIZINHOS = queen | rook
#
# CUSTO. No SAR tudo o que depende de rho cabe em tres escalares por modelo. No
# SEM a projecao muda a cada lambda, entao ha um sistema k x k por ponto da
# grade. Os blocos de produto cruzado sao pre-computados uma vez, o que reduz o
# custo de uma decomposicao QR de n x k para um solve de k x k.

source(if (file.exists("../caminhos.R")) "../caminhos.R" else "caminhos.R")

suppressPackageStartupMessages(library(spdep))

AQUI <- ESPACIAL
RAIZ <- ESTIMACOES
hyp <- readRDS(file.path(AQUI, "saida", "hyperg_aux.rds"))
tb  <- readRDS(file.path(AQUI, "saida", "tabela_f21_estendida.rds"))

DESFECHO <- Sys.getenv("DESFECHO", "casos")
VIZINHOS <- Sys.getenv("VIZINHOS", "queen")
N_MC <- as.integer(Sys.getenv("MC3_SORTEIOS", "1000000"))
BURN <- as.integer(Sys.getenv("MC3_BURN", "200000"))

COLUNA <- c(casos = "z_log_casos100k", obitos = "z_log_obitos100k",
            letalidade = "z_logit_letal")[[DESFECHO]]

dados <- read.csv(file.path(RAIZ, "resultados/estimacoes/matrizes_estimacao/principal",
                            paste0(DESFECHO, ".csv")),
                  stringsAsFactors = FALSE, colClasses = c(codigo_ibge = "character"))
covariaveis <- setdiff(names(dados), c("codigo_ibge", "municipio", COLUNA))
y <- as.numeric(dados[[COLUNA]])
X <- as.matrix(dados[, covariaveis])
n <- length(y); K <- length(covariaveis)
# Centragem de X: exigida pela prescricao da pagina 248 do artigo, ver 07.
X <- sweep(X, 2, colMeans(X), "-")

listw <- readRDS(file.path(RAIZ, "resultados/diagnostico_espacial",
                           paste0("pesos_", VIZINHOS, ".rds")))
W <- spdep::listw2mat(listw); dimnames(W) <- NULL
stopifnot(nrow(W) == n, max(abs(rowSums(W) - 1)) < 1e-9, max(abs(diag(W))) < 1e-12)

f21a <- hyp$f21a; a_par <- hyp$a_par; L <- hyp$L; J_MAX <- hyp$J_MAX
TAB <- tb$TAB; passo <- tb$passo; N_TAB <- tb$N_TAB; Z_MAX <- tb$Z_MAX

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

# ----------------------------------------------- grade de lambda
mu <- Re(eigen(W, only.values = TRUE)$values)
# Resolucao da grade de lambda. O SAR usa 201 pontos, mas la o custo por ponto e
# escalar. Aqui cada ponto custa um solve k x k, entao vale medir: contra uma
# grade de 801 pontos, o erro maximo da log-marginal e 1,1e-10 com 101 pontos e
# 1,1e-10 com 201; com 51 sobe para 2,3e-05. Fica em 101, com metade do custo.
ng <- as.integer(Sys.getenv("N_GRADE_LAMBDA", "101"))
# Simpson exige numero par de intervalos, ou seja, numero impar de pontos. Com
# ng par os pesos abaixo saem errados (terminam em 2 em vez de 4) e a integracao
# fica silenciosamente viesada. Como ng vem de variavel de ambiente documentada,
# a checagem evita o caso.
stopifnot(ng >= 3L, ng %% 2L == 1L)
lambda <- seq(-0.999, 0.999, length.out = ng)
ld <- as.numeric(outer(lambda, mu, function(r, m) log1p(-r * m)) %*% rep(1, length(mu)))
a0 <- 1.01
lpr <- (a0 - 1) * log1p(lambda) + (a0 - 1) * log1p(-lambda) -
  (2 * a0 - 1) * log(2) - lbeta(a0, a0)
ld_mais_prior <- ld + lpr
pesos <- c(1, rep(c(4, 2), length.out = ng - 2), 1); h <- lambda[2] - lambda[1]

# ----------------------------------------------- blocos pre-computados
Wy <- as.numeric(W %*% y); WX <- W %*% X
A00 <- crossprod(X); A01 <- crossprod(X, WX); A11 <- crossprod(WX)
b0 <- as.numeric(crossprod(X, y));  b1 <- as.numeric(crossprod(X, Wy))
c0 <- as.numeric(crossprod(WX, y)); c1 <- as.numeric(crossprod(WX, Wy))
s0 <- colSums(X); s1 <- colSums(WX)
yy <- sum(y^2); yWy <- sum(y * Wy); WyWy <- sum(Wy^2)
sy <- sum(y); sWy <- sum(Wy)

ys_ys <- yy - 2 * lambda * yWy + lambda^2 * WyWy
soma_ys <- sy - lambda * sWy
alfa <- soma_ys / n
tss_s <- ys_ys - n * alfa^2

# CENTRAGEM DEPOIS DA FILTRAGEM. X chega centrada, mas X*(lambda) = X - lambda WX
# nao e: a media de cada coluna de WX nao e zero. Entao o intercepto e beta nao
# sao ortogonais em lambda diferente de zero, e resolver X*'X* beta = X*'y*
# devolve um beta que nao e o do ajuste com intercepto. A verificacao em
# lambda = 0 nao pega o erro, porque nesse ponto X* = X ja esta centrada; so a
# comparacao com o ajuste explicito em lambda diferente de zero pega, e foi
# assim que apareceu, com diferenca de RSS de 3,7e-2. O mesmo defeito estava em
# 07_mc3_sem.R, cujos resultados foram descartados.
#
# A correcao usa os produtos cruzados centrados, obtidos dos nao centrados sem
# refazer nenhuma multiplicacao de matriz:
#     Xc*'Xc* = X*'X* - (X*'i)(X*'i)'/n      Xc*'yc* = X*'y* - (X*'i)(i'y*)/n
# e entao RSS = ||yc*||^2 - beta' Xc*'yc*.
# CENTRAGEM = "sim" (padrao) ou "nao". Com "nao" o script reproduz a forma da
# pagina 248 de LeSage e Parent (2007), que estima alfa por media(y*) e beta por
# (X*'X*)^-1 X*'y*. Isso e exato quando X*'i = 0, o que vale no SAR, onde X e
# fixa e estudentizada, e nao vale no SEM, onde X*(lambda) = X - lambda WX tem
# medias de coluna diferentes de zero. A opcao existe para medir o tamanho da
# diferenca, nao para ser usada.
CENTRAGEM <- Sys.getenv("CENTRAGEM", "sim")

rss_lambda <- function(idx) {
  k <- length(idx)
  if (k == 0L) return(tss_s)
  if (CENTRAGEM == "nao") return(rss_lambda_sem_centrar(idx))
  A00i <- A00[idx, idx, drop = FALSE]; A01i <- A01[idx, idx, drop = FALSE]
  A11i <- A11[idx, idx, drop = FALSE]
  b0i <- b0[idx]; b1i <- b1[idx]; c0i <- c0[idx]; c1i <- c1[idx]
  s0i <- s0[idx]; s1i <- s1[idx]
  sim <- A01i + t(A01i); soma_bc <- b1i + c0i
  out <- numeric(ng)
  for (j in seq_len(ng)) {
    Lj <- lambda[j]
    Xti <- s0i - Lj * s1i                                    # X*'iota
    XtX <- A00i - Lj * sim + Lj^2 * A11i - tcrossprod(Xti) / n
    Xty <- b0i - Lj * soma_bc + Lj^2 * c1i - Xti * (soma_ys[j] / n)
    beta <- tryCatch(solve(XtX, Xty), error = function(e) NULL)
    if (is.null(beta)) { out[j] <- NA_real_; next }
    out[j] <- tss_s[j] - sum(beta * Xty)
  }
  out
}

# A forma da pagina 248, reproduzida para comparacao.
rss_lambda_sem_centrar <- function(idx) {
  A00i <- A00[idx, idx, drop = FALSE]; A01i <- A01[idx, idx, drop = FALSE]
  A11i <- A11[idx, idx, drop = FALSE]
  b0i <- b0[idx]; b1i <- b1[idx]; c0i <- c0[idx]; c1i <- c1[idx]
  s0i <- s0[idx]; s1i <- s1[idx]
  sim <- A01i + t(A01i); soma_bc <- b1i + c0i
  out <- numeric(ng)
  for (j in seq_len(ng)) {
    Lj <- lambda[j]
    XtX <- A00i - Lj * sim + Lj^2 * A11i
    Xty <- b0i - Lj * soma_bc + Lj^2 * c1i
    Xti <- s0i - Lj * s1i
    beta <- tryCatch(solve(XtX, Xty), error = function(e) NULL)
    if (is.null(beta)) { out[j] <- NA_real_; next }
    out[j] <- ys_ys[j] - 2 * alfa[j] * soma_ys[j] + n * alfa[j]^2 -
      2 * sum(beta * Xty) + 2 * alfa[j] * sum(beta * Xti) +
      as.numeric(crossprod(beta, XtX %*% beta))
  }
  out
}

lm_sem <- function(idx) {
  k <- length(idx)
  rss <- rss_lambda(idx)
  if (anyNA(rss) || any(rss <= 0) || any(tss_s <= 0)) return(-Inf)
  lf <- ld_mais_prior + (1 - n) / 2 * log(tss_s)
  if (k > 0L) lf <- lf + log(f21a - 2) - log(k + f21a - 2) +
      log_f21(pmax(0, pmin(1 - rss / tss_s, 1 - 1e-12)), k)
  m <- max(lf); if (!is.finite(m)) return(-Inf)
  m + log(h / 3 * sum(pesos * exp(lf - m)))
}

# Verificacao em toda a grade, nao so em lambda = 0: a forma pre-computada tem
# de reproduzir o RSS de um ajuste explicito de y*(lambda) sobre X*(lambda) com
# intercepto.
cat("=== Verificacao: RSS pre-computado contra ajuste explicito ===\n")
set.seed(20200731L)
am <- replicate(12, sort(sample(K, sample(1:10, 1))), simplify = FALSE)
js <- unique(c(1L, round(ng * c(0.1, 0.3)), which.min(abs(lambda)),
               round(ng * c(0.7, 0.9)), ng))
pior <- 0
for (idx in am) {
  r <- rss_lambda(idx)
  for (j in js) {
    ys <- y - lambda[j] * as.numeric(W %*% y)
    Xs <- X - lambda[j] * (W %*% X)
    e <- lm.fit(cbind(1, Xs[, idx, drop = FALSE]), ys)$residuals
    pior <- max(pior, abs(r[j] - sum(e^2)))
  }
}
cat("  modelos:", length(am), "| pontos da grade:", length(js), "\n")
cat("  maior diferenca de RSS:", format(pior, scientific = TRUE), "\n")
if (CENTRAGEM == "nao") {
  cat("  CENTRAGEM = nao: a forma da pagina 248 nao minimiza a soma de\n")
  cat("  quadrados com intercepto, e a diferenca acima e o objeto do exercicio.\n\n")
} else {
  if (pior > 1e-6) stop("verificacao do RSS reprovada")
  cat("  veredito: APROVADO\n\n")
}

# ------------------------------------------------------------ amostrador
lpm <- vapply(0:K, function(k) lgamma(1 + k) + lgamma(1 + K - k), numeric(1))
tipos_disp <- function(k) { t <- character(0)
  if (k < K) t <- c(t, "birth"); if (k > 0L) t <- c(t, "death")
  if (k > 0L && k < K) t <- c(t, "move"); t }
log_q <- function(k, tipo) -log(length(tipos_disp(k))) -
  log(switch(tipo, birth = K - k, death = k, move = k * (K - k)))

mc3 <- function(semente) {
  set.seed(semente)
  atual <- sort(sample(seq_len(K), sample(0:K, 1)))
  post <- lm_sem(atual) + lpm[length(atual) + 1]
  inc <- numeric(K); g <- 0L; sk <- 0; ace <- 0L
  for (it in seq_len(N_MC)) {
    k <- length(atual); disp <- tipos_disp(k)
    tp <- disp[sample.int(length(disp), 1L)]
    fora <- setdiff(seq_len(K), atual)
    prop <- switch(tp,
                   birth = sort(c(atual, fora[sample.int(length(fora), 1L)])),
                   death = atual[-sample.int(k, 1L)],
                   move  = sort(c(atual[-sample.int(k, 1L)],
                                  fora[sample.int(length(fora), 1L)])))
    pp <- lm_sem(prop) + lpm[length(prop) + 1]
    inv <- switch(tp, birth = "death", death = "birth", move = "move")
    la <- pp - post + log_q(length(prop), inv) - log_q(k, tp)
    if (is.finite(la) && log(runif(1)) < la) {
      atual <- prop; post <- pp; ace <- ace + 1L
    }
    if (it > BURN) { g <- g + 1L
      if (length(atual) > 0L) inc[atual] <- inc[atual] + 1
      sk <- sk + length(atual) }
    if (it %% 25000L == 0L) {
      cat(sprintf("    cadeia %d: %d de %d | k = %d | %.0f s\n", semente, it,
                  N_MC, length(atual), proc.time()[3] - T0))
      flush.console()
    }
  }
  list(pip = inc / g, tam = sk / g, taxa = ace / N_MC)
}

cat("=== SEM sob a priori do TCC ===\n")
cat("  desfecho:", DESFECHO, "| vizinhanca:", VIZINHOS, "\n")
cat("  sorteios:", N_MC, "| descarte:", BURN, "\n\n")
T0 <- proc.time()[3]
a <- mc3(20200731L); b <- mc3(987654L)
cat("\n  tempo total:", round(proc.time()[3] - T0, 1), "s\n")
cat("  aceitacao:", round(a$taxa, 3), "e", round(b$taxa, 3), "\n")
cat("  tamanho medio:", round(a$tam, 2), "e", round(b$tam, 2), "\n")
cat("  maior diferenca de PIP entre cadeias:", signif(max(abs(a$pip - b$pip)), 3), "\n")
cat("  avaliacoes fora da tabela de interpolacao:", fora_tab, "\n\n")

pip <- (a$pip + b$pip) / 2
Dst <- function(p) { m <- mean(p); Kk <- length(p)
  sd(p) * sqrt((Kk - 1) / Kk) / sqrt(m * (1 - m)) }

dir.create(file.path(AQUI, "saida", "grade"), showWarnings = FALSE)
saveRDS(list(desfecho = DESFECHO, vizinhos = VIZINHOS, covariaveis = covariaveis,
             a = a, b = b, pip = pip, tam = (a$tam + b$tam) / 2, D = Dst(pip)),
        file.path(AQUI, "saida", "grade",
                  sprintf("sem_%s_%s.rds", DESFECHO, VIZINHOS)))

sar <- readRDS(file.path(AQUI, "saida", "grade",
                         sprintf("%s_%s.rds", DESFECHO, VIZINHOS)))
pip_mqo <- sar$res$MQO$pip; pip_sar <- sar$res$SAR$pip
tab <- data.frame(covariavel = covariaveis, mqo = pip_mqo,
                  sar = pip_sar, sem = pip, stringsAsFactors = FALSE)
tab <- tab[order(-tab$mqo), ]
print(transform(tab, mqo = round(mqo, 4), sar = round(sar, 4), sem = round(sem, 4)),
      row.names = FALSE)
cat("\n  D: MQO", round(Dst(pip_mqo), 4), "| SAR", round(Dst(pip_sar), 4),
    "| SEM", round(Dst(pip), 4), "\n")
cat("  tamanho medio: MQO", round(sar$res$MQO$tam, 2), "| SAR",
    round(sar$res$SAR$tam, 2), "| SEM", round((a$tam + b$tam) / 2, 2), "\n")
cat("  acima de 0,50: MQO", sum(pip_mqo >= 0.5), "| SAR", sum(pip_sar >= 0.5),
    "| SEM", sum(pip >= 0.5), "\n")
cat("  correlacao de Spearman SAR-SEM:",
    round(cor(pip_sar, pip, method = "spearman"), 4), "\n")
cat("  correlacao de Spearman MQO-SEM:",
    round(cor(pip_mqo, pip, method = "spearman"), 4), "\n")
write.csv(tab, file.path(AQUI, "saida",
                         sprintf("sem_hyperg_%s_%s.csv", DESFECHO, VIZINHOS)),
          row.names = FALSE)
cat("\nGravado em saida/\n")
