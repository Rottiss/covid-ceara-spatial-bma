# MC3 over covariate subsets and OLS (MQO), SAR, and SEM specifications.
# Scripts 72 and 73 validate marginal-likelihood scaling and transitions
# against enumeration. This script uses the same functions and repeats the
# scale check on five randomly selected models before sampling.
#
# Priors: hyper-g with UIP calibration; binomial-beta with expected size K/2;
# a0=1.01 for rho and lambda; uniform prior over the three specifications.
# Integration uses 201 rho points and 101 lambda points.
#
# Marginal likelihoods are cached by specification. Chain A starts in SAR,
# chain B in SEM. Posterior rho/lambda distributions are not accumulated here;
# separate estimation scripts 33 and 31 provide those outputs.
#
# Writes to saida/tres_vias/ or SAIDA_TRES; existing outputs are not overwritten.
#
#   DESFECHO     = casos | obitos | letalidade
#   VIZINHOS     = queen | rook
#   MC3_SORTEIOS   default 1,000,000
#   MC3_BURN       default 200,000
#   SAIDA_TRES     default saida/tres_vias

source(if (file.exists("../caminhos.R")) "../caminhos.R" else "caminhos.R")

suppressPackageStartupMessages(library(spdep))
AQUI <- ESPACIAL
RAIZ <- ESTIMACOES
DESFECHO <- Sys.getenv("DESFECHO", "casos")
VIZINHOS <- Sys.getenv("VIZINHOS", "queen")
N_MC <- as.integer(Sys.getenv("MC3_SORTEIOS", "1000000"))
BURN <- as.integer(Sys.getenv("MC3_BURN", "200000"))
DEST <- Sys.getenv("SAIDA_TRES", file.path(AQUI, "saida", "tres_vias"))
stopifnot(BURN < N_MC)
stopifnot(DESFECHO %in% c("casos", "obitos", "letalidade"),
          VIZINHOS %in% c("queen", "rook"),
          is.finite(N_MC), is.finite(BURN), N_MC > 0, BURN >= 0)
arq <- file.path(DEST, sprintf("%s_%s.rds", DESFECHO, VIZINHOS))
if (file.exists(arq)) stop("Output exists; choose another SAIDA_TRES directory: ", arq)

hyp <- readRDS(file.path(AQUI, "saida", "hyperg_aux.rds"))
tb  <- readRDS(file.path(AQUI, "saida", "tabela_f21_estendida.rds"))
f21a <- hyp$f21a; a_par <- hyp$a_par; L <- hyp$L; J_MAX <- hyp$J_MAX
TAB <- tb$TAB; passo <- tb$passo; N_TAB <- tb$N_TAB; Z_MAX <- tb$Z_MAX

COLUNA <- c(casos = "z_log_casos100k", obitos = "z_log_obitos100k",
            letalidade = "z_logit_letal")[[DESFECHO]]
dados <- read.csv(file.path(RAIZ, "resultados/estimacoes/matrizes_estimacao/principal",
                            paste0(DESFECHO, ".csv")),
                  stringsAsFactors = FALSE, colClasses = c(codigo_ibge = "character"))
covariaveis <- setdiff(names(dados), c("codigo_ibge", "municipio", COLUNA))
y <- as.numeric(dados[[COLUNA]])
X <- as.matrix(dados[, covariaveis]); n <- length(y); K <- ncol(X)
Xc <- sweep(X, 2, colMeans(X), "-")
listw <- readRDS(file.path(RAIZ, "resultados/diagnostico_espacial",
                           paste0("pesos_", VIZINHOS, ".rds")))
W <- spdep::listw2mat(listw); dimnames(W) <- NULL
stopifnot(nrow(W) == n, max(abs(rowSums(W) - 1)) < 1e-9, max(abs(diag(W))) < 1e-12)

# --------------------------------------------------------- Hypergeometric function
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
  list(v = v, ldp = ld + lpr, lpr = lpr, j0 = which(abs(v) < 1e-12),
       pes = c(1, rep(c(4, 2), length.out = ngp - 2), 1), h = v[2] - v[1])
}
gr <- grade(201L); gl <- grade(101L)
simpson <- function(lf, g) {
  m <- max(lf); if (!is.finite(m)) return(-Inf)
  m + log(g$h / 3 * sum(g$pes * exp(lf - m)))
}

# --------------------------------------------------------- OLS (MQO) and SAR
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
lm_ols <- function(idx) {
  e <- estat(idx); k <- unname(e["k"]); nl <- (1 - n) / 2 * log(ta)
  if (k == 0L) return(nl)
  unname(nl + log(f21a - 2) - log(k + f21a - 2) + log_f21(1 - e["a"] / ta, k))
}
lf_sar <- function(idx) {
  e <- estat(idx); k <- unname(e["k"])
  rss <- e["a"] - 2 * gr$v * e["b"] + gr$v^2 * e["c"]
  tss <- ta - 2 * gr$v * tbb + gr$v^2 * tcc
  if (any(rss <= 0) || any(tss <= 0)) return(rep(-Inf, length(gr$v)))
  lf <- gr$ldp + (1 - n) / 2 * log(tss)
  if (k > 0L) lf <- lf + log(f21a - 2) - log(k + f21a - 2) +
      log_f21(pmax(0, 1 - rss / tss), k)
  unname(lf)
}

# --------------------------------------------------------- SEM with centering
WX <- W %*% X
A00 <- crossprod(X); A01 <- crossprod(X, WX); A11 <- crossprod(WX)
b0 <- as.numeric(crossprod(X, y));  b1 <- as.numeric(crossprod(X, Wy))
c0 <- as.numeric(crossprod(WX, y)); c1 <- as.numeric(crossprod(WX, Wy))
s0 <- colSums(X); s1 <- colSums(WX)
yy <- sum(y^2); yWy <- sum(y * Wy); WyWy <- sum(Wy^2); sy <- sum(y); sWy <- sum(Wy)
ys_ys <- yy - 2 * gl$v * yWy + gl$v^2 * WyWy
soma_ys <- sy - gl$v * sWy
tss_s <- ys_ys - n * (soma_ys / n)^2
lf_sem <- function(idx) {
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
  if (anyNA(rss) || any(rss <= 0)) return(rep(-Inf, length(gl$v)))
  lf <- gl$ldp + (1 - n) / 2 * log(tss_s)
  if (k > 0L) lf <- lf + log(f21a - 2) - log(k + f21a - 2) +
      log_f21(pmax(0, pmin(1 - rss / tss_s, 1 - 1e-12)), k)
  lf
}

# --------------------------------------------------------- Scale check
set.seed(1L)
for (t in 1:5) {
  idx <- sort(sample(seq_len(K), sample(0:10, 1)))
  ref <- lm_ols(idx)
  d1 <- lf_sar(idx)[gr$j0] - gr$lpr[gr$j0] - ref
  d2 <- lf_sem(idx)[gl$j0] - gl$lpr[gl$j0] - ref
  if (!(abs(d1) < 1e-8 && abs(d2) < 1e-8))
    stop(sprintf("Common-scale check failed: k=%d, SAR %.3e, SEM %.3e", length(idx), d1, d2))
}

# --------------------------------------------------------- Cached marginal likelihoods
ESP <- c("MQO", "SAR", "SEM")
cache <- list(MQO = new.env(hash = TRUE), SAR = new.env(hash = TRUE),
              SEM = new.env(hash = TRUE))
marg <- function(idx, esp) {
  chave <- if (length(idx)) paste(idx, collapse = ",") else "vazio"
  v <- cache[[esp]][[chave]]
  if (!is.null(v)) return(v)
  v <- switch(esp, MQO = lm_ols(idx), SAR = simpson(lf_sar(idx), gr),
              SEM = simpson(lf_sem(idx), gl))
  assign(chave, v, envir = cache[[esp]])
  v
}
lpm <- vapply(0:K, function(k) lgamma(1 + k) + lgamma(1 + K - k), numeric(1))

tipos_disp <- function(k) {
  t <- "troca"
  if (k < K) t <- c(t, "birth"); if (k > 0L) t <- c(t, "death")
  if (k > 0L && k < K) t <- c(t, "move"); t
}
# The switching term cancels in the ratio (see 73), but is explicit here.
log_q <- function(k, tipo)
  -log(length(tipos_disp(k))) -
  log(switch(tipo, troca = length(ESP) - 1, birth = K - k, death = k, move = k * (K - k)))

T0 <- proc.time()[3]
mc3 <- function(semente, esp_inicial) {
  set.seed(semente)
  atual <- sort(sample(seq_len(K), sample(0:K, 1))); esp <- esp_inicial
  post <- marg(atual, esp) + lpm[length(atual) + 1]
  inc <- matrix(0, K, 3, dimnames = list(covariaveis, ESP))
  ce <- setNames(numeric(3), ESP); tam <- setNames(numeric(3), ESP)
  prop_tr <- 0L; ace_tr <- 0L; ace <- 0L; g <- 0L
  traco <- integer(0)
  for (it in seq_len(N_MC)) {
    k <- length(atual); disp <- tipos_disp(k)
    tp <- disp[sample.int(length(disp), 1L)]
    if (tp == "troca") {
      prop <- atual; alt <- setdiff(ESP, esp)
      esp_p <- alt[sample.int(2L, 1L)]; prop_tr <- prop_tr + 1L
    } else {
      esp_p <- esp; fora <- setdiff(seq_len(K), atual)
      prop <- switch(tp,
                     birth = sort(c(atual, fora[sample.int(length(fora), 1L)])),
                     death = atual[-sample.int(k, 1L)],
                     move  = sort(c(atual[-sample.int(k, 1L)],
                                    fora[sample.int(length(fora), 1L)])))
    }
    pp <- marg(prop, esp_p) + lpm[length(prop) + 1]
    inv <- switch(tp, troca = "troca", birth = "death", death = "birth", move = "move")
    la <- pp - post + log_q(length(prop), inv) - log_q(k, tp)
    if (is.finite(la) && log(runif(1)) < la) {
      atual <- prop; esp <- esp_p; post <- pp; ace <- ace + 1L
      if (tp == "troca") ace_tr <- ace_tr + 1L
    }
    if (it > BURN) {
      g <- g + 1L; ce[esp] <- ce[esp] + 1; tam[esp] <- tam[esp] + length(atual)
      if (length(atual)) inc[atual, esp] <- inc[atual, esp] + 1
    }
    if (it %% 1000L == 0L) traco <- c(traco, match(esp, ESP))
    if (it %% 25000L == 0L) {
      cat(sprintf("    chain %d: %d of %d | %s | k = %d | cache %d/%d/%d | %.0f s\n",
                  semente, it, N_MC, esp, length(atual),
                  length(ls(cache$MQO)), length(ls(cache$SAR)), length(ls(cache$SEM)),
                  proc.time()[3] - T0))
      flush.console()
    }
  }
  list(inc = inc, ce = ce, tam = tam, g = g, taxa = ace / N_MC,
       taxa_troca = if (prop_tr) ace_tr / prop_tr else NA, traco = traco)
}

cat("=== MC3 with OLS (MQO), SAR, and SEM ===\n")
cat(sprintf("  outcome: %s | contiguity: %s | K = %d | iterations: %d | burn-in: %d\n\n",
            DESFECHO, VIZINHOS, K, N_MC, BURN))
a <- mc3(20200731L, "SAR")
b <- mc3(987654L, "SEM")

pa <- a$ce / a$g; pb <- b$ce / b$g; p <- (pa + pb) / 2
cat("\n=== posterior specification probabilities ===\n")
for (e in ESP)
  cat(sprintf("  P(%s|y)  chain A %.5f | chain B %.5f | mean %.5f\n", e, pa[e], pb[e], p[e]))
cat(sprintf("  maximum difference between chains: %.5f\n", max(abs(pa - pb))))
cat(sprintf("  overall acceptance %.3f/%.3f | switch acceptance %.3f/%.3f\n",
            a$taxa, b$taxa, a$taxa_troca, b$taxa_troca))

# Conditional PIPs by specification, preserving the original covariate order.
pip <- matrix(NA_real_, K, 3, dimnames = list(covariaveis, ESP))
for (e in ESP) {
  n_e <- a$ce[e] + b$ce[e]
  if (n_e > 0) pip[, e] <- (a$inc[, e] + b$inc[, e]) / n_e
}
tabela <- data.frame(covariavel = covariaveis, pip, check.names = FALSE,
                     stringsAsFactors = FALSE)
stopifnot(identical(tabela$covariavel, rownames(a$inc)),
          identical(rownames(a$inc), rownames(b$inc)))
cat("\n=== conditional PIPs, visited specifications ===\n")
for (e in ESP) {
  if (all(is.na(pip[, e]))) next
  sel <- covariaveis[!is.na(pip[, e]) & pip[, e] >= 0.5]
  cat(sprintf("  %s: above 0.50 = %s | mean size %.2f\n", e,
              if (length(sel)) paste(sel, collapse = ", ") else "nenhuma",
              (a$tam[e] + b$tam[e]) / max(1, a$ce[e] + b$ce[e])))
}
cat(sprintf("\n  evaluations outside the f21 table: %d\n", fora_tab))
cat(sprintf("  total time: %.1f s\n", proc.time()[3] - T0))

dir.create(DEST, showWarnings = FALSE, recursive = TRUE)
arq <- file.path(DEST, sprintf("%s_%s.rds", DESFECHO, VIZINHOS))
stopifnot(!file.exists(arq))
saveRDS(list(desfecho = DESFECHO, vizinhos = VIZINHOS, n_mc = N_MC, burn = BURN,
             p_esp = p, p_esp_a = pa, p_esp_b = pb, tabela = tabela,
             a = a, b = b, fora_tab = fora_tab), arq)
cat("Saved to", arq, "\n")
