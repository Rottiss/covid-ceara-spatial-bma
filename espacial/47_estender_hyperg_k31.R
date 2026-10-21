# Estende os auxiliares da hipergeométrica de K = 30 para K = 31.
#
# A sensibilidade de tempo de exposição do TCC acrescenta uma covariável, os dias
# desde o primeiro registro municipal. O vetor L do 11_hyperg_sar.R e a tabela do
# 15_otimizar_hyperg.R cobrem k de 0 a 30; com 31 covariáveis, o modelo cheio
# leria fora da tabela sem erro nenhum.
#
# Mesmas fórmulas, mesmo n = 184 (f21a e a_par não mudam), mesma grade de z.
# As colunas de 0 a 30 são copiadas e conferidas como idênticas; a nova é
# conferida contra a série exata. Grava arquivos novos e não toca nos antigos,
# que sustentam as saídas já verificadas.

source(if (file.exists("../caminhos.R")) "../caminhos.R" else "caminhos.R")

AQUI <- ESPACIAL
K_NOVO <- as.integer(Sys.getenv("K_NOVO", "31"))

hyp <- readRDS(file.path(AQUI, "saida", "hyperg_aux.rds"))
tb  <- readRDS(file.path(AQUI, "saida", "tabela_f21_estendida.rds"))
f21a <- hyp$f21a; a_par <- hyp$a_par; J_MAX <- hyp$J_MAX
K_ANT <- length(hyp$L) - 1L
TAB <- matrix(tb$TAB, nrow = tb$N_TAB)
stopifnot(ncol(TAB) == K_ANT + 1L, K_NOVO > K_ANT)
cat("auxiliares atuais: k de 0 a", K_ANT, "| tabela", tb$N_TAB, "pontos até z =", tb$Z_MAX, "\n")

# Mesma fórmula do 11_hyperg_sar.R.
L <- hyp$L
for (k in (K_ANT + 1L):K_NOVO) {
  cc <- (f21a + k) / 2
  L[[k + 1L]] <- cumsum(log(a_par + 0:(J_MAX - 1)) - log(cc + 0:(J_MAX - 1)))
}
L_ref <- lapply(0:K_ANT, function(k) {
  cc <- (f21a + k) / 2
  cumsum(log(a_par + 0:(J_MAX - 1)) - log(cc + 0:(J_MAX - 1)))
})
stopifnot(isTRUE(all.equal(L_ref, hyp$L, tolerance = 0)))
cat("L: fórmula reproduz as", K_ANT + 1L, "entradas gravadas\n")

# Mesma série do 15_otimizar_hyperg.R.
log_f21_serie <- function(z, k) {
  cc <- (f21a + k) / 2
  vapply(z, function(zi) {
    j_pico <- max(0, ceiling((zi * a_par - cc) / max(1 - zi, 1e-12)))
    J <- min(J_MAX, max(300L, as.integer(j_pico * 4)))
    lt <- c(0, log(zi) * seq_len(J) + L[[k + 1]][seq_len(J)])
    m <- max(lt)
    m + log(sum(exp(lt - m)))
  }, numeric(1))
}
set.seed(45)
z_conf <- runif(200, 0, tb$Z_MAX)
dif_ant <- max(abs(TAB[, K_ANT + 1L][as.integer(z_conf / tb$passo) + 1L] -
                   log_f21_serie(tb$z_tab[as.integer(z_conf / tb$passo) + 1L], K_ANT)))
cat("série contra coluna gravada k =", K_ANT, "nos pontos da grade:", format(dif_ant, digits = 3), "\n")
stopifnot(dif_ant < 1e-10)

t0 <- proc.time()[3]
novas <- vapply((K_ANT + 1L):K_NOVO, function(k) log_f21_serie(tb$z_tab, k), numeric(tb$N_TAB))
cat("colunas novas em", round(proc.time()[3] - t0, 1), "s\n")
TAB2 <- cbind(TAB, novas)

log_f21_tab <- function(z, k) {
  pos <- z / tb$passo; i <- as.integer(pos) + 1L
  i[i >= tb$N_TAB] <- tb$N_TAB - 1L
  frac <- pos - (i - 1L); col <- TAB2[, k + 1L]
  col[i] * (1 - frac) + col[i + 1L] * frac
}
set.seed(1)
for (k in (K_ANT + 1L):K_NOVO) {
  z_teste <- runif(400, 0, tb$Z_MAX)
  e <- max(abs(log_f21_tab(z_teste, k) - log_f21_serie(z_teste, k)))
  cat(sprintf("k = %d: maior |erro| de interpolação em log %.3e\n", k, e))
  if (e > 1e-7) stop("erro de interpolação alto demais")
}
stopifnot(identical(TAB2[, seq_len(K_ANT + 1L)], TAB))

saveRDS(list(f21a = f21a, a_par = a_par, L = L, J_MAX = J_MAX),
        file.path(AQUI, "saida", sprintf("hyperg_aux_k%d.rds", K_NOVO)))
saveRDS(list(TAB = as.numeric(TAB2), z_tab = tb$z_tab, passo = tb$passo,
             N_TAB = tb$N_TAB, Z_MAX = tb$Z_MAX),
        file.path(AQUI, "saida", sprintf("tabela_f21_k%d.rds", K_NOVO)))
cat("gravados: hyperg_aux_k", K_NOVO, ".rds e tabela_f21_k", K_NOVO, ".rds\n", sep = "")
