# Auditoria da escolha MQO/SAR e da sensibilidade a priori de rho.
#
# Distingue duas quantidades:
# 1. fator de Bayes condicional ao conjunto M, p(y|M,SAR)/p(y|M,MQO);
# 2. fator de Bayes agregado, soma_M p(y|M,SAR)P(M) dividida pela soma analoga.
#
# A segunda e estimada por Monte Carlo direto da priori de modelos. O ESS dos
# pesos de evidencia e reportado porque a estimativa pode ser instavel quando a
# evidencia se concentra em poucos modelos.

source(if (file.exists("../caminhos.R")) "../caminhos.R" else "caminhos.R")

AQUI <- ESPACIAL
arq33 <- file.path(AQUI, "33_probabilidade_espacial.R")

# Carrega dados e funcoes do script 33 sem executar o amostrador longo.
expr <- parse(arq33)
for (e in expr) {
  txt <- paste(deparse(e, width.cutoff = 500L), collapse = " ")
  if (grepl("^mc3 <- function", txt)) break
  eval(e, envir = .GlobalEnv)
}

logsumexp <- function(x) {
  m <- max(x)
  m + log(sum(exp(x - m)))
}
logmeanexp <- function(x) logsumexp(x) - log(length(x))
ess_logpesos <- function(x) exp(2 * logsumexp(x) - logsumexp(2 * x))

a0_grade <- c(1.01, 1.1, 1.5, 2, 5, 10)
# Padrão igual ao usado nos CSVs citados no artigo e conferidos pelo 38.
N_MODELOS <- as.integer(Sys.getenv("N_MODELOS", "20000"))
SEMENTE_MODELOS <- as.integer(Sys.getenv("SEMENTE_MODELOS", "20260910"))
set.seed(SEMENTE_MODELOS)
modelos <- replicate(N_MODELOS, {
  k <- sample(0:K, 1L)
  if (k == 0L) integer(0) else sort(sample.int(K, k))
}, simplify = FALSE)
est <- lapply(modelos, estat)
lm_mqo <- vapply(est, lm_ols, numeric(1))

saida <- lapply(a0_grade, function(a0_i) {
  lpr <<- (a0_i - 1) * log1p(rho) + (a0_i - 1) * log1p(-rho) -
    (2 * a0_i - 1) * log(2) - lbeta(a0_i, a0_i)
  ld_mais_prior <<- ld + lpr
  lm_sar_i <- vapply(est, lm_sar, numeric(1))
  dif <- lm_sar_i - lm_mqo
  bloco <- rep(seq_len(10L), length.out = N_MODELOS)
  bf_blocos <- vapply(seq_len(10L), function(b) {
    ii <- bloco == b
    logmeanexp(lm_sar_i[ii]) - logmeanexp(lm_mqo[ii])
  }, numeric(1))
  log_bf_agregado <- logmeanexp(lm_sar_i) - logmeanexp(lm_mqo)
  data.frame(
    desfecho = DESFECHO,
    a0 = a0_i,
    n_modelos = N_MODELOS,
    log_bf_cond_mediano = median(dif),
    log_bf_cond_min = min(dif),
    log_bf_cond_max = max(dif),
    log_bf_agregado_mc = log_bf_agregado,
    p_sar_agregada_mc = plogis(log_bf_agregado),
    dp_log_bf_blocos = sd(bf_blocos),
    ess_sar = ess_logpesos(lm_sar_i),
    ess_mqo = ess_logpesos(lm_mqo)
  )
})
saida <- do.call(rbind, saida)

print(transform(saida,
                log_bf_cond_mediano = round(log_bf_cond_mediano, 3),
                log_bf_cond_min = round(log_bf_cond_min, 3),
                log_bf_cond_max = round(log_bf_cond_max, 3),
                log_bf_agregado_mc = round(log_bf_agregado_mc, 3),
                ess_sar = round(ess_sar, 1), ess_mqo = round(ess_mqo, 1)),
      row.names = FALSE)

arq_saida <- file.path(AQUI, "saida",
                       sprintf("sensibilidade_prior_rho_%s.csv", DESFECHO))
write.csv(saida, arq_saida, row.names = FALSE)
cat("\nGravado em", arq_saida, "\n")
cat("O BF agregado por Monte Carlo direto so deve ser usado se os ESS forem\n")
cat("adequados e a repeticao com outra semente for estavel.\n")
