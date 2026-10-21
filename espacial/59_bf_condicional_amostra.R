# Recalcula e grava os fatores de Bayes condicionais SAR contra MQO em 200
# conjuntos de covariáveis sorteados, que o 33 só imprime na tela quando a cadeia
# não visita o MQO (casos e óbitos). O relatório e o artigo citam a mediana, o
# mínimo e o fator condicional mediano; sem arquivo, o verificador não os confere.
#
# Reaproveita as definições do 33 (dados, W, log-marginais) sem rodar as cadeias.
# O sorteio usa set.seed(1L), o mesmo do 33, e independe das cadeias.
#
# Saída: saida/bf_condicional_amostra.csv

source(if (file.exists("../caminhos.R")) "../caminhos.R" else "caminhos.R")

AQUI <- ESPACIAL
src <- readLines(file.path(AQUI, "33_probabilidade_espacial.R"), encoding = "UTF-8", warn = FALSE)
corte <- grep('cat("=== Especificacao dentro do espaco amostral ===\\n")', src, fixed = TRUE)
stopifnot(length(corte) == 1L)
defs <- parse(text = src[seq_len(corte - 1L)])

VIZ <- Sys.getenv("VIZINHOS", "queen")
res <- list()
for (d in c("casos", "obitos")) {
  Sys.setenv(DESFECHO = d, VIZINHOS = VIZ)
  env <- new.env(parent = globalenv())
  eval(defs, envir = env)
  local({
    set.seed(1L)
    amostra <- replicate(200, sort(sample(K, sample(0:K, 1))), simplify = FALSE)
    dif <- vapply(amostra, function(idx) lm_sar(estat(idx)) - lm_ols(estat(idx)), numeric(1))
    res[[d]] <<- data.frame(desfecho = d, vizinhos = VIZ, n_conjuntos = length(dif),
                            log_bf_mediana = median(dif), log_bf_minimo = min(dif),
                            log_bf_maximo = max(dif), bf_mediano = exp(median(dif)),
                            fora_tabela = fora_tab)
  }, envir = env)
}
tab <- do.call(rbind, res)
arq <- file.path(AQUI, "saida", "bf_condicional_amostra.csv")
write.csv(tab, arq, row.names = FALSE)
print(tab, digits = 4, row.names = FALSE)
cat("Gravado em", arq, "\n")
