# Fator de Bayes agregado e sensibilidade à priori de rho na sensibilidade de
# tempo de exposição.
#
# Roda o 43_sensibilidade_prior_rho.R, que carrega as funções do 33 sem executar a
# cadeia, com o 33 já apontado para a matriz de 31 covariáveis. Mesma grade de a0,
# N_MODELOS = 20.000 e semente dos CSVs citados no artigo.
#
#   DESFECHO = casos | obitos | letalidade    VIZINHOS = queen | rook
#
# Saída: saida/tempo_exposicao/sensibilidade_prior_rho_<desfecho>.csv

source(if (file.exists("../caminhos.R")) "../caminhos.R" else "caminhos.R")

source(file.path(ESPACIAL, "tempo_exposicao_comum.R"))

SRC33_TEMPO <- fonte_tempo("33_probabilidade_espacial.R")
src <- paste(readLines(file.path(AQUI, "43_sensibilidade_prior_rho.R"), encoding = "UTF-8",
                       warn = FALSE), collapse = "\n")
src <- trocar1(src, "expr <- parse(arq33)", "expr <- parse(text = SRC33_TEMPO)")
src <- trocar1(src, 'arq_saida <- file.path(AQUI, "saida",', sprintf('arq_saida <- file.path("%s",', DEST))

padrao_env("DESFECHO", "casos")
padrao_env("VIZINHOS", "queen")
cat("=== tempo de exposição: fator de Bayes agregado e priori de rho ===\n")
eval(parse(text = src), envir = globalenv())
