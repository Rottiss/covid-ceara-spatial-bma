# Fator de Bayes agregado e sensibilidade à priori de rho, especificação com os
# dias, com a matriz Rook.
#
# Reproduz as trocas do 54_tempo_exposicao_prior_rho.R sobre o 43, mais uma:
# o 43 grava sensibilidade_prior_rho_<desfecho>.csv, sem a vizinhança no nome,
# e o 54 só redireciona a pasta (para DEST), não o nome do arquivo. Rodar com
# VIZINHOS=rook sem essa troca extra sobrescreveria os arquivos da Queen já
# citados no artigo e conferidos pelo 38. Nem o 43 nem o 54 são alterados.
#
#   DESFECHO = casos | obitos | letalidade    VIZINHOS = rook (obrigatório)
#
# Mesma grade de a0, N_MODELOS = 20.000 e semente do 54, para ficar comparável
# ao par Queen/Rook que a Seção 4.4 já reporta na especificação principal.
#
# Saída: saida/tempo_exposicao/sensibilidade_prior_rho_<desfecho>_rook.csv

source(if (file.exists("../caminhos.R")) "../caminhos.R" else "caminhos.R")

source(file.path(ESPACIAL, "tempo_exposicao_comum.R"))
stopifnot(Sys.getenv("VIZINHOS") == "rook")

SRC33_TEMPO <- fonte_tempo("33_probabilidade_espacial.R")
src <- paste(readLines(file.path(AQUI, "43_sensibilidade_prior_rho.R"), encoding = "UTF-8",
                       warn = FALSE), collapse = "\n")
src <- trocar1(src, "expr <- parse(arq33)", "expr <- parse(text = SRC33_TEMPO)")
src <- trocar1(src, 'arq_saida <- file.path(AQUI, "saida",', sprintf('arq_saida <- file.path("%s",', DEST))
src <- trocar1(src, 'sprintf("sensibilidade_prior_rho_%s.csv", DESFECHO)',
               'sprintf("sensibilidade_prior_rho_%s_rook.csv", DESFECHO)')

padrao_env("DESFECHO", "casos")
cat(sprintf("=== tempo de exposicao, priori de rho com Rook: %s ===\n", Sys.getenv("DESFECHO")))
t0 <- proc.time()[3]
eval(parse(text = src), envir = globalenv())
cat(sprintf("tempo: %.1f s\n", proc.time()[3] - t0))
