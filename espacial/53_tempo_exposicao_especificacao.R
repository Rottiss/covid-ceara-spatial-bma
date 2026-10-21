# Especificação dentro do espaço amostral na sensibilidade de tempo de exposição.
#
# Roda o 33_probabilidade_espacial.R com a matriz de 31 covariáveis (dias desde o
# primeiro registro municipal) e a mesma configuração usada nas estimações
# principais da Tabela 4: duas cadeias de 300.000 iterações, 60.000 de burn-in.
# Dá a posteriori de rho sob BMA e a frequência do SAR contra o MQO, para saber se
# o rho de casos cai quando os dias entram no modelo.
#
#   DESFECHO = casos | obitos | letalidade    VIZINHOS = queen | rook
#
# Saída: saida/tempo_exposicao/especificacao_<desfecho>_<vizinhos>.rds

source(if (file.exists("../caminhos.R")) "../caminhos.R" else "caminhos.R")

source(file.path(ESPACIAL, "tempo_exposicao_comum.R"))

src <- fonte_tempo("33_probabilidade_espacial.R")
src <- trocar1(src, 'file.path(AQUI, "saida", sprintf("especificacao_%s_%s.rds", DESFECHO, VIZINHOS))',
               sprintf('file.path("%s", sprintf("especificacao_%%s_%%s.rds", DESFECHO, VIZINHOS))', DEST))
stopifnot(!grepl('file.path(AQUI, "saida", sprintf("especificacao', src, fixed = TRUE))

padrao_env("DESFECHO", "casos")
padrao_env("VIZINHOS", "queen")
padrao_env("MC3_SORTEIOS", "300000")
padrao_env("MC3_BURN", "60000")
cat("=== tempo de exposição: especificação dentro do espaço amostral ===\n")
eval(parse(text = src), envir = globalenv())
