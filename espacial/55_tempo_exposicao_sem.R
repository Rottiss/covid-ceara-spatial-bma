# SEM na sensibilidade de tempo de exposição.
#
# Roda o 31_sem_hyperg.R com a matriz de 31 covariáveis, na mesma configuração
# das estimações principais: duas cadeias de 1.000.000 de iterações, 200.000 de
# burn-in, grade de 101 pontos em lambda. A comparação final do 31 lê o MQO e o
# SAR da mesma pasta, que é onde o 48 gravou as cadeias com os dias.
#
#   DESFECHO = casos | obitos | letalidade    VIZINHOS = queen | rook
#
# Saída: saida/tempo_exposicao/sem_<desfecho>_<vizinhos>.rds e sem_hyperg_*.csv.
# Exige o 48 para a mesma combinação. Cerca de 3 horas por combinação.

source(if (file.exists("../caminhos.R")) "../caminhos.R" else "caminhos.R")

source(file.path(ESPACIAL, "tempo_exposicao_comum.R"))

trocar_n <- function(txt, a, b, n) {
  m <- gregexpr(a, txt, fixed = TRUE)[[1]]
  achou <- if (m[1] == -1L) 0L else length(m)
  if (achou != n) stop(sprintf("esperava %d ocorrência(s), achei %d: %s", n, achou, a))
  gsub(a, b, txt, fixed = TRUE)
}

arq_sar <- file.path(DEST, sprintf("%s_%s.rds", Sys.getenv("DESFECHO", "casos"), Sys.getenv("VIZINHOS", "queen")))
if (!file.exists(arq_sar)) stop("falta o SAR desta combinação (rodar o 48): ", arq_sar)

src <- fonte_tempo("31_sem_hyperg.R")
src <- trocar_n(src, 'dir.create(file.path(AQUI, "saida", "grade"), showWarnings = FALSE)',
                sprintf('dir.create("%s", showWarnings = FALSE)', DEST), 1L)
# Duas ocorrências: a gravação do SEM e a leitura do MQO/SAR para a comparação.
src <- trocar_n(src, 'file.path(AQUI, "saida", "grade",', sprintf('file.path("%s",', DEST), 2L)
src <- trocar_n(src, 'write.csv(tab, file.path(AQUI, "saida",', sprintf('write.csv(tab, file.path("%s",', DEST), 1L)
stopifnot(!grepl('"saida", "grade"', src, fixed = TRUE))

padrao_env("DESFECHO", "casos")
padrao_env("VIZINHOS", "queen")
cat("=== tempo de exposição: SEM ===\n")
eval(parse(text = src), envir = globalenv())
