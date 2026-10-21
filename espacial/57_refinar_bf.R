# Refina o fator de Bayes agregado SAR contra MQO por Monte Carlo da priori de
# modelos, com mais modelos sorteados.
#
# Em óbitos, o tamanho efetivo dos pesos do MQO ficou entre 82 e 121 com 20.000
# modelos na especificação principal e em 66 com os dias desde o primeiro
# registro. Aqui o 43 roda com N_MODELOS maior e sementes diferentes, para medir
# a estabilidade da estimativa. Não sobrescreve os CSVs citados no artigo.
#
#   ESPEC = principal | tempo      DESFECHO = obitos (padrão)
#   N_MODELOS (padrão 200000)      SEMENTE_MODELOS (padrão 20260910)
#   SAIDA_REFINO muda a pasta de saída (testes rápidos)
#
# Saída: saida/refino_bf/<espec>_<desfecho>_n<N>_s<semente>.csv

source(if (file.exists("../caminhos.R")) "../caminhos.R" else "caminhos.R")

AQUI <- ESPACIAL
ESPEC <- Sys.getenv("ESPEC", "principal")
stopifnot(ESPEC %in% c("principal", "tempo"))
if (Sys.getenv("DESFECHO") == "") Sys.setenv(DESFECHO = "obitos")
if (Sys.getenv("VIZINHOS") == "") Sys.setenv(VIZINHOS = "queen")
if (Sys.getenv("N_MODELOS") == "") Sys.setenv(N_MODELOS = "200000")

trocar1 <- function(txt, a, b) {
  m <- gregexpr(a, txt, fixed = TRUE)[[1]]
  if (m[1] == -1L || length(m) != 1L) stop("esperava 1 ocorrência de: ", a)
  sub(a, b, txt, fixed = TRUE)
}

src <- paste(readLines(file.path(AQUI, "43_sensibilidade_prior_rho.R"), encoding = "UTF-8",
                       warn = FALSE), collapse = "\n")
if (ESPEC == "tempo") {
  source(file.path(AQUI, "tempo_exposicao_comum.R"))
  SRC33_TEMPO <- fonte_tempo("33_probabilidade_espacial.R")
  src <- trocar1(src, "expr <- parse(arq33)", "expr <- parse(text = SRC33_TEMPO)")
}
DEST_REFINO <- Sys.getenv("SAIDA_REFINO", file.path(AQUI, "saida", "refino_bf"))
dir.create(DEST_REFINO, showWarnings = FALSE, recursive = TRUE)
src <- trocar1(src, 'arq_saida <- file.path(AQUI, "saida",', sprintf('arq_saida <- file.path("%s",', DEST_REFINO))
src <- trocar1(src, 'sprintf("sensibilidade_prior_rho_%s.csv", DESFECHO)',
               sprintf('sprintf("%s_%%s_n%%d_s%%d.csv", DESFECHO, N_MODELOS, SEMENTE_MODELOS)', ESPEC))

cat(sprintf("=== refino do fator de Bayes agregado: %s, %s, N_MODELOS = %s ===\n",
            ESPEC, Sys.getenv("DESFECHO"), Sys.getenv("N_MODELOS")))
t0 <- proc.time()[3]
eval(parse(text = src), envir = globalenv())
cat(sprintf("tempo: %.1f s\n", proc.time()[3] - t0))
