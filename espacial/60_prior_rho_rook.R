# Fator de Bayes agregado e sensibilidade à priori de rho com a matriz Rook.
#
# O 43 grava sensibilidade_prior_rho_<desfecho>.csv, sem a vizinhança no nome.
# Rodá-lo com VIZINHOS=rook sobrescreveria os arquivos da Queen citados no artigo
# e conferidos pelo 38. Aqui o nome de saída recebe o sufixo _rook, por
# substituição de texto, sem alterar o 43.
#
#   DESFECHO = casos | obitos | letalidade     VIZINHOS = rook (obrigatório)
#
# Saída: saida/sensibilidade_prior_rho_<desfecho>_rook.csv

source(if (file.exists("../caminhos.R")) "../caminhos.R" else "caminhos.R")

AQUI <- ESPACIAL
stopifnot(Sys.getenv("VIZINHOS") == "rook")

trocar1 <- function(txt, a, b) {
  m <- gregexpr(a, txt, fixed = TRUE)[[1]]
  if (m[1] == -1L || length(m) != 1L) stop("esperava 1 ocorrência de: ", a)
  sub(a, b, txt, fixed = TRUE)
}
src <- paste(readLines(file.path(AQUI, "43_sensibilidade_prior_rho.R"),
                       encoding = "UTF-8", warn = FALSE), collapse = "\n")
src <- trocar1(src, 'sprintf("sensibilidade_prior_rho_%s.csv", DESFECHO)',
               'sprintf("sensibilidade_prior_rho_%s_rook.csv", DESFECHO)')
cat(sprintf("=== priori de rho com Rook: %s ===\n", Sys.getenv("DESFECHO")))
t0 <- proc.time()[3]
eval(parse(text = src), envir = globalenv())
cat(sprintf("tempo: %.1f s\n", proc.time()[3] - t0))
