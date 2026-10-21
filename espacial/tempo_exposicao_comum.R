# Trocas de texto comuns aos scripts que rodam a sensibilidade de tempo de
# exposição (dias desde o primeiro registro municipal, 31 covariáveis).
#
# Nenhum script original é alterado: o texto é lido, trocado e avaliado. Cada
# troca exige exatamente uma ocorrência, para que uma mudança futura no script de
# origem faça a execução parar em vez de rodar a especificação errada.
#
# SAIDA_TEMPO muda a pasta de saída; serve para testes rápidos fora de saida/.

source(if (file.exists("../caminhos.R")) "../caminhos.R" else "caminhos.R")

AQUI <- ESPACIAL
DEST <- Sys.getenv("SAIDA_TEMPO", file.path(AQUI, "saida", "tempo_exposicao"))
dir.create(DEST, showWarnings = FALSE, recursive = TRUE)
for (f in c("hyperg_aux_k31.rds", "tabela_f21_k31.rds"))
  if (!file.exists(file.path(AQUI, "saida", f))) stop("rodar antes o 47_estender_hyperg_k31.R: falta ", f)

trocar1 <- function(txt, a, b) {
  m <- gregexpr(a, txt, fixed = TRUE)[[1]]
  if (m[1] == -1L || length(m) != 1L) stop("esperava 1 ocorrência de: ", a)
  sub(a, b, txt, fixed = TRUE)
}

# Lê um script do projeto já apontado para a matriz de tempo de exposição e para
# os auxiliares da hipergeométrica com K = 31.
fonte_tempo <- function(arquivo) {
  src <- paste(readLines(file.path(AQUI, arquivo), encoding = "UTF-8", warn = FALSE),
               collapse = "\n")
  src <- trocar1(src, '"hyperg_aux.rds"', '"hyperg_aux_k31.rds"')
  src <- trocar1(src, '"tabela_f21_estendida.rds"', '"tabela_f21_k31.rds"')
  src <- trocar1(src, '"resultados/estimacoes/matrizes_estimacao/principal"',
                 '"resultados/estimacoes/matrizes_estimacao/sensibilidade_tempo_exposicao"')
  src <- trocar1(src, "n <- length(y); K <- length(covariaveis)",
                 paste("n <- length(y); K <- length(covariaveis)",
                       'stopifnot(K == 31L, "z_dias.primeiro.registro" %in% covariaveis)',
                       sep = "\n"))
  src
}

padrao_env <- function(nome, valor) if (Sys.getenv(nome) == "") do.call(Sys.setenv, setNames(list(valor), nome))
