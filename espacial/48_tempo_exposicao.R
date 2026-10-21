# MQO e SAR na sensibilidade de tempo de exposição, com a configuração principal.
#
# No TCC, as sensibilidades rodaram com uma cadeia de 100.000 de burn-in e
# 500.000 iterações; o modelo principal, com duas cadeias de 200.000 e 1.000.000.
# Aqui a especificação com os dias desde o primeiro registro municipal (31
# covariáveis) roda como o principal: duas cadeias, 200.000 de burn-in e
# 1.000.000 de iterações, MQO e SAR com a mesma priori. A cadeia do MQO refaz a
# sensibilidade com a configuração principal; a do SAR é a versão espacial.
#
# Executa o 20_grade_empirica.R com trocas de texto, sem alterar o arquivo: a
# matriz de dados, os auxiliares da hipergeométrica para K = 31 (47) e a saída em
# saida/tempo_exposicao/, nunca em saida/grade/.
#
#   DESFECHO = casos | obitos | letalidade    VIZINHOS = queen | rook

source(if (file.exists("../caminhos.R")) "../caminhos.R" else "caminhos.R")

AQUI <- ESPACIAL
DEST <- file.path(AQUI, "saida", "tempo_exposicao")
dir.create(DEST, showWarnings = FALSE, recursive = TRUE)
for (f in c("hyperg_aux_k31.rds", "tabela_f21_k31.rds"))
  if (!file.exists(file.path(AQUI, "saida", f))) stop("rodar antes o 47_estender_hyperg_k31.R: falta ", f)

src <- paste(readLines(file.path(AQUI, "20_grade_empirica.R"), encoding = "UTF-8",
                       warn = FALSE), collapse = "\n")
trocar1 <- function(txt, a, b) {
  m <- gregexpr(a, txt, fixed = TRUE)[[1]]
  if (m[1] == -1L || length(m) != 1L) stop("esperava 1 ocorrência de: ", a)
  sub(a, b, txt, fixed = TRUE)
}
src <- trocar1(src, '"hyperg_aux.rds"', '"hyperg_aux_k31.rds"')
src <- trocar1(src, '"tabela_f21_estendida.rds"', '"tabela_f21_k31.rds"')
src <- trocar1(src, '"resultados/estimacoes/matrizes_estimacao/principal"',
               '"resultados/estimacoes/matrizes_estimacao/sensibilidade_tempo_exposicao"')
src <- trocar1(src, 'dir.create(file.path(AQUI, "saida", "grade"), showWarnings = FALSE, recursive = TRUE)',
               sprintf('dir.create("%s", showWarnings = FALSE, recursive = TRUE)', DEST))
src <- trocar1(src, 'file.path(AQUI, "saida", "grade", paste0(DESFECHO, "_", VIZINHOS, ".rds"))',
               sprintf('file.path("%s", paste0(DESFECHO, "_", VIZINHOS, ".rds"))', DEST))
src <- trocar1(src, "n <- length(y); K <- length(covariaveis)",
               paste("n <- length(y); K <- length(covariaveis)",
                     'stopifnot(K == 31L, "z_dias.primeiro.registro" %in% covariaveis)', sep = "\n"))
stopifnot(!grepl('"saida", "grade"', src, fixed = TRUE))

if (Sys.getenv("DESFECHO") == "") Sys.setenv(DESFECHO = "casos")
if (Sys.getenv("VIZINHOS") == "") Sys.setenv(VIZINHOS = "queen")
if (Sys.getenv("MC3_SORTEIOS") == "") Sys.setenv(MC3_SORTEIOS = "1000000")
if (Sys.getenv("MC3_BURN") == "") Sys.setenv(MC3_BURN = "200000")
cat("=== sensibilidade de tempo de exposição, configuração principal ===\n")
eval(parse(text = src), envir = globalenv())
