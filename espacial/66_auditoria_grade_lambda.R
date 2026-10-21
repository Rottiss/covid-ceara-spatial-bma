# Auditoria: o erro de quadratura da grade de lambda do SEM.
#
# O 31_sem_hyperg.R integra lambda por Simpson em 101 pontos e o comentario do
# script afirma erro maximo de 1,1e-10 contra 801 pontos. Esta e a peca menos
# coberta por teste do projeto: a verificacao interna do 31 confere o RSS
# pre-computado contra ajuste explicito, nao a integracao.
#
# Aqui a afirmacao e medida de forma independente: carrega as definicoes do 31
# ate antes da cadeia, avalia lm_sem em modelos sorteados com ng = 51, 101, 201
# e 801, e compara. Nenhuma cadeia e executada e nada e gravado.
#
#   DESFECHO = casos (padrao)    VIZINHOS = queen (padrao)

source(if (file.exists("../caminhos.R")) "../caminhos.R" else "caminhos.R")

AQUI <- ESPACIAL
if (Sys.getenv("DESFECHO") == "") Sys.setenv(DESFECHO = "casos")
if (Sys.getenv("VIZINHOS") == "") Sys.setenv(VIZINHOS = "queen")

src <- readLines(file.path(AQUI, "31_sem_hyperg.R"), encoding = "UTF-8", warn = FALSE)
corte <- grep('cat("=== Verificacao: RSS pre-computado contra ajuste explicito ===', src, fixed = TRUE)
stopifnot(length(corte) == 1L)

avaliar <- function(ng, modelos) {
  Sys.setenv(N_GRADE_LAMBDA = as.character(ng))
  e <- new.env(parent = globalenv())
  eval(parse(text = paste(src[seq_len(corte - 1L)], collapse = "\n")), envir = e)
  stopifnot(length(e$lambda) == ng)
  vapply(modelos, e$lm_sem, numeric(1))
}

set.seed(20260914)
K <- 30
modelos <- c(list(integer(0)), replicate(60, sort(sample(K, sample(1:K, 1))), simplify = FALSE))
cat(sprintf("=== %s, %s: %d modelos, grade de lambda ===\n\n",
            Sys.getenv("DESFECHO"), Sys.getenv("VIZINHOS"), length(modelos)))

ref <- avaliar(801, modelos)
for (ng in c(51, 101, 201, 401)) {
  v <- avaliar(ng, modelos)
  d <- abs(v - ref)
  cat(sprintf("  ng = %3d | erro maximo %.2e | erro mediano %.2e\n", ng, max(d), median(d)))
}

# Simpson exige numero impar de pontos. O script nao valida isso.
cat("\n=== Pesos de Simpson para ng par (armadilha, nao afeta os resultados atuais) ===\n")
for (ng in c(100, 101)) {
  pesos <- c(1, rep(c(4, 2), length.out = ng - 2), 1)
  cat(sprintf("  ng = %d | primeiros %s | ultimos %s | soma %d | esperado %d\n",
              ng, paste(head(pesos, 4), collapse = ","),
              paste(tail(pesos, 4), collapse = ","), sum(pesos), 3 * (ng - 1)))
}
