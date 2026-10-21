# SAR da grade sob outra priori Beta(a0, a0) para rho.
#
# A auditoria de 10/09/2026 perguntou quanto a seleção sob SAR em letalidade
# depende da priori de rho. A grade usa a0 = 1,01, quase uniforme em (-1, 1).
# Aqui a mesma cadeia do 20_grade_empirica.R roda com outro a0.
#
# Executa o 20 com três trocas de texto, sem alterar o arquivo: o valor de a0;
# só a cadeia do SAR, porque o MQO não depende de rho e é lido da grade gravada;
# e a saída em saida/sensibilidade_a0/, nunca em saida/grade/.
#
#   DESFECHO = casos | obitos | letalidade    VIZINHOS = queen | rook
#   A0_RHO   = parâmetro da Beta(a0, a0), padrão 10

source(if (file.exists("../caminhos.R")) "../caminhos.R" else "caminhos.R")

AQUI <- ESPACIAL
A0 <- as.numeric(Sys.getenv("A0_RHO", "10"))
DEST <- file.path(AQUI, "saida", "sensibilidade_a0")
dir.create(DEST, showWarnings = FALSE, recursive = TRUE)

src <- paste(readLines(file.path(AQUI, "20_grade_empirica.R"), encoding = "UTF-8",
                       warn = FALSE), collapse = "\n")
trocar1 <- function(txt, a, b) {
  m <- gregexpr(a, txt, fixed = TRUE)[[1]]
  if (m[1] == -1L || length(m) != 1L) stop("esperava 1 ocorrência de: ", a)
  sub(a, b, txt, fixed = TRUE)
}
alvo_grade <- 'file.path(AQUI, "saida", "grade", paste0(DESFECHO, "_", VIZINHOS, ".rds"))'
src <- trocar1(src, paste0(alvo_grade, ")"),
               sprintf('file.path("%s", paste0(DESFECHO, "_", VIZINHOS, "_a0_%s.rds")))',
                       DEST, format(A0)))
src <- trocar1(src, "a0 <- 1.01", sprintf("a0 <- %s", format(A0)))
src <- trocar1(src, "res <- list()", paste0("res <- list(MQO = readRDS(", alvo_grade, ")$res$MQO)"))
src <- trocar1(src, 'for (rot in c("MQO", "SAR")) {', 'for (rot in c("SAR")) {')
ocorr <- gregexpr(alvo_grade, src, fixed = TRUE)[[1]]
stopifnot(length(ocorr) == 1L, ocorr[1] != -1L,
          grepl(paste0("readRDS(", alvo_grade), src, fixed = TRUE))

if (Sys.getenv("DESFECHO") == "") Sys.setenv(DESFECHO = "letalidade")
if (Sys.getenv("VIZINHOS") == "") Sys.setenv(VIZINHOS = "queen")
if (Sys.getenv("MC3_SORTEIOS") == "") Sys.setenv(MC3_SORTEIOS = "1000000")
if (Sys.getenv("MC3_BURN") == "") Sys.setenv(MC3_BURN = "200000")
cat(sprintf("=== SAR com a0 = %s ===\n", format(A0)))
eval(parse(text = src), envir = globalenv())
