# Auditoria: o 33 ainda reproduz os .rds gravados em 09/09?
#
# Os logs de saida/log_especificacao_*.txt sao de 09/09 21:05 e o
# 33_probabilidade_espacial.R foi alterado em 10/09 02:32. A mensagem impressa
# mudou (o script passou a dizer que zero visitas NAO autoriza o limite 1/N).
# Falta saber se a estimacao mudou junto, o que tornaria os .rds desatualizados.
#
# Roda o 33 atual com a saida desviada para saida/auditoria/, sem tocar nos
# arquivos originais, e compara.
#
#   DESFECHO = casos (padrao)    VIZINHOS = queen (padrao)
#
# Saida: saida/auditoria/especificacao_<desfecho>_<viz>.rds

source(if (file.exists("../caminhos.R")) "../caminhos.R" else "caminhos.R")

AQUI <- ESPACIAL
DEST <- file.path(AQUI, "saida", "auditoria")
dir.create(DEST, showWarnings = FALSE, recursive = TRUE)

if (Sys.getenv("DESFECHO") == "") Sys.setenv(DESFECHO = "casos")
if (Sys.getenv("VIZINHOS") == "") Sys.setenv(VIZINHOS = "queen")

trocar1 <- function(txt, a, b) {
  m <- gregexpr(a, txt, fixed = TRUE)[[1]]
  if (m[1] == -1L || length(m) != 1L) stop("esperava 1 ocorrencia de: ", a)
  sub(a, b, txt, fixed = TRUE)
}

src <- paste(readLines(file.path(AQUI, "33_probabilidade_espacial.R"),
                       encoding = "UTF-8", warn = FALSE), collapse = "\n")
src <- trocar1(src, 'file.path(AQUI, "saida", sprintf("especificacao_%s_%s.rds", DESFECHO, VIZINHOS))',
               sprintf('file.path("%s", sprintf("especificacao_%%s_%%s.rds", DESFECHO, VIZINHOS))', DEST))
stopifnot(!grepl('file.path(AQUI, "saida", sprintf("especificacao', src, fixed = TRUE))

eval(parse(text = src), envir = globalenv())

# --- comparacao com o arquivo original
d <- Sys.getenv("DESFECHO"); v <- Sys.getenv("VIZINHOS")
novo <- readRDS(file.path(DEST, sprintf("especificacao_%s_%s.rds", d, v)))
velho <- readRDS(file.path(AQUI, "saida", sprintf("especificacao_%s_%s.rds", d, v)))

rho_medio <- function(e) {
  w <- (e$a$w_rho + e$b$w_rho) / 2; w <- w / sum(w); sum(w * e$rho)
}

cat("\n=== AUDITORIA: novo contra gravado ===\n")
cat(sprintf("  p_sar cadeia a: %.6f -> %.6f\n", velho$a$p_sar, novo$a$p_sar))
cat(sprintf("  p_sar cadeia b: %.6f -> %.6f\n", velho$b$p_sar, novo$b$p_sar))
cat(sprintf("  visitas ao MQO: %d / %d -> %d / %d\n",
            velho$a$n_ols, velho$b$n_ols, novo$a$n_ols, novo$b$n_ols))
cat(sprintf("  rho medio: %.6f -> %.6f\n", rho_medio(velho), rho_medio(novo)))
cat(sprintf("  maior diferenca de PIP (SAR): %.3e\n",
            max(abs(velho$a$pip_sar - novo$a$pip_sar))))
identicos <- isTRUE(all.equal(velho, novo, tolerance = 1e-12))
cat("  objetos identicos:", identicos, "\n")
if (!identicos) cat("  ATENCAO: os .rds gravados nao correspondem ao script atual.\n")
