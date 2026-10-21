# Erro das duas tabelas de interpolação de log 2F1, medido como no 15.
#
# O 15_otimizar_hyperg.R constrói a tabela e mede o erro na mesma execução, mas
# não grava a medição. O artigo citava o erro da tabela preliminar (80.001
# pontos) como se fosse o da tabela usada na grade (300.001 pontos). A auditoria
# de 10/09/2026 apontou a troca.
#
# Este script reaproveita as seções 1 a 4 do 15, com as mesmas sementes e os
# mesmos modelos, carrega as tabelas já gravadas em vez de reconstruí-las e grava
# o resultado em saida/erro_interpolacao.csv. Não grava mais nada.

source(if (file.exists("../caminhos.R")) "../caminhos.R" else "caminhos.R")

AQUI <- ESPACIAL
src <- readLines(file.path(AQUI, "15_otimizar_hyperg.R"), encoding = "UTF-8", warn = FALSE)
corte <- grep("^# -+ 5\\. ganho", src)
ini <- grep("^N_TAB <- as.integer", src)
fim <- grep("^cat\\(\"  construida em", src)
stopifnot(length(corte) == 1L, length(ini) == 1L, length(fim) == 1L, ini < fim, fim < corte)

# O 15 usa "tb" para o escalar sum(yc * Wyc); a tabela carregada leva outro nome.
carga <- c(
  'tab_carregada <- readRDS(file.path(AQUI, "saida", ARQ_TABELA))',
  'TAB <- tab_carregada$TAB',
  'N_TAB <- tab_carregada$N_TAB',
  'Z_MAX <- tab_carregada$Z_MAX',
  'z_tab <- seq(0, Z_MAX, length.out = N_TAB)')
codigo <- c(src[1:(ini - 1L)], carga, src[(fim + 1L):(corte - 1L)])

linhas <- list()
for (arq in c("tabela_f21.rds", "tabela_f21_estendida.rds")) {
  env <- new.env()
  assign("ARQ_TABELA", arq, envir = env)
  assign("saveRDS", function(...) invisible(NULL), envir = env)
  assign("write.csv", function(...) invisible(NULL), envir = env)
  cat("\n########", arq, "########\n")
  eval(parse(text = codigo), envir = env)
  linhas[[length(linhas) + 1L]] <- data.frame(
    tabela = arq, pontos = env$N_TAB, teto_r2 = env$Z_MAX,
    erro_log_max = env$pior, erro_logmarginal_max = max(abs(env$d)),
    modelos = length(env$d), stringsAsFactors = FALSE)
}
tab <- do.call(rbind, linhas)
print(tab, row.names = FALSE)
write.csv(tab, file.path(AQUI, "saida", "erro_interpolacao.csv"), row.names = FALSE)
cat("\nGravado em saida/erro_interpolacao.csv\n")
