# Consolida o refino do fator de Bayes agregado de óbitos (57): duas sementes
# com 200.000 modelos cada, por especificação (principal e com os dias).
#
# Por valor de a0: combinação das duas estimativas do log do fator de Bayes por
# log-sum-exp, média das duas para a probabilidade do SAR e para a mediana
# condicional, mínimo e máximo condicionais sobre as duas, tamanho efetivo mínimo
# e máximo entre sementes e diferença entre as sementes como erro de Monte Carlo.
#
# Duas observações sobre a combinação, levantadas na auditoria de 11/09/2026:
#
#   1. O log-sum-exp e a média dos logs diferem em 0,0001 nesta grade, contra
#      0,016 a 0,030 de diferença entre as sementes. O estimador exato juntaria
#      os 400.000 modelos no numerador e no denominador, o que exigiria os
#      valores por modelo; o 43 grava só o agregado de cada rodada.
#   2. O tamanho efetivo NÃO é somado. A soma valeria para uma estimativa única
#      sobre as duas amostras; aqui há duas estimativas combinadas, e o que
#      descreve a precisão de cada uma é o tamanho efetivo dela. O artigo cita o
#      menor dos dois.
#
# As colunas do 43 são mantidas com o mesmo nome, para que 50 e 56 leiam o
# consolidado no lugar de sensibilidade_prior_rho_obitos.csv.
#
# Saída: saida/refino_bf/consolidado_<espec>_<desfecho>.csv
#
# DESFECHO = obitos (padrão) | casos | letalidade. O padrão reproduz exatamente
# o arquivo de óbitos lido por 38, 50 e 56.

source(if (file.exists("../caminhos.R")) "../caminhos.R" else "caminhos.R")

AQUI <- ESPACIAL
RB <- file.path(AQUI, "saida", "refino_bf")
SEMENTES <- c(20260910, 20260911)
DESFECHO <- Sys.getenv("DESFECHO", "obitos")
log_media <- function(x, y) { m <- pmax(x, y); m + log((exp(x - m) + exp(y - m)) / 2) }

for (espec in c("principal", "tempo")) {
  s <- lapply(SEMENTES, function(x)
    read.csv(file.path(RB, sprintf("%s_%s_n200000_s%d.csv", espec, DESFECHO, x))))
  stopifnot(identical(s[[1]]$a0, s[[2]]$a0))
  a <- s[[1]]; b <- s[[2]]
  out <- data.frame(
    desfecho = DESFECHO, a0 = a$a0, n_modelos = a$n_modelos + b$n_modelos,
    log_bf_cond_mediano = (a$log_bf_cond_mediano + b$log_bf_cond_mediano) / 2,
    log_bf_cond_min = pmin(a$log_bf_cond_min, b$log_bf_cond_min),
    log_bf_cond_max = pmax(a$log_bf_cond_max, b$log_bf_cond_max),
    log_bf_agregado_mc = log_media(a$log_bf_agregado_mc, b$log_bf_agregado_mc),
    p_sar_agregada_mc = (a$p_sar_agregada_mc + b$p_sar_agregada_mc) / 2,
    dp_log_bf_blocos = pmax(a$dp_log_bf_blocos, b$dp_log_bf_blocos),
    ess_sar = pmin(a$ess_sar, b$ess_sar), ess_mqo = pmin(a$ess_mqo, b$ess_mqo),
    ess_sar_max = pmax(a$ess_sar, b$ess_sar), ess_mqo_max = pmax(a$ess_mqo, b$ess_mqo),
    dif_log_bf_sementes = abs(a$log_bf_agregado_mc - b$log_bf_agregado_mc))
  arq <- file.path(RB, sprintf("consolidado_%s_%s.csv", espec, DESFECHO))
  write.csv(out, arq, row.names = FALSE)
  cat("==", espec, "\n"); print(out, digits = 4, row.names = FALSE)
  cat("Gravado em", arq, "\n\n")
}
