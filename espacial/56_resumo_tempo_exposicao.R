# Resumo da sensibilidade de tempo de exposição (dias desde o primeiro registro)
# contra as estimações principais, na mesma configuração de cadeia.
#
# Lê as saídas do 48 (MQO e SAR), 53 (especificação dentro do espaço amostral),
# 54 (fator de Bayes agregado) e 55 (SEM), e os arquivos equivalentes das
# estimações principais (20, 33, 43, 31). Vizinhança Queen.
#
# Saída: saida/tempo_exposicao/resumo_sensibilidade.csv

source(if (file.exists("../caminhos.R")) "../caminhos.R" else "caminhos.R")

AQUI <- ESPACIAL
TE <- file.path(AQUI, "saida", "tempo_exposicao")
DIAS <- "z_dias.primeiro.registro"
D <- function(p) { m <- mean(p); K <- length(p)
  sd(p) * sqrt((K - 1) / K) / sqrt(m * (1 - m)) }

rho_resumo <- function(e) {
  w <- (e$a$w_rho + e$b$w_rho) / 2; w <- w / sum(w); ac <- cumsum(w)
  c(media = sum(w * e$rho), inf = e$rho[which.max(ac >= 0.025)], sup = e$rho[which.max(ac >= 0.975)],
    p_sar = (e$a$p_sar + e$b$p_sar) / 2)
}
linha <- function(versao, d, g, e, s, pri) {
  pip_mqo <- setNames(g$res$MQO$pip, NULL); pip_sar <- g$res$SAR$pip
  cov_g <- g$tabela$covariavel[order(match(g$tabela$covariavel, g$tabela$covariavel))]
  r <- rho_resumo(e)
  dias_idx <- function(covs) match(DIAS, covs)
  data.frame(
    versao = versao, desfecho = d, K = g$K,
    rho_media = r[["media"]], rho_inf = r[["inf"]], rho_sup = r[["sup"]], freq_sar = r[["p_sar"]],
    log_bf_agregado = pri$log_bf_agregado_mc[pri$a0 == 1.01],
    p_sar_a0_1.01 = pri$p_sar_agregada_mc[pri$a0 == 1.01], p_sar_a0_10 = pri$p_sar_agregada_mc[pri$a0 == 10],
    k_mqo = g$res$MQO$tam, k_sar = g$res$SAR$tam, k_sem = s$tam,
    D_mqo = D(pip_mqo), D_sar = D(pip_sar), D_sem = D(s$pip),
    razao_sar = D(pip_sar) / D(pip_mqo), razao_sem = D(s$pip) / D(pip_mqo),
    acima_mqo = sum(pip_mqo >= 0.5), acima_sar = sum(pip_sar >= 0.5), acima_sem = sum(s$pip >= 0.5),
    ref_sar = g$res$SAR$tam / g$K, ref_sem = s$tam / g$K,
    conv_sem = max(abs(s$a$pip - s$b$pip)),
    pip_dias_mqo = if (g$K == 31) g$tabela$mqo[g$tabela$covariavel == DIAS] else NA,
    pip_dias_sar = if (g$K == 31) g$tabela$sar[g$tabela$covariavel == DIAS] else NA,
    pip_dias_sem = if (g$K == 31) s$pip[s$covariaveis == DIAS] else NA,
    stringsAsFactors = FALSE)
}

# Fator de Bayes agregado: casos e óbitos usam o refino consolidado (57 e
# 58), com tamanho efetivo maior; letalidade fica nos 20.000 modelos padrao
# do 43/54, pois seu tamanho efetivo ja e alto (Seção 4.4).
RB <- file.path(AQUI, "saida", "refino_bf")
REFINADOS <- c("casos", "obitos")
pri_principal <- function(d) read.csv(if (d %in% REFINADOS) file.path(RB, sprintf("consolidado_principal_%s.csv", d))
                                      else file.path(AQUI, "saida", sprintf("sensibilidade_prior_rho_%s.csv", d)))
pri_tempo <- function(d) read.csv(if (d %in% REFINADOS) file.path(RB, sprintf("consolidado_tempo_%s.csv", d))
                                  else file.path(TE, sprintf("sensibilidade_prior_rho_%s.csv", d)))

L <- list()
for (d in c("casos", "obitos", "letalidade")) {
  L[[length(L) + 1]] <- linha("principal", d,
    readRDS(file.path(AQUI, "saida", "grade", sprintf("%s_queen.rds", d))),
    readRDS(file.path(AQUI, "saida", sprintf("especificacao_%s_queen.rds", d))),
    readRDS(file.path(AQUI, "saida", "grade", sprintf("sem_%s_queen.rds", d))),
    pri_principal(d))
  L[[length(L) + 1]] <- linha("com dias", d,
    readRDS(file.path(TE, sprintf("%s_queen.rds", d))),
    readRDS(file.path(TE, sprintf("especificacao_%s_queen.rds", d))),
    readRDS(file.path(TE, sprintf("sem_%s_queen.rds", d))),
    pri_tempo(d))
}
tab <- do.call(rbind, L)
write.csv(tab, file.path(TE, "resumo_sensibilidade.csv"), row.names = FALSE)
op <- options(width = 250)
print(tab[, c("versao", "desfecho", "rho_media", "rho_inf", "rho_sup", "freq_sar", "log_bf_agregado",
              "p_sar_a0_1.01", "p_sar_a0_10")], digits = 3, row.names = FALSE)
cat("\n")
print(tab[, c("versao", "desfecho", "k_mqo", "k_sar", "k_sem", "D_mqo", "D_sar", "D_sem", "razao_sar",
              "razao_sem", "acima_mqo", "acima_sar", "acima_sem", "ref_sar", "ref_sem", "conv_sem",
              "pip_dias_mqo", "pip_dias_sar", "pip_dias_sem")], digits = 3, row.names = FALSE)
cat("\n-- SEM com dias, covariaveis acima de 0,50:\n")
for (d in c("casos", "obitos", "letalidade")) {
  s <- readRDS(file.path(TE, sprintf("sem_%s_queen.rds", d)))
  cat(" ", d, ":", paste(s$covariaveis[s$pip >= 0.5], collapse = ", "), "\n")
}
