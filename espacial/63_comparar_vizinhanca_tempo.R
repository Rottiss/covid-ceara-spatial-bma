# Queen contra Rook na Seção 4.8 (sensibilidade ao tempo desde o primeiro
# registro). O 56 fez essa tabela só com a Queen. Aqui a mesma conta roda com
# a Rook, para conferir se a conclusão sobrevive à troca de matriz, do mesmo
# jeito que o 61 fez para a especificação principal e o 62 para o SEM.
#
# Reusa quatro peças da especificação "com dias": 48 (MQO e SAR), 53
# (especificação dentro do espaço amostral), 55 (SEM) e 54/64 (fator de Bayes
# agregado por Monte Carlo direto, a0 = 1,01, 20.000 modelos nas duas
# matrizes, para comparar com o mesmo desenho da Seção 4.4).
#
# Exige: saida/tempo_exposicao/{<d>,especificacao_<d>,sem_<d>}_<viz>.rds e
# saida/tempo_exposicao/sensibilidade_prior_rho_<d>[_rook].csv
# Saída: saida/tempo_exposicao/comparacao_vizinhanca_tempo.csv

source(if (file.exists("../caminhos.R")) "../caminhos.R" else "caminhos.R")

AQUI <- ESPACIAL
TE <- file.path(AQUI, "saida", "tempo_exposicao")
DIAS <- "z_dias.primeiro.registro"
D <- function(p) { m <- mean(p); K <- length(p)
  sd(p) * sqrt((K - 1) / K) / sqrt(m * (1 - m)) }

arq_pri <- function(d, viz) file.path(TE, sprintf("sensibilidade_prior_rho_%s%s.csv",
                                                    d, if (viz == "rook") "_rook" else ""))

linha <- function(d, viz) {
  g <- readRDS(file.path(TE, sprintf("%s_%s.rds", d, viz)))
  s <- readRDS(file.path(TE, sprintf("sem_%s_%s.rds", d, viz)))
  pri <- read.csv(arq_pri(d, viz))
  pip_mqo <- setNames(g$res$MQO$pip, NULL); pip_sar <- g$res$SAR$pip
  data.frame(
    desfecho = d, vizinhos = viz, K = g$K,
    D_mqo = D(pip_mqo), D_sar = D(pip_sar), D_sem = D(s$pip),
    razao_sar = D(pip_sar) / D(pip_mqo), razao_sem = D(s$pip) / D(pip_mqo),
    acima_mqo = sum(pip_mqo >= 0.5), acima_sar = sum(pip_sar >= 0.5), acima_sem = sum(s$pip >= 0.5),
    pip_dias_mqo = if (g$K == 31) g$tabela$mqo[g$tabela$covariavel == DIAS] else NA,
    pip_dias_sar = if (g$K == 31) g$tabela$sar[g$tabela$covariavel == DIAS] else NA,
    pip_dias_sem = if (g$K == 31) s$pip[s$covariaveis == DIAS] else NA,
    log_bf_agregado = pri$log_bf_agregado_mc[pri$a0 == 1.01],
    p_sar_agregada = pri$p_sar_agregada_mc[pri$a0 == 1.01],
    stringsAsFactors = FALSE)
}

L <- list()
for (d in c("casos", "obitos", "letalidade"))
  for (viz in c("queen", "rook"))
    L[[length(L) + 1]] <- linha(d, viz)
tab <- do.call(rbind, L)
write.csv(tab, file.path(TE, "comparacao_vizinhanca_tempo.csv"), row.names = FALSE)

op <- options(width = 200)
print(tab, digits = 3, row.names = FALSE)

cat("\n=== Queen contra Rook, especificacao com dias ===\n")
for (d in c("casos", "obitos", "letalidade")) {
  q <- tab[tab$desfecho == d & tab$vizinhos == "queen", ]
  r <- tab[tab$desfecho == d & tab$vizinhos == "rook", ]
  cat(sprintf("  %-11s D_sem %.4f -> %.4f | acima_sem %d -> %d | PIP dias SEM %.3f -> %.3f | log BF %.2f -> %.2f\n",
              d, q$D_sem, r$D_sem, q$acima_sem, r$acima_sem, q$pip_dias_sem, r$pip_dias_sem,
              q$log_bf_agregado, r$log_bf_agregado))
}
cat("\n  maior diferenca de log BF entre matrizes, com dias:",
    round(max(abs(tab$log_bf_agregado[tab$vizinhos == "queen"] -
                  tab$log_bf_agregado[tab$vizinhos == "rook"])), 3), "\n")
options(op)
cat("\nGravado em saida/tempo_exposicao/comparacao_vizinhanca_tempo.csv\n")
