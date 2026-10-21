# Consolida a grade de sensibilidade a priori sob SEM (75) ao lado de MQO e SAR
# (32, e a linha de base em saida/grade).
#
# Tudo e alinhado pelo nome da covariavel, nunca por posicao.
#
# GRAVA saida/priori_sem/consolidado_sem.csv, regeneravel a partir dos .rds.
# Nenhum .rds e alterado.

source(if (file.exists("../caminhos.R")) "../caminhos.R" else "caminhos.R")

S <- file.path(ESPACIAL, "saida")
COMB <- c("binbeta_hyperg", "binbeta_uip", "binbeta5_hyperg", "binomial5_hyperg", "uniforme_hyperg")

Dst <- function(p) { m <- mean(p); Kk <- length(p)
  sd(p) * sqrt((Kk - 1) / Kk) / sqrt(m * (1 - m)) }

linhas <- list()
for (d in c("casos", "obitos", "letalidade")) for (v in c("queen", "rook")) for (cb in COMB) {
  sem <- readRDS(file.path(S, "priori_sem", sprintf("%s_%s_%s.rds", d, v, cb)))
  stopifnot(sem$n_mc == 1e6, sem$burn == 2e5, !anyDuplicated(names(sem$pip)))
  tab <- if (cb == "binbeta_hyperg") readRDS(file.path(S, "grade", sprintf("%s_%s.rds", d, v)))$tabela else
    readRDS(file.path(S, "priori", sprintf("%s_%s_%s.rds", d, v, cb)))$tabela
  stopifnot(!anyDuplicated(tab$covariavel), setequal(tab$covariavel, names(sem$pip)))
  covs <- sort(tab$covariavel)
  mqo <- setNames(tab$mqo, tab$covariavel)[covs]
  sar <- setNames(tab$sar, tab$covariavel)[covs]
  pse <- sem$pip[covs]
  topo <- function(p) names(p)[which.max(p)]
  linhas[[length(linhas) + 1]] <- data.frame(
    desfecho = d, vizinhos = v, combinacao = cb,
    tam_sem = sem$tam, conv_sem = sem$conv, fora_tab = sem$fora_tab,
    reproducao_31 = sem$reproducao_31,
    acima_mqo = sum(mqo >= 0.5), acima_sar = sum(sar >= 0.5), acima_sem = sum(pse >= 0.5),
    D_mqo = Dst(mqo), D_sar = Dst(sar), D_sem = Dst(pse),
    razao_sar = Dst(sar) / Dst(mqo), razao_sem = Dst(pse) / Dst(mqo),
    media_sem = mean(pse),
    topo_mqo = topo(mqo), topo_sar = topo(sar), topo_sem = topo(pse),
    pip_topo_sem = max(pse),
    spearman_sar_sem = cor(sar, pse, method = "spearman"),
    acima_sem_lista = paste(covs[pse >= 0.5], collapse = " "),
    stringsAsFactors = FALSE)
}
res <- do.call(rbind, linhas)

# a linha de base do 75 deve coincidir com a do D do SEM ja calculado no 31
base <- res[res$combinacao == "binbeta_hyperg", ]
stopifnot(all(base$reproducao_31 < 1e-12))

write.csv(res, file.path(S, "priori_sem", "consolidado_sem.csv"), row.names = FALSE)
op <- options(width = 200)
print(transform(res[, c("desfecho", "vizinhos", "combinacao", "tam_sem", "acima_mqo", "acima_sar",
                        "acima_sem", "D_mqo", "D_sem", "razao_sar", "razao_sem", "media_sem",
                        "topo_sar", "topo_sem", "pip_topo_sem", "conv_sem", "fora_tab")],
                tam_sem = round(tam_sem, 2), D_mqo = round(D_mqo, 3), D_sem = round(D_sem, 3),
                razao_sar = round(razao_sar, 3), razao_sem = round(razao_sem, 3),
                media_sem = round(media_sem, 3), pip_topo_sem = round(pip_topo_sem, 3),
                conv_sem = round(conv_sem, 4)), row.names = FALSE)
options(op)
cat("\nrazao de D do SEM < 1 em", sum(res$razao_sem < 1), "de", nrow(res), "\n")
cat("faixa da razao SEM:", round(range(res$razao_sem), 3), "\n")
cat("maior discordancia entre cadeias SEM:", round(max(res$conv_sem), 4), "\n")
cat("avaliacoes fora da tabela:", sum(res$fora_tab), "\n")
for (d in unique(res$desfecho)) {
  r <- res[res$desfecho == d, ]
  cat(sprintf("%s: covariavel de maior PIP sob SEM nas 10 combinacoes: %s\n", d,
              paste(names(table(r$topo_sem)), table(r$topo_sem), sep = " x", collapse = "; ")))
}
cat("\nGravado: saida/priori_sem/consolidado_sem.csv\n")
