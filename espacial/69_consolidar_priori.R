# Consolida a grade de sensibilidade a priori sob Queen, nas tres variaveis
# dependentes, em duas tabelas: uma por eixo do desenho em estrela.
#
#   Eixo 1, priori de modelo, com hyper-g mantida:
#     binbeta (linha de base) | binbeta5 | binomial5 | uniforme
#   Eixo 2, priori de g, com Binomial-Beta mantida:
#     hyperg (linha de base) | uip
#
# A linha de base vem de saida/grade/<desfecho>_queen.rds, gravada pelo 20; as
# variacoes vem de saida/priori/, gravadas pelo 32. Todas com a mesma
# configuracao: 2 cadeias, 1.000.000 de sorteios, 200.000 de descarte.
#
# GRAVA saida/priori/consolidado.csv, regeneravel a partir dos .rds, e imprime
# um resumo no console. Nao altera nenhum .rds.

source(if (file.exists("../caminhos.R")) "../caminhos.R" else "caminhos.R")

AQUI <- ESPACIAL
SAIDA <- file.path(AQUI, "saida")

Dst <- function(p) { m <- mean(p); K <- length(p)
  sd(p) * sqrt((K - 1) / K) / sqrt(m * (1 - m)) }

# rotulo -> caminho. A linha de base nao esta em saida/priori.
caminho <- function(desfecho, combo) {
  if (combo == "binbeta_hyperg")
    file.path(SAIDA, "grade", sprintf("%s_queen.rds", desfecho))
  else
    file.path(SAIDA, "priori", sprintf("%s_queen_%s.rds", desfecho, combo))
}

COMBOS <- c("binbeta_hyperg", "binbeta_uip", "binbeta5_hyperg",
            "binomial5_hyperg", "uniforme_hyperg")
NOMES <- c(binbeta_hyperg = "Binomial-Beta K/2 + hyper-g (TCC)",
           binbeta_uip = "Binomial-Beta K/2 + g fixo UIP",
           binbeta5_hyperg = "Binomial-Beta E[k]=5 + hyper-g",
           binomial5_hyperg = "Binomial 5/K + hyper-g",
           uniforme_hyperg = "Uniforme sobre 2^K + hyper-g")
DESFECHOS <- c("casos", "obitos", "letalidade")

linhas <- list(); pips <- list()
for (d in DESFECHOS) {
  for (cb in COMBOS) {
    arq <- caminho(d, cb)
    if (!file.exists(arq)) { message("ausente: ", arq); next }
    o <- readRDS(arq); tb <- o$tabela
    p <- setNames(tb$mqo, tb$covariavel)
    pips[[paste(d, cb)]] <- p
    linhas[[length(linhas) + 1]] <- data.frame(
      desfecho = d, combinacao = cb, nome = unname(NOMES[cb]),
      tam_mqo = o$res$MQO$tam, tam_sar = o$res$SAR$tam,
      acima_mqo = sum(tb$mqo >= 0.5), acima_sar = sum(tb$sar >= 0.5),
      D_mqo = Dst(tb$mqo), D_sar = Dst(tb$sar),
      razao = Dst(tb$sar) / Dst(tb$mqo),
      conv = max(o$res$MQO$conv, o$res$SAR$conv),
      stringsAsFactors = FALSE)
  }
}
tab <- do.call(rbind, linhas)
write.csv(tab, file.path(SAIDA, "priori", "consolidado.csv"), row.names = FALSE)

cat("\n=== Contagem acima de 0,50 e dispersao, por combinacao ===\n\n")
for (d in DESFECHOS) {
  s <- tab[tab$desfecho == d, ]
  if (!nrow(s)) next
  cat(sprintf("--- %s ---\n", d))
  cat(sprintf("  %-34s %5s %5s %7s %7s %7s %8s\n",
              "combinacao", ">=.50", "SAR", "D MQO", "D SAR", "razao", "tam MQO"))
  for (i in seq_len(nrow(s)))
    cat(sprintf("  %-34s %5d %5d %7.4f %7.4f %7.3f %8.2f\n",
                s$nome[i], s$acima_mqo[i], s$acima_sar[i],
                s$D_mqo[i], s$D_sar[i], s$razao[i], s$tam_mqo[i]))
  cat(sprintf("  razao de D entre %.3f e %.3f | discordancia maxima entre cadeias: %.4f\n\n",
              min(s$razao), max(s$razao), max(s$conv)))
}

# As tabelas de cada .rds vem ordenadas por PIP, e a ordem muda entre
# combinacoes. Empilhar os vetores com sapply alinha por POSICAO e troca as
# covariaveis de linha sem aviso. A juncao abaixo e por NOME, e o teste de
# permutacao confirma que o resultado
# nao depende da ordem das linhas.
alinhar <- function(lista) {
  covs <- sort(names(lista[[1]]))
  stopifnot(!anyDuplicated(covs))
  M <- vapply(lista, function(p) {
    stopifnot(!anyDuplicated(names(p)), setequal(names(p), covs))
    unname(p[covs])
  }, numeric(length(covs)))
  rownames(M) <- covs
  M
}
set.seed(16092026L)
.teste <- list(c(a = 1, b = 2, c = 3), c(c = 30, a = 10, b = 20))
stopifnot(identical(alinhar(.teste)["a", ], c(1, 10)),
          identical(alinhar(lapply(.teste, function(p) p[sample(length(p))])),
                    alinhar(.teste)))

cat("=== Quantas combinacoes cada covariavel supera 0,50, sob MQO ===\n")
for (d in DESFECHOS) {
  cbs <- COMBOS[paste(d, COMBOS) %in% names(pips)]
  if (!length(cbs)) next
  M <- alinhar(lapply(cbs, function(cb) pips[[paste(d, cb)]]))
  n50 <- rowSums(M >= 0.5)
  amp <- apply(M, 1, max) - apply(M, 1, min)
  r <- data.frame(covariavel = rownames(M), n = n50, min = apply(M, 1, min),
                  max = apply(M, 1, max), amplitude = amp)
  r <- r[r$n > 0, ]; r <- r[order(-r$n, -r$min), ]
  cat(sprintf("\n--- %s (%d combinacoes disponiveis) ---\n", d, length(cbs)))
  for (i in seq_len(nrow(r)))
    cat(sprintf("  %-16s %d de %d   min %.4f  max %.4f  amplitude %.4f\n",
                r$covariavel[i], r$n[i], length(cbs), r$min[i], r$max[i], r$amplitude[i]))
}
cat("\nGravado: saida/priori/consolidado.csv\n")
