# Compara a grade de sensibilidade a priori entre as duas matrizes de
# vizinhanca. Queen e Rook, cinco combinacoes, tres variaveis dependentes:
# 30 estimacoes ao todo.
#
# GRAVA saida/priori/consolidado_vizinhanca.csv, regeneravel a partir dos .rds,
# e imprime um resumo no console. Nao altera nenhum .rds.
#
# A secao da covariavel que atravessa as cinco prioris usa as PIP sob MQO, que
# nao dependem de W: o resultado e identico nas duas matrizes por construcao e
# nao informa nada sobre a vizinhanca.

source(if (file.exists("../caminhos.R")) "../caminhos.R" else "caminhos.R")

S <- file.path(ESPACIAL, "saida")
Dst <- function(p) { m <- mean(p); K <- length(p)
  sd(p) * sqrt((K - 1) / K) / sqrt(m * (1 - m)) }
cam <- function(d, viz, cb) if (cb == "binbeta_hyperg")
  file.path(S, "grade", sprintf("%s_%s.rds", d, viz)) else
  file.path(S, "priori", sprintf("%s_%s_%s.rds", d, viz, cb))

COMBOS <- c("binbeta_hyperg", "binbeta_uip", "binbeta5_hyperg",
            "binomial5_hyperg", "uniforme_hyperg")
NOMES <- c(binbeta_hyperg = "Binomial-Beta K/2 + hyper-g (TCC)",
           binbeta_uip = "Binomial-Beta K/2 + g fixo",
           binbeta5_hyperg = "Binomial-Beta E[k]=5",
           binomial5_hyperg = "Binomial 5/K",
           uniforme_hyperg = "Uniforme sobre 2^K")
DESF <- c("casos", "obitos", "letalidade")

linhas <- list(); pip <- list()
for (d in DESF) for (viz in c("queen", "rook")) for (cb in COMBOS) {
  a <- cam(d, viz, cb)
  if (!file.exists(a)) { message("ausente: ", a); next }
  o <- readRDS(a); tb <- o$tabela
  pip[[paste(d, viz, cb)]] <- setNames(tb$mqo, tb$covariavel)
  linhas[[length(linhas) + 1]] <- data.frame(
    desfecho = d, vizinhos = viz, combinacao = cb,
    acima_mqo = sum(tb$mqo >= .5), acima_sar = sum(tb$sar >= .5),
    ref_sar = o$res$SAR$tam / o$K,
    D_mqo = Dst(tb$mqo), D_sar = Dst(tb$sar), razao = Dst(tb$sar) / Dst(tb$mqo),
    stringsAsFactors = FALSE)
}
tab <- do.call(rbind, linhas)
write.csv(tab, file.path(S, "priori", "consolidado_vizinhanca.csv"), row.names = FALSE)

cat("\n=== Acima de 0,50 sob SAR e razao de D, por matriz ===\n\n")
for (d in DESF) {
  cat(sprintf("--- %s ---\n", d))
  cat(sprintf("  %-34s %13s %15s\n", "combinacao", "SAR Q / R", "razao D Q / R"))
  for (cb in COMBOS) {
    q <- tab[tab$desfecho == d & tab$vizinhos == "queen" & tab$combinacao == cb, ]
    r <- tab[tab$desfecho == d & tab$vizinhos == "rook" & tab$combinacao == cb, ]
    if (!nrow(q) || !nrow(r)) next
    marca <- if (q$acima_sar != r$acima_sar) "  <-- muda" else ""
    cat(sprintf("  %-34s %6d %6d %7.3f %7.3f%s\n", NOMES[cb],
                q$acima_sar, r$acima_sar, q$razao, r$razao, marca))
  }
  s <- tab[tab$desfecho == d, ]
  cat(sprintf("  razao de D: Queen [%.3f; %.3f] | Rook [%.3f; %.3f]\n\n",
              min(s$razao[s$vizinhos == "queen"]), max(s$razao[s$vizinhos == "queen"]),
              min(s$razao[s$vizinhos == "rook"]), max(s$razao[s$vizinhos == "rook"])))
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

cat("=== A covariavel que atravessa as cinco prioris, por matriz ===\n")
for (d in DESF) for (viz in c("queen", "rook")) {
  chaves <- paste(d, viz, COMBOS)
  chaves <- chaves[chaves %in% names(pip)]
  M <- alinhar(lapply(chaves, function(k) pip[[k]]))
  n50 <- rowSums(M >= 0.5)
  todas <- names(n50)[n50 == length(chaves)]
  amp <- apply(M, 1, max) - apply(M, 1, min)
  cat(sprintf("  %-11s %-6s %s (amplitude %.4f)\n", d, viz,
              paste(todas, collapse = ", "),
              if (length(todas)) max(amp[todas]) else NA))
}
cat("\nGravado: saida/priori/consolidado_vizinhanca.csv\n")
