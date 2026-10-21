# I de Moran do y simulado nos tres cenarios da primeira rodada.
#
# O 18_experimento_simulacao.R nao guardou esse numero, nem no .rds nem no log.
# O texto do artigo citava valores cuja origem nao era rastreavel a nenhum
# arquivo de saida. Esta rotina reproduz o DGP com as mesmas sementes e os
# mesmos betas calibrados, e recalcula.

source(if (file.exists("../caminhos.R")) "../caminhos.R" else "caminhos.R")

suppressPackageStartupMessages(library(spdep))
AQUI <- ESPACIAL
base <- readRDS(file.path(AQUI, "saida", "base_sar.rds"))
des  <- readRDS(file.path(AQUI, "saida", "desenho_simulacao.rds"))
cal  <- readRDS(file.path(AQUI, "saida", "calibracao_sinal.rds"))

X <- base$X; W <- base$W; n <- base$n; listw <- base$listw
Xc <- sweep(X, 2, colMeans(X), "-")
RHO_CTRL <- 0.6
S_inv <- solve(diag(n) - RHO_CTRL * W)

cenarios <- list(
  suaves   = list(idx = des$config$suaves$idx,
                  beta = rep(cal$suaves$beta_escolhido, 4),  rho = 0),
  rugosas  = list(idx = des$config$rugosas$idx,
                  beta = rep(cal$rugosas$beta_escolhido, 4), rho = 0),
  controle = list(idx = des$config$suaves$idx,
                  beta = rep(cal$suaves$beta_escolhido, 4),  rho = RHO_CTRL))
stopifnot(all(vapply(cenarios, function(z) length(z$idx) == 4L, logical(1))))

cat("=== I de Moran do y simulado, media de 10 replicacoes ===\n")
saida <- list()
for (nome in names(cenarios)) {
  cn <- cenarios[[nome]]
  mis <- numeric(10)
  for (r in seq_len(10)) {
    set.seed(1000L * which(names(cenarios) == nome) + r)
    sinal <- as.numeric(Xc[, cn$idx, drop = FALSE] %*% cn$beta)
    y <- if (cn$rho == 0) sinal + rnorm(n) else as.numeric(S_inv %*% (sinal + rnorm(n)))
    mis[r] <- spdep::moran.test(y, listw)$estimate[["Moran I statistic"]]
  }
  cat(sprintf("  %-9s rho %.1f | I de Moran medio %.4f | faixa [%.4f, %.4f]\n",
              nome, cn$rho, mean(mis), min(mis), max(mis)))
  saida[[length(saida) + 1]] <- data.frame(cenario = nome, rho = cn$rho,
                                           mi = mean(mis), mi_min = min(mis),
                                           mi_max = max(mis))
}
write.csv(do.call(rbind, saida), file.path(AQUI, "saida", "moran_cenarios.csv"),
          row.names = FALSE)
cat("\nGravado em saida/moran_cenarios.csv\n")
