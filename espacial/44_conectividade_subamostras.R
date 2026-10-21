# Conectividade das subamostras usadas no teste de escala.
#
# O 30_escala_n.R descarta subamostras com município sem vizinho (ilha). Isso
# não garante que a malha restante seja conexa: ela pode se dividir em blocos
# que não se tocam. A auditoria de 10/09/2026 apontou a diferença.
#
# O 30 não permite reconstruir as subamostras exatas que aceitou, porque cada
# execução do MC3 chama set.seed e reinicia o gerador. Por isso a propriedade
# do critério é medida aqui por sorteio independente: para cada tamanho, quantas
# subamostras aleatórias não têm ilha, e quantas dessas são de fato conexas.

source(if (file.exists("../caminhos.R")) "../caminhos.R" else "caminhos.R")

suppressPackageStartupMessages(library(spdep))

AQUI <- ESPACIAL
base <- readRDS(file.path(AQUI, "saida", "base_sar.rds"))
B <- (base$W > 0) * 1

componentes <- function(Bs) {
  nb <- lapply(seq_len(nrow(Bs)), function(i) {
    v <- which(Bs[i, ] > 0)
    if (length(v)) as.integer(v) else 0L
  })
  class(nb) <- "nb"
  spdep::n.comp.nb(nb)$nc
}

R <- as.integer(Sys.getenv("SORTEIOS_SUBAMOSTRA", "2000"))
tamanhos <- c(80L, 120L, 160L)

cat("=== Subamostras sem ilhas contra subamostras conexas ===\n")
cat(sprintf("  malha completa, n = 184: %d componente(s)\n", componentes(B)))
cat(sprintf("  sorteios por tamanho: %d\n\n", R))

set.seed(20260910L)
linhas <- list()
for (nn in tamanhos) {
  sem_ilha <- 0L; nc <- integer(0)
  for (r in seq_len(R)) {
    idx <- sort(sample(184L, nn)); Bs <- B[idx, idx]
    if (any(rowSums(Bs) == 0)) next
    sem_ilha <- sem_ilha + 1L
    nc <- c(nc, componentes(Bs))
  }
  conexas <- sum(nc == 1L)
  cat(sprintf("  n = %3d | sem ilhas: %4d de %d | conexas entre essas: %4d (%s)\n",
              nn, sem_ilha, R, conexas,
              if (sem_ilha > 0L) sprintf("%.1f%%", 100 * conexas / sem_ilha) else "sem casos"))
  linhas[[length(linhas) + 1L]] <- data.frame(
    n = nn, sorteios = R, sem_ilha = sem_ilha, conexas = conexas,
    prop_conexas = if (sem_ilha > 0L) conexas / sem_ilha else NA_real_,
    componentes_mediana = if (length(nc)) median(nc) else NA_real_,
    componentes_max = if (length(nc)) max(nc) else NA_integer_)
}
tab <- do.call(rbind, linhas)
write.csv(tab, file.path(AQUI, "saida", "conectividade_subamostras.csv"), row.names = FALSE)
cat("\nGravado em saida/conectividade_subamostras.csv\n")
