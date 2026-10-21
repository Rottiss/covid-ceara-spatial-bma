source(if (file.exists("../caminhos.R")) "../caminhos.R" else "caminhos.R")

base <- readRDS(file.path(ESPACIAL, "saida/base_sar.rds"))
X <- base$X
medias <- colMeans(X)
nz <- which(abs(medias) > 1e-10)
cat("Colunas com media diferente de zero:\n")
for (i in nz) cat(sprintf("  %-16s media = %.4f\n", colnames(X)[i], medias[i]))
cat("\ntotal:", length(nz), "de", ncol(X), "\n")
