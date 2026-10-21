source(if (file.exists("../caminhos.R")) "../caminhos.R" else "caminhos.R")

cr <- readRDS(file.path(ESPACIAL, "saida/cenario_realista.rds"))
cat("conjuntos salvos:", paste(names(cr$resultados), collapse = ", "), "\n")
for (nm in names(cr$resultados)) {
  r <- cr$resultados[[nm]]
  cat(sprintf("  %-8s replicacoes %d | verdadeiras %d | NA em ols: %d\n",
              nm, nrow(r$ols), length(r$idx), sum(is.na(r$ols))))
}
