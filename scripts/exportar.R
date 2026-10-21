source("../caminhos.R")
dest <- file.path(RAIZ, "execucoes_locais", "resumos")
dir.create(dest, recursive = TRUE, showWarnings = FALSE)
prob <- list(); pip <- list()
for (d in c("casos", "obitos", "letalidade")) for (w in c("queen", "rook")) {
  x <- readRDS(file.path(ESPACIAL, "saida/tres_vias", paste0(d, "_", w, ".rds")))
  prob[[length(prob) + 1L]] <- data.frame(variavel_dependente = d, vizinhanca = w,
    especificacao = names(x$p_esp), probabilidade = as.numeric(x$p_esp),
    cadeia_a = as.numeric(x$p_esp_a), cadeia_b = as.numeric(x$p_esp_b))
  pip[[length(pip) + 1L]] <- cbind(variavel_dependente = d, vizinhanca = w, x$tabela)
}
write.csv(do.call(rbind, prob), file.path(dest, "probabilidades.csv"), row.names = FALSE)
write.csv(do.call(rbind, pip), file.path(dest, "pips_condicionais.csv"), row.names = FALSE)
cat("Summaries exported to ", dest, "\n", sep = "")
