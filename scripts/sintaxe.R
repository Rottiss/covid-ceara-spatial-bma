source("../caminhos.R")
arquivos <- c(list.files(RAIZ, pattern = "[.]R$", full.names = TRUE),
              list.files(file.path(RAIZ, "espacial"), pattern = "[.]R$", full.names = TRUE),
              list.files(ESTIMACOES, pattern = "[.]R$", full.names = TRUE),
              list.files(file.path(RAIZ, "scripts"), pattern = "[.]R$", full.names = TRUE))
for (f in arquivos) parse(file = f, encoding = "UTF-8")
cat("R syntax checked: ", length(arquivos), " files.\n", sep = "")
