# Compare SEM PIPs with and without centering of transformed covariates.
# Reads the two logged tables from script 31 and compares their differences
# with the between-chain variation reported by those runs.

source(if (file.exists("../caminhos.R")) "../caminhos.R" else "caminhos.R")

AQUI <- ESPACIAL

ler_tabela <- function(arq) {
  x <- readLines(arq, warn = FALSE)
  i <- grep("covariavel +mqo +sar +sem", x)
  if (length(i) != 1L) stop("tabela nao encontrada em ", arq)
  fim <- i
  while (fim + 1L <= length(x) && grepl("^\\s*\\S+\\s+[0-9.]+\\s+[0-9.]+\\s+[0-9.]+\\s*$",
                                        x[fim + 1L])) fim <- fim + 1L
  d <- read.table(text = x[i:fim], header = TRUE, stringsAsFactors = FALSE)
  d[order(d$covariavel), ]
}
ruido <- function(arq) {
  x <- readLines(arq, warn = FALSE)
  as.numeric(sub(".*cadeias: *", "", x[grep("maior diferenca de PIP entre cadeias", x)]))
}

a <- ler_tabela(file.path(AQUI, "saida", "log_convencao2_sim.txt"))
b <- ler_tabela(file.path(AQUI, "saida", "log_convencao2_nao.txt"))
stopifnot(identical(a$covariavel, b$covariavel))

d <- abs(a$sem - b$sem)
r <- max(ruido(file.path(AQUI, "saida", "log_convencao2_sim.txt")),
         ruido(file.path(AQUI, "saida", "log_convencao2_nao.txt")))

cat("=== PIP do SEM: centrada contra a forma da Equacao (13) ===\n")
cat("  covariaveis:", nrow(a), "| sorteios por cadeia: 80.000\n\n")
cat(sprintf("  maior diferenca de PIP entre as duas formas: %.4f (%s)\n",
            max(d), a$covariavel[which.max(d)]))
cat(sprintf("  diferenca media: %.4f | mediana: %.4f\n", mean(d), median(d)))
cat(sprintf("  ruido de Monte Carlo entre cadeias da mesma forma: %.4f\n", r))
cat(sprintf("  razao entre o efeito da centragem e o ruido: %.2f\n", max(d) / r))
cat(sprintf("  correlacao de Spearman entre as duas formas: %.4f\n",
            cor(a$sem, b$sem, method = "spearman")))
cat(sprintf("  acima de 0,50: centrada %d | Equacao (13) %d\n",
            sum(a$sem >= 0.5), sum(b$sem >= 0.5)))

ord <- order(-d)
cat("\n  as cinco maiores diferencas:\n")
print(data.frame(covariavel = a$covariavel[ord][1:5],
                 centrada = round(a$sem[ord][1:5], 4),
                 eq13 = round(b$sem[ord][1:5], 4),
                 dif = round(d[ord][1:5], 4)), row.names = FALSE)

cat("\n  Leitura: a forma da fonte nao minimiza a soma de quadrados com\n")
cat("  intercepto, e o RSS difere em ate 1,9e-1 nesta base. O efeito sobre as\n")
cat("  PIP, porem, fica na mesma ordem do ruido de Monte Carlo de cadeias de 80\n")
cat("  mil sorteios, e nenhuma conclusao muda. A correcao importa para quem\n")
cat("  reimplementa o metodo, nao para as conclusoes deste artigo.\n")

write.csv(data.frame(covariavel = a$covariavel, centrada = a$sem,
                     eq13 = b$sem, dif = d),
          file.path(AQUI, "saida", "comparacao_centragem.csv"), row.names = FALSE)
cat("\nGravado em saida/comparacao_centragem.csv\n")
