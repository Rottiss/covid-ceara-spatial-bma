# Queen contra Rook no modelo de erro espacial.
#
# POR QUE ISTO EXISTE. O 61 fechou a comparacao de vizinhanca para a escolha da
# especificacao e para rho. As PIP do SEM, citadas no artigo e na Figura 3,
# vinham so da Queen. Aqui as duas matrizes ficam lado a lado, nas tres
# variaveis dependentes.
#
# Exige as seis rodadas do 31: saida/grade/sem_<desfecho>_<vizinhos>.rds

source(if (file.exists("../caminhos.R")) "../caminhos.R" else "caminhos.R")

AQUI <- ESPACIAL

Dst <- function(p) { m <- mean(p); K <- length(p)
  sd(p) * sqrt((K - 1) / K) / sqrt(m * (1 - m)) }

desfechos <- c("casos", "obitos", "letalidade")
linhas <- list(); detalhe <- list()

for (d in desfechos) {
  q <- readRDS(file.path(AQUI, "saida", "grade", sprintf("sem_%s_queen.rds", d)))
  r <- readRDS(file.path(AQUI, "saida", "grade", sprintf("sem_%s_rook.rds", d)))
  stopifnot(identical(q$covariaveis, r$covariaveis))

  dif <- abs(q$pip - r$pip)
  # Ruido de Monte Carlo de cada configuracao: as duas cadeias da mesma rodada.
  ruido <- max(max(abs(q$a$pip - q$b$pip)), max(abs(r$a$pip - r$b$pip)))
  pior <- which.max(dif)

  linhas[[d]] <- data.frame(
    desfecho = d,
    pip_max_queen = max(q$pip), pip_max_rook = max(r$pip),
    acima_50_queen = sum(q$pip >= 0.5), acima_50_rook = sum(r$pip >= 0.5),
    d_queen = Dst(q$pip), d_rook = Dst(r$pip),
    tam_queen = q$tam, tam_rook = r$tam,
    dif_pip_max = max(dif), covariavel_dif_max = q$covariaveis[pior],
    ruido_entre_cadeias = ruido,
    spearman_queen_rook = cor(q$pip, r$pip, method = "spearman"),
    stringsAsFactors = FALSE)

  detalhe[[d]] <- data.frame(desfecho = d, covariavel = q$covariaveis,
                             pip_queen = q$pip, pip_rook = r$pip,
                             diferenca = q$pip - r$pip, stringsAsFactors = FALSE)
}

res <- do.call(rbind, linhas); rownames(res) <- NULL
det <- do.call(rbind, detalhe); rownames(det) <- NULL

write.csv(res, file.path(AQUI, "saida", "comparacao_vizinhanca_sem.csv"), row.names = FALSE)
write.csv(det, file.path(AQUI, "saida", "comparacao_vizinhanca_sem_pip.csv"), row.names = FALSE)

cat("=== SEM: Queen contra Rook ===\n\n")
print(transform(res, pip_max_queen = round(pip_max_queen, 4),
                pip_max_rook = round(pip_max_rook, 4),
                d_queen = round(d_queen, 4), d_rook = round(d_rook, 4),
                tam_queen = round(tam_queen, 2), tam_rook = round(tam_rook, 2),
                dif_pip_max = round(dif_pip_max, 4),
                ruido_entre_cadeias = round(ruido_entre_cadeias, 4),
                spearman_queen_rook = round(spearman_queen_rook, 4)),
      row.names = FALSE)

cat("\n  maior diferenca de PIP entre matrizes:", round(max(res$dif_pip_max), 4), "\n")
cat("  maior ruido entre cadeias da mesma configuracao:",
    round(max(res$ruido_entre_cadeias), 4), "\n")
cat("  a diferenca entre matrizes excede o ruido de Monte Carlo em",
    sum(res$dif_pip_max > res$ruido_entre_cadeias), "de 3 variaveis\n")

cat("\n=== Cinco maiores diferencas em cada variavel ===\n")
for (d in desfechos) {
  x <- detalhe[[d]]; x <- x[order(-abs(x$diferenca)), ][1:5, ]
  cat("\n", d, "\n")
  print(transform(x[, c("covariavel", "pip_queen", "pip_rook", "diferenca")],
                  pip_queen = round(pip_queen, 4), pip_rook = round(pip_rook, 4),
                  diferenca = round(diferenca, 4)), row.names = FALSE)
}

cat("\nGravado em saida/comparacao_vizinhanca_sem.csv e _pip.csv\n")
