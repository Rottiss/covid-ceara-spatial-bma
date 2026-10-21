source(if (file.exists("../caminhos.R")) "../caminhos.R" else "caminhos.R")

tb <- readRDS(file.path(ESPACIAL, "saida/tabela_f21_estendida.rds"))
TAB <- tb$TAB; N <- tb$N_TAB
set.seed(1)
i <- sample(N - 1L, 201L); k <- 7L

t0 <- proc.time()[3]
for (r in 1:2000) { col <- TAB[, k + 1L]; v <- col[i] * 0.5 + col[i + 1L] * 0.5 }
t_col <- (proc.time()[3] - t0) / 2000

off <- k * N
t0 <- proc.time()[3]
for (r in 1:2000) { v2 <- TAB[i + off] * 0.5 + TAB[i + 1L + off] * 0.5 }
t_dir <- (proc.time()[3] - t0) / 2000

cat(sprintf("extrair a coluna inteira: %.4f ms\n", t_col * 1000))
cat(sprintf("indexar direto          : %.4f ms\n", t_dir * 1000))
cat(sprintf("ganho: %.0f vezes\n", t_col / t_dir))
cat("resultados identicos:", isTRUE(all.equal(v, v2)), "\n")
