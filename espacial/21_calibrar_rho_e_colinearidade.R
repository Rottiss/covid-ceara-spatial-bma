# Duas coisas que faltavam para o experimento ficar completo.
#
# (A) CALIBRAR rho PARA REPRODUZIR A ESTRUTURA ESPACIAL REAL. Nenhum cenario
# rodado ate agora atingiu o I de Moran do desfecho real (0,4834). Os cenarios
# com rho = 0 pararam em 0,13. O controle usou rho = 0,6 sem verificar onde isso
# cai. O teste mais proximo do caso real exige um DGP com I de Moran comparavel.
#
# (B) MEDIR A COLINEARIDADE. Na calibracao do sinal apareceu que covariaveis
# espacialmente suaves precisam do dobro de R2 para a mesma detectabilidade.
# Isso foi subproduto; aqui vira medida.

source(if (file.exists("../caminhos.R")) "../caminhos.R" else "caminhos.R")

suppressPackageStartupMessages(library(spdep))

AQUI <- ESPACIAL
base <- readRDS(file.path(AQUI, "saida", "base_sar.rds"))
des  <- readRDS(file.path(AQUI, "saida", "desenho_simulacao.rds"))
cal  <- readRDS(file.path(AQUI, "saida", "calibracao_sinal.rds"))

X <- base$X; W <- base$W; n <- base$n; covariaveis <- base$covariaveis
listw <- base$listw; y_real <- base$y
Xc <- sweep(X, 2, colMeans(X), "-")
mi_real <- spdep::moran.test(y_real, listw)$estimate[["Moran I statistic"]]

# =================== (A) rho que reproduz o I de Moran real ===================
cat("=== (A) Calibracao de rho ===\n")
cat("  I de Moran do desfecho real:", round(mi_real, 4), "\n\n")

idx_s <- des$config$suaves$idx
b_s <- cal$suaves$beta_escolhido
beta_s <- rep(b_s, length(idx_s))

mi_para_rho <- function(rho, reps = 30L, semente = 4242L) {
  S_inv <- solve(diag(n) - rho * W)
  set.seed(semente)
  mean(replicate(reps, {
    y <- as.numeric(S_inv %*% (as.numeric(Xc[, idx_s, drop = FALSE] %*% beta_s) + rnorm(n)))
    spdep::moran.test(y, listw)$estimate[["Moran I statistic"]]
  }))
}

cat("  rho     I de Moran de y   R2 do modelo verdadeiro\n")
r2_para_rho <- function(rho, reps = 30L, semente = 4242L) {
  S_inv <- solve(diag(n) - rho * W)
  set.seed(semente)
  mean(replicate(reps, {
    y <- as.numeric(S_inv %*% (as.numeric(Xc[, idx_s, drop = FALSE] %*% beta_s) + rnorm(n)))
    yc <- y - mean(y)
    Q <- qr.Q(qr(Xc[, idx_s, drop = FALSE]))
    sum((crossprod(Q, yc))^2) / sum(yc^2)
  }))
}

grade_rho <- c(0, 0.2, 0.4, 0.5, 0.6, 0.7, 0.8)
mis <- numeric(length(grade_rho)); r2s <- numeric(length(grade_rho))
for (i in seq_along(grade_rho)) {
  mis[i] <- mi_para_rho(grade_rho[i]); r2s[i] <- r2_para_rho(grade_rho[i])
  cat(sprintf("  %.2f    %13.4f   %19.4f\n", grade_rho[i], mis[i], r2s[i]))
}

# Interpola para achar o rho que da o I de Moran real.
rho_alvo <- approx(mis, grade_rho, xout = mi_real)$y
cat("\n  rho que reproduz I de Moran de", round(mi_real, 4), ":",
    round(rho_alvo, 4), "\n")
cat("  conferido:", round(mi_para_rho(rho_alvo), 4), "\n\n")

# O beta precisa ser reduzido: o multiplicador espacial amplifica o sinal e o R2
# do modelo verdadeiro sobe. Recalibra para manter R2 comparavel ao caso rho = 0.
r2_no_alvo <- r2_para_rho(rho_alvo)
cat("  R2 do modelo verdadeiro em rho =", round(rho_alvo, 3), ":", round(r2_no_alvo, 4), "\n")
cat("  (em rho = 0 com o mesmo beta era", round(r2s[1], 4), ")\n\n")

# =========================== (B) colinearidade ================================
cat("=== (B) Colinearidade dos subconjuntos ===\n")

diagnostico <- function(vars, rotulo) {
  idx <- match(vars, covariaveis)
  M <- Xc[, idx, drop = FALSE]
  R <- cor(M)
  fora <- R[upper.tri(R)]
  cond <- kappa(scale(M), exact = TRUE)
  # VIF de cada uma contra as demais do proprio subconjunto
  vifs <- vapply(seq_along(idx), function(j) {
    r2 <- summary(lm(M[, j] ~ M[, -j, drop = FALSE]))$r.squared
    1 / (1 - r2)
  }, numeric(1))
  cat(sprintf("  %-8s |r| medio entre elas: %.4f | maximo: %.4f\n",
              rotulo, mean(abs(fora)), max(abs(fora))))
  cat(sprintf("           numero de condicao: %.3f | VIF maximo: %.3f\n", cond, max(vifs)))
  cat(sprintf("           I de Moran medio: %.4f\n", mean(des$moran[vars])))
  list(r_medio = mean(abs(fora)), r_max = max(abs(fora)),
       cond = cond, vif_max = max(vifs), mi = mean(des$moran[vars]))
}

d_s <- diagnostico(des$config$suaves$vars, "suaves")
d_r <- diagnostico(des$config$rugosas$vars, "rugosas")

cat("\n  As suaves sao mais colineares entre si? ",
    ifelse(d_s$r_medio > d_r$r_medio, "sim", "nao"), "\n", sep = "")

# A hipotese mais geral: correlacao entre covariaveis cresce com a semelhanca
# de estrutura espacial. Testa sobre todos os pares.
cat("\n=== Correlacao entre covariaveis vs. semelhanca de I de Moran ===\n")
cont <- covariaveis[grepl("^z_", covariaveis)]
ic <- match(cont, covariaveis)
R <- cor(Xc[, ic, drop = FALSE])
mi_c <- des$moran[cont]
pares <- which(upper.tri(R), arr.ind = TRUE)
tab <- data.frame(
  r_abs = abs(R[upper.tri(R)]),
  mi_min = pmin(mi_c[pares[, 1]], mi_c[pares[, 2]]),
  mi_media = (mi_c[pares[, 1]] + mi_c[pares[, 2]]) / 2
)
cat("  pares:", nrow(tab), "\n")
cat("  correlacao entre |r| do par e o I de Moran medio do par:",
    round(cor(tab$r_abs, tab$mi_media), 4), "\n")
cat("  |r| medio quando ambos tem I de Moran > 0,30:",
    round(mean(tab$r_abs[tab$mi_min > 0.30]), 4),
    "(", sum(tab$mi_min > 0.30), "pares )\n")
cat("  |r| medio quando ambos tem I de Moran < 0,10:",
    round(mean(tab$r_abs[tab$mi_min < 0.10]), 4),
    "(", sum(tab$mi_min < 0.10), "pares )\n")

saveRDS(list(rho_alvo = rho_alvo, mi_real = mi_real, grade_rho = grade_rho,
             mis = mis, r2s = r2s, r2_no_alvo = r2_no_alvo,
             colinearidade = list(suaves = d_s, rugosas = d_r), pares = tab),
        file.path(AQUI, "saida", "calibracao_rho.rds"))
cat("\nGravado em saida/calibracao_rho.rds\n")
