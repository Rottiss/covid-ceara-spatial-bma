# PIP com a especificacao integrada, e a relacao entre rho posterior e severidade.
#
# Em letalidade a posteriori da especificacao nao e degenerada: o SAR fica com
# dois tercos da massa. Nesse caso nem a PIP do MQO nem a do SAR e a resposta
# bayesiana; a resposta e a media das duas ponderada pela massa de cada
# especificacao, que e o que a cadeia de 33 produz diretamente.

source(if (file.exists("../caminhos.R")) "../caminhos.R" else "caminhos.R")

AQUI <- ESPACIAL
D <- function(p) { m <- mean(p); K <- length(p)
  sd(p) * sqrt((K - 1) / K) / sqrt(m * (1 - m)) }

desfechos <- c("casos", "obitos", "letalidade")
linhas <- list()

for (d in desfechos) {
  arq <- file.path(AQUI, "saida", sprintf("especificacao_%s_queen.rds", d))
  if (!file.exists(arq)) { cat("faltando:", d, "\n"); next }
  e <- readRDS(arq)
  g <- readRDS(file.path(AQUI, "saida", "grade", paste0(d, "_queen.rds")))

  p_sar <- (e$a$p_sar + e$b$p_sar) / 2
  w <- (e$a$w_rho + e$b$w_rho) / 2; w <- w / sum(w)
  rho_med <- sum(w * e$rho)
  ac <- cumsum(w)
  q <- e$rho[c(which.max(ac >= 0.025), which.max(ac >= 0.975))]

  # PIP integrando a especificacao. Quando a cadeia nunca visita o MQO, as
  # visitas ao MQO sao zero e a media coincide com a PIP do SAR.
  n_s <- e$a$n_sar + e$b$n_sar; n_o <- e$a$n_ols + e$b$n_ols
  pip_s <- if (n_s > 0) (e$a$pip_sar * e$a$n_sar + e$b$pip_sar * e$b$n_sar) / n_s else NA
  pip_o <- if (n_o > 0) (e$a$pip_ols * e$a$n_ols + e$b$pip_ols * e$b$n_ols) / n_o else NA
  pip_int <- if (n_o > 0) (pip_s * n_s + pip_o * n_o) / (n_s + n_o) else pip_s

  cat(sprintf("\n=== %s ===\n", d))
  cat(sprintf("  P(SAR | y) = %.4f | rho posterior %.3f [%.3f, %.3f]\n",
              p_sar, rho_med, q[1], q[2]))
  cat(sprintf("  D: MQO da grade %.4f | SAR da grade %.4f | razao %.3f\n",
              D(g$res$MQO$pip), D(g$res$SAR$pip),
              D(g$res$SAR$pip) / D(g$res$MQO$pip)))
  cat(sprintf("  acima de 0,50: MQO %d | SAR %d | especificacao integrada %d\n",
              sum(g$res$MQO$pip >= 0.5), sum(g$res$SAR$pip >= 0.5),
              sum(pip_int >= 0.5)))
  if (n_o > 0) {
    sel_mqo <- g$covariaveis <- e$covariaveis
    acima <- e$covariaveis[pip_int >= 0.5]
    cat("  quais, com a especificacao integrada:",
        if (length(acima)) paste(acima, collapse = ", ") else "nenhuma", "\n")
    cat(sprintf("  maior |PIP integrada - PIP do SAR|: %.4f\n",
                max(abs(pip_int - pip_s))))
  }
  linhas[[length(linhas) + 1]] <- data.frame(
    desfecho = d, mi_y = g$mi_y, p_sar = p_sar, rho = rho_med,
    rho_inf = q[1], rho_sup = q[2],
    D_mqo = D(g$res$MQO$pip), D_sar = D(g$res$SAR$pip),
    razao = D(g$res$SAR$pip) / D(g$res$MQO$pip),
    acima_mqo = sum(g$res$MQO$pip >= 0.5), acima_sar = sum(g$res$SAR$pip >= 0.5),
    acima_int = sum(pip_int >= 0.5), stringsAsFactors = FALSE)
}

tab <- do.call(rbind, linhas)
cat("\n=== O que governa a severidade ===\n")
print(transform(tab[, c("desfecho", "mi_y", "rho", "p_sar", "razao")],
                mi_y = round(mi_y, 4), rho = round(rho, 4),
                p_sar = round(p_sar, 4), razao = round(razao, 3)),
      row.names = FALSE)
if (nrow(tab) == 3) {
  cat("\n  Com tres pontos nao ha teste, so ordenacao. Ainda assim:\n")
  cat(sprintf("  I de Moran de y: %.3f, %.3f, %.3f (casos e obitos quase iguais)\n",
              tab$mi_y[1], tab$mi_y[2], tab$mi_y[3]))
  cat(sprintf("  rho posterior:   %.3f, %.3f, %.3f (casos e obitos bem diferentes)\n",
              tab$rho[1], tab$rho[2], tab$rho[3]))
  cat(sprintf("  razao de D:      %.3f, %.3f, %.3f\n",
              tab$razao[1], tab$razao[2], tab$razao[3]))
  cat("\n  O I de Moran de y mistura a suavidade que as covariaveis explicam com a\n")
  cat("  que sobra; rho e so a que sobra, e e a que compete com as covariaveis.\n")
}
write.csv(tab, file.path(AQUI, "saida", "severidade_rho.csv"), row.names = FALSE)
cat("\nGravado em saida/severidade_rho.csv\n")
