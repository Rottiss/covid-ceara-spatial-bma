# Figuras do artigo. Graficos de base, saida em PDF vetorial em saida/figuras/.
#
# Figura 1  PIP sob MQO contra PIP sob SAR, um painel por desfecho. E a figura
#           que carrega o resultado central: os tres modos de falha aparecem na
#           forma da nuvem, nao na tabela.
# Figura 2  Razao entre o D do SAR e o do MQO nos cenarios com processo gerador conhecido
#           e nos dados reais. Separa os dois mundos numa linha so.
# Figura 3  Posteriori de rho sob media de modelos (quando 33 ja tiver rodado).

source(if (file.exists("../caminhos.R")) "../caminhos.R" else "caminhos.R")

AQUI <- ESPACIAL
FIG <- file.path(AQUI, "saida", "figuras")
dir.create(FIG, showWarnings = FALSE, recursive = TRUE)

D <- function(p) { m <- mean(p); K <- length(p)
  sd(p) * sqrt((K - 1) / K) / sqrt(m * (1 - m)) }

cinza <- "grey45"; escuro <- "grey15"

# ------------------------------------------------------------------ figura 1
desfechos <- c(casos = "casos", obitos = "óbitos", letalidade = "letalidade")
dados <- lapply(names(desfechos), function(d)
  readRDS(file.path(AQUI, "saida", "grade", paste0(d, "_queen.rds"))))
names(dados) <- names(desfechos)

pdf(file.path(FIG, "fig1_pip_mqo_sar.pdf"), width = 9, height = 3.4,
    pointsize = 10)
par(mfrow = c(1, 3), mar = c(4, 4, 3, 1), mgp = c(2.4, 0.8, 0), las = 1)
for (d in names(desfechos)) {
  z <- dados[[d]]
  x <- z$res$MQO$pip; y <- z$res$SAR$pip
  plot(x, y, xlim = c(0, 1), ylim = c(0, 1), pch = 21, bg = "white",
       col = escuro, cex = 1.1, xlab = "PIP sob MQO", ylab = "PIP sob SAR",
       panel.first = {
         abline(0, 1, col = "grey80")
         abline(h = 0.5, v = 0.5, col = "grey85", lty = 3)
         abline(h = z$res$SAR$tam / z$K, col = "grey60", lty = 2)
       })
  # A linha tracejada e a PIP media sob SAR: o nivel que uma covariavel sem
  # informacao alguma alcanca so por causa do tamanho medio do modelo.
  sobe <- which(y > x + 0.05)
  if (length(sobe)) points(x[sobe], y[sobe], pch = 19, col = escuro, cex = 1.1)
  title(main = sprintf("%s  (I de Moran %.3f)", desfechos[[d]], z$mi_y),
        font.main = 1)
  mtext(sprintf("D: %.3f para %.3f", D(x), D(y)), side = 3, line = -1.2,
        adj = 0.04, cex = 0.75, col = cinza)
}
dev.off()
cat("figura 1 gravada\n")

# ------------------------------------------------------------------ figura 2
linhas <- list()
sim <- readRDS(file.path(AQUI, "saida", "experimento_simulacao.rds"))
for (nm in names(sim)) {
  r <- sim[[nm]]
  linhas[[length(linhas) + 1]] <- data.frame(
    rotulo = sprintf("simulação: %s", nm), tipo = "simulação",
    razao = mean(apply(r$sar, 1, D)) / mean(apply(r$ols, 1, D)),
    stringsAsFactors = FALSE)
}
for (arq in c("cenario_realista.rds", "cenario_realista_oito.rds",
              "cenario_rho070_quatro.rds", "cenario_rho070_oito.rds")) {
  cam <- file.path(AQUI, "saida", arq)
  if (!file.exists(cam)) next
  cr <- readRDS(cam)
  for (nm in names(cr$resultados)) {
    r <- cr$resultados[[nm]]
    linhas[[length(linhas) + 1]] <- data.frame(
      rotulo = sprintf("realista: %s verd., rho %.2f", nm, cr$rho),
      tipo = "simulação",
      razao = mean(apply(r$sar, 1, D)) / mean(apply(r$ols, 1, D)),
      stringsAsFactors = FALSE)
  }
}
for (d in names(desfechos)) {
  z <- dados[[d]]
  linhas[[length(linhas) + 1]] <- data.frame(
    rotulo = sprintf("dados reais: %s", desfechos[[d]]), tipo = "dados",
    razao = D(z$res$SAR$pip) / D(z$res$MQO$pip), stringsAsFactors = FALSE)
}
tab <- do.call(rbind, linhas)
tab <- tab[order(tab$tipo, tab$razao), ]

pdf(file.path(FIG, "fig2_razao_D.pdf"), width = 7, height = 0.32 * nrow(tab) + 1.2,
    pointsize = 10)
par(mar = c(4, 15, 1, 1), mgp = c(2.4, 0.8, 0), las = 1)
plot(tab$razao, seq_len(nrow(tab)), xlim = c(0, 1.1), yaxt = "n", bty = "n",
     xlab = "razão entre o D do SAR e o do MQO", ylab = "",
     pch = ifelse(tab$tipo == "dados", 19, 21),
     bg = "white", col = escuro, cex = 1.2,
     panel.first = {
       abline(v = seq(0, 1, 0.25), col = "grey90")
       segments(0, seq_len(nrow(tab)), tab$razao, seq_len(nrow(tab)),
                col = "grey80")
     })
axis(2, at = seq_len(nrow(tab)), labels = tab$rotulo, tick = FALSE,
     cex.axis = 0.8)
dev.off()
cat("figura 2 gravada com", nrow(tab), "linhas\n")
write.csv(tab, file.path(AQUI, "saida", "razoes_D.csv"), row.names = FALSE)

# ------------------------------------------------------------------ figura 3
# Os três desfechos no mesmo painel: é a figura que sustenta a Seção 4.5, porque
# mostra que casos e óbitos, de I de Moran quase igual, têm dependência residual
# bem diferente.
arqs <- file.path(AQUI, "saida",
                  sprintf("especificacao_%s_queen.rds", names(desfechos)))
tem <- file.exists(arqs)
if (any(tem)) {
  es <- lapply(arqs[tem], readRDS)
  nomes <- unlist(desfechos)[tem]
  ws <- lapply(es, function(e) { w <- (e$a$w_rho + e$b$w_rho) / 2; w / sum(w) })
  passo <- es[[1]]$rho[2] - es[[1]]$rho[1]
  ymax <- max(vapply(ws, function(w) max(w / passo), numeric(1)))
  ltys <- c(1, 2, 4)

  pdf(file.path(FIG, "fig3_posteriori_rho.pdf"), width = 6, height = 3.8,
      pointsize = 10)
  par(mar = c(4, 4, 1, 1), mgp = c(2.4, 0.8, 0), las = 1)
  plot(NA, xlim = c(0, 1), ylim = c(0, ymax * 1.08), bty = "n",
       xlab = expression(rho), ylab = "densidade a posteriori")
  legenda <- character(0)
  for (i in seq_along(es)) {
    lines(es[[i]]$rho, ws[[i]] / passo, col = escuro, lwd = 2, lty = ltys[i])
    m <- sum(ws[[i]] * es[[i]]$rho)
    p_sar <- (es[[i]]$a$p_sar + es[[i]]$b$p_sar) / 2
    legenda[i] <- sprintf("%s: média %.3f, P(SAR | y) = %.3f",
                          nomes[i], m, p_sar)
  }
  legend("topleft", legend = legenda, lty = ltys[seq_along(es)], lwd = 2,
         col = escuro, bty = "n", cex = 0.8)
  dev.off()
  cat("figura 3 gravada com", length(es), "desfechos\n")
} else cat("figura 3 pendente: rode 33_probabilidade_espacial.R\n")

cat("\nSaida em", FIG, "\n")
