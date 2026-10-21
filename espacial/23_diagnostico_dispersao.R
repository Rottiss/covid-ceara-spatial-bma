# Summarize normalized PIP dispersion in observed and simulated data.
# For mean PIP p, D divides the population standard deviation of the PIPs
# by sqrt(p * (1-p)). Equal PIPs give D=0; PIPs restricted to 0/1 give D=1
# when 0<p<1. The sum of PIPs is reported separately as mean model size.
# D is descriptive; it does not by itself establish selection accuracy,
# information loss, or a failure of the sampler.

source(if (file.exists("../caminhos.R")) "../caminhos.R" else "caminhos.R")

AQUI <- ESPACIAL

D <- function(pip) {
  p <- mean(pip)
  if (p <= 0 || p >= 1) return(NA_real_)
  sd(pip) * sqrt((length(pip) - 1) / length(pip)) / sqrt(p * (1 - p))
}

resumo <- function(pip, rotulo) {
  p <- mean(pip); K <- length(pip)
  cat(sprintf("  %-28s k_medio %6.3f | p %.4f | dp %.4f | dp_max %.4f | D %.4f\n",
              rotulo, sum(pip), p, sd(pip) * sqrt((K - 1) / K), sqrt(p * (1 - p)), D(pip)))
  D(pip)
}

cat("=== Dados reais, casos com Queen, 1.000.000 de iteracoes ===\n")
real <- read.csv(file.path(AQUI, "saida", "comparacao_1milhao.csv"),
                 stringsAsFactors = FALSE)
d_pub <- resumo(real$publicada, "TCC publicado (MQO)")
d_mqo <- resumo(real$mqo, "MQO reimplementado")
d_sar <- resumo(real$sar, "SAR")
cat(sprintf("\n  queda de D entre MQO e SAR: %.4f para %.4f, ou %.0f%% do valor original\n\n",
            d_mqo, d_sar, 100 * d_sar / d_mqo))

cat("=== Simulacao com processo gerador conhecido ===\n")
cat("  Referencia: nos cenarios simulados sabemos que 4 de 30 sao verdadeiras.\n")
cat("  D descreve a dispersao das PIP; nao mede, isoladamente, a selecao correta\n")
cat("  das covariaveis com coeficientes diferentes de zero.\n\n")

sim <- readRDS(file.path(AQUI, "saida", "experimento_simulacao.rds"))
for (nome in names(sim)) {
  rr <- sim[[nome]]
  d_o <- mean(apply(rr$ols, 1, D)); d_s <- mean(apply(rr$sar, 1, D))
  cat(sprintf("  %-10s (rho do DGP = %.1f)  D no MQO: %.4f | D no SAR: %.4f | razao: %.3f\n",
              nome, rr$rho, d_o, d_s, d_s / d_o))
}

arq_real <- file.path(AQUI, "saida", "cenario_realista.rds")
if (file.exists(arq_real)) {
  cat("\n=== Cenario realista (rho calibrado ao I de Moran observado) ===\n")
  cr <- readRDS(arq_real)
  for (nome in names(cr$resultados)) {
    rr <- cr$resultados[[nome]]
    d_o <- mean(apply(rr$ols, 1, D)); d_s <- mean(apply(rr$sar, 1, D))
    cat(sprintf("  %-8s (%d verdadeiras)  D no MQO: %.4f | D no SAR: %.4f | razao: %.3f\n",
                nome, length(rr$idx), d_o, d_s, d_s / d_o))
  }
}

pasta <- file.path(AQUI, "saida", "grade")
arqs <- list.files(pasta, pattern = "\\.rds$", full.names = TRUE)
if (length(arqs) > 0) {
  cat("\n=== Grade empirica ===\n")
  linhas <- lapply(arqs, function(a) {
    g <- readRDS(a)
    data.frame(desfecho = g$desfecho, vizinhos = g$vizinhos,
               mi_y = g$mi_y,
               k_mqo = sum(g$res$MQO$pip), k_sar = sum(g$res$SAR$pip),
               D_mqo = D(g$res$MQO$pip), D_sar = D(g$res$SAR$pip),
               acima_mqo = sum(g$res$MQO$pip >= 0.5),
               acima_sar = sum(g$res$SAR$pip >= 0.5),
               stringsAsFactors = FALSE)
  })
  tabela <- do.call(rbind, linhas)
  tabela$razao_D <- tabela$D_sar / tabela$D_mqo
  print(transform(tabela, mi_y = round(mi_y, 4),
                  k_mqo = round(k_mqo, 2), k_sar = round(k_sar, 2),
                  D_mqo = round(D_mqo, 4), D_sar = round(D_sar, 4),
                  razao_D = round(razao_D, 3)), row.names = FALSE)
  write.csv(tabela, file.path(AQUI, "saida", "diagnostico_dispersao.csv"),
            row.names = FALSE)
}

cat("\nD proximo de zero indica PIP semelhantes entre as covariaveis.\n")
cat("Isso nao implica concentracao posterior do tamanho do modelo.\n")
