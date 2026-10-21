# Confere cada numero do artigo contra o arquivo de saida que o produziu.
#
# Os valores esperados sao transcritos do ARTIGO_RASCUNHO.md a mao. O script
# recalcula tudo dos .rds e .csv e aponta divergencia acima da tolerancia. E a
# defesa contra o erro mais provavel num texto com muitos numeros: o valor certo
# em execucao antiga que ficou no texto depois de a execucao ser refeita.

source(if (file.exists("../caminhos.R")) "../caminhos.R" else "caminhos.R")

AQUI <- ESPACIAL
D <- function(p) { m <- mean(p); K <- length(p)
  sd(p) * sqrt((K - 1) / K) / sqrt(m * (1 - m)) }

falhas <- 0L; testes <- 0L
conf <- function(rotulo, obtido, esperado, tol = 5e-3) {
  testes <<- testes + 1L
  ok <- is.finite(obtido) && abs(obtido - esperado) <= tol
  if (!ok) falhas <<- falhas + 1L
  cat(sprintf("  [%s] %-46s artigo %10.4f | recalculado %10.4f\n",
              if (ok) "ok" else "XX", rotulo, esperado, obtido))
}

cat("=== 4.1 Grade empirica ===\n")
grade <- list(
  list("casos", "queen", 13.10, 7.19, 0.3235, 0.0973, 6, 0, 0.759),
  list("casos", "rook", 13.10, 8.26, 0.3235, 0.1058, 6, 0, 0.778),
  list("obitos", "queen", 5.83, 6.10, 0.4652, 0.3299, 2, 2, 0.822),
  list("obitos", "rook", 5.83, 6.04, 0.4652, 0.3340, 2, 2, 0.828),
  list("letalidade", "queen", 15.05, 12.29, 0.2000, 0.1413, 12, 3, 0.979),
  list("letalidade", "rook", 15.05, 12.77, 0.2000, 0.1494, 12, 5, 0.982))
for (g in grade) {
  z <- readRDS(file.path(AQUI, "saida", "grade",
                         sprintf("%s_%s.rds", g[[1]], g[[2]])))
  pre <- sprintf("%s/%s", g[[1]], substr(g[[2]], 1, 1))
  conf(paste(pre, "k MQO"), z$res$MQO$tam, g[[3]], 0.01)
  conf(paste(pre, "k SAR"), z$res$SAR$tam, g[[4]], 0.01)
  conf(paste(pre, "D MQO"), D(z$res$MQO$pip), g[[5]], 1e-4)
  conf(paste(pre, "D SAR"), D(z$res$SAR$pip), g[[6]], 1e-4)
  conf(paste(pre, "n >= 0,50 MQO"), sum(z$res$MQO$pip >= 0.5), g[[7]], 0)
  conf(paste(pre, "n >= 0,50 SAR"), sum(z$res$SAR$pip >= 0.5), g[[8]], 0)
  conf(paste(pre, "Spearman"),
       cor(z$res$MQO$pip, z$res$SAR$pip, method = "spearman"), g[[9]], 1e-3)
}

cat("\n=== 4.2 PIP citadas no texto ===\n")
pip <- function(d, v) {
  z <- readRDS(file.path(AQUI, "saida", "grade", sprintf("%s_queen.rds", d)))
  cv <- z$tabela$covariavel
  c(mqo = z$tabela$mqo[cv == v], sar = z$tabela$sar[cv == v])
}
esperado_42 <- list(
  list("casos", "z_idosos", 0.9451, 0.3002), list("casos", "z_leitos1k", 0.8269, 0.3800),
  list("casos", "metrop", 0.6810, 0.2933), list("casos", "semiarido", 0.6287, 0.2337),
  list("casos", "z_ideb.fund", 0.6223, 0.2147), list("casos", "z_int.circ", 0.5741, 0.2245),
  list("obitos", "z_idosos", 0.9754, 0.4221), list("obitos", "z_ex.pobr", 0.5055, 0.2898),
  list("obitos", "z_pop.rural", 0.4112, 0.7243), list("obitos", "z_eleitas.fem", 0.3800, 0.5602))
for (e in esperado_42) {
  p <- pip(e[[1]], e[[2]])
  conf(sprintf("%s/%s MQO", e[[1]], e[[2]]), p[["mqo"]], e[[3]], 1e-4)
  conf(sprintf("%s/%s SAR", e[[1]], e[[2]]), p[["sar"]], e[[4]], 1e-4)
}

cat("\n=== 4.4 e 4.5 Especificacao e rho ===\n")
esp <- list(list("casos", 1.0000, 0.708), list("obitos", 1.0000, 0.555),
            list("letalidade", 0.665, 0.254))
for (e in esp) {
  z <- readRDS(file.path(AQUI, "saida", sprintf("especificacao_%s_queen.rds", e[[1]])))
  w <- (z$a$w_rho + z$b$w_rho) / 2; w <- w / sum(w)
  conf(sprintf("%s P(SAR | y)", e[[1]]), (z$a$p_sar + z$b$p_sar) / 2, e[[2]], 2e-3)
  conf(sprintf("%s rho posterior", e[[1]]), sum(w * z$rho), e[[3]], 1e-3)
}

# Fatores de Bayes condicionais em 200 conjuntos sorteados (59, antes só na tela do 33).
bfc <- read.csv(file.path(AQUI, "saida", "bf_condicional_amostra.csv"))
b1 <- function(d, col) bfc[bfc$desfecho == d, col]
conf("casos: 200 conjuntos, mediana do log BF", b1("casos", "log_bf_mediana"), 33.4, 0.05)
conf("casos: 200 conjuntos, minimo do log BF", b1("casos", "log_bf_minimo"), 29.5, 0.05)
conf("casos: BF condicional mediano, log10", log10(b1("casos", "bf_mediano")), log10(3.3e14), 0.007)
conf("obitos: BF condicional mediano, log10", log10(b1("obitos", "bf_mediano")), log10(1.3e8), 0.017)
conf("200 conjuntos: log BF > 0 em todos", as.numeric(min(bfc$log_bf_minimo) > 0), 1, 0)
conf("200 conjuntos: valores fora da tabela", sum(bfc$fora_tabela), 0, 0)

# 4.4: a mesma comparacao com a matriz Rook (33 e 60, consolidados pelo 61).
cv <- read.csv(file.path(AQUI, "saida", "comparacao_vizinhanca_especificacao.csv"),
               stringsAsFactors = FALSE)
v1 <- function(d, viz, col) cv[cv$desfecho == d & cv$vizinhos == viz, col]
for (e in list(list("casos", 0.679, 30.03), list("obitos", 0.542, 16.33),
               list("letalidade", 0.245, 0.74))) {
  d <- e[[1]]
  conf(paste(d, "rho posterior Rook"), v1(d, "rook", "rho_media"), e[[2]], 6e-4)
  conf(paste(d, "log BF agregado Rook"), v1(d, "rook", "log_bf_agregado"), e[[3]], 6e-3)
}
conf("casos e obitos: zero visitas ao MQO com Rook",
     sum(v1("casos", "rook", "visitas_mqo"), v1("obitos", "rook", "visitas_mqo")), 0, 0)
conf("letalidade: frequencia SAR com Rook", v1("letalidade", "rook", "p_sar_cadeia"), 0.674, 6e-4)
conf("maior diferenca de log BF entre matrizes",
     max(abs(cv$log_bf_agregado[cv$vizinhos == "queen"] - cv$log_bf_agregado[cv$vizinhos == "rook"])), 1.913, 6e-3)

cat("\n=== 4.6 Sensibilidade a priori ===\n")
pri <- list(
  list("uniforme_hyperg", 14.68, 14.80, 0.3667, 0.2074, 9, 8),
  list("binbeta5_hyperg", 6.75, 1.31, 0.4253, 0.0828, 2, 0),
  list("binomial5_hyperg", 5.41, 4.59, 0.4894, 0.2058, 2, 0),
  list("binbeta_uip", 1.88, 0.70, 0.7007, 0.1996, 1, 0))
for (p in pri) {
  z <- readRDS(file.path(AQUI, "saida", "priori",
                         sprintf("casos_queen_%s.rds", p[[1]])))
  conf(paste(p[[1]], "k MQO"), z$res$MQO$tam, p[[2]], 0.01)
  conf(paste(p[[1]], "k SAR"), z$res$SAR$tam, p[[3]], 0.01)
  conf(paste(p[[1]], "D MQO"), z$D_mqo, p[[4]], 1e-4)
  conf(paste(p[[1]], "D SAR"), z$D_sar, p[[5]], 1e-4)
  conf(paste(p[[1]], "n >= 0,50 MQO"), sum(z$res$MQO$pip >= 0.5), p[[6]], 0)
  conf(paste(p[[1]], "n >= 0,50 SAR"), sum(z$res$SAR$pip >= 0.5), p[[7]], 0)
}

cat("\n=== 4.7 SEM ===\n")
sem <- list(list("casos", 0.0988, 0, 4.79), list("obitos", 0.2343, 1, 5.14),
            list("letalidade", 0.1418, 2, 10.90))
for (s in sem) {
  z <- readRDS(file.path(AQUI, "saida", "grade", sprintf("sem_%s_queen.rds", s[[1]])))
  conf(sprintf("%s D SEM", s[[1]]), D(z$pip), s[[2]], 1e-4)
  conf(sprintf("%s n >= 0,50 SEM", s[[1]]), sum(z$pip >= 0.5), s[[3]], 0)
  conf(sprintf("%s k SEM", s[[1]]), z$tam, s[[4]], 0.01)
}
z <- readRDS(file.path(AQUI, "saida", "grade", "sem_casos_queen.rds"))
conf("z_idosos SEM casos", z$pip[z$covariaveis == "z_idosos"], 0.1408, 1e-4)
conf("linha de nao informacao SEM casos", z$tam / 30, 0.160, 1e-3)
z <- readRDS(file.path(AQUI, "saida", "grade", "sem_obitos_queen.rds"))
conf("z_idosos SEM obitos", z$pip[z$covariaveis == "z_idosos"], 0.1610, 1e-4)
conf("z_pop.rural SEM obitos", z$pip[z$covariaveis == "z_pop.rural"], 0.5224, 1e-4)
conf("z_eleitas.fem SEM obitos", z$pip[z$covariaveis == "z_eleitas.fem"], 0.3700, 1e-4)
conf("linha de nao informacao SEM obitos", z$tam / 30, 0.171, 1e-3)

# Queen contra Rook no SEM, do 62_comparar_vizinhanca_sem.R
cs <- read.csv(file.path(AQUI, "saida", "comparacao_vizinhanca_sem.csv"),
               stringsAsFactors = FALSE)
w1 <- function(d, col) cs[cs$desfecho == d, col]
for (e in list(list("casos", 0.1328, 0), list("obitos", 0.2408, 1),
               list("letalidade", 0.1472, 2))) {
  d <- e[[1]]
  conf(paste(d, "D SEM Rook"), w1(d, "d_rook"), e[[2]], 1e-4)
  conf(paste(d, "n >= 0,50 SEM Rook"), w1(d, "acima_50_rook"), e[[3]], 0)
}
conf("SEM: maior diferenca de PIP entre matrizes", max(cs$dif_pip_max), 0.127, 6e-4)
conf("SEM: segunda maior diferenca, casos", w1("casos", "dif_pip_max"), 0.107, 6e-4)
conf("SEM: maior ruido entre cadeias", max(cs$ruido_entre_cadeias), 0.023, 6e-4)
conf("SEM: menor Spearman Queen-Rook", min(cs$spearman_queen_rook), 0.960, 6e-4)
conf("SEM: maior Spearman Queen-Rook", max(cs$spearman_queen_rook), 0.987, 6e-4)
conf("SEM casos: k sob Rook", w1("casos", "tam_rook"), 5.79, 0.01)
conf("SEM obitos: pop rural sob Rook", w1("obitos", "pip_max_rook"), 0.5002, 1e-4)
conf("SEM obitos: ex.pobr Queen", {
  z <- read.csv(file.path(AQUI, "saida", "comparacao_vizinhanca_sem_pip.csv"),
                stringsAsFactors = FALSE)
  z$pip_queen[z$desfecho == "obitos" & z$covariavel == "z_ex.pobr"] }, 0.3528, 1e-4)
conf("SEM obitos: ex.pobr Rook", {
  z <- read.csv(file.path(AQUI, "saida", "comparacao_vizinhanca_sem_pip.csv"),
                stringsAsFactors = FALSE)
  z$pip_rook[z$desfecho == "obitos" & z$covariavel == "z_ex.pobr"] }, 0.4797, 1e-4)

cat("\n=== 5.3 Simulacoes ===\n")
razoes <- read.csv(file.path(AQUI, "saida", "razoes_D.csv"), stringsAsFactors = FALSE)
esperado_53 <- c("simulação: suaves" = 0.99, "simulação: rugosas" = 1.00,
                 "simulação: controle" = 0.88,
                 "realista: quatro verd., rho 0.53" = 0.86,
                 "realista: quatro verd., rho 0.70" = 0.75,
                 "realista: oito verd., rho 0.70" = 0.57)
for (nm in names(esperado_53)) {
  i <- which(razoes$rotulo == nm)
  if (length(i) != 1L) { cat("  [XX] rotulo ausente ou duplicado:", nm, "\n")
    falhas <- falhas + 1L; testes <- testes + 1L; next }
  conf(nm, razoes$razao[i], esperado_53[[nm]], 6e-3)
}
oito53 <- sort(razoes$razao[razoes$rotulo == "realista: oito verd., rho 0.53"])
conf("realista oito rho 0,53 (menor)", oito53[1], 0.74, 6e-3)
conf("realista oito rho 0,53 (maior)", oito53[2], 0.78, 6e-3)
conf("dados reais casos", razoes$razao[razoes$rotulo == "dados reais: casos"], 0.30, 6e-3)

cat("\n=== 5.6 Escala ===\n")
esc <- read.csv(file.path(AQUI, "saida", "escala_n.csv"), stringsAsFactors = FALSE)
conf("n=120 razao", esc$razao[esc$n == 120], 0.269, 1e-3)
conf("n=120 desvio", esc$dp_razao[esc$n == 120], 0.126, 1e-3)
conf("n=160 razao", esc$razao[esc$n == 160], 0.320, 1e-3)
conf("n=160 desvio", esc$dp_razao[esc$n == 160], 0.224, 1e-3)
conf("n=184 razao", esc$razao[esc$n == 184], 0.294, 1e-3)
conf("ligacoes n=184", esc$ligacoes[esc$n == 184], 5.27, 5e-3)
conf("ligacoes n=120", esc$ligacoes[esc$n == 120], 3.54, 5e-3)
r160 <- as.numeric(strsplit(esc$razoes[esc$n == 160], ";")[[1]])
conf("n=160 menor razao individual", min(r160), 0.128, 1e-3)
conf("n=160 maior razao individual", max(r160), 0.732, 1e-3)

cat("\n=== 5.3 D das simulacoes citado no texto ===\n")
ts <- read.csv(file.path(AQUI, "saida", "tabela_simulacoes.csv"),
               stringsAsFactors = FALSE)
esp53 <- list(
  list("rugosas", 0.000, 0.307, 0.307), list("suaves", 0.000, 0.583, 0.575),
  list("controle", 0.600, 0.612, 0.538))
for (e in esp53) {
  i <- which(ts$cenario == e[[1]])
  conf(paste(e[[1]], "D MQO"), ts$D_ols[i], e[[3]], 1e-3)
  conf(paste(e[[1]], "D SAR"), ts$D_sar[i], e[[4]], 1e-3)
}
r70 <- ts[abs(ts$rho - 0.704) < 1e-6, ]
r70 <- r70[order(-r70$D_ols), ]
conf("realista 4 verd. rho 0,704 D MQO", r70$D_ols[1], 0.584, 1e-3)
conf("realista 4 verd. rho 0,704 D SAR", r70$D_sar[1], 0.438, 1e-3)
conf("realista 8 verd. rho 0,704 D MQO", r70$D_ols[2], 0.464, 1e-3)
conf("realista 8 verd. rho 0,704 D SAR", r70$D_sar[2], 0.264, 1e-3)

cat("\n=== I de Moran dos cenarios ===\n")
mc <- read.csv(file.path(AQUI, "saida", "moran_cenarios.csv"),
               stringsAsFactors = FALSE)
conf("Moran suaves", mc$mi[mc$cenario == "suaves"], 0.125, 1e-3)
conf("Moran rugosas", mc$mi[mc$cenario == "rugosas"], 0.026, 1e-3)
conf("Moran controle", mc$mi[mc$cenario == "controle"], 0.551, 1e-3)

cat("\n=== 2.1 Fator de Ley e Steel ===\n")
K <- 30; m <- 15
conf("fator (K/m + K)/(K/m + 1)", (K / m + K) / (K / m + 1), 10.7, 0.05)

# 2.1: a_g da hyper-g, padrao hyper=UIP do BMS (Zeugner e Feldkircher, 2015).
conf("a_g da hyper-g com n = 184",
     readRDS(file.path(AQUI, "saida", "hyperg_aux.rds"))$f21a, 2.0109, 6e-5)

cat("\n=== 2.2 Erro das tabelas de interpolacao (45) ===\n")
ei <- read.csv(file.path(AQUI, "saida", "erro_interpolacao.csv"), stringsAsFactors = FALSE)
eg <- ei[ei$tabela == "tabela_f21_estendida.rds", ]; ep <- ei[ei$tabela == "tabela_f21.rds", ]
conf("tabela da grade: pontos", eg$pontos, 300001, 0)
conf("tabela da grade: erro no log (e-9)", eg$erro_log_max * 1e9, 9.5, 0.05)
conf("tabela da grade: erro na log-marginal (e-10)", eg$erro_logmarginal_max * 1e10, 5.0, 0.05)
conf("tabela preliminar: erro no log (e-9)", ep$erro_log_max * 1e9, 4.3, 0.05)
conf("tabela preliminar: erro na log-marginal (e-9)", ep$erro_logmarginal_max * 1e9, 2.5, 0.05)

cat("\n=== 4.2 Letalidade: faixa das PIP sob MQO ===\n")
lt <- readRDS(file.path(AQUI, "saida", "grade", "letalidade_queen.rds"))$tabela
lsel <- lt[lt$mqo >= 0.5, ]
conf("acima de 0,50 sob MQO", nrow(lsel), 12, 0)
conf("acima de 0,66", sum(lsel$mqo > 0.66), 2, 0)
conf("entre 0,50 e 0,66", sum(lsel$mqo <= 0.66), 10, 0)
lcai <- lsel[lsel$sar < 0.5, ]
conf("caem abaixo de 0,50 sob SAR", nrow(lcai), 9, 0)
conf("menor queda", min(lcai$mqo - lcai$sar), 0.08, 0.005)
conf("maior queda", max(lcai$mqo - lcai$sar), 0.15, 0.005)

cat("\n=== 4.4 Letalidade: faixa do log BF condicional ===\n")
slt <- read.csv(file.path(AQUI, "saida", "sensibilidade_prior_rho_letalidade.csv"))
conf("minimo com a0 = 1,01", slt$log_bf_cond_min[slt$a0 == 1.01], -0.69, 0.005)
conf("maximo com a0 = 1,01", slt$log_bf_cond_max[slt$a0 == 1.01], 2.55, 0.005)

cat("\n=== 5.1 Calibracao do sinal ===\n")
cs <- readRDS(file.path(AQUI, "saida", "calibracao_sinal.rds"))
pip_escolhido <- function(z) z$pips[z$candidatos == z$beta_escolhido]
conf("calibracao, PIP das verdadeiras, suaves", pip_escolhido(cs$suaves), 0.81, 0.005)
conf("calibracao, PIP das verdadeiras, rugosas", pip_escolhido(cs$rugosas), 0.76, 0.005)
tsim <- read.csv(file.path(AQUI, "saida", "tabela_simulacoes.csv"), stringsAsFactors = FALSE)
conf("experimento, PIP das verdadeiras, suaves", tsim$verd_ols[tsim$cenario == "suaves"], 0.84, 0.005)
conf("experimento, PIP das verdadeiras, rugosas", tsim$verd_ols[tsim$cenario == "rugosas"], 0.64, 0.005)

cat("\n=== 5.6 Conectividade das subamostras (44) ===\n")
cn <- read.csv(file.path(AQUI, "saida", "conectividade_subamostras.csv"))
conf("n = 120: fracao conexa entre as sem ilha (%)", 100 * cn$prop_conexas[cn$n == 120], 19.4, 0.05)
conf("n = 160: fracao conexa entre as sem ilha (%)", 100 * cn$prop_conexas[cn$n == 160], 90.8, 0.05)
conf("sorteios por tamanho", cn$sorteios[1], 2000, 0)

cat("\n=== 2.4 Referencia nula simplificada de D ===\n")
nd <- readRDS(file.path(AQUI, "saida", "nulo_de_D.rds"))
conf("p95 de D sob referencia MQO", unname(nd$MQO$q["95%"]), 0.0155, 1e-4)
conf("p95 de D sob referencia SAR", unname(nd$SAR$q["95%"]), 0.0103, 1e-4)

cat("\n=== 4.4 e 4.5 Sensibilidade a priori de rho (43) ===\n")
sp <- function(d) read.csv(file.path(AQUI, "saida",
                                     sprintf("sensibilidade_prior_rho_%s.csv", d)))
sc20 <- sp("casos"); sl <- sp("letalidade")
so20 <- sp("obitos")
# Casos e obitos no texto: refino consolidado (57 e 58), 2 sementes x 200.000
# modelos cada. O 20.000-modelos original tinha tamanho efetivo baixo do lado
# do MQO nas duas variaveis.
RB <- file.path(AQUI, "saida", "refino_bf")
sc <- read.csv(file.path(RB, "consolidado_principal_casos.csv"))
so <- read.csv(file.path(RB, "consolidado_principal_obitos.csv"))
conf("letalidade: modelos sorteados", sl$n_modelos[1], 20000, 0)
conf("casos log BF agregado, a0 = 1,01", sc$log_bf_agregado_mc[sc$a0 == 1.01], 31.92, 0.006)
conf("casos log BF agregado, a0 = 10", sc$log_bf_agregado_mc[sc$a0 == 10], 27.77, 0.006)
conf("obitos log BF agregado, a0 = 1,01", so$log_bf_agregado_mc[so$a0 == 1.01], 15.85, 0.006)
conf("obitos log BF agregado, a0 = 10", so$log_bf_agregado_mc[so$a0 == 10], 14.20, 0.006)
conf("letalidade P(SAR) agregada, a0 = 1,01", sl$p_sar_agregada_mc[sl$a0 == 1.01], 0.677, 6e-4)
conf("letalidade P(SAR) agregada, a0 = 10", sl$p_sar_agregada_mc[sl$a0 == 10], 0.799, 6e-4)
conf("casos, 20.000 modelos: ESS do MQO", sc20$ess_mqo[1], 284, 0.5)
conf("obitos, 20.000 modelos: menor ESS dos pesos", min(so20$ess_sar, so20$ess_mqo), 82, 0.5)
conf("obitos, 20.000 modelos: maior ESS dos pesos", max(so20$ess_sar, so20$ess_mqo), 136, 0.5)
conf("casos refino: modelos por semente", sc$n_modelos[1] / 2, 200000, 0)
conf("casos refino: menor ESS MQO", min(sc$ess_mqo), 1397, 0.5)
conf("casos refino: menor ESS SAR", min(sc$ess_sar), 24575, 0.5)
conf("casos refino: diferenca entre sementes <= 0,03", as.numeric(max(sc$dif_log_bf_sementes) <= 0.03), 1, 0)
conf("obitos refino: modelos por semente", so$n_modelos[1] / 2, 200000, 0)
conf("obitos refino: menor ESS MQO", min(so$ess_mqo), 501, 0.5)
conf("obitos refino: menor ESS SAR", min(so$ess_sar), 1189, 0.5)
conf("obitos refino: diferenca entre sementes <= 0,03", as.numeric(max(so$dif_log_bf_sementes) <= 0.03), 1, 0)
# "favorecem o SAR em todos os conjuntos examinados": o menor log-BF condicional
# precisa ser positivo em toda a grade de a0.
conf("casos: log BF condicional > 0 em todos os conjuntos",
     as.numeric(min(sc$log_bf_cond_min) > 0), 1, 0)
conf("obitos: log BF condicional > 0 em todos os conjuntos",
     as.numeric(min(so$log_bf_cond_min) > 0), 1, 0)
# Leitura na escala de Kass e Raftery (1995, p. 777), em 2 log B.
conf("casos: menor log BF condicional na grade de a0", min(sc$log_bf_cond_min), 22.95, 0.006)
conf("obitos: menor log BF condicional na grade de a0", min(so$log_bf_cond_min), 10.06, 0.006)
conf("casos: menor 2 log B", 2 * min(sc$log_bf_cond_min), 45.9, 0.06)
conf("obitos: menor 2 log B", 2 * min(so$log_bf_cond_min), 20.1, 0.06)
conf("letalidade: menor 2 log B, a0 = 1,01", 2 * sl$log_bf_cond_min[sl$a0 == 1.01], -1.4, 0.06)
conf("letalidade: maior 2 log B, a0 = 1,01", 2 * sl$log_bf_cond_max[sl$a0 == 1.01], 5.1, 0.06)
# Contagem integrada em letalidade, mantidas as PIP condicionais da cadeia.
el <- readRDS(file.path(AQUI, "saida", "especificacao_letalidade_queen.rds"))
ns <- el$a$n_sar + el$b$n_sar; no <- el$a$n_ols + el$b$n_ols
ps <- (el$a$pip_sar * el$a$n_sar + el$b$pip_sar * el$b$n_sar) / ns
po <- (el$a$pip_ols * el$a$n_ols + el$b$pip_ols * el$b$n_ols) / no
pint <- function(w) w * ps + (1 - w) * po
w0 <- ns / (ns + no)
conf("letalidade frequencia SAR na cadeia", w0, 0.665, 6e-4)
# A contagem integrada depende da fonte das PIP condicionais e da priori de rho
# (auditoria do Antigravity, 10/09/2026). Cada fonte é conferida separadamente.
conf("letalidade, cadeia conjunta: acima de 0,50 sob MQO", sum(po >= 0.5), 14, 0)
conf("letalidade, cadeia conjunta: acima de 0,50 sob SAR", sum(ps >= 0.5), 8, 0)
conf("letalidade, cadeia conjunta: integrada, peso da cadeia", sum(pint(w0) >= 0.5), 8, 0)
lg <- readRDS(file.path(AQUI, "saida", "grade", "letalidade_queen.rds"))$tabela
lg <- lg[match(el$covariaveis, lg$covariavel), ]
conf("letalidade: maior |conjunta - grade| no SAR", max(abs(ps - lg$sar)), 0.071, 6e-4)
conf("letalidade, grade: integrada, peso da cadeia",
     sum(w0 * lg$sar + (1 - w0) * lg$mqo >= 0.5), 6, 0)
z10 <- readRDS(file.path(AQUI, "saida", "sensibilidade_a0", "letalidade_queen_a0_10.rds"))$tabela
z10 <- z10[match(el$covariaveis, z10$covariavel), ]
w10 <- sl$p_sar_agregada_mc[sl$a0 == 10]
conf("letalidade, a0 = 10: acima de 0,50 sob SAR", sum(z10$sar >= 0.5), 7, 0)
conf("letalidade, a0 = 10: razao de D", D(z10$sar) / D(lg$mqo), 0.773, 6e-4)
conf("letalidade, a0 = 10: integrada, peso agregado",
     sum(w10 * z10$sar + (1 - w10) * lg$mqo >= 0.5), 7, 0)

cat("\n=== 1 e 5.2 Suavidade contra correlacao com o desfecho ===\n")
dd <- read.csv(file.path(ESTIMACOES,
                         "resultados/estimacoes/matrizes_estimacao/principal/casos.csv"),
               colClasses = c(codigo_ibge = "character"), stringsAsFactors = FALSE)
de <- readRDS(file.path(AQUI, "saida", "desenho_simulacao.rds"))
cvs <- setdiff(names(dd), c("codigo_ibge", "municipio", "z_log_casos100k"))
ar <- vapply(cvs, function(v) abs(cor(dd$z_log_casos100k, dd[[v]])), numeric(1))
conf("cor(|r| com o desfecho, I de Moran)",
     cor(ar, as.numeric(de$moran[cvs])), 0.69, 6e-3)
conf("r de z_idosos com o desfecho", cor(dd$z_log_casos100k, dd$z_idosos), -0.38, 6e-3)
conf("I de Moran de z_idosos", de$moran[["z_idosos"]], 0.662, 1e-3)
conf("I de Moran de z_tax.agua", de$moran[["z_tax.agua"]], 0.000, 1e-3)
# Figura 1: escolha das covariáveis pelos extremos do I de Moran.
mor <- as.numeric(de$moran[cvs]); names(mor) <- cvs
indic <- vapply(cvs, function(v) length(unique(dd[[v]])) == 2L, logical(1))
conf("Figura 1: covariaveis indicadoras", sum(indic), 3, 0)
conf("Figura 1: idosos e a continua de maior I de Moran",
     as.numeric(names(which.max(mor[!indic])) == "z_idosos"), 1, 0)
conf("Figura 1: agua e a de menor I de Moran",
     as.numeric(names(which.min(mor)) == "z_tax.agua"), 1, 0)
conf("Figura 1: indicadoras com I de Moran acima de idosos",
     sum(mor[indic] > mor[["z_idosos"]]), 2, 0)
conf("Figura 1: I de Moran de casos",
     readRDS(file.path(AQUI, "saida", "grade", "casos_queen.rds"))$mi_y, 0.483, 6e-4)
# Discussão, confundimento espacial: escala das covariáveis selecionadas em casos.
conf("I de Moran de metrop", de$moran[["metrop"]], 0.789, 1e-3)
conf("I de Moran de semiarido", de$moran[["semiarido"]], 0.715, 1e-3)
conf("I de Moran de z_leitos1k", de$moran[["z_leitos1k"]], 0.106, 1e-3)
gq <- readRDS(file.path(AQUI, "saida", "grade", "casos_queen.rds"))$tabela
conf("casos: leitos PIP MQO", gq$mqo[gq$covariavel == "z_leitos1k"], 0.8269, 6e-5)
conf("casos: leitos PIP SAR", gq$sar[gq$covariavel == "z_leitos1k"], 0.3800, 6e-5)
conf("casos: selecionadas com I de Moran > 0,66",
     sum(de$moran[gq$covariavel[gq$mqo >= 0.5]] > 0.66), 3, 0)

cat("\n=== 4.8 Tempo desde o primeiro registro ===\n")
TE <- file.path(AQUI, "saida", "tempo_exposicao")
rs <- read.csv(file.path(TE, "resumo_sensibilidade.csv"), stringsAsFactors = FALSE)
cpt <- read.csv(file.path(TE, "comparacao.csv"), stringsAsFactors = FALSE)
r1 <- function(v, d, col) rs[rs$versao == v & rs$desfecho == d, col]
conf("tempo: maior |MQO - TCC| entre variaveis", max(cpt$dif_max_mqo_tcc[cpt$vizinhos == "queen"]), 0.010, 6e-4)
esp_te <- list(
  list("casos", 0.708, 0.694, 0.301, 0.820, 0.305, 0.681, 0.968, 0.954, 0.833, 7, 2, 1),
  list("obitos", 0.555, 0.553, 0.709, 0.695, 0.504, 0.516, 0.342, 0.271, 0.290, 1, 2, 1),
  list("letalidade", 0.254, 0.245, 0.706, 0.807, 0.709, 0.752, 0.810, 0.734, 0.629, 31, 29, 10))
for (e in esp_te) {
  d <- e[[1]]
  conf(paste(d, "rho principal"), r1("principal", d, "rho_media"), e[[2]], 6e-4)
  conf(paste(d, "rho com dias"), r1("com dias", d, "rho_media"), e[[3]], 6e-4)
  conf(paste(d, "razao D SAR principal"), r1("principal", d, "razao_sar"), e[[4]], 6e-4)
  conf(paste(d, "razao D SAR com dias"), r1("com dias", d, "razao_sar"), e[[5]], 6e-4)
  conf(paste(d, "razao D SEM principal"), r1("principal", d, "razao_sem"), e[[6]], 6e-4)
  conf(paste(d, "razao D SEM com dias"), r1("com dias", d, "razao_sem"), e[[7]], 6e-4)
  conf(paste(d, "PIP dias MQO"), r1("com dias", d, "pip_dias_mqo"), e[[8]], 6e-4)
  conf(paste(d, "PIP dias SAR"), r1("com dias", d, "pip_dias_sar"), e[[9]], 6e-4)
  conf(paste(d, "PIP dias SEM"), r1("com dias", d, "pip_dias_sem"), e[[10]], 6e-4)
  conf(paste(d, ">= 0,50 MQO com dias"), r1("com dias", d, "acima_mqo"), e[[11]], 0)
  conf(paste(d, ">= 0,50 SAR com dias"), r1("com dias", d, "acima_sar"), e[[12]], 0)
  conf(paste(d, ">= 0,50 SEM com dias"), r1("com dias", d, "acima_sem"), e[[13]], 0)
}
conf("casos com dias: rho, limite inferior", r1("com dias", "casos", "rho_inf"), 0.569, 6e-4)
conf("casos com dias: rho, limite superior", r1("com dias", "casos", "rho_sup"), 0.809, 6e-4)
conf("casos: log BF agregado principal", r1("principal", "casos", "log_bf_agregado"), 31.92, 6e-3)
conf("casos: log BF agregado com dias", r1("com dias", "casos", "log_bf_agregado"), 31.89, 6e-3)
conf("obitos: log BF agregado principal", r1("principal", "obitos", "log_bf_agregado"), 15.85, 6e-3)
conf("obitos: log BF agregado com dias", r1("com dias", "obitos", "log_bf_agregado"), 15.81, 6e-3)
conf("letalidade P(SAR) a0 1,01 principal", r1("principal", "letalidade", "p_sar_a0_1.01"), 0.677, 6e-4)
conf("letalidade P(SAR) a0 1,01 com dias", r1("com dias", "letalidade", "p_sar_a0_1.01"), 0.627, 6e-4)
conf("letalidade P(SAR) a0 10 principal", r1("principal", "letalidade", "p_sar_a0_10"), 0.799, 6e-4)
conf("letalidade P(SAR) a0 10 com dias", r1("com dias", "letalidade", "p_sar_a0_10"), 0.770, 6e-4)
conf("casos/obitos: freq SAR com dias = 1", min(r1("com dias", "casos", "freq_sar"), r1("com dias", "obitos", "freq_sar")), 1, 0)
tc <- readRDS(file.path(TE, "casos_queen.rds")); sc <- readRDS(file.path(TE, "sem_casos_queen.rds"))
tt <- tc$tabela; tt$sem <- sc$pip[match(tt$covariavel, sc$covariaveis)]
seis <- tt[tt$covariavel %in% c("z_idosos", "z_leitos1k", "metrop", "semiarido", "z_ideb.fund", "z_int.circ"), ]
conf("casos com dias: seis, menor PIP MQO", min(seis$mqo), 0.574, 6e-4)
conf("casos com dias: seis, maior PIP MQO", max(seis$mqo), 0.889, 6e-4)
conf("casos com dias: seis, menor PIP SAR", min(seis$sar), 0.219, 6e-4)
conf("casos com dias: seis, maior PIP SAR", max(seis$sar), 0.381, 6e-4)
conf("casos com dias: seis, menor PIP SEM", min(seis$sem), 0.225, 6e-4)
conf("casos com dias: seis, maior PIP SEM", max(seis$sem), 0.405, 6e-4)
conf("casos com dias: imunizacao sob SAR", tt$sar[tt$covariavel == "z_imuni"], 0.516, 6e-4)
conf("casos com dias: discordancia SAR", tc$res$SAR$conv, 0.024, 6e-4)
to <- readRDS(file.path(TE, "obitos_queen.rds")); so2 <- readRDS(file.path(TE, "sem_obitos_queen.rds"))
ot <- to$tabela; ot$sem <- so2$pip[match(ot$covariavel, so2$covariaveis)]
conf("obitos com dias: idosos MQO", ot$mqo[ot$covariavel == "z_idosos"], 0.9669, 6e-5)
conf("obitos com dias: idosos SAR", ot$sar[ot$covariavel == "z_idosos"], 0.3995, 6e-5)
conf("obitos com dias: idosos SEM", ot$sem[ot$covariavel == "z_idosos"], 0.1577, 6e-5)
conf("obitos com dias: pop rural SAR", ot$sar[ot$covariavel == "z_pop.rural"], 0.6818, 6e-5)
conf("obitos com dias: pop rural SEM", ot$sem[ot$covariavel == "z_pop.rural"], 0.5072, 6e-5)
conf("obitos com dias: mulheres eleitas SAR", ot$sar[ot$covariavel == "z_eleitas.fem"], 0.5229, 6e-5)
conf("letalidade com dias: k MQO", r1("com dias", "letalidade", "k_mqo"), 19.62, 6e-3)
conf("letalidade com dias: referencia MQO", r1("com dias", "letalidade", "k_mqo") / 31, 0.633, 6e-4)
conf("letalidade com dias: referencia SAR", r1("com dias", "letalidade", "ref_sar"), 0.579, 6e-4)
conf("letalidade com dias: referencia SEM", r1("com dias", "letalidade", "ref_sem"), 0.508, 6e-4)

# Queen contra Rook na especificacao "com dias" (48, 53 e 55), do
# 63_comparar_vizinhanca_tempo.R
ct <- read.csv(file.path(TE, "comparacao_vizinhanca_tempo.csv"), stringsAsFactors = FALSE)
w2 <- function(d, viz, col) ct[ct$desfecho == d & ct$vizinhos == viz, col]
conf("com dias, Rook: nenhuma visita ao MQO em casos e obitos",
     sum(w2("casos", "rook", "acima_mqo") >= 0, w2("obitos", "rook", "acima_mqo") >= 0), 2, 0)
conf("com dias casos Rook: D SEM", w2("casos", "rook", "D_sem"), 0.2892, 1e-4)
conf("com dias casos Rook: acima 0,50 SEM", w2("casos", "rook", "acima_sem"), 2, 0)
conf("com dias casos Rook: PIP dias SEM", w2("casos", "rook", "pip_dias_sem"), 0.906, 6e-4)
conf("com dias letalidade Rook: acima 0,50 SAR", w2("letalidade", "rook", "acima_sar"), 31, 0)
conf("com dias letalidade Rook: acima 0,50 SEM", w2("letalidade", "rook", "acima_sem"), 13, 0)
conf("com dias obitos Rook: acima 0,50 SEM", w2("obitos", "rook", "acima_sem"), 1, 0)

conf("com dias casos Rook: log BF agregado", w2("casos", "rook", "log_bf_agregado"), 30.48, 6e-3)
conf("com dias obitos Rook: log BF agregado", w2("obitos", "rook", "log_bf_agregado"), 16.04, 6e-3)
conf("com dias letalidade Rook: log BF agregado", w2("letalidade", "rook", "log_bf_agregado"), 0.62, 6e-3)
conf("com dias: maior diferenca de log BF entre matrizes",
     max(abs(ct$log_bf_agregado[ct$vizinhos == "queen"] - ct$log_bf_agregado[ct$vizinhos == "rook"])),
     1.454, 6e-3)

cat(sprintf("\n=== %d testes, %d divergencias ===\n", testes, falhas))
if (falhas > 0L) cat("Corrigir o texto ou a transcricao antes de submeter.\n")
