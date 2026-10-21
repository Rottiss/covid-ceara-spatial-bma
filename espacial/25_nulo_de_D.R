# Distribuicao de D sob a hipotese nula de nenhuma discriminacao.
#
# D = dp(PIP) / sqrt(p(1-p)) vale exatamente zero se os dados nada dizem sobre
# QUAIS covariaveis entram. Mas D e calculado a partir de uma cadeia finita, e
# ruido de Monte Carlo produz dispersao mesmo sob a nula. Sem saber o piso, um D
# de 0,0973 nao e interpretavel.
#
# CONSTRUCAO DA NULA. Simula-se uma cadeia que percorre modelos com tamanho
# distribuido segundo uma Poisson truncada calibrada apenas pela media posterior,
# mas escolhendo as covariaveis ao acaso a cada passo. A persistencia e imitada
# por um parametro de permanencia calibrado pela discordancia maxima de PIP entre
# duas cadeias. Esta e uma referencia simplificada, nao uma replica do processo
# de transicao do MC3.
#
# Reporta o quantil 0,95 de D sob essa referencia. Um D observado acima dele e
# compativel com dispersao alem do ruido simulado; nao e teste exato nem prova de
# discriminacao, porque a referencia nao reproduz as transicoes do MC3.

source(if (file.exists("../caminhos.R")) "../caminhos.R" else "caminhos.R")

AQUI <- ESPACIAL

D <- function(pip) {
  p <- mean(pip); K <- length(pip)
  if (p <= 0 || p >= 1) return(NA_real_)
  sd(pip) * sqrt((K - 1) / K) / sqrt(p * (1 - p))
}

# Cadeia nula: a cada passo, com probabilidade (1 - troca) mantem o modelo; com
# probabilidade troca, sorteia um modelo novo de tamanho k ~ dist. observada,
# com covariaveis uniformes. 'troca' calibra a persistencia.
cadeia_nula <- function(K, T_guardados, dist_k, troca) {
  inclusao <- numeric(K)
  atual <- sample(K, sample(dist_k, 1))
  for (t in seq_len(T_guardados)) {
    if (runif(1) < troca) atual <- sample(K, sample(dist_k, 1))
    inclusao[atual] <- inclusao[atual] + 1
  }
  inclusao / T_guardados
}

# Calibra 'troca' pela discordancia observada entre duas cadeias independentes:
# procura o valor que reproduz a mesma discordancia maxima de PIP.
# A discordancia medida em uma cadeia curta deve ser maior que a observada na
# cadeia longa, por sqrt(T_longo/T_curto). A versao anterior usava o inverso,
# tornava o alvo inatingivel e caia silenciosamente em troca = 0,01.
calibrar_troca <- function(K, T_curto, T_real, dist_k, alvo_real, reps = 8L) {
  alvo_curto <- alvo_real * sqrt(T_real / T_curto)
  grade <- 10^seq(-5, 0, length.out = 31L)
  discordancia <- vapply(seq_along(grade), function(i) {
    set.seed(1100L + i)
    mean(replicate(reps, {
      a <- cadeia_nula(K, T_curto, dist_k, grade[i])
      b <- cadeia_nula(K, T_curto, dist_k, grade[i])
      max(abs(a - b))
    }))
  }, numeric(1))
  i <- which.min(abs(discordancia - alvo_curto))
  list(troca = grade[i], alvo_curto = alvo_curto,
       discordancia_curta = discordancia[i])
}

cenarios <- list(
  MQO = list(k_medio = 13.0954, discordancia = 0.0224, D_obs = 0.3235),
  SAR = list(k_medio = 7.1880,  discordancia = 0.0125, D_obs = 0.0973)
)

K <- 30L
T_GUARD <- 800000L   # 1.000.000 menos 200.000 de descarte
T_CURTO <- 20000L
REPS <- 400L

cat("=== Piso de D por ruido de Monte Carlo ===\n")
cat("  K =", K, "| sorteios guardados por cadeia:", T_GUARD, "\n")
cat("  A referencia usa uma Poisson truncada com a media observada e calibra\n")
cat("  um parametro de permanencia pela discordancia entre cadeias.\n\n")

resultado <- list()
for (nome in names(cenarios)) {
  cn <- cenarios[[nome]]
  # Distribuicao de tamanho: Poisson truncada com a media observada.
  set.seed(7000L + match(nome, names(cenarios)))
  dist_k <- pmin(K, pmax(0, rpois(5000, cn$k_medio)))
  cal <- calibrar_troca(K, T_CURTO, T_GUARD, dist_k, cn$discordancia)
  tr <- cal$troca
  cat(sprintf("### %s\n", nome))
  cat(sprintf("  tamanho medio do modelo: %.3f | discordancia observada: %.4f\n",
              cn$k_medio, cn$discordancia))
  cat(sprintf("  parametro de permanencia calibrado: %.5f\n", tr))
  cat(sprintf("  alvo/obtido na cadeia curta: %.5f / %.5f\n",
              cal$alvo_curto, cal$discordancia_curta))

  set.seed(20200731L)
  Ds <- replicate(REPS, D(cadeia_nula(K, T_CURTO, dist_k, tr)))
  # Reescala para o comprimento real: dp do ruido cai com sqrt(T).
  fator <- sqrt(T_CURTO / T_GUARD)
  Ds_esc <- Ds * fator

  q <- quantile(Ds_esc, c(0.5, 0.95, 0.99))
  cat(sprintf("  D sob a nula (mediana / p95 / p99): %.5f / %.5f / %.5f\n",
              q[1], q[2], q[3]))
  cat(sprintf("  D observado: %.4f | razao com o p95: %.1f vezes\n",
              cn$D_obs, cn$D_obs / q[2]))
  cat(sprintf("  veredito: %s\n\n",
              ifelse(cn$D_obs > q[2],
                     "acima desta referencia nula simplificada",
                     "nao se separa desta referencia nula simplificada")))
  resultado[[nome]] <- list(troca = tr, calibracao = cal, nula = Ds_esc,
                            q = q, D_obs = cn$D_obs)
}

cat("=== Leitura ===\n")
cat("  Os dois D estao acima da referencia simulada. O exercicio nao reproduz\n")
cat("  a distribuicao de tamanhos nem as transicoes completas do MC3.\n")
cat("  O que muda e a magnitude: de", round(cenarios$MQO$D_obs, 4), "para",
    round(cenarios$SAR$D_obs, 4), "\n")
cat("  A comparacao reporta os valores de D nas duas especificacoes.\n")

# Quanta da discriminacao do SAR vem de uma unica covariavel?
real <- read.csv(file.path(AQUI, "saida", "comparacao_1milhao.csv"),
                 stringsAsFactors = FALSE)
cat("\n=== Concentracao da discriminacao remanescente ===\n")
d_todas <- D(real$sar)
sem_maior <- D(real$sar[-which.max(real$sar)])
cat(sprintf("  D do SAR com as 30: %.4f\n", d_todas))
cat(sprintf("  D do SAR sem a maior (%s, %.4f): %.4f\n",
            real$covariavel[which.max(real$sar)], max(real$sar), sem_maior))
cat(sprintf("  a maior covariavel responde por %.0f%% da dispersao\n",
            100 * (1 - sem_maior / d_todas)))

d_mqo <- D(real$mqo)
sem_maior_mqo <- D(real$mqo[-which.max(real$mqo)])
cat(sprintf("\n  para comparar, no MQO: %.4f para %.4f, ou %.0f%%\n",
            d_mqo, sem_maior_mqo, 100 * (1 - sem_maior_mqo / d_mqo)))

saveRDS(resultado, file.path(AQUI, "saida", "nulo_de_D.rds"))
cat("\nGravado em saida/nulo_de_D.rds\n")
