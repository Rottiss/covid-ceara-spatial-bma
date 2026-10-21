# Export consolidated CSV inputs for Python plotting scripts.

source(if (file.exists("../caminhos.R")) "../caminhos.R" else "caminhos.R")

AQUI <- ESPACIAL
SAIDA <- file.path(AQUI, "saida")
EXP <- file.path(SAIDA, "export_python")
dir.create(EXP, showWarnings = FALSE, recursive = TRUE)

source(file.path(AQUI, "rotulos.R"))
base <- readRDS(file.path(SAIDA, "base_sar.rds"))
des  <- readRDS(file.path(SAIDA, "desenho_simulacao.rds"))

D_calc <- function(p) {
  m <- mean(p); K <- length(p)
  sd(p) * sqrt((K - 1) / K) / sqrt(m * (1 - m))
}

# -------------------------------------------------------------
# 1. Tabela Completa de PIPs e Transicoes
# -------------------------------------------------------------
desfechos <- c("casos", "obitos", "letalidade")
covs <- base$covariaveis

# Definir grupos tematicos
grupos <- list(
  Demografia = c("z_pop", "z_area", "z_densidade", "z_idosos", "z_pop.rural"),
  `Vulnerabilidade e Renda` = c("z_pibpc", "z_desp.educ", "z_ideb.fund", "z_ideb3", 
                                "z_bolsaf", "z_ex.pobr", "z_emprego", "z_ivs.infra"),
  `Saúde e Saneamento` = c("z_tax.agua", "z_energia100k", "z_tax.hom", "z_sus1k", 
                           "z_leitos1k", "z_prof1k", "z_desp.saude", "z_imuni", 
                           "z_int.circ", "z_int.resp", "z_int.diab", "z_int.asma"),
  `Território e Política` = c("metrop", "semiarido", "aliadogov", "z_votos", "z_eleitas.fem")
)

grupo_map <- character(length(covs))
names(grupo_map) <- covs
for (g_nome in names(grupos)) {
  for (v in grupos[[g_nome]]) {
    if (v %in% covs) grupo_map[v] <- g_nome
  }
}

lista_pips <- list()

for (d in desfechos) {
  g_queen <- readRDS(file.path(SAIDA, "grade", paste0(d, "_queen.rds")))
  g_rook  <- readRDS(file.path(SAIDA, "grade", paste0(d, "_rook.rds")))
  f_sem   <- file.path(SAIDA, "grade", paste0("sem_", d, "_queen.rds"))
  
  t_q <- g_queen$tabela
  t_r <- g_rook$tabela
  
  df <- data.frame(
    covariavel = t_q$covariavel,
    rotulo = rotular(t_q$covariavel, "pt"),
    grupo = grupo_map[t_q$covariavel],
    desfecho = d,
    moran_cov = as.numeric(des$moran[t_q$covariavel]),
    pip_mqo = t_q$mqo,
    pip_sar_queen = t_q$sar,
    pip_sar_rook = t_r$sar[match(t_q$covariavel, t_r$covariavel)],
    stringsAsFactors = FALSE
  )
  
  if (file.exists(f_sem)) {
    g_sem <- readRDS(f_sem)
    df$pip_sem_queen <- g_sem$pip[match(t_q$covariavel, g_sem$covariaveis)]
  } else {
    df$pip_sem_queen <- NA_real_
  }
  
  # Categorizacao do destino MQO -> SAR Queen
  df$status_sar <- ifelse(df$pip_mqo >= 0.5 & df$pip_sar_queen < 0.5, "Sai da lista (≥ 0,50 → < 0,50)",
                   ifelse(df$pip_mqo < 0.5 & df$pip_sar_queen >= 0.5, "Entra na lista (< 0,50 → ≥ 0,50)",
                   ifelse(df$pip_mqo >= 0.5 & df$pip_sar_queen >= 0.5, "Permanece na lista (≥ 0,50 → ≥ 0,50)",
                          "Abaixo do limiar (< 0,50)")))
  
  # Linhas de base
  k_mqo <- sum(df$pip_mqo)
  k_sar <- sum(df$pip_sar_queen)
  df$linha_base_mqo <- k_mqo / 30
  df$linha_base_sar <- k_sar / 30
  
  lista_pips[[d]] <- df
}

pips_total <- do.call(rbind, lista_pips)
write.csv(pips_total, file.path(EXP, "pips_completas.csv"), row.names = FALSE)
cat("Exportado: pips_completas.csv\n")

# -------------------------------------------------------------
# 2. Densidades Posteriores de Rho
# -------------------------------------------------------------
lista_rho <- list()
moran_ys <- c(casos = 0.483366, obitos = 0.441355, letalidade = 0.138410)

for (d in desfechos) {
  esp <- readRDS(file.path(SAIDA, paste0("especificacao_", d, "_queen.rds")))
  rho_grid <- esp$rho
  h <- rho_grid[2] - rho_grid[1]
  w_media <- (esp$a$w_rho + esp$b$w_rho) / 2
  dens <- w_media / h
  
  df_rho <- data.frame(
    desfecho = d,
    rho = rho_grid,
    densidade = dens,
    moran_y = moran_ys[d],
    stringsAsFactors = FALSE
  )
  lista_rho[[d]] <- df_rho
}
rho_total <- do.call(rbind, lista_rho)
write.csv(rho_total, file.path(EXP, "densidades_rho.csv"), row.names = FALSE)
cat("Exportado: densidades_rho.csv\n")

# -------------------------------------------------------------
# 3. Replicações Detalhadas de Simulação
# -------------------------------------------------------------
sim_det <- list()

# 1. Experimentos da primeira rodada
exp_sim <- readRDS(file.path(SAIDA, "experimento_simulacao.rds"))
for (nm in names(exp_sim)) {
  obj <- exp_sim[[nm]]
  reps <- nrow(obj$ols)
  d_ols <- apply(obj$ols, 1, D_calc)
  d_sar <- apply(obj$sar, 1, D_calc)
  verd_ols <- rowMeans(obj$ols[, obj$idx, drop = FALSE])
  verd_sar <- rowMeans(obj$sar[, obj$idx, drop = FALSE])
  
  df <- data.frame(
    cenario = nm,
    categoria = "Sintético (processo gerador conhecido)",
    rho_dgp = obj$rho,
    n_verdadeiras = length(obj$idx),
    replicacao = 1:reps,
    D_mqo = d_ols,
    D_sar = d_sar,
    razao_D = d_sar / d_ols,
    verd_mqo = verd_ols,
    verd_sar = verd_sar,
    stringsAsFactors = FALSE
  )
  sim_det[[length(sim_det) + 1]] <- df
}

# 2. Cenários realistas
arqs_real <- c("cenario_realista.rds", "cenario_realista_oito.rds",
               "cenario_rho070_quatro.rds", "cenario_rho070_oito.rds")
for (ar in arqs_real) {
  f <- file.path(SAIDA, ar)
  if (!file.exists(f)) next
  cr <- readRDS(f)
  for (nm in names(cr$resultados)) {
    obj <- cr$resultados[[nm]]
    reps <- nrow(obj$ols)
    d_ols <- apply(obj$ols, 1, D_calc)
    d_sar <- apply(obj$sar, 1, D_calc)
    verd_ols <- rowMeans(obj$ols[, obj$idx, drop = FALSE])
    verd_sar <- rowMeans(obj$sar[, obj$idx, drop = FALSE])
    
    rotulo <- sprintf("Realista %s verd. (rho = %.3f)", nm, cr$rho)
    df <- data.frame(
      cenario = rotulo,
      categoria = "Realista (calibrado)",
      rho_dgp = cr$rho,
      n_verdadeiras = length(obj$idx),
      replicacao = 1:reps,
      D_mqo = d_ols,
      D_sar = d_sar,
      razao_D = d_sar / d_ols,
      verd_mqo = verd_ols,
      verd_sar = verd_sar,
      stringsAsFactors = FALSE
    )
    sim_det[[length(sim_det) + 1]] <- df
  }
}

sim_total <- do.call(rbind, sim_det)
write.csv(sim_total, file.path(EXP, "simulacoes_detalhadas.csv"), row.names = FALSE)
cat("Exportado: simulacoes_detalhadas.csv\n")

# -------------------------------------------------------------
# 4. Dados do Mecanismo (Colinearidade x Moran x Queda)
# -------------------------------------------------------------
mec <- read.csv(file.path(SAIDA, "mecanismo_dados_reais.csv"), stringsAsFactors = FALSE)
mec$rotulo <- rotular(mec$covariavel, "pt")
mec$grupo <- grupo_map[mec$covariavel]
write.csv(mec, file.path(EXP, "mecanismo_completo.csv"), row.names = FALSE)
cat("Exportado: mecanismo_completo.csv\n")

cat("Todos os dados foram preparados com sucesso em:", EXP, "\n")
