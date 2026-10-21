# Export CSV/GeoJSON inputs for figure scripts 51 and 52.
# Calculations follow scripts 20, 31, 39, 41, and gerar_tabelas.R.
#
#   pips_queen.csv          Figure 3: PIPs under OLS (MQO), SAR, and SEM
#   malha_ce.geojson        Figure 1: boundaries with standardized variables
#   suavidade.csv           Figure 2: Moran's I, case correlation, OLS PIP
#   simulacao_reps.csv      Figure 4: D ratio by replication
#   simulacao_cenarios.csv  Figure 4: D ratios by scenario and observed data
#   rho_densidade.csv       Figure 5: posterior density of rho
#   rho_resumo.csv          Figure 5: mean, interval, Moran's I, D ratio
#   priori_rho.csv          Figure 6: Bayes factor and P(SAR | y) against a0
# Exported identifiers and Portuguese source labels are kept for compatibility.
# The Python display layer selects English or Portuguese labels.

source(if (file.exists("../caminhos.R")) "../caminhos.R" else "caminhos.R")

suppressPackageStartupMessages(library(sf))

AQUI <- ESPACIAL
RAIZ <- ESTIMACOES
DEST <- Sys.getenv("SAIDA_DADOS_FIGURAS", file.path(AQUI, "saida", "figuras_v3", "dados"))
dir.create(DEST, showWarnings = FALSE, recursive = TRUE)
source(file.path(AQUI, "rotulos.R"))

D <- function(p) { m <- mean(p); K <- length(p)
  sd(p) * sqrt((K - 1) / K) / sqrt(m * (1 - m)) }
gravar_csv <- function(tab, nome) {
  stopifnot(!anyNA(tab))
  old_opts <- options(digits = 15)
  con <- file(file.path(DEST, nome), open = "w", encoding = "UTF-8")
  write.csv(tab, con, row.names = FALSE); close(con)
  options(old_opts)
  cat(sprintf("  %-24s %d rows\n", nome, nrow(tab)))
}
DESF <- c("casos", "obitos", "letalidade")
grade <- lapply(setNames(DESF, DESF), function(d)
  readRDS(file.path(AQUI, "saida", "grade", sprintf("%s_queen.rds", d))))

# ------------------------------------------------------------ Figure 3
linhas <- list()
for (d in DESF) {
  g <- grade[[d]]
  s <- readRDS(file.path(AQUI, "saida", "grade", sprintf("sem_%s_queen.rds", d)))
  stopifnot(length(s$pip) == length(s$covariaveis), setequal(s$covariaveis, g$tabela$covariavel))
  t <- g$tabela
  t$sem <- s$pip[match(t$covariavel, s$covariaveis)]
  linhas[[d]] <- data.frame(
    desfecho = d, covariavel = t$covariavel, rotulo = inicial_maiuscula(rotular(t$covariavel, "pt")),
    pip_mqo = t$mqo, pip_sar = t$sar, pip_sem = t$sem,
    ref_mqo = g$res$MQO$tam / g$K, ref_sar = g$res$SAR$tam / g$K, ref_sem = s$tam / g$K,
    stringsAsFactors = FALSE)
}
pips <- do.call(rbind, linhas)
gravar_csv(pips, "pips_queen.csv")

# ------------------------------------------------------- Figures 1 and 2
dados <- read.csv(file.path(RAIZ, "resultados/estimacoes/matrizes_estimacao/principal/casos.csv"),
                  colClasses = c(codigo_ibge = "character"), stringsAsFactors = FALSE)
des <- readRDS(file.path(AQUI, "saida", "desenho_simulacao.rds"))

malha <- readRDS(file.path(AQUI, "saida", "malha_ce.rds"))
malha$codigo_ibge <- as.character(malha$code_muni)
malha <- merge(malha, dados[, c("codigo_ibge", "z_log_casos100k", "z_idosos", "z_tax.agua")],
               by = "codigo_ibge")
stopifnot(nrow(malha) == 184, !anyNA(st_drop_geometry(malha)[, c("z_log_casos100k", "z_idosos", "z_tax.agua")]))
names(malha)[names(malha) == "z_tax.agua"] <- "z_tax_agua"
arq_geo <- file.path(DEST, "malha_ce.geojson")
if (file.exists(arq_geo)) file.remove(arq_geo)
st_write(malha[, c("codigo_ibge", "z_log_casos100k", "z_idosos", "z_tax_agua")], arq_geo,
         driver = "GeoJSON", quiet = TRUE)
cat(sprintf("  %-24s %d municipalities | Moran's I: cases %.3f, aged 60+ %.3f, water %.3f\n",
            "malha_ce.geojson", nrow(malha), des$mi_y_real, des$moran[["z_idosos"]], des$moran[["z_tax.agua"]]))

cv <- setdiff(names(dados), c("codigo_ibge", "municipio", "z_log_casos100k"))
tc <- grade$casos$tabela
suav <- data.frame(
  covariavel = cv, rotulo = inicial_maiuscula(rotular(cv, "pt")),
  moran = as.numeric(des$moran[cv]),
  r_casos = vapply(cv, function(v) cor(dados$z_log_casos100k, dados[[v]]), numeric(1)),
  pip_mqo = tc$mqo[match(cv, tc$covariavel)], stringsAsFactors = FALSE)
gravar_csv(suav, "suavidade.csv")
cat(sprintf("    correlation between Moran's I and |r|: %.3f\n", cor(suav$moran, abs(suav$r_casos))))

# ------------------------------------------------------------ Figure 4
nome_cen <- c(rugosas = "Baixa autocorrelação", suaves = "Alta autocorrelação", controle = "Controle")
cen <- list()
sim <- readRDS(file.path(AQUI, "saida", "experimento_simulacao.rds"))
for (nm in names(sim)) cen[[length(cen) + 1]] <- list(nome = nome_cen[[nm]], r = sim[[nm]], rho = sim[[nm]]$rho)
for (a in c("cenario_realista.rds", "cenario_realista_oito.rds",
            "cenario_rho070_quatro.rds", "cenario_rho070_oito.rds")) {
  cr <- readRDS(file.path(AQUI, "saida", a))
  for (nm in names(cr$resultados))
    cen[[length(cen) + 1]] <- list(nome = sprintf("Realista, %s verdadeiras", nm),
                                  r = cr$resultados[[nm]], rho = cr$rho)
}
resumo_cen <- do.call(rbind, lapply(seq_along(cen), function(i) {
  z <- cen[[i]]
  data.frame(id = i, nome = z$nome, rho = z$rho,
             razao = mean(apply(z$r$sar, 1, D)) / mean(apply(z$r$ols, 1, D)),
             stringsAsFactors = FALSE)
}))
resumo_cen <- resumo_cen[order(-resumo_cen$razao), ]
# Match gerar_tabelas.R: repeated scenarios receive (1)/(2), ordered by ratio.
chave <- paste(resumo_cen$nome, round(resumo_cen$rho, 2))
resumo_cen$sufixo <- ave(seq_along(chave), chave, FUN = function(i) if (length(i) > 1) seq_along(i) else 0L)
resumo_cen$rotulo <- sprintf("%s%s, ρ = %s", resumo_cen$nome,
                             ifelse(resumo_cen$sufixo > 0, sprintf(" (%d)", resumo_cen$sufixo), ""),
                             gsub(".", ",", sprintf("%.2f", resumo_cen$rho), fixed = TRUE))
tsim <- read.csv(file.path(AQUI, "saida", "tabela_simulacoes.csv"))
stopifnot(isTRUE(all.equal(sort(resumo_cen$razao), sort(tsim$razao), tolerance = 1e-10)))
reps <- do.call(rbind, lapply(seq_len(nrow(resumo_cen)), function(j) {
  z <- cen[[resumo_cen$id[j]]]
  data.frame(rotulo = resumo_cen$rotulo[j], replicacao = seq_len(nrow(z$r$ols)),
             razao = apply(z$r$sar, 1, D) / apply(z$r$ols, 1, D), stringsAsFactors = FALSE)
}))
reais <- data.frame(rotulo = c("Casos", "Óbitos", "Letalidade"), tipo = "observado",
                    razao = vapply(DESF, function(d) D(grade[[d]]$res$SAR$pip) / D(grade[[d]]$res$MQO$pip), numeric(1)),
                    stringsAsFactors = FALSE)
cenarios <- rbind(data.frame(rotulo = resumo_cen$rotulo, tipo = "simulado", razao = resumo_cen$razao,
                             stringsAsFactors = FALSE), reais)
gravar_csv(reps, "simulacao_reps.csv")
gravar_csv(cenarios, "simulacao_cenarios.csv")

# ------------------------------------------------------------ Figure 5
pesos_simpson <- function(ng) c(1, rep(c(4, 2), length.out = ng - 2), 1)
dens <- list(); res_rho <- list()
for (d in DESF) {
  e <- readRDS(file.path(AQUI, "saida", sprintf("especificacao_%s_queen.rds", d)))
  w <- (e$a$w_rho + e$b$w_rho) / 2
  ps <- pesos_simpson(length(e$rho)); h <- e$rho[2] - e$rho[1]
  ex <- w / ps
  dens[[d]] <- data.frame(desfecho = d, rho = e$rho, densidade = ex / (h / 3 * sum(ps * ex)))
  wn <- w / sum(w); ac <- cumsum(wn)
  res_rho[[d]] <- data.frame(
    desfecho = d, media = sum(wn * e$rho),
    inf = e$rho[which.max(ac >= 0.025)], sup = e$rho[which.max(ac >= 0.975)],
    moran_y = grade[[d]]$mi_y,
    razao_D = D(grade[[d]]$res$SAR$pip) / D(grade[[d]]$res$MQO$pip))
}
gravar_csv(do.call(rbind, dens), "rho_densidade.csv")
rr <- do.call(rbind, res_rho)
gravar_csv(rr, "rho_resumo.csv")
print(rr, digits = 3, row.names = FALSE)

# ------------------------------------------------------------ Figure 6
# Cases/deaths use the consolidated refinement (57/58: 2 x 200,000 models).
# The 20,000-model run in script 43 had low effective size for OLS weights.
REFINADOS_FIG6 <- c("casos", "obitos")
pri <- do.call(rbind, lapply(DESF, function(d) {
  arq <- if (d %in% REFINADOS_FIG6) file.path(AQUI, "saida", "refino_bf", sprintf("consolidado_principal_%s.csv", d))
         else file.path(AQUI, "saida", sprintf("sensibilidade_prior_rho_%s.csv", d))
  read.csv(arq)[, c("desfecho", "a0", "log_bf_cond_min", "log_bf_cond_mediano", "log_bf_cond_max",
                    "log_bf_agregado_mc", "p_sar_agregada_mc")]
}))
gravar_csv(pri[, c("desfecho", "a0", "log_bf_cond_min", "log_bf_cond_mediano", "log_bf_cond_max",
                   "log_bf_agregado_mc", "p_sar_agregada_mc")], "priori_rho.csv")
