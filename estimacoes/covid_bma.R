# Data preparation and non-spatial BMA estimation.
# Municipal data for Ceara; cumulative COVID-19 records through 2020-07-31.

pacotes_necessarios <- c("readxl", "dplyr", "BMS")
pacotes_ausentes <- pacotes_necessarios[
  !vapply(pacotes_necessarios, requireNamespace, logical(1), quietly = TRUE)
]

if (length(pacotes_ausentes) > 0L) {
  stop(
    "Install missing packages: ",
    paste(pacotes_ausentes, collapse = ", "),
    call. = FALSE
  )
}

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(BMS)
})

# 1. Read the dataset -----------------------------------------------------

arquivo_dados <- "dadosv6.xlsx"

if (!file.exists(arquivo_dados)) {
  stop(
    "dadosv6.xlsx was not found in the project directory.",
    call. = FALSE
  )
}

dados <- read_excel(
  path = arquivo_dados,
  sheet = "raw",
  na = c("", "NA")
)

# The original population column header includes whitespace.
names(dados) <- trimws(names(dados))

colunas_obrigatorias <- c(
  "Municipios", "codigo_ibge", "covtotal", "cov100k", "obito",
  "obito100k", "letal", "pop", "recursofed", "auxem", "area",
  "densidade", "pibpc", "idhm", "idm", "desp.saude", "desp.educ",
  "ideb5", "ideb9", "ideb3", "tax.agua", "con.energia", "tax.hom",
  "eleitas.fem", "bolsaf", "sus1k", "leitos1k", "prof1k", "metrop",
  "semiarido", "pop.rural", "aliadogov", "votos", "ivs", "ivs.infra",
  "ivs.capital", "ivs.renda", "emprego", "ex.pobr", "tmi", "imuni",
  "idosos", "int.circ", "int.resp", "int.diab", "int.asma", "ideb.fund",
  "dias.primeiro.registro"
)

colunas_ausentes <- setdiff(colunas_obrigatorias, names(dados))

if (length(colunas_ausentes) > 0L) {
  stop(
    "Required columns missing from the workbook: ",
    paste(colunas_ausentes, collapse = ", "),
    call. = FALSE
  )
}

colunas_numericas <- setdiff(colunas_obrigatorias, "Municipios")

dados <- dados |>
  mutate(across(all_of(colunas_numericas), as.numeric))

# 2. Validate the dataset -------------------------------------------------

if (nrow(dados) != 184L) {
  stop(
    "Expected 184 municipalities; found ",
    nrow(dados),
    ".",
    call. = FALSE
  )
}

if (anyNA(dados$codigo_ibge) || anyDuplicated(dados$codigo_ibge)) {
  stop(
    "IBGE codes must be present and unique.",
    call. = FALSE
  )
}

if (anyNA(dados$Municipios) || anyDuplicated(dados$Municipios)) {
  stop(
    "Municipality names must be present and unique.",
    call. = FALSE
  )
}

if (anyNA(dados$pop) || any(dados$pop <= 0)) {
  stop(
    "Population values must be present and positive.",
    call. = FALSE
  )
}

if (anyNA(dados$covtotal) || any(dados$covtotal <= 0)) {
  stop(
    "Cumulative case counts must be present and positive.",
    call. = FALSE
  )
}

if (anyNA(dados$obito) || any(dados$obito < 0)) {
  stop(
    "Death counts must be present and nonnegative.",
    call. = FALSE
  )
}

if (any(dados$obito > dados$covtotal)) {
  stop(
    "A municipality has more deaths than confirmed cases.",
    call. = FALSE
  )
}

validar_calculo <- function(observado, esperado, nome, tolerancia = 1e-7) {
  erro <- abs(observado - esperado)
  limite <- tolerancia * pmax(1, abs(esperado))

  if (anyNA(erro) || any(erro > limite)) {
    stop(
      "Column ", nome, " does not match its source values.",
      call. = FALSE
    )
  }
}

validar_calculo(
  observado = dados$cov100k,
  esperado = dados$covtotal / dados$pop * 100000,
  nome = "cov100k"
)

validar_calculo(
  observado = dados$obito100k,
  esperado = dados$obito / dados$pop * 100000,
  nome = "obito100k"
)

validar_calculo(
  observado = dados$letal,
  esperado = dados$obito / dados$covtotal,
  nome = "letal"
)

validar_calculo(
  observado = dados$ideb.fund,
  esperado = (dados$ideb5 + dados$ideb9) / 2,
  nome = "ideb.fund"
)

colunas_finitas <- setdiff(colunas_numericas, "tmi")
if (any(!is.finite(as.matrix(dados[colunas_finitas])))) {
  stop(
    "A numeric analysis column contains a non-finite value.",
    call. = FALSE
  )
}

variaveis_binarias <- c("metrop", "semiarido", "aliadogov")

for (variavel in variaveis_binarias) {
  valor <- dados[[variavel]]
  invalido <- !is.na(valor) & !(valor %in% c(0, 1))

  if (any(invalido)) {
    stop(
      "Dummy variable ", variavel,
      " contains a value other than 0 or 1.",
      call. = FALSE
    )
  }
}

# 3. Transform outcomes and electricity consumption ----------------------

dados <- dados |>
  mutate(
    log_casos100k = log(cov100k),
    log_obitos100k = log((obito + 0.5) / pop * 100000),
    logit_letal = log(
      (obito + 0.5) /
        (covtotal - obito + 0.5)
    ),
    log1p_obitos100k = log1p(obito100k),
    energia100k = con.energia / pop * 100000
  )

if (anyNA(dados$energia100k) || any(dados$energia100k <= 0)) {
  stop(
    "Electricity consumption per 100,000 inhabitants must be positive.",
    call. = FALSE
  )
}

padronizar <- function(x, nome) {
  if (anyNA(x)) {
    stop(
      "Variable ", nome, " contains missing values.",
      call. = FALSE
    )
  }

  desvio <- sd(x)

  if (!is.finite(desvio) || desvio == 0) {
    stop(
      "Variable ", nome, " cannot be standardized.",
      call. = FALSE
    )
  }

  as.numeric((x - mean(x)) / desvio)
}

desfechos_principais <- c(
  casos = "log_casos100k",
  obitos = "log_obitos100k",
  letalidade = "logit_letal"
)

dados <- dados |>
  mutate(
    across(
      all_of(unname(desfechos_principais)),
      ~ padronizar(.x, cur_column()),
      .names = "z_{.col}"
    )
  )

# 4. Model covariates -----------------------------------------------------

# The main set uses IVS Urban Infrastructure and the mean of the grade 5/9
# IDEB indicators. Aggregate IVS, IDHM, IDM, the other two IVS components,
# and TMI are retained in the source data but excluded from estimation.
# Fiscal covariates enter only the expanded model. Days since the first
# municipal record enters a separate sensitivity analysis of the main set.
covariaveis_principais <- c(
  "pop", "area", "densidade", "pibpc", "desp.saude", "desp.educ",
  "ideb.fund", "ideb3", "tax.agua", "energia100k",
  "tax.hom", "eleitas.fem", "bolsaf", "sus1k", "leitos1k", "prof1k",
  "metrop", "semiarido", "pop.rural", "aliadogov", "votos", "ivs.infra",
  "emprego", "ex.pobr", "imuni", "idosos", "int.circ", "int.resp",
  "int.diab", "int.asma"
)

covariaveis_fiscais <- c("recursofed", "auxem")
covariaveis_temporais <- "dias.primeiro.registro"
covariaveis_ampliadas <- c(covariaveis_principais, covariaveis_fiscais)
covariaveis_com_tempo <- c(covariaveis_principais, covariaveis_temporais)
covariaveis_continuas <- setdiff(
  unique(c(covariaveis_ampliadas, covariaveis_temporais)),
  variaveis_binarias
)

variaveis_retiradas <- c("idhm", "idm", "ivs.capital", "ivs.renda", "ideb5", "ideb9")

if (length(covariaveis_principais) != 30L) {
  stop("The main model must contain exactly 30 covariates.", call. = FALSE)
}

if (length(covariaveis_ampliadas) != 32L) {
  stop("The expanded model must contain exactly 32 covariates.", call. = FALSE)
}

if (length(covariaveis_com_tempo) != 31L) {
  stop("The temporal sensitivity model must contain exactly 31 covariates.", call. = FALSE)
}

if (any(variaveis_retiradas %in% covariaveis_ampliadas)) {
  stop("An excluded variable is still present in the model.", call. = FALSE)
}

covariaveis_usadas <- unique(c(covariaveis_ampliadas, covariaveis_temporais))
faltantes_por_variavel <- colSums(is.na(dados[covariaveis_usadas]))
faltantes_por_variavel <- faltantes_por_variavel[faltantes_por_variavel > 0L]

if (length(faltantes_por_variavel) > 0L) {
  stop(
    "Model covariates contain NA: ",
    paste(names(faltantes_por_variavel), collapse = ", "),
    call. = FALSE
  )
}

dados <- dados |>
  mutate(
    across(
      all_of(covariaveis_continuas),
      ~ padronizar(.x, cur_column()),
      .names = "z_{.col}"
    )
  )

nomes_para_modelo <- function(variaveis) {
  ifelse(
    variaveis %in% variaveis_binarias,
    variaveis,
    paste0("z_", variaveis)
  )
}

covariaveis_principais_modelo <- nomes_para_modelo(covariaveis_principais)
covariaveis_ampliadas_modelo <- nomes_para_modelo(covariaveis_ampliadas)
covariaveis_com_tempo_modelo <- nomes_para_modelo(covariaveis_com_tempo)

# 5. Prepare BMS inputs ---------------------------------------------------

montar_dados_bms <- function(base, desfecho, covariaveis) {
  dados_bms <- base |>
    select(all_of(c(desfecho, covariaveis)))

  if (anyNA(dados_bms)) {
    stop(
      "The BMS input contains missing values.",
      call. = FALSE
    )
  }

  if (any(!is.finite(as.matrix(dados_bms)))) {
    stop(
      "The BMS input contains an infinite or nonnumeric value.",
      call. = FALSE
    )
  }

  dados_bms
}

preparar_modelos <- function(base, desfechos, covariaveis) {
  modelos <- lapply(
    unname(desfechos),
    function(desfecho) {
      montar_dados_bms(
        base = base,
        desfecho = paste0("z_", desfecho),
        covariaveis = covariaveis
      )
    }
  )

  names(modelos) <- names(desfechos)
  modelos
}

# Each input has the dependent variable first, followed by covariates.
modelos_principais <- preparar_modelos(
  base = dados,
  desfechos = desfechos_principais,
  covariaveis = covariaveis_principais_modelo
)

modelos_ampliados <- preparar_modelos(
  base = dados,
  desfechos = desfechos_principais,
  covariaveis = covariaveis_ampliadas_modelo
)

modelos_com_tempo <- preparar_modelos(
  base = dados,
  desfechos = desfechos_principais,
  covariaveis = covariaveis_com_tempo_modelo
)

# 6. Sensitivity analyses -------------------------------------------------

# Alternative scales for the three dependent variables.
dados <- dados |>
  mutate(
    z_casos100k_nivel = padronizar(cov100k, "cov100k"),
    z_obitos100k_nivel = padronizar(obito100k, "obito100k"),
    z_log1p_obitos100k = padronizar(log1p_obitos100k, "log1p_obitos100k"),
    z_letal_nivel = padronizar(letal, "letal")
  )

desfechos_sensibilidade <- c(
  casos_nivel = "z_casos100k_nivel",
  obitos_nivel = "z_obitos100k_nivel",
  obitos_log1p = "z_log1p_obitos100k",
  letalidade_nivel = "z_letal_nivel"
)

modelos_sensibilidade_escala <- lapply(
  desfechos_sensibilidade,
  function(desfecho) {
    montar_dados_bms(
      base = dados,
      desfecho = desfecho,
      covariaveis = covariaveis_principais_modelo
    )
  }
)

# Exclude Fortaleza by IBGE code rather than name spelling.
dados_sem_fortaleza <- dados |>
  filter(codigo_ibge != 2304400)

modelos_sem_fortaleza <- preparar_modelos(
  base = dados_sem_fortaleza,
  desfechos = desfechos_principais,
  covariaveis = covariaveis_principais_modelo
)

# Population-cutoff sensitivity samples retain full-sample standardization
# so coefficient scales remain comparable across specifications.
dados_pop_10mil <- dados |>
  filter(pop >= 10000)

dados_pop_20mil <- dados |>
  filter(pop >= 20000)

if (nrow(dados_pop_10mil) != 164L || nrow(dados_pop_20mil) != 101L) {
  stop(
    "Population cutoffs did not yield the expected samples (164 and 101 municipalities).",
    call. = FALSE
  )
}

modelos_pop_10mil <- preparar_modelos(
  base = dados_pop_10mil,
  desfechos = desfechos_principais,
  covariaveis = covariaveis_principais_modelo
)

modelos_pop_20mil <- preparar_modelos(
  base = dados_pop_20mil,
  desfechos = desfechos_principais,
  covariaveis = covariaveis_principais_modelo
)

cat(
  "Preparation complete.\n",
  "Municipalities in the main model:", nrow(dados), "\n",
  "Covariates in the main model:", length(covariaveis_principais), "\n",
  "Covariates in the expanded model:", length(covariaveis_ampliadas), "\n",
  "Covariates in the temporal sensitivity model:", length(covariaveis_com_tempo), "\n",
  "Municipalities at the 10,000 population cutoff:", nrow(dados_pop_10mil), "\n",
  "Municipalities at the 20,000 population cutoff:", nrow(dados_pop_20mil), "\n"
)

# 7. Estimation settings --------------------------------------------------

modo_teste <- identical(Sys.getenv("COVID_BMA_MODO_TESTE"), "1")
reprocessar_existentes <- identical(Sys.getenv("COVID_BMA_REPROCESSAR"), "1")
somente_preparar <- identical(Sys.getenv("COVID_BMA_SOMENTE_PREPARAR"), "1")

configuracao_principal <- list(
  burn = if (modo_teste) 100 else 200000,
  iter = if (modo_teste) 500 else 1000000,
  nmodel = if (modo_teste) 100 else 2000,
  mcmc = "bd",
  g = "hyper=UIP",
  mprior = "random",
  force.full.ols = TRUE
)

configuracao_sensibilidade <- configuracao_principal
configuracao_sensibilidade$burn <- if (modo_teste) 100 else 100000
configuracao_sensibilidade$iter <- if (modo_teste) 500 else 500000

sementes_principais <- c(20200731, 20200801)
sementes_sensibilidade <- 20200731

pasta_raiz_resultados <- Sys.getenv("COVID_BMA_SAIDA", file.path(
  "resultados", if (modo_teste) "teste" else "estimacoes"
))
dir.create(pasta_raiz_resultados, recursive = TRUE, showWarnings = FALSE)

resultado_compativel <- function(objeto, dados_modelo, configuracao) {
  x_salvo <- tryCatch(
    objeto$arguments$X.data,
    error = function(e) NULL
  )

  if (is.null(x_salvo)) {
    return(FALSE)
  }

  if (!identical(dim(x_salvo), dim(dados_modelo)) ||
      !identical(colnames(x_salvo), colnames(dados_modelo))) {
    return(FALSE)
  }

  diferenca <- max(abs(as.matrix(x_salvo) - as.matrix(dados_modelo)))
  argumentos <- as.list(objeto$arguments)
  configuracao_igual <- identical(as.numeric(argumentos$burn), as.numeric(configuracao$burn)) &&
    identical(as.numeric(argumentos$iter), as.numeric(configuracao$iter)) &&
    identical(as.numeric(argumentos$nmodel), as.numeric(configuracao$nmodel)) &&
    identical(as.character(argumentos$mcmc), as.character(configuracao$mcmc)) &&
    identical(as.character(argumentos$g), as.character(configuracao$g)) &&
    identical(as.character(argumentos$mprior), as.character(configuracao$mprior)) &&
    identical(isTRUE(argumentos$force.full.ols), isTRUE(configuracao$force.full.ols))

  is.finite(diferenca) && diferenca <= 1e-10 && configuracao_igual
}

estimar_cadeias <- function(dados_modelo, pasta, configuracao, sementes) {
  dir.create(pasta, recursive = TRUE, showWarnings = FALSE)
  numero_covariaveis <- ncol(dados_modelo) - 1L
  cadeias <- vector("list", length(sementes))
  names(cadeias) <- paste0("cadeia_", seq_along(sementes))
  valores_iniciais <- rep(list(NA), length(sementes))
  if (length(sementes) > 1L) {
    valores_iniciais[[2]] <- 0
  }

  for (i in seq_along(sementes)) {
    arquivo <- file.path(pasta, paste0("cadeia_", i, ".rds"))

    objeto_salvo <- NULL
    cache_compativel <- FALSE

    if (file.exists(arquivo) && file.mtime(arquivo) >= file.mtime(arquivo_dados)) {
      objeto_salvo <- tryCatch(readRDS(arquivo), error = function(e) NULL)
      cache_compativel <- !is.null(objeto_salvo) &&
        resultado_compativel(objeto_salvo, dados_modelo, configuracao)
    }

    if (cache_compativel && !reprocessar_existentes) {
      message("Loading existing result: ", arquivo)
      cadeias[[i]] <- objeto_salvo
      next
    }

    if (file.exists(arquivo) && !cache_compativel) {
      message(
        "The saved result does not match the current matrix. ",
        "The chain will be reestimated: ", arquivo
      )
    }

    set.seed(sementes[i])
    message(
      "Estimating ", basename(pasta), ", chain ", i,
      " of ", length(sementes), "."
    )

    cadeias[[i]] <- BMS::bms(
      X.data = dados_modelo,
      burn = configuracao$burn,
      iter = configuracao$iter,
      nmodel = configuracao$nmodel,
      mcmc = configuracao$mcmc,
      g = configuracao$g,
      mprior = configuracao$mprior,
      mprior.size = numero_covariaveis / 2,
      start.value = valores_iniciais[[i]],
      user.int = FALSE,
      g.stats = TRUE,
      logfile = file.path(pasta, paste0("cadeia_", i, ".log")),
      logstep = if (modo_teste) 100 else 100000,
      force.full.ols = configuracao$force.full.ols
    )

    saveRDS(cadeias[[i]], arquivo)
  }

  saveRDS(cadeias, file.path(pasta, "todas_as_cadeias.rds"))
  cadeias
}

extrair_coeficientes <- function(cadeias, coluna) {
  tabelas <- lapply(
    cadeias,
    coef,
    exact = FALSE,
    order.by.pip = FALSE
  )
  variaveis <- rownames(tabelas[[1]])
  valores <- sapply(
    tabelas,
    function(tabela) tabela[variaveis, coluna]
  )
  if (is.null(dim(valores))) {
    valores <- matrix(valores, ncol = 1)
  }
  rownames(valores) <- variaveis
  valores
}

comparar_cadeias <- function(cadeias) {
  pips <- extrair_coeficientes(cadeias, "PIP")
  medias <- extrair_coeficientes(cadeias, "Post Mean")
  sinais <- extrair_coeficientes(cadeias, "Cond.Pos.Sign")
  colnames(pips) <- paste0("PIP_", names(cadeias))

  data.frame(
    variavel = rownames(pips),
    pips,
    PIP_media = rowMeans(pips),
    PIP_minima = apply(pips, 1, min),
    PIP_maxima = apply(pips, 1, max),
    diferenca_maxima = apply(pips, 1, function(x) max(x) - min(x)),
    Post_Mean_media = rowMeans(medias),
    Cond_Pos_Sign_media = rowMeans(sinais),
    row.names = NULL
  ) |>
    dplyr::arrange(dplyr::desc(PIP_media))
}

consolidar_cadeias <- function(cadeias) {
  pips <- extrair_coeficientes(cadeias, "PIP")
  medias <- extrair_coeficientes(cadeias, "Post Mean")
  desvios <- extrair_coeficientes(cadeias, "Post SD")
  sinais <- extrair_coeficientes(cadeias, "Cond.Pos.Sign")

  draws <- vapply(
    cadeias,
    function(x) as.numeric(summary(x)[["Draws"]]),
    numeric(1)
  )
  pesos <- draws / sum(draws)
  pesos_matriz <- matrix(
    pesos,
    nrow = nrow(pips),
    ncol = length(pesos),
    byrow = TRUE
  )

  pip_final <- as.numeric(pips %*% pesos)
  media_final <- as.numeric(medias %*% pesos)
  segundo_momento <- as.numeric((desvios^2 + medias^2) %*% pesos)
  peso_inclusao <- rowSums(pips * pesos_matriz)

  resultado <- data.frame(
    variavel = rownames(pips),
    PIP = pip_final,
    post_mean = media_final,
    post_sd = sqrt(pmax(segundo_momento - media_final^2, 0)),
    cond_pos_sign = ifelse(
      peso_inclusao > 0,
      rowSums(pips * sinais * pesos_matriz) / peso_inclusao,
      NA_real_
    ),
    PIP_maior_igual_05 = pip_final >= 0.5,
    check.names = FALSE
  )
  names(resultado)[3:5] <- c("Post Mean", "Post SD", "Cond.Pos.Sign")

  resultado |>
    dplyr::arrange(dplyr::desc(PIP))
}

diagnosticar_cadeias <- function(cadeias, especificacao, desfecho, sementes) {
  dplyr::bind_rows(lapply(seq_along(cadeias), function(i) {
    resumo <- summary(cadeias[[i]])
    data.frame(
      especificacao = especificacao,
      desfecho = desfecho,
      cadeia = names(cadeias)[i],
      semente = sementes[i],
      tamanho_medio_modelo = resumo[["Mean no. regressors"]],
      iteracoes = resumo[["Draws"]],
      burn_in = resumo[["Burnins"]],
      tempo = resumo[["Time"]],
      modelos_visitados = resumo[["No. models visited"]],
      percentual_visitado = resumo[["% visited"]],
      corr_PMP = resumo[["Corr PMP"]],
      prior_modelo = resumo[["Model Prior"]],
      g_prior = resumo[["g-Prior"]],
      shrinkage = resumo[["Shrinkage-Stats"]]
    )
  }))
}

avaliar_diagnostico <- function(diagnosticos, comparacao) {
  if (modo_teste) {
    return("teste")
  }
  corr_minima <- min(as.numeric(diagnosticos$corr_PMP), na.rm = TRUE)
  criterio_pip <- TRUE
  if (nrow(diagnosticos) > 1L) {
    diferenca_pip <- max(comparacao$diferenca_maxima, na.rm = TRUE)
    criterio_pip <- diferenca_pip <= 0.05
  }
  if (corr_minima >= 0.95 && criterio_pip) {
    "aprovado"
  } else {
    "revisar"
  }
}

executar_grupo <- function(modelos, nome, configuracao, sementes) {
  pasta <- file.path(pasta_raiz_resultados, nome)
  pasta_tabelas <- file.path(pasta, "tabelas_tecnicas")
  dir.create(pasta_tabelas, recursive = TRUE, showWarnings = FALSE)

  objetos <- lapply(names(modelos), function(desfecho) {
    estimar_cadeias(
      modelos[[desfecho]],
      file.path(pasta, desfecho),
      configuracao,
      sementes
    )
  })
  names(objetos) <- names(modelos)

  comparacoes <- lapply(objetos, comparar_cadeias)
  consolidados <- lapply(objetos, consolidar_cadeias)
  diagnosticos <- dplyr::bind_rows(
    lapply(names(objetos), function(desfecho) {
      diagnosticar_cadeias(
        objetos[[desfecho]],
        nome,
        desfecho,
        sementes
      )
    })
  )

  resumos <- dplyr::bind_rows(lapply(names(objetos), function(desfecho) {
    diagnostico <- diagnosticos[
      diagnosticos$desfecho == desfecho,
      ,
      drop = FALSE
    ]
    data.frame(
      especificacao = nome,
      desfecho = desfecho,
      cadeias = nrow(diagnostico),
      iteracoes_retidas = sum(as.numeric(diagnostico$iteracoes)),
      burn_in_total = sum(as.numeric(diagnostico$burn_in)),
      corr_PMP_minima = min(as.numeric(diagnostico$corr_PMP)),
      corr_PMP_maxima = max(as.numeric(diagnostico$corr_PMP)),
      maior_diferenca_PIP = if (nrow(diagnostico) > 1L) {
        max(comparacoes[[desfecho]]$diferenca_maxima)
      } else {
        NA_real_
      },
      modelos_visitados_por_cadeia = paste(
        diagnostico$modelos_visitados,
        collapse = ";"
      ),
      status = avaliar_diagnostico(
        diagnostico,
        comparacoes[[desfecho]]
      )
    )
  }))

  for (desfecho in names(objetos)) {
    write.csv(
      comparacoes[[desfecho]],
      file.path(pasta_tabelas, paste0("comparacao_", desfecho, ".csv")),
      row.names = FALSE,
      fileEncoding = "UTF-8"
    )
    write.csv(
      consolidados[[desfecho]],
      file.path(pasta_tabelas, paste0("resultado_", desfecho, ".csv")),
      row.names = FALSE,
      fileEncoding = "UTF-8"
    )
  }

  write.csv(
    diagnosticos,
    file.path(pasta_tabelas, "diagnosticos_cadeias.csv"),
    row.names = FALSE,
    fileEncoding = "UTF-8"
  )
  write.csv(
    resumos,
    file.path(pasta_tabelas, "resumo_diagnosticos.csv"),
    row.names = FALSE,
    fileEncoding = "UTF-8"
  )
  saveRDS(objetos, file.path(pasta, "resultados.rds"))

  list(
    objetos = objetos,
    comparacoes = comparacoes,
    consolidados = consolidados,
    diagnosticos = diagnosticos,
    resumo = resumos
  )
}

# 9. Inventory and exact estimation matrices -----------------------------

exportar_matrizes <- function(modelos, base, nome_grupo) {
  pasta <- file.path(pasta_raiz_resultados, "matrizes_estimacao", nome_grupo)
  dir.create(pasta, recursive = TRUE, showWarnings = FALSE)

  if (any(vapply(modelos, nrow, integer(1)) != nrow(base))) {
    stop("The matrix and identifier table have different sizes.", call. = FALSE)
  }

  for (desfecho in names(modelos)) {
    matriz <- data.frame(
      codigo_ibge = sprintf("%.0f", base$codigo_ibge),
      municipio = base$Municipios,
      modelos[[desfecho]],
      check.names = FALSE
    )
    write.csv(
      matriz,
      file.path(pasta, paste0(desfecho, ".csv")),
      row.names = FALSE,
      fileEncoding = "UTF-8"
    )
  }

  saveRDS(modelos, file.path(pasta, "matrizes_bms.rds"))
}

inventariar_grupo <- function(modelos, nome, configuracao, sementes) {
  dplyr::bind_rows(lapply(names(modelos), function(desfecho) {
    k <- ncol(modelos[[desfecho]]) - 1L
    data.frame(
      especificacao = nome,
      desfecho = desfecho,
      municipios = nrow(modelos[[desfecho]]),
      covariaveis = k,
      tamanho_esperado_prior = k / 2,
      cadeias = length(sementes),
      sementes = paste(sementes, collapse = ";"),
      burn_in_por_cadeia = configuracao$burn,
      iteracoes_por_cadeia = configuracao$iter,
      g_prior = configuracao$g,
      prior_modelo = configuracao$mprior,
      stringsAsFactors = FALSE
    )
  }))
}

grupos_programados <- list(
  principal = list(modelos_principais, dados, configuracao_principal, sementes_principais),
  ampliado = list(modelos_ampliados, dados, configuracao_principal, sementes_principais),
  sensibilidade_escalas = list(
    modelos_sensibilidade_escala, dados, configuracao_sensibilidade, sementes_sensibilidade
  ),
  sensibilidade_sem_fortaleza = list(
    modelos_sem_fortaleza, dados_sem_fortaleza,
    configuracao_sensibilidade, sementes_sensibilidade
  ),
  sensibilidade_pop_10mil = list(
    modelos_pop_10mil, dados_pop_10mil,
    configuracao_sensibilidade, sementes_sensibilidade
  ),
  sensibilidade_pop_20mil = list(
    modelos_pop_20mil, dados_pop_20mil,
    configuracao_sensibilidade, sementes_sensibilidade
  ),
  sensibilidade_tempo_exposicao = list(
    modelos_com_tempo, dados,
    configuracao_sensibilidade, sementes_sensibilidade
  )
)

for (nome_grupo in names(grupos_programados)) {
  grupo <- grupos_programados[[nome_grupo]]
  exportar_matrizes(grupo[[1]], grupo[[2]], nome_grupo)
}

inventario_estimacoes <- dplyr::bind_rows(lapply(
  names(grupos_programados),
  function(nome_grupo) {
    grupo <- grupos_programados[[nome_grupo]]
    inventariar_grupo(grupo[[1]], nome_grupo, grupo[[3]], grupo[[4]])
  }
))

write.csv(
  inventario_estimacoes,
  file.path(pasta_raiz_resultados, "inventario_estimacoes.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

# 10. Run the models ------------------------------------------------------

if (somente_preparar) {
  message(
    "Preparation complete without estimation. Check inventario_estimacoes.csv and ",
    "the matrizes_estimacao directory before running full estimation."
  )
} else {

resultado_principal <- executar_grupo(
  modelos_principais,
  "principal",
  configuracao_principal,
  sementes_principais
)

resultado_ampliado <- executar_grupo(
  modelos_ampliados,
  "ampliado",
  configuracao_principal,
  sementes_principais
)

resultado_escalas <- executar_grupo(
  modelos_sensibilidade_escala,
  "sensibilidade_escalas",
  configuracao_sensibilidade,
  sementes_sensibilidade
)

resultado_sem_fortaleza <- executar_grupo(
  modelos_sem_fortaleza,
  "sensibilidade_sem_fortaleza",
  configuracao_sensibilidade,
  sementes_sensibilidade
)

resultado_pop_10mil <- executar_grupo(
  modelos_pop_10mil,
  "sensibilidade_pop_10mil",
  configuracao_sensibilidade,
  sementes_sensibilidade
)

resultado_pop_20mil <- executar_grupo(
  modelos_pop_20mil,
  "sensibilidade_pop_20mil",
  configuracao_sensibilidade,
  sementes_sensibilidade
)

resultado_tempo_exposicao <- executar_grupo(
  modelos_com_tempo,
  "sensibilidade_tempo_exposicao",
  configuracao_sensibilidade,
  sementes_sensibilidade
)

resumos_gerais <- dplyr::bind_rows(
  resultado_principal$resumo,
  resultado_ampliado$resumo,
  resultado_escalas$resumo,
  resultado_sem_fortaleza$resumo,
  resultado_pop_10mil$resumo,
  resultado_pop_20mil$resumo,
  resultado_tempo_exposicao$resumo
)

write.csv(
  resumos_gerais,
  file.path(pasta_raiz_resultados, "resumo_geral_diagnosticos.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

saveRDS(
  list(
    principal = resultado_principal,
    ampliado = resultado_ampliado,
    escalas = resultado_escalas,
    sem_fortaleza = resultado_sem_fortaleza,
    pop_10mil = resultado_pop_10mil,
    pop_20mil = resultado_pop_20mil,
    tempo_exposicao = resultado_tempo_exposicao
  ),
  file.path(pasta_raiz_resultados, "todos_os_resultados.rds")
)

if (interactive()) {
  View(resumos_gerais, title = "Resumo dos diagnosticos")
}

if (!modo_teste && any(resumos_gerais$status == "revisar")) {
  warning(
    "Uma ou mais estimacoes exigem revisao dos diagnosticos. ",
    "Consulte resumo_geral_diagnosticos.csv.",
    call. = FALSE
  )
} else {
  message("All scheduled estimations have finished.")
}

}

writeLines(
  capture.output(sessionInfo()),
  file.path(pasta_raiz_resultados, "sessionInfo.txt")
)
