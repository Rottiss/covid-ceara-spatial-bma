# Gera o apendice do TCC com as PIP das 30 covariaveis sob SAR e SEM, nas duas
# matrizes de vizinhanca. Le as saidas e GRAVA um unico .tex na copia do TCC.
#
# O MQO tem coluna unica porque nao usa W: as duas execucoes rodam o mesmo
# modelo com as mesmas sementes e devolvem PIP identicas. O script confere isso
# em vez de supor.
#
# Os rotulos sao os do proprio TCC, e o script verifica que cada um aparece no
# Apendice B antes de gravar, para nao criar nome novo para covariavel que ja
# tem nome no documento.
#
# Saida: manuscrito/3-pos-textuais/apendices/apendice-g.tex

source(if (file.exists("../caminhos.R")) "../caminhos.R" else "caminhos.R")

S <- file.path(ESPACIAL, "saida")
TCC <- MANUSCRITO
APB <- file.path(TCC, "3-pos-textuais/apendices/apendice-b.tex")
DEST <- file.path(TCC, "3-pos-textuais/apendices/apendice-g.tex")

ROTULO <- c(
  z_idosos = "Proporção da população com 60 anos ou mais",
  z_leitos1k = "Leitos vinculados ao SUS",
  metrop = "Município integrante da Região Metropolitana de Fortaleza",
  semiarido = "Município integrante do semiárido",
  z_ideb.fund = "IDEB do ensino fundamental",
  z_int.circ = "Internações por doenças do aparelho circulatório",
  z_bolsaf = "População beneficiária do Programa Bolsa Família",
  z_pop.rural = "Proporção da população rural",
  z_int.resp = "Internações por doenças do aparelho respiratório",
  z_prof1k = "Profissionais de saúde vinculados ao SUS",
  z_energia100k = "Consumo de energia elétrica por 100 mil habitantes",
  z_ex.pobr = "Proporção da população em extrema pobreza",
  z_tax.hom = "Taxa de homicídios",
  z_desp.educ = "Despesa municipal com educação e cultura",
  z_int.asma = "Internações por asma",
  z_ivs.infra = "IVS Infraestrutura Urbana",
  z_imuni = "Cobertura média de imunização em menores de um ano",
  z_desp.saude = "Despesa municipal com saúde e saneamento",
  z_pibpc = "Produto Interno Bruto per capita",
  z_tax.agua = "Cobertura urbana de abastecimento de água",
  z_sus1k = "Unidades de saúde vinculadas ao SUS",
  z_int.diab = "Internações por diabetes mellitus",
  z_eleitas.fem = "Proporção de mulheres eleitas",
  z_emprego = "Proporção da população em empregos formais",
  z_pop = "População estimada",
  aliadogov = "Alinhamento partidário do prefeito com o governo federal",
  z_densidade = "Densidade demográfica",
  z_area = "Área territorial",
  z_ideb3 = "IDEB do 3º ano do ensino médio",
  z_votos = "Proporção de votos válidos do prefeito eleito")

apb <- readLines(APB, encoding = "UTF-8", warn = FALSE)
apb <- paste(apb, collapse = "\n")
falta <- ROTULO[!vapply(ROTULO, function(x) grepl(x, apb, fixed = TRUE), logical(1))]
stopifnot("rotulo ausente do Apendice B" = length(falta) == 0)

vir <- function(x) formatC(x, format = "f", digits = 4, decimal.mark = ",")

DESF <- c(casos = "casos de COVID-19 por 100 mil habitantes",
          obitos = "óbitos por COVID-19 por 100 mil habitantes",
          letalidade = "letalidade por COVID-19")

linhas <- character(0)
add <- function(...) linhas <<- c(linhas, sprintf(...))

add("\\apendice{Probabilidades de inclusão sob SAR e SEM nas duas matrizes de vizinhança}")
add("\\label{ap:vizinhanca-espacial}")
add("")
add(paste("As Tabelas~\\ref{tab:apG-casos} a~\\ref{tab:apG-letalidade} reúnem as PIP das 30 covariáveis sob as especificações",
          "espaciais da Seção~\\ref{sec:mc3-espacial}, com as matrizes Queen e Rook. A coluna MQO",
          "é única porque o modelo de Mínimos Quadrados Ordinários não utiliza a matriz de pesos:",
          "as duas execuções percorrem o mesmo espaço de modelos com as mesmas sementes e devolvem",
          "PIP idênticas. Esse MQO é a reestimação pela cadeia própria descrita na",
          "Seção~\\ref{sec:mc3-espacial}, não o resultado do pacote BMS apresentado no",
          "Apêndice~\\ref{ap:resultados-completos}. As covariáveis aparecem em ordem decrescente",
          "de PIP sob MQO."))
add("")

for (d in names(DESF)) {
  q <- readRDS(file.path(S, "grade", sprintf("%s_queen.rds", d)))$tabela
  r <- readRDS(file.path(S, "grade", sprintf("%s_rook.rds", d)))$tabela
  r <- r[match(q$covariavel, r$covariavel), ]
  sq <- read.csv(file.path(S, sprintf("sem_hyperg_%s_queen.csv", d)), stringsAsFactors = FALSE)
  sr <- read.csv(file.path(S, sprintf("sem_hyperg_%s_rook.csv", d)), stringsAsFactors = FALSE)
  sq <- sq[match(q$covariavel, sq$covariavel), ]
  sr <- sr[match(q$covariavel, sr$covariavel), ]
  stopifnot(identical(q$covariavel, r$covariavel), identical(q$covariavel, sq$covariavel),
            identical(q$covariavel, sr$covariavel))
  # o MQO nao pode depender da matriz
  stopifnot("MQO difere entre matrizes" = max(abs(q$mqo - r$mqo)) < 1e-12)

  o <- order(-q$mqo)
  add("\\scriptsize")
  # Guia de normalizacao da UFC (2022), 4.10 e: nas paginas seguintes o titulo e
  # repetido com "(continuacao)" ou "(conclusao)" pela macro \tabelacontinuada
  # (manuscrito/lib/preambulo.tex); o pe das paginas que continuam traz "(continua)",
  # sem traco de fechamento; o rotulo "-fim" marca a ultima linha da tabela.
  add("\\tabelacontinuada{Probabilidades de inclusão sob SAR e SEM, Queen e Rook, %s}{tab:apG-%s-fim}",
      DESF[[d]], d)
  add("\\begin{longtable}{p{5.6cm}rrrrr}")
  add("    \\caption{Probabilidades de inclusão sob SAR e SEM, Queen e Rook, %s}\\label{tab:apG-%s} \\\\",
      DESF[[d]], d)
  cab <- paste("    Covariável & MQO & SAR Queen & SAR Rook & SEM Queen & SEM Rook \\\\")
  add("    \\toprule"); add(cab); add("    \\midrule"); add("    \\endfirsthead")
  add("    \\toprule"); add(cab); add("    \\midrule"); add("    \\endhead")
  add("    \\multicolumn{6}{r}{\\small (continua)} \\\\")
  add("    \\endfoot"); add("    \\bottomrule"); add("    \\endlastfoot")
  add("")
  for (i in o)
    add("    %s%s & %s & %s & %s & %s & %s \\\\",
        ROTULO[[q$covariavel[i]]],
        if (i == o[length(o)]) sprintf("\\label{tab:apG-%s-fim}", d) else "",
        vir(q$mqo[i]), vir(q$sar[i]), vir(r$sar[i]),
        vir(sq$sem[i]), vir(sr$sem[i]))
  add("\\end{longtable}")
  add("\\Fonte{Elaboração própria com base nos resultados da pesquisa.}")
  add("\\normalsize")
  add("")
}

writeLines(linhas, DEST, useBytes = FALSE)
cat("gravado:", DEST, "\n")
cat("linhas:", length(linhas), "| 3 tabelas de 30 covariaveis\n")
