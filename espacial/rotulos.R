# Nomes das covariáveis para uso nas figuras, em português e em inglês.
#
# Os nomes de coluna da base são abreviações de trabalho e não servem para
# figura nenhuma. As descrições vêm de COVID_BMA/DICIONARIO_DADOS_FASE_1.md,
# encurtadas para caber num rótulo. O prefixo z_ marca padronização e não
# aparece nos rótulos, porque todas as contínuas estão padronizadas.

ROTULOS <- data.frame(
  variavel = c(
    "z_pop", "z_area", "z_densidade", "z_pibpc", "z_desp.saude", "z_desp.educ",
    "z_ideb.fund", "z_ideb3", "z_tax.agua", "z_energia100k", "z_tax.hom",
    "z_eleitas.fem", "z_bolsaf", "z_sus1k", "z_leitos1k", "z_prof1k",
    "metrop", "semiarido", "z_pop.rural", "aliadogov", "z_votos",
    "z_ivs.infra", "z_emprego", "z_ex.pobr", "z_imuni", "z_idosos",
    "z_int.circ", "z_int.resp", "z_int.diab", "z_int.asma"),
  pt = c(
    "população", "área", "densidade demográfica", "PIB per capita",
    "despesa com saúde", "despesa com educação", "IDEB fundamental",
    "IDEB ensino médio", "cobertura de água", "consumo de energia",
    "taxa de homicídios", "mulheres eleitas", "Bolsa Família",
    "unidades do SUS", "leitos do SUS", "profissionais de saúde",
    "Região Metropolitana", "semiárido", "população rural",
    "prefeito aliado ao governo", "votos do prefeito", "IVS infraestrutura",
    "emprego formal", "extrema pobreza", "cobertura de imunização",
    "proporção de idosos", "internações circulatórias",
    "internações respiratórias", "internações por diabetes",
    "internações por asma"),
  # Nomes próprios de instituições e índices brasileiros não se traduzem: SUS,
  # IDEB, IVS, RMF e Bolsa Família aparecem assim na literatura em inglês sobre
  # o Brasil. Traduzi-los apaga a referência e inventa siglas que já pertencem a
  # outros índices, como SVI, que é o da agência de saúde dos Estados Unidos.
  en = c(
    "population", "area", "population density", "GDP per capita",
    "health spending", "education spending", "primary-school IDEB",
    "secondary-school IDEB", "water coverage", "electricity use",
    "homicide rate", "women elected", "Bolsa Família coverage",
    "SUS facilities", "SUS hospital beds", "SUS health workers",
    "Fortaleza metro region (RMF)", "semi-arid region", "rural population",
    "mayor aligned with federal govt", "mayor's vote share",
    "IVS urban infrastructure", "formal employment", "extreme poverty",
    "immunization coverage", "share aged 60+", "circulatory admissions",
    "respiratory admissions", "diabetes admissions", "asthma admissions"),
  stringsAsFactors = FALSE)

rotular <- function(v, idioma = "pt") {
  i <- match(v, ROTULOS$variavel)
  if (anyNA(i)) stop("sem rótulo para: ", paste(v[is.na(i)], collapse = ", "))
  ROTULOS[[idioma]][i]
}

# Rótulo isolado (eixo de figura, linha de tabela) começa com maiúscula, como nas
# figuras do TCC; o resto fica como está, para preservar SUS, IDEB, IVS e nomes
# próprios. No meio de frase, use rotular() sem esta função.
inicial_maiuscula <- function(s) paste0(toupper(substr(s, 1, 1)), substring(s, 2))
