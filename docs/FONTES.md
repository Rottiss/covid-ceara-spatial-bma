# Fontes e transformações

[English](SOURCES.md) | Português | [README](../README.pt-BR.md)

A base consolidada está em `estimacoes/dadosv6.xlsx`. O script `covid_bma.R` lê exclusivamente a aba `raw`. As demais abas auxiliam a consulta e não entram como dados de estimação.

As instituições e os períodos abaixo correspondem à documentação da base utilizada na pesquisa. A coleta nos portais de origem não é automatizada neste repositório. A reprodução computacional começa na planilha consolidada; não há um comando que reconstrua essa planilha a partir de todos os arquivos administrativos originais.

| Variáveis | Referência | Fonte registrada |
|---|---|---|
| Casos e óbitos | Acumulado até 31/07/2020 | Secretarias Estaduais de Saúde, dados consolidados pelo Brasil.IO |
| População; proporção de pessoas com 60 anos ou mais | 2019 | DATASUS |
| Área; densidade demográfica | 2019 | IBGE |
| PIB per capita | 2018 | IBGE |
| Proporção rural | 2010 | IBGE |
| Empregos formais | 2019 | IPECE |
| Extrema pobreza | 2010 | IPECE |
| Beneficiários do Bolsa Família | 2020 | Ministério da Cidadania, SENARC |
| IVS Infraestrutura Urbana | 2010 | IPEA, Atlas da Vulnerabilidade Social |
| Região Metropolitana de Fortaleza | Delimitação de 2018 | IPECE |
| Semiárido | 2020 | IPECE |
| Alinhamento partidário | Classificação de 2019 | TRE-CE e TSE |
| Votos válidos e proporção de mulheres eleitas | Eleição de 2016 | TRE-CE e TSE |
| Homicídios | 2019 | SIM, DATASUS |
| Despesas com saúde e saneamento; educação e cultura | 2018 | IPEAdata, com dados da STN |
| IDEB do ensino fundamental e do ensino médio | 2019 | INEP |
| Cobertura urbana de água | 2019 | IPECE |
| Consumo de energia elétrica | 2019 | IPECE e Coelce |
| Imunização em menores de um ano | 2019 | IPECE e SESA |
| Unidades de saúde vinculadas ao SUS | 2019 | SESA |
| Leitos e profissionais vinculados ao SUS | 2019 | IPECE e SESA |
| Internações por doenças circulatórias, respiratórias, diabetes e asma | 2019 | SIH, DATASUS |
| Transferências federais para combate à COVID-19 | Março a julho de 2020 | Portal da Transparência, CGU |
| Auxílio emergencial | Abril a julho de 2020 | Portal da Transparência, CGU |
| Dias desde o primeiro registro municipal | Primeiro registro até 31/07/2020 | Secretarias Estaduais de Saúde, dados consolidados pelo Brasil.IO |
| Geometria dos municípios | 2020 | Malha do IBGE obtida por meio do pacote geobr |

## Transformações usadas nas estimações

O código em `estimacoes/covid_bma.R` define e verifica as transformações. As principais são:

- Casos: `log(covtotal / pop * 100000)`.
- Óbitos: `log((obito + 0.5) / pop * 100000)`.
- Letalidade: `log((obito + 0.5) / (covtotal - obito + 0.5))`.
- IDEB do ensino fundamental: média simples de `ideb5` e `ideb9`.
- Energia: `con.energia / pop * 100000`.
- Padronização: `(x - mean(x)) / sd(x)`, aplicada às variáveis contínuas conforme cada especificação. As três dummies permanecem em 0/1.

As covariáveis do modelo principal estão listadas em `covariaveis_principais` no script. A planilha preserva campos que não entram nesse modelo. A dummy `metrop` identifica a Região Metropolitana de Fortaleza; os demais municípios formam a categoria de referência.

O comando `Rscript reproduzir.R preparar` reconstrói as matrizes e as compara com as referências incluídas. O comando `Rscript reproduzir.R geometria` alinha a malha por código IBGE, reconstrói as duas vizinhanças e compara os pesos com os originais.
