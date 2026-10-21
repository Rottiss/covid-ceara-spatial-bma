# Inputs and results

English | [Português](DADOS.md) | [README](../README.md)

## Estimation matrices

`casos.csv`, `obitos.csv`, and `letalidade.csv` in `estimacoes/resultados/estimacoes/matrizes_estimacao/principal/` contain the main analysis inputs. Each row is a municipality.

| Field | Meaning |
|---|---|
| `codigo_ibge` | Municipal identifier; read as text |
| `municipio` | Municipality name |
| `z_log_casos100k`, `z_log_obitos100k`, or `z_logit_letal` | Transformed and standardized dependent variable, depending on the file |
| Other columns | The 30 covariates used in estimation |

Continuous variables are standardized using z-scores. The dummy variables `metrop`, `semiarido`, and `aliadogov` retain their 0/1 coding. These matrices are not raw data: their columns already incorporate the estimation transformations.

The consolidated dataset is in `estimacoes/dadosv6.xlsx`; formulas are in `estimacoes/covid_bma.R`. See [SOURCES.md](SOURCES.md) for provenance and transformations. The `preparar` command rebuilds and compares the 22 included matrices covering the main model and sensitivity analyses.

## Variable glossary

Names are preserved to maintain compatibility with the code. In the estimation matrices, a `z_` prefix denotes z-score standardization. This glossary lists the identifiers used in the main matrices without that prefix.

| Identifier | English description |
|---|---|
| `pop` | Population |
| `area` | Municipal land area |
| `densidade` | Population density |
| `pibpc` | GDP per capita |
| `desp.saude` | Municipal spending on health and sanitation |
| `desp.educ` | Municipal spending on education and culture |
| `ideb.fund` | Basic Education Development Index (IDEB), elementary and lower-secondary education; mean of `ideb5` and `ideb9` |
| `ideb3` | IDEB, upper-secondary education |
| `tax.agua` | Urban water-supply coverage |
| `energia100k` | Electricity consumption per 100,000 inhabitants |
| `tax.hom` | Homicide rate |
| `eleitas.fem` | Proportion of elected women |
| `bolsaf` | Population receiving Bolsa Família benefits |
| `sus1k` | SUS-affiliated health facilities per 1,000 inhabitants |
| `leitos1k` | SUS-affiliated hospital beds per 1,000 inhabitants |
| `prof1k` | SUS-affiliated health professionals per 1,000 inhabitants |
| `metrop` | Membership in the Fortaleza Metropolitan Region (dummy) |
| `semiarido` | Membership in the semi-arid region (dummy) |
| `pop.rural` | Rural population share |
| `aliadogov` | Mayor's party alignment with the federal government (dummy) |
| `votos` | Elected mayor's share of valid votes |
| `ivs.infra` | Urban Infrastructure component of the Social Vulnerability Index (IVS) |
| `emprego` | Share of the population in formal employment |
| `ex.pobr` | Share of the population in extreme poverty |
| `imuni` | Immunization coverage among children under one year of age |
| `idosos` | Share of the population aged 60 or older |
| `int.circ` | Hospital admissions for circulatory diseases |
| `int.resp` | Hospital admissions for respiratory diseases |
| `int.diab` | Hospital admissions for diabetes mellitus |
| `int.asma` | Hospital admissions for asthma |

SUS is Brazil's Unified Health System. The `metrop` reference category includes every municipality outside the Fortaleza Metropolitan Region, including municipalities in other metropolitan regions. Detailed source units and periods are recorded in the workbook dictionary; the estimation script defines any rescaling.

## Spatial weights

`pesos_queen.rds` and `pesos_rook.rds` are `listw` objects read by `spdep`. The estimation script converts each object to a matrix and checks its dimensions, zero diagonal, and row sums of one. Preserve municipal order: do not sort a CSV without also reordering the weights.

## Joint-comparison results

| Field | Meaning |
|---|---|
| `desfecho`, `vizinhos` | Dependent variable and contiguity rule |
| `n_mc`, `burn` | Number of iterations and burn-in |
| `p_esp` | Estimated probabilities of the three specifications |
| `p_esp_a`, `p_esp_b` | Separate estimates for each chain |
| `tabela` | Covariate PIPs conditional on each specification |
| `a`, `b` | Summaries for each chain |
| `fora_tab` | Number of evaluations outside the interpolation table's range |

These files contain estimation summaries, not a complete history of visited states. An estimated probability of zero due to no visits should not be read as proof that a specification is impossible. `MQO` denotes OLS in the code and saved results.

## Auxiliary results

`espacial/saida/` also contains outputs from separate estimations by specification, sensitivity analyses, and simulations. Files under `grade`, `priori`, `priori_sem`, `refino_bf`, and `tempo_exposicao` feed post-processing and figure scripts. They do not represent the same posterior distribution as the joint three-way chain. Reference outputs are preserved; new runs made through the documented commands are written to `execucoes_locais/`.
