# Sources and transformations

English | [Português](FONTES.md) | [README](../README.md)

The consolidated dataset is in `estimacoes/dadosv6.xlsx`. `covid_bma.R` reads only the `raw` sheet. Other sheets support inspection and are not estimation inputs.

The institutions and periods below follow the documentation of the dataset used in the research. Collection from the original portals is not automated here. Computational reproduction starts from the consolidated workbook; no command rebuilds it from all original administrative files.

| Variables | Reference period | Recorded source |
|---|---|---|
| Cases and deaths | Cumulative through July 31, 2020 | State Health Departments, consolidated by Brasil.IO |
| Population; share aged 60 or older | 2019 | DATASUS |
| Area; population density | 2019 | IBGE |
| GDP per capita | 2018 | IBGE |
| Rural population share | 2010 | IBGE |
| Formal employment | 2019 | IPECE |
| Extreme poverty | 2010 | IPECE |
| Bolsa Família beneficiaries | 2020 | Ministry of Citizenship, SENARC |
| IVS Urban Infrastructure | 2010 | IPEA, Atlas da Vulnerabilidade Social |
| Fortaleza Metropolitan Region | 2018 boundaries | IPECE |
| Semi-arid region | 2020 | IPECE |
| Party alignment | 2019 classification | TRE-CE and TSE |
| Valid votes and proportion of elected women | 2016 election | TRE-CE and TSE |
| Homicides | 2019 | SIM, DATASUS |
| Health and sanitation spending; education and culture spending | 2018 | IPEAdata, using STN data |
| Elementary, lower-secondary, and upper-secondary IDEB | 2019 | INEP |
| Urban water-supply coverage | 2019 | IPECE |
| Electricity consumption | 2019 | IPECE and Coelce |
| Immunization among children under one year | 2019 | IPECE and SESA |
| SUS-affiliated health facilities | 2019 | SESA |
| SUS-affiliated hospital beds and health professionals | 2019 | IPECE and SESA |
| Admissions for circulatory diseases, respiratory diseases, diabetes, and asthma | 2019 | SIH, DATASUS |
| Federal transfers to address COVID-19 | March to July 2020 | Portal da Transparência, CGU |
| Emergency cash assistance | April to July 2020 | Portal da Transparência, CGU |
| Days since the first municipal record | First record through July 31, 2020 | State Health Departments, consolidated by Brasil.IO |
| Municipal boundaries | 2020 | IBGE boundaries obtained through the geobr package |

Institutional acronyms and source names retain their Brazilian forms so they can be matched to the workbook documentation.

## Estimation transformations

`estimacoes/covid_bma.R` defines and checks the transformations. The main ones are:

- Cases: `log(covtotal / pop * 100000)`.
- Deaths: `log((obito + 0.5) / pop * 100000)`.
- Case fatality: `log((obito + 0.5) / (covtotal - obito + 0.5))`.
- Elementary and lower-secondary IDEB: simple mean of `ideb5` and `ideb9`.
- Electricity consumption: `con.energia / pop * 100000`.
- Standardization: `(x - mean(x)) / sd(x)`, applied to continuous variables within each specification. The three dummy variables remain coded 0/1.

The main covariates are listed in `covariaveis_principais` in the script. The workbook retains fields not used in that model. The dummy `metrop` identifies the Fortaleza Metropolitan Region; all other municipalities form the reference category.

`Rscript reproduzir.R preparar` rebuilds the matrices and compares them with the included references. `Rscript reproduzir.R geometria` aligns the boundaries by IBGE code, rebuilds both contiguity rules, and compares the weights with the originals.
