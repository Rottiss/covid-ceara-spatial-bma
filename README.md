# COVID-19 in Ceará: spatial Bayesian model averaging

English | [Português](README.pt-BR.md)

This project studies associations between municipal characteristics and recorded COVID-19 cases, deaths, and case fatality in Ceará, Brazil. Bayesian model averaging (BMA) accounts for uncertainty about covariate inclusion and spatial specification by comparing ordinary least squares (OLS), spatial autoregressive (SAR), and spatial error (SEM) models.

## Analysis

The main estimation matrices contain 184 municipalities and 30 candidate covariates. Cases, deaths, and case fatality are analyzed separately, using Queen and Rook contiguity weights.

| Specification | Spatial component |
|---|---|
| OLS (`MQO` in the code and saved results) | No spatial lag |
| SAR | Spatial lag of the dependent variable |
| SEM | Spatial dependence in the errors |

The MC³ sampler explores models that differ in their covariates and specification. Results describe municipal associations, not causal effects or individual-level relationships. Posterior inclusion probabilities (PIPs) depend on the model space and priors.

## Available results

The six files in [espacial/saida/tres_vias](espacial/saida/tres_vias) contain estimated specification probabilities, conditional PIPs, and summaries of two chains for each outcome and contiguity rule. Each chain uses 1,000,000 iterations, including a burn-in of 200,000.

You can inspect these results without rerunning the sampler. From the repository root in R:

```r
result <- readRDS("espacial/saida/tres_vias/casos_queen.rds")
result$p_esp  # estimated probabilities of OLS (MQO), SAR, and SEM
result$tabela # PIPs conditional on each specification
```

No visits to a specification does not establish that its posterior probability is exactly zero.

## Code and data

| Path | Contents |
|---|---|
| [espacial/](espacial) | Estimation, numerical tests, sensitivity analyses, and post-processing |
| [espacial/saida/tres_vias/](espacial/saida/tres_vias) | Joint OLS/SAR/SEM comparison results |
| [Estimation matrices](estimacoes/resultados/estimacoes/matrizes_estimacao/principal) | Main analysis inputs |
| [Spatial weights](estimacoes/resultados/diagnostico_espacial) | Queen and Rook weights and original environment records |
| [Execution guide](docs/EXECUTION.md) | Dependencies, commands, and reproduction limits |
| [Data description](docs/DATA.md) | Input and output structure, with an English variable glossary |
| [Sources and transformations](docs/SOURCES.md) | Data provenance and variable construction |
| [Reproduction checks](docs/VERIFICATION.md) | Tests performed and their scope |

`espacial/74_probabilidade_tres_vias.R` is the entry point for the joint comparison. Scripts 72 and 73 test the scale of marginal likelihoods and transitions between specifications. File numbers do not define an execution order.

## Run

```sh
Rscript instalar.R
Rscript reproduzir.R verificar
Rscript reproduzir.R preparar
Rscript reproduzir.R geometria
```

These commands check saved results and rebuild estimation matrices and spatial weights from the included workbook and municipal boundaries. To generate the seven spatial figures:

```sh
python -m pip install -r requirements.txt
Rscript reproduzir.R dados-figuras
python figuras.py
```

Generated files are written to `execucoes_locais/`. The [execution guide](docs/EXECUTION.md) covers interpolation, numerical validation, and full estimation. Preparation and artifact generation were checked in copies containing only tracked files. Full estimation chains were not rerun as part of that check.

Reproduction starts from the consolidated dataset. Collection from the original administrative portals is not automated. See [verification coverage](docs/VERIFICATION.md).

File names, command names, and data fields retain their original Portuguese identifiers to avoid breaking dependencies. For example, `casos` means cases, `obitos` means deaths, and `letalidade` means case fatality. The seven main figures default to English and are saved under `execucoes_locais/figuras/en/`. Set `FIGURE_LANGUAGE=pt` to generate Portuguese versions in a separate folder. Supplementary plotting scripts and some internal comments and messages remain in Portuguese.

## Author

Lucas Siqueira de Castro. Developed from undergraduate research in Economics at the Federal University of Ceará, Brazil.

## License

No license for use or redistribution of the code has been selected yet. Third-party data remain subject to the terms of their respective sources.
