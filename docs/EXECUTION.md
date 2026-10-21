# Execution

English | [Português](EXECUCAO.md) | [README](../README.md)

Run the commands below from the repository root. `reproduzir.R` also accepts invocation by absolute path and resolves directories from its own location. No machine-specific path edits are required.

Command names remain in Portuguese. Each command is explained below.

## Installation

The checked environment uses R 4.6.1 and Python 3.10.7. R package versions are recorded in [ambiente_testado.csv](ambiente_testado.csv). The installer obtains missing packages from the configured R repository, defaulting to CRAN. It does not pin transitive dependencies to the versions in that record.

```sh
Rscript instalar.R
python -m pip install -r requirements.txt
```

`sf` and `spdep` require geospatial libraries. On Windows, use CRAN binaries where available. On Ubuntu/Debian, install the system libraries before installing R packages:

```sh
sudo apt-get update
sudo apt-get install -y libudunits2-dev libgdal-dev libgeos-dev libproj-dev
```

To record the local package environment:

```sh
Rscript reproduzir.R ambiente
```

## Checks and data preparation

```sh
python scripts/auditar_codigo.py
Rscript reproduzir.R codigo
Rscript reproduzir.R verificar
Rscript reproduzir.R preparar
Rscript reproduzir.R geometria
```

- `auditar_codigo.py` checks Python syntax, UTF-8 encoding, and absolute personal paths, and inventories file reads in R scripts. This is static inspection, not execution of every script.
- `codigo` parses every R file to check syntax.
- `verificar` checks the main matrices, dummy variables, weights, and summaries of the six saved estimation runs.
- `preparar` reads the `raw` sheet in `dadosv6.xlsx`, rebuilds 22 matrices across seven groups, and compares them with the references at a tolerance of 1e-12.
- `geometria` aligns municipal boundaries by IBGE code, rebuilds Queen and Rook weights, and compares them with the reference weights.

Rebuilt matrices are written to `execucoes_locais/preparacao/`; rebuilt weights to `execucoes_locais/geometria/`. Reference inputs remain in `estimacoes/resultados/`.

## Auxiliary objects and interpolation

```sh
Rscript reproduzir.R auxiliares
Rscript reproduzir.R interpolacao
```

`auxiliares` runs scripts 03, 04, and 11 in a separate directory, `execucoes_locais/auxiliares/`, and compares the three generated objects with the included references. It refuses to reuse an existing directory. This checks the construction of `base_sar.rds`, `grade_sar.rds`, and `hyperg_aux.rds`.

`interpolacao` uses these three reference files to build `espacial/saida/tabela_f21_estendida.rds`, with 300,001 grid points, 31 columns, and an upper bound of 0.9. The approximately 66 MB file is generated locally and ignored by Git. Construction took about six minutes in the test environment. The command refuses to overwrite an existing table.

The regenerated table matched the original with a maximum numerical difference of zero. It does not need to be regenerated before each run.

## Numerical validation and smoke test

Once the interpolation table is available:

```sh
Rscript reproduzir.R validar
Rscript reproduzir.R teste
```

`validar` runs scripts 72 and 73. The first tests marginal likelihoods over 1,024 models. The second compares the sampler with enumeration over 3,072 states and includes a negative control with an intentionally incorrect transition. With the fixed test seeds, the correct chains must meet a tolerance of 0.03 for the maximum PIP error and specification-probability error. This tolerance checks the numerical test; it is not a convergence criterion for full estimation runs.

`teste` runs all six outcome/contiguity combinations, using two chains of 30 iterations each and discarding the first 10. These results only test inputs, computation, and output writing. They must not be used for inference.

## Full estimation

```sh
Rscript reproduzir.R estimar
```

This runs cases, deaths, and case fatality under Queen and Rook weights. Each combination uses two chains of 1,000,000 iterations, including a burn-in of 200,000, with seeds defined in script 74. Results are written to a new directory under `execucoes_locais/estimar_*/`. Full estimation is time-consuming and is not part of the quick checks.

To rerun non-spatial BMA and its sensitivity analyses:

```sh
Rscript reproduzir.R bma
```

Outputs go to `execucoes_locais/bma_*/`. This command runs the full settings in `estimacoes/covid_bma.R`, including the expanded model, alternative outcome scales, exclusion of Fortaleza, population cutoffs, and time since the first recorded case.

To run a single spatial combination, change to `espacial/` and set the environment variables `DESFECHO`, `VIZINHOS`, `MC3_SORTEIOS`, `MC3_BURN`, and `SAIDA_TRES` before invoking `Rscript 74_probabilidade_tres_vias.R`. The script checks for an existing output file before starting and refuses to overwrite it. Outcome identifiers are `casos`, `obitos`, and `letalidade`; contiguity identifiers are `queen` and `rook`.

## Tables and figures

```sh
Rscript reproduzir.R exportar
Rscript reproduzir.R dados-figuras
python figuras.py
```

`exportar` reads the six joint-comparison results and exports specification probabilities and conditional PIPs as CSV files to `execucoes_locais/resumos/`.

`dados-figuras` rebuilds eight CSV/GeoJSON input files from the included estimation and simulation outputs. `figuras.py` uses these tables to generate seven figures as PDF and PNG files in `execucoes_locais/figuras/en/`, with English labels and decimal points. The source tables retain their existing identifiers and labels.

Set `FIGURE_LANGUAGE=pt` to generate Portuguese labels and decimal commas in `execucoes_locais/figuras/pt/`, without replacing the English figures. In PowerShell:

```powershell
$env:FIGURE_LANGUAGE = "pt"
python figuras.py
$env:FIGURE_LANGUAGE = "en"
```

On Linux/macOS:

```sh
FIGURE_LANGUAGE=pt python figuras.py
```

In figure 4, "4 active" and "8 active" refer to the number of covariates with nonzero coefficients in the data-generating process. Figure 3 uses "IDEB (grades 5/9)" for the mean of the grade 5 and grade 9 indicators. OLS labels correspond to the `MQO` identifiers in the source tables. SUS, IDEB, IVS, and Bolsa Família retain their Brazilian names.

Figures 3 and 7 compare PIPs and priors from estimations run separately for each specification. They are not calculated from the joint-chain visit frequencies in script 74. Figure 6 concerns the OLS/SAR comparison. Technical sources are listed in script 50 and in the CSV that consolidates the prior grid.

## Coverage

Data preparation, geometry and auxiliary-object reconstruction, interpolation, the smoke test, exports, and generation of the seven figures were executed in copies containing only tracked files. These copies used packages already installed in the test environment; this does not constitute a fresh installation test of every package.

Full chains were not rerun during this check. Previously estimated results are included for inspection and comparison. Other numbered scripts record validations, sensitivity analyses, and experiments; they are not intended to run in numerical order. `scripts/auditar_codigo.py` generates the file-read inventory.

Reproduction starts from the included consolidated dataset and municipal boundaries. Collection of all original administrative data and manuscript compilation are outside this workflow. See [SOURCES.md](SOURCES.md) and [VERIFICATION.md](VERIFICATION.md).

## Supplementary scripts

The Python generators named `gerar_figuras_*.py` provide additional plot designs and require packages from `requirements-extra.txt`. R plotting scripts 34, 41, and 42 also use `ggplot2`, `ggrepel`, `patchwork`, and `tidyr`. These generators are not called by `figuras.py` and are outside the seven-figure check described above.

The language selector applies to the seven main figures, not to these supplementary generators. Some internal comments and messages also remain in Portuguese.

Scripts such as `38_verificar_artigo.R` and appendix generators also refer to article text or manuscript LaTeX files. Those documents are not included in this code package, so these scripts are not standalone reproduction entry points. The documented estimation and figure workflow does not depend on them.
