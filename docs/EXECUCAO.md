# Execução

[English](EXECUTION.md) | Português | [README](../README.pt-BR.md)

Execute os comandos abaixo a partir da raiz do repositório. O arquivo `reproduzir.R` também aceita chamada por caminho absoluto e resolve os diretórios a partir da própria localização. Não é necessário editar caminhos para usar outro computador.

## Instalação

O ambiente conferido usa R 4.6.1 e Python 3.10.7. As versões dos pacotes R estão em [ambiente_testado.csv](ambiente_testado.csv). O instalador abaixo busca os pacotes ausentes no repositório R configurado, com CRAN como padrão; não congela as dependências transitivas nas versões desse registro.

```sh
Rscript instalar.R
python -m pip install -r requirements.txt
```

`sf` e `spdep` dependem de bibliotecas geoespaciais. No Windows, prefira os binários do CRAN. Instalações Linux também exigem as bibliotecas de sistema correspondentes.

Em Ubuntu/Debian, instale também as bibliotecas de sistema antes dos pacotes R:

```sh
sudo apt-get update
sudo apt-get install -y libudunits2-dev libgdal-dev libgeos-dev libproj-dev
```

O registro do ambiente é gerado com:

```sh
Rscript reproduzir.R ambiente
```

## Verificações e preparação

```sh
python scripts/auditar_codigo.py
Rscript reproduzir.R codigo
Rscript reproduzir.R verificar
Rscript reproduzir.R preparar
Rscript reproduzir.R geometria
```

- `auditar_codigo.py` verifica a sintaxe Python, UTF-8 e caminhos pessoais absolutos e inventaria as leituras nos scripts R. É uma inspeção estática, não uma execução de todos os scripts.
- `codigo` analisa a sintaxe de todos os arquivos R.
- `verificar` confere as matrizes principais, as dummies, os pesos e os resumos das seis cadeias salvas.
- `preparar` lê a aba `raw` de `dadosv6.xlsx`, reconstrói 22 matrizes em sete grupos e compara todas com as referências, com tolerância de 1e-12.
- `geometria` alinha a malha municipal por código IBGE, reconstrói Queen e Rook e compara os pesos com os arquivos de referência.

As matrizes reconstruídas são gravadas em `execucoes_locais/preparacao/`; os pesos reconstruídos, em `execucoes_locais/geometria/`. As entradas de referência permanecem em `estimacoes/resultados/`.

## Auxiliares e interpolação

```sh
Rscript reproduzir.R auxiliares
Rscript reproduzir.R interpolacao
```

`auxiliares` executa os scripts 03, 04 e 11 em uma pasta separada, `execucoes_locais/auxiliares/`, e compara os três objetos gerados com os auxiliares incluídos. Recusa reutilizar uma pasta existente. Essa etapa permite verificar a origem de `base_sar.rds`, `grade_sar.rds` e `hyperg_aux.rds`.

`interpolacao` usa esses três arquivos de referência para construir `espacial/saida/tabela_f21_estendida.rds`, com 300.001 pontos, 31 colunas e limite superior 0,9. O arquivo, de aproximadamente 66 MB, é gerado localmente e ignorado pelo Git. Sua construção levou cerca de seis minutos no teste. O comando recusa sobrescrever uma tabela existente.

A tabela regenerada foi comparada à original e apresentou diferença numérica máxima igual a zero. Não é necessário regenerá-la antes de cada execução.

## Validação e teste de execução

Com a tabela disponível:

```sh
Rscript reproduzir.R validar
Rscript reproduzir.R teste
```

`validar` executa os scripts 72 e 73. O primeiro testa as marginais em 1.024 modelos. O segundo compara o amostrador com enumeração em 3.072 estados e inclui um controle negativo com transição deliberadamente incorreta. As cadeias corretas devem respeitar a margem de 0,03 para o erro máximo de PIP e de probabilidade de especificação, sob as sementes fixas desse teste. Essa margem verifica o teste numérico, não a convergência das estimações completas.

`teste` executa as seis combinações de variável dependente e vizinhança, com duas cadeias de 30 iterações e descarte de 10. Esses resultados servem exclusivamente para verificar entradas, cálculo e gravação. Não devem ser usados para inferência.

## Estimações completas

```sh
Rscript reproduzir.R estimar
```

Esse comando executa casos, óbitos e letalidade sob Queen e Rook. Cada combinação usa duas cadeias de 1.000.000 de iterações, burn-in de 200.000 e as sementes definidas no script 74. Os resultados são gravados em uma pasta nova em `execucoes_locais/estimar_*/`. A execução completa é demorada e não integra os testes rápidos.

Para refazer também o BMA não espacial e suas sensibilidades:

```sh
Rscript reproduzir.R bma
```

As saídas ficam em `execucoes_locais/bma_*/`. Esse comando executa as configurações completas de `estimacoes/covid_bma.R`, incluindo o modelo ampliado, alterações de escala, exclusão de Fortaleza, cortes populacionais e tempo desde o primeiro registro.

Para executar somente uma combinação espacial, entre em `espacial/` e use as variáveis de ambiente `DESFECHO`, `VIZINHOS`, `MC3_SORTEIOS`, `MC3_BURN` e `SAIDA_TRES` antes de chamar `Rscript 74_probabilidade_tres_vias.R`. O script verifica a existência da saída antes de iniciar e recusa sobrescrevê-la.

## Tabelas e figuras

```sh
Rscript reproduzir.R exportar
Rscript reproduzir.R dados-figuras
python figuras.py
```

`exportar` lê os seis resultados da comparação conjunta e grava probabilidades de especificação e PIPs condicionais em CSV em `execucoes_locais/resumos/`.

`dados-figuras` reconstrói oito arquivos CSV/GeoJSON a partir das saídas de estimação e simulação incluídas. `figuras.py` usa essas tabelas para gerar sete figuras em PDF e PNG, em inglês por padrão, em `execucoes_locais/figuras/en/`.

Para gerar a versão em português sem substituir a inglesa, defina `FIGURE_LANGUAGE=pt`. No PowerShell:

```powershell
$env:FIGURE_LANGUAGE = "pt"
python figuras.py
$env:FIGURE_LANGUAGE = "en"
```

As versões em português ficam em `execucoes_locais/figuras/pt/`. O seletor se aplica às sete figuras principais, não aos geradores complementares `gerar_figuras_*.py`. Na Figura 4, os rótulos ingleses "4 active" e "8 active" indicam o número de covariáveis com coeficientes diferentes de zero no processo gerador dos dados.

As figuras 3 e 7 comparam PIPs e prioris de estimações feitas separadamente por especificação. Elas não são calculadas a partir das frequências da cadeia conjunta do script 74. A figura 6 refere-se à comparação MQO/SAR. As fontes técnicas estão no script 50 e no CSV de consolidação da grade de prioris.

## Cobertura

Os comandos de preparação, reconstrução da geometria e dos auxiliares, geração da interpolação, teste curto, exportação e geração das sete figuras foram executados em cópias extraídas apenas dos arquivos versionados. A cópia usa os pacotes instalados no ambiente de teste; isso não equivale a testar uma instalação nova de todos os pacotes.

As cadeias completas não foram reestimadas durante essa verificação. Os resultados anteriores estão incluídos para consulta e comparação. Os demais scripts numerados registram validações, sensibilidades e experimentos; sua presença não significa que devam ser executados em ordem numérica. O inventário de leituras é gerado por `scripts/auditar_codigo.py`.

A reprodução começa na base consolidada e na malha incluídas. A coleta de todos os dados administrativos de origem e a compilação do manuscrito não fazem parte deste fluxo. Consulte [FONTES.md](FONTES.md).

## Scripts complementares

Os geradores Python `gerar_figuras_*.py` incluem outros desenhos gráficos e exigem os pacotes de `requirements-extra.txt`. As versões gráficas em R, nos scripts 34, 41 e 42, usam também `ggplot2`, `ggrepel`, `patchwork` e `tidyr`. Esses geradores não integram o comando `figuras.py` nem a verificação das sete figuras descrita acima.

Scripts como `38_verificar_artigo.R` e os geradores de apêndices também se referem ao texto do artigo ou a arquivos LaTeX do manuscrito. Como esses documentos não fazem parte do pacote de código, esses comandos não são entradas autônomas de reprodução. O fluxo documentado de estimação e figuras independe deles.
