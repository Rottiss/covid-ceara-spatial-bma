# COVID-19 no Ceará: Bayesian Model Averaging espacial

[English](README.md) | Português

Este projeto investiga as associações entre características municipais e os registros de casos, óbitos e letalidade por COVID-19 no Ceará. A análise utiliza Bayesian Model Averaging (BMA) para considerar a incerteza sobre a inclusão de covariáveis e a especificação espacial, comparando MQO, SAR e SEM.

## Análise

As matrizes de estimação contêm 184 municípios e 30 covariáveis. A comparação é realizada separadamente para casos, óbitos e letalidade, com matrizes de vizinhança Queen e Rook.

| Especificação | Componente espacial |
|---|---|
| MQO | Sem defasagem espacial |
| SAR | Defasagem espacial da variável dependente |
| SEM | Dependência espacial no erro |

O amostrador MC³ percorre modelos que diferem no conjunto de covariáveis e na especificação. Os resultados descrevem associações entre municípios, não efeitos causais ou relações individuais. As probabilidades de inclusão dependem do conjunto de modelos e das prioris adotadas.

## Resultados disponíveis

Os seis arquivos em [espacial/saida/tres_vias](espacial/saida/tres_vias) contêm as probabilidades estimadas das especificações e as PIPs, além de resumos das duas cadeias de cada combinação de variável dependente e vizinhança. Cada execução registra 1.000.000 de iterações e burn-in de 200.000.

É possível consultar esses resultados sem reestimar os modelos. Em R, a partir da raiz:

```r
resultado <- readRDS("espacial/saida/tres_vias/casos_queen.rds")
resultado$p_esp   # probabilidades estimadas de MQO, SAR e SEM
resultado$tabela # PIPs condicionais a cada especificação
```

Uma frequência zero de visitas a uma especificação não demonstra que sua probabilidade posterior seja exatamente zero.

## Código e dados

| Caminho | Conteúdo |
|---|---|
| [espacial/](espacial) | Estimação, testes numéricos, sensibilidades e pós-processamento |
| [espacial/saida/tres_vias/](espacial/saida/tres_vias) | Resultados da comparação MQO/SAR/SEM |
| [Matrizes de estimação](estimacoes/resultados/estimacoes/matrizes_estimacao/principal) | Entradas da análise principal |
| [Pesos espaciais](estimacoes/resultados/diagnostico_espacial) | Pesos Queen e Rook e registro do ambiente de origem |
| [Guia de execução](docs/EXECUCAO.md) | Dependências, comandos e limites de reprodução |
| [Descrição dos dados](docs/DADOS.md) | Estrutura das entradas e saídas |
| [Fontes e transformações](docs/FONTES.md) | Origem da base e construção das variáveis |
| [Verificação de reprodução](docs/VERIFICACAO.md) | Testes executados e seu alcance |

O ponto de entrada da comparação conjunta é `74_probabilidade_tres_vias.R`. Os scripts `72_validar_escala_tres.R` e `73_validar_troca_tres_vias.R` testam a escala das verossimilhanças marginais e a transição entre especificações. A numeração dos arquivos não define uma sequência de execução.

## Executar

```sh
Rscript instalar.R
Rscript reproduzir.R verificar
Rscript reproduzir.R preparar
Rscript reproduzir.R geometria
```

Esses comandos verificam os resultados salvos e reconstroem as matrizes e os pesos a partir da base consolidada e da malha incluídas. Para gerar as sete figuras espaciais:

```sh
python -m pip install -r requirements.txt
Rscript reproduzir.R dados-figuras
python figuras.py
```

Os arquivos gerados ficam em `execucoes_locais/`. O [guia de execução](docs/EXECUCAO.md) também descreve a geração da tabela de interpolação, os testes numéricos e as estimações completas. As etapas de preparação e geração de artefatos foram conferidas em cópias contendo apenas os arquivos versionados; as cadeias completas não foram reestimadas nessa verificação.

As sete figuras principais são geradas em inglês por padrão, em `execucoes_locais/figuras/en/`. Defina `FIGURE_LANGUAGE=pt` para gerar a versão em português em uma pasta separada. Os identificadores dos dados e os comandos não foram traduzidos.

O fluxo começa na base consolidada. A coleta em cada portal administrativo de origem não é automatizada. O escopo e os limites dos testes estão em [VERIFICACAO.md](docs/VERIFICACAO.md).

## Autoria

Lucas Siqueira de Castro. Projeto desenvolvido a partir de pesquisa de graduação em Ciências Econômicas na Universidade Federal do Ceará.

## Licença

Ainda não há licença de uso e redistribuição definida para o código. Os dados de terceiros permanecem sujeitos aos termos de suas respectivas fontes.
