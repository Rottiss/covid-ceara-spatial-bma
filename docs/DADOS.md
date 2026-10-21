# Entradas e resultados

[English](DATA.md) | Português | [README](../README.pt-BR.md)

## Matrizes de estimação

Os arquivos `casos.csv`, `obitos.csv` e `letalidade.csv`, em `estimacoes/resultados/estimacoes/matrizes_estimacao/principal/`, contêm as entradas da análise principal. Cada linha corresponde a um município.

| Campo | Conteúdo |
|---|---|
| `codigo_ibge` | Identificador municipal; deve ser lido como texto |
| `municipio` | Nome do município |
| `z_log_casos100k`, `z_log_obitos100k` ou `z_logit_letal` | Variável dependente transformada e padronizada, conforme o arquivo |
| Demais colunas | 30 covariáveis usadas na estimação |

As variáveis contínuas entram padronizadas por z-score. As dummies `metrop`, `semiarido` e `aliadogov` permanecem na escala 0/1. Essas matrizes não são dados brutos: suas colunas já incorporam as transformações anteriores à estimação.

A base consolidada está em `estimacoes/dadosv6.xlsx`, e as fórmulas estão em `estimacoes/covid_bma.R`. As fontes e transformações principais são descritas em [FONTES.md](FONTES.md). O comando `preparar` reconstrói e compara as 22 matrizes incluídas, que cobrem o modelo principal e as sensibilidades.

## Vizinhança

`pesos_queen.rds` e `pesos_rook.rds` são objetos `listw`, lidos pelo pacote `spdep`. O script de estimação converte cada objeto em matriz e verifica a dimensão, a diagonal nula e a soma unitária das linhas. A ordem dos municípios deve ser preservada: não ordene o CSV sem reordenar também os pesos.

## Resultados da comparação conjunta

| Campo | Conteúdo |
|---|---|
| `desfecho`, `vizinhos` | Variável dependente e matriz de vizinhança |
| `n_mc`, `burn` | Número de iterações e burn-in |
| `p_esp` | Probabilidades estimadas das três especificações |
| `p_esp_a`, `p_esp_b` | Estimativas separadas por cadeia |
| `tabela` | PIPs das covariáveis condicionais a cada especificação |
| `a`, `b` | Resumos de cada cadeia |
| `fora_tab` | Contagem de avaliações fora do intervalo da tabela de interpolação |

Os arquivos são resumos de estimação, não um histórico completo de todos os estados visitados. Uma probabilidade estimada igual a zero por ausência de visitas não deve ser interpretada como impossibilidade da especificação.

## Resultados auxiliares

`espacial/saida/` contém também saídas de estimações separadas por especificação, sensibilidades e simulações. Os arquivos nas pastas `grade`, `priori`, `priori_sem`, `refino_bf` e `tempo_exposicao` alimentam os scripts de pós-processamento e figuras. Eles não representam a mesma distribuição posterior da cadeia conjunta de três vias. As saídas de referência são preservadas; novas execuções feitas pelos comandos documentados ficam em `execucoes_locais/`.
