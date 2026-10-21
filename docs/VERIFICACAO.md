# Verificação de reprodução

[English](VERIFICATION.md) | Português | [README](../README.pt-BR.md)

Verificações realizadas em 20/10/2026, em cópias extraídas do índice do Git. Os arquivos locais ignorados pelo Git não foram incorporados à extração.

| Etapa | Verificação |
|---|---|
| Preparação | 22 matrizes reconstruídas a partir da planilha; todas concordaram com as referências com tolerância de 1e-12 |
| Geometria | Pesos Queen e Rook reconstruídos a partir da malha; diferença máxima igual a zero |
| Auxiliares | Os três objetos gerados pelos scripts 03, 04 e 11 concordaram com os arquivos de referência com tolerância de 1e-12 |
| Interpolação | Tabela gerada pelo script 15; diferença máxima de seus valores em relação ao auxiliar original igual a zero |
| Escala das marginais | Erro máximo de 5,684e-14 no problema com 1.024 modelos |
| Transição entre especificações | Erro máximo de PIP de 0,00419 e erro de probabilidade de especificação de 0,00036 no teste de três vias contra enumeração |
| Teste curto | Seis combinações executadas e gravadas, com duas cadeias de 30 iterações por combinação |
| Figuras | Oito tabelas de entrada reconstruídas; sete figuras geradas em PNG e PDF e inspecionadas visualmente |
| Idiomas das figuras | Sete figuras principais geradas em inglês e português; os arrays dos dados desenhados concordaram entre os idiomas |

O teste de transição do script 73 compara a cadeia com enumeração em um problema reduzido. Seus resultados devem ser lidos separadamente dos testes de execução curta e das estimações completas.

A validação de caminhos foi realizada a partir de diretórios diferentes do diretório de trabalho. A configuração de locale UTF-8 é feita pelo comando principal antes de ler os scripts. Os arquivos Python e R também passam por verificação de sintaxe.

O registro das versões dos pacotes está em [ambiente_testado.csv](ambiente_testado.csv); Python 3.10.7 foi usado com as versões de [requirements.txt](../requirements.txt).

A [execução independente no GitHub Actions](https://github.com/Rottiss/covid-ceara-spatial-bma/actions/runs/39516352017) foi concluída com sucesso em uma instalação nova de R no Ubuntu 24.04. Ela verificou a sintaxe R, os contratos das entradas e resultados, a reconstrução das 22 matrizes e a reconstrução dos pesos Queen e Rook. O workflow fixa essa versão do Ubuntu.

Essa verificação não reexecutou as cadeias completas nem reconstruiu a coleta das fontes administrativas.

A tradução também preservou os cálculos dos arquivos R, conferidos pela estrutura sintática após excluir mensagens de exibição e comentários. A preparação dos dados, a checagem dos resultados de referência, a reconstrução dos pesos e a exportação dos resumos foram executadas novamente após essa alteração.
