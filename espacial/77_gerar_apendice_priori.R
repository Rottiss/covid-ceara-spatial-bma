# RASCUNHO do apendice com a grade completa de sensibilidade a priori:
# 5 combinacoes x 3 variaveis dependentes x 2 matrizes, sob MQO, SAR e SEM.
#
# Le saida/priori_sem/consolidado_sem.csv (76), que ja alinha tudo por nome.
# GRAVA apenas em rascunho/apendice-h.tex, fora da copia do TCC. Copiado para
# 3-pos-textuais/apendices/apendice-h.tex em 17/09; ao regenerar, copiar de novo.

source(if (file.exists("../caminhos.R")) "../caminhos.R" else "caminhos.R")

S <- file.path(ESPACIAL, "saida")
DEST <- file.path(ESPACIAL, "rascunho/apendice-h.tex")
res <- read.csv(file.path(S, "priori_sem", "consolidado_sem.csv"), stringsAsFactors = FALSE)
stopifnot(nrow(res) == 30, all(res$fora_tab == 0))

COMB <- c(binbeta_hyperg   = "Binomial-Beta $K/2$ e hiper-$g$",
          binbeta_uip      = "Binomial-Beta $K/2$ e $g = n$",
          binbeta5_hyperg  = "Binomial-Beta $E(k) = 5$ e hiper-$g$",
          binomial5_hyperg = "Binomial $5/K$ e hiper-$g$",
          uniforme_hyperg  = "Uniforme sobre $2^K$ e hiper-$g$")
DESF <- c(casos = "Casos", obitos = "Óbitos", letalidade = "Letalidade")
vir <- function(x) formatC(x, format = "f", digits = 3, decimal.mark = ",")

L <- character(0); add <- function(...) L <<- c(L, sprintf(...))
add("\\apendice{Sensibilidade às distribuições a priori sob SAR e SEM nas duas matrizes de vizinhança}")
add("\\label{ap:priori-completa}")
add("")
add("\\begin{table}[htbp]")
add("    \\centering")
add("    \\Caption{\\label{tab:priori-completa} Covariáveis com PIP de pelo menos 0,50 e razão $D$, por combinação de prioris, especificação e matriz}")
add("    \\UFCtab{}{")
add("        \\scriptsize")
add("        \\begin{adjustbox}{max width=0.98\\textwidth}")
add("        \\begin{tabular}{lrrrrrrrrr}")
add("            \\toprule")
add("            & & \\multicolumn{4}{c}{Queen} & \\multicolumn{4}{c}{Rook} \\\\")
add("            \\cmidrule(lr){3-6} \\cmidrule(lr){7-10}")
add("            Combinação de prioris & MQO & SAR & SEM & Razão SAR & Razão SEM & SAR & SEM & Razão SAR & Razão SEM \\\\")
add("            \\midrule \\midrule")
for (d in names(DESF)) {
  add("            \\multicolumn{10}{l}{\\textit{%s}} \\\\", DESF[[d]])
  for (cb in names(COMB)) {
    q <- res[res$desfecho == d & res$vizinhos == "queen" & res$combinacao == cb, ]
    r <- res[res$desfecho == d & res$vizinhos == "rook" & res$combinacao == cb, ]
    stopifnot(nrow(q) == 1, nrow(r) == 1, q$acima_mqo == r$acima_mqo)
    add("            %s & %d & %d & %d & %s & %s & %d & %d & %s & %s \\\\", COMB[[cb]], q$acima_mqo,
        q$acima_sar, q$acima_sem, vir(q$razao_sar), vir(q$razao_sem),
        r$acima_sar, r$acima_sem, vir(r$razao_sar), vir(r$razao_sem))
  }
  if (d != "letalidade") add("            \\midrule")
}
add("            \\bottomrule")
add("        \\end{tabular}")
add("        \\end{adjustbox}")
add("    }{")
add("        \\Fonte{Elaboração própria com base nos resultados da pesquisa.}")
add("        \\Nota{As colunas MQO, SAR e SEM trazem o número de covariáveis, entre as 30, com PIP igual ou superior a 0,50. A coluna MQO é única porque o modelo não utiliza a matriz de pesos. A razão $D$ compara a dispersão das PIP sob a especificação espacial com a observada sob MQO na mesma combinação de prioris, como na Tabela~\\ref{tab:d-espacial}. A primeira linha de cada bloco é a combinação adotada no restante do trabalho. Maior diferença de PIP entre as duas cadeias de cada estimação sob SEM: %s.}",
    formatC(max(res$conv_sem), format = "f", digits = 4, decimal.mark = ","))
add("    }")
add("\\end{table}")

dir.create(dirname(DEST), showWarnings = FALSE, recursive = TRUE)
writeLines(L, DEST, useBytes = FALSE)
cat("gravado:", DEST, "\n")
