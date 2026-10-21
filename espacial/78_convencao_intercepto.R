# Diagnostico: efeito da convencao do intercepto nas marginais espaciais, no
# problema reduzido do 73 (10 covariaveis, enumeracao exata). Nao toca em
# nenhuma saida existente; grava saida/convencao_intercepto.csv.
#
#   filtrada : priori plana no intercepto da equacao filtrada (SAR e SEM), a implementada
#   literal  : priori plana no alpha como escrito no TCC; SAR sem fator, SEM com 1/|1-lambda|
#   nivel    : priori plana no nivel medio mu nos dois; SAR com 1/|1-rho|, SEM com 1/|1-lambda|
source(if (file.exists("../caminhos.R")) "../caminhos.R" else "caminhos.R")

B <- ESPACIAL
s0 <- readLines(file.path(B, "73_validar_troca_tres_vias.R"), encoding = "UTF-8")
s0 <- s0[1:(grep("^# -+ amostrador", s0)[1] - 1)]
COL <- c(casos = "z_log_casos100k", obitos = "z_log_obitos100k", letalidade = "z_logit_letal")

saida <- list()

for (d in names(COL)) {
  s <- sub('DESFECHO <- "casos"; COLUNA <- "z_log_casos100k"',
           sprintf('DESFECHO <- "%s"; COLUNA <- "%s"', d, COL[[d]]), s0, fixed = TRUE)
  stopifnot(sum(grepl(sprintf('DESFECHO <- "%s"', d), s)) == 1)
  en <- new.env()
  invisible(capture.output(eval(parse(text = s), en)))
  base <- en$M
  # versoes com fator de jacobiano dentro da integral
  patch <- function(fname, grade) {
    fn <- deparse(en[[fname]])
    alvo <- sprintf("lf <- %s$ldp +", grade)
    stopifnot(sum(grepl(alvo, fn, fixed = TRUE)) == 1)
    fn <- sub(alvo, sprintf("lf <- %s$ldp - log1p(-%s$v) +", grade, grade), fn, fixed = TRUE)
    f <- eval(parse(text = fn)); environment(f) <- en; f
  }
  sar_j <- patch("lm_sar", "gr"); sem_j <- patch("lm_sem", "gl")
  Msar <- vapply(en$modelos, function(idx) sar_j(en$estat(idx)), numeric(1))
  Msem <- vapply(en$modelos, function(idx) sem_j(idx), numeric(1))
  conv <- list(filtrada = base,
               literal = cbind(MQO = base[, "MQO"], SAR = base[, "SAR"], SEM = Msem),
               nivel = cbind(MQO = base[, "MQO"], SAR = Msar, SEM = Msem))
  cat("\n==", d, "(10 covariaveis) ==\n")
  ref <- NULL
  for (nm in names(conv)) {
    en$POST <- conv[[nm]] + en$lpm[en$tam + 1]
    tot <- en$exato(en$ESP)
    # PIP condicional a cada especificacao
    pc <- sapply(c("SAR", "SEM"), function(e) en$exato(e)$pip)
    lbf <- function(e) { p <- conv[[nm]][, e] + en$lpm[en$tam + 1]; q <- conv[[nm]][, "MQO"] + en$lpm[en$tam + 1]
      (max(p) + log(sum(exp(p - max(p))))) - (max(q) + log(sum(exp(q - max(q))))) }
    if (is.null(ref)) ref <- pc
    saida[[length(saida) + 1]] <- data.frame(
      desfecho = d, convencao = nm,
      p_mqo = tot$esp[1], p_sar = tot$esp[2], p_sem = tot$esp[3],
      lbf_sar = lbf("SAR"), lbf_sem = lbf("SEM"),
      dpip_sar = max(abs(pc[, "SAR"] - ref[, "SAR"])),
      dpip_sem = max(abs(pc[, "SEM"] - ref[, "SEM"])),
      stringsAsFactors = FALSE)
    cat(sprintf("  %-9s P(MQO,SAR,SEM) = %.4f %.4f %.4f | log FB SAR/MQO %.2f | SEM/MQO %.2f | dPIP SAR %.4f SEM %.4f\n",
        nm, tot$esp[1], tot$esp[2], tot$esp[3], lbf("SAR"), lbf("SEM"),
        max(abs(pc[, "SAR"] - ref[, "SAR"])), max(abs(pc[, "SEM"] - ref[, "SEM"]))))
  }
}

saida <- do.call(rbind, saida)
arq <- file.path(B, "saida", "convencao_intercepto.csv")
write.csv(saida, arq, row.names = FALSE)
cat("\ngravado:", arq, "\n")
