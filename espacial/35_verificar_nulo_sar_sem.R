# O modelo nulo tem a mesma log-marginal no SAR e no SEM: com k = 0 as duas
# especificacoes filtram so o y, e sao a mesma coisa. As implementacoes sao
# independentes uma da outra, entao a igualdade e uma verificacao cruzada.
source(if (file.exists("../caminhos.R")) "../caminhos.R" else "caminhos.R")

AQUI <- ESPACIAL
Sys.setenv(N_GRADE_LAMBDA = "201")   # igualar a grade do SAR

src <- readLines(file.path(AQUI, "31_sem_hyperg.R"))
corte <- grep("^# Verificacao em toda a grade", src)
eval(parse(text = paste(src[1:(corte - 1)], collapse = "\n")))
v_sem <- lm_sem(integer(0))

# SAR, reimplementado aqui a partir das mesmas pecas.
yc <- y - mean(y)
Wyc <- as.numeric(W %*% y); Wyc <- Wyc - mean(Wyc)
ta <- sum(yc^2); tbb <- sum(yc * Wyc); tcc <- sum(Wyc^2)
rho <- lambda
tss <- ta - 2 * rho * tbb + rho^2 * tcc
lf <- ld_mais_prior + (1 - n) / 2 * log(tss)
m <- max(lf)
v_sar <- m + log(h / 3 * sum(pesos * exp(lf - m)))

cat(sprintf("  log-marginal do modelo nulo | SEM %.10f | SAR %.10f | dif %.3e\n",
            v_sem, v_sar, v_sem - v_sar))
