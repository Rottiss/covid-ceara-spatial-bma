# Explore associations between covariate Moran's I and OLS/SAR PIP differences.
# The regressions are descriptive and do not identify a causal mechanism.

source(if (file.exists("../caminhos.R")) "../caminhos.R" else "caminhos.R")

suppressPackageStartupMessages(library(spdep))

AQUI <- ESPACIAL
base <- readRDS(file.path(AQUI, "saida", "base_sar.rds"))
des  <- readRDS(file.path(AQUI, "saida", "desenho_simulacao.rds"))
real <- read.csv(file.path(AQUI, "saida", "comparacao_1milhao.csv"),
                 stringsAsFactors = FALSE)

X <- base$X; covariaveis <- base$covariaveis; listw <- base$listw
moran <- des$moran

d <- merge(real, data.frame(covariavel = names(moran), mi = as.numeric(moran),
                            stringsAsFactors = FALSE), by = "covariavel")
d$queda <- d$mqo - d$sar
d$queda_rel <- d$queda / d$mqo

cat("=== Queda da PIP contra I de Moran da covariavel ===\n")
cat("  dados reais, casos com Queen, 1.000.000 de iteracoes\n\n")

fit_abs <- lm(queda ~ mi, data = d)
fit_rel <- lm(queda_rel ~ mi, data = d)
cat(sprintf("  queda absoluta ~ I de Moran: coeficiente %.4f (ep %.4f), p = %.5f, R2 = %.4f\n",
            coef(fit_abs)[2], summary(fit_abs)$coefficients[2, 2],
            summary(fit_abs)$coefficients[2, 4], summary(fit_abs)$r.squared))
cat(sprintf("  queda relativa ~ I de Moran: coeficiente %.4f (ep %.4f), p = %.5f, R2 = %.4f\n",
            coef(fit_rel)[2], summary(fit_rel)$coefficients[2, 2],
            summary(fit_rel)$coefficients[2, 4], summary(fit_rel)$r.squared))
cat(sprintf("\n  correlacao de Spearman entre queda e I de Moran: %.4f\n",
            cor(d$queda, d$mi, method = "spearman")))

# Larger initial OLS PIPs allow larger absolute decreases.
# Controla por isso.
fit_ctl <- lm(queda ~ mi + mqo, data = d)
cat("\n=== Controlando pela PIP no MQO ===\n")
cat("  (covariavel com PIP alta tem mais de onde cair; isso precisa sair)\n\n")
print(round(summary(fit_ctl)$coefficients, 5))
cat(sprintf("\n  R2: %.4f\n", summary(fit_ctl)$r.squared))

cat("\n=== Terços de I de Moran ===\n")
d$tercil <- cut(d$mi, breaks = quantile(d$mi, c(0, 1/3, 2/3, 1)),
                include.lowest = TRUE, labels = c("baixo", "medio", "alto"))
resumo <- aggregate(cbind(mi, mqo, sar, queda, queda_rel) ~ tercil, data = d, FUN = mean)
print(transform(resumo, mi = round(mi, 4), mqo = round(mqo, 4), sar = round(sar, 4),
                queda = round(queda, 4), queda_rel = round(queda_rel, 4)),
      row.names = FALSE)

# A previsao alternativa: talvez a queda dependa so da PIP inicial, e o I de
# Moran nao acrescente nada. Compara os dois modelos.
cat("\n=== Comparacao de modelos para a queda ===\n")
m1 <- lm(queda ~ mqo, data = d)
m2 <- lm(queda ~ mqo + mi, data = d)
cat(sprintf("  so PIP do MQO      : R2 = %.4f\n", summary(m1)$r.squared))
cat(sprintf("  PIP do MQO + Moran : R2 = %.4f\n", summary(m2)$r.squared))
print(anova(m1, m2))

# Colinearidade de cada covariavel com as demais, que e o mecanismo proposto.
Xc <- sweep(X, 2, colMeans(X), "-")
r2_contra_resto <- vapply(seq_along(covariaveis), function(j) {
  summary(lm(Xc[, j] ~ Xc[, -j]))$r.squared
}, numeric(1))
names(r2_contra_resto) <- covariaveis
d <- merge(d, data.frame(covariavel = names(r2_contra_resto),
                         r2_resto = as.numeric(r2_contra_resto),
                         stringsAsFactors = FALSE), by = "covariavel")

cat("\n=== Colinearidade de cada covariavel com as outras 29 ===\n")
cat(sprintf("  correlacao entre I de Moran e R2 contra o resto: %.4f\n",
            cor(d$mi, d$r2_resto)))
cat(sprintf("  correlacao entre queda da PIP e R2 contra o resto: %.4f\n",
            cor(d$queda, d$r2_resto)))

m3 <- lm(queda ~ mqo + r2_resto, data = d)
cat("\n  queda ~ PIP do MQO + colinearidade:\n")
print(round(summary(m3)$coefficients, 5))
cat(sprintf("  R2: %.4f\n", summary(m3)$r.squared))

d <- d[order(-d$queda), ]
write.csv(d, file.path(AQUI, "saida", "mecanismo_dados_reais.csv"), row.names = FALSE)
cat("\nGravado em saida/mecanismo_dados_reais.csv\n")
