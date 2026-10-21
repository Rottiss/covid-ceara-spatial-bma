# Extrai a definicao exata da log-verossimilhanca marginal com priori hyper-g do
# BMS, e da priori Binomial-Beta, para reimplementa-las sem adivinhar.

suppressPackageStartupMessages(library(BMS))
ns <- getNamespace("BMS")

cat("################ .lprob.hyperg.init ################\n")
print(get(".lprob.hyperg.init", envir = ns))

cat("\n################ .f21_4hyperg ################\n")
print(get(".f21_4hyperg", envir = ns))

cat("\n################ .gprior.hyperg.init ################\n")
print(get(".gprior.hyperg.init", envir = ns))

cat("\n################ .mprior.randomt.init ################\n")
print(get(".mprior.randomt.init", envir = ns))
