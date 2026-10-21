# Inspect the BMS definitions of the hyper-g UIP and binomial-beta priors
# used in the OLS/SAR comparison.

suppressPackageStartupMessages(library(BMS))

cat("=== BMS", as.character(packageVersion("BMS")), "===\n\n")

objetos <- ls(getNamespace("BMS"), all.names = TRUE)
cat("objetos com 'hyper':", paste(grep("hyper", objetos, value = TRUE, ignore.case = TRUE),
                                  collapse = ", "), "\n")
cat("objetos com 'gprior':", paste(grep("gprior", objetos, value = TRUE, ignore.case = TRUE),
                                   collapse = ", "), "\n")
cat("objetos com 'mprior':", paste(grep("mprior", objetos, value = TRUE, ignore.case = TRUE),
                                   collapse = ", "), "\n\n")

mostrar <- function(nome) {
  if (!exists(nome, envir = getNamespace("BMS"))) {
    cat("### ", nome, ": nao existe\n\n", sep = "")
    return(invisible(NULL))
  }
  obj <- get(nome, envir = getNamespace("BMS"))
  cat("### ", nome, "\n", sep = "")
  print(obj)
  cat("\n")
}

for (nome in c(".fixedg", ".hyperg", "gprior.hyperg", ".mprior.randomfixed.init",
               ".f21hyper", "f21hyper")) {
  mostrar(nome)
}

cat("=== Estrutura do objeto de priori g, para a configuracao do TCC ===\n")
set.seed(1)
n <- 60; k <- 5
Xf <- matrix(rnorm(n * k), n, k)
yf <- Xf[, 1] + rnorm(n)
aj <- BMS::bms(cbind(yf, Xf), mcmc = "enumerate", g = "hyper=UIP",
               mprior = "random", mprior.size = k / 2, nmodel = 2^k,
               user.int = FALSE)
cat("\ngprior.info:\n"); print(aj$gprior.info)
cat("\nmprior.info:\n"); print(aj$mprior.info)
cat("\ninfo:\n"); print(aj$info[1:8])
