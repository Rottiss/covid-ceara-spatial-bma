# English figures, leaner than the Portuguese set: short subtitles, no captions.
# Meant to be readable on their own, without the paper.
#
# Output in saida/figures_en/.

source(if (file.exists("../caminhos.R")) "../caminhos.R" else "caminhos.R")

suppressPackageStartupMessages({
  library(ggplot2); library(sf); library(dplyr); library(tidyr)
  library(patchwork); library(ggrepel)
})

AQUI <- ESPACIAL
RAIZ <- ESTIMACOES
FIG  <- file.path(AQUI, "saida", "figures_en")
dir.create(FIG, showWarnings = FALSE, recursive = TRUE)
source(file.path(AQUI, "rotulos.R"))

INK   <- "#1C1C1C"
GREY  <- "#9A9A9A"
WARM  <- "#B4462F"
AMBER <- "#C98A3C"
COOL  <- "#2E5E7E"
PAPER <- "#FBFAF8"

tema <- function(base = 11) {
  theme_minimal(base_size = base) +
    theme(
      plot.background  = element_rect(fill = PAPER, colour = NA),
      panel.background = element_rect(fill = PAPER, colour = NA),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(colour = "#E8E5E0", linewidth = 0.3),
      plot.title    = element_text(face = "bold", size = base * 1.15, colour = INK),
      plot.subtitle = element_text(size = base * 0.92, colour = "#5A5A5A",
                                   margin = margin(b = 10)),
      axis.title    = element_text(size = base * 0.9, colour = "#4A4A4A"),
      axis.text     = element_text(colour = "#4A4A4A"),
      strip.text    = element_text(face = "bold", size = base * 0.95, colour = INK),
      legend.title  = element_text(size = base * 0.85, colour = "#4A4A4A"),
      legend.text   = element_text(size = base * 0.82, colour = "#4A4A4A"),
      plot.title.position = "plot")
}
save_fig <- function(p, name, w, h) {
  ggsave(file.path(FIG, paste0(name, ".pdf")), p, width = w, height = h,
         device = cairo_pdf, bg = PAPER)
  ggsave(file.path(FIG, paste0(name, ".png")), p, width = w, height = h,
         dpi = 200, bg = PAPER)
  cat("  saved:", name, "\n")
}
D <- function(p) { m <- mean(p); K <- length(p)
  sd(p) * sqrt((K - 1) / K) / sqrt(m * (1 - m)) }
en <- function(v) rotular(v, "en")

OUT <- c(casos = "cases", obitos = "deaths", letalidade = "case fatality")
grade <- lapply(names(OUT), function(d)
  readRDS(file.path(AQUI, "saida", "grade", paste0(d, "_queen.rds"))))
names(grade) <- names(OUT)
des <- readRDS(file.path(AQUI, "saida", "desenho_simulacao.rds"))
dados <- read.csv(file.path(RAIZ, "resultados/estimacoes/matrizes_estimacao/principal/casos.csv"),
                  colClasses = c(codigo_ibge = "character"), stringsAsFactors = FALSE)

# ------------------------------------------------------------------- fig 1
cat("figure 1: maps\n")
malha <- readRDS(file.path(AQUI, "saida", "malha_ce.rds"))
mapa <- malha %>%
  mutate(codigo_ibge = as.character(code_muni)) %>%
  left_join(dados[, c("codigo_ibge", "z_log_casos100k", "z_idosos", "z_tax.agua")],
            by = "codigo_ibge")
r_y <- function(v) cor(dados$z_log_casos100k, dados[[v]])
lim <- c(-1, 1) * max(abs(range(c(mapa$z_log_casos100k, mapa$z_idosos,
                                  mapa$`z_tax.agua`), na.rm = TRUE)))
painel <- function(col, title, sub) {
  ggplot(mapa) +
    geom_sf(aes(fill = .data[[col]]), colour = "white", linewidth = 0.07) +
    coord_sf(expand = FALSE, datum = NA) +
    scale_fill_gradient2(low = COOL, mid = "#F4F1EC", high = WARM, midpoint = 0,
                         limits = lim, name = "standard deviations",
                         guide = guide_colourbar(barheight = 0.45, barwidth = 11,
                                                 title.position = "top",
                                                 title.hjust = 0.5)) +
    labs(title = title, subtitle = sub) +
    tema() +
    theme(axis.text = element_blank(), axis.title = element_blank(),
          axis.ticks = element_blank(), panel.grid = element_blank(),
          legend.position = "bottom", plot.margin = margin(2, 2, 2, 2),
          plot.title = element_text(size = 11.5, face = "bold"),
          plot.subtitle = element_text(size = 8.6, margin = margin(b = 4)))
}
f1 <- (painel("z_log_casos100k", "The outcome",
              sprintf("COVID-19 cases per 100k\nMoran's I %.2f", des$mi_y_real)) |
       painel("z_idosos", "A smooth covariate",
              sprintf("%s\nMoran's I %.2f  |  r with outcome %+.2f",
                      en("z_idosos"), des$moran[["z_idosos"]], r_y("z_idosos"))) |
       painel("z_tax.agua", "A rough covariate",
              sprintf("%s\nMoran's I %.2f  |  r with outcome %+.2f",
                      en("z_tax.agua"), des$moran[["z_tax.agua"]], r_y("z_tax.agua")))) +
  plot_layout(guides = "collect") +
  plot_annotation(
    title = "COVID-19 cases and two covariates, 184 municipalities of Ceará, Brazil",
    subtitle = paste("Common scale across panels, in standard deviations. Moran's I",
                     "and correlation with the outcome shown per panel."),
    theme = tema()) &
  theme(legend.position = "bottom")
save_fig(f1, "fig1_maps", 11, 5.4)

# ------------------------------------------------------------------- fig 2
cat("figure 2: smoothness\n")
cv <- setdiff(names(dados), c("codigo_ibge", "municipio", "z_log_casos100k"))
sv <- data.frame(
  variavel = cv,
  nome = en(cv),
  r = vapply(cv, function(v) cor(dados$z_log_casos100k, dados[[v]]), numeric(1)),
  mi = as.numeric(des$moran[cv]), stringsAsFactors = FALSE) %>%
  mutate(ar = abs(r),
         pip = grade$casos$tabela$mqo[match(variavel, grade$casos$tabela$covariavel)],
         picked = ifelse(pip >= 0.5, "selected by OLS BMA", "not selected"))
rr <- cor(sv$ar, sv$mi)

f2 <- ggplot(sv, aes(mi, ar)) +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE, colour = GREY,
              fill = GREY, alpha = 0.12, linewidth = 0.5, linetype = "22") +
  geom_point(aes(colour = picked, size = picked), alpha = 0.9) +
  geom_text_repel(data = filter(sv, picked == "selected by OLS BMA" | mi > 0.55 | ar > 0.3),
                  aes(label = nome, colour = picked), size = 2.8, box.padding = 0.45,
                  min.segment.length = 0, segment.colour = GREY, segment.size = 0.2,
                  show.legend = FALSE) +
  annotate("text", x = 0.02, y = 0.40, hjust = 0, size = 3.2, colour = INK,
           label = sprintf("correlation between axes: %.2f", rr)) +
  scale_colour_manual(values = c("selected by OLS BMA" = WARM,
                                 "not selected" = GREY), name = NULL) +
  scale_size_manual(values = c("selected by OLS BMA" = 3.2, "not selected" = 2.2),
                    guide = "none") +
  labs(title = "Moran's I of each covariate against its correlation with the outcome",
       subtitle = paste("30 candidate covariates. In red, the six with PIP above 0.50",
                        "in the non-spatial BMA. Outcome: log cases per 100k."),
       x = "Moran's I of the covariate", y = "|correlation| with the outcome") +
  tema()
save_fig(f2, "fig2_smoothness", 9, 6.2)

# ------------------------------------------------------------------- fig 3
cat("figure 3: slopegraph\n")
sl <- bind_rows(lapply(names(OUT), function(d) {
  z <- grade[[d]]
  data.frame(outcome = OUT[[d]], nome = en(z$tabela$covariavel),
             mqo = z$tabela$mqo, sar = z$tabela$sar,
             base_mqo = z$res$MQO$tam / z$K, base_sar = z$res$SAR$tam / z$K,
             stringsAsFactors = FALSE)
})) %>%
  mutate(outcome = factor(outcome, levels = unname(OUT)),
         fate = factor(case_when(
           mqo >= 0.5 & sar < 0.5 ~ "drops out",
           mqo < 0.5 & sar >= 0.5 ~ "comes in",
           mqo >= 0.5 & sar >= 0.5 ~ "stays in",
           TRUE ~ "does not cross"),
           levels = c("drops out", "comes in", "stays in", "does not cross")))
lg <- sl %>% pivot_longer(c(mqo, sar), names_to = "model", values_to = "pip") %>%
  mutate(model = factor(ifelse(model == "mqo", "OLS", "SAR"), levels = c("OLS", "SAR")))
bs <- sl %>% distinct(outcome, base_mqo, base_sar) %>%
  pivot_longer(c(base_mqo, base_sar), names_to = "model", values_to = "base") %>%
  mutate(x = ifelse(model == "base_mqo", 1, 2))
lab <- lg %>% filter(fate != "does not cross", model == "SAR")

f3 <- ggplot(lg, aes(x = model, y = pip, group = nome)) +
  geom_hline(yintercept = 0.5, colour = INK, linewidth = 0.4, linetype = "22") +
  geom_segment(data = bs, aes(x = x - 0.32, xend = x + 0.32, y = base, yend = base),
               inherit.aes = FALSE, colour = GREY, linewidth = 0.5, linetype = "11") +
  geom_line(aes(colour = fate, linewidth = fate, alpha = fate)) +
  geom_point(aes(colour = fate, alpha = fate), size = 1.5) +
  geom_text_repel(data = lab, aes(label = nome, colour = fate), nudge_x = 0.26,
                  direction = "y", hjust = 0, size = 2.6, segment.size = 0.2,
                  segment.colour = GREY, min.segment.length = 0, box.padding = 0.28,
                  force = 3, max.overlaps = 30, seed = 1, show.legend = FALSE) +
  facet_wrap(~ outcome, nrow = 1) +
  scale_colour_manual(values = c("drops out" = WARM, "comes in" = COOL,
                                 "stays in" = INK, "does not cross" = GREY), name = NULL) +
  scale_linewidth_manual(values = c(0.85, 0.85, 0.85, 0.35), guide = "none") +
  scale_alpha_manual(values = c(0.95, 0.95, 0.95, 0.45), guide = "none") +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
  scale_x_discrete(expand = expansion(add = c(0.35, 1.5))) +
  labs(title = "PIP under OLS and under SAR, by covariate and outcome",
       subtitle = paste("Dashed: the 0.50 threshold. Dotted: the no-information line,",
                        "equal to average model size divided by K.\nTwo chains of",
                        "1,000,000 draws, 200,000 burn-in, Queen contiguity."),
       x = NULL, y = "posterior inclusion probability") +
  tema() + theme(legend.position = "top", panel.grid.major.x = element_blank())
save_fig(f3, "fig3_slopegraph", 12, 6)

# ------------------------------------------------------------------- fig 4
cat("figure 4: two worlds\n")
ts <- read.csv(file.path(AQUI, "saida", "tabela_simulacoes.csv"), stringsAsFactors = FALSE)
nomes_en <- c(rugosas = "rough covariates", suaves = "smooth covariates",
              controle = "control")
sim <- ts %>%
  mutate(nome = ifelse(cenario %in% names(nomes_en), nomes_en[cenario],
                       sub("realista, (\\w+) verd\\.",
                           "realistic, \\1 true covariates", cenario)),
         nome = sub("quatro", "four", nome), nome = sub("oito", "eight", nome),
         label = sprintf("%s, rho = %.2f", nome, rho),
         world = "simulation, truth known") %>%
  select(label, razao, world) %>%
  group_by(label) %>%
  mutate(label = if (n() > 1) sprintf("%s (%d)", label, row_number()) else label) %>%
  ungroup()
real <- data.frame(
  label = sprintf("real data: %s", unname(OUT)),
  razao = vapply(names(OUT), function(d)
    D(grade[[d]]$res$SAR$pip) / D(grade[[d]]$res$MQO$pip), numeric(1)),
  world = "real data", stringsAsFactors = FALSE)
dm <- bind_rows(sim, real)
dm$label <- factor(dm$label, levels = dm$label[order(dm$razao)])
band <- range(sim$razao)

f4 <- ggplot(dm, aes(x = razao, y = label)) +
  annotate("rect", xmin = band[1], xmax = band[2], ymin = -Inf, ymax = Inf,
           fill = COOL, alpha = 0.07) +
  geom_segment(aes(x = 0, xend = razao, yend = label, colour = world),
               linewidth = 0.4, alpha = 0.45) +
  geom_point(aes(colour = world, shape = world), size = 3.4) +
  geom_text(aes(label = sprintf("%.2f", razao)), hjust = -0.55, size = 2.9,
            colour = INK) +
  scale_colour_manual(values = c("real data" = WARM,
                                 "simulation, truth known" = COOL), name = NULL) +
  scale_shape_manual(values = c("real data" = 18,
                                "simulation, truth known" = 16), name = NULL) +
  scale_x_continuous(limits = c(0, 1.12), breaks = seq(0, 1, 0.25),
                     expand = expansion(mult = c(0, 0.02))) +
  labs(title = "Ratio of D, SAR to OLS: simulated scenarios and observed outcomes",
       subtitle = sprintf(paste("Shaded: the range across the eight known-truth",
                                "scenarios, %.2f to %.2f. D is the standard deviation",
                                "of the PIPs divided by\nits theoretical maximum given",
                                "average model size."), band[1], band[2]),
       x = "ratio of D, SAR to OLS", y = NULL) +
  tema() + theme(legend.position = "top", panel.grid.major.y = element_blank())
save_fig(f4, "fig4_two_worlds", 10, 6)

# ------------------------------------------------------------------- fig 5
cat("figure 5: rho\n")
esp <- lapply(names(OUT), function(d)
  readRDS(file.path(AQUI, "saida", sprintf("especificacao_%s_queen.rds", d))))
names(esp) <- names(OUT)
ng <- length(esp[[1]]$rho)
wS <- c(1, rep(c(4, 2), length.out = ng - 2), 1)
dens_of <- function(e) {
  w <- (e$a$w_rho + e$b$w_rho) / 2; ex <- w / wS; h <- e$rho[2] - e$rho[1]
  ex / (h / 3 * sum(wS * ex))
}
dens <- bind_rows(lapply(names(OUT), function(d)
  data.frame(outcome = OUT[[d]], rho = esp[[d]]$rho, dens = dens_of(esp[[d]]),
             stringsAsFactors = FALSE))) %>%
  mutate(outcome = factor(outcome, levels = unname(OUT)))
sm <- bind_rows(lapply(names(OUT), function(d) {
  e <- esp[[d]]; w <- (e$a$w_rho + e$b$w_rho) / 2; w <- w / sum(w)
  data.frame(outcome = OUT[[d]], rho = sum(w * e$rho), mi = grade[[d]]$mi_y,
             razao = D(grade[[d]]$res$SAR$pip) / D(grade[[d]]$res$MQO$pip),
             stringsAsFactors = FALSE)
})) %>% mutate(outcome = factor(outcome, levels = unname(OUT)))
cores <- setNames(c(WARM, AMBER, COOL), unname(OUT))

p5a <- ggplot(dens, aes(rho, dens, fill = outcome, colour = outcome)) +
  geom_area(position = "identity", alpha = 0.16, linewidth = 0.8) +
  geom_vline(data = sm, aes(xintercept = rho, colour = outcome),
             linetype = "22", linewidth = 0.5, show.legend = FALSE) +
  geom_text(data = sm, aes(x = rho, y = 0, colour = outcome,
                           label = sprintf("%.2f", rho)),
            vjust = 1.6, size = 3, show.legend = FALSE) +
  scale_fill_manual(values = cores, name = NULL) +
  scale_colour_manual(values = cores, name = NULL) +
  scale_x_continuous(limits = c(0, 1)) +
  coord_cartesian(clip = "off") +
  labs(title = "Model-averaged posterior of rho",
       subtitle = "Dashed vertical: posterior mean.", x = expression(rho),
       y = "density") +
  tema() + theme(legend.position = "top", plot.margin = margin(5, 5, 16, 5))

cmp <- sm %>% select(outcome, razao, `Moran's I of the outcome` = mi,
                     `posterior rho` = rho) %>%
  pivot_longer(c(`Moran's I of the outcome`, `posterior rho`),
               names_to = "measure", values_to = "value")
p5b <- ggplot(cmp, aes(value, razao, colour = outcome)) +
  geom_point(size = 4) +
  geom_text_repel(aes(label = outcome), size = 2.8, box.padding = 0.6,
                  show.legend = FALSE, segment.colour = GREY, segment.size = 0.2) +
  facet_wrap(~ measure, nrow = 1, scales = "free_x") +
  scale_colour_manual(values = cores, guide = "none") +
  scale_y_continuous(limits = c(0.2, 0.85)) +
  labs(title = "Ratio of D against Moran's I and against rho",
       subtitle = "Three outcomes, Queen contiguity.",
       x = NULL, y = "ratio of D, SAR to OLS") +
  tema()
f5 <- (p5a | p5b) + plot_layout(widths = c(1, 1.25))
save_fig(f5, "fig5_rho", 12, 5.2)

# ------------------------------------------------------------------- fig 6
cat("figure 6: priors\n")
rot <- c(uniforme_hyperg = "uniform over 2^K",
         binbeta5_hyperg = "Binomial-Beta, m = 5",
         binomial5_hyperg = "Binomial, theta = 5/K",
         binbeta_uip = "Binomial-Beta m = K/2, g = n")
pri <- bind_rows(lapply(names(rot), function(k) {
  z <- readRDS(file.path(AQUI, "saida", "priori", paste0("casos_queen_", k, ".rds")))
  data.frame(prior = rot[[k]], k_ols = z$res$MQO$tam, razao = z$D_sar / z$D_mqo,
             above = sum(z$res$SAR$pip >= 0.5), stringsAsFactors = FALSE)
}))
bg <- grade$casos
pri <- bind_rows(pri, data.frame(
  prior = "Binomial-Beta m = K/2, hyper-g (the paper's)",
  k_ols = bg$res$MQO$tam, razao = D(bg$res$SAR$pip) / D(bg$res$MQO$pip),
  above = sum(bg$res$SAR$pip >= 0.5), stringsAsFactors = FALSE)) %>%
  mutate(tag = sprintf("%s\n%s above 0.50 under SAR", prior,
                       ifelse(above == 0, "none", as.character(above))))

f6 <- ggplot(pri, aes(k_ols, razao)) +
  geom_hline(yintercept = 1, colour = GREY, linetype = "22", linewidth = 0.4) +
  annotate("text", x = max(pri$k_ols), y = 1.03, hjust = 1, size = 2.9,
           colour = GREY, label = "ratio = 1") +
  geom_point(colour = WARM, size = 3.6, alpha = 0.9) +
  geom_text_repel(aes(label = tag), size = 2.8, colour = INK, box.padding = 0.9,
                  min.segment.length = 0, lineheight = 1.05,
                  segment.colour = GREY, segment.size = 0.2) +
  scale_y_continuous(limits = c(0, 1.1), breaks = seq(0, 1, 0.25)) +
  scale_x_continuous(limits = c(0, 17)) +
  labs(title = "Ratio of D under five combinations of model prior and g prior",
       subtitle = paste("Cases, Queen contiguity, 1,000,000 draws per chain. Label:",
                        "covariates above 0.50 under SAR.\nUnder the uniform prior the",
                        "no-information line rises to 0.489, so the 0.50 threshold no",
                        "longer separates anything."),
       x = "average model size under OLS", y = "ratio of D, SAR to OLS") +
  tema()
save_fig(f6, "fig6_priors", 9.5, 5.8)

cat("\nAll in", FIG, "\n")
