# Figuras do artigo, segunda geração: ggplot2 + sf, com identidade visual única.
#
# As figuras de base (34_figuras.R) documentam. Estas explicam. Cada uma carrega
# um argumento inteiro do artigo.
#
#   fig1_mapas          por que o problema existe, para quem nunca viu o método
#   fig2_suavidade      a versão quantitativa do mesmo problema
#   fig3_slopegraph     o que acontece com a seleção, covariável por covariável
#   fig4_dois_mundos    processo gerador conhecido contra dados reais
#   fig5_rho            o que governa a severidade
#   fig6_prioris        não é escolha de priori
#
# Saída em saida/figuras/ nos dois formatos: PDF para o texto, PNG para tela.
# Números com vírgula decimal, como nas tabelas; a versão em inglês é o 42.

source(if (file.exists("../caminhos.R")) "../caminhos.R" else "caminhos.R")

suppressPackageStartupMessages({
  library(ggplot2); library(sf); library(dplyr); library(tidyr)
  library(patchwork); library(ggrepel)
})

AQUI <- ESPACIAL
RAIZ <- ESTIMACOES
FIG  <- file.path(AQUI, "saida", "figuras")
dir.create(FIG, showWarnings = FALSE, recursive = TRUE)

# ---------------------------------------------------------------- identidade
TINTA  <- "#1C1C1C"
CINZA  <- "#9A9A9A"
QUENTE <- "#B4462F"
AMBAR  <- "#C98A3C"
FRIO   <- "#2E5E7E"
PAPEL  <- "#FBFAF8"

tema <- function(base = 11) {
  theme_minimal(base_size = base) +
    theme(
      plot.background  = element_rect(fill = PAPEL, colour = NA),
      panel.background = element_rect(fill = PAPEL, colour = NA),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(colour = "#E8E5E0", linewidth = 0.3),
      plot.title    = element_text(face = "bold", size = base * 1.15, colour = TINTA),
      plot.subtitle = element_text(size = base * 0.92, colour = "#5A5A5A",
                                   margin = margin(b = 10)),
      plot.caption  = element_text(size = base * 0.78, colour = "#7A7A7A",
                                   hjust = 0, margin = margin(t = 12)),
      axis.title    = element_text(size = base * 0.9, colour = "#4A4A4A"),
      axis.text     = element_text(colour = "#4A4A4A"),
      strip.text    = element_text(face = "bold", size = base * 0.95, colour = TINTA),
      legend.title  = element_text(size = base * 0.85, colour = "#4A4A4A"),
      legend.text   = element_text(size = base * 0.82, colour = "#4A4A4A"),
      plot.title.position = "plot",
      plot.caption.position = "plot")
}
# SEM_TITULO=1 grava sem título, subtítulo e nota de rodapé. Serve ao documento
# em LaTeX, onde essa informação vai na legenda e repeti-la dentro da figura é
# ruído. Fora dele a figura precisa se explicar sozinha.
SEM_TITULO <- nzchar(Sys.getenv("SEM_TITULO"))
if (SEM_TITULO) FIG <- file.path(AQUI, "saida", "figuras_sem_titulo")
dir.create(FIG, showWarnings = FALSE, recursive = TRUE)

gravar <- function(p, nome, largura, altura) {
  if (SEM_TITULO) {
    # Num patchwork, só o título geral sai: os títulos dos painéis internos são
    # o que identifica cada painel e precisam ficar. Numa figura simples, o
    # título é o que a legenda do LaTeX repete, e é ele que sai.
    if (inherits(p, "patchwork")) {
      p <- p + plot_annotation(title = NULL, subtitle = NULL, caption = NULL)
      altura <- altura - 0.6
    } else {
      p <- p + labs(title = NULL, subtitle = NULL, caption = NULL)
      altura <- altura - 1.1
    }
  }
  ggsave(file.path(FIG, paste0(nome, ".pdf")), p, width = largura, height = altura,
         device = cairo_pdf, bg = PAPEL)
  if (!SEM_TITULO)
    ggsave(file.path(FIG, paste0(nome, ".png")), p, width = largura, height = altura,
           dpi = 200, bg = PAPEL)
  cat("  gravada:", nome, "\n")
}
D <- function(p) { m <- mean(p); K <- length(p)
  sd(p) * sqrt((K - 1) / K) / sqrt(m * (1 - m)) }

# Vírgula decimal: `virg` para números em texto, `eixo` para rótulos de escala.
virg <- function(x, d = 2, sinal = FALSE)
  gsub(".", ",", sprintf(paste0("%", if (sinal) "+" else "", ".", d, "f"), x), fixed = TRUE)
eixo <- function(x) ifelse(is.na(x), NA_character_,
  gsub(".", ",", format(x, trim = TRUE, drop0trailing = TRUE, scientific = FALSE),
       fixed = TRUE))

DESF <- c(casos = "casos", obitos = "óbitos", letalidade = "letalidade")
grade <- lapply(names(DESF), function(d)
  readRDS(file.path(AQUI, "saida", "grade", paste0(d, "_queen.rds"))))
names(grade) <- names(DESF)
des <- readRDS(file.path(AQUI, "saida", "desenho_simulacao.rds"))
source(file.path(AQUI, "rotulos.R"))
pt <- function(v) rotular(v, "pt")
dados <- read.csv(file.path(RAIZ, "resultados/estimacoes/matrizes_estimacao/principal/casos.csv"),
                  colClasses = c(codigo_ibge = "character"), stringsAsFactors = FALSE)

# ============================================================ fig 1: os mapas
cat("figura 1: mapas\n")
malha <- readRDS(file.path(AQUI, "saida", "malha_ce.rds"))
mapa <- malha %>%
  mutate(codigo_ibge = as.character(code_muni)) %>%
  left_join(dados[, c("codigo_ibge", "z_log_casos100k", "z_idosos", "z_tax.agua")],
            by = "codigo_ibge")
stopifnot(sum(is.na(mapa$z_idosos)) == 0)

r_com_y <- function(v) cor(dados$z_log_casos100k, dados[[v]])
paineis <- list(
  list("z_log_casos100k", "A variável dependente", "casos por 100 mil, padronizado",
       sprintf("I de Moran %s", virg(des$mi_y_real, 3))),
  list("z_idosos", "Alta autocorrelação espacial", pt("z_idosos"),
       sprintf("I de Moran %s\nr com a variável dependente %s",
               virg(des$moran[["z_idosos"]], 3), virg(r_com_y("z_idosos"), 2, TRUE))),
  list("z_tax.agua", "Baixa autocorrelação espacial", pt("z_tax.agua"),
       sprintf("I de Moran %s\nr com a variável dependente %s",
               virg(des$moran[["z_tax.agua"]], 3), virg(r_com_y("z_tax.agua"), 2, TRUE))))

lim <- c(-1, 1) * max(abs(range(c(mapa$z_log_casos100k, mapa$z_idosos,
                                  mapa$`z_tax.agua`), na.rm = TRUE)))
fazer_mapa <- function(cfg) {
  ggplot(mapa) +
    geom_sf(aes(fill = .data[[cfg[[1]]]]), colour = "white", linewidth = 0.07) +
    coord_sf(expand = FALSE, datum = NA) +
    scale_fill_gradient2(low = FRIO, mid = "#F4F1EC", high = QUENTE, midpoint = 0,
                         limits = lim, name = "desvios-padrão", labels = eixo,
                         guide = guide_colourbar(barheight = 0.45, barwidth = 11,
                                                 title.position = "top",
                                                 title.hjust = 0.5)) +
    labs(title = cfg[[2]], subtitle = paste0(cfg[[3]], "\n", cfg[[4]])) +
    tema() +
    theme(axis.text = element_blank(), axis.title = element_blank(),
          axis.ticks = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), legend.position = "bottom",
          plot.margin = margin(2, 2, 2, 2),
          plot.title = element_text(size = 11.5, face = "bold"),
          plot.subtitle = element_text(size = 8.6, margin = margin(b = 4)))
}
f1 <- (fazer_mapa(paineis[[1]]) | fazer_mapa(paineis[[2]]) | fazer_mapa(paineis[[3]])) +
  plot_layout(guides = "collect") +
  plot_annotation(
    title = "Casos por 100 mil e duas covariáveis, 184 municípios do Ceará",
    subtitle = paste("Escala comum aos três painéis, em desvios-padrão. I de Moran e",
                     "correlação com a variável dependente indicados em cada painel."),
    caption = paste("Malha municipal do IBGE, 2020, via geobr. Vizinhança Queen,",
                    "normalizada por linha. O I de Moran vai de 0, sem padrão espacial,",
                    "a 1."),
    theme = tema()) &
  theme(legend.position = "bottom")
gravar(f1, "fig1_mapas", 11, 5.6)

# ================================================== fig 2: suavidade e sinal
cat("figura 2: suavidade\n")
cv <- setdiff(names(dados), c("codigo_ibge", "municipio", "z_log_casos100k"))
sv <- data.frame(
  covariavel = cv,
  nome = pt(cv),
  r = vapply(cv, function(v) cor(dados$z_log_casos100k, dados[[v]]), numeric(1)),
  mi = as.numeric(des$moran[cv]),
  pip_mqo = grade$casos$res$MQO$pip[match(cv, grade$casos$covariaveis)],
  stringsAsFactors = FALSE) %>%
  mutate(ar = abs(r),
         selecionada = ifelse(pip_mqo >= 0.5, "selecionada pelo MQO",
                              "não selecionada"))
if (all(is.na(sv$pip_mqo))) {
  sv$pip_mqo <- grade$casos$tabela$mqo[match(sv$covariavel, grade$casos$tabela$covariavel)]
  sv$selecionada <- ifelse(sv$pip_mqo >= 0.5, "selecionada pelo MQO", "não selecionada")
}
rho_sv <- cor(sv$ar, sv$mi)

f2 <- ggplot(sv, aes(mi, ar)) +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE, colour = CINZA,
              fill = CINZA, alpha = 0.12, linewidth = 0.5, linetype = "22") +
  geom_point(aes(colour = selecionada, size = selecionada), alpha = 0.9) +
  geom_text_repel(data = filter(sv, selecionada == "selecionada pelo MQO" | mi > 0.55 | ar > 0.3),
                  aes(label = nome, colour = selecionada),
                  size = 2.8, box.padding = 0.45, min.segment.length = 0,
                  segment.colour = CINZA, segment.size = 0.2, show.legend = FALSE) +
  annotate("text", x = 0.02, y = 0.40, hjust = 0, size = 3.2, colour = TINTA,
           label = sprintf("correlação entre os dois eixos: %s", virg(rho_sv))) +
  scale_colour_manual(values = c("selecionada pelo MQO" = QUENTE,
                                 "não selecionada" = CINZA), name = NULL) +
  scale_size_manual(values = c("selecionada pelo MQO" = 3.2,
                               "não selecionada" = 2.2), guide = "none") +
  scale_x_continuous(labels = eixo) + scale_y_continuous(labels = eixo) +
  labs(title = "I de Moran da covariável contra sua correlação com a variável dependente",
       subtitle = paste("30 covariáveis candidatas. Em vermelho, as seis com PIP",
                        "acima de 0,50 no BMA sem termo espacial."),
       x = "I de Moran da covariável", y = "|correlação| com a variável dependente",
       caption = paste("Variável dependente: log de casos por 100 mil, padronizado. Reta por",
                       "mínimos quadrados, faixa de 95%.")) +
  tema()
gravar(f2, "fig2_suavidade", 9, 6.4)

# ==================================================== fig 3: o slopegraph
cat("figura 3: slopegraph\n")
sl <- bind_rows(lapply(names(DESF), function(d) {
  z <- grade[[d]]
  data.frame(desfecho = DESF[[d]], covariavel = z$tabela$covariavel,
             nome = pt(z$tabela$covariavel),
             mqo = z$tabela$mqo, sar = z$tabela$sar,
             base_mqo = z$res$MQO$tam / z$K, base_sar = z$res$SAR$tam / z$K,
             stringsAsFactors = FALSE)
})) %>%
  mutate(desfecho = factor(desfecho, levels = unname(DESF)),
         destino = factor(case_when(
           mqo >= 0.5 & sar < 0.5 ~ "sai da lista",
           mqo < 0.5 & sar >= 0.5 ~ "entra na lista",
           mqo >= 0.5 & sar >= 0.5 ~ "permanece na lista",
           TRUE ~ "não cruza o limiar"),
           levels = c("sai da lista", "entra na lista", "permanece na lista",
                      "não cruza o limiar")))
longo <- sl %>%
  pivot_longer(c(mqo, sar), names_to = "modelo", values_to = "pip") %>%
  mutate(modelo = factor(ifelse(modelo == "mqo", "MQO", "SAR"),
                         levels = c("MQO", "SAR")))
bases <- sl %>% distinct(desfecho, base_mqo, base_sar) %>%
  pivot_longer(c(base_mqo, base_sar), names_to = "modelo", values_to = "base") %>%
  mutate(x = ifelse(modelo == "base_mqo", 1, 2))
rot <- longo %>% filter(destino != "não cruza o limiar", modelo == "SAR")

f3 <- ggplot(longo, aes(x = modelo, y = pip, group = covariavel)) +
  geom_hline(yintercept = 0.5, colour = TINTA, linewidth = 0.4, linetype = "22") +
  geom_segment(data = bases,
               aes(x = x - 0.32, xend = x + 0.32, y = base, yend = base),
               inherit.aes = FALSE, colour = CINZA, linewidth = 0.5, linetype = "11") +
  geom_line(aes(colour = destino, linewidth = destino, alpha = destino)) +
  geom_point(aes(colour = destino, alpha = destino), size = 1.5) +
  geom_text_repel(data = rot, aes(label = nome, colour = destino),
                  nudge_x = 0.26, direction = "y", hjust = 0, size = 2.6,
                  segment.size = 0.2, segment.colour = CINZA,
                  min.segment.length = 0, box.padding = 0.28, force = 3,
                  max.overlaps = 30, seed = 1, show.legend = FALSE) +
  facet_wrap(~ desfecho, nrow = 1) +
  scale_colour_manual(values = c("sai da lista" = QUENTE, "entra na lista" = FRIO,
                                 "permanece na lista" = TINTA, "não cruza o limiar" = CINZA),
                      name = NULL) +
  scale_linewidth_manual(values = c(0.85, 0.85, 0.85, 0.35), guide = "none") +
  scale_alpha_manual(values = c(0.95, 0.95, 0.95, 0.45), guide = "none") +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25), labels = eixo) +
  scale_x_discrete(expand = expansion(add = c(0.35, 1.5))) +
  labs(title = "PIP sob MQO e sob SAR, por covariável e variável dependente",
       subtitle = paste("Tracejado: limiar de 0,50. Pontilhado: PIP de referência,",
                        "igual ao tamanho médio do modelo dividido por K."),
       x = NULL, y = "probabilidade posterior de inclusão",
       caption = paste("Duas cadeias de 1.000.000 de iterações, 200.000 de burn-in,",
                       "vizinhança Queen. Discordância máxima entre cadeias: 0,022.")) +
  tema() + theme(legend.position = "top", panel.grid.major.x = element_blank())
gravar(f3, "fig3_slopegraph", 12, 6)

# ================================================== fig 4: os dois mundos
cat("figura 4: dois mundos\n")
ts <- read.csv(file.path(AQUI, "saida", "tabela_simulacoes.csv"), stringsAsFactors = FALSE)
sim <- ts %>%
  mutate(cenario = recode(cenario, rugosas = "baixa autocorrelação",
                          suaves = "alta autocorrelação")) %>%
  mutate(rotulo = sprintf("%s, \u03c1 = %s", cenario, virg(rho)),
         mundo = "simulação, processo gerador conhecido") %>%
  select(rotulo, razao, mundo) %>%
  # O cenário de oito verdadeiras em rho = 0,53 foi executado duas vezes de forma
  # independente; sem sufixo os dois rótulos colidiriam no eixo.
  group_by(rotulo) %>%
  mutate(rotulo = if (n() > 1) sprintf("%s (%d)", rotulo, row_number()) else rotulo) %>%
  ungroup()
real <- data.frame(
  rotulo = sprintf("dados reais: %s", unname(DESF)),
  razao = vapply(names(DESF), function(d)
    D(grade[[d]]$res$SAR$pip) / D(grade[[d]]$res$MQO$pip), numeric(1)),
  mundo = "dados reais", stringsAsFactors = FALSE)
dm <- bind_rows(sim, real)
dm$rotulo <- factor(dm$rotulo, levels = dm$rotulo[order(dm$razao)])
faixa <- range(sim$razao)

f4 <- ggplot(dm, aes(x = razao, y = rotulo)) +
  annotate("rect", xmin = faixa[1], xmax = faixa[2], ymin = -Inf, ymax = Inf,
           fill = FRIO, alpha = 0.07) +
  geom_segment(aes(x = 0, xend = razao, yend = rotulo, colour = mundo),
               linewidth = 0.4, alpha = 0.45) +
  geom_point(aes(colour = mundo, shape = mundo), size = 3.4) +
  geom_text(aes(label = virg(razao)), hjust = -0.55, size = 2.9,
            colour = TINTA) +
  scale_colour_manual(values = c("dados reais" = QUENTE,
                                 "simulação, processo gerador conhecido" = FRIO), name = NULL) +
  scale_shape_manual(values = c("dados reais" = 18,
                                "simulação, processo gerador conhecido" = 16), name = NULL) +
  scale_x_continuous(limits = c(0, 1.12), breaks = seq(0, 1, 0.25), labels = eixo,
                     expand = expansion(mult = c(0, 0.02))) +
  labs(title = "Razão entre o D do SAR e o do MQO: cenários simulados e variáveis dependentes observadas",
       subtitle = sprintf(paste("Faixa sombreada: intervalo dos oito cenários com",
                                "processo gerador conhecido, de %s a %s."),
                          virg(faixa[1]), virg(faixa[2])),
       x = "razão de D entre SAR e MQO", y = NULL,
       caption = paste("D é o desvio-padrão das PIP dividido pelo seu máximo teórico",
                       "dado o tamanho médio do modelo. Dez replicações por cenário,",
                       "120.000 iterações\npor cadeia; nas variáveis dependentes observadas,",
                       "1.000.000 de iterações.")) +
  tema() + theme(legend.position = "top", panel.grid.major.y = element_blank())
gravar(f4, "fig4_dois_mundos", 9.8, 6.2)

# =========================================== fig 5: rho e a severidade
cat("figura 5: rho\n")
esp <- lapply(names(DESF), function(d)
  readRDS(file.path(AQUI, "saida", sprintf("especificacao_%s_queen.rds", d))))
names(esp) <- names(DESF)

# O w guardado pelo amostrador traz embutidos os pesos de Simpson (1,4,2,...,1).
# Eles servem à quadratura, não à densidade: plotados crus produzem um serrilhado
# de período dois. Dividir por eles devolve a densidade.
ng <- length(esp[[1]]$rho)
pesos_simpson <- c(1, rep(c(4, 2), length.out = ng - 2), 1)
densidade <- function(e) {
  w <- (e$a$w_rho + e$b$w_rho) / 2
  ex <- w / pesos_simpson
  h <- e$rho[2] - e$rho[1]
  ex / (h / 3 * sum(pesos_simpson * ex))
}
dens <- bind_rows(lapply(names(DESF), function(d) {
  data.frame(desfecho = DESF[[d]], rho = esp[[d]]$rho, dens = densidade(esp[[d]]),
             stringsAsFactors = FALSE)
})) %>% mutate(desfecho = factor(desfecho, levels = unname(DESF)))
resumo <- bind_rows(lapply(names(DESF), function(d) {
  e <- esp[[d]]; w <- (e$a$w_rho + e$b$w_rho) / 2; w <- w / sum(w)
  data.frame(desfecho = DESF[[d]], rho = sum(w * e$rho),
             pico = max(dens$dens[dens$desfecho == DESF[[d]]]),
             p_sar = (e$a$p_sar + e$b$p_sar) / 2, mi = grade[[d]]$mi_y,
             razao = D(grade[[d]]$res$SAR$pip) / D(grade[[d]]$res$MQO$pip),
             stringsAsFactors = FALSE)
})) %>% mutate(desfecho = factor(desfecho, levels = unname(DESF)))
cores <- setNames(c(QUENTE, AMBAR, FRIO), unname(DESF))

# A média vai no alto de cada curva: embaixo, colidia com os rótulos do eixo.
p5a <- ggplot(dens, aes(rho, dens, fill = desfecho, colour = desfecho)) +
  geom_area(position = "identity", alpha = 0.16, linewidth = 0.8) +
  geom_vline(data = resumo, aes(xintercept = rho, colour = desfecho),
             linetype = "22", linewidth = 0.5, show.legend = FALSE) +
  geom_text(data = resumo, aes(x = rho, y = pico, colour = desfecho,
                               label = virg(rho)),
            hjust = -0.2, vjust = -0.5, size = 3, show.legend = FALSE) +
  scale_fill_manual(values = cores, name = NULL) +
  scale_colour_manual(values = cores, name = NULL) +
  scale_x_continuous(limits = c(0, 1), labels = eixo) +
  scale_y_continuous(labels = eixo, expand = expansion(mult = c(0, 0.1))) +
  labs(title = "Posteriori de \u03c1 sob BMA",
       subtitle = "Tracejado vertical: média a posteriori.",
       x = expression(rho), y = "densidade") +
  tema() + theme(legend.position = "top")

comp <- resumo %>%
  select(desfecho, razao, mi, rho) %>%
  pivot_longer(c(mi, rho), names_to = "medida", values_to = "valor") %>%
  mutate(medida = factor(ifelse(medida == "mi", "I de Moran da vari\u00e1vel dependente",
                                "\u03c1 a posteriori"),
                         levels = c("I de Moran da vari\u00e1vel dependente", "\u03c1 a posteriori")))

p5b <- ggplot(comp, aes(valor, razao, colour = desfecho)) +
  geom_point(size = 4) +
  geom_text_repel(aes(label = desfecho), size = 2.8, box.padding = 0.6,
                  show.legend = FALSE, segment.colour = CINZA, segment.size = 0.2) +
  facet_wrap(~ medida, nrow = 1, scales = "free_x") +
  scale_colour_manual(values = cores, guide = "none") +
  scale_x_continuous(labels = eixo) +
  scale_y_continuous(limits = c(0.2, 0.85), labels = eixo) +
  labs(title = "Razão de D contra I de Moran e contra \u03c1",
       subtitle = "Três variáveis dependentes, vizinhança Queen.",
       x = NULL, y = "razão de D entre SAR e MQO") +
  tema()

f5 <- (p5a | p5b) + plot_layout(widths = c(1, 1.25)) +
  plot_annotation(
    caption = paste("\u03c1 estimado com a especificação dentro do espaço amostral, em",
                    "300.000 iterações por cadeia. Três pontos não constituem teste."),
    theme = tema())
gravar(f5, "fig5_rho", 12, 5.4)

# ================================================ fig 6: prioris
cat("figura 6: prioris\n")
rot_pri <- c(uniforme_hyperg = "uniforme sobre todos os modelos, hyper-g",
             binbeta5_hyperg = "Binomial-Beta, m = 5, hyper-g",
             binomial5_hyperg = "Binomial, \u03b8 = 5/K, hyper-g",
             binbeta_uip = "Binomial-Beta, m = K/2, g = n")
pri <- bind_rows(lapply(names(rot_pri), function(k) {
  z <- readRDS(file.path(AQUI, "saida", "priori", paste0("casos_queen_", k, ".rds")))
  data.frame(priori = rot_pri[[k]], k_mqo = z$res$MQO$tam,
             razao = z$D_sar / z$D_mqo, acima = sum(z$res$SAR$pip >= 0.5),
             base = z$res$SAR$tam / z$K, stringsAsFactors = FALSE)
}))
bg <- grade$casos
pri <- bind_rows(pri, data.frame(
  priori = "Binomial-Beta, m = K/2, hyper-g (principal)",
  k_mqo = bg$res$MQO$tam, razao = D(bg$res$SAR$pip) / D(bg$res$MQO$pip),
  acima = sum(bg$res$SAR$pip >= 0.5), base = bg$res$SAR$tam / bg$K,
  stringsAsFactors = FALSE)) %>%
  mutate(etiqueta = sprintf("%s\n%s sob SAR", priori,
                            ifelse(acima == 0, "nenhuma covariável acima de 0,50",
                                   sprintf("%d covariáveis acima de 0,50", acima))))

f6 <- ggplot(pri, aes(k_mqo, razao)) +
  geom_hline(yintercept = 1, colour = CINZA, linetype = "22", linewidth = 0.4) +
  annotate("text", x = max(pri$k_mqo), y = 1.03, hjust = 1, size = 2.9,
           colour = CINZA, label = "razão = 1") +
  geom_point(colour = QUENTE, size = 3.6, alpha = 0.9) +
  geom_text_repel(aes(label = etiqueta), size = 2.8, colour = TINTA,
                  box.padding = 0.9, min.segment.length = 0, lineheight = 1.05,
                  segment.colour = CINZA, segment.size = 0.2, seed = 1) +
  scale_y_continuous(limits = c(0, 1.1), breaks = seq(0, 1, 0.25), labels = eixo) +
  scale_x_continuous(limits = c(0, 17), labels = eixo) +
  labs(title = "Razão de D sob cinco combinações de priori de modelo e de g",
       subtitle = paste("Casos com vizinhança Queen, 1.000.000 de iterações por cadeia.",
                        "Rótulo: covariáveis acima de 0,50 sob SAR."),
       x = "tamanho médio do modelo sob MQO", y = "razão de D entre SAR e MQO",
       caption = paste("Sob priori uniforme a PIP de referência sobe para 0,489,",
                       "de modo que o limiar de 0,50 não separa covariáveis informativas",
                       "das demais.")) +
  tema()
gravar(f6, "fig6_prioris", 9.5, 6)

cat("\nTodas em", FIG, "\n")
