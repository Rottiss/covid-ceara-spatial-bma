# -*- coding: utf-8 -*-
"""Figure 3: conditional PIPs under OLS, SAR, and SEM for 30 covariates.

Input: pips_queen.csv exported by script 50. Labels default to English;
set FIGURE_LANGUAGE=pt for Portuguese. Data identifiers are unchanged.
"""

import os
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import FixedLocator, FuncFormatter
import pandas as pd
from figure_labels import OUTCOMES, covariate_label, text

AQUI = Path(__file__).resolve().parent
DADOS = Path(os.environ.get("DADOS_FIGURAS", AQUI / "saida" / "figuras_v3" / "dados")) / "pips_queen.csv"
SAIDA = Path(os.environ.get("FIGURAS_SAIDA", AQUI / "saida" / "figuras_v3"))

CM = 1 / 2.54
LARGURA = 16 * CM
ALTURA = 19.5 * CM

TINTA = "#1a252f"
MOLDURA = "#2c3e50"
GRADE = "#e2e8f0"
CINZA = "#9aa5b1"
LIMIAR = "#9b2c26"
CORES = {"casos": "#1b4f72", "obitos": "#78281f", "letalidade": "#4a235a"}
TITULOS = OUTCOMES

# Display groups only; these do not affect estimation.
GRUPOS = [
    (text("Demography and geography", "Demografia e território"),
     ["z_pop", "z_area", "z_densidade", "z_pop.rural", "z_idosos", "metrop", "semiarido"]),
    (text("Health and sanitation", "Saúde e saneamento"),
     ["z_sus1k", "z_leitos1k", "z_prof1k", "z_imuni", "z_int.circ", "z_int.resp",
      "z_int.diab", "z_int.asma", "z_tax.agua", "z_desp.saude"]),
    (text("Economy and education", "Socioeconomia e educação"),
     ["z_pibpc", "z_emprego", "z_ex.pobr", "z_bolsaf", "z_ivs.infra", "z_ideb.fund",
      "z_ideb3", "z_desp.educ", "z_energia100k", "z_tax.hom"]),
    (text("Local politics", "Política local"), ["aliadogov", "z_votos", "z_eleitas.fem"]),
]


def clarear(hex_cor, fracao):
    """Mix with white: fracao=0 keeps the color; fracao=1 returns white."""
    h = hex_cor.lstrip("#")
    r, g, b = (int(h[i:i + 2], 16) for i in (0, 2, 4))
    mix = lambda c: round(c + (255 - c) * fracao)
    return "#{:02x}{:02x}{:02x}".format(mix(r), mix(g), mix(b))


def virgula(x, _pos=None):
    return text(f"{x:g}", f"{x:g}".replace(".", ","))


def main():
    df = pd.read_csv(DADOS, encoding="utf-8")
    todas = [c for _, cs in GRUPOS for c in cs]
    assert sorted(todas) == sorted(df["covariavel"].unique()), "Groups must cover all 30 covariates"
    assert len(todas) == len(set(todas)) == 30

    # Within each group, order by the maximum PIP across the nine
    # outcome/specification combinations.
    maxima = df.set_index("covariavel")[["pip_mqo", "pip_sar", "pip_sem"]].max(axis=1).groupby(level=0).max()
    linhas, cabecalhos, y = [], [], 0.0
    for nome, cs in GRUPOS:
        cabecalhos.append((y, nome))
        y += 1.0
        for c in sorted(cs, key=lambda c: -maxima[c]):
            linhas.append((y, c))
            y += 1.0
        y += 0.45
    y_max = y - 0.45
    pos = {c: yy for yy, c in linhas}
    rotulo = df.drop_duplicates("covariavel").set_index("covariavel")["rotulo"]

    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["DejaVu Sans", "Arial", "Helvetica"],
        "font.size": 8,
        "axes.edgecolor": MOLDURA,
        "axes.linewidth": 0.6,
        "axes.titlesize": 9,
        "axes.titleweight": "bold",
        "xtick.labelsize": 7.5,
        "ytick.labelsize": 7.5,
        "xtick.color": TINTA,
        "ytick.color": TINTA,
        "legend.fontsize": 7.5,
        "pdf.fonttype": 3,  # Embed only the glyphs used.
    })

    fig, eixos = plt.subplots(1, 3, figsize=(LARGURA, ALTURA), sharey=True, facecolor="white")
    fig.subplots_adjust(left=0.34, right=0.985, top=0.935, bottom=0.055, wspace=0.07)

    for ax, desfecho in zip(eixos, ["casos", "obitos", "letalidade"]):
        cor = CORES[desfecho]
        sub = df[df["desfecho"] == desfecho].set_index("covariavel")
        ref_sar = sub["ref_sar"].iloc[0]

        ax.set_facecolor("white")
        ax.grid(axis="x", color=GRADE, linewidth=0.5, linestyle="-", zorder=0)
        ax.axvline(0.5, color=LIMIAR, linestyle=(0, (4, 2)), linewidth=0.8, zorder=1)
        ax.axvline(ref_sar, color=CINZA, linestyle=(0, (1, 1.5)), linewidth=0.8, zorder=1)

        for c, r in sub.iterrows():
            yy = pos[c]
            # Highlight covariates with PIP >= 0.50 in any specification.
            destaque = max(r["pip_mqo"], r["pip_sar"], r["pip_sem"]) >= 0.5
            cm = cor if destaque else CINZA

            ax.plot([r["pip_mqo"], r["pip_sar"]], [yy, yy], color=cm,
                    linewidth=1.1 if destaque else 0.7, alpha=1 if destaque else 0.8,
                    solid_capstyle="round", zorder=2)
            ax.scatter(r["pip_sem"], yy, marker="^", s=14, color=clarear(cm, 0.35),
                       edgecolors="white", linewidths=0.4, zorder=3)
            ax.scatter(r["pip_mqo"], yy, marker="o", s=16, facecolors="white",
                       edgecolors=cm, linewidths=0.9, zorder=4)
            ax.scatter(r["pip_sar"], yy, marker="o", s=16, color=cm,
                       edgecolors="white", linewidths=0.4, zorder=5)

        ax.set_title(TITULOS[desfecho], color=TINTA, pad=4)
        ax.set_xlim(-0.02, 1.02)
        ax.xaxis.set_major_locator(FixedLocator([0, 0.25, 0.5, 0.75, 1]))
        ax.xaxis.set_major_formatter(FuncFormatter(virgula))
        ax.tick_params(axis="x", length=2, width=0.5, pad=2)
        ax.tick_params(axis="y", length=0, pad=3)
        for lado in ("top", "right"):
            ax.spines[lado].set_visible(True)

    eixos[0].set_ylim(y_max - 0.3, -0.7)
    eixos[0].set_yticks([yy for yy, _ in linhas])
    eixos[0].set_yticklabels([covariate_label(c, rotulo[c]) for _, c in linhas])
    for yy, nome in cabecalhos:
        eixos[0].text(-0.03, yy, nome, transform=eixos[0].get_yaxis_transform(),
                      ha="right", va="center", fontsize=7.5, fontweight="bold", color=TINTA)
    fig.text(0.66, 0.012, text("Posterior inclusion probability", "Probabilidade posterior de inclusão"), ha="center",
             fontsize=8, color=TINTA)

    chaves = [
        Line2D([], [], marker="o", linestyle="", markersize=4.2, markerfacecolor="white",
               markeredgecolor=TINTA, markeredgewidth=0.9, label=text("OLS", "MQO")),
        Line2D([], [], marker="o", linestyle="", markersize=4.2, color=TINTA, label="SAR"),
        Line2D([], [], marker="^", linestyle="", markersize=4.2, color=clarear(TINTA, 0.35), label="SEM"),
        Line2D([], [], color=LIMIAR, linestyle=(0, (4, 2)), linewidth=0.8, label=text("0.50", "0,50")),
        Line2D([], [], color=CINZA, linestyle=(0, (1, 1.5)), linewidth=0.8, label=text("Reference PIP (SAR)", "PIP de referência (SAR)")),
    ]
    fig.legend(handles=chaves, loc="upper center", ncol=5, frameon=False,
               bbox_to_anchor=(0.625, 1.0), handletextpad=0.4, columnspacing=1.3)

    SAIDA.mkdir(parents=True, exist_ok=True)
    for ext in ("pdf", "png"):
        fig.savefig(SAIDA / f"fig3_pip_especificacoes.{ext}", dpi=300, facecolor="white")
    print("Saved:", SAIDA / "fig3_pip_especificacoes.pdf")


if __name__ == "__main__":
    main()
