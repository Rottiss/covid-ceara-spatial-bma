# -*- coding: utf-8 -*-
"""Figure 7: D ratios for the 60 spatial fits in the prior grid.

Input: saida/priori_sem/consolidado_sem.csv, consolidated by script 76.
Uses the styles and language setting of script 52.
"""

import importlib.util
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.ticker import FixedLocator
from figure_labels import text

AQUI = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location("figuras_v3", AQUI / "52_figuras_v3.py")
v3 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(v3)

CM, CORES, NOMES = v3.CM, v3.CORES, v3.NOMES
TINTA, CINZA, LIMIAR = v3.TINTA, v3.CINZA, v3.LIMIAR

ROTULOS = {
    "binbeta_hyperg": text("Beta-binomial, $m=K/2$\nhyper-g", "Binomial-Beta, $m=K/2$\nhyper-g"),
    "binbeta_uip": text("Beta-binomial, $m=K/2$\n$g=n$", "Binomial-Beta, $m=K/2$\n$g=n$"),
    "binbeta5_hyperg": text("Beta-binomial, $m=5$\nhyper-g", "Binomial-Beta, $m=5$\nhyper-g"),
    "binomial5_hyperg": "Binomial, $\\theta=5/K$\nhyper-g",
    "uniforme_hyperg": text("Uniform over $2^K$\nhyper-g", "Uniforme sobre $2^K$\nhyper-g"),
}


def fig7_priori_sem():
    df = pd.read_csv(AQUI / "saida" / "priori_sem" / "consolidado_sem.csv",
                     encoding="utf-8")
    ordem = list(ROTULOS)
    fig, eixos = plt.subplots(1, 3, figsize=(16 * CM, 7.4 * CM), facecolor="white",
                              sharey=True)
    fig.subplots_adjust(left=0.24, right=0.99, top=0.86, bottom=0.17, wspace=0.10)

    for ax, d in zip(eixos, ("casos", "obitos", "letalidade")):
        v3.estilo_eixo(ax, grade_y=False)
        cor = CORES[d]
        ax.axvline(1, color=LIMIAR, linewidth=0.8, linestyle=(0, (4, 3)), zorder=1)
        for j, comb in enumerate(ordem):
            y = len(ordem) - 1 - j
            ax.axhline(y, color=v3.GRADE, linewidth=0.5, zorder=0)
            for w, desloc in (("queen", 0.16), ("rook", -0.16)):
                r = df[(df.desfecho == d) & (df.combinacao == comb) & (df.vizinhos == w)]
                assert len(r) == 1, (d, comb, w)
                cheio = w == "queen"
                for col, marca in (("razao_sar", "o"), ("razao_sem", "^")):
                    ax.plot(r[col].iloc[0], y + desloc, marker=marca, markersize=4.6,
                            markerfacecolor=cor if cheio else "white",
                            markeredgecolor=cor, markeredgewidth=0.9, zorder=3,
                            linestyle="none")
        ax.set_title(NOMES[d], color=TINTA, pad=4)
        ax.set_xlim(0.1, 1.09)
        ax.xaxis.set_major_locator(FixedLocator([0.2, 0.4, 0.6, 0.8, 1]))
        ax.set_ylim(-0.6, len(ordem) - 0.4)
        ax.set_xlabel(text("$D$ ratio relative to OLS", "Razão de $D$ contra o MQO"))

    eixos[0].yaxis.set_major_locator(FixedLocator(range(len(ordem))))
    eixos[0].set_yticklabels([ROTULOS[c] for c in reversed(ordem)])
    eixos[0].tick_params(axis="y", length=0)

    chaves = [
        Line2D([], [], marker="o", linestyle="none", markersize=4.6, color=TINTA,
               markerfacecolor=TINTA, label="SAR"),
        Line2D([], [], marker="^", linestyle="none", markersize=4.6, color=TINTA,
               markerfacecolor=TINTA, label="SEM"),
        Line2D([], [], marker="s", linestyle="none", markersize=4.6, color=TINTA,
               markerfacecolor=TINTA, label="Queen"),
        Line2D([], [], marker="s", linestyle="none", markersize=4.6, color=TINTA,
               markerfacecolor="white", markeredgewidth=0.9, label="Rook"),
        Line2D([], [], color=LIMIAR, linewidth=0.8, linestyle=(0, (4, 3)),
               label=text("Ratio equal to 1", "Razão igual a 1")),
    ]
    fig.legend(handles=chaves, loc="lower center", ncol=5, frameon=False,
               bbox_to_anchor=(0.58, -0.01), handletextpad=0.35, columnspacing=1.4)
    v3.gravar(fig, "fig7_priori_sem")


if __name__ == "__main__":
    fig7_priori_sem()
