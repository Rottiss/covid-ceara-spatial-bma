# -*- coding: utf-8 -*-
"""Figures 1, 2, 4, 5, and 6 from the tables exported by script 50.

Labels default to English; set FIGURE_LANGUAGE=pt for Portuguese.
Dependencies: numpy, pandas, and matplotlib.
"""

import json
import os
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.cm import ScalarMappable
from matplotlib.collections import PolyCollection
from matplotlib.colors import Normalize
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib.ticker import FixedLocator, FuncFormatter, NullLocator
from figure_labels import OUTCOMES, covariate_label, scenario_label, text

AQUI = Path(__file__).resolve().parent
DADOS = Path(os.environ.get("DADOS_FIGURAS", AQUI / "saida" / "figuras_v3" / "dados"))
SAIDA = Path(os.environ.get("FIGURAS_SAIDA", AQUI / "saida" / "figuras_v3"))
CM = 1 / 2.54

TINTA = "#1a252f"
MOLDURA = "#2c3e50"
GRADE = "#e2e8f0"
CINZA = "#9aa5b1"
LIMIAR = "#9b2c26"
CORES = {"casos": "#1b4f72", "obitos": "#78281f", "letalidade": "#4a235a"}
NOMES = OUTCOMES

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["DejaVu Sans", "Arial", "Helvetica"],
    "font.size": 8,
    "axes.edgecolor": MOLDURA,
    "axes.linewidth": 0.6,
    "axes.labelcolor": TINTA,
    "axes.labelsize": 8,
    "axes.titlesize": 8.5,
    "axes.titleweight": "bold",
    "xtick.labelsize": 7.5,
    "ytick.labelsize": 7.5,
    "xtick.color": TINTA,
    "ytick.color": TINTA,
    "legend.fontsize": 7.5,
    "axes.unicode_minus": True,
    "pdf.fonttype": 3,  # Embed only the glyphs used.
})


def virgula(x, _pos=None):
    return text(f"{x:g}", f"{x:g}".replace(".", ",")).replace("-", "\u2212")


def estilo_eixo(ax, grade_x=True, grade_y=True):
    ax.set_facecolor("white")
    if grade_x:
        ax.grid(axis="x", color=GRADE, linewidth=0.5, linestyle="-", zorder=0)
    if grade_y:
        ax.grid(axis="y", color=GRADE, linewidth=0.5, linestyle="-", zorder=0)
    ax.set_axisbelow(True)
    ax.tick_params(length=2, width=0.5, pad=2)
    ax.xaxis.set_major_formatter(FuncFormatter(virgula))
    ax.yaxis.set_major_formatter(FuncFormatter(virgula))


def gravar(fig, nome):
    SAIDA.mkdir(parents=True, exist_ok=True)
    for ext in ("pdf", "png"):
        fig.savefig(SAIDA / f"{nome}.{ext}", dpi=300, facecolor="white")
    plt.close(fig)
    print("Saved:", nome)


# ------------------------------------------------------------------ Figure 1
def fig1_mapas():
    with open(DADOS / "malha_ce.geojson", encoding="utf-8") as f:
        geo = json.load(f)
    colunas = [("z_log_casos100k", text("Log cases per 100,000", "Log de casos por 100 mil")),
               ("z_idosos", text("Population aged 60+", "Proporção de idosos")),
               ("z_tax_agua", text("Water-supply coverage", "Cobertura de água"))]
    poligonos, valores = [], {c: [] for c, _ in colunas}
    for feat in geo["features"]:
        geom = feat["geometry"]
        partes = [geom["coordinates"]] if geom["type"] == "Polygon" else geom["coordinates"]
        for parte in partes:
            poligonos.append(np.asarray(parte[0])[:, :2])
            for c, _ in colunas:
                valores[c].append(feat["properties"][c])
    lim = max(np.abs(np.asarray(v)).max() for v in valores.values())
    norma = Normalize(-lim, lim)
    cmap = plt.get_cmap("RdBu_r")

    fig, eixos = plt.subplots(1, 3, figsize=(16 * CM, 7.4 * CM), facecolor="white")
    fig.subplots_adjust(left=0.01, right=0.99, top=0.93, bottom=0.2, wspace=0.04)
    todos = np.vstack(poligonos)
    for ax, (c, titulo) in zip(eixos, colunas):
        col = PolyCollection(poligonos, array=np.asarray(valores[c]), cmap=cmap, norm=norma,
                             edgecolors=MOLDURA, linewidths=0.12)
        ax.add_collection(col)
        ax.set_xlim(todos[:, 0].min() - 0.05, todos[:, 0].max() + 0.05)
        ax.set_ylim(todos[:, 1].min() - 0.05, todos[:, 1].max() + 0.05)
        ax.set_aspect("equal")
        ax.set_axis_off()
        ax.set_title(titulo, color=TINTA, pad=3)
    cax = fig.add_axes([0.3, 0.1, 0.4, 0.035])
    barra = fig.colorbar(ScalarMappable(norm=norma, cmap=cmap), cax=cax, orientation="horizontal")
    barra.set_ticks([-2, -1, 0, 1, 2])
    barra.ax.xaxis.set_major_formatter(FuncFormatter(virgula))
    barra.ax.tick_params(length=2, width=0.5, labelsize=7.5, pad=2)
    barra.outline.set_linewidth(0.5)
    barra.outline.set_edgecolor(MOLDURA)
    barra.set_label(text("Standard deviations", "Desvios-padrão"), fontsize=7.5, color=TINTA, labelpad=2)
    gravar(fig, "fig1_mapas")


# ------------------------------------------------------------------ Figure 2
def fig2_suavidade():
    df = pd.read_csv(DADOS / "suavidade.csv", encoding="utf-8")
    df["abs_r"] = df["r_casos"].abs()
    sel = df["pip_mqo"] >= 0.5

    fig, ax = plt.subplots(figsize=(12 * CM, 8.5 * CM), facecolor="white")
    fig.subplots_adjust(left=0.12, right=0.97, top=0.9, bottom=0.14)
    estilo_eixo(ax)
    b, a = np.polyfit(df["moran"], df["abs_r"], 1)
    xs = np.array([0, 0.82])
    ax.plot(xs, a + b * xs, color=CINZA, linewidth=0.9, zorder=1)
    ax.scatter(df.loc[~sel, "moran"], df.loc[~sel, "abs_r"], s=20, color=CINZA,
               edgecolors="white", linewidths=0.5, zorder=2)
    ax.scatter(df.loc[sel, "moran"], df.loc[sel, "abs_r"], s=30, color=CORES["casos"],
               edgecolors="white", linewidths=0.6, zorder=3)
    posicao = {  # Label offsets in data units and text alignment.
        "z_idosos": (-0.015, 0.0, "right", "center"),
        "semiarido": (0.0, -0.016, "center", "top"),
        "metrop": (-0.018, 0.0, "right", "center"),
        "z_leitos1k": (0.0, 0.013, "center", "bottom"),
        "z_int.circ": (0.013, 0.0, "left", "center"),
        "z_ideb.fund": (-0.015, 0.0, "right", "center"),
    }
    for _, r in df[sel].iterrows():
        dx, dy, ha, va = posicao[r["covariavel"]]
        ax.text(r["moran"] + dx, r["abs_r"] + dy, covariate_label(r["covariavel"], r["rotulo"]), ha=ha, va=va,
                fontsize=7, color=TINTA, zorder=4,
                bbox=dict(facecolor="white", edgecolor="none", alpha=0.85, pad=0.5))
    ax.set_xlim(-0.02, 0.84)
    ax.set_ylim(-0.01, 0.42)
    ax.xaxis.set_major_locator(FixedLocator([0, 0.2, 0.4, 0.6, 0.8]))
    ax.yaxis.set_major_locator(FixedLocator([0, 0.1, 0.2, 0.3, 0.4]))
    ax.set_xlabel(text("Covariate Moran's I", "I de Moran da covariável"))
    ax.set_ylabel(text("Absolute correlation with\nlog cases per 100,000", "Módulo da correlação com\nlog de casos por 100 mil"))
    chaves = [Line2D([], [], marker="o", linestyle="", markersize=5, color=CORES["casos"],
                     markeredgecolor="white", label=text("OLS PIP ≥ 0.50", "PIP ≥ 0,50 no MQO")),
              Line2D([], [], marker="o", linestyle="", markersize=4.2, color=CINZA,
                     markeredgecolor="white", label=text("OLS PIP < 0.50", "PIP < 0,50 no MQO")),
              Line2D([], [], color=CINZA, linewidth=0.9, label=text("Least-squares fit", "Mínimos quadrados"))]
    fig.legend(handles=chaves, loc="upper center", ncol=3, frameon=False,
               bbox_to_anchor=(0.545, 1.0), handletextpad=0.4, columnspacing=1.2)
    gravar(fig, "fig2_suavidade")


# ------------------------------------------------------------------ Figure 4
def fig4_razao_d():
    reps = pd.read_csv(DADOS / "simulacao_reps.csv", encoding="utf-8")
    cen = pd.read_csv(DADOS / "simulacao_cenarios.csv", encoding="utf-8")
    sim = cen[cen["tipo"] == "simulado"].reset_index(drop=True)
    obs = cen[cen["tipo"] == "observado"].reset_index(drop=True)
    chave_desf = {"Casos": "casos", "Óbitos": "obitos", "Letalidade": "letalidade"}

    y_sim = {r: i for i, r in enumerate(sim["rotulo"])}
    y0 = len(sim) + 0.7
    y_obs = {r: y0 + i for i, r in enumerate(obs["rotulo"])}

    fig, ax = plt.subplots(figsize=(16 * CM, 9.2 * CM), facecolor="white")
    fig.subplots_adjust(left=0.33, right=0.97, top=0.9, bottom=0.12)
    estilo_eixo(ax, grade_y=False)
    ax.axvline(1, color=CINZA, linewidth=0.8, linestyle=(0, (1, 1.5)), zorder=1)
    ax.axhline(len(sim) - 0.15, color=GRADE, linewidth=0.8, zorder=0)

    rng = np.random.default_rng(20260911)
    for rot, g in reps.groupby("rotulo", sort=False):
        yy = y_sim[rot] + rng.uniform(-0.18, 0.18, len(g))
        ax.scatter(g["razao"], yy, s=9, color=CINZA, alpha=0.9, edgecolors="none", zorder=2)
    ax.scatter(sim["razao"], [y_sim[r] for r in sim["rotulo"]], s=34, color=TINTA,
               edgecolors="white", linewidths=0.6, zorder=3)
    for _, r in obs.iterrows():
        ax.scatter(r["razao"], y_obs[r["rotulo"]], s=40, marker="D",
                   color=CORES[chave_desf[r["rotulo"]]], edgecolors="white", linewidths=0.6, zorder=3)

    ax.set_yticks(list(y_sim.values()) + list(y_obs.values()))
    ax.set_yticklabels([scenario_label(r) for r in y_sim] +
                       [text(f"Observed: {NOMES[chave_desf[r]].lower()}", f"Dados observados: {r.lower()}") for r in y_obs])
    ax.yaxis.set_minor_locator(NullLocator())
    ax.tick_params(axis="y", length=0)
    ax.set_ylim(y0 + len(obs) - 0.4, -0.6)
    ax.set_xlim(0, 1.2)
    ax.xaxis.set_major_locator(FixedLocator([0, 0.25, 0.5, 0.75, 1]))
    ax.set_xlabel(text("D ratio (SAR/OLS)", "Razão entre o D do SAR e o do MQO"))
    chaves = [Line2D([], [], marker="o", linestyle="", markersize=3.2, color=CINZA, label=text("Replication", "Replicação")),
              Line2D([], [], marker="o", linestyle="", markersize=5.2, color=TINTA,
                     markeredgecolor="white", label=text("Scenario", "Cenário")),
              Line2D([], [], marker="D", linestyle="", markersize=4.6, color=TINTA,
                     markeredgecolor="white", label=text("Observed data", "Dados observados"))]
    fig.legend(handles=chaves, loc="upper center", ncol=3, frameon=False,
               bbox_to_anchor=(0.65, 1.0), handletextpad=0.4, columnspacing=1.4)
    gravar(fig, "fig4_razao_d")


# ------------------------------------------------------------------ Figure 5
def fig5_rho():
    dens = pd.read_csv(DADOS / "rho_densidade.csv", encoding="utf-8")
    res = pd.read_csv(DADOS / "rho_resumo.csv", encoding="utf-8").set_index("desfecho")

    fig = plt.figure(figsize=(16 * CM, 6.6 * CM), facecolor="white")
    grade = fig.add_gridspec(1, 3, width_ratios=[1.5, 1, 1], left=0.07, right=0.985,
                             top=0.86, bottom=0.18, wspace=0.42)
    ax_d = fig.add_subplot(grade[0])
    ax_m = fig.add_subplot(grade[1])
    ax_r = fig.add_subplot(grade[2], sharey=ax_m)

    estilo_eixo(ax_d)
    for d in ("casos", "obitos", "letalidade"):
        sub = dens[(dens["desfecho"] == d) & (dens["rho"] >= -0.2)]
        ax_d.fill_between(sub["rho"], 0, sub["densidade"], color=CORES[d], alpha=0.12, linewidth=0, zorder=1)
        ax_d.plot(sub["rho"], sub["densidade"], color=CORES[d], linewidth=1.2, zorder=2)
        m = res.loc[d, "media"]
        ax_d.plot([m, m], [0, np.interp(m, sub["rho"], sub["densidade"])], color=CORES[d],
                  linewidth=0.8, zorder=2)
    ax_d.set_xlim(-0.2, 1)
    ax_d.set_ylim(0, None)
    ax_d.xaxis.set_major_locator(FixedLocator([0, 0.25, 0.5, 0.75, 1]))
    ax_d.set_xlabel("ρ")
    ax_d.set_ylabel(text("Posterior density", "Densidade a posteriori"))

    for ax, col, rotulo, lim in ((ax_m, "moran_y", text("Outcome Moran's I", "I de Moran da variável dependente"), (0.05, 0.6)),
                                 (ax_r, "media", text("Posterior mean of ρ", "ρ a posteriori"), (0.15, 0.8))):
        estilo_eixo(ax)
        for d in ("casos", "obitos", "letalidade"):
            ax.scatter(res.loc[d, col], res.loc[d, "razao_D"], s=36, color=CORES[d],
                       edgecolors="white", linewidths=0.6, zorder=3)
        ax.set_xlim(*lim)
        ax.set_xlabel(rotulo)
    ax_m.set_ylim(0.2, 0.8)
    ax_m.yaxis.set_major_locator(FixedLocator([0.2, 0.4, 0.6, 0.8]))
    ax_m.set_ylabel(text("D ratio (SAR/OLS)", "Razão de D (SAR/MQO)"))
    plt.setp(ax_r.get_yticklabels(), visible=False)
    ax_m.xaxis.set_major_locator(FixedLocator([0.1, 0.3, 0.5]))
    ax_r.xaxis.set_major_locator(FixedLocator([0.2, 0.4, 0.6, 0.8]))

    chaves = [Line2D([], [], color=CORES[d], linewidth=1.6, label=NOMES[d]) for d in CORES]
    fig.legend(handles=chaves, loc="upper center", ncol=3, frameon=False,
               bbox_to_anchor=(0.53, 1.0), handletextpad=0.5, columnspacing=1.6)
    gravar(fig, "fig5_rho")


# ------------------------------------------------------------------ Figure 6
def fig6_priori_rho():
    df = pd.read_csv(DADOS / "priori_rho.csv", encoding="utf-8")
    ticks = [1.01, 1.5, 2, 5, 10]

    fig, eixos = plt.subplots(2, 2, figsize=(16 * CM, 11 * CM), facecolor="white")
    fig.subplots_adjust(left=0.09, right=0.98, top=0.9, bottom=0.11, hspace=0.42, wspace=0.22)
    paineis = [("casos", eixos[0, 0]), ("obitos", eixos[0, 1]), ("letalidade", eixos[1, 0])]
    for d, ax in paineis:
        sub = df[df["desfecho"] == d].sort_values("a0")
        estilo_eixo(ax)
        ax.set_xscale("log")
        if d == "letalidade":
            ax.axhline(0, color=CINZA, linewidth=0.8, zorder=1)
        ax.fill_between(sub["a0"], sub["log_bf_cond_min"], sub["log_bf_cond_max"],
                        color=CORES[d], alpha=0.14, linewidth=0, zorder=1)
        ax.plot(sub["a0"], sub["log_bf_agregado_mc"], color=CORES[d], linewidth=1.3,
                marker="o", markersize=3.2, markeredgecolor="white", markeredgewidth=0.5, zorder=3)
        ax.set_title(NOMES[d], color=TINTA, pad=3)
        ax.set_ylabel(text("Log Bayes factor", "Log do fator de Bayes"))

    ax = eixos[1, 1]
    estilo_eixo(ax)
    ax.set_xscale("log")
    for d in ("casos", "obitos", "letalidade"):
        sub = df[df["desfecho"] == d].sort_values("a0")
        ax.plot(sub["a0"], sub["p_sar_agregada_mc"], color=CORES[d], linewidth=1.3,
                marker="o", markersize=3.2, markeredgecolor="white", markeredgewidth=0.5, zorder=3)
    ax.set_ylim(0.6, 1.03)
    ax.yaxis.set_major_locator(FixedLocator([0.6, 0.7, 0.8, 0.9, 1]))
    ax.set_title(text("Posterior probability of SAR", "Probabilidade a posteriori do SAR"), color=TINTA, pad=3)
    ax.text(10, 0.975, text("Cases and deaths", "Casos e óbitos"), ha="right", va="top", fontsize=7.5, color=TINTA)
    ax.text(10, 0.775, NOMES["letalidade"], ha="right", va="top", fontsize=7.5, color=TINTA)

    for ax in eixos.flat:
        ax.xaxis.set_major_locator(FixedLocator(ticks))
        ax.xaxis.set_minor_locator(NullLocator())
        ax.xaxis.set_major_formatter(FuncFormatter(virgula))
        ax.set_xlim(0.95, 11)
    for ax in eixos[1]:
        ax.set_xlabel(text("$a_0$ of the Beta$(a_0, a_0)$ prior for ρ", "$a_0$ da priori Beta$(a_0, a_0)$ para ρ"))

    chaves = [Line2D([], [], color=TINTA, linewidth=1.3, marker="o", markersize=3.2,
                     markeredgecolor="white", label=text("Aggregate estimate", "Estimativa agregada")),
              Patch(facecolor=TINTA, alpha=0.14, label=text("Range across covariate sets", "Faixa entre conjuntos de covariáveis"))]
    fig.legend(handles=chaves, loc="upper center", ncol=2, frameon=False,
               bbox_to_anchor=(0.53, 1.0), handletextpad=0.4, columnspacing=1.1)
    gravar(fig, "fig6_priori_rho")


if __name__ == "__main__":
    fig1_mapas()
    fig2_suavidade()
    fig4_razao_d()
    fig5_rho()
    fig6_priori_rho()
