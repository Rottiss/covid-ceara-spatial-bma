#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
gerar_figuras_simulacao_mecanismo.py
-----------------------------------
Gera duas figuras para comparar resultados observados e simulados:
1. fig_dois_mundos_executivo (.png / .pdf):
   Comparação dos 8 cenários simulados de Monte Carlo (distribuições empíricas
   via boxenplot/stripplot) contra os 3 desfechos reais (Casos, Óbitos, Letalidade),
   com faixas de razão de D abaixo de 0,50 e acima de 0,55.
2. fig_mecanismo_jointplot (.png / .pdf):
   Scatterplot com densidades marginais KDE:
   Eixo X: I de Moran da covariável
   Eixo Y: |r| com o log de casos por 100 mil
   Tamanho: R² contra as outras 29 covariáveis (colinearidade interna)
   Cor: Queda absoluta de PIP sob o modelo SAR (gradiente verde -> coral)
   Anotações detalhadas nas 6 covariáveis centrais.

"""

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from caminhos import ESTIMACOES

import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.ticker as ticker
import seaborn as sns

# ==============================================================================
# CONFIGURAÇÕES GERAIS E CAMINHOS
# ==============================================================================
BASE_DIR = str(Path(__file__).resolve().parent)
EXP_DIR = os.path.join(BASE_DIR, "saida", "export_python")
SAIDA_DIR = os.path.join(BASE_DIR, "saida", "figuras_novas")
CASOS_CSV = str(ESTIMACOES / "resultados/estimacoes/matrizes_estimacao/principal/casos.csv")
RESUMO_GRADE_CSV = os.path.join(BASE_DIR, "saida", "resumo_grade.csv")

os.makedirs(SAIDA_DIR, exist_ok=True)

# Tipografia e estética global
plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Segoe UI", "Arial", "Calibri", "DejaVu Sans"]
plt.rcParams["mathtext.fontset"] = "dejavusans"
plt.rcParams["text.color"] = "#1C1C1C"
plt.rcParams["axes.labelcolor"] = "#2B2D42"
plt.rcParams["xtick.color"] = "#4A5568"
plt.rcParams["ytick.color"] = "#4A5568"
plt.rcParams["axes.edgecolor"] = "#CBD5E0"
plt.rcParams["axes.linewidth"] = 0.8

# Função auxiliar para formatação em português (vírgula decimal)
def pt_num(val, decimals=2):
    return f"{val:.{decimals}f}".replace(".", ",")


# ==============================================================================
# FIGURA 1: COMPARAÇÃO DE DADOS OBSERVADOS E SIMULAÇÕES
# ==============================================================================
def gerar_figura_dois_mundos():
    print("\n[1/2] Gerando 'fig_dois_mundos_executivo'...")
    
    # 1. Carregar dados de simulação
    p_sim = os.path.join(EXP_DIR, "simulacoes_detalhadas.csv")
    df_sim = pd.read_csv(p_sim, encoding="utf-8")
    
    cenario_labels = {
        0: "Sintético: Covariáveis Suaves (ρ = 0,00, k = 4)",
        1: "Sintético: Covariáveis Rugosas (ρ = 0,00, k = 4)",
        2: "Sintético: Controle (ρ = 0,60, k = 4)",
        3: "Realista: 4 verdadeiras (ρ = 0,53, k = 4)",
        4: "Realista: 8 verdadeiras (ρ = 0,53, k = 8) — Lote A",
        5: "Realista: 8 verdadeiras (ρ = 0,53, k = 8) — Lote B",
        6: "Realista: 4 verdadeiras (ρ = 0,70, k = 4)",
        7: "Realista: 8 verdadeiras (ρ = 0,70, k = 8)",
    }
    
    df_sim["bloco_id"] = df_sim.index // 10
    df_sim["cenario_formatado"] = df_sim["bloco_id"].map(cenario_labels)
    
    # Calcular médias para ordenação
    medias = df_sim.groupby("cenario_formatado")["razao_D"].mean()
    ordem_sim = medias.sort_values(ascending=True).index.tolist()
    
    # 2. Carregar dados reais (resumo_grade.csv)
    df_grade = pd.read_csv(RESUMO_GRADE_CSV)
    df_real = df_grade[df_grade["vizinhos"] == "queen"].copy()
    
    real_map = {
        "casos": ("Casos COVID-19 (Taxa / 100k)", 0.300672, "#C0392B"),
        "letalidade": ("Letalidade COVID-19", 0.706365, "#2E5E7E"),
        "obitos": ("Óbitos COVID-19 (Taxa / 100k)", 0.709281, "#2E5E7E"),
    }
    
    # Organizar os paineis em duas linhas
    fig = plt.figure(figsize=(14.0, 10.8), dpi=300)
    gs = fig.add_gridspec(2, 1, height_ratios=[3.0, 7.0], hspace=0.32)
    ax_real = fig.add_subplot(gs[0])
    ax_sim = fig.add_subplot(gs[1], sharex=ax_real)
    
    # Configurar faixas conceituais sombreadas
    for ax in [ax_real, ax_sim]:
        ax.axvspan(0.0, 0.50, color="#FDE8E8", alpha=0.55, zorder=0)
        ax.axvspan(0.50, 0.55, color="#F8F9FA", alpha=0.45, zorder=0)
        ax.axvspan(0.55, 1.25, color="#EBF5FB", alpha=0.50, zorder=0)
        
        # Reference and parity lines
        ax.axvline(0.50, color="#E53E3E", linestyle="--", linewidth=1.3, alpha=0.85, zorder=1)
        ax.axvline(0.55, color="#718096", linestyle=":", linewidth=1.1, alpha=0.7, zorder=1)
        ax.axvline(1.00, color="#4A5568", linestyle="-.", linewidth=1.2, alpha=0.8, zorder=1)
        
        ax.grid(axis="x", color="#CBD5E0", linestyle=":", linewidth=0.7, alpha=0.7, zorder=0)
        ax.set_facecolor("#FFFFFF")
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.spines["left"].set_color("#CBD5E0")
        ax.spines["bottom"].set_color("#CBD5E0")

    # ==================== PAINEL REAL (TOP) ====================
    desfechos_ordem = [
        ("casos", 2),
        ("letalidade", 1),
        ("obitos", 0)
    ]
    
    for nome, y_idx in desfechos_ordem:
        rotulo, val, cor = real_map[nome]
        
        # Linha horizontal de suporte com degradê ou alpha
        ax_real.hlines(y_idx, xmin=0, xmax=val, color=cor, linestyle="-", linewidth=2.4, alpha=0.35, zorder=2)
        
        # Ponto em destaque (Losango/Estrela)
        marker = "D" if nome != "casos" else "*"
        msize = 15 if nome != "casos" else 24
        ax_real.plot(val, y_idx, marker=marker, markersize=msize, color=cor, 
                     markeredgecolor="white", markeredgewidth=2.2, zorder=4)
        
        # Texto com o valor exato posicionado estrategicamente
        if nome == "casos":
            ax_real.text(val + 0.02, y_idx - 0.28, f"Razão D = {pt_num(val, 3)}", 
                         va="center", ha="left", fontsize=10.5, fontweight="bold", color=cor,
                         bbox=dict(boxstyle="round,pad=0.25", fc="#FFFFFF", ec=cor, lw=1.2, alpha=0.95))
        elif nome == "letalidade":
            ax_real.text(val + 0.025, y_idx + 0.12, f"Razão D = {pt_num(val, 3)}", 
                         va="center", ha="left", fontsize=10.0, fontweight="bold", color=cor,
                         bbox=dict(boxstyle="round,pad=0.25", fc="#FFFFFF", ec=cor, lw=1.2, alpha=0.95))
        else: # obitos
            ax_real.text(val + 0.025, y_idx - 0.12, f"Razão D = {pt_num(val, 3)}", 
                         va="center", ha="left", fontsize=10.0, fontweight="bold", color=cor,
                         bbox=dict(boxstyle="round,pad=0.25", fc="#FFFFFF", ec=cor, lw=1.2, alpha=0.95))

    ax_real.set_yticks([0, 1, 2])
    ax_real.set_yticklabels(["Óbitos COVID-19", "Letalidade COVID-19", "Casos COVID-19"], 
                            fontsize=11, fontweight="bold")
    ax_real.set_ylim(-0.7, 2.8)
    
    # Annotations for observed data
    ax_real.annotate(
        "CASOS: RAZÃO DE D = 0,301\n"
        "D no SAR corresponde a cerca de 30% de D no MQO.\n"
        "A razão observada fica abaixo das médias dos cenários mostrados.",
        xy=(0.301, 2.0), xytext=(0.04, 1.32),
        arrowprops=dict(arrowstyle="->", color="#C0392B", lw=1.8, connectionstyle="arc3,rad=-0.15"),
        fontsize=9.2, fontweight="bold", color="#900C3F",
        bbox=dict(boxstyle="round,pad=0.45", fc="#FDEDEC", ec="#E6B0AA", lw=1.4, alpha=0.96),
        zorder=5
    )
    
    ax_real.annotate(
        "ÓBITOS E LETALIDADE\n"
        "Razão de D aproximadamente igual a 0,71.\n"
        "Valores dentro da faixa dos cenários mostrados.",
        xy=(0.708, 1.0), xytext=(0.73, 1.82),
        arrowprops=dict(arrowstyle="->", color="#2E5E7E", lw=1.6, connectionstyle="arc3,rad=0.15"),
        fontsize=9.0, fontweight="semibold", color="#1F3A52",
        bbox=dict(boxstyle="round,pad=0.45", fc="#EBF5FB", ec="#AED6F1", lw=1.3, alpha=0.96),
        zorder=5
    )

    ax_real.set_title("MUNDO REAL: Desfechos Epidemiológicos Observados (Ceará, 184 Municípios, Vizinhança Queen)", 
                      fontsize=12, fontweight="bold", loc="left", pad=16, color="#1C1C1C")

    # ==================== PAINEL SIMULAÇÕES (BOTTOM) ====================
    palette_box = ["#4A6984" if "Realista" in s else "#5C8374" for s in ordem_sim]
    
    # Distribuicao das replicacoes
    sns.boxenplot(
        data=df_sim,
        y="cenario_formatado",
        x="razao_D",
        hue="cenario_formatado",
        order=ordem_sim,
        ax=ax_sim,
        palette=palette_box,
        legend=False,
        box_kws=dict(alpha=0.38),
        k_depth=4,
        zorder=2
    )
    
    # Stripplot com transparência para exibir as 10 replicações exatas de cada cenário
    sns.stripplot(
        data=df_sim,
        y="cenario_formatado",
        x="razao_D",
        order=ordem_sim,
        ax=ax_sim,
        size=7.5,
        color="#1F3A52",
        alpha=0.75,
        jitter=0.22,
        edgecolor="white",
        linewidth=0.7,
        zorder=3
    )
    
    # Adicionar marcador de média de cada cenário
    for i, cenario in enumerate(ordem_sim):
        m = medias[cenario]
        ax_sim.plot(m, i, marker="o", markersize=9, color="#E65100", 
                    markeredgecolor="white", markeredgewidth=1.8, zorder=5)
        ax_sim.text(1.17, i, f"μ = {pt_num(m, 2)}", va="center", ha="left", 
                    fontsize=9.0, fontweight="bold", color="#2D3748")

    ax_sim.set_yticks(range(len(ordem_sim)))
    ax_sim.set_yticklabels(ordem_sim, fontsize=10)
    ax_sim.set_ylabel("")
    ax_sim.set_xlabel("Razão de Dispersão das Inclusões Posteriores:  D(SAR) / D(MQO)", 
                      fontsize=11.5, fontweight="bold", labelpad=12, color="#1A202C")
    
    ax_sim.set_xlim(0.0, 1.25)
    ax_sim.xaxis.set_major_locator(ticker.MultipleLocator(0.10))
    ax_sim.xaxis.set_major_formatter(ticker.FuncFormatter(lambda x, pos: pt_num(x, 1)))

    ax_sim.set_title("MUNDO SIMULADO: 8 Cenários de Monte Carlo (Distribuição Empírica de 10 Replicações Independentes)", 
                     fontsize=12, fontweight="bold", loc="left", pad=12, color="#1C1C1C")

    # Anotação de destaque no cenário mais severo (coordenadas seguras)
    ax_sim.annotate(
        "Cenário com ρ = 0,70 e 8 variáveis: \n"
        "Média da razão = 0,59 | Mínimo pontual = 0,22.\n"
        "A média e o mínimo descrevem aspectos distintos\n"
        "da distribuição das replicações.",
        xy=(0.588, 0.0), xytext=(0.08, 1.85),
        arrowprops=dict(arrowstyle="->", color="#D97706", lw=1.6, connectionstyle="arc3,rad=-0.18"),
        fontsize=8.8, fontweight="semibold", color="#78350F",
        bbox=dict(boxstyle="round,pad=0.45", fc="#FEF3C7", ec="#FDE68A", lw=1.3, alpha=0.96),
        zorder=6
    )

    # Rótulos nas zonas conceituais sombreadas
    ax_real.text(0.25, 2.78, "RAZÃO DE D < 0,50", ha="center", va="center",
                 fontsize=9.0, fontweight="heavy", color="#C0392B", alpha=0.9)
    ax_real.text(0.88, 2.78, "RAZÃO DE D > 0,55", ha="center", va="center",
                 fontsize=9.0, fontweight="heavy", color="#1F3A52", alpha=0.9)
    
    ax_real.text(1.005, -0.65, "Paridade: SAR = MQO", ha="left", va="bottom",
                 fontsize=8.5, fontstyle="italic", color="#4A5568", rotation=90)

    # Legenda customizada
    legend_elements = [
        mpatches.Patch(facecolor="#FDE8E8", edgecolor="#E53E3E", label="Razão de D < 0,50"),
        mpatches.Patch(facecolor="#EBF5FB", edgecolor="#2E5E7E", label="Razão de D > 0,55"),
        plt.Line2D([0], [0], marker="*", color="w", markerfacecolor="#C0392B", markersize=14, label="Dados observados: casos"),
        plt.Line2D([0], [0], marker="D", color="w", markerfacecolor="#2E5E7E", markersize=9, label="Desfechos Reais: Óbitos / Letalidade"),
        plt.Line2D([0], [0], marker="o", color="w", markerfacecolor="#1F3A52", markersize=8, label="Replicação Simulada de Monte Carlo"),
        plt.Line2D([0], [0], marker="o", color="w", markerfacecolor="#E65100", markersize=9, label="Média Empírica do Cenário"),
    ]
    fig.legend(handles=legend_elements, loc="upper center", bbox_to_anchor=(0.5, 0.985), 
               ncol=3, frameon=True, facecolor="#FAFAFA", edgecolor="#E2E8F0", fontsize=9.0)

    # Título Geral
    fig.suptitle("COMPARAÇÃO ENTRE DADOS OBSERVADOS E SIMULAÇÕES",
                 fontsize=15, fontweight="bold", y=1.025, color="#0F172A")
    
    plt.figtext(0.5, 0.995, 
                "Razão de D entre SAR e MQO nas três variáveis dependentes e nos cenários simulados",
                ha="center", fontsize=10.5, color="#475569", fontstyle="italic")

    # Rodapé Metodológico
    nota = (
        "Nota: D é a métrica de dispersão normalizada das PIPs: D = sd(PIP) * sqrt((K-1)/K) / sqrt(m*(1-m)). "
        "Razões abaixo de 0,50 indicam que D no SAR corresponde a menos da metade de D no MQO.\n"
        "Simulações baseadas em 10 replicações independentes por cenário com 120.000 iterações de MC³. "
        "Desfechos reais baseados em 1.000.000 de iterações em malha Queen de 184 municípios."
    )
    plt.figtext(0.08, -0.015, nota, ha="left", fontsize=8.2, color="#64748B")

    # Salvar nos dois formatos
    png_path = os.path.join(SAIDA_DIR, "fig_dois_mundos_executivo.png")
    pdf_path = os.path.join(SAIDA_DIR, "fig_dois_mundos_executivo.pdf")
    
    plt.savefig(png_path, dpi=300, bbox_inches="tight")
    plt.savefig(pdf_path, bbox_inches="tight")
    plt.close()
    
    print(f"  -> Gerado com sucesso: {png_path}")
    print(f"  -> Gerado com sucesso: {pdf_path}")


# ==============================================================================
# FIGURA 2: MECANISMO JOINTPLOT
# ==============================================================================
def gerar_figura_mecanismo():
    print("\n[2/2] Gerando 'fig_mecanismo_jointplot'...")
    
    # 1. Carregar dados do mecanismo
    p_mec = os.path.join(EXP_DIR, "mecanismo_completo.csv")
    df_mec = pd.read_csv(p_mec, encoding="utf-8")
    
    # 2. Carregar matriz de casos para calcular |r| com casos
    df_casos = pd.read_csv(CASOS_CSV)
    df_mec["r_casos"] = df_mec["covariavel"].map(lambda c: df_casos["z_log_casos100k"].corr(df_casos[c]))
    df_mec["abs_r"] = df_mec["r_casos"].abs()
    
    # Colormap contínuo
    cmap_cores = ["#2A9D8F", "#8AB17D", "#E9C46A", "#F4A261", "#E76F51", "#B4462F", "#800F2F"]
    cmap_mecanismo = LinearSegmentedColormap.from_list("queda_pip", cmap_cores, N=256)
    
    # Tamanho proporcional ao R² contra o resto das covariáveis
    min_r2, max_r2 = df_mec["r2_resto"].min(), df_mec["r2_resto"].max()
    tamanho_pontos = 90 + (df_mec["r2_resto"] - min_r2) / (max_r2 - min_r2) * 430
    
    # GridSpec
    fig = plt.figure(figsize=(11.5, 11.5), dpi=300)
    gs = fig.add_gridspec(2, 2, width_ratios=[8.2, 1.8], height_ratios=[1.8, 8.2],
                           wspace=0.04, hspace=0.04)
    
    ax_top = fig.add_subplot(gs[0, 0])
    ax_main = fig.add_subplot(gs[1, 0])
    ax_right = fig.add_subplot(gs[1, 1])
    ax_corner = fig.add_subplot(gs[0, 1])
    ax_corner.axis("off")
    
    # ==================== DISPERSÃO CENTRAL ====================
    ax_main.set_facecolor("#FCFCFD")
    ax_main.grid(True, linestyle=":", color="#E2E8F0", linewidth=0.8, zorder=0)
    
    # Linha de regressão (OLS)
    sns.regplot(
        data=df_mec,
        x="mi",
        y="abs_r",
        ax=ax_main,
        scatter=False,
        color="#64748B",
        line_kws=dict(linestyle="--", linewidth=1.7, alpha=0.85, zorder=1),
        ci=95,
        truncate=False
    )
    
    corr_val = df_mec["mi"].corr(df_mec["abs_r"])
    ax_main.text(
        0.02, 0.408, 
        f"Regressão OLS: |r| vs. I de Moran (r = {pt_num(corr_val, 2)}, p < 0,001)\n"
        "Correlação entre o I de Moran e |r| com o log da taxa de casos",
        fontsize=8.6, color="#334155", fontstyle="italic",
        bbox=dict(boxstyle="round,pad=0.35", fc="#FFFFFF", ec="#CBD5E0", lw=0.9, alpha=0.95),
        zorder=5
    )
    
    scatter = ax_main.scatter(
        df_mec["mi"],
        df_mec["abs_r"],
        s=tamanho_pontos,
        c=df_mec["queda"],
        cmap=cmap_mecanismo,
        edgecolor="#1E293B",
        linewidth=1.1,
        alpha=0.92,
        zorder=4
    )
    
    # Annotations for six selected covariates
    anotacoes = {
        "z_idosos": {
            "texto": "Idosos (z_idosos)\nI = 0,66 | |r| = 0,38 | R² = 0,56\nDiferença de PIP (SAR − MQO): −0,64",
            "xy": (0.6619, 0.3830),
            "xytext": (0.45, 0.405),
            "rad": -0.08,
            "cor": "#800F2F"
        },
        "semiarido": {
            "texto": "Semiárido (semiarido)\nI = 0,71 | |r| = 0,34 | R² = 0,45\nDiferença de PIP (SAR − MQO): −0,39",
            "xy": (0.7150, 0.3422),
            "xytext": (0.50, 0.285),
            "rad": -0.15,
            "cor": "#B4462F"
        },
        "z_leitos1k": {
            "texto": "Leitos SUS (z_leitos1k)\nI = 0,11 | |r| = 0,25 | R² = 0,56\nDiferença de PIP (SAR − MQO): −0,45",
            "xy": (0.1059, 0.2507),
            "xytext": (0.17, 0.265),
            "rad": -0.15,
            "cor": "#900C3F"
        },
        "z_ex.pobr": {
            "texto": "Extrema Pobreza (z_ex.pobr)\nI = 0,48 | |r| = 0,19 | R² = 0,71\nDiferença de PIP (SAR − MQO): −0,19",
            "xy": (0.4776, 0.1907),
            "xytext": (0.31, 0.125),
            "rad": -0.15,
            "cor": "#D97706"
        },
        "metrop": {
            "texto": "Região Metropolitana (metrop)\nI = 0,79 | |r| = 0,11 | R² = 0,61\nDiferença de PIP (SAR − MQO): −0,39",
            "xy": (0.7886, 0.1065),
            "xytext": (0.55, 0.170),
            "rad": 0.15,
            "cor": "#B4462F"
        },
        "z_tax.agua": {
            "texto": "Cobertura de Água (z_tax.agua)\nI = 0,00 | |r| = 0,12 | R² = 0,16\nDiferença de PIP (SAR − MQO): −0,12",
            "xy": (0.0002, 0.1227),
            "xytext": (0.02, 0.035),
            "rad": -0.15,
            "cor": "#2A9D8F"
        },
    }
    
    for cov_id, cfg in anotacoes.items():
        ax_main.annotate(
            cfg["texto"],
            xy=cfg["xy"],
            xytext=cfg["xytext"],
            arrowprops=dict(
                arrowstyle="->",
                color=cfg["cor"],
                lw=1.5,
                shrinkA=3,
                shrinkB=6,
                connectionstyle=f"arc3,rad={cfg['rad']}"
            ),
            fontsize=8.3,
            fontweight="bold",
            color="#1E293B",
            bbox=dict(boxstyle="round,pad=0.35", fc="#FFFFFF", ec=cfg["cor"], lw=1.3, alpha=0.96),
            zorder=6
        )

    ax_main.set_xlim(-0.04, 0.85)
    ax_main.set_ylim(-0.02, 0.44)
    ax_main.xaxis.set_major_locator(ticker.MultipleLocator(0.10))
    ax_main.xaxis.set_major_formatter(ticker.FuncFormatter(lambda x, pos: pt_num(x, 1)))
    ax_main.yaxis.set_major_locator(ticker.MultipleLocator(0.05))
    ax_main.yaxis.set_major_formatter(ticker.FuncFormatter(lambda y, pos: pt_num(y, 2)))
    
    ax_main.set_xlabel("I de Moran da Covariável (Autocorrelação Espacial / Suavidade)", 
                       fontsize=11.5, fontweight="bold", labelpad=10, color="#1A202C")
    ax_main.set_ylabel("Módulo da Correlação com o Desfecho:  |r(x, casos)|", 
                       fontsize=11.5, fontweight="bold", labelpad=10, color="#1A202C")

    # ==================== MARGINAIS (KDE) ====================
    ax_top.set_facecolor("#FCFCFD")
    sns.kdeplot(data=df_mec["mi"], ax=ax_top, color="#2E5E7E", fill=True, alpha=0.35, linewidth=1.5)
    sns.rugplot(data=df_mec["mi"], ax=ax_top, color="#2E5E7E", height=0.2, alpha=0.8)
    ax_top.set_xlim(ax_main.get_xlim())
    ax_top.set_ylabel("Densidade", fontsize=9.0, color="#4A5568")
    ax_top.set_xlabel("")
    ax_top.tick_params(labelbottom=False, labelsize=8.5)
    ax_top.spines["top"].set_visible(False)
    ax_top.spines["right"].set_visible(False)
    ax_top.grid(axis="x", linestyle=":", color="#E2E8F0", linewidth=0.7)

    ax_right.set_facecolor("#FCFCFD")
    sns.kdeplot(y=df_mec["abs_r"], ax=ax_right, color="#B4462F", fill=True, alpha=0.35, linewidth=1.5)
    sns.rugplot(y=df_mec["abs_r"], ax=ax_right, color="#B4462F", height=0.2, alpha=0.8)
    ax_right.set_ylim(ax_main.get_ylim())
    ax_right.set_xlabel("Densidade", fontsize=9.0, color="#4A5568")
    ax_right.set_ylabel("")
    ax_right.tick_params(labelleft=False, labelsize=8.5)
    ax_right.spines["top"].set_visible(False)
    ax_right.spines["right"].set_visible(False)
    ax_right.grid(axis="y", linestyle=":", color="#E2E8F0", linewidth=0.7)

    # ==================== LEGENDAS E COLORBAR ====================
    # Horizontal colorbar in the upper-left panel
    # Inset axes dentro de ax_main
    from mpl_toolkits.axes_grid1.inset_locator import inset_axes
    cax = inset_axes(ax_main, width="32%", height="3.2%", loc="upper left",
                     bbox_to_anchor=(0.04, -0.15, 1, 1), bbox_transform=ax_main.transAxes)
    cbar = fig.colorbar(scatter, cax=cax, orientation="horizontal")
    cbar.set_label("Queda de PIP:  ΔPIP = PIP(MQO) − PIP(SAR)", fontsize=8.5, fontweight="bold", color="#1F2937", labelpad=5)
    cbar.ax.xaxis.set_major_formatter(ticker.FuncFormatter(lambda v, p: pt_num(v, 2)))
    cbar.ax.tick_params(labelsize=8)

    # Legenda de tamanhos (R² contra as 29 covariáveis) posicionada no quadrante inferior direito livre
    r2_vals = [0.20, 0.50, 0.80]
    legend_sizes = [90 + (rv - min_r2) / (max_r2 - min_r2) * 430 for rv in r2_vals]
    size_handles = [
        plt.scatter([], [], s=sz, c="#E2E8F0", edgecolors="#1E293B", linewidth=1.1, alpha=0.9)
        for sz in legend_sizes
    ]
    labels_sz = [
        f"R² = {pt_num(r2_vals[0], 2)} (baixa colinearidade)", 
        f"R² = {pt_num(r2_vals[1], 2)} (moderada)", 
        f"R² = {pt_num(r2_vals[2], 2)} (alta colinearidade)"
    ]
    leg = ax_main.legend(
        size_handles, labels_sz,
        title="Colinearidade Interna (R² contra 29 covariáveis)",
        loc="lower right",
        bbox_to_anchor=(0.98, 0.02),
        frameon=True,
        facecolor="#FFFFFF",
        edgecolor="#CBD5E0",
        fontsize=8.3,
        title_fontsize=8.8
    )
    leg.get_title().set_fontweight("bold")

    fig.suptitle("ASSOCIAÇÕES ENTRE COLINEARIDADE, MORAN E MUDANÇAS DE PIP",
                 fontsize=15, fontweight="bold", y=0.965, color="#0F172A")
    
    plt.figtext(
        0.5, 0.940,
        "Comparação exploratória entre autocorrelação espacial, correlação com casos, R² interno e diferença de PIP",
        ha="center", fontsize=10.5, color="#475569", fontstyle="italic"
    )

    nota_mec = (
        "Nota: I de Moran calculado na malha Queen (184 municípios). |r| é a correlação linear bivariada com o log da taxa de casos padronizada. "
        "R² interno obtido pela regressão de cada covariável contra todas as demais 29 candidatas.\n"
        "A cor representa PIP(MQO) − PIP(SAR); as anotações mostram a diferença inversa, PIP(SAR) − PIP(MQO)."
    )
    plt.figtext(0.12, 0.035, nota_mec, ha="left", fontsize=8.2, color="#64748B")

    png_path = os.path.join(SAIDA_DIR, "fig_mecanismo_jointplot.png")
    pdf_path = os.path.join(SAIDA_DIR, "fig_mecanismo_jointplot.pdf")
    
    plt.savefig(png_path, dpi=300, bbox_inches="tight")
    plt.savefig(pdf_path, bbox_inches="tight")
    plt.close()
    
    print(f"  -> Gerado com sucesso: {png_path}")
    print(f"  -> Gerado com sucesso: {pdf_path}")


# ==============================================================================
# EXECUÇÃO PRINCIPAL
# ==============================================================================
if __name__ == "__main__":
    print("=" * 70)
    print("Iniciando geração de figuras científicas para o artigo do MC³ espacial...")
    print(f"Diretório de saída: {SAIDA_DIR}")
    print("=" * 70)
    
    gerar_figura_dois_mundos()
    gerar_figura_mecanismo()
    
    print("\n" + "=" * 70)
    print("Todas as figuras foram geradas com sucesso em PNG (300 DPI) e PDF vetorial!")
    print("=" * 70)
