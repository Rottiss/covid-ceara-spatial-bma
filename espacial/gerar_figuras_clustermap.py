# -*- coding: utf-8 -*-
"""
gerar_figuras_clustermap.py
---------------------------
Graficos de PIPs por covariavel e especificacao:
1. fig_clustermap_pips: Clustermap hierárquico das PIPs (MQO, SAR, SEM para Casos, Óbitos e Letalidade).
2. fig_heatmap_delta: Heatmap divergente de transição e mudança líquida de PIPs (ΔSAR e ΔSEM).
"""

import os
from pathlib import Path
import sys
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns

if hasattr(sys.stdout, 'reconfigure'):
    sys.stdout.reconfigure(encoding='utf-8')

plt.style.use('seaborn-v0_8-white' if 'seaborn-v0_8-white' in plt.style.available else 'default')
mpl.rcParams['font.sans-serif'] = ['DejaVu Sans', 'Arial', 'Helvetica']
mpl.rcParams['font.size'] = 9.5
mpl.rcParams['axes.edgecolor'] = '#333333'
mpl.rcParams['axes.linewidth'] = 0.8
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42

DIR_BASE = str(Path(__file__).resolve().parent)
ARQ_CSV = os.path.join(DIR_BASE, "saida", "export_python", "pips_completas.csv")
DIR_SAIDA = os.path.join(DIR_BASE, "saida", "figuras_novas")
os.makedirs(DIR_SAIDA, exist_ok=True)

ROTULOS_PT = {
    "z_pop": "População",
    "z_area": "Área territorial",
    "z_densidade": "Densidade demográfica",
    "z_idosos": "Proporção de idosos (60+)",
    "z_pop.rural": "População rural",
    "z_pibpc": "PIB per capita",
    "z_desp.educ": "Despesa com educação",
    "z_ideb.fund": "IDEB fundamental",
    "z_ideb3": "IDEB ensino médio",
    "z_bolsaf": "Bolsa Família (cobertura)",
    "z_ex.pobr": "Extrema pobreza",
    "z_emprego": "Emprego formal",
    "z_ivs.infra": "IVS infraestrutura",
    "z_tax.agua": "Cobertura de água",
    "z_energia100k": "Consumo de energia",
    "z_tax.hom": "Taxa de homicídios",
    "z_sus1k": "Unidades do SUS",
    "z_leitos1k": "Leitos hospitalares SUS",
    "z_prof1k": "Profissionais de saúde",
    "z_desp.saude": "Despesa com saúde",
    "z_imuni": "Cobertura de imunização",
    "z_int.circ": "Internações circulatórias",
    "z_int.resp": "Internações respiratórias",
    "z_int.diab": "Internações por diabetes",
    "z_int.asma": "Internações por asma",
    "metrop": "Região Metropolitana (RMF)",
    "semiarido": "Semiárido",
    "aliadogov": "Prefeito aliado ao governo",
    "z_votos": "Votação do prefeito",
    "z_eleitas.fem": "Mulheres eleitas",
}

GRUPOS = {
    "z_pop": "Demografia",
    "z_area": "Demografia",
    "z_densidade": "Demografia",
    "z_idosos": "Demografia",
    "z_pop.rural": "Demografia",
    "z_pibpc": "Vulnerabilidade e Renda",
    "z_desp.educ": "Vulnerabilidade e Renda",
    "z_ideb.fund": "Vulnerabilidade e Renda",
    "z_ideb3": "Vulnerabilidade e Renda",
    "z_bolsaf": "Vulnerabilidade e Renda",
    "z_ex.pobr": "Vulnerabilidade e Renda",
    "z_emprego": "Vulnerabilidade e Renda",
    "z_ivs.infra": "Vulnerabilidade e Renda",
    "z_tax.agua": "Saúde e Saneamento",
    "z_energia100k": "Saúde e Saneamento",
    "z_tax.hom": "Saúde e Saneamento",
    "z_sus1k": "Saúde e Saneamento",
    "z_leitos1k": "Saúde e Saneamento",
    "z_prof1k": "Saúde e Saneamento",
    "z_desp.saude": "Saúde e Saneamento",
    "z_imuni": "Saúde e Saneamento",
    "z_int.circ": "Saúde e Saneamento",
    "z_int.resp": "Saúde e Saneamento",
    "z_int.diab": "Saúde e Saneamento",
    "z_int.asma": "Saúde e Saneamento",
    "metrop": "Território e Política",
    "semiarido": "Território e Política",
    "aliadogov": "Território e Política",
    "z_votos": "Território e Política",
    "z_eleitas.fem": "Território e Política",
}

CORES_GRUPOS = {
    "Demografia": "#2B5C8F",              # Azul ardósia
    "Vulnerabilidade e Renda": "#D96B27", # Terracota/âmbar quente
    "Saúde e Saneamento": "#2A9D8F",      # Green
    "Território e Política": "#7D53A2",   # Violeta nobre
}

def carregar_dados():
    df = pd.read_csv(ARQ_CSV, encoding='latin1')
    piv = df.pivot(index='covariavel', columns='desfecho', values=['pip_mqo', 'pip_sar_queen', 'pip_sem_queen'])
    
    col_order = [
        ('pip_mqo', 'casos'), ('pip_sar_queen', 'casos'), ('pip_sem_queen', 'casos'),
        ('pip_mqo', 'obitos'), ('pip_sar_queen', 'obitos'), ('pip_sem_queen', 'obitos'),
        ('pip_mqo', 'letalidade'), ('pip_sar_queen', 'letalidade'), ('pip_sem_queen', 'letalidade'),
    ]
    
    col_names = ['MQO', 'SAR', 'SEM', 'MQO', 'SAR', 'SEM', 'MQO', 'SAR', 'SEM']
    
    mat_pips = piv[col_order].copy()
    mat_pips.columns = col_names
    
    df_c = df[df['desfecho'] == 'casos'].set_index('covariavel')
    df_o = df[df['desfecho'] == 'obitos'].set_index('covariavel')
    df_l = df[df['desfecho'] == 'letalidade'].set_index('covariavel')
    
    mat_deltas = pd.DataFrame(index=mat_pips.index)
    mat_deltas['ΔSAR Casos'] = df_c['pip_sar_queen'] - df_c['pip_mqo']
    mat_deltas['ΔSEM Casos'] = df_c['pip_sem_queen'] - df_c['pip_mqo']
    mat_deltas['ΔSAR Óbitos'] = df_o['pip_sar_queen'] - df_o['pip_mqo']
    mat_deltas['ΔSEM Óbitos'] = df_o['pip_sem_queen'] - df_o['pip_mqo']
    mat_deltas['ΔSAR Letalidade'] = df_l['pip_sar_queen'] - df_l['pip_mqo']
    mat_deltas['ΔSEM Letalidade'] = df_l['pip_sem_queen'] - df_l['pip_mqo']
    
    return df, mat_pips, mat_deltas

def gerar_clustermap_pips(mat_pips):
    """
    Figura 1: Clustermap hierárquico das PIPs.
    """
    print("Gerando Figura 1: Clustermap Hierárquico de PIPs...")
    
    grupos_series = pd.Series([GRUPOS[c] for c in mat_pips.index], index=mat_pips.index)
    row_colors = grupos_series.map(CORES_GRUPOS)
    row_colors.name = ""
    
    mat_plot = mat_pips.copy()
    rotulos_indices = [ROTULOS_PT[c] for c in mat_plot.index]
    mat_plot.index = rotulos_indices
    row_colors.index = rotulos_indices
    
    # Anotações: apenas células com PIP >= 0,50 recebem texto
    annot_mat = mat_plot.map(lambda v: f"{v:.2f}".replace('.', ',') if v >= 0.50 else "")
    
    cmap = sns.color_palette("crest", as_cmap=True)
    
    g = sns.clustermap(
        mat_plot,
        method="ward",
        metric="euclidean",
        row_cluster=True,
        col_cluster=False,
        row_colors=row_colors,
        cmap=cmap,
        vmin=0.0,
        vmax=1.0,
        annot=annot_mat,
        fmt="",
        annot_kws={"size": 8.5, "va": "center", "ha": "center"},
        figsize=(12.8, 12.0),
        dendrogram_ratio=(0.13, 0.02),
        colors_ratio=0.024,
        cbar_pos=(0.025, 0.88, 0.16, 0.022),
        cbar_kws={"orientation": "horizontal"},
        linewidths=0.6,
        linecolor="#EAEAEA",
    )
    
    ax = g.ax_heatmap
    
    # Contraste de texto fino
    for text in ax.texts:
        txt = text.get_text().replace(',', '.')
        if txt != "":
            try:
                val = float(txt)
                if val >= 0.70:
                    text.set_color('#FFFFFF')
                    text.set_weight('bold')
                else:
                    text.set_color('#062436')
                    text.set_weight('bold')
            except ValueError:
                pass
    
    g.ax_row_colors.set_xlabel("")
    g.ax_row_colors.set_xticks([])
    
    # Separadores verticais entre desfechos
    ax.axvline(3, color='#FFFFFF', linewidth=2.8, zorder=5)
    ax.axvline(6, color='#FFFFFF', linewidth=2.8, zorder=5)
    
    # Rótulos dos eixos
    ax.set_xticklabels(ax.get_xticklabels(), rotation=0, fontsize=10.0, fontweight='bold', color='#222222')
    ax.set_yticklabels(ax.get_yticklabels(), rotation=0, fontsize=9.2, color='#1A1A1A')
    ax.set_xlabel("")
    ax.set_ylabel("")
    
    # Cabeçalhos temáticos
    ax.text(1.5, -0.65, "CASOS", ha="center", va="bottom", fontsize=11.5, fontweight="bold", color="#1B3B6F",
            bbox=dict(boxstyle="round,pad=0.28", facecolor="#EAF2F8", edgecolor="#B0C4DE", linewidth=1.2))
    ax.text(4.5, -0.65, "ÓBITOS", ha="center", va="bottom", fontsize=11.5, fontweight="bold", color="#6E2C00",
            bbox=dict(boxstyle="round,pad=0.28", facecolor="#FBEEE6", edgecolor="#EDBB99", linewidth=1.2))
    ax.text(7.5, -0.65, "LETALIDADE", ha="center", va="bottom", fontsize=11.5, fontweight="bold", color="#145A32",
            bbox=dict(boxstyle="round,pad=0.28", facecolor="#E8F8F5", edgecolor="#A3E4D7", linewidth=1.2))
    
    # Colorbar Horizontal
    cbar = g.ax_cbar
    cbar.set_title("PIP (Probabilidade a Posteriori)", fontsize=9.0, fontweight="bold", pad=4)
    cbar.set_xticks([0.0, 0.50, 1.0])
    cbar.set_xticklabels(["0,0", "0,5", "1,0"], fontsize=8.5)
    cbar.tick_params(top=False, bottom=True, labeltop=False, labelbottom=True, size=3.5, pad=3)
    cbar.axvline(0.50, color="#D90429", linestyle="--", linewidth=1.6)
    cbar.text(0.50, -1.9, "Limiar ≥ 0,50", transform=cbar.get_xaxis_transform(),
              fontsize=7.8, color="#D90429", fontweight="bold", ha="center", va="top")
    
    # Legenda dos Grupos Temáticos
    patches_legenda = [
        mpatches.Patch(facecolor=cor, edgecolor="#444444", linewidth=0.5, label=grp)
        for grp, cor in CORES_GRUPOS.items()
    ]
    g.fig.legend(
        handles=patches_legenda,
        title="Grupo Temático",
        title_fontsize=9.0,
        fontsize=8.2,
        loc="upper left",
        bbox_to_anchor=(0.024, 0.81),
        frameon=True,
        framealpha=0.96,
        edgecolor="#D0D0D0",
        facecolor="#FAFAFA"
    )
    
    # Título principal e nota de rodapé
    g.fig.suptitle(
        "Clustermap Hierárquico das Probabilidades de Inclusão a Posteriori (PIPs)",
        fontsize=14.5, fontweight="bold", y=0.992, color="#111111"
    )
    g.fig.text(
        0.53, 0.015,
        "Nota: Rótulos numéricos exibidos exclusivamente nas células com PIP ≥ 0,50 (Modelo de Probabilidade Mediana de Barbieri & Berger).\n"
        "O dendrograma hierárquico nas linhas (ligação de Ward, distância Euclidiana) evidencia blocos de covariáveis com trajetória similar de resposta espacial.",
        ha="center", fontsize=8.5, color="#4F4F4F", style="italic"
    )
    
    plt.subplots_adjust(top=0.92, bottom=0.07, left=0.18, right=0.81)
    
    out_png = os.path.join(DIR_SAIDA, "fig_clustermap_pips.png")
    out_pdf = os.path.join(DIR_SAIDA, "fig_clustermap_pips.pdf")
    g.savefig(out_png, dpi=300, bbox_inches="tight")
    g.savefig(out_pdf, format="pdf", bbox_inches="tight")
    plt.close()
    
    print(f"  -> Concluído: {out_png}")
    print(f"  -> Concluído: {out_pdf}")
    
    row_order_idx = g.dendrogram_row.reordered_ind
    ordered_covs = [mat_pips.index[i] for i in row_order_idx]
    return ordered_covs

def gerar_heatmap_delta(mat_deltas, ordered_covs):
    """
    Figura 2: Heatmap de Mudança Líquida (Delta PIP = SAR - MQO e SEM - MQO).
    """
    print("Gerando Figura 2: Heatmap de Mudança Líquida (Delta PIP)...")
    
    mat_deltas_ord = mat_deltas.loc[ordered_covs].copy()
    rotulos_indices = [ROTULOS_PT[c] for c in mat_deltas_ord.index]
    mat_deltas_ord.index = rotulos_indices
    
    col_labels = [
        "ΔSAR\n(SAR − MQO)", "ΔSEM\n(SEM − MQO)",
        "ΔSAR\n(SAR − MQO)", "ΔSEM\n(SEM − MQO)",
        "ΔSAR\n(SAR − MQO)", "ΔSEM\n(SEM − MQO)"
    ]
    mat_deltas_ord.columns = col_labels
    
    def formatar_delta(v):
        if abs(v) < 0.005:
            return "0,00"
        return f"{v:+.2f}".replace('.', ',')
        
    annot_mat = mat_deltas_ord.map(formatar_delta)
    
    fig = plt.figure(figsize=(12.0, 12.0))
    gs = fig.add_gridspec(1, 2, width_ratios=[0.024, 0.976], wspace=0.015,
                           left=0.06, right=0.74, top=0.89, bottom=0.08)
    
    ax_bar = fig.add_subplot(gs[0, 0])
    ax_heat = fig.add_subplot(gs[0, 1])
    
    lim_max = 0.85
    cmap = sns.color_palette("vlag", as_cmap=True)
    
    # 1. Barra lateral de grupos temáticos
    grupos_ord = [GRUPOS[c] for c in ordered_covs]
    grp_unique = ["Demografia", "Vulnerabilidade e Renda", "Saúde e Saneamento", "Território e Política"]
    grp_map_idx = {g: i for i, g in enumerate(grp_unique)}
    color_matrix = np.array([[grp_map_idx[g]] for g in grupos_ord])
    
    cmap_grp = mpl.colors.ListedColormap([CORES_GRUPOS[g] for g in grp_unique])
    bounds = np.arange(len(grp_unique) + 1) - 0.5
    norm_grp = mpl.colors.BoundaryNorm(bounds, cmap_grp.N)
    
    ax_bar.imshow(color_matrix, aspect='auto', cmap=cmap_grp, norm=norm_grp)
    ax_bar.set_xticks([])
    ax_bar.set_yticks([])
    for sp in ax_bar.spines.values():
        sp.set_visible(False)
        
    # 2. Heatmap de Deltas
    cbar_ax = fig.add_axes([0.94, 0.35, 0.018, 0.30])
    
    sns.heatmap(
        mat_deltas_ord,
        ax=ax_heat,
        cmap=cmap,
        center=0.0,
        vmin=-lim_max,
        vmax=lim_max,
        annot=annot_mat,
        fmt="",
        annot_kws={"size": 8.2, "va": "center", "ha": "center"},
        linewidths=0.6,
        linecolor="#EAEAEA",
        cbar_ax=cbar_ax,
        cbar_kws={
            "ticks": [-0.80, -0.60, -0.40, -0.20, 0.00, 0.20, 0.40],
        }
    )
    
    for text in ax_heat.texts:
        txt = text.get_text().replace(',', '.')
        try:
            val = float(txt)
            if abs(val) >= 0.45:
                text.set_color('#FFFFFF')
                text.set_weight('bold')
            elif abs(val) >= 0.12:
                text.set_color('#111111')
                text.set_weight('semibold')
            else:
                text.set_color('#555555')
                text.set_weight('normal')
        except ValueError:
            pass
            
    # Linhas divisórias brancas entre blocos de desfecho
    ax_heat.axvline(2, color='#FFFFFF', linewidth=2.8, zorder=5)
    ax_heat.axvline(4, color='#FFFFFF', linewidth=2.8, zorder=5)
    
    # Rótulos no eixo Y à direita
    ax_heat.yaxis.tick_right()
    ax_heat.set_yticklabels(mat_deltas_ord.index, rotation=0, fontsize=9.2, color='#1A1A1A')
    ax_heat.tick_params(axis='y', length=0, pad=8)
    
    # Rótulos no eixo X
    ax_heat.set_xticklabels(ax_heat.get_xticklabels(), rotation=0, fontsize=9.0, fontweight='semibold', color='#222222')
    ax_heat.set_ylabel("")
    ax_heat.set_xlabel("")
    
    # Títulos superiores dos desfechos
    ax_heat.text(1.0, -0.72, "CASOS", ha="center", va="bottom",
                 fontsize=10.5, fontweight="bold", color="#1B3B6F",
                 bbox=dict(boxstyle="round,pad=0.3", facecolor="#EAF2F8", edgecolor="#B0C4DE", linewidth=1.2))
    ax_heat.text(3.0, -0.72, "ÓBITOS", ha="center", va="bottom",
                 fontsize=10.5, fontweight="bold", color="#6E2C00",
                 bbox=dict(boxstyle="round,pad=0.3", facecolor="#FBEEE6", edgecolor="#EDBB99", linewidth=1.2))
    ax_heat.text(5.0, -0.72, "LETALIDADE", ha="center", va="bottom",
                 fontsize=10.5, fontweight="bold", color="#145A32",
                 bbox=dict(boxstyle="round,pad=0.3", facecolor="#E8F8F5", edgecolor="#A3E4D7", linewidth=1.2))
    
    # Barra de cores
    cbar_ax.set_title("Δ PIP\n", fontsize=9.0, fontweight="bold", pad=4)
    yticklabels = [f"{t:+.2f}".replace('.', ',') for t in [-0.80, -0.60, -0.40, -0.20, 0.00, 0.20, 0.40]]
    cbar_ax.set_yticklabels(yticklabels, fontsize=8.0)
    cbar_ax.axhline(0.00, color="#111111", linewidth=1.2)
    
    # Textos explicativos na colorbar
    cbar_ax.text(0.5, 1.14, "Aumento de PIP", ha="center", va="bottom",
                transform=cbar_ax.transAxes, fontsize=7.8, color="#B2182B", fontweight="bold")
    cbar_ax.text(0.5, -0.14, "Redução de PIP", ha="center", va="top",
                transform=cbar_ax.transAxes, fontsize=7.8, color="#2166AC", fontweight="bold")
    
    # Legenda dos Grupos Temáticos na base
    patches_legenda = [
        mpatches.Patch(facecolor=CORES_GRUPOS[grp], edgecolor="#333333", linewidth=0.5, label=grp)
        for grp in grp_unique
    ]
    fig.legend(
        handles=patches_legenda,
        title="Grupo Temático",
        title_fontsize=9.0,
        fontsize=8.5,
        loc="lower center",
        bbox_to_anchor=(0.42, 0.016),
        ncol=4,
        frameon=True,
        framealpha=0.96,
        edgecolor="#D0D0D0",
        facecolor="#FAFAFA"
    )
    
    fig.suptitle(
        "Matriz de Transição e Mudança Líquida de PIP (Δ PIP = Modelo Espacial − MQO)",
        fontsize=14.5, fontweight="bold", y=0.988, color="#111111"
    )
    fig.text(
        0.50, 0.002,
        "Nota: ΔSAR = PIP(SAR Queen) − PIP(MQO); ΔSEM = PIP(SEM Queen) − PIP(MQO). Vermelho indica aumento de PIP; azul indica redução.\n"
        "As mudanças de PIP são condicionais à especificação. As covariáveis seguem a ordem do agrupamento hierárquico.",
        ha="center", fontsize=8.5, color="#4F4F4F", style="italic"
    )
    
    out_png = os.path.join(DIR_SAIDA, "fig_heatmap_delta.png")
    out_pdf = os.path.join(DIR_SAIDA, "fig_heatmap_delta.pdf")
    fig.savefig(out_png, dpi=300, bbox_inches="tight")
    fig.savefig(out_pdf, format="pdf", bbox_inches="tight")
    plt.close()
    
    print(f"  -> Concluído: {out_png}")
    print(f"  -> Concluído: {out_pdf}")

def main():
    print("=" * 65)
    print("Iniciando geração de figuras científicas para MC³ Espacial...")
    print(f"Ambiente: Python {sys.version.split()[0]} | Seaborn {sns.__version__}")
    print(f"Diretório de saída: {DIR_SAIDA}")
    print("=" * 65)
    
    df, mat_pips, mat_deltas = carregar_dados()
    
    ordered_covs = gerar_clustermap_pips(mat_pips)
    gerar_heatmap_delta(mat_deltas, ordered_covs)
    
    print("=" * 65)
    print("Todas as figuras foram geradas com sucesso em PNG (300 DPI) e PDF vetorial!")
    print("=" * 65)

if __name__ == "__main__":
    main()
