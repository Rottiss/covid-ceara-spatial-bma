"""
=============================================================================
gerar_figuras_severidade_rho.py
=============================================================================
Densidades posteriores de rho, estatísticas de Moran e sensibilidade às prioris.

Gera:
1. fig_densidades_rho_moran.png e fig_densidades_rho_moran.pdf:
   - Curvas de densidade posterior refinadas de rho para Casos, Óbitos e Letalidade
     (corrigindo os pesos de quadratura de Simpson para devolver a densidade real suave).
   - Preenchimento sob a curva em gradiente translúcido (alpha gradient).
   - Linha vertical marcando a média posterior de rho e o intervalo HPD de 95% sombreado.
   - Linha pontilhada contrastante marcando o I de Moran observado de y para cada desfecho.
   - Comparação entre Moran, média posterior de rho e razão de D entre SAR e MQO.
2. fig_bayes_factor_sensibilidade.png e fig_bayes_factor_sensibilidade.pdf:
   - Painéis com log do Fator de Bayes e probabilidade posterior P(SAR | y) em função de a0.
   - Destaque claro para a estabilidade de Casos e Óbitos (sempre favoráveis ao SAR)
     contrastando com a sensibilidade de Letalidade.
=============================================================================
"""

import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.patches import Polygon, Rectangle, FancyBboxPatch
from matplotlib.lines import Line2D
import matplotlib.gridspec as gridspec
from scipy.interpolate import interp1d

# Diretórios do projeto
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
DADOS_DIR = os.path.join(BASE_DIR, "saida")
EXPORT_DIR = os.path.join(DADOS_DIR, "export_python")
SAIDA_FIGURAS = os.path.join(DADOS_DIR, "figuras_novas")
os.makedirs(SAIDA_FIGURAS, exist_ok=True)

# ---------------------------------------------------------------------------
# Cores e tipografia
# ---------------------------------------------------------------------------
CORES = {
    "casos": "#B4462F",        # Terracota
    "obitos": "#C98A3C",       # Amber
    "letalidade": "#2E5E7E",   # Azul petróleo / Ardósia (Frio)
    "fundo": "#FBFAF8",        # Papel editorial suave
    "tinta": "#1C1C1C",        # Grafite escuro / Texto principal
    "cinza": "#5A5A5A",        # Texto secundário
    "grade": "#EAE7E2",        # Linhas de grade suaves
    "destaque": "#0F52BA",     # Azul safira
    "moran": "#3E3E3E",        # Grafite neutro para Moran
    "alerta": "#C0392B",       # Red
    "seguro": "#27AE60"        # Verde para estabilidade
}

ROTULOS = {
    "casos": "Casos Confirmados",
    "obitos": "Óbitos Acumulados",
    "letalidade": "Taxa de Letalidade"
}

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Segoe UI", "Arial", "Helvetica", "Calibri", "DejaVu Sans"],
    "text.color": CORES["tinta"],
    "axes.labelcolor": CORES["tinta"],
    "xtick.color": CORES["tinta"],
    "ytick.color": CORES["tinta"],
    "figure.facecolor": CORES["fundo"],
    "axes.facecolor": CORES["fundo"],
    "axes.edgecolor": "#D0CDC7",
    "axes.linewidth": 0.8,
    "grid.color": CORES["grade"],
    "grid.linewidth": 0.6,
    "pdf.fonttype": 42,
    "ps.fonttype": 42
})


# ---------------------------------------------------------------------------
# Funções Auxiliares: Desempacotamento de Simpson e Estatísticas Posteriores
# ---------------------------------------------------------------------------
def desempacotar_densidade_simpson(sub_df):
    """
    Remove os pesos de quadratura de Simpson embutidos na exportação
    e normaliza a densidade para que a integral sobre rho seja 1.
    """
    sub = sub_df.sort_values("rho").reset_index(drop=True)
    ng = len(sub)
    pesos = np.ones(ng)
    pesos[1:-1:2] = 4.0
    pesos[2:-1:2] = 2.0
    
    rho_grid = sub["rho"].values
    dens_bruta = sub["densidade"].values
    dens_unweighted = dens_bruta / pesos
    
    # Integral numérica pela regra do trapézio
    integral = np.trapezoid(dens_unweighted, rho_grid)
    dens_norm = dens_unweighted / integral
    
    # Interpolação fina para renderização ultrassuave (1000 pontos)
    rho_fino = np.linspace(rho_grid.min(), rho_grid.max(), 1000)
    f_interp = interp1d(rho_grid, dens_norm, kind="cubic", fill_value="extrapolate")
    dens_fina = np.maximum(0.0, f_interp(rho_fino))
    integral_fina = np.trapezoid(dens_fina, rho_fino)
    dens_fina = dens_fina / integral_fina
    
    # Média posterior exata
    media_rho = np.trapezoid(rho_fino * dens_fina, rho_fino)
    
    # Modo posterior
    modo_rho = rho_fino[np.argmax(dens_fina)]
    
    # Intervalo HPD de 95%
    h = rho_fino[1] - rho_fino[0]
    idx_sorted = np.argsort(dens_fina)[::-1]
    cum_mass = np.cumsum(dens_fina[idx_sorted] * h)
    idx_cutoff = np.searchsorted(cum_mass, 0.95)
    cutoff = dens_fina[idx_sorted[min(idx_cutoff, len(dens_fina) - 1)]]
    in_hpd = dens_fina >= cutoff
    hpd_inf = rho_fino[in_hpd].min()
    hpd_sup = rho_fino[in_hpd].max()
    
    moran = sub["moran_y"].iloc[0]
    
    return {
        "rho": rho_fino,
        "densidade": dens_fina,
        "media": media_rho,
        "modo": modo_rho,
        "hpd_inf": hpd_inf,
        "hpd_sup": hpd_sup,
        "moran": moran,
        "cutoff": cutoff
    }


def criar_gradiente_vertical(ax, x, y, cor_hex, y_max_ref, alpha_min=0.04, alpha_max=0.52, n_passos=120):
    """
    Cria um gradiente vertical suave translúcido preenchendo a área sob a curva y = f(x).
    """
    r, g, b = plt.matplotlib.colors.to_rgb(cor_hex)
    x_min, x_max = x.min(), x.max()
    y_max = max(y.max(), y_max_ref)
    
    # Matriz RGBA de gradiente
    grad_array = np.zeros((n_passos, 1, 4))
    grad_array[:, 0, 0] = r
    grad_array[:, 0, 1] = g
    grad_array[:, 0, 2] = b
    grad_array[:, 0, 3] = np.linspace(alpha_min, alpha_max, n_passos)
    
    # Imagem de gradiente no plano de dados
    im = ax.imshow(
        grad_array,
        aspect="auto",
        extent=[x_min, x_max, 0, y_max],
        origin="lower",
        zorder=2
    )
    
    # Polígono para clipar a imagem exatamente sob a curva
    poly_pts = np.column_stack([
        np.concatenate(([x_min], x, [x_max])),
        np.concatenate(([0.0], y, [0.0]))
    ])
    clip_polygon = Polygon(poly_pts, closed=True, transform=ax.transData)
    im.set_clip_path(clip_polygon)
    
    return im


# ---------------------------------------------------------------------------
# FIGURA 1: DENSIDADES DE RHO, MORAN E RAZÕES DE D
# ---------------------------------------------------------------------------
def gerar_figura_1():
    print("Gerando Figura 1: Densidades Posteriores de Rho e Razões de D...")
    
    # 1. Carregar dados
    caminho_dens = os.path.join(EXPORT_DIR, "densidades_rho.csv")
    df_dens = pd.read_csv(caminho_dens)
    
    # Carregar dados de severidade para obter D_mqo, D_sar e razão
    caminho_sev = os.path.join(DADOS_DIR, "severidade_rho.csv")
    if os.path.exists(caminho_sev):
        df_sev = pd.read_csv(caminho_sev).set_index("desfecho")
    else:
        df_sev = pd.DataFrame({
            "mi_y": [0.4834, 0.4414, 0.1384],
            "rho": [0.7079, 0.5549, 0.2544],
            "D_mqo": [0.3235, 0.4652, 0.2000],
            "D_sar": [0.0973, 0.3299, 0.1413],
            "razao": [0.3007, 0.7093, 0.7064]
        }, index=["casos", "obitos", "letalidade"])
    
    desfechos = ["casos", "obitos", "letalidade"]
    stats = {}
    for d in desfechos:
        sub = df_dens[df_dens["desfecho"] == d]
        stats[d] = desempacotar_densidade_simpson(sub)
    
    # 2. Configurar Grid de Subplots
    fig = plt.figure(figsize=(15.8, 9.6), dpi=300)
    
    # Espaco entre cabecalho e paineis
    gs_master = gridspec.GridSpec(1, 2, width_ratios=[1.35, 1.0], wspace=0.30,
                                  left=0.07, right=0.96, top=0.845, bottom=0.085)
    
    # Coluna Esquerda: 3 densidades empilhadas
    gs_left = gridspec.GridSpecFromSubplotSpec(3, 1, subplot_spec=gs_master[0], hspace=0.36)
    axes_dens = [fig.add_subplot(gs_left[i]) for i in range(3)]
    
    # Coluna Direita: 2 diagnósticos complementares da tese
    gs_right = gridspec.GridSpecFromSubplotSpec(2, 1, subplot_spec=gs_master[1], hspace=0.38)
    ax_comp = fig.add_subplot(gs_right[0])
    ax_sev = fig.add_subplot(gs_right[1])
    
    # -------------------------------------------------------------
    # Renderização da Coluna Esquerda: Curvas de Rho
    # -------------------------------------------------------------
    y_max_global = 7.2
    
    for i, d in enumerate(desfechos):
        ax = axes_dens[i]
        st = stats[d]
        cor = CORES[d]
        
        ax.grid(True, linestyle=":", alpha=0.5, color=CORES["grade"], zorder=1)
        ax.set_axisbelow(True)
        
        # Gradiente suave translúcido
        criar_gradiente_vertical(ax, st["rho"], st["densidade"], cor, y_max_ref=y_max_global,
                                 alpha_min=0.04, alpha_max=0.50)
        
        # Curva de densidade
        ax.plot(st["rho"], st["densidade"], color=cor, linewidth=2.2, zorder=4)
        
        # Intervalo HPD de 95% sombreado
        mask_hpd = (st["rho"] >= st["hpd_inf"]) & (st["rho"] <= st["hpd_sup"])
        ax.fill_between(
            st["rho"][mask_hpd],
            0,
            st["densidade"][mask_hpd],
            color=cor,
            alpha=0.20,
            zorder=3
        )
        ax.vlines([st["hpd_inf"], st["hpd_sup"]], ymin=0,
                  ymax=[st["cutoff"], st["cutoff"]],
                  color=cor, linestyle="--", linewidth=1.0, alpha=0.7, zorder=5)
        
        # Média Posterior de rho
        pico_y = st["densidade"].max()
        ax.vlines(st["media"], ymin=0, ymax=pico_y + 0.35, color=cor,
                  linestyle="-", linewidth=2.0, zorder=6)
        ax.plot(st["media"], pico_y + 0.35, marker="v", markersize=6, color=cor, zorder=6)
        
        # Moran I de y
        ax.vlines(st["moran"], ymin=0, ymax=y_max_global * 0.82, color=CORES["moran"],
                  linestyle=(0, (3, 2)), linewidth=1.8, zorder=5)
        ax.plot(st["moran"], y_max_global * 0.82, marker="o", markersize=4.5,
                color=CORES["moran"], zorder=5)
        
        # Badge do desfecho
        badge_box = dict(boxstyle="square,pad=0.35", facecolor=CORES["fundo"],
                         edgecolor=cor, linewidth=1.4)
        ax.text(0.02, 0.84, ROTULOS[d], transform=ax.transAxes,
                fontsize=11.0, fontweight="bold", color=cor, bbox=badge_box, va="top", zorder=7)
        
        # Badges com valores
        ax.text(st["media"] + 0.018, pico_y * 0.90,
                f"$\\bar{{\\rho}}_{{\\mathrm{{post}}}} = {st['media']:.2f}$\n(HPD: [{st['hpd_inf']:.2f}; {st['hpd_sup']:.2f}])",
                fontsize=9.2, fontweight="bold", color=cor, va="center",
                bbox=dict(boxstyle="round,pad=0.25", facecolor="#FFFFFF", edgecolor=cor, alpha=0.92, linewidth=0.8),
                zorder=7)
        
        pos_y_moran = pico_y * 0.45 if d != "letalidade" else pico_y * 0.75
        align_moran = "right" if st["moran"] < st["media"] else "left"
        offset_moran = -0.018 if align_moran == "right" else 0.018
        ax.text(st["moran"] + offset_moran, pos_y_moran,
                f"Moran $I(y) = {st['moran']:.2f}$",
                fontsize=8.8, fontweight="bold", color=CORES["moran"], va="center", ha=align_moran,
                bbox=dict(boxstyle="round,pad=0.2", facecolor="#FFFFFF", edgecolor=CORES["moran"], alpha=0.88, linewidth=0.7),
                zorder=7)
        
        # Seta Moran -> Rho
        delta = st["media"] - st["moran"]
        y_arrow = pico_y * 0.22
        ax.annotate(
            "",
            xy=(st["media"], y_arrow),
            xytext=(st["moran"], y_arrow),
            arrowprops=dict(arrowstyle="->", color=cor, lw=1.3, shrinkA=3, shrinkB=3),
            zorder=6
        )
        ax.text((st["moran"] + st["media"]) / 2, y_arrow + 0.30,
                f"$\\Delta = {delta:+.2f}$",
                ha="center", va="bottom", fontsize=8.0, fontweight="bold", color=cor,
                bbox=dict(boxstyle="round,pad=0.15", facecolor="#FFFFFF", edgecolor="none", alpha=0.85),
                zorder=7)
        
        # Ajustes de Eixos
        ax.set_xlim(-0.05, 1.0)
        ax.set_ylim(0, y_max_global)
        ax.xaxis.set_major_locator(ticker.MultipleLocator(0.2))
        ax.xaxis.set_minor_locator(ticker.MultipleLocator(0.1))
        ax.yaxis.set_major_locator(ticker.MultipleLocator(2.0))
        
        if i < 2:
            ax.set_xticklabels([])
        else:
            ax.set_xlabel("Coeficiente de Autocorrelação Espacial ($\\rho$ no SAR e $I$ de Moran)",
                          fontsize=10.0, fontweight="bold", labelpad=7)
            
        ax.set_ylabel("Densidade", fontsize=9.2, labelpad=5)
        ax.axhline(0, color="#888888", linewidth=0.8, zorder=3)
    
    axes_dens[0].set_title(
        "(A) Curvas Posteriores $p(\\rho \\mid y)$ sob BMA vs. Autocorrelação Bruta $I(y)$",
        fontsize=11.5, fontweight="bold", pad=10, loc="left", color=CORES["tinta"]
    )
    
    # -------------------------------------------------------------
    # Coluna Direita - Painel Superior: Moran vs Rho
    # -------------------------------------------------------------
    ax_comp.set_title(
        "(B) I de Moran e média posterior de rho",
        fontsize=11.2, fontweight="bold", pad=10, loc="left", color=CORES["tinta"]
    )
    ax_comp.grid(True, linestyle=":", alpha=0.5, color=CORES["grade"], zorder=1)
    ax_comp.set_axisbelow(True)
    
    y_pos = [2, 1, 0]  # Casos, Óbitos, Letalidade
    
    for idx, d in enumerate(desfechos):
        y_val = y_pos[idx]
        cor = CORES[d]
        m_val = stats[d]["moran"]
        r_val = stats[d]["media"]
        
        ax_comp.plot([m_val, r_val], [y_val, y_val], color=cor, linewidth=3.5,
                     solid_capstyle="round", alpha=0.85, zorder=3)
        
        ax_comp.scatter(m_val, y_val, color="#FFFFFF", edgecolors=CORES["moran"],
                        s=110, linewidth=2.2, zorder=5)
        ax_comp.scatter(m_val, y_val, color=CORES["moran"], s=28, zorder=6)
        
        ax_comp.scatter(r_val, y_val, color=cor, edgecolors="#FFFFFF",
                        s=140, linewidth=1.5, zorder=5)
        
        ax_comp.text(m_val - 0.025, y_val, f"{m_val:.2f}", ha="right", va="center",
                     fontsize=9.2, fontweight="bold", color=CORES["moran"])
        ax_comp.text(r_val + 0.025, y_val, f"{r_val:.2f}", ha="left", va="center",
                     fontsize=9.2, fontweight="bold", color=cor)
        
        ax_comp.text(-0.02, y_val + 0.28, ROTULOS[d], fontsize=9.5, fontweight="bold",
                     color=cor, ha="left")
    
    # Destacar a comparacao entre casos e obitos
    rect_contraste = FancyBboxPatch(
        (0.38, 0.65), 0.39, 1.68,
        boxstyle="round,pad=0.03",
        linewidth=1.2, edgecolor="#C0392B", facecolor="#FDEDEC",
        alpha=0.45, zorder=2, linestyle="--"
    )
    ax_comp.add_patch(rect_contraste)
    
    ax_comp.text(
        0.83, 1.5,
        "Diferença no Moran:\n$|0{,}48 - 0{,}44| = 0{,}04$\n\nDiferença em $\\rho$:\n$|0{,}71 - 0{,}55| = 0{,}16$",
        fontsize=8.3, fontweight="bold", color="#78281F", va="center",
        bbox=dict(boxstyle="round,pad=0.3", facecolor="#FFFFFF", edgecolor="#C0392B", alpha=0.92, lw=0.9),
        zorder=6
    )
    
    ax_comp.set_xlim(-0.05, 1.05)
    ax_comp.set_ylim(-0.6, 2.7)
    ax_comp.set_yticks([])
    ax_comp.set_xlabel("Magnitude da Autocorrelação", fontsize=9.2, fontweight="bold", labelpad=5)
    
    legend_elements = [
        Line2D([0], [0], marker="o", color="w", markerfacecolor=CORES["moran"],
               markeredgecolor="#FFFFFF", markersize=8, label="Moran $I(y)$ bruto"),
        Line2D([0], [0], marker="o", color="w", markerfacecolor="#444444",
               markeredgecolor="#FFFFFF", markersize=9, label="Média $\\rho$ posterior (SAR)")
    ]
    ax_comp.legend(handles=legend_elements, loc="lower right", fontsize=8.2,
                   frameon=True, facecolor="#FFFFFF", framealpha=0.92)
    
    # -------------------------------------------------------------
    # Lower-right panel: D ratios against posterior rho
    # -------------------------------------------------------------
    ax_sev.set_title(
        "(C) Razão $D_{\\mathrm{SAR}} / D_{\\mathrm{MQO}}$",
        fontsize=11.2, fontweight="bold", pad=10, loc="left", color=CORES["tinta"]
    )
    ax_sev.grid(True, linestyle=":", alpha=0.5, color=CORES["grade"], zorder=1)
    ax_sev.set_axisbelow(True)
    
    ax_sev.axhspan(0.0, 0.50, color="#FADBD8", alpha=0.55, zorder=2)
    ax_sev.axhspan(0.50, 1.0, color="#D5F5E3", alpha=0.35, zorder=2)
    
    ax_sev.axhline(0.50, color="#C0392B", linestyle="--", linewidth=1.4, zorder=3)
    ax_sev.text(0.08, 0.465, "Referência visual: $D_{\\mathrm{SAR}}/D_{\\mathrm{MQO}} = 0{,}50$",
                color="#922B21", fontsize=8.2, fontweight="bold", zorder=6)
    
    ax_sev.text(0.08, 0.08, "$D_{\\mathrm{SAR}}/D_{\\mathrm{MQO}} < 0{,}50$",
                color="#78281F", fontsize=8.0, fontweight="bold", zorder=6)
    ax_sev.text(0.08, 0.91, "$D_{\\mathrm{SAR}}/D_{\\mathrm{MQO}} > 0{,}50$",
                color="#1E8449", fontsize=8.0, fontweight="bold", zorder=6)
    
    # Plotar os pontos com posições de anotações desempacotadas
    for d in desfechos:
        cor = CORES[d]
        rho_val = float(df_sev.loc[d, "rho"])
        razao_val = float(df_sev.loc[d, "razao"])
        
        ax_sev.scatter(rho_val, razao_val, color=cor, edgecolors="#FFFFFF",
                       s=160, linewidth=2.0, zorder=7)
        
        if d == "casos":
            # Ponto em (0.71, 0.30): anotação à esquerda no quadrante vermelho
            ax_sev.annotate(
                f"{ROTULOS[d]}\n$\\bar{{\\rho}} = {rho_val:.2f}$  |  Razão $D = {razao_val:.2f}$\nD no SAR corresponde a cerca de 30% do D no MQO",
                xy=(rho_val, razao_val),
                xytext=(rho_val - 0.32, razao_val - 0.14),
                bbox=dict(boxstyle="round,pad=0.35", facecolor="#FFFFFF", edgecolor=cor, lw=1.2),
                arrowprops=dict(arrowstyle="->", color=cor, lw=1.4, shrinkB=6),
                fontsize=8.3, fontweight="bold", color=cor, zorder=8
            )
        elif d == "obitos":
            # Ponto em (0.55, 0.71): anotação acima à direita na zona verde
            ax_sev.annotate(
                f"{ROTULOS[d]}\n$\\bar{{\\rho}} = {rho_val:.2f}$  |  Razão $D = {razao_val:.2f}$\nD no SAR corresponde a cerca de 71% do D no MQO",
                xy=(rho_val, razao_val),
                xytext=(rho_val + 0.04, razao_val + 0.12),
                bbox=dict(boxstyle="round,pad=0.35", facecolor="#FFFFFF", edgecolor=cor, lw=1.2),
                arrowprops=dict(arrowstyle="->", color=cor, lw=1.4, shrinkB=6),
                fontsize=8.3, fontweight="bold", color=cor, zorder=8
            )
        else:  # letalidade
            # Ponto em (0.25, 0.71): anotação logo abaixo à esquerda na faixa neutra (y=0.53)
            ax_sev.annotate(
                f"{ROTULOS[d]}\n$\\bar{{\\rho}} = {rho_val:.2f}$  |  Razão $D = {razao_val:.2f}$\nComparação descritiva",
                xy=(rho_val, razao_val),
                xytext=(rho_val - 0.17, razao_val - 0.18),
                bbox=dict(boxstyle="round,pad=0.35", facecolor="#FFFFFF", edgecolor=cor, lw=1.2),
                arrowprops=dict(arrowstyle="->", color=cor, lw=1.4, shrinkB=6),
                fontsize=8.3, fontweight="bold", color=cor, zorder=8
            )
            
    ax_sev.set_xlim(0.05, 0.85)
    ax_sev.set_ylim(0.0, 1.0)
    ax_sev.set_xlabel("$\\bar{\\rho}$ Posterior (SAR)", fontsize=9.2, fontweight="bold", labelpad=5)
    ax_sev.set_ylabel("Razão de D ($D_{\\mathrm{SAR}} / D_{\\mathrm{MQO}}$)",
                      fontsize=9.2, fontweight="bold", labelpad=5)
    
    # -------------------------------------------------------------
    # Titulo e nota explicativa
    # -------------------------------------------------------------
    fig.suptitle(
        "Densidades posteriores de $\\rho$, Moran e dispersão das PIPs",
        fontsize=15.0, fontweight="heavy", y=0.965, color=CORES["tinta"]
    )
    fig.text(
        0.5, 0.905,
        "Casos ($I=0{,}48$) e óbitos ($I=0{,}44$) apresentam valores próximos de Moran,\n"
        "mas diferem na média posterior de $\\rho$ e na razão de D entre SAR e MQO. A comparação é descritiva.",
        ha="center", fontsize=9.6, color="#353535", linespacing=1.35
    )
    
    png_path = os.path.join(SAIDA_FIGURAS, "fig_densidades_rho_moran.png")
    pdf_path = os.path.join(SAIDA_FIGURAS, "fig_densidades_rho_moran.pdf")
    
    plt.savefig(png_path, dpi=300, bbox_inches="tight", facecolor=CORES["fundo"])
    plt.savefig(pdf_path, dpi=300, bbox_inches="tight", facecolor=CORES["fundo"])
    plt.close()
    
    print(f"  -> Salvo: {png_path}")
    print(f"  -> Salvo: {pdf_path}")


# ---------------------------------------------------------------------------
# FIGURA 2: SENSIBILIDADE DO FATOR DE BAYES E PROBABILIDADE POSTERIOR A a0
# ---------------------------------------------------------------------------
def gerar_figura_2():
    print("Gerando Figura 2: Sensibilidade do Fator de Bayes e P(SAR | y)...")
    
    # 1. Carregar os CSVs de sensibilidade
    df_casos = pd.read_csv(os.path.join(DADOS_DIR, "sensibilidade_prior_rho_casos.csv"))
    df_obitos = pd.read_csv(os.path.join(DADOS_DIR, "sensibilidade_prior_rho_obitos.csv"))
    df_let = pd.read_csv(os.path.join(DADOS_DIR, "sensibilidade_prior_rho_letalidade.csv"))
    
    dados = {
        "casos": df_casos.sort_values("a0").reset_index(drop=True),
        "obitos": df_obitos.sort_values("a0").reset_index(drop=True),
        "letalidade": df_let.sort_values("a0").reset_index(drop=True)
    }
    
    # 2. Configurar Grid de Subplots
    # Layout limpo:
    # Coluna Esquerda dividida em dois painéis empilhados:
    #   (A1) Log BF Casos & Óbitos (evidência massiva, escala 10 a 45)
    #   (A2) Log BF Letalidade (evidência tênue, escala -1 a 4 com zoom analítico)
    # Coluna Direita:
    #   (B) Probabilidade Posterior P(SAR | y) para todos os três desfechos (50% a 104%)
    fig = plt.figure(figsize=(15.2, 7.8), dpi=300)
    
    gs_master = gridspec.GridSpec(1, 2, width_ratios=[1.15, 1.0], wspace=0.28,
                                  left=0.07, right=0.96, top=0.84, bottom=0.095)
    
    gs_left = gridspec.GridSpecFromSubplotSpec(2, 1, subplot_spec=gs_master[0], height_ratios=[1.15, 1.0], hspace=0.34)
    ax_bf1 = fig.add_subplot(gs_left[0])
    ax_bf2 = fig.add_subplot(gs_left[1])
    ax_psar = fig.add_subplot(gs_master[1])
    
    # -------------------------------------------------------------
    # Panel A1: log Bayes factors for cases and deaths
    # -------------------------------------------------------------
    ax_bf1.set_title(
        "(A1) Log do Fator de Bayes: Casos & Óbitos",
        fontsize=11.2, fontweight="bold", pad=9, loc="left"
    )
    ax_bf1.grid(True, linestyle=":", alpha=0.5, color=CORES["grade"], zorder=1)
    ax_bf1.set_axisbelow(True)
    
    # Reference region: log BF > 5
    ax_bf1.axhspan(5.0, 45.0, color="#E8F8F5", alpha=0.55, zorder=2)
    ax_bf1.axhline(5.0, color="#16A085", linestyle="--", linewidth=1.0, zorder=3)
    ax_bf1.text(1.05, 7.0, "Referência: $\\log \\mathrm{BF} = 5$",
                fontsize=8.0, fontweight="bold", color="#117864", zorder=4)
    
    for d in ["casos", "obitos"]:
        df = dados[d]
        cor = CORES[d]
        
        # Envelope Min-Max
        ax_bf1.fill_between(
            df["a0"],
            df["log_bf_cond_min"],
            df["log_bf_cond_max"],
            color=cor,
            alpha=0.18,
            zorder=4
        )
        
        # BF Agregado MC
        ax_bf1.plot(
            df["a0"],
            df["log_bf_agregado_mc"],
            marker="o",
            markersize=6.5,
            linewidth=2.4,
            color=cor,
            label=f"{ROTULOS[d]} (Agregado)",
            zorder=6
        )
        
        # Mediana Condicional
        ax_bf1.plot(
            df["a0"],
            df["log_bf_cond_mediano"],
            marker="s",
            markersize=4.0,
            linestyle=":",
            linewidth=1.3,
            color=cor,
            alpha=0.8,
            zorder=5
        )
    
    ax_bf1.text(2.2, 38.0, "Casos: $\\log \\mathrm{BF} \\approx 32 \\rightarrow 28$\n(Favorável ao SAR)",
                fontsize=8.2, fontweight="bold", color=CORES["casos"], ha="center",
                bbox=dict(boxstyle="round,pad=0.25", facecolor="#FFFFFF", edgecolor=CORES["casos"], alpha=0.92))
    
    ax_bf1.text(2.2, 21.5, "Óbitos: $\\log \\mathrm{BF} \\approx 16 \\rightarrow 14{,}5$\n(Favorável ao SAR)",
                fontsize=8.2, fontweight="bold", color=CORES["obitos"], ha="center",
                bbox=dict(boxstyle="round,pad=0.25", facecolor="#FFFFFF", edgecolor=CORES["obitos"], alpha=0.92))
    
    ax_bf1.set_xscale("log")
    ax_bf1.set_xlim(0.95, 11.5)
    ax_bf1.set_ylim(4.0, 44.0)
    ax_bf1.set_xticks([1.01, 1.1, 1.5, 2.0, 5.0, 10.0])
    ax_bf1.get_xaxis().set_major_formatter(ticker.ScalarFormatter())
    ax_bf1.set_xticklabels([])
    ax_bf1.set_ylabel("$\\log \\mathrm{BF}(\\mathrm{SAR} : \\mathrm{MQO})$",
                      fontsize=9.2, fontweight="bold", labelpad=6)
    
    # -------------------------------------------------------------
    # Painel A2: Log BF Letalidade (Zoom Analítico e Sensibilidade)
    # -------------------------------------------------------------
    ax_bf2.set_title(
        "(A2) Log do Fator de Bayes: Letalidade (Zoom Detalhado na Transição)",
        fontsize=11.2, fontweight="bold", pad=9, loc="left"
    )
    ax_bf2.grid(True, linestyle=":", alpha=0.5, color=CORES["grade"], zorder=1)
    ax_bf2.set_axisbelow(True)
    
    # Faixas de interpretação
    ax_bf2.axhspan(1.0, 4.0, color="#FEF9E7", alpha=0.50, zorder=2)
    ax_bf2.axhspan(0.0, 1.0, color="#FDFEFE", alpha=0.30, zorder=2)
    ax_bf2.axhspan(-1.5, 0.0, color="#FDEDEC", alpha=0.45, zorder=2)
    
    ax_bf2.axhline(1.0, color="#B7950B", linestyle="--", linewidth=1.0, zorder=3)
    ax_bf2.axhline(0.0, color="#7F8C8D", linestyle="-", linewidth=1.2, zorder=3)
    
    ax_bf2.text(1.04, 3.2, "Evidência Positiva / Moderada ($1 < \\log \\mathrm{BF} \\leq 3$)",
                fontsize=8.0, fontweight="bold", color="#9A7D0A", zorder=4)
    ax_bf2.text(1.04, 0.35, "Evidência Fraca / Anedótica ($0 < \\log \\mathrm{BF} \\leq 1$)",
                fontsize=8.0, fontweight="bold", color="#7D6608", zorder=4)
    ax_bf2.text(1.04, -1.1, "Evidência Favorável ao MQO ($\\log \\mathrm{BF} < 0$)",
                fontsize=8.0, fontweight="bold", color="#922B21", zorder=4)
    
    df_l = dados["letalidade"]
    cor_l = CORES["letalidade"]
    
    # Envelope Min-Max de Letalidade
    ax_bf2.fill_between(
        df_l["a0"],
        df_l["log_bf_cond_min"],
        df_l["log_bf_cond_max"],
        color=cor_l,
        alpha=0.20,
        zorder=4
    )
    
    ax_bf2.plot(
        df_l["a0"],
        df_l["log_bf_agregado_mc"],
        marker="o",
        markersize=6.5,
        linewidth=2.4,
        color=cor_l,
        label=f"{ROTULOS['letalidade']} (Agregado)",
        zorder=6
    )
    
    ax_bf2.plot(
        df_l["a0"],
        df_l["log_bf_cond_mediano"],
        marker="s",
        markersize=4.0,
        linestyle=":",
        linewidth=1.3,
        color=cor_l,
        alpha=0.8,
        label="Mediana Condicional",
        zorder=5
    )
    
    # Nota explicativa sobre a travessia do zero
    ax_bf2.text(
        4.0, -0.65,
        "Min. condicional $< 0$ para $a_0 \\leq 2$:\nEm submodelos com ruído,\no MQO supera o SAR!",
        fontsize=8.0, fontweight="bold", color="#78281F",
        bbox=dict(boxstyle="round,pad=0.2", facecolor="#FFFFFF", edgecolor="#C0392B", lw=0.9),
        zorder=7
    )
    
    ax_bf2.set_xscale("log")
    ax_bf2.set_xlim(0.95, 11.5)
    ax_bf2.set_ylim(-1.4, 3.8)
    ax_bf2.set_xticks([1.01, 1.1, 1.5, 2.0, 5.0, 10.0])
    ax_bf2.get_xaxis().set_major_formatter(ticker.ScalarFormatter())
    ax_bf2.set_xlabel("Hiperparâmetro $a_0$ da Priori $\\mathrm{Beta}(a_0, a_0)$ sobre $\\rho$ (Escala Log)",
                      fontsize=9.2, fontweight="bold", labelpad=6)
    ax_bf2.set_ylabel("$\\log \\mathrm{BF}(\\mathrm{SAR} : \\mathrm{MQO})$",
                      fontsize=9.2, fontweight="bold", labelpad=6)
    
    # -------------------------------------------------------------
    # Painel B: Probabilidade Posterior P(SAR | y)
    # -------------------------------------------------------------
    ax_psar.set_title(
        "(B) Probabilidade Posterior Agregada $P(\\mathrm{SAR} \\mid y)$ vs. $a_0$",
        fontsize=11.8, fontweight="bold", pad=10, loc="left"
    )
    ax_psar.grid(True, linestyle=":", alpha=0.5, color=CORES["grade"], zorder=1)
    ax_psar.set_axisbelow(True)
    
    for d in ["casos", "obitos", "letalidade"]:
        df = dados[d]
        cor = CORES[d]
        
        ax_psar.plot(
            df["a0"],
            df["p_sar_agregada_mc"],
            marker="o",
            markersize=7.5,
            linewidth=2.6,
            color=cor,
            label=ROTULOS[d],
            zorder=5
        )
    
    # Probabilidades estimadas nos valores de a0 avaliados
    ax_psar.annotate(
        "Casos e óbitos\n$P(\\mathrm{SAR} \\mid y) \\approx 1$\n"
        "Estimativas próximas de 1 nos valores\nde $a_0$ avaliados.",
        xy=(2.0, 1.0),
        xytext=(1.2, 0.89),
        bbox=dict(boxstyle="round,pad=0.4", facecolor="#FFFFFF", edgecolor="#333333", lw=1.2),
        arrowprops=dict(arrowstyle="->", color="#333333", lw=1.4, shrinkB=6),
        fontsize=8.8, fontweight="bold", color="#1C1C1C", zorder=7
    )
    
    # Destacar a sensibilidade de Letalidade
    p_min_let = df_let["p_sar_agregada_mc"].iloc[0] * 100
    p_max_let = df_let["p_sar_agregada_mc"].iloc[-1] * 100
    ax_psar.annotate(
        f"Letalidade: Sensibilidade Moderada\n$P(\\mathrm{{SAR}} \\mid y)$ varia de {p_min_let:.1f}% ($a_0=1{{,}}01$) a {p_max_let:.1f}% ($a_0=10$)\n"
        "Com sinal espacial fraco, a priori em $\\rho$\nmodula a partição posterior entre SAR e MQO.",
        xy=(5.0, df_let.loc[df_let["a0"] == 5.0, "p_sar_agregada_mc"].values[0]),
        xytext=(1.35, 0.58),
        bbox=dict(boxstyle="round,pad=0.4", facecolor="#FFFFFF", edgecolor=CORES["letalidade"], lw=1.2),
        arrowprops=dict(arrowstyle="->", color=CORES["letalidade"], lw=1.4, shrinkB=6),
        fontsize=8.8, fontweight="bold", color=CORES["letalidade"], zorder=7
    )
    
    ax_psar.set_xscale("log")
    ax_psar.set_xlim(0.95, 11.5)
    ax_psar.set_ylim(0.50, 1.04)
    ax_psar.set_xticks([1.01, 1.1, 1.5, 2.0, 5.0, 10.0])
    ax_psar.get_xaxis().set_major_formatter(ticker.ScalarFormatter())
    ax_psar.yaxis.set_major_formatter(ticker.PercentFormatter(xmax=1.0, decimals=0))
    ax_psar.set_xlabel("Hiperparâmetro $a_0$ da Priori $\\mathrm{Beta}(a_0, a_0)$ sobre $\\rho$ (Escala Log)",
                       fontsize=9.5, fontweight="bold", labelpad=7)
    ax_psar.set_ylabel("Probabilidade Posterior $P(\\mathrm{SAR} \\mid y)$",
                       fontsize=9.5, fontweight="bold", labelpad=6)
    
    ax_psar.axhline(1.0, color="#27AE60", linestyle="--", linewidth=1.0, alpha=0.7)
    ax_psar.axhline(0.5, color="#E74C3C", linestyle=":", linewidth=1.0, alpha=0.7)
    
    ax_psar.legend(loc="lower right", fontsize=9.2, frameon=True, facecolor="#FFFFFF", framealpha=0.92)
    
    # -------------------------------------------------------------
    # Título Geral e Rodapé Informativo
    # -------------------------------------------------------------
    fig.suptitle(
        "Sensibilidade dos fatores de Bayes e das probabilidades à priori de $\\rho$",
        fontsize=15.0, fontweight="heavy", y=0.965, color=CORES["tinta"]
    )
    fig.text(
        0.5, 0.905,
        "Avaliação do Fator de Bayes e de $P(\\mathrm{SAR} \\mid y)$ sob variação da concentração $a_0 \\in [1{,}01; 10{,}0]$ na priori simétrica de $\\rho$.\n"
        "$a_0 = 1{,}01$ reflete indiferença quase uniforme; $a_0 = 10{,}0$ concentra densidade em $\\rho=0$, favorecendo fortemente o modelo clássico (MQO).",
        ha="center", fontsize=9.4, color="#353535", linespacing=1.35
    )
    
    png_path = os.path.join(SAIDA_FIGURAS, "fig_bayes_factor_sensibilidade.png")
    pdf_path = os.path.join(SAIDA_FIGURAS, "fig_bayes_factor_sensibilidade.pdf")
    
    plt.savefig(png_path, dpi=300, bbox_inches="tight", facecolor=CORES["fundo"])
    plt.savefig(pdf_path, dpi=300, bbox_inches="tight", facecolor=CORES["fundo"])
    plt.close()
    
    print(f"  -> Salvo: {png_path}")
    print(f"  -> Salvo: {pdf_path}")


# ---------------------------------------------------------------------------
# Execução Principal
# ---------------------------------------------------------------------------
if __name__ == "__main__":
    print("=" * 70)
    print("Executando geração de figuras de severidade e sensibilidade de rho...")
    print("=" * 70)
    
    gerar_figura_1()
    gerar_figura_2()
    
    print("=" * 70)
    print("Processo concluído com sucesso!")
    print(f"Arquivos salvos no diretório: {SAIDA_FIGURAS}")
    print("=" * 70)
