"""
Mapas municipais e diagramas de Moran para a análise de COVID-19 no Ceará.
Usa a malha e as tabelas de resultados incluídas no projeto.
"""

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from caminhos import ESTIMACOES

import os
import numpy as np
import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.patches import Patch
import libpysal as lp
from esda.moran import Moran
import seaborn as sns

# Caminhos
AQUI = str(Path(__file__).resolve().parent)
RAIZ = str(ESTIMACOES)
DIR_SAIDA = os.path.join(AQUI, "saida", "figuras_novas")
os.makedirs(DIR_SAIDA, exist_ok=True)

# Estilo dos graficos
plt.rcParams['font.sans-serif'] = ['Segoe UI', 'DejaVu Sans', 'Helvetica', 'Arial']
plt.rcParams['font.family'] = 'sans-serif'

# Carregar malha
print("Carregando malha e bases analíticas...")
malha_path = os.path.join(AQUI, "saida", "export_python", "malha_ce.geojson")
gdf = gpd.read_file(malha_path)
gdf['code_muni'] = gdf['code_muni'].astype(int)

# Dados de casos, óbitos e letalidade
df_casos = pd.read_csv(os.path.join(RAIZ, "resultados", "estimacoes", "matrizes_estimacao", "principal", "casos.csv"))
df_obitos = pd.read_csv(os.path.join(RAIZ, "resultados", "estimacoes", "matrizes_estimacao", "principal", "obitos.csv"))
df_letal = pd.read_csv(os.path.join(RAIZ, "resultados", "estimacoes", "matrizes_estimacao", "principal", "letalidade.csv"))

df_casos['codigo_ibge'] = df_casos['codigo_ibge'].astype(int)
df_obitos['codigo_ibge'] = df_obitos['codigo_ibge'].astype(int)
df_letal['codigo_ibge'] = df_letal['codigo_ibge'].astype(int)

gdf = gdf.merge(df_casos, left_on='code_muni', right_on='codigo_ibge')
gdf['z_log_obitos100k'] = df_obitos['z_log_obitos100k'].values
gdf['z_logit_letal'] = df_letal['z_logit_letal'].values

# Matriz Queen
w_queen = lp.weights.Queen.from_dataframe(gdf, use_index=False)
w_queen.transform = 'R'

polos = {
    'Fortaleza': (-38.5267, -3.7319),
    'Sobral': (-40.3497, -3.6894),
    'Juazeiro do Norte': (-39.3150, -7.2131)
}

# ==============================================================================
# FIGURA 1: MAPAS COMPARATIVOS
# ==============================================================================
print("Gerando Figura 1: mapas comparativos...")

fig = plt.figure(figsize=(17, 6.5), facecolor='#FAFAFA')

# Grid com 4 colunas para os 4 mapas
grid = plt.GridSpec(1, 4, figure=fig, left=0.03, right=0.97, top=0.76, bottom=0.20, wspace=0.10)

configs = [
    {
        'col': 'z_log_casos100k',
        'letra': 'A',
        'titulo': '(A) Variável Dependente',
        'subtitulo': 'Casos por 100 mil hab. (z-score)',
        'moran': 0.483,
        'tag': 'I de Moran = 0,483',
        'tag_cor': '#1E40AF',
        'is_dummy': False
    },
    {
        'col': 'z_idosos',
        'letra': 'B',
        'titulo': '(B) Covariável Suave',
        'subtitulo': 'Proporção de Idosos (z-score)',
        'moran': 0.662,
        'tag': 'PIP: MQO 0,95 → SAR 0,30',
        'tag_cor': '#B91C1C',
        'is_dummy': False
    },
    {
        'col': 'metrop',
        'letra': 'C',
        'titulo': '(C) Estrutura Polarizada',
        'subtitulo': 'Região Metropolitana de Fortaleza',
        'moran': 0.789,
        'tag': 'I de Moran = 0,789',
        'tag_cor': '#6D28D9',
        'is_dummy': True
    },
    {
        'col': 'z_tax.agua',
        'letra': 'D',
        'titulo': '(D) Covariável Rugosa',
        'subtitulo': 'Cobertura de Água (z-score)',
        'moran': 0.000,
        'tag': 'I de Moran ≈ 0,000',
        'tag_cor': '#047857',
        'is_dummy': False
    }
]

cmap_div = sns.color_palette("vlag", as_cmap=True)
norm_div = mcolors.Normalize(vmin=-2.6, vmax=2.6)

for idx, cfg in enumerate(configs):
    ax = fig.add_subplot(grid[0, idx])
    ax.set_facecolor('#FAFAFA')
    ax.axis('off')
    
    if not cfg['is_dummy']:
        gdf.plot(
            column=cfg['col'],
            cmap=cmap_div,
            norm=norm_div,
            linewidth=0.3,
            edgecolor='#FFFFFF',
            ax=ax
        )
    else:
        # Plotagem distinta para a dummy binária da RMF
        gdf[gdf['metrop'] == 0].plot(
            color='#E2E8F0',
            linewidth=0.3,
            edgecolor='#CBD5E1',
            ax=ax
        )
        gdf[gdf['metrop'] == 1].plot(
            color='#6D28D9',
            linewidth=0.4,
            edgecolor='#4C1D95',
            ax=ax
        )
        leg_handles = [
            Patch(facecolor='#6D28D9', edgecolor='#4C1D95', label='RMF (19 mun.)'),
            Patch(facecolor='#E2E8F0', edgecolor='#CBD5E1', label='Interior (165 mun.)')
        ]
        ax.legend(handles=leg_handles, loc='lower right', fontsize=7.5, 
                  frameon=True, facecolor='#FFFFFF', edgecolor='#E2E8F0', framealpha=0.95)

    # Anotação de polos no Painel A
    if idx == 0:
        for cidade, (px, py) in polos.items():
            ax.scatter(px, py, color='#0F172A', s=14, zorder=5, edgecolor='#FFFFFF', linewidth=0.5)
            offset = (6, 3) if cidade != 'Juazeiro do Norte' else (6, -6)
            ax.annotate(cidade, xy=(px, py), xytext=offset, textcoords='offset points',
                        fontsize=7.0, fontweight='bold', color='#0F172A',
                        bbox=dict(boxstyle='round,pad=0.15', facecolor='#FFFFFF', alpha=0.85, edgecolor='none'))

    # Títulos e Subtítulos estruturados NO TOPO DO AXIS
    ax.set_title(f"{cfg['titulo']}\n{cfg['subtitulo']}", 
                 fontsize=10.0, fontweight='bold', color='#1E293B', pad=10, linespacing=1.2)
    
    # Badge de status no rodapé do painel
    badge_box = dict(boxstyle='round,pad=0.35', facecolor='#FFFFFF', edgecolor=cfg['tag_cor'], alpha=0.95, linewidth=1.0)
    moran_fmt = f"{cfg['moran']:.3f}".replace('.', ',')
    ax.text(0.5, -0.06, f"I de Moran = {moran_fmt} | {cfg['tag']}", 
            transform=ax.transAxes, fontsize=7.6, fontweight='bold', 
            color=cfg['tag_cor'], ha='center', va='top', bbox=badge_box)

# Barra de cores horizontal para as variáveis contínuas
cbar_ax = fig.add_axes([0.25, 0.05, 0.50, 0.022])
sm = plt.cm.ScalarMappable(cmap=cmap_div, norm=norm_div)
sm.set_array([])
cbar = fig.colorbar(sm, cax=cbar_ax, orientation='horizontal')
cbar.set_label('Escala Padronizada (z-score com média 0 e desvio-padrão 1)', fontsize=9.0, color='#334155', labelpad=4)
cbar.set_ticks([-2.5, -1.5, 0, 1.5, 2.5])
cbar.set_ticklabels(['-2,5 (Muito Baixo)', '-1,5', '0 (Média)', '+1,5', '+2,5 (Muito Alto)'])
cbar.ax.tick_params(labelsize=8.0, color='#64748B')
cbar.outline.set_edgecolor('#CBD5E1')
cbar.outline.set_linewidth(0.8)

# Cabeçalho da figura bem espaçado
fig.suptitle('A Geografia da Dependência Espacial: Variável Dependente vs. Estrutura das Covariáveis',
             fontsize=14.0, fontweight='bold', color='#0F172A', y=0.97)
fig.text(0.5, 0.92,
         '184 municípios do Ceará: casos, proporção de idosos, região metropolitana e cobertura de água.',
         ha='center', fontsize=9.5, color='#475569')

out_png1 = os.path.join(DIR_SAIDA, "fig1_mapas_python.png")
out_pdf1 = os.path.join(DIR_SAIDA, "fig1_mapas_python.pdf")
plt.savefig(out_png1, dpi=300, facecolor='#FAFAFA', bbox_inches='tight')
plt.savefig(out_pdf1, facecolor='#FAFAFA', bbox_inches='tight')
plt.close()
print(f"Salvo com sucesso: {out_png1}")


# ==============================================================================
# FIGURA 2: OS TRÊS DESFECHOS E DIAGRAMAS DE MORAN (LAYOUT 2 LINHAS x 3 COLUNAS)
# ==============================================================================
print("Gerando Figura 2: Tríade dos Desfechos e Moran em 2 Linhas...")

fig = plt.figure(figsize=(16, 11.0), facecolor='#FAFAFA')

# Grid 2 linhas x 3 colunas: Linha 0 = Mapas, Linha 1 = Scatter de Moran
grid2 = plt.GridSpec(2, 3, figure=fig, height_ratios=[1.2, 0.9], 
                     left=0.06, right=0.94, top=0.89, bottom=0.08, hspace=0.32, wspace=0.15)

desfechos_cfg = [
    {
        'col': 'z_log_casos100k',
        'letra': 'A',
        'nome': 'Casos de COVID-19',
        'sub': 'Log de casos por 100 mil hab. (z-score)',
        'moran': 0.483,
        'rho_sar': 0.708,
        'razao_d': 0.301,
        'cor': '#1D4ED8',
        'status_curto': 'Razão de D = 0,30 | ρ = 0,71'
    },
    {
        'col': 'z_log_obitos100k',
        'letra': 'B',
        'nome': 'Óbitos de COVID-19',
        'sub': 'Log de óbitos por 100 mil hab. (z-score)',
        'moran': 0.441,
        'rho_sar': 0.555,
        'razao_d': 0.709,
        'cor': '#D97706',
        'status_curto': 'Inversão Seletiva (Razão = 0,71 | ρ = 0,55)'
    },
    {
        'col': 'z_logit_letal',
        'letra': 'C',
        'nome': 'Taxa de Letalidade',
        'sub': 'Logit empírico óbitos/casos (z-score)',
        'moran': 0.138,
        'rho_sar': 0.254,
        'razao_d': 0.706,
        'cor': '#059669',
        'status_curto': 'Fragilidade de Limiar (Razão = 0,71 | ρ = 0,25)'
    }
]

cmap_seq = sns.color_palette("mako_r", as_cmap=True)
norm_seq = mcolors.Normalize(vmin=-2.5, vmax=2.5)

for idx, dcfg in enumerate(desfechos_cfg):
    # ----------------------------------------------------
    # LINHA 0: MAPA COROPLÉTICO
    # ----------------------------------------------------
    ax_map = fig.add_subplot(grid2[0, idx])
    ax_map.set_facecolor('#FAFAFA')
    ax_map.axis('off')
    
    gdf.plot(
        column=dcfg['col'],
        cmap=cmap_seq,
        norm=norm_seq,
        linewidth=0.3,
        edgecolor='#FFFFFF',
        ax=ax_map
    )
    
    # Título do mapa acima do painel
    ax_map.set_title(f"({dcfg['letra']}) {dcfg['nome']}\n{dcfg['sub']}", 
                     fontsize=10.5, fontweight='bold', color='#0F172A', pad=10, linespacing=1.2)
    
    # Badge do diagnóstico espacial posicionado logo abaixo do mapa
    box_diag = dict(boxstyle='round,pad=0.35', facecolor='#FFFFFF', edgecolor=dcfg['cor'], alpha=0.95, linewidth=1.1)
    ax_map.text(0.5, -0.05, dcfg['status_curto'], transform=ax_map.transAxes, 
                fontsize=7.8, fontweight='bold', color=dcfg['cor'], 
                ha='center', va='top', bbox=box_diag)
    
    # ----------------------------------------------------
    # LINHA 1: DIAGRAMA DE MORAN (MORAN SCATTER PLOT)
    # ----------------------------------------------------
    ax_moran = fig.add_subplot(grid2[1, idx])
    ax_moran.set_facecolor('#FFFFFF')
    
    y_val = gdf[dcfg['col']].values
    w_y_val = lp.weights.lag_spatial(w_queen, y_val)
    moran_obj = Moran(y_val, w_queen)
    
    # Linhas de quadrante
    ax_moran.axvline(0, color='#E2E8F0', linestyle='--', linewidth=0.8, zorder=1)
    ax_moran.axhline(0, color='#E2E8F0', linestyle='--', linewidth=0.8, zorder=1)
    
    # Scatter points com cor temática
    ax_moran.scatter(y_val, w_y_val, s=18, color=dcfg['cor'], alpha=0.70, edgecolors='none', zorder=2)
    
    # Linha de regressão do Moran
    x_line = np.linspace(-3.2, 3.2, 50)
    ax_moran.plot(x_line, moran_obj.I * x_line, color='#0F172A', linewidth=1.6, zorder=3,
                  label=f"Inclinação = {moran_obj.I:.3f}".replace('.', ','))
    
    ax_moran.set_xlim(-3.2, 3.2)
    ax_moran.set_ylim(-3.2, 3.2)
    ax_moran.set_xlabel('Valor do Município (z-score)', fontsize=8.2, color='#334155')
    ax_moran.set_ylabel('Defasagem Espacial Vizinhos (Wy)', fontsize=8.2, color='#334155')
    ax_moran.tick_params(labelsize=7.5, colors='#64748B')
    
    # Quadrantes descritivos sutis
    ax_moran.text(0.95, 0.93, 'Alto-Alto', transform=ax_moran.transAxes, fontsize=7.0, color='#94A3B8', ha='right', va='top')
    ax_moran.text(0.05, 0.07, 'Baixo-Baixo', transform=ax_moran.transAxes, fontsize=7.0, color='#94A3B8', ha='left', va='bottom')
    ax_moran.text(0.05, 0.93, 'Baixo-Alto', transform=ax_moran.transAxes, fontsize=7.0, color='#94A3B8', ha='left', va='top')
    ax_moran.text(0.95, 0.07, 'Alto-Baixo', transform=ax_moran.transAxes, fontsize=7.0, color='#94A3B8', ha='right', va='bottom')
    
    moran_str = f"{moran_obj.I:.3f}".replace('.', ',')
    p_val_str = "< 0,001" if moran_obj.p_sim < 0.001 else f"= {moran_obj.p_sim:.3f}".replace('.', ',')
    ax_moran.set_title(f"Diagrama de Dispersão de Moran: I = {moran_str} (p {p_val_str})", 
                       fontsize=8.8, fontweight='bold', color='#1E293B', pad=8)
    
    for spine in ax_moran.spines.values():
        spine.set_color('#CBD5E1')
        spine.set_linewidth(0.8)

# Barra de cores horizontal na base
cbar_ax2 = fig.add_axes([0.30, 0.025, 0.40, 0.016])
sm2 = plt.cm.ScalarMappable(cmap=cmap_seq, norm=norm_seq)
sm2.set_array([])
cbar2 = fig.colorbar(sm2, cax=cbar_ax2, orientation='horizontal')
cbar2.set_label('Taxa Padronizada do Desfecho (z-score municipal)', fontsize=8.2, color='#334155', labelpad=3)
cbar2.set_ticks([-2, -1, 0, 1, 2])
cbar2.set_ticklabels(['-2,0 (Subnotificado / Baixo)', '-1,0', '0 (Média)', '+1,0', '+2,0 (Incidência Severa)'])
cbar2.ax.tick_params(labelsize=7.5, color='#64748B')
cbar2.outline.set_edgecolor('#CBD5E1')
cbar2.outline.set_linewidth(0.8)

# Cabeçalho geral da Figura 2 bem espaçado
fig.suptitle('A Geografia da Pandemia no Ceará e os Diagnósticos de Autocorrelação Global',
             fontsize=14.5, fontweight='bold', color='#0F172A', y=0.975)
fig.text(0.5, 0.945,
         'Casos, óbitos e letalidade: distribuição municipal e estatísticas de autocorrelação espacial.',
         ha='center', fontsize=9.2, color='#475569')

out_png2 = os.path.join(DIR_SAIDA, "fig_mapas_desfechos_moran.png")
out_pdf2 = os.path.join(DIR_SAIDA, "fig_mapas_desfechos_moran.pdf")
plt.savefig(out_png2, dpi=300, facecolor='#FAFAFA', bbox_inches='tight')
plt.savefig(out_pdf2, facecolor='#FAFAFA', bbox_inches='tight')
plt.close()
print(f"Salvo com sucesso: {out_png2}")

print("Mapas e diagramas de Moran gerados.")

