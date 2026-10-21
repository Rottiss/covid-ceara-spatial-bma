# -*- coding: utf-8 -*-
"""
===============================================================================
SCRIPT: gerar_figuras_dumbbell_pips.py
PROJETO: MC3 Espacial - Seleção Bayesiana de Modelos sob Dependência Espacial
DATA: 2026-09-11
OBJETIVO:
  1. fig_dumbbell_pips.png / .pdf:
     Dumbbell Plot de PIPs (MQO -> SAR -> SEM) com Agrupamento Temático e Status.
  2. fig_bump_ranking.png / .pdf:
     Bump Chart de reorganização dos rankings das covariáveis (Top-10).
===============================================================================
"""

import os
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D
import matplotlib.ticker as ticker
from matplotlib.patches import FancyBboxPatch, Rectangle

# -----------------------------------------------------------------------------
# 1. Configurações Globais e Paleta Semântica
# -----------------------------------------------------------------------------
BASE_DIR = str(Path(__file__).resolve().parent)
CSV_PATH = os.path.join(BASE_DIR, "saida", "export_python", "pips_completas.csv")
OUT_DIR = os.path.join(BASE_DIR, "saida", "figuras_novas")
os.makedirs(OUT_DIR, exist_ok=True)

# Tipografia e Estilo Nature / The Economist
plt.rcParams['font.sans-serif'] = ['Arial', 'Helvetica', 'DejaVu Sans']
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['figure.dpi'] = 300
plt.rcParams['savefig.dpi'] = 300
plt.rcParams['savefig.bbox'] = 'tight'
plt.rcParams['axes.edgecolor'] = '#D1D5DB'
plt.rcParams['axes.linewidth'] = 0.8
plt.rcParams['figure.facecolor'] = '#FFFFFF'
plt.rcParams['axes.facecolor'] = '#FFFFFF'
plt.rcParams['text.color'] = '#0F172A'
plt.rcParams['axes.labelcolor'] = '#1E293B'
plt.rcParams['xtick.color'] = '#475569'
plt.rcParams['ytick.color'] = '#1E293B'

# Colors for threshold-crossing categories
CORES_STATUS = {
    'Sai da lista (≥ 0,50 → < 0,50)': '#DC2626',       # Red
    'Entra na lista (< 0,50 → ≥ 0,50)': '#0D9488',      # Teal
    'Permanece na lista (≥ 0,50 → ≥ 0,50)': '#16A34A',  # Green
    'Abaixo do limiar (< 0,50)': '#94A3B8'             # Cinza ardósia neutro
}

# Cores dos Grupos Temáticos
CORES_GRUPO = {
    'Território e Política': '#7C3AED',   # Violeta
    'Demografia': '#2563EB',              # Azul real
    'Vulnerabilidade e Renda': '#D97706', # Âmbar sóbrio
    'Saúde e Saneamento': '#0891B2'       # Cyan
}

ORDEM_GRUPOS = [
    'Território e Política',
    'Demografia',
    'Vulnerabilidade e Renda',
    'Saúde e Saneamento'
]

CONFIG_DESFECHOS = {
    'casos': {
        'titulo': 'CASOS POR 100 MIL HAB.',
        'subtitulo': 'Nenhuma PIP ≥ 0,50\nno SAR e no SEM',
        'moran_y': 0.483,
        'diagnostico': 'Defasagem Wy absorve o sinal de todas as covariáveis; nenhuma atinge o limiar MPM (0,50).'
    },
    'obitos': {
        'titulo': 'ÓBITOS POR 100 MIL HAB.',
        'subtitulo': 'Inversão de Liderança\n(População rural e mulheres sobem)',
        'moran_y': 0.441,
        'diagnostico': 'Idosos e pobreza colapsam; população rural e mulheres eleitas assumem o topo.'
    },
    'letalidade': {
        'titulo': 'TAXA DE LETALIDADE (%)',
        'subtitulo': 'Fragilidade de Limiar\n(9 covariáveis na margem 0,40–0,48)',
        'moran_y': 0.138,
        'diagnostico': '9 das 12 covariáveis selecionadas no MQO caem para a faixa marginal sob SAR/SEM.'
    }
}

# -----------------------------------------------------------------------------
# 2. Carga e Preparação dos Dados
# -----------------------------------------------------------------------------
def carregar_dados():
    df = pd.read_csv(CSV_PATH, encoding='utf-8')
    
    # Padronização de rótulos de grupo
    df['grupo'] = df['grupo'].replace({
        'Sade e Saneamento': 'Saúde e Saneamento',
        'Territrio e Poltica': 'Território e Política'
    })
    
    # Ordenação das covariáveis dentro de cada grupo temático:
    max_pip = df.groupby('covariavel')[['pip_mqo', 'pip_sar_queen', 'pip_sem_queen']].max().max(axis=1)
    df['max_pip_global'] = df['covariavel'].map(max_pip)
    
    return df

# -----------------------------------------------------------------------------
# 3. FIGURA 1: Dumbbell Plot de PIPs (MQO -> SAR -> SEM)
# -----------------------------------------------------------------------------
def gerar_figura_dumbbell(df):
    print("Gerando Figura 1: Dumbbell Plot de PIPs...")
    
    # Estrutura de ordenação temática (de cima para baixo na leitura)
    cov_info = df[['covariavel', 'rotulo', 'grupo', 'max_pip_global']].drop_duplicates()
    cov_info['grupo_order'] = cov_info['grupo'].map(lambda g: ORDEM_GRUPOS.index(g) if g in ORDEM_GRUPOS else 99)
    # Inverter para que o primeiro grupo apareça no topo
    cov_info = cov_info.sort_values(by=['grupo_order', 'max_pip_global'], ascending=[False, True]).reset_index(drop=True)
    
    # Inserir espaçamento visual entre grupos
    y_positions = {}
    y_current = 0.0
    last_group = None
    group_boundaries = {}
    
    for idx, row in cov_info.iterrows():
        g = row['grupo']
        if last_group is not None and g != last_group:
            y_current += 1.4  # Espaço extra entre grupos temáticos
        if g not in group_boundaries:
            group_boundaries[g] = {'start': y_current}
        group_boundaries[g]['end'] = y_current
        y_positions[row['covariavel']] = y_current
        y_current += 1.0
        last_group = g
        
    cov_info['y_pos'] = cov_info['covariavel'].map(y_positions)
    total_y_max = y_current
    
    # Margem esquerda de 0.32 garante que rótulos longos não colidam com a barra lateral
    fig, axes = plt.subplots(1, 3, figsize=(18.8, 13.0), sharey=True)
    fig.subplots_adjust(left=0.32, right=0.97, top=0.76, bottom=0.07, wspace=0.10)
    
    desfechos = ['casos', 'obitos', 'letalidade']
    
    for ax_idx, d in enumerate(desfechos):
        ax = axes[ax_idx]
        sub = df[df['desfecho'] == d].copy()
        sub['y_pos'] = sub['covariavel'].map(y_positions)
        
        cfg = CONFIG_DESFECHOS[d]
        
        # 1. Shading de fundo para grupos alternados
        for g_idx, g in enumerate(ORDEM_GRUPOS):
            if g in group_boundaries:
                b = group_boundaries[g]
                fill_col = '#F8FAFC' if g_idx % 2 == 1 else '#FFFFFF'
                ax.axhspan(b['start'] - 0.5, b['end'] + 0.5, color=fill_col, zorder=0, alpha=0.9)
            
        # 2. Linha de base não-informativa (k/K) e Limiar MPM (0,50)
        base_mqo = sub['linha_base_mqo'].iloc[0]
        base_sar = sub['linha_base_sar'].iloc[0]
        base_sem = sub['pip_sem_queen'].mean()
        b_min = min(base_mqo, base_sar, base_sem)
        b_max = max(base_mqo, base_sar, base_sem)
        
        # Baseline band
        ax.axvspan(b_min, b_max, color='#E2E8F0', alpha=0.55, zorder=1)
        ax.axvline(base_sar, color='#94A3B8', linestyle=':', linewidth=1.1, zorder=2)
        
        # Limiar de 0,50 (Median Probability Model)
        ax.axvline(0.50, color='#334155', linestyle='--', linewidth=1.3, zorder=3)
        
        # 3. Guias horizontais sutis para cada linha
        for _, r in sub.iterrows():
            ax.axhline(r['y_pos'], color='#E2E8F0', linewidth=0.5, zorder=1)
            
        # 4. Desenho dos Dumbbells (MQO -> SAR -> SEM)
        for _, r in sub.iterrows():
            y = r['y_pos']
            p_mqo = r['pip_mqo']
            p_sar = r['pip_sar_queen']
            p_sem = r['pip_sem_queen']
            status = r['status_sar']
            cor = CORES_STATUS.get(status, '#94A3B8')
            
            p_min = min(p_mqo, p_sar, p_sem)
            p_max = max(p_mqo, p_sar, p_sem)
            
            is_active = (status != 'Abaixo do limiar (< 0,50)')
            lw = 2.4 if is_active else 1.0
            alpha_line = 0.95 if is_active else 0.45
            
            # Segmento principal do Dumbbell (extensão completa)
            ax.plot([p_min, p_max], [y, y], color=cor, linewidth=lw, alpha=alpha_line, zorder=4, solid_capstyle='round')
            
            # Segmento conectando SAR ao SEM (destaque de modelo espacial)
            ax.plot([p_sar, p_sem], [y, y], color=cor, linewidth=lw*1.15, alpha=alpha_line, zorder=4)
            
            # Marcadores dos 3 Modelos:
            # a) MQO: Círculo aberto com borda da cor de status
            ax.scatter(p_mqo, y, marker='o', facecolor='#FFFFFF', edgecolor=cor, s=65, linewidth=1.8, zorder=5)
            
            # b) SAR (Queen): Diamante sólido preenchido
            ax.scatter(p_sar, y, marker='D', facecolor=cor, edgecolor='#FFFFFF', s=55, linewidth=0.8, zorder=6)
            
            # c) SEM (Queen): Triângulo sólido complementar
            sem_face = cor if is_active else '#64748B'
            ax.scatter(p_sem, y, marker='^', facecolor=sem_face, edgecolor='#FFFFFF', s=45, linewidth=0.7, zorder=6)

        # Configurações dos Eixos
        ax.set_xlim(-0.02, 1.02)
        ax.set_ylim(-0.8, total_y_max + 0.2)
        ax.set_xticks([0.0, 0.25, 0.50, 0.75, 1.0])
        ax.set_xticklabels(['0,00', '0,25', '0,50', '0,75', '1,00'], fontsize=9.5)
        ax.set_xlabel('Probabilidade Posterior de Inclusão (PIP)', fontsize=10.0, fontweight='bold', labelpad=8)
        
        # Título do Painel com diagnóstico integrado
        tit_txt = f"{cfg['titulo']}\nI de Moran y = {cfg['moran_y']:.3f}\n{cfg['subtitulo']}"
        ax.set_title(tit_txt, fontsize=10.8, fontweight='bold', pad=14, color='#0F172A', linespacing=1.25)
                
        # Grade vertical suave
        ax.grid(axis='x', color='#E2E8F0', linestyle='-', linewidth=0.6, zorder=1)
        ax.tick_params(axis='both', which='both', length=3.5, color='#94A3B8')

    # Configuração dos Y-Ticks no painel da esquerda
    axes[0].set_yticks([cov_info.loc[i, 'y_pos'] for i in range(len(cov_info))])
    axes[0].set_yticklabels([cov_info.loc[i, 'rotulo'] for i in range(len(cov_info))], fontsize=9.5, color='#1E293B')
    
    # Rotulos dos grupos tematicos na margem esquerda
    for g in ORDEM_GRUPOS:
        if g in group_boundaries:
            b = group_boundaries[g]
            y_mid = (b['start'] + b['end']) / 2.0
            cor_g = CORES_GRUPO[g]
            
            # Barra lateral colorida de grupo bem à esquerda
            axes[0].plot([-0.56, -0.56], [b['start'] - 0.4, b['end'] + 0.4], 
                         transform=axes[0].get_yaxis_transform(), clip_on=False, 
                         color=cor_g, linewidth=4.0, solid_capstyle='round')
            
            # Nome do grupo temático em destaque rotacionado
            axes[0].text(-0.60, y_mid, g.upper(), transform=axes[0].get_yaxis_transform(),
                         fontsize=8.5, fontweight='bold', color=cor_g, rotation=90, 
                         va='center', ha='right', clip_on=False)
                     
    # Título Geral da Figura
    fig.suptitle('Transição de Probabilidades Posteriores de Inclusão (PIPs): MQO → SAR → SEM',
                 fontsize=15.5, fontweight='bold', color='#0F172A', y=0.97)
    fig.text(0.5, 0.940, 
             'Comparação entre seleção não-espacial (MQO), defasagem espacial (SAR Queen) e erro autorregressivo (SEM Queen) para 30 covariáveis',
             fontsize=10.2, ha='center', color='#475569')

    # Legenda Unificada no Topo (posicionamento sem sobreposição)
    h_mqo = Line2D([0], [0], marker='o', color='w', markeredgecolor='#0F172A', markeredgewidth=1.6, markersize=7.5, label='MQO (Não-espacial)')
    h_sar = Line2D([0], [0], marker='D', color='w', markerfacecolor='#0F172A', markeredgecolor='w', markersize=7.5, label='SAR Queen (Defasagem)')
    h_sem = Line2D([0], [0], marker='^', color='w', markerfacecolor='#0F172A', markeredgecolor='w', markersize=7.5, label='SEM Queen (Erro Espacial)')
    
    h_sai = mpatches.Patch(color=CORES_STATUS['Sai da lista (≥ 0,50 → < 0,50)'], label='Sai da lista (MQO ≥ 0,50 → SAR < 0,50)')
    h_entra = mpatches.Patch(color=CORES_STATUS['Entra na lista (< 0,50 → ≥ 0,50)'], label='Entra na lista (MQO < 0,50 → SAR ≥ 0,50)')
    h_perm = mpatches.Patch(color=CORES_STATUS['Permanece na lista (≥ 0,50 → ≥ 0,50)'], label='Permanece na lista (MQO ≥ 0,50 e SAR ≥ 0,50)')
    h_abaixo = mpatches.Patch(color=CORES_STATUS['Abaixo do limiar (< 0,50)'], label='Abaixo do limiar (< 0,50 em ambos)')
    
    h_lim = Line2D([0], [0], color='#334155', linestyle='--', linewidth=1.3, label='Limiar MPM (0,50)')
    h_base = mpatches.Patch(facecolor='#E2E8F0', edgecolor='#CBD5E1', label='Faixa de não-informação (k/K)')

    leg1 = fig.legend(handles=[h_mqo, h_sar, h_sem, h_lim, h_base],
                      loc='upper center', bbox_to_anchor=(0.645, 0.895),
                      ncol=5, frameon=True, facecolor='#FAFAFA', edgecolor='#E2E8F0',
                      fontsize=8.5, handletextpad=0.5, columnspacing=1.4)
    leg1.get_frame().set_linewidth(0.7)
    
    leg2 = fig.legend(handles=[h_sai, h_entra, h_perm, h_abaixo],
                      loc='upper center', bbox_to_anchor=(0.645, 0.860),
                      ncol=4, frameon=True, facecolor='#FAFAFA', edgecolor='#E2E8F0',
                      fontsize=8.5, handletextpad=0.5, columnspacing=1.4)
    leg2.get_frame().set_linewidth(0.7)

    # Salvar em alta resolução
    png_path = os.path.join(OUT_DIR, 'fig_dumbbell_pips.png')
    pdf_path = os.path.join(OUT_DIR, 'fig_dumbbell_pips.pdf')
    plt.savefig(png_path, dpi=300)
    plt.savefig(pdf_path)
    plt.close()
    print(f"Salvo: {png_path}")
    print(f"Salvo: {pdf_path}")

# -----------------------------------------------------------------------------
# 4. FIGURA 2: Bump Chart de Rankings (Top-10 Covariáveis)
# -----------------------------------------------------------------------------
def gerar_figura_bump_ranking(df):
    print("Gerando Figura 2: Bump Chart de Rankings (Top-10)...")
    
    # Dimensões e espaçamentos equilibrados
    fig, axes = plt.subplots(1, 3, figsize=(19.5, 9.8))
    fig.subplots_adjust(left=0.13, right=0.87, top=0.76, bottom=0.09, wspace=0.64)
    
    desfechos = ['casos', 'obitos', 'letalidade']
    
    for ax_idx, d in enumerate(desfechos):
        ax = axes[ax_idx]
        sub = df[df['desfecho'] == d].copy()
        
        # Calcular rankings exatos (menor rank = maior PIP)
        sub['rank_mqo'] = sub['pip_mqo'].rank(ascending=False, method='min').astype(int)
        sub['rank_sar'] = sub['pip_sar_queen'].rank(ascending=False, method='min').astype(int)
        sub['rank_sem'] = sub['pip_sem_queen'].rank(ascending=False, method='min').astype(int)
        
        # Selecionar variáveis que atingem Top 10 em MQO, SAR ou SEM, ou que tiveram PIP >= 0.50
        mask_top = (sub['rank_mqo'] <= 10) | (sub['rank_sar'] <= 10) | (sub['rank_sem'] <= 10) | (sub['pip_mqo'] >= 0.50)
        top_df = sub[mask_top].copy()
        
        cfg = CONFIG_DESFECHOS[d]
        
        # Contagem MPM (PIP >= 0.50)
        n_mpm_mqo = (sub['pip_mqo'] >= 0.50).sum()
        n_mpm_sar = (sub['pip_sar_queen'] >= 0.50).sum()
        n_mpm_sem = (sub['pip_sem_queen'] >= 0.50).sum()
        
        # Faixas de seleção MPM (PIP >= 0.50) por coluna
        if n_mpm_mqo > 0:
            ax.add_patch(Rectangle((0.68, 0.5), 0.64, n_mpm_mqo, color='#DCFCE7', alpha=0.7, zorder=0))
            ax.text(1.0, 0.22, f"MPM: {n_mpm_mqo} vars", fontsize=8.0, ha='center', color='#166534', fontweight='bold')
        else:
            ax.text(1.0, 0.22, "Nenhuma ≥ 0,50", fontsize=7.8, ha='center', color='#991B1B', fontweight='bold', style='italic')

        if n_mpm_sar > 0:
            ax.add_patch(Rectangle((1.68, 0.5), 0.64, n_mpm_sar, color='#DCFCE7', alpha=0.7, zorder=0))
            ax.text(2.0, 0.22, f"MPM: {n_mpm_sar} vars", fontsize=8.0, ha='center', color='#166534', fontweight='bold')
        else:
            ax.text(2.0, 0.22, "Nenhuma ≥ 0,50", fontsize=7.8, ha='center', color='#991B1B', fontweight='bold', style='italic')
            
        if n_mpm_sem > 0:
            ax.add_patch(Rectangle((2.68, 0.5), 0.64, n_mpm_sem, color='#DCFCE7', alpha=0.7, zorder=0))
            ax.text(3.0, 0.22, f"MPM: {n_mpm_sem} vars", fontsize=8.0, ha='center', color='#166534', fontweight='bold')
        else:
            ax.text(3.0, 0.22, "Nenhuma ≥ 0,50", fontsize=7.8, ha='center', color='#991B1B', fontweight='bold', style='italic')
            
        # Linhas de grade horizontais de rank
        for rk in range(1, 11):
            ax.axhline(rk, color='#F1F5F9', linewidth=0.8, zorder=1)
            
        # Linha limite para além do Top-10
        ax.axhline(11.0, color='#E2E8F0', linestyle=':', linewidth=0.8, zorder=1)

        # Curvas de transição Bezier para cada covariável
        for _, r in top_df.iterrows():
            r1_raw, r2_raw, r3_raw = r['rank_mqo'], r['rank_sar'], r['rank_sem']
            p1, p2, p3 = r['pip_mqo'], r['pip_sar_queen'], r['pip_sem_queen']
            status = r['status_sar']
            
            # Limitar a 11 para visualização limpa
            r1 = min(r1_raw, 11)
            r2 = min(r2_raw, 11)
            r3 = min(r3_raw, 11)
            
            # Destaque de cor conforme o papel no artigo
            if status == 'Entra na lista (< 0,50 → ≥ 0,50)':
                cor = '#0D9488'  # Teal
                lw = 3.2
                alpha = 1.0
                z_order = 10
            elif status == 'Sai da lista (≥ 0,50 → < 0,50)':
                cor = '#DC2626'  # Carmesim
                lw = 2.8
                alpha = 0.95
                z_order = 9
            elif status == 'Permanece na lista (≥ 0,50 → ≥ 0,50)':
                cor = '#16A34A'  # Verde
                lw = 2.8
                alpha = 0.95
                z_order = 8
            else:
                cor = '#94A3B8'  # Cinza ardósia
                lw = 1.1
                alpha = 0.40
                z_order = 3
                
            # Curva suave cúbica spline entre 1 e 2, e 2 e 3
            t = np.linspace(0, 1, 40)
            x_seg1 = (1 - t)**3 * 1.0 + 3 * (1 - t)**2 * t * 1.5 + 3 * (1 - t) * t**2 * 1.5 + t**3 * 2.0
            y_seg1 = (1 - t)**3 * r1 + 3 * (1 - t)**2 * t * r1 + 3 * (1 - t) * t**2 * r2 + t**3 * r2
            
            x_seg2 = (1 - t)**3 * 2.0 + 3 * (1 - t)**2 * t * 2.5 + 3 * (1 - t) * t**2 * 2.5 + t**3 * 3.0
            y_seg2 = (1 - t)**3 * r2 + 3 * (1 - t)**2 * t * r2 + 3 * (1 - t) * t**2 * r3 + t**3 * r3
            
            ax.plot(x_seg1, y_seg1, color=cor, linewidth=lw, alpha=alpha, zorder=z_order)
            ax.plot(x_seg2, y_seg2, color=cor, linewidth=lw, alpha=alpha, zorder=z_order)
            
            # Nós circulares em cada posição (tamanho proporcional à PIP)
            s1 = 45 + p1 * 150
            s2 = 45 + p2 * 150
            s3 = 45 + p3 * 150
            
            ax.scatter(1.0, r1, s=s1, facecolor=cor, edgecolor='#FFFFFF', linewidth=1.2, zorder=z_order+2)
            ax.scatter(2.0, r2, s=s2, facecolor=cor, edgecolor='#FFFFFF', linewidth=1.2, zorder=z_order+2)
            ax.scatter(3.0, r3, s=s3, facecolor=cor, edgecolor='#FFFFFF', linewidth=1.2, zorder=z_order+2)
            
            # Rótulos no eixo esquerdo (MQO) para ranks <= 10: prefixados pelo rank
            if r1_raw <= 10:
                lbl_left = f"#{r1_raw} {r['rotulo']} ({p1:.2f})"
                ax.text(0.90, r1, lbl_left, ha='right', va='center', fontsize=7.8, 
                        fontweight='bold' if p1 >= 0.50 else 'normal',
                        color=cor if alpha >= 0.8 else '#64748B')
                        
            # Rótulos no eixo direito (SEM) para ranks <= 10: prefixados pelo rank
            if r3_raw <= 10:
                lbl_right = f"#{r3_raw} {r['rotulo']} ({p3:.2f})"
                ax.text(3.10, r3, lbl_right, ha='left', va='center', fontsize=7.8,
                        fontweight='bold' if p3 >= 0.50 else 'normal',
                        color=cor if alpha >= 0.8 else '#64748B')
            elif r1_raw == 1 and r3_raw > 10:
                # Se era líder (#1) no MQO e despencou para além do Top-10 no SEM (caso emblemático dos idosos em Casos)
                lbl_out = f"#{r3_raw} {r['rotulo']} ({p3:.2f})"
                ax.text(3.10, 11.0, lbl_out, ha='left', va='center', fontsize=7.5,
                        fontweight='bold', color='#991B1B', style='italic')

        # Configurações de Eixos
        ax.set_xlim(0.38, 3.62)
        ax.set_ylim(11.4, -0.1)  # Invertido: Rank 1 no topo
        ax.set_yticks([])  # Rótulos integrados diretamente com os nomes das covariáveis
        ax.set_xticks([1.0, 2.0, 3.0])
        ax.set_xticklabels(['MQO\n(Não-espacial)', 'SAR Queen\n(Defasagem)', 'SEM Queen\n(Erro Espacial)'], 
                           fontsize=9.8, fontweight='bold')
                           
        # Título do Painel com diagnóstico
        tit_txt = f"{cfg['titulo']}\n{cfg['subtitulo']}"
        ax.set_title(tit_txt, fontsize=11.0, fontweight='bold', pad=14, color='#0F172A', linespacing=1.22)
                     
        ax.tick_params(axis='x', which='both', length=3.5, color='#94A3B8')

    # Título Geral da Figura
    fig.suptitle('Comparação dos Rankings de PIP: MQO, SAR e SEM',
                 fontsize=15.5, fontweight='bold', color='#0F172A', y=0.97)
    fig.text(0.5, 0.938, 
             'Trajetória das covariáveis do Top-10 entre os modelos. O tamanho dos círculos é proporcional à PIP. Faixa verde: Zona de Inclusão MPM (PIP ≥ 0,50).',
             fontsize=10.2, ha='center', color='#475569')

    # Legenda Superior (posicionamento sem sobreposição)
    h_sai = Line2D([0], [0], color=CORES_STATUS['Sai da lista (≥ 0,50 → < 0,50)'], linewidth=3.0, label='Sai da lista (MQO ≥ 0,50 → SAR < 0,50)')
    h_entra = Line2D([0], [0], color=CORES_STATUS['Entra na lista (< 0,50 → ≥ 0,50)'], linewidth=3.0, label='Entra na lista (MQO < 0,50 → SAR ≥ 0,50)')
    h_perm = Line2D([0], [0], color=CORES_STATUS['Permanece na lista (≥ 0,50 → ≥ 0,50)'], linewidth=3.0, label='Permanece na lista (MQO ≥ 0,50 e SAR ≥ 0,50)')
    h_abaixo = Line2D([0], [0], color=CORES_STATUS['Abaixo do limiar (< 0,50)'], linewidth=1.5, label='Abaixo do limiar (< 0,50 em ambos)')
    h_zone = mpatches.Patch(facecolor='#DCFCE7', edgecolor='#86EFAC', label='Zona de Inclusão MPM (PIP ≥ 0,50)')

    fig.legend(handles=[h_sai, h_entra, h_perm, h_abaixo, h_zone],
               loc='upper center', bbox_to_anchor=(0.50, 0.880),
               ncol=5, frameon=True, facecolor='#FAFAFA', edgecolor='#E2E8F0',
               fontsize=8.8, handletextpad=0.6, columnspacing=1.4)

    # Salvar em alta resolução
    png_path = os.path.join(OUT_DIR, 'fig_bump_ranking.png')
    pdf_path = os.path.join(OUT_DIR, 'fig_bump_ranking.pdf')
    plt.savefig(png_path, dpi=300)
    plt.savefig(pdf_path)
    plt.close()
    print(f"Salvo: {png_path}")
    print(f"Salvo: {pdf_path}")

# -----------------------------------------------------------------------------
# 5. Execução Principal
# -----------------------------------------------------------------------------
if __name__ == '__main__':
    print("===============================================================")
    print("Iniciando geração de figuras para o artigo MC3 Espacial...")
    print("===============================================================")
    df = carregar_dados()
    gerar_figura_dumbbell(df)
    gerar_figura_bump_ranking(df)
    print("===============================================================")
    print("Processo concluído com sucesso!")
    print("===============================================================")

