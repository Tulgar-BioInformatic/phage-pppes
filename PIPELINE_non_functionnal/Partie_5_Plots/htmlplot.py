#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
htmlplot.py
-----------
But : lire des coordonnées UMAP (TSV) + des métadonnées (annotations.tsv),
puis générer une page HTML Plotly interactive avec :
  - 4×2 sous-graphiques (brin, clusters, outil, MOG, ratio, histogrammes, longueur),
  - sliders pour la taille et l'opacité des points,
  - légende “intelligente” (clic = masquer/afficher, double-clic = isoler),
  - synchro des zooms entre sous-graphiques,
  - export SVG et export des IDs sélectionnés (lasso/box).

Usage :
  python htmlplot.py -i Partie_4_UMAP/output/umap.tsv -o umap_interactive_search.html
"""

import argparse, numpy as np, pandas as pd, matplotlib.cm as cm
from pathlib import Path
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import plotly.io as pio

# ────────────────────────── CLI ──────────────────────────
ap = argparse.ArgumentParser()
ap.add_argument("-i","--input",  type=Path, required=True, help="TSV avec colonnes ID, x, y")
ap.add_argument("-o","--output", type=Path, default="umap_interactive_search.html", help="Fichier HTML de sortie")
args = ap.parse_args()

# ──────────── Données UMAP + métadonnées ────────────
umap_df = pd.read_csv(args.input, sep="\t")  # attend des colonnes ID, x, y

# annotations.tsv (index = ID des séquences)
meta_df = pd.read_csv(
    "Partie_2_Annotations_DataFrame/output/annotations.tsv",
    sep="\t", index_col=0
)
meta_dic = meta_df.to_dict("index")

# Récupérer quelques colonnes utiles des annotations (si absentes → NaN)
for col in ["strand","ratio","predictionTool","clusterID","mogID","mogNumber","mogText","length"]:
    umap_df[col] = umap_df["ID"].map(lambda i: meta_dic.get(i, {}).get(col, np.nan))

# Valeurs catégorielles : remplacer NaN par "NA" pour l’affichage
for col in ["strand","predictionTool","clusterID","mogNumber"]:
    umap_df[col] = umap_df[col].fillna("NA")

# Valeurs numériques : forcer en numérique + défaut si NaN
umap_df["ratio"]   = pd.to_numeric(umap_df["ratio"], errors="coerce").fillna(0)
umap_df["length"]  = pd.to_numeric(umap_df["length"], errors="coerce")  # on laisse NaN si inconnu

# ─────────── Palettes de couleurs ───────────
# brin / outil : palettes fixes
strand2col = {"-":"blue", "+":"red", "NA":"green"}
tool2col   = {"Getorf":"purple","Prodigal":"orange","Expe":"green"}

# mogNumber / clusterID : palettes automatiques (tab20)
base = cm.get_cmap("tab20").colors
mog2col = {m: f"rgb{tuple(int(c*255) for c in base[i%20])}"
           for i, m in enumerate(sorted(x for x in umap_df["mogNumber"].unique() if x!="NA"))}
mog2col["NA"] = "gray"
clu2col = {c: f"rgb{tuple(int(c*255) for c in base[i%20])}"
           for i, c in enumerate(sorted(x for x in umap_df["clusterID"].unique() if x!="NA"))}
clu2col["NA"] = "gray"

# Paramètres par défaut des points
PS_DEFAULT      = 14    # taille
OPACITY_DEFAULT = 0.8   # opacité

# ─────────── Helpers pour ajouter des traces ───────────
def add_by_cat(fig, df, col, palette, row, col_nb, prefix="", lg_title=None, order=None):
    """Ajouter des nuages par catégorie (ex. brin, outil, cluster, MOG)."""
    cats = order if order else palette.keys()
    for cat in cats:
        if cat not in palette:
            continue
        sub = df[df[col] == cat]
        if sub.empty:
            continue
        fig.add_trace(go.Scatter(
            x=sub["x"], y=sub["y"], mode="markers",
            marker=dict(color=palette[cat], size=PS_DEFAULT, opacity=OPACITY_DEFAULT),
            name=f"{prefix}{cat} ({len(sub)})",
            legendgroup=lg_title or col,
            legendgrouptitle=dict(text=lg_title) if lg_title else None,
            showlegend=True,
            customdata=sub[["ID","mogNumber","mogText"]].values,
            hovertemplate="<b>%{customdata[0]}</b><br>MOG #: %{customdata[1]}<br>%{customdata[2]}<extra></extra>"
        ), row=row, col=col_nb)

def add_ratio(fig, df, row, col_nb):
    """Nuage coloré par le ratio (Viridis)."""
    fig.add_trace(go.Scatter(
        x=df["x"], y=df["y"], mode="markers",
        marker=dict(color=df["ratio"], colorscale="Viridis",
                    size=PS_DEFAULT, opacity=OPACITY_DEFAULT, showscale=False),
        hoverinfo="skip", showlegend=False,
        customdata=df[["ID","mogNumber","mogText"]].values,
        hovertemplate="<b>%{customdata[0]}</b><br>MOG #: %{customdata[1]}<br>%{customdata[2]}<extra></extra>"
    ), row=row, col=col_nb)

def add_length(fig, df, row, col_nb, colorscale="Turbo"):
    """Nuage coloré par la longueur (si disponible)."""
    sub = df.copy()
    cds = sub[["ID","mogNumber","mogText","length"]].values  # pour le hover
    fig.add_trace(go.Scatter(
        x=sub["x"], y=sub["y"], mode="markers",
        marker=dict(color=sub["length"], colorscale=colorscale,
                    size=PS_DEFAULT, opacity=OPACITY_DEFAULT, showscale=False),
        hoverinfo="skip", showlegend=False,
        customdata=cds,
        hovertemplate="<b>%{customdata[0]}</b><br>MOG #: %{customdata[1]}<br>%{customdata[2]}<br>Longueur: %{customdata[3]}<extra></extra>"
    ), row=row, col=col_nb)

# ─────────── Figure avec 4 lignes × 2 colonnes ───────────
fig = make_subplots(
    rows=4, cols=2,
    specs=[[{}, {}], [{}, {}], [{}, {}], [{}, {}]],
    subplot_titles=(
        "Brin", "Clusters",
        "Outil de prédiction", "MOG ID",
        "Conservation / Ratio", "Histogramme ratio",
        "Longueur (nt)", "Histogramme longueur"
    ),
    horizontal_spacing=0.03, vertical_spacing=0.05
)

# Ligne 1 : brin, clusters
add_by_cat(fig, umap_df, "strand",    strand2col, 1, 1, prefix="Brin ", lg_title="Brin")
add_by_cat(fig, umap_df, "clusterID", clu2col,    1, 2, prefix="Cluster ", lg_title="Clusters")

# Ligne 2 : outil, MOG
tool_order = ["Getorf","Prodigal","Expe"]
add_by_cat(fig, umap_df, "predictionTool", tool2col, 2, 1, lg_title="Outil de prédiction", order=tool_order)
mog_vals  = [m for m in sorted(umap_df["mogNumber"].unique()) if m != "NA"]
mog_order = ["NA"] + mog_vals
add_by_cat(fig, umap_df, "mogNumber", mog2col, 2, 2, prefix="MOG ", lg_title="MOG ID", order=mog_order)

# Ligne 3 : ratio (nuage + histo)
add_ratio(fig, umap_df, 3, 1)

# Ligne 4 : longueur (nuage + histo)
add_length(fig, umap_df, 4, 1)

# ==== Paramètres pour les “bandes d’échelle” colorées sous les histogrammes ====
BAND_OFFSET_FRAC    = 0.12  # distance sous l’axe X (en fraction de la hauteur max)
BAND_THICKNESS_FRAC = 0.03  # épaisseur de la bande (fraction hauteur max)

# ---------- Histogramme ratio (3,2) + bande dégradée ----------
cnt, bins = np.histogram(umap_df["ratio"].values, bins=50)
cent = 0.5 * (bins[:-1] + bins[1:])
ymax = cnt.max() if len(cnt) else 1

offset    = BAND_OFFSET_FRAC    * max(1, ymax)
thickness = BAND_THICKNESS_FRAC * max(1, ymax)
band_top, band_bottom = -offset, -(offset + thickness)

# Bande de couleur (Heatmap) placée sous l’axe (y négatif)
fig.add_trace(
    go.Heatmap(
        x=cent, y=[band_bottom, band_top],
        z=np.tile(cent, (2, 1)),
        colorscale="Viridis", showscale=False,
        hoverinfo="skip", xgap=0, ygap=0
    ),
    row=3, col=2
)
# Barres de l’histogramme
fig.add_trace(
    go.Bar(x=cent, y=cnt, marker=dict(color=cent, colorscale="Viridis"), showlegend=False),
    row=3, col=2
)
# Étendre l’axe Y pour voir la bande sous 0
fig.update_yaxes(range=[band_bottom, ymax*1.05], row=3, col=2)

# ---------- Histogramme longueurs (4,2) + bande dégradée ----------
len_vals = pd.to_numeric(umap_df["length"], errors="coerce").dropna().values
cntL, binsL = np.histogram(len_vals, bins=50)
centL = 0.5 * (binsL[:-1] + binsL[1:])
ymaxL = cntL.max() if len(cntL) else 1

offsetL    = BAND_OFFSET_FRAC    * max(1, ymaxL)
thicknessL = BAND_THICKNESS_FRAC * max(1, ymaxL)
band_topL, band_bottomL = -offsetL, -(offsetL + thicknessL)

fig.add_trace(
    go.Heatmap(
        x=centL, y=[band_bottomL, band_topL],
        z=np.tile(centL, (2, 1)),
        colorscale="Turbo", showscale=False,
        hoverinfo="skip", xgap=0, ygap=0
    ),
    row=4, col=2
)
fig.add_trace(
    go.Bar(x=centL, y=cntL, marker=dict(color=centL, colorscale="Turbo"), showlegend=False),
    row=4, col=2
)
fig.update_yaxes(range=[band_bottomL, ymaxL*1.05], row=4, col=2)

# Mise en page générale : grande légende à droite, drag = zoom
fig.update_layout(
    height=1950,
    margin=dict(t=40, b=30, l=30, r=300),
    dragmode="zoom",
    legend=dict(
        x=1.02, xanchor="left", y=0.5, yanchor="middle",
        bgcolor="rgba(255,255,255,0.65)", bordercolor="black", borderwidth=1,
        itemsizing="constant"
    )
)

# ─────────── Sliders (taille & opacité) ───────────
# On repère les traces Scatter (les histos/heatmaps ne doivent pas bouger)
scatter_idxs = [i for i, tr in enumerate(fig.data) if isinstance(tr, go.Scatter)]

size_slider = dict(
    active=(PS_DEFAULT - 4)//2,
    currentvalue={"prefix": "Taille : "},
    pad={"t": 40},
    steps=[
        dict(
            label=str(sz),
            method="restyle",
            args=[{"marker.size": [sz] * len(scatter_idxs)}, scatter_idxs]
        )
        for sz in range(4, 21, 2)
    ]
)
opacity_slider = dict(
    active=int((OPACITY_DEFAULT - 0.2)/0.1),
    pad={"t": 90},
    currentvalue={"prefix": "Opacité : "},
    steps=[
        dict(
            label=f"{op:.1f}",
            method="restyle",
            args=[{"marker.opacity": [op] * len(scatter_idxs)}, scatter_idxs]
        )
        for op in np.arange(0.2, 1.01, 0.1)
    ]
)
fig.update_layout(sliders=[size_slider, opacity_slider])

# ─────────── JavaScript embarqué pour l’interactivité avancée ───────────
#  - Toolbar perso (boutons Export SVG / Export IDs sélectionnés)
#  - Légende : clic (toggle), double-clic (isoler), clic droit (amener au-dessus)
#  - Synchronisation des zooms/pans entre plusieurs sous-graphiques
post_js = r"""
document.addEventListener('DOMContentLoaded', () => {
  const gd = document.getElementById('plotly-graph');

  /* ───────── Toolbar propre au-dessus du graphe ───────── */
  const container = gd.parentElement;
  const bar = document.createElement('div');
  bar.className = 'umap-toolbar';
  bar.style.cssText = [
    'position:sticky','top:0',
    'display:flex','gap:8px','align-items:center','flex-wrap:wrap',
    'margin:0 0 8px 0','padding:6px 8px',
    'border:1px solid #ddd','border-radius:8px',
    'background:rgba(255,255,255,0.9)','backdrop-filter:blur(2px)',
    'box-shadow:0 1px 6px rgba(0,0,0,.08)','z-index:10'
  ].join(';');
  container.insertBefore(bar, gd);

  const mkBtn = (label) => {
    const b = document.createElement('button');
    b.textContent = label;
    b.style.cssText = [
      'padding:6px 10px','border:1px solid #ccc','border-radius:8px',
      'background:#f7f7f7','cursor:pointer','font-size:13px'
    ].join(';');
    b.onmouseenter = () => (b.style.background = '#eee');
    b.onmouseleave = () => (b.style.background = '#f7f7f7');
    bar.appendChild(b);
    return b;
  };

  const exportSvgBtn = mkBtn('Exporter SVG');
  const exportIdsBtn = mkBtn('Exporter IDs sélectionnés');
  exportIdsBtn.disabled = true;

  /* ───────── Export SVG groupé par sous-plot et légende ───────── */
  exportSvgBtn.onclick = () => {
    const svgOrig = gd.querySelector('svg');
    if(!svgOrig){ alert('SVG introuvable'); return; }
    const svg = svgOrig.cloneNode(true);

    svg.querySelectorAll('g.subplot').forEach((sp,i)=>{
      const wrap = document.createElementNS(svg.namespaceURI,'g');
      wrap.setAttribute('id',`subplot_${i+1}`);
      sp.parentNode.insertBefore(wrap, sp);
      wrap.appendChild(sp);
    });
    svg.querySelectorAll('g.legend').forEach((leg,i)=>{
      const wrap = document.createElementNS(svg.namespaceURI,'g');
      wrap.setAttribute('id',`legend_${i+1}`);
      leg.parentNode.insertBefore(wrap, leg);
      wrap.appendChild(leg);
    });

    const serializer = new XMLSerializer();
    const svgText = serializer.serializeToString(svg);
    const blob = new Blob([svgText], {type:'image/svg+xml'});
    const url  = URL.createObjectURL(blob);
    const a    = document.createElement('a');
    a.href = url; a.download = 'umap_grouped.svg'; a.click();
    setTimeout(()=>URL.revokeObjectURL(url),100000);
  };

  /* ───────── Sélection lasso/box → export d’IDs ───────── */
  let selIDs = [];
  gd.on('plotly_selected', ev => {
    if(!ev || !ev.points){ selIDs = []; exportIdsBtn.disabled = true; return; }
    selIDs = [...new Set(ev.points.map(pt=>{
      const cd = gd.data[pt.curveNumber].customdata;
      return cd ? cd[pt.pointNumber][0] : null;
    }).filter(Boolean))];
    exportIdsBtn.disabled = selIDs.length === 0;
  });

  exportIdsBtn.onclick = () => {
    if(!selIDs.length){ alert('Aucun point sélectionné'); return; }
    const blob = new Blob([selIDs.join('\n')], {type:'text/plain'});
    const url  = URL.createObjectURL(blob);
    const a    = document.createElement('a');
    a.href = url; a.download = 'selected_ids.txt'; a.click();
    setTimeout(()=>URL.revokeObjectURL(url),1000);
  };

  /* ───────── Légende : clic / double-clic / clic droit ───────── */
  gd.on('plotly_legendclick',ev=>{
    const grp  = gd.data[ev.curveNumber].legendgroup;
    const idxs = gd.data.map((t,i)=>t.legendgroup===grp?i:-1).filter(i=>i>=0);
    const cur  = ev.curveNumber, clicks = ev.event.detail;
    if(clicks===2){
      const visCount = idxs.filter(i=>gd.data[i].visible===true||gd.data[i].visible===undefined).length;
      if(visCount===1&&(gd.data[cur].visible===true||gd.data[cur].visible===undefined))
        Plotly.restyle(gd,{visible:true},idxs);
      else{
        const vis = idxs.map(i=>i===cur?true:'legendonly');
        Plotly.restyle(gd,{visible:vis},idxs);
      }
      return false;
    }
    if(clicks===1){
      const curVis = gd.data[cur].visible;
      Plotly.restyle(gd,{visible:(curVis==='legendonly')?true:'legendonly'},[cur]);
      return false;
    }
  });

  gd.once('plotly_afterplot',()=>{
    gd.querySelectorAll('.legendtoggle').forEach(el=>{
      el.addEventListener('contextmenu',e=>{
        e.preventDefault();
        const curve = el.__data__.trace.index;
        Plotly.moveTraces(gd,curve,gd.data.length-1);
        if(gd.data[curve].visible==='legendonly')
          Plotly.restyle(gd,{visible:true},[curve]);
      });
    });
  });

  /* zoom pan sync */
  const xa=['xaxis','xaxis2','xaxis3','xaxis4','xaxis5','xaxis7'],
        ya=['yaxis','yaxis2','yaxis3','yaxis4','yaxis5','yaxis7'];

  let syncing=false;
  gd.on('plotly_relayout',ev=>{
    if(syncing) return;
    const xr0=ev[Object.keys(ev).find(k=>/^xaxis\d*\.range\[0]/.test(k))];
    if(xr0===undefined) return;
    const xr1=ev[Object.keys(ev).find(k=>/^xaxis\d*\.range\[1]/.test(k))],
          yr0=ev[Object.keys(ev).find(k=>/^yaxis\d*\.range\[0]/.test(k))],
          yr1=ev[Object.keys(ev).find(k=>/^yaxis\d*\.range\[1]/.test(k))];
    const up={}; xa.forEach(a=>up[`${a}.range`]=[xr0,xr1]);
                 ya.forEach(a=>up[`${a}.range`]=[yr0,yr1]);
    syncing=true; Plotly.relayout(gd,up).then(()=>{syncing=false;});
  });
});
"""

# Config Plotly (logo masqué, responsive, outils lasso/box ajoutés)
config = {
    "displaylogo": False,
    "responsive": True,
    "modeBarButtonsToAdd": ["lasso2d", "select2d"]
}

# Générer le fragment HTML du graphe (avec JS embarqué)
html = pio.to_html(
    fig, full_html=False, include_plotlyjs="embed",
    div_id="plotly-graph", post_script=post_js, config=config
)

# Un peu de CSS pour la page
page_css = """
<style>
  html,body{height:100%; margin:0; background:#f5f7fb;}
  .modebar{ gap:8px; }
  .modebar-btn{ padding:10px !important; }
  .modebar-btn svg{ width:28px !important; height:28px !important; }
  .umap-toolbar button{ padding:10px 14px; font-size:14px; }
</style>
"""

# Page complète (doctype minimal + encodage)
full_html = (
    "<!DOCTYPE html><html><head><meta charset='utf-8'>"
    + page_css +
    "</head><body>" + html + "</body></html>"
)

# Écriture du fichier final
Path(args.output).write_text(full_html, encoding="utf8")
print(f"✅ Fichier HTML créé : {args.output}")
