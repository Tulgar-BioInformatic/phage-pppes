#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
gpumap.py
---------
But : charger un ou plusieurs fichiers de vecteurs (TSV "ID  v0 v1 v2 ..."),
ne garder que les IDs listés dans un fichier texte, puis calculer une
projection 2D avec UMAP (implémentation GPU de cuML). Sauvegarde les
coordonnées dans un TSV "ID  x  y". Optionnellement, affiche un scatterplot.

Entrées (exemples) :
  -i Partie_3_Embedding/output/xxx.tsv           (répéter -i pour plusieurs sources)
  -IDs Partie_4_UMAP/input/IDs.txt               (1 ID par ligne)
  -o Partie_4_UMAP/output/umap.tsv               (fichier TSV de sortie)
  -n 15 --min-dist 0.01 --metric cosine --n-epochs 5000 --random-state 42

Sortie :
  - TSV avec colonnes : ID, x, y
"""

import argparse
from pathlib import Path
import numpy as np
import pandas as pd
from cuml.manifold import UMAP          # UMAP GPU (RAPIDS cuML)
import matplotlib.pyplot as plt

# ---------------------------------------------------------------------
# 1) Arguments CLI
# ---------------------------------------------------------------------
def parse_args():
    p = argparse.ArgumentParser(description="Exécute UMAP-GPU et retourne les coordonnées.")
    p.add_argument("-i", "--input", action="append", required=True,
                   help="Chemin(s) vers fichier(s) de vecteurs TSV (répéter -i si plusieurs).")
    p.add_argument("-IDs", "--inputIDs", required=True,
                   help="Chemin vers fichier texte listant les IDs à retenir (1 par ligne).")
    p.add_argument("-n", "--n-neighbors", type=int, default=150,
                   help="Nombre de voisins UMAP (ex. 15, 50, 150...).")
    p.add_argument("--min-dist", type=float, default=0.01, help="Paramètre min_dist d'UMAP.")
    p.add_argument("--metric", default="cosine", help="Métrique (cosine, euclidean, ...).")
    p.add_argument("--random-state", type=int, default=42, help="Graine aléatoire.")
    p.add_argument("--n-epochs", type=int, default=5000, help="Nombre d’époques d’optimisation.")
    p.add_argument("-p", "--plot", action="store_true", help="Afficher le scatterplot (matplotlib).")
    p.add_argument("-o", "--output", type=Path, required=True,
                   help="Fichier TSV de sortie (sera créé si nécessaire).")
    return p.parse_args()

# ---------------------------------------------------------------------
# 2) Utilitaires d’E/S
# ---------------------------------------------------------------------
def load_ids_list(path: str) -> list[str]:
    """Lit un fichier texte et renvoie la liste d'IDs (1 par ligne)."""
    ids = []
    with open(path) as f:
        for line in f:
            ids.append(line.strip())
    return ids

def load_single_source(path: str, keep_ids: list[str]) -> tuple[np.ndarray, np.ndarray]:
    """
    Lit un TSV de vecteurs. Formats acceptés :
      ID\t"v0 v1 v2 ..."  (une seule colonne vecteur en chaîne)
      ID\tv0\tv1\tv2 ...  (colonnes séparées)
    Ne conserve que les lignes dont l'ID est dans keep_ids.
    Renvoie X (np.ndarray) et ids (np.ndarray).
    """
    data, ids = [], []
    keep_set = set(keep_ids)
    with open(path) as f:
        for line in f:
            parts = line.rstrip().split("\t")
            if not parts:
                continue
            uid = parts[0]
            if uid not in keep_set:
                continue
            vec_parts = parts[1:]
            if len(vec_parts) == 1 and " " in vec_parts[0]:
                arr = np.fromstring(vec_parts[0], sep=" ", dtype=float)
            else:
                arr = np.array(vec_parts, dtype=float)
            ids.append(uid)
            data.append(arr)
    if not data:
        return np.empty((0, 0)), np.empty((0,), dtype=object)
    return np.vstack(data), np.array(ids, dtype=object)

def load_multiple_sources(paths: list[str], keep_ids: list[str]) -> tuple[np.ndarray, np.ndarray]:
    """
    Concatène verticalement les points provenant de plusieurs sources
    (les IDs sont répétés si la même protéine a plusieurs représentations).
    """
    X_list, id_list = [], []
    for p in paths:
        X, ids = load_single_source(p, keep_ids)
        if X.size == 0:
            continue
        X_list.append(X)
        id_list.append(ids)
    if not X_list:
        return np.empty((0, 0)), np.empty((0,), dtype=object)
    return np.vstack(X_list), np.concatenate(id_list)

def plot_scatter(df: pd.DataFrame) -> None:
    """Affiche un nuage de points (x,y)."""
    plt.figure(figsize=(9, 9))
    plt.scatter(df["x"], df["y"], s=8, alpha=0.7)
    plt.xlabel("UMAP 1")
    plt.ylabel("UMAP 2")
    plt.title("Projection UMAP (GPU)")
    plt.grid(True, linewidth=0.3)
    plt.tight_layout()
    plt.show()

# ---------------------------------------------------------------------
# 3) Main
# ---------------------------------------------------------------------
def main():
    args = parse_args()

    # 3.1 Charger la liste d’IDs à conserver
    keep_ids = load_ids_list(args.inputIDs)

    # 3.2 Charger/concaténer les vecteurs depuis 1+ fichiers
    X, ids = load_multiple_sources(args.input, keep_ids)
    if X.size == 0:
        raise SystemExit("Aucun vecteur chargé (vérifie -i et -IDs).")

    # 3.3 Instancier UMAP (cuML)
    umap = UMAP(
        n_neighbors=args.n_neighbors,
        min_dist=args.min_dist,
        metric=args.metric,
        random_state=args.random_state,
        n_epochs=args.n_epochs,
        verbose=True,
        build_algo="brute_force_knn",   # explicite pour reproductibilité
    )

    # 3.4 Apprendre et projeter
    coords = umap.fit_transform(X)     # shape: (n_points, 2)

    # 3.5 Préparer sortie
    out_path = Path(args.output).expanduser().resolve()
    out_path.parent.mkdir(parents=True, exist_ok=True)

    df = pd.DataFrame({
        "ID": ids,                     # on garde l’ID original tel quel
        "x": coords[:, 0],
        "y": coords[:, 1]
    })
    df.to_csv(out_path, sep="\t", index=False)

    # 3.6 Plot optionnel
    if args.plot:
        plot_scatter(df)

if __name__ == "__main__":
    main()
