#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
MOGs.py
-------
But : enrichir la table `annotations.tsv` avec les informations MOG
(Microvirus Orthologous Groups).

Entrées :
  - output/annotations.tsv                                  (table à enrichir)
  - input/all_annotated_microvirus_proteins.faa             (définitions MOG)
      >MOG_ID MOG_NUMBER MOG_DESCRIPTION...
  - input/allseqs_vs_mogs.m8                                (meilleurs matches)
      <seqID>\t<MOG_ID>\t...

Sortie :
  - output/annotations.tsv  (réécrit avec colonnes : mogID, mogNumber, mogText)
"""

from pathlib import Path
import pandas as pd

# --- 1) Charger la table d’annotations en dict --------------------------------
data_df = pd.read_csv(
    "Partie_2_Annotations_DataFrame/output/annotations.tsv",
    sep="\t",
    index_col=0
)
data_dict = data_df.to_dict("index")

# --- 2) Lire les définitions MOG depuis le FASTA annoté -----------------------
#     Chaque en-tête : >MOG_ID MOG_NUMBER MOG_DESCRIPTION...
mogNumber_dict = {}
mogText_dict = {}
with open("Partie_2_Annotations_DataFrame/input/all_annotated_microvirus_proteins.faa") as f:
    for line in f:
        if not line.startswith(">"):
            continue
        header = line.strip(">").strip()
        parts = header.split()
        if len(parts) < 2:
            # si l'en-tête n'a pas au moins ID et NUMBER, on ignore
            continue
        mogID = parts[0]
        mogNumber = parts[1]
        mogText = " ".join(parts[2:]) if len(parts) > 2 else ""
        mogNumber_dict[mogID] = mogNumber
        mogText_dict[mogID] = mogText

# --- 3) Lire les correspondances seqID -> MOG_ID (résultat d'alignement) ------
match_dict = {}
with open("Partie_2_Annotations_DataFrame/input/allseqs_vs_mogs.m8") as f:
    for line in f:
        parts = line.split()
        if len(parts) < 2:
            continue
        seqID, mogID = parts[0], parts[1]
        # on garde le meilleur match déjà préparé dans le .m8 (supposé trié)
        if seqID not in match_dict:
            match_dict[seqID] = mogID

# --- 4) Injecter les infos MOG dans la table ----------------------------------
for seqID, mogID in match_dict.items():
    if seqID not in data_dict:
        continue  # si la séquence n'est pas dans annotations.tsv
    data_dict[seqID]["mogID"] = mogID
    data_dict[seqID]["mogNumber"] = mogNumber_dict.get(mogID)
    data_dict[seqID]["mogText"] = mogText_dict.get(mogID)

# --- 5) Réécrire annotations.tsv ----------------------------------------------
out_df = pd.DataFrame.from_dict(data_dict, orient="index")
out_df.to_csv("Partie_2_Annotations_DataFrame/output/annotations.tsv", sep="\t")
