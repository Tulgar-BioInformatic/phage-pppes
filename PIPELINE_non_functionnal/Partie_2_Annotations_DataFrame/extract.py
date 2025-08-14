#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
extract.py
----------
But : extraire du FASTA `allseqs.faa` les séquences dont l'ID figure dans
la table d'annotations et qui respectent ces critères :
  - pas de protéine correspondante détectée dans le même génome
    (matchedProtSameGenome est NaN),
  - ET (le génome est le représentant de son cluster) OU (la source est "Ref")
  - ET séquence extraite d'un génome appartenant à un cluster de taille supérieure à 60.

Sortie : écrit un FASTA `extraction_clusterRep_n0_or_ref.faa` avec les séquences retenues.
"""

import pandas as pd

# --- 1) Charger la table d’annotations comme dict --------------------------------
df_annotations = pd.read_csv(
    "Partie_2_Annotations_DataFrame/output/annotations.tsv",
    sep="\t",
    index_col=0,
    low_memory=False,
)
dict_annotations = df_annotations.to_dict(orient="index")

# --- 2) Charger le FASTA en dictionnaire {seqID: sequence} -----------------------
def read_fasta_to_dict(path_faa: str) -> dict[str, str]:
    seqs = {}
    header = None
    chunks = []
    with open(path_faa) as f:
        for line in f:
            if line.startswith(">"):
                # on sauvegarde l'entrée précédente si elle existe
                if header is not None:
                    seqs[header] = "".join(chunks)
                # nouvel en-tête (ID = premier mot après '>')
                header = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line.strip())
        # ne pas oublier la dernière séquence
        if header is not None:
            seqs[header] = "".join(chunks)
    return seqs

dict_allseqs = read_fasta_to_dict("Partie_2_Annotations_DataFrame/input/allseqs.faa")

# --- 3) Sélectionner les IDs selon les critères ----------------------------------
ids = [
    seq_id for seq_id, v in dict_annotations.items()
    
    if (
        pd.isna(v.get("matchedProtSameGenome"))                             # pas de match même génome
        and (
            v["nRelatives"] >= 0 )
        and (
            (v.get("genomeID") == v.get("clusterID"))                       # génome représentant du cluster
            or (v.get("source") == "Ref")                                   # ou séquence de référence
        )
    )
]

# --- 4) Écrire les séquences retenues dans un nouveau FASTA ----------------------
outPath = "Partie_2_Annotations_DataFrame/output/extraction_clusterRep_n0_or_ref.faa"

with open(outPath, "w") as f:
    for seq_id in ids:
        seq = dict_allseqs.get(seq_id)
        if not seq:
            # ID absent du FASTA → on ignore proprement
            continue
        # écrire en FASTA (en-tête minimal + séquence)
        f.write(f">{seq_id}\n{seq}\n")
