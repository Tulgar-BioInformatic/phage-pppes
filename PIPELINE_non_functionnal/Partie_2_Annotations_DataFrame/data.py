#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
data.py
-------
But : construire la table d'annotations finale (annotations.tsv) en combinant :
  • les comptages de similarités « all » et « near » (df_match_vs_all / df_match_vs_near),
  • les liens ORF/prot du même génome (liens_orf_prot.tsv),
  • l’appartenance des génomes à des clusters,
  • les méta-infos encodées dans l’ID des séquences (source, outil de prédiction, brin),
  • les en-têtes FASTA (pour garder une trace).

Entrées attendues :
  - tmp/df_match_vs_near.tsv
  - tmp/df_match_vs_all.tsv
  - tmp/liens_orf_prot.tsv
  - tmp/clusterized_capsid_prodigal_allseqs_cluster.tsv
  - tmp/conservations/input/allseqs.faa

Sortie :
  - output/annotations.tsv
"""

from collections import defaultdict, Counter
import re
import pandas as pd

# --------- petites aides ------------------------------------------------------

def suffix_abbrev(abbr: str) -> str:
    """Traduire les abréviations trouvées dans l'ID."""
    mapping = {
        "Prod": "Prodigal",
        "Geto": "Getorf",
        "Posi": "+",
        "Nega": "-",
    }
    return mapping.get(abbr, abbr)  # si inconnu, on renvoie tel quel

# --------- charger les tables intermédiaires ---------------------------------

df_match_vs_near = pd.read_csv(
    "Partie_2_Annotations_DataFrame/tmp/df_match_vs_near.tsv",
    sep="\t", index_col=0
)
dict_match_vs_near = df_match_vs_near.to_dict(orient="index")

df_match_vs_all = pd.read_csv(
    "Partie_2_Annotations_DataFrame/tmp/df_match_vs_all.tsv",
    sep="\t", index_col=0
)
dict_match_vs_all = df_match_vs_all.to_dict(orient="index")

df_liens_orf_prot = pd.read_csv(
    "Partie_2_Annotations_DataFrame/tmp/liens_orf_prot.tsv",
    sep="\t", index_col=0
)
dict_liens_orf_prot = df_liens_orf_prot.to_dict(orient="index")

# --------- dictionnaire cluster (genomeID -> clusterID) -----------------------

cluster_dict: dict[str, str] = {}
with open("Partie_2_Annotations_DataFrame/tmp/clusterized_capsid_prodigal_allseqs_cluster.tsv") as f:
    for line in f:
        ls = line.split()
        # On retire le suffixe "_<nombre>_<Mot>" présent en fin d'ID
        genomeID  = re.sub(r"_[\d]+_[a-zA-Z]+$", "", ls[1])
        clusterID = re.sub(r"_[\d]+_[a-zA-Z]+$", "", ls[0])
        cluster_dict[genomeID] = clusterID

# Taille de chaque cluster (pour nRelatives)
cluster_sizes = Counter(cluster_dict.values())

# --------- parcourir le FASTA pour remplir la table --------------------------

rows = defaultdict(dict)  # seqID -> dict(col -> valeur)

with open("Partie_2_Annotations_DataFrame/tmp/conservations/input/allseqs.faa") as f:
    for header in f:
        if not header.startswith(">"):
            continue

        # ID brut (premier mot après '>')
        seqID = header.split()[0].lstrip(">")

        # Extraire les suffixes encodés (ex. "..._RefProdPosi")
        suffixe  = seqID.rsplit("_", 1)[-1]
        suffixes = re.findall(r"[A-Z][a-z]+", suffixe)  # ex. ["Ref", "Prod", "Posi"]

        # Retirer le suffixe final pour retrouver genomeID
        genomeID = re.sub(r"_[\d]+_[a-zA-Z]+$", "", seqID)

        # Champs de base
        rows[seqID]["fastaHeader"] = header.rstrip("\n")
        rows[seqID]["genomeID"] = genomeID
        rows[seqID]["clusterID"] = cluster_dict.get(genomeID)

        # Comptages « relatives / matchs »
        cid = rows[seqID]["clusterID"]
        rows[seqID]["nRelatives"] = cluster_sizes.get(cid, 0)

        # nMatchs parmi les génomes « proches » (même cluster)
        if seqID in dict_match_vs_near:
            rows[seqID]["nMatchsToRelatives"] = dict_match_vs_near[seqID].get("nGenomes", 0)
        else:
            rows[seqID]["nMatchsToRelatives"] = 0

        # nMatchs parmi « all »
        if seqID in dict_match_vs_all:
            rows[seqID]["nMatchsToAll"] = dict_match_vs_all[seqID].get("nMatchs", 0)
        else:
            rows[seqID]["nMatchsToAll"] = 0

        # ratio = nMatchsToRelatives / nRelatives (si possible)
        nrel = rows[seqID]["nRelatives"]
        rows[seqID]["ratio"] = (rows[seqID]["nMatchsToRelatives"] / nrel) if nrel else 0

        # Décoder la source / outil / brin à partir des suffixes
        # On tolère les cas incomplets (longueur < 3)
        src  = suffixes[0] if len(suffixes) >= 1 else None
        tool = suffix_abbrev(suffixes[1]) if len(suffixes) >= 2 else None
        strand = suffix_abbrev(suffixes[2]) if len(suffixes) >= 3 else None

        rows[seqID]["source"] = src
        rows[seqID]["predictionTool"] = tool
        rows[seqID]["strand"] = strand

        # Ajouter les infos ORF↔prot si présentes
        if seqID in dict_liens_orf_prot:
            r = dict_liens_orf_prot[seqID]
            rows[seqID]["matchedProtSameGenome"] = r.get("matchedProtSameGenome")
            rows[seqID]["nMatchesOnProtOthersGenomes"] = r.get("nMatchesOnProtOthersGenomes")
            rows[seqID]["protMatchedID"] = r.get("protID")
        else:
            rows[seqID]["matchedProtSameGenome"] = None
            rows[seqID]["nMatchesOnProtOthersGenomes"] = 0
            rows[seqID]["protMatchedID"] = None

# --------- DataFrame et export -----------------------------------------------

df = pd.DataFrame.from_dict(rows, orient="index")
df.to_csv("Partie_2_Annotations_DataFrame/output/annotations.tsv", sep="\t")
