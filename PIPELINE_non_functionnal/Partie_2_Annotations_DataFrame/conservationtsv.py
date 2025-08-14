#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
conservationtsv.py
------------------
But : à partir d’un fichier d’alignements (allseqs_vs_self.m8), calculer pour
chaque séquence :
  - le nombre total de matchs ("all") et le nombre de génomes distincts touchés,
  - le nombre de matchs restreints au même cluster ("near") et le nombre de
    génomes distincts touchés dans ce cluster,
  - la taille (nb de génomes) du cluster de la séquence (nRelatives).

Sorties :
  - tmp/df_match_vs_all.tsv   (colonnes : nMatchs, nRelatives, nGenomes)
  - tmp/df_match_vs_near.tsv  (idem, mais restreint au même cluster)
"""

from collections import defaultdict, Counter
import re
import pandas as pd

# --- 1) Charger le mapping genomeID -> clusterID ------------------------------
cluster_dict = {}
with open("Partie_2_Annotations_DataFrame/tmp/clusterized_capsid_prodigal_allseqs_cluster.tsv") as f:
    for l in f:
        ls = l.split()
        genomeID  = re.sub(r"_[\d]+_[a-zA-Z]+$", "", ls[1])  # nettoie le suffixe de fin
        clusterID = re.sub(r"_[\d]+_[a-zA-Z]+$", "", ls[0])
        cluster_dict[genomeID] = clusterID

# Taille de chaque cluster : clusterID -> nb de génomes
cluster_sizes = Counter(cluster_dict.values())

# --- 2) Préparer compteurs et ensembles d'accumulation ------------------------
# Compteurs par séquence
all_counts  = defaultdict(lambda: {"nMatchs": 0, "nRelatives": 0})
near_counts = defaultdict(lambda: {"nMatchs": 0, "nRelatives": 0})

# Ensembles des génomes atteints (pour compter des génomes distincts)
all_genome_sets  = defaultdict(set)   # seqID -> {genomeID, ...}
near_genome_sets = defaultdict(set)

# --- 3) Parcourir les alignements --------------------------------------------
with open("Partie_2_Annotations_DataFrame/tmp/conservations/out/allseqs_vs_self.m8") as f:
    for l in f:
        ls = l.split()
        query, target = ls[0], ls[1]

        # On ignore l'auto-alignement
        if query == target:
            continue

        # Récupérer les genomeID bruts (sans suffixe final)
        queryGenomeID  = re.sub(r"_[\d]+_[a-zA-Z]+$", "", query)
        targetGenomeID = re.sub(r"_[\d]+_[a-zA-Z]+$", "", target)

        # --- Comptage "all" : tous les matchs, tous génomes confondus
        all_counts[query]["nMatchs"]  += 1
        all_counts[target]["nMatchs"] += 1
        all_genome_sets[query].add(targetGenomeID)
        all_genome_sets[target].add(queryGenomeID)

        # --- Comptage "near" : uniquement si même cluster
        if cluster_dict.get(queryGenomeID) == cluster_dict.get(targetGenomeID):
            near_counts[query]["nMatchs"]  += 1
            near_counts[target]["nMatchs"] += 1
            near_genome_sets[query].add(targetGenomeID)
            near_genome_sets[target].add(queryGenomeID)

# --- 4) Ajouter nRelatives et nGenomes pour chaque séquence -------------------
for seqID in all_counts:
    currentGenomeID  = re.sub(r"_[\d]+_[a-zA-Z]+$", "", seqID)
    currentClusterID = cluster_dict[currentGenomeID]
    all_counts[seqID]["nRelatives"] = cluster_sizes[currentClusterID]      # taille du cluster
    all_counts[seqID]["nGenomes"]   = len(all_genome_sets.get(seqID, set()))

for seqID in near_counts:
    currentGenomeID  = re.sub(r"_[\d]+_[a-zA-Z]+$", "", seqID)
    currentClusterID = cluster_dict[currentGenomeID]
    near_counts[seqID]["nRelatives"] = cluster_sizes[currentClusterID]
    near_counts[seqID]["nGenomes"]   = len(near_genome_sets.get(seqID, set()))

# --- 5) Sauvegarde en TSV -----------------------------------------------------
df_all  = pd.DataFrame.from_dict(all_counts,  orient="index")
df_near = pd.DataFrame.from_dict(near_counts, orient="index")

df_all.to_csv("Partie_2_Annotations_DataFrame/tmp/df_match_vs_all.tsv",  sep="\t")
df_near.to_csv("Partie_2_Annotations_DataFrame/tmp/df_match_vs_near.tsv", sep="\t")
