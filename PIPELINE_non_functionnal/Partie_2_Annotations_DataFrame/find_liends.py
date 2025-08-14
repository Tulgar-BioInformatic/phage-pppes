#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
find_liends.py
--------------
Ce script lit le résultat d’un alignement entre ORFs et protéines
(`result_getorf_vs_prodigal_FE.m8`) et construit un tableau indiquant,
pour chaque ORF :
    - l’ID de la protéine correspondante (si trouvée dans le même génome),
    - si cette correspondance est bien dans le même génome (booléen 0/1),
    - le nombre total de protéines alignées (tous génomes confondus).

Sortie : fichier TSV `liens_orf_prot.tsv`
"""

from collections import defaultdict
import sys

# Dictionnaires pour stocker les correspondances
dorfProt     = defaultdict(str)   # ORF → protéine correspondante (même génome)
dorfNbProt   = defaultdict(int)   # ORF → nombre total de protéines alignées
dorfProtBool = defaultdict(int)   # ORF → 1 si match même génome & identité > 95%
sOrf         = set()               # Ensemble des ORFs rencontrés

# Lecture du fichier m8 (tabulé : query, subject, %identity, ...)
with open("Partie_2_Annotations_DataFrame/tmp/getorf_vs_prodigal/result_getorf_vs_prodigal_FE.m8",'r') as f1:
    for lig in f1:
        li  = lig.rstrip("\n").split("\t")
        orf = li[0]
        prot = li[1]
        v1  = "_".join(orf.split("_")[:-2])   # ID du génome de l’ORF
        v2  = "_".join(prot.split("_")[:-2])  # ID du génome de la protéine
        idp = float(li[2])                    # identité en %

        # Si c’est le même génome et identité > 95% → correspondance principale
        if (v1 == v2) and idp > 0.95: 
            dorfProt[orf] = prot
            dorfProtBool[orf] = 1
            sOrf.add(orf)

        # Incrémente le nombre total de protéines matchées par cet ORF
        dorfNbProt[orf] += 1
        sOrf.add(orf)

# Écriture du fichier TSV avec en-tête
with open("Partie_2_Annotations_DataFrame/tmp/liens_orf_prot.tsv", "w") as fe1:
    print("", "protID", "matchedProtSameGenome", "nMatchesOnProtOthersGenomes", sep="\t", file=fe1)

    # Chaque ligne = ORF + infos associées
    for o in sOrf:
        print(o, dorfProt[o], dorfProtBool[o], dorfNbProt[o], sep="\t", file=fe1)
