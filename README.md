# Phage-PPPES : Phage Protein Prediction by Projection of Embedded Sequences

Ce projet explore l’utilisation de **modèles de langage protéiques** (ProtT5) couplés à une **projection UMAP** pour discriminer les protéines codantes des séquences non codantes dans les génomes de **Microviridae**.

Les **UMAPs interactifs** sont consultables ici :  
🔗 **[Visualisations interactives](https://tulgar-bioinformatic.github.io/phage-pppes/)**

## Table des principaux outils utilisés

| Catégorie | Outil | Version | Paramètres / Remarques | Lien |
|-----------|-------|---------|------------------------|------|
| **Prédiction et recherche de gènes** | [Getorf](http://emboss.sourceforge.net/apps/release/6.6/emboss/apps/getorf.html) (suite EMBOSS) | - | Détection des ORFs à partir de codons start/stop | [Site officiel](http://emboss.sourceforge.net/) |
| | [Prodigal](https://github.com/hyattpd/Prodigal) | 2.6.3 | `-p single` ou `-p meta` | [GitHub](https://github.com/hyattpd/Prodigal) |
| | [MMseqs2](https://github.com/soedinglab/MMseqs2) | - | `easy-search`, `easy-cluster` | [GitHub](https://github.com/soedinglab/MMseqs2) |
| **Modèles de langage et embeddings** | [ProtT5-XL-UniRef50](https://github.com/agemagician/ProtTrans) | - | Modèle de langage protéique pour vectorisation des séquences | [GitHub](https://github.com/agemagician/ProtTrans) |
| **Réduction de dimension et clustering** | [UMAP-learn](https://umap-learn.readthedocs.io/) | - | `n_neighbors` / `min_dist` ajustés selon la vue | [Documentation](https://umap-learn.readthedocs.io/) |
| | [cuML-UMAP](https://docs.rapids.ai/api/cuml/stable/) | - | Implémentation GPU de UMAP, `metric=cosine`, `brute_force_knn` | [RAPIDS cuML](https://docs.rapids.ai/api/cuml/stable/) |
| | [t-SNE](https://scikit-learn.org/stable/modules/generated/sklearn.manifold.TSNE.html) | - | Mentionné en perspective | [scikit-learn](https://scikit-learn.org/stable/) |
| | [MMseqs2](https://github.com/soedinglab/MMseqs2) | - | Utilisé aussi pour clustering et recherche | [GitHub](https://github.com/soedinglab/MMseqs2) |
| **Outils et langages généraux** | [Python](https://www.python.org/) | - | Langage principal pour scripts et analyses | [Site officiel](https://www.python.org/) |
| | [Conda](https://docs.conda.io/) | - | Gestion des environnements (`env_cpu.yml`, `env_gpu.yml`) | [Documentation](https://docs.conda.io/) |
| | GPU NVIDIA RTX 4070 | - | Accélération pour ProtT5 & UMAP cuML | [Fiche technique](https://www.nvidia.com/) |

## Structure du dépôt

- `env_cpu.yml` : environnement Conda CPU (analyses légères)
- `env_gpu.yml` : environnement Conda GPU (embedding ProtT5 + UMAP cuML)
- `scripts/` : scripts Python pour le prétraitement, l’embedding, et la visualisation
- `data/` : exemples de données d’entrée (ou liens vers les données)
- `umap_global.html` et `umap_central.html` : visualisations interactives exportées
- `README.md` : ce fichier

## Référence

Rapport de stage M1 Bioinformatique : *Exploitation des modèles de langage protéique et de la projection UMAP pour améliorer la prédiction de gènes codants chez les phages : application aux Microviridae*  
Auteur : **Youlen Iglesias** – Encadrant : **François Enault** – 2025
