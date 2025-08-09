# Phage-PPPES : Phage Protein Prediction by Projection of Embedded Sequences

Ce projet explore l’utilisation de **modèles de langage protéiques** (ProtT5) couplés à une **projection UMAP** pour discriminer les protéines codantes des séquences non codantes dans les génomes de **Microviridae**.

Les **UMAPs interactifs** sont consultables ici :  
🔗 **[Visualisations interactives](https://tulgar-bioinformatic.github.io/phage-pppes/)**

## Table des principaux outils utilisés

| Outil | Version | Paramètres principaux | Fonction |
|-------|---------|-----------------------|----------|
| [singularity](https://github.com/apptainer/singularity) | 4.1.1 | `exec --bind` | Exécuter des outils via des containers/images |
| [Prodigal](https://github.com/hyattpd/Prodigal) | 2.6.3 | `-p single` ou `-p meta` | Prédiction des gènes codants |
| [MMseqs2](https://github.com/soedinglab/MMseqs2) | v15-6f452+ds-2 | `easy-search` | Recherche de similarité entre séquences (~BLAST) |
| [UMAP-learn](https://umap-learn.readthedocs.io/) | 0.5.7 | `n_neighbors` / `min_dist` ajustés | Réduction de dimension pour la visualisation |
| [cuML-UMAP (RAPIDS)](https://docs.rapids.ai/api/cuml/stable/) | 25.04 | `brute_force_knn` / `metric=cosine` | UMAP accéléré GPU |
| [Biopython](https://biopython.org/) | 1.85 | - | Manipulation des séquences biologiques |
| [Plotly](https://plotly.com/python/) | 6.1.2 | - | Visualisation interactive des résultats |

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
