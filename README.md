# Phage-PPPES : Phage Protein Prediction by Projection of Embedded Sequences

Ce projet explore l’utilisation de **modèles de langage protéiques** (ProtT5) couplés à une **projection UMAP** pour discriminer les séquences codantes des séquences non codantes dans les génomes de **Microviridae**.

Les **UMAPs interactifs** sont consultables ici :  
 **[Visualisations interactives](https://tulgar-bioinformatic.github.io/phage-pppes/)**

## Outils utilisés

L’ensemble des principaux outils et bibliothèques utilisés dans ce projet est listé ci-dessous.  
Chaque entrée inclut la version (si connue), les paramètres principaux employés, ainsi qu’une brève description de sa fonction dans le pipeline.

| OUTIL | VERSION | PARAMÈTRES | FONCTION |
|-------|---------|------------|----------|
| [MEGAHIT](https://github.com/voutcn/megahit) | 1.2.9 | `--presets meta-large` | Assemblage de génomes à partir de métagénomes |
| [Getorf](http://emboss.sourceforge.net/apps/release/6.6/emboss/apps/getorf.html) (suite EMBOSS) | - | Codons start alternatifs, circulaire, taille min 25 aa | Détection des ORFs à partir de codons start/stop |
| [Prodigal](https://github.com/hyattpd/Prodigal) | 2.6.3 | `-p single` (génomes entiers) / `-p meta` (contigs courts) | Prédiction des gènes codants |
| [MMseqs2](https://github.com/soedinglab/MMseqs2) | - | `easy-search`, `easy-cluster`, `--min-seq-id 0.95 -c 0.95` (retrait redondances) | Recherche de similarité et clustering de séquences |
| [RNAcode](https://github.com/s-will/RNAcode) | - | Alignements multiples + détection sélection négative | Détection de régions codantes par analyse évolutive |
| [ProtT5-XL-UniRef50](https://github.com/agemagician/ProtTrans) | - | - | Vectorisation (embedding) de séquences protéiques |
| [UMAP-learn](https://umap-learn.readthedocs.io/) | - | `n_neighbors`, `min_dist` ajustés selon vue globale/locale | Réduction de dimension |
| [cuML-UMAP](https://docs.rapids.ai/api/cuml/stable/) | 25.04 | `metric=cosine`, `brute_force_knn`, `n_epochs=5000` | UMAP accéléré sur GPU |
| [DNAViewer](https://dnacanvas.com/) | - | - | Visualisation de cartes génomiques |
| [Plotly](https://plotly.com/python/) | 6.1.2 | - | Visualisation interactive des résultats |
| [Python](https://www.python.org/) | - | - | Langage principal pour scripts et analyses |
| [Conda](https://docs.conda.io/) | - | - | Gestion des environnements (`env_cpu.yml`, `env_gpu.yml`) |
| GPU NVIDIA RTX 4070 | - | - | Accélération des calculs ProtT5 et UMAP cuML |

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
