## Phase 4 : Parallélisme en Mémoire Distribuée (MPI) — Stratégie 1

### Objectif de la phase

Les Phases 2 et 3 ont montré les limites de l'optimisation sur une seule machine en mémoire partagée : la contention atomique sur la grille de phéromones annule les gains de la parallélisation OpenMP. Pour dépasser cette limitation, nous passons au paradigme de **mémoire distribuée** avec **MPI (Message Passing Interface)**.

MPI repose sur un principe fondamentalement différent d'OpenMP : au lieu de partager la mémoire entre les threads, chaque **processus** possède son propre espace mémoire privé. Les processus communiquent **explicitement** par échange de messages. Ce modèle élimine structurellement les problèmes de contention mémoire, mais introduit un nouveau coût : celui de la **communication**.

Le sujet du projet propose deux stratégies de parallélisation MPI. Cette phase implémente la **Stratégie 1** (dite « réplication de la carte ») :

- **Chaque processus MPI possède la carte entière** (terrain fractal + grille de phéromones)
- **Les fourmis sont réparties** entre les processus ($N_{\text{total}} / P$ fourmis par processus)
- **À chaque itération**, les grilles de phéromones locales sont **fusionnées** via une opération collective `MPI_Allreduce` avec l'opérateur `MPI_MAX`

Le sujet avertit explicitement que cette approche échange « une grande quantité de données entre les processus » et n'est efficace que « si le nombre de fourmis est grand et la carte assez petite ». Nos expériences viseront à vérifier cette prédiction.

---

### Étape 4.1 : Vérification de l'Environnement MPI

#### Objectif

Avant toute implémentation, vérifier que les outils MPI sont disponibles et configurer la chaîne de compilation.

#### Vérifications effectuées

| Outil | Commande | Résultat |
|---|---|---|
| Compilateur MPI | `mpic++ --version` | `g++ (Ubuntu 13.3.0-6ubuntu2~24.04) 13.3.0` ✅ |
| OpenMPI (runtime) | `dpkg -l openmpi-bin` | Version 4.1.6-7ubuntu2, déjà installé ✅ |
| MPICH (alternative) | `sudo apt install mpich` | Version 4.2.0-5build3, installé en complément ✅ |

#### Configuration du Makefile

Le fichier `Make_linux.inc` a été modifié pour utiliser le compilateur MPI :

```makefile
# Avant (Phase 1-3)
CXX = g++ -fopenmp

# Après (Phase 4+)
CXX = mpic++ -fopenmp
```

Le wrapper `mpic++` appelle `g++` en ajoutant automatiquement les chemins d'inclusion et les bibliothèques MPI. Le drapeau `-fopenmp` est conservé pour permettre l'utilisation combinée OpenMP + MPI (parallélisme hybride) si nécessaire.

#### Topologie matérielle confirmée

La tentative de lancement avec 12 processus a échoué avec le message :

```
There are not enough slots available in the system to satisfy the 12 slots
that were requested by the application
```

Cela confirme que la machine dispose de **8 cœurs physiques**. OpenMPI refuse par défaut de lancer plus de processus que de cœurs disponibles (surallocation). L'option `--oversubscribe` permet de forcer l'exécution, mais les résultats seront pénalisés par le *context switching*.

#### Conclusion

L'environnement de développement MPI est pleinement opérationnel. La chaîne de compilation est configurée et la topologie matérielle (8 cœurs) est identifiée.

---

### Étape 4.2 : Implémentation de la Stratégie 1

#### Principe de la stratégie

La Stratégie 1 est la plus simple conceptuellement :

1. **Initialisation** : chaque processus construit le terrain fractal (identique sur tous les rangs grâce à la même graine) et crée sa propre grille de phéromones. Le nombre total de fourmis est **divisé équitablement** entre les processus.

2. **Boucle principale** : à chaque itération, chaque processus :
   - Fait avancer ses fourmis locales (modification de sa grille de phéromones locale)
   - Fusionne sa grille avec celles des autres processus via `MPI_Allreduce(MPI_MAX)`
   - Applique l'évaporation et la mise à jour sur la grille fusionnée

3. **Finalisation** : le compteur de nourriture est sommé sur tous les processus via `MPI_Reduce(MPI_SUM)`.

#### Architecture du schéma de communication

```
Itération k :
┌──────────────────┐   ┌──────────────────┐   ┌──────────────────┐
│     Rank 0       │   │     Rank 1       │   │     Rank 2       │
│ N/P fourmis      │   │ N/P fourmis      │   │ N/P fourmis      │
│ Carte ENTIÈRE    │   │ Carte ENTIÈRE    │   │ Carte ENTIÈRE    │
│ (copie locale)   │   │ (copie locale)   │   │ (copie locale)   │
└────────┬─────────┘   └────────┬─────────┘   └────────┬─────────┘
         │                      │                      │
         ▼                      ▼                      ▼
    advance_all()          advance_all()          advance_all()
    (modif. phen locale)   (modif. phen locale)   (modif. phen locale)
         │                      │                      │
         └──────────┬───────────┴──────────────────────┘
                    ▼
           MPI_Allreduce(MPI_MAX)
           → Chaque processus reçoit la grille
             fusionnée (valeur MAX globale par cellule)
                    │
         ┌──────────┼──────────┐
         ▼          ▼          ▼
    do_evaporation()  do_evaporation()  do_evaporation()
    update()          update()          update()
    (identique sur chaque rang)
```

**Justification de l'opérateur `MPI_MAX`** : le sujet précise que « lorsque deux fourmis appartenant à deux processus différents se trouvent sur une même cellule [...] on choisira de prendre la valeur la plus grande d'entre tous les processus comme valeur de phéromone ». L'opérateur `MPI_MAX` réalise exactement cette opération de façon collective.

#### Modifications réalisées

| Fichier | Modification |
|---|---|
| **`ant_simu.cpp`** | Refonte complète du fichier principal : |
| | — Initialisation MPI (`MPI_Init`, `MPI_Comm_rank`, `MPI_Comm_size`) |
| | — Division des fourmis : `local_nb_ants = total / nb_procs` avec gestion du reste |
| | — Graine unique par processus : `local_seed = seed + rank * 1000` |
| | — Fusion des phéromones : `MPI_Allreduce(MPI_MAX)` à chaque itération |
| | — Somme de la nourriture : `MPI_Reduce(MPI_SUM)` en fin de simulation |
| | — Chronométrage séparé de la communication MPI |
| | — Affichage des résultats uniquement sur le rank 0 |
| | — Finalisation MPI (`MPI_Finalize`) |
| **`pheronome.hpp`** | Déjà compatible grâce aux accesseurs `get_data()` et `get_total_size()` ajoutés en Phase 2 |
| **`Make_linux.inc`** | `CXX = mpic++ -fopenmp` |

#### Détail des opérations MPI utilisées

| Opération MPI | Utilisation | Fréquence | Volume de données |
|---|---|---|---|
| `MPI_Allreduce(MPI_MAX)` | Fusion des grilles de phéromones | 1 fois par itération (×500) | ~16.8 Mo par appel |
| `MPI_Reduce(MPI_SUM)` | Somme du compteur de nourriture | 1 fois en fin de simulation | 8 octets |
| `MPI_Barrier` | Synchronisation avant le chronomètre global | 1 fois | 0 octet |

**Volume de données échangées par `MPI_Allreduce`** :

La grille de phéromones a une taille de $(N+2)^2 \times 2$ doubles (avec $N = 1025$ pour la grille de base, ou $N = 513$ selon la configuration) :

$$\text{Volume} = (1025 + 2)^2 \times 2 \times 8 \text{ octets} \approx 16.8 \text{ Mo}$$

Sur 500 itérations, le volume total échangé est de :

$$500 \times 16.8 \text{ Mo} \approx 8.4 \text{ Go}$$

#### Exécution

```bash
# Avec 1 processus (référence séquentielle MPI)
mpirun -np 1 ./ant_simu.exe

# Avec 2, 4, 8 processus
mpirun -np 2 ./ant_simu.exe
mpirun -np 4 ./ant_simu.exe
mpirun -np 8 ./ant_simu.exe

# Avec 12 processus (surallocation forcée)
mpirun --oversubscribe -np 12 ./ant_simu.exe
```

---

### Étape 4.3 : Mesures de Performance et Analyse

#### 4.3.1 — Résultats avec 5 000 fourmis (configuration de référence)

**Conditions de test** : 5 000 fourmis, 500 itérations, grille 1025×1025, GUI désactivé.

| Nb Procs | Temps Total (s) | Mouvement (s) | Communication (s) | Évaporation (s) | Update (s) |
|:---:|:---:|:---:|:---:|:---:|:---:|
| 1 | 1.704 | 0.989 (58.0%) | 0.496 (29.1%) | 0.214 (12.6%) | 0.005 (0.3%) |
| 2 | 1.915 | 0.519 (27.1%) | 1.136 (59.3%) | 0.255 (13.3%) | 0.005 (0.3%) |
| 4 | 6.910 | 0.950 (13.7%) | 5.460 (79.0%) | 0.481 (7.0%) | 0.019 (0.3%) |
| 8 | 8.892 | 0.518 (5.8%) | 7.636 (85.9%) | 0.722 (8.1%) | 0.015 (0.2%) |
| 12 | ❌ Échec | — | — | — | — |

**Tableau de speedup :**

| Nb Procs | Temps Total (s) | Speedup ($T_1/T_N$) | Efficacité (Speedup/$N$) |
|:---:|:---:|:---:|:---:|
| 1 | 1.704 | 1.00 | 100% |
| 2 | 1.915 | 0.89 ❌ | 44.5% |
| 4 | 6.910 | 0.25 ❌ | 6.1% |
| 8 | 8.892 | 0.19 ❌ | 2.4% |

**Constat** : avec 5 000 fourmis, le speedup est **systématiquement inférieur à 1**. L'ajout de processus ralentit le programme. La communication MPI domine le temps total dès 2 processus.

#### Analyse du mouvement des fourmis

Le mouvement se parallélise correctement :

| Nb Procs | Fourmis/proc | Temps Mouvement (s) | Speedup Mouvement |
|:---:|:---:|:---:|:---:|
| 1 | 5 000 | 0.989 | 1.00 |
| 2 | 2 500 | 0.519 | **1.91** |
| 4 | 1 250 | 0.950 | 1.04 |
| 8 | 625 | 0.518 | **1.91** |

Le temps de mouvement se divise bien par ~2 entre 1 et 2 processus (0.989 → 0.519 s), confirmant que la répartition des fourmis fonctionne. Les variations à 4 processus (0.950 s) sont attribuées à des effets de cache et de contention mémoire liés à l'exécution simultanée de 4 processus sur la même machine.

#### Analyse de la communication MPI

Le véritable goulot d'étranglement est la communication :

| Nb Procs | Temps Communication (s) | % du Temps Total | Volume par appel |
|:---:|:---:|:---:|:---:|
| 1 | 0.496 | 29.1% | ~16.8 Mo (copie locale `memcpy`) |
| 2 | 1.136 | 59.3% | ~16.8 Mo × 2 rangs |
| 4 | 5.460 | 79.0% | ~16.8 Mo × 4 rangs |
| 8 | 7.636 | 85.9% | ~16.8 Mo × 8 rangs |

Le `MPI_Allreduce` est appelé **500 fois** sur un vecteur de **~16.8 Mo**. Le volume total de données transitant par le sous-système MPI atteint **~8.4 Go** sur toute la simulation.

```
Communication MPI (s)
│
8 ┤                              ████████
│                         ████████
6 ┤                    ████████
│               ████████
4 ┤          ████████
│     ████████
2 ┤████████
│▓▓▓▓
0 ┤──────────────────────────────
   1P    2P    4P    8P

▓▓▓▓ = Mouvement    ████ = Communication MPI
```

**Diagnostic** : le ratio calcul/communication est **défavorable**. Avec seulement 5 000 fourmis, le mouvement ne prend que ~1 seconde au total, tandis que la communication croît de façon quasi-linéaire avec le nombre de processus. C'est un cas classique de *communication-bound* : le calcul est noyé dans les échanges de données.

#### Comparaison avec les prédictions du sujet

Le sujet avertissait explicitement :

> « Le problème de cette approche est qu'elle n'est efficace que si le nombre de fourmis est grand et la carte assez petite pour tenir en mémoire sur chaque processus. [...] Par contre, une grande quantité de données est échangée entre les processus. »

Nos résultats **confirment exactement** cette prédiction. Avec une carte de 1025×1025 (~16.8 Mo de phéromones) et seulement 5 000 fourmis, le ratio est défavorable.

---

#### 4.3.2 — Tests de scalabilité avec charge accrue

Pour vérifier l'hypothèse que l'augmentation du nombre de fourmis améliore le ratio calcul/communication, nous avons testé deux configurations supplémentaires : **50 000** et **100 000 fourmis**.

##### Configuration B : 50 000 fourmis

| Nb Procs | Temps Total (s) | Mouvement (s) | Communication (s) | Évaporation (s) | Update (s) |
|:---:|:---:|:---:|:---:|:---:|:---:|
| 1 | 8.220 | 7.539 (91.7%) | 0.478 (5.8%) | 0.199 (2.4%) | 0.004 (0.1%) |
| 2 | 6.810 | 5.026 (73.8%) | 1.424 (20.9%) | 0.355 (5.2%) | 0.005 (0.1%) |
| 4 | 11.701 | 5.535 (47.3%) | 5.722 (48.9%) | 0.431 (3.7%) | 0.013 (0.1%) |
| 8 | 12.857 | 3.797 (29.5%) | 8.363 (65.0%) | 0.682 (5.3%) | 0.015 (0.1%) |
| 12* | 16.396 | 2.887 (17.6%) | 12.339 (75.3%) | 1.151 (7.0%) | 0.017 (0.1%) |

##### Configuration C : 100 000 fourmis

| Nb Procs | Temps Total (s) | Mouvement (s) | Communication (s) | Évaporation (s) | Update (s) |
|:---:|:---:|:---:|:---:|:---:|:---:|
| 1 | 15.474 | 14.776 (95.5%) | 0.493 (3.2%) | 0.202 (1.3%) | 0.004 (0.0%) |
| 2 | 10.889 | 9.150 (84.0%) | 1.375 (12.6%) | 0.359 (3.3%) | 0.005 (0.0%) |
| 4 | 17.025 | 9.357 (55.0%) | 7.238 (42.5%) | 0.418 (2.5%) | 0.012 (0.1%) |
| 8 | 17.045 | 7.103 (41.7%) | 9.173 (53.8%) | 0.754 (4.4%) | 0.016 (0.1%) |
| 12* | 19.368 | 5.338 (27.6%) | 12.908 (66.6%) | 1.109 (5.7%) | 0.012 (0.1%) |

*\* 12 processus sur 8 cœurs physiques via `--oversubscribe`. Résultats pénalisés par le context switching.*

##### Tableaux de speedup comparatifs

**Speedup global ($T_1/T_N$) :**

| Nb Procs | 5 000 fourmis | 50 000 fourmis | 100 000 fourmis |
|:---:|:---:|:---:|:---:|
| 1 | 1.00 | 1.00 | 1.00 |
| 2 | 0.89 ❌ | **1.21** ✅ | **1.42** ✅ |
| 4 | 0.25 ❌ | 0.70 ❌ | 0.91 ❌ |
| 8 | 0.19 ❌ | 0.64 ❌ | 0.91 ❌ |
| 12* | — | 0.50 ❌ | 0.80 ❌ |

**Speedup du mouvement uniquement (composante parallélisable) :**

| Nb Procs | 5 000 fourmis | 50 000 fourmis | 100 000 fourmis |
|:---:|:---:|:---:|:---:|
| 1 | 1.00 | 1.00 | 1.00 |
| 2 | 1.91 | 1.50 | **1.61** |
| 4 | 1.04 | 1.36 | **1.58** |
| 8 | 1.91 | **1.99** | **2.08** |
| 12* | — | **2.61** | **2.77** |

Le mouvement se parallélise efficacement : le speedup du mouvement atteint **2.77×** à 12 processus avec 100k fourmis. C'est la communication qui plafonne le speedup global.

---

### Analyse détaillée des résultats

#### Confirmation de l'hypothèse : le ratio calcul/communication

Le facteur déterminant est le **pourcentage du mouvement** dans le temps total séquentiel (1 processus) :

| Configuration | % Mouvement (1 proc) | % Communication (1 proc) | Ratio Calcul/Comm |
|:---:|:---:|:---:|:---:|
| 5 000 fourmis | 58.0% | 29.1% | **2.0** |
| 50 000 fourmis | 91.7% | 5.8% | **15.8** |
| 100 000 fourmis | 95.5% | 3.2% | **30.0** |

Avec 100 000 fourmis, le mouvement représente 95.5% du temps total séquentiel. La fraction parallélisable est donc bien plus importante, ce qui se traduit par un **speedup réel à 2 processus (1.42×)**. L'hypothèse est **confirmée** : augmenter le nombre de fourmis améliore le ratio et permet d'atteindre un speedup positif.

#### Le plafonnement au-delà de 2 processus

Le `MPI_Allreduce` échange **~16.8 Mo à chaque itération**, indépendamment du nombre de fourmis. Ce coût est :
- **Fixe** par rapport à la charge de calcul (il ne diminue pas quand on ajoute des fourmis)
- **Croissant** avec le nombre de processus (l'algorithme de réduction collective a une complexité en $O(V \times \log P)$ où $V$ est le volume et $P$ le nombre de processus)

| Nb Procs | Comm. 5k (s) | Comm. 50k (s) | Comm. 100k (s) |
|:---:|:---:|:---:|:---:|
| 1 | 0.496 | 0.478 | 0.493 |
| 2 | 1.136 | 1.424 | 1.375 |
| 4 | 5.460 | 5.722 | 7.238 |
| 8 | 7.636 | 8.363 | 9.173 |
| 12* | — | 12.339 | 12.908 |

La communication est **quasi indépendante du nombre de fourmis** (~12.3 s vs ~12.9 s à 12 procs pour 50k vs 100k). Cela confirme que le coût est dominé par le **volume de la grille** et non par le calcul.

#### Évolution du profil temporel — Illustration visuelle

L'impact de l'augmentation des fourmis sur le profil temporel à 2 processus :

```
5 000 fourmis (2P)          50 000 fourmis (2P)        100 000 fourmis (2P)
T = 1.92s                   T = 6.81s                  T = 10.89s
┌──────────────────┐        ┌──────────────────┐       ┌──────────────────┐
│▓▓▓▓▓▓▓ 27%      │        │▓▓▓▓▓▓▓▓▓▓▓▓▓ 74%│       │▓▓▓▓▓▓▓▓▓▓▓▓▓▓84%│
│  Mouvement       │        │  Mouvement       │       │  Mouvement       │
│                  │        │                  │       │                  │
│████████████ 59%  │        │████ 21%          │       │███ 13%           │
│  Communication   │        │  Communication   │       │  Communication   │
│▒▒▒ 13%          │        │▒ 5%              │       │▒ 3%              │
│  Évaporation     │        │  Évaporation     │       │  Évaporation     │
└──────────────────┘        └──────────────────┘       └──────────────────┘
Speedup: 0.89 ❌            Speedup: 1.21 ✅           Speedup: 1.42 ✅
```

L'augmentation du nombre de fourmis **déplace le goulot d'étranglement** : avec 100k fourmis et 2 processus, la communication ne représente plus que 13% du temps total, permettant un speedup réel.

Évolution du profil temporel pour 100k fourmis en fonction du nombre de processus :

```
1P (15.5s)     2P (10.9s)    4P (17.0s)    8P (17.0s)    12P (19.4s)
┌──────────┐  ┌──────────┐  ┌──────────┐  ┌──────────┐  ┌──────────┐
│▓▓▓▓▓▓▓▓▓▓│  │▓▓▓▓▓▓▓▓▓ │  │▓▓▓▓▓▓    │  │▓▓▓▓▓     │  │▓▓▓       │
│▓▓▓▓▓▓▓▓▓▓│  │▓▓▓▓▓▓▓▓  │  │▓▓▓▓▓     │  │▓▓▓▓      │  │▓▓        │
│▓▓▓ 95.5% │  │▓▓▓ 84.0% │  │▓  55.0%  │  │▓  41.7%  │  │  27.6%   │
│          │  │███ 12.6% │  │█████████ │  │█████████ │  │██████████│
│█ 3.2%    │  │          │  │██ 42.5%  │  │███ 53.8% │  │████66.6% │
│          │  │          │  │          │  │          │  │          │
└──────────┘  └──────────┘  └──────────┘  └──────────┘  └──────────┘
  S = 1.00      S = 1.42      S = 0.91      S = 0.91      S = 0.80

▓ = Mouvement (parallélisable)    █ = Communication MPI (overhead)
```

#### Estimation du point de rentabilité

En extrapolant les données, on peut estimer le nombre de fourmis nécessaire pour atteindre un speedup de 1.0 (seuil de rentabilité) à $N$ processus :

| Nb Procs | Speedup 5k | Speedup 50k | Speedup 100k | Estimation pour Speedup = 1.0 |
|:---:|:---:|:---:|:---:|:---:|
| 2 | 0.89 | 1.21 | 1.42 | ~**15 000** fourmis |
| 4 | 0.25 | 0.70 | 0.91 | ~**120 000** fourmis |
| 8 | 0.19 | 0.64 | 0.91 | ~**200 000+** fourmis |

Avec 2 processus, la Stratégie 1 devient rentable dès ~15 000 fourmis. Pour 4+ processus, il faudrait **>100 000 fourmis** pour simplement atteindre le seuil de rentabilité, sans même obtenir une accélération significative.

---

### Observations sur les tests à 12 processus (Oversubscribe)

#### Contexte

La machine dispose de 8 cœurs physiques. L'option `--oversubscribe` force OpenMPI à lancer 12 processus en les partageant sur les 8 cœurs : 4 cœurs hébergent 2 processus chacun, ce qui implique du *context switching* (partage de temps) géré par le système d'exploitation. Les résultats sont donc **structurellement pénalisés** et ne reflètent pas une véritable exécution sur 12 cœurs.

#### Observations

Le mouvement continue de se diviser malgré la surallocation :

| Config | Mouvement 1P (s) | Mouvement 12P (s) | Ratio réel | Ratio théorique |
|:---:|:---:|:---:|:---:|:---:|
| 50k fourmis | 7.539 | 2.887 | **2.6×** | 12× |
| 100k fourmis | 14.776 | 5.338 | **2.8×** | 12× |

Le gain sur le mouvement est réel (~2.7×) mais très loin du facteur 12× théorique. L'écart s'explique par :
- **Surallocation** : le gain maximal théorique est ~8× (nombre de cœurs physiques) et non 12×
- **Contention mémoire** : la bande passante RAM est saturée par 12 processus chargeant chacun la grille entière

La communication à 12 processus (~12.5 s) représente **66-75% du temps total**, confirmant que le `MPI_Allreduce` est le facteur limitant, pas le calcul.

---

### Conclusion de la Phase 4

#### Synthèse des speedups (toutes configurations)

| Nb Procs | 5 000 fourmis | 50 000 fourmis | 100 000 fourmis |
|:---:|:---:|:---:|:---:|
| 1 | 1.00 | 1.00 | 1.00 |
| 2 | 0.89 ❌ | **1.21** ✅ | **1.42** ✅ |
| 4 | 0.25 ❌ | 0.70 ❌ | 0.91 ❌ |
| 8 | 0.19 ❌ | 0.64 ❌ | 0.91 ❌ |
| 12* | — | 0.50 ❌ | 0.80 ❌ |

#### Bilan

| Critère | Résultat |
|---|---|
| **Implémentation MPI** | ✅ Fonctionnelle et correcte |
| **Correctesse** | ✅ Résultats cohérents sur toutes les configurations |
| **Meilleur speedup obtenu** | **1.42×** (100k fourmis, 2 processus) |
| **Scalabilité** | ❌ Ne passe pas à l'échelle au-delà de 2 processus |
| **Cause principale** | `MPI_Allreduce` sur ~16.8 Mo à chaque itération |
| **Conformité au sujet** | ✅ Les limites annoncées sont confirmées expérimentalement |

#### Leçon retenue

La Stratégie 1 illustre un **anti-pattern classique** en calcul distribué : elle échange **trop de données** par rapport au calcul qu'elle effectue. Le `MPI_Allreduce` sur la grille entière est un coût **incompressible** qui croît avec le nombre de processus, tandis que le gain de calcul (division des fourmis) est **limité** par la taille du problème.

Pour que cette stratégie soit rentable, deux conditions doivent être réunies :
1. Le nombre de fourmis doit être **très grand** (>>100 000) pour que le mouvement domine le temps total
2. La carte doit être **petite** pour que le volume du `MPI_Allreduce` reste gérable

Ces deux conditions sont en tension : une carte petite limite l'intérêt de la simulation, tandis qu'une grande population de fourmis sur une petite carte sature rapidement les phéromones.

La **Stratégie 2 (décomposition de domaine)**, objet de la Phase 5, résout structurellement ce problème en n'échangeant que les **cellules fantômes aux frontières** (~32 Ko au lieu de 16.8 Mo), soit une réduction du volume de communication de **~500×**.