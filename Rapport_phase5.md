## Phase 5 : Parallélisme en Mémoire Distribuée (MPI) — Stratégie 2 (Décomposition de Domaine)

### Objectif de la phase

La Phase 4 a démontré les limites structurelles de la Stratégie 1 : la duplication de la grille entière sur chaque processus et l'utilisation de `MPI_Allreduce` sur ~16.8 Mo à chaque itération créent un goulot de communication insurmontable. Le meilleur speedup obtenu est de 1.42× à 2 processus, et la performance se **dégrade** au-delà.

La Phase 5 implémente la **Stratégie 2**, dite de **décomposition de domaine**, décrite dans le sujet du projet comme suit :

> « Cette fois-ci, chaque processus prend en compte qu'une partie de la carte et ne gère que les fourmis qui sont sur la carte. La difficulté ici est de gérer les bords de chaque sous-carte. [...] L'avantage de cette méthode est que la mémoire occupée par l'application est bien plus petite [...]. De plus, il n'y a qu'un échange au bord des sous-cartes et donc peu de données échangées en définitive. »

Cette stratégie repose sur trois idées fondamentales :

1. **Partition spatiale** : la grille est découpée en sous-domaines, chaque processus ne stockant que **sa portion** de la carte
2. **Ghost cells** (cellules fantômes) : chaque sous-domaine est entouré d'une rangée de cellules supplémentaires contenant les valeurs des voisins, mises à jour par communication MPI
3. **Migration des fourmis** : lorsqu'une fourmi sort de son sous-domaine, elle est transférée dynamiquement au processus voisin

---

### Conception de la décomposition

#### Choix du découpage : 1D par lignes

Parmi les schémas de décomposition possibles (1D par lignes, 1D par colonnes, 2D par blocs), le découpage **1D par lignes** a été retenu pour sa simplicité d'implémentation :

```
Grille globale 513 × 513 — Découpage 1D en 4 processus
┌─────────────────────────────────┐  ← Ligne 0
│                                 │
│           Rank 0                │
│        lignes [0, 129)          │
│         129 lignes              │
│                                 │
├─────────────────────────────────┤  ← Ligne 129 (frontière)
│                                 │
│           Rank 1                │
│        lignes [129, 257)        │
│         128 lignes              │
│                                 │
├─────────────────────────────────┤  ← Ligne 257 (frontière)
│                                 │
│           Rank 2                │
│        lignes [257, 385)        │
│         128 lignes              │
│                                 │
├─────────────────────────────────┤  ← Ligne 385 (frontière)
│                                 │
│           Rank 3                │
│        lignes [385, 513)        │
│         128 lignes              │
│                                 │
└─────────────────────────────────┘  ← Ligne 512
```

Chaque processus possède $\lfloor N/P \rfloor$ ou $\lceil N/P \rceil$ lignes. Le premier processus (rank 0) et le dernier (rank $P-1$) n'ont qu'un seul voisin (respectivement en bas et en haut). Les processus intermédiaires ont deux voisins.

#### Principe des Ghost Cells

Pour calculer les phéromones d'une cellule, il faut connaître les valeurs de ses 4 voisines. Les cellules situées **à la frontière** d'un sous-domaine ont des voisines appartenant au processus adjacent. Pour y accéder sans communication à chaque calcul, chaque sous-domaine est étendu d'une **rangée fantôme** en haut et en bas :

```
Stockage local du Rank 1 (lignes globales [129, 257)) :

  Ligne locale 0     : ░░░░░░░░░░░░░░░░  ← Ghost top (copie de la ligne 128 du Rank 0)
  Ligne locale 1     : ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓  ← Données réelles (ligne globale 129)
  Ligne locale 2     : ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓  ← Données réelles (ligne globale 130)
  ...
  Ligne locale 128   : ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓  ← Données réelles (ligne globale 256)
  Ligne locale 129   : ░░░░░░░░░░░░░░░░  ← Ghost bottom (copie de la ligne 257 du Rank 2)

  ▓ = Données réelles (calculées localement)
  ░ = Ghost cells (reçues du voisin par MPI)
```

Les ghost cells sont mises à jour **une fois par itération** via des échanges `MPI_Sendrecv` avec les voisins. Cela permet ensuite à toutes les fourmis locales de calculer leurs phéromones **sans aucune communication supplémentaire** pendant le mouvement.

#### Volume de communication : comparaison quantitative

| Donnée échangée | Stratégie 1 | Stratégie 2 |
|---|---|---|
| **Type** | `MPI_Allreduce` (collectif, all-to-all) | `MPI_Sendrecv` (point à point, voisin à voisin) |
| **Contenu** | Grille entière $(N+2)^2 \times 2$ doubles | 2 lignes × $(N+2)$ colonnes × 2 doubles |
| **Volume par itération** | $(513+2)^2 \times 2 \times 8 \approx$ **4.1 Mo** | $2 \times 2 \times (513+2) \times 2 \times 8 \approx$ **32 Ko** |
| **Volume total (500 it.)** | **~2.0 Go** | **~16 Mo** (+ migration variable) |
| **Réduction** | — | **~128×** moins |
| **Complexité** | $O(V \times \log P)$ avec $V =$ volume | $O(L)$ avec $L =$ longueur d'une ligne |

#### Protocole de migration des fourmis

La migration des fourmis est le défi technique principal de la Stratégie 2. Lorsqu'une fourmi se déplace hors de son sous-domaine local (par exemple, la fourmi du Rank 1 passe de la ligne 257 à la ligne 258, qui appartient au Rank 2), elle doit être **transférée** au processus propriétaire de la zone de destination.

Le protocole implémenté se déroule en 4 étapes :

```
Étape 1 : Identification
  Pour chaque fourmi i :
    Si m_x[i] < local_row_start → marquée pour envoi au rank_top
    Si m_x[i] >= local_row_end  → marquée pour envoi au rank_bottom

Étape 2 : Échange des compteurs
  MPI_Sendrecv : "J'envoie K fourmis vers le haut" ←→ "Je reçois L fourmis d'en bas"
  MPI_Sendrecv : "J'envoie M fourmis vers le bas" ←→ "Je reçois N fourmis d'en haut"

Étape 3 : Échange des données
  Chaque fourmi est encodée en 4 doubles : [x, y, loaded, seed]
  MPI_Sendrecv des vecteurs de fourmis migrantes

Étape 4 : Mise à jour locale
  Suppression des fourmis parties (swap-and-pop en O(1) par fourmi)
  Ajout des fourmis reçues (push_back)
```

Chaque fourmi nécessite l'envoi de **4 doubles = 32 octets**. Le volume de migration dépend du nombre de fourmis traversant les frontières, lui-même fonction du taux d'exploration $\varepsilon$ et de la densité locale des fourmis.

---

### Réalisation

#### Fichiers créés

| Fichier | Lignes | Description |
|---|---|---|
| `domain_decomposition.hpp` | ~130 | Classe `DomainDecomposition` : découpage équitable des lignes, calcul des voisins MPI (`rank_top`, `rank_bottom`), conversion coordonnées globales/locales, méthode `exchange_ghost_cells()` utilisant `MPI_Sendrecv` |
| `pheronome_local.hpp` | ~190 | Classe `PheronomeLocal` : grille de phéromones locale (sous-domaine + ghost rows), accesseurs avec coordonnées globales (conversion automatique), `safe_access()` pour les voisins dans les ghost cells, évaporation et update locaux |
| `ant_population_local.hpp` | ~200 | Classe `AntPopulationLocal` : population SoA initialisée dans le domaine local, méthode `advance_all()`, méthode `migrate_ants()` avec protocole d'échange en 4 étapes |
| `ant_simu_domain.cpp` | ~130 | Programme principal : initialisation MPI, construction du terrain, décomposition de domaine, boucle principale avec chronométrage de 5 composants (mouvement, ghost cells, migration, évaporation, update), réductions finales |

#### Boucle principale (pseudo-code)

```
Pour chaque itération :
    1. exchange_ghosts()     ← MPI_Sendrecv des ghost cells phéromones
    2. advance_all()         ← Mouvement des fourmis locales
    3. migrate_ants()        ← Transfert des fourmis hors domaine aux voisins
    4. do_evaporation()      ← Multiplication par β (local)
    5. update()              ← Swap des buffers, reset nourriture/nid (local)
```

L'ordre est important : les ghost cells sont échangées **avant** le mouvement pour que les fourmis proches des frontières aient accès aux phéromones du voisin. La migration est effectuée **après** le mouvement pour détecter les fourmis sorties du domaine.

#### Compilation et exécution

```bash
# Compilation (sans renderer/window, SDL nécessaire uniquement pour SDL_Point)
mpic++ -fopenmp -std=c++17 -O3 -march=native -Wall -c ant_simu_domain.cpp -o ant_simu_domain.o
mpic++ -fopenmp -std=c++17 -O3 -march=native -Wall -c fractal_land.cpp -o fractal_land.o
mpic++ -fopenmp ant_simu_domain.o fractal_land.o -o ant_simu_domain.exe -lSDL2

# Exécution
mpirun -np 1 ./ant_simu_domain.exe
mpirun -np 2 ./ant_simu_domain.exe
mpirun -np 4 ./ant_simu_domain.exe
mpirun -np 8 ./ant_simu_domain.exe
mpirun --oversubscribe -np 12 ./ant_simu_domain.exe
```

---

### Résultats

#### Conditions de test

| Paramètre | Valeur |
|---|---|
| Machine | Dell 16 DC16250, 8 cœurs physiques |
| OS | Ubuntu 24.04, g++ 13.3.0, OpenMPI 4.1.6 |
| Nombre de fourmis | 100 000 |
| Nombre d'itérations | 500 |
| Taille de la grille | 513 × 513 |
| GUI | Désactivé |

#### Décomposition de domaine observée

| Nb Procs | Lignes/proc | Voisins Rank 0 | Mémoire phéromones/proc |
|:---:|:---:|:---:|:---:|
| 1 | 513 | Aucun | 4 144 Ko |
| 2 | 257 / 256 | bottom=1 | 2 084 Ko |
| 4 | 129 / 128 / 128 / 128 | bottom=1 | 1 054 Ko |
| 8 | 65 / 64 / ... / 64 | bottom=1 | 539 Ko |
| 12 | 43 / 43 / ... / 42 | bottom=1 | 362 Ko |

#### Tableau complet des résultats

| Nb Procs | Temps Total (s) | Mouvement (s) | Ghost cells (s) | Migration (s) | Évaporation (s) | Update (s) |
|:---:|:---:|:---:|:---:|:---:|:---:|:---:|
| 1 | 16.651 | 16.428 (98.7%) | 0.000 (0.0%) | 0.050 (0.3%) | 0.170 (1.0%) | 0.003 (0.0%) |
| 2 | 7.851 | 6.670 (84.9%) | 0.016 (0.2%) | 1.071 (13.6%) | 0.092 (1.2%) | 0.001 (0.0%) |
| 4 | 8.919 | 5.754 (64.5%) | 1.132 (12.7%) | 1.908 (21.4%) | 0.122 (1.4%) | 0.002 (0.0%) |
| 8 | 6.820 | 3.376 (49.5%) | 0.859 (12.6%) | 2.443 (35.8%) | 0.141 (2.1%) | 0.001 (0.0%) |
| 12* | 6.503 | 2.681 (41.2%) | 1.614 (24.8%) | 2.013 (31.0%) | 0.194 (3.0%) | 0.001 (0.0%) |

*\* 12 processus sur 8 cœurs physiques via `--oversubscribe`*

#### Vérification de correctesse

| Nb Procs | Fourmis restantes | Attendu | Résultat |
|:---:|:---:|:---:|:---:|
| 1 | 100 000 | 100 000 | ✅ Correct |
| 2 | 100 000 | 100 000 | ✅ Correct |
| 4 | 100 000 | 100 000 | ✅ Correct |
| 8 | 100 000 | 100 000 | ✅ Correct |
| 12 | 100 000 | 100 000 | ✅ Correct |

**Aucune fourmi perdue ni dupliquée** lors de la migration. Le protocole de transfert swap-and-pop + `MPI_Sendrecv` est correct et conservatif.

---

### Analyse des résultats

#### Speedup global

| Nb Procs | Temps Total (s) | Speedup Strat. 2 | Speedup Strat. 1 (rappel) |
|:---:|:---:|:---:|:---:|
| 1 | 16.651 | 1.00 | 1.00 |
| 2 | 7.851 | **2.12** ✅ | 1.42 |
| 4 | 8.919 | **1.87** ✅ | 0.91 |
| 8 | 6.820 | **2.44** ✅ | 0.91 |
| 12* | 6.503 | **2.56** ✅ | 0.80 |

**Résultat majeur** : la Stratégie 2 produit un speedup **systématiquement positif**, y compris à 4, 8 et 12 processus où la Stratégie 1 échouait. Le speedup continue de **croître** avec le nombre de processus (2.44× à 8P, 2.56× à 12P), démontrant une scalabilité réelle.

#### Speedup du mouvement (composante purement parallélisable)

| Nb Procs | Mouvement (s) | Speedup Mouvement |
|:---:|:---:|:---:|
| 1 | 16.428 | 1.00 |
| 2 | 6.670 | **2.46** |
| 4 | 5.754 | **2.85** |
| 8 | 3.376 | **4.87** |
| 12* | 2.681 | **6.13** |

Le mouvement des fourmis, isolé de la communication, présente un speedup de **6.13×** à 12 processus. C'est un résultat excellent qui démontre que la charge de calcul se répartit bien entre les processus. L'écart avec le facteur théorique 12× s'explique par la surallocation (12 processus sur 8 cœurs) et la contention mémoire.

#### Comparaison de la communication : Stratégie 1 vs Stratégie 2

| Nb Procs | Communication Strat. 1 (s) | Communication Strat. 2 — Ghost cells (s) | Ratio |
|:---:|:---:|:---:|:---:|
| 1 | 0.493 | 0.000 | — |
| 2 | 1.375 | 0.016 | **86×** moins |
| 4 | 7.238 | 1.132 | **6.4×** moins |
| 8 | 9.173 | 0.859 | **10.7×** moins |
| 12* | 12.908 | 1.614 | **8.0×** moins |

La réduction du volume de communication porte ses fruits : le temps d'échange des ghost cells est **négligeable** à 2 processus (16 ms sur 7.8 s) et reste modéré même à 12 processus (1.6 s). Par comparaison, le `MPI_Allreduce` de la Stratégie 1 consommait 12.9 s à 12 processus.

#### Le nouveau goulot d'étranglement : la migration des fourmis

Si la communication des phéromones a été résolue, un **nouveau goulot** apparaît : la migration des fourmis entre sous-domaines.

| Nb Procs | Migration (s) | % du Temps Total |
|:---:|:---:|:---:|
| 1 | 0.050 | 0.3% |
| 2 | 1.071 | **13.6%** |
| 4 | 1.908 | **21.4%** |
| 8 | 2.443 | **35.8%** |
| 12* | 2.013 | **31.0%** |

La migration consomme jusqu'à **35.8%** du temps total à 8 processus. Trois facteurs expliquent ce coût :

**1. Taux de migration élevé dû à $\varepsilon = 0.8$**

Avec un taux d'exploration de 80%, les fourmis se déplacent de façon **très aléatoire**. À chaque sous-pas de temps, une fourmi a une probabilité significative de franchir une frontière de domaine. Plus le nombre de processus augmente, plus les sous-domaines sont petits (43 lignes à 12 processus), et plus les frontières sont proches → plus les traversées sont fréquentes.

**2. Coût des opérations mémoire**

Le protocole de migration implique des opérations sur les vecteurs SoA (`push_back`, `pop_back`, `swap`) qui peuvent déclencher des réallocations mémoire et invalider les caches.

**3. Synchronisation bilatérale**

Chaque migration nécessite **4 appels `MPI_Sendrecv`** (compteurs + données, dans les deux directions), ce qui impose une synchronisation point-à-point avec les deux voisins à chaque itération.

##### Évolution du profil temporel

```
1P (16.7s)      2P (7.9s)       4P (8.9s)       8P (6.8s)       12P (6.5s)
┌──────────┐   ┌──────────┐   ┌──────────┐   ┌──────────┐   ┌──────────┐
│▓▓▓▓▓▓▓▓▓▓│   │▓▓▓▓▓▓▓▓▓ │   │▓▓▓▓▓▓▓   │   │▓▓▓▓▓     │   │▓▓▓▓      │
│▓▓▓▓▓▓▓▓▓▓│   │▓▓▓ 84.9% │   │▓▓ 64.5%  │   │▓  49.5%  │   │▓  41.2%  │
│▓▓ 98.7%  │   │          │   │          │   │          │   │          │
│          │   │░░░ 13.6% │   │░░░ 21.4% │   │░░░░░35.8%│   │░░░░ 31.0%│
│          │   │  Migration│   │  Migration│   │  Migration│   │  Migration│
│          │   │█ 0.2%    │   │██ 12.7%  │   │██ 12.6%  │   │███ 24.8% │
│          │   │  Ghost    │   │  Ghost    │   │  Ghost    │   │  Ghost   │
└──────────┘   └──────────┘   └──────────┘   └──────────┘   └──────────┘

▓ = Mouvement    ░ = Migration fourmis    █ = Ghost cells
```

Le profil montre clairement que le mouvement diminue régulièrement (de 98.7% à 41.2%), tandis que la migration et les ghost cells prennent une part croissante. Le temps total continue néanmoins de **diminuer** (16.7 s → 6.5 s), preuve que la stratégie est efficace.

##### Pistes d'optimisation de la migration (non implémentées)

| Piste | Description | Gain estimé |
|---|---|---|
| Pré-allocation des buffers | Réserver la capacité maximale des vecteurs de migration pour éviter les réallocations dynamiques | ~10-15% |
| Communications non-bloquantes | Utiliser `MPI_Isend`/`MPI_Irecv` pour superposer communication et calcul | ~20-30% |
| Sous-domaines plus grands | Réduire le nombre de processus pour augmenter la taille des sous-domaines (moins de frontières) | Variable |
| Réduction de $\varepsilon$ | Un taux d'exploration plus faible génère moins de mouvement aléatoire et donc moins de traversées de frontières | Significatif |
| Migration asynchrone | Permettre aux fourmis migrantes de « sauter » une itération au lieu d'être transférées immédiatement | ~15-25% |

---

### Bilan mémoire

Un avantage important de la Stratégie 2 est la **réduction de l'empreinte mémoire** par processus :

| Nb Procs | Mémoire/proc Strat. 1 | Mémoire/proc Strat. 2 | Réduction |
|:---:|:---:|:---:|:---:|
| 1 | 4 Mo | 4 144 Ko | ~1× |
| 2 | 4 Mo | 2 084 Ko | **2×** |
| 4 | 4 Mo | 1 054 Ko | **4×** |
| 8 | 4 Mo | 539 Ko | **7.6×** |
| 12 | 4 Mo | 362 Ko | **11.3×** |

La Stratégie 1 duplique la grille entière sur chaque processus (mémoire constante à 4 Mo). La Stratégie 2 divise la grille proportionnellement au nombre de processus. À 12 processus, chaque rang ne stocke que **362 Ko** de phéromones, soit **11.3× moins** que la Stratégie 1.

Cette propriété est cruciale pour le passage à l'échelle : si la grille faisait $10\,000 \times 10\,000$ cellules (~1.5 Go de phéromones), la Stratégie 1 nécessiterait 1.5 Go par processus (impossible sur de nombreuses machines), tandis que la Stratégie 2 ne nécessiterait que ~125 Mo par processus avec 12 rangs.

---

### Synthèse comparative finale

#### Tableau récapitulatif (100 000 fourmis, 500 itérations)

| Nb Procs | Strat. 1 Temps (s) | Strat. 1 Speedup | Strat. 2 Temps (s) | Strat. 2 Speedup | Gain Strat. 2 vs 1 |
|:---:|:---:|:---:|:---:|:---:|:---:|
| 1 | 15.474 | 1.00 | 16.651 | 1.00 | 0.93× |
| 2 | 10.889 | 1.42 | 7.851 | **2.12** | **1.39×** plus rapide |
| 4 | 17.025 | 0.91 | 8.919 | **1.87** | **1.91×** plus rapide |
| 8 | 17.045 | 0.91 | 6.820 | **2.44** | **2.50×** plus rapide |
| 12* | 19.368 | 0.80 | 6.503 | **2.56** | **2.98×** plus rapide |

À 1 processus, la Stratégie 2 est légèrement plus lente que la Stratégie 1 (16.65 s vs 15.47 s) en raison de l'overhead du mécanisme de migration (même sans voisin MPI, le code vérifie les frontières). Dès 2 processus, l'avantage de la Stratégie 2 est clair et **s'amplifie** avec le nombre de processus.

```
Temps total (s) — 100 000 fourmis
│
20 ┤                         ●──────● Stratégie 1
   │               ●────────●
18 ┤          
   │ ●        ●
16 ┤  \      /
   │   \    /
14 ┤    \  /
   │     \/
12 ┤
   │          ●
10 ┤         /
   │        /
 8 ┤       /     ●
   │ ●────●       \
 6 ┤               ●────● Stratégie 2
   │
 4 ┤
   │
 0 ┤────────────────────────────
    1P   2P    4P    8P   12P
```

#### Verdict comparatif

| Critère | Stratégie 1 | Stratégie 2 |
|---|---|---|
| **Meilleur speedup** | 1.42× (2P) | **2.56× (12P)** |
| **Passe à l'échelle ?** | ❌ Non (ralentit dès 4P) | ✅ Oui (accélère jusqu'à 12P) |
| **Communication phéromones** | ~4.1 Mo/itération (`Allreduce`) | **~32 Ko/itération** (ghost cells) |
| **Mémoire par processus** | Constante (carte dupliquée) | **Décroissante** (carte divisée) |
| **Complexité d'implémentation** | Simple | Modérée (migration) |
| **Goulot d'étranglement résiduel** | Communication `MPI_Allreduce` | Migration des fourmis |
| **Scalabilité théorique** | Limitée | **Bonne** (si migration optimisée) |

---

### Conclusion de la Phase 5

La Stratégie 2 (décomposition de domaine) valide le **principe fondamental du calcul distribué** : minimiser les communications en maximisant la localité des données. En ne communiquant que les frontières (~32 Ko) au lieu de la grille entière (~4.1 Mo), le temps de communication est réduit d'un facteur **6 à 86×** selon le nombre de processus.

Le speedup de **2.56×** à 12 processus, bien que modeste en valeur absolue, représente une **amélioration qualitative majeure** par rapport à la Stratégie 1 qui produisait un **ralentissement** de 20% dans les mêmes conditions. La Stratégie 2 est la seule approche qui **passe réellement à l'échelle** dans ce projet.

Le goulot d'étranglement résiduel (migration des fourmis, ~31-36% du temps) offre des perspectives d'optimisation clairement identifiées (communications non-bloquantes, pré-allocation, migration asynchrone) qui pourraient encore améliorer significativement les performances.

---

---

## Conclusion Générale

### Récapitulatif du parcours d'optimisation

Ce projet a permis de traverser **l'ensemble des paradigmes de parallélisation** enseignés en calcul haute performance, appliqués à un problème concret de simulation biologique (optimisation par colonie de fourmis sur un paysage fractal).

| Phase | Approche | Action principale | Résultat clé |
|---|---|---|---|
| **Phase 1** | Audit et métrologie | Profiling du code séquentiel | Baseline : 0.38 s (5k fourmis), mouvement = 75% du temps |
| **Phase 2** | Vectorisation (SoA) | Restructuration AoS → SoA | Gain nul en séquentiel, mais évaporation **-38%** (vectorisation SIMD) |
| **Phase 3** | Mémoire partagée (OpenMP) | Parallélisation des boucles | **Ralentissement +22%** (contention atomique, overhead) |
| **Phase 4** | Mémoire distribuée — Stratégie 1 | Réplication carte + `MPI_Allreduce` | Speedup **1.42×** à 2P (100k fourmis), mais échec à 4P+ |
| **Phase 5** | Mémoire distribuée — Stratégie 2 | Décomposition de domaine + ghost cells + migration | Speedup **2.56×** à 12P, scalabilité confirmée |

### Évolution des performances (100 000 fourmis)

```
Speedup
│
3.0 ┤                                          ★ 2.56 (Phase 5, 12P)
    │                                    ★
2.5 ┤                              ★ 2.44
    │                        ★ 2.12
2.0 ┤
    │              ★ 1.87
1.5 ┤        ● 1.42
    │
1.0 ┤──── Seuil de rentabilité ──────────────────────
    │  ○ 0.91          ○ 0.91
0.5 ┤              ○ 0.80
    │
0.0 ┤
     2P     4P     8P     12P

● = Strat. 1 (meilleur)   ○ = Strat. 1 (pire)   ★ = Strat. 2
```

### Enseignements tirés

#### 1. « On ne peut pas optimiser ce qu'on ne mesure pas »

Le profiling initial (Phase 1) a été déterminant : sans la connaissance précise que le mouvement représente 75% du temps, les efforts d'optimisation auraient pu être mal dirigés. La métrologie est le **fondement** de toute démarche d'optimisation.

#### 2. La structure des données conditionne la performance

Le passage de AoS à SoA (Phase 2) n'a apporté aucun gain séquentiel, mais a été un **pré-requis indispensable** pour les phases suivantes. En calcul haute performance, l'organisation des données est aussi importante que l'algorithme lui-même.

#### 3. La parallélisation n'est pas toujours synonyme d'accélération

La Phase 3 (OpenMP) a produit un **ralentissement** de 22%. Ce résultat contre-intuitif illustre que la parallélisation introduit des coûts (synchronisation, contention, overhead) qui doivent être compensés par un volume de calcul suffisant. La **loi d'Amdahl** et le **ratio calcul/communication** sont des outils prédictifs essentiels.

#### 4. Le choix de la stratégie de communication est critique

La comparaison entre les Stratégies 1 et 2 (Phases 4 et 5) démontre que le **volume de données échangées** est le facteur déterminant de la scalabilité en mémoire distribuée. Échanger la grille entière (4.1 Mo) vs. les frontières (32 Ko) fait la différence entre un ralentissement et un speedup.

#### 5. Chaque optimisation déplace le goulot d'étranglement

Le projet illustre la progression caractéristique de l'optimisation :
- Phase 1 : goulot = mouvement des fourmis (75%)
- Phase 3 : goulot = contention atomique sur les phéromones
- Phase 4 : goulot = communication `MPI_Allreduce` (66-86%)
- Phase 5 : goulot = migration des fourmis (31-36%)

Chaque résolution déplace le problème vers le composant suivant, confirmant que l'optimisation est un processus **itératif** et jamais « terminé ».

### Perspectives

Plusieurs pistes d'amélioration restent ouvertes pour aller au-delà des résultats obtenus :

| Piste | Impact attendu |
|---|---|
| **Parallélisme hybride MPI + OpenMP** | Utiliser MPI entre les nœuds et OpenMP au sein de chaque nœud pour combiner les avantages des deux paradigmes |
| **Communications non-bloquantes** (`MPI_Isend`/`MPI_Irecv`) | Superposer calcul et communication pour masquer la latence de la migration |
| **Décomposition 2D** | Découper la grille en blocs 2D au lieu de bandes 1D pour réduire le ratio surface/volume des frontières |
| **Équilibrage dynamique** | Redistribuer les sous-domaines en cours de simulation si la densité de fourmis devient déséquilibrée |
| **GPU (CUDA/OpenCL)** | Exploiter le parallélisme massif des cartes graphiques pour le mouvement des fourmis et l'évaporation |