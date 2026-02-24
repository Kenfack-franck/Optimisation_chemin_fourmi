# Optimisation Par Colonie De Fourmis D'un Chemin Sur Un Paysage Fractal

## Rapport de Projet — Calcul Haute Performance

**Auteur :** KENFACK Franck  
**Formation :** ENSTA — S2 Systèmes Distribués  
**Date :** Février 2026

---

## Introduction Générale

Ce projet porte sur l'optimisation et la parallélisation d'une simulation biologique inspirée du comportement des colonies de fourmis. L'algorithme simulé, appelé **Optimisation par Colonie de Fourmis (ACO — Ant Colony Optimization)**, résout le problème de *fourragement* : trouver le chemin le plus court entre une fourmilière et une source de nourriture sur un terrain accidenté.

Le code fourni initialement est **séquentiel** et repose sur une interface graphique SDL2. L'objectif du projet est de le transformer progressivement pour qu'il exploite les architectures modernes multi-cœurs et multi-nœuds, en suivant trois paradigmes :

1. **Vectorisation** (réorganisation des données pour favoriser le cache et le SIMD)
2. **Parallélisme en mémoire partagée** (OpenMP)
3. **Parallélisme en mémoire distribuée** (MPI)

Le présent rapport documente de manière détaillée chaque phase du travail réalisé, des mesures initiales jusqu'aux résultats de scalabilité.

---

## Rappel du Modèle Simulé

### Principe de la simulation

La simulation met en jeu les éléments suivants :

- **Un terrain** : une grille cartésienne 2D de taille $N \times N$ (par défaut $N = 513$). Chaque cellule contient une valeur entre 0 et 1 représentant le coût en temps pour la traverser. Ce terrain est généré de façon fractale (algorithme de type *plasma*) à partir d'une graine aléatoire, garantissant la reproductibilité.

- **$m$ fourmis artificielles** : des agents simples qui se déplacent sur la grille. Chaque fourmi peut être dans l'un de deux états :
  - **Non chargée** : elle cherche la nourriture
  - **Chargée** : elle a trouvé la nourriture et cherche à retourner au nid

- **Deux types de phéromones** stockés par cellule :
  - $V_1(s)$ : phéromone « nourriture » — guide les fourmis non chargées vers la source de nourriture. Vaut 1 à la position de la nourriture.
  - $V_2(s)$ : phéromone « nid » — guide les fourmis chargées vers la fourmilière. Vaut 1 à la position du nid.

### Dynamique à chaque pas de temps

À chaque itération, pour chaque fourmi :

1. **Mise à jour des phéromones** de la cellule courante $s$ en fonction des 4 voisines $N(s)$ :

$$V_i(s) \leftarrow \alpha \cdot \max_{s' \in N(s)} V_i(s') + (1-\alpha) \cdot \frac{1}{4}\sum_{s' \in N(s)} V_i(s')$$

   où $\alpha \in [0,1]$ est le paramètre de bruit (pondération entre max et moyenne).

2. **Déplacement** :
   - Avec probabilité $\varepsilon$ (taux d'exploration) : mouvement **aléatoire** vers une cellule voisine valide
   - Avec probabilité $1-\varepsilon$ : mouvement **guidé** vers la cellule voisine ayant le phéromone $V_1$ le plus fort (si non chargée) ou $V_2$ le plus fort (si chargée)

3. **Gestion de l'énergie** : la fourmi consomme un budget de mouvement de 1 par pas de temps. Chaque déplacement coûte la valeur de la cellule de destination. La fourmi continue de se déplacer tant que son budget n'est pas épuisé.

4. **Changement d'état** :
   - Arrivée sur la nourriture → état « chargée »
   - Arrivée au nid en état « chargée » → incrémentation du compteur de nourriture + état « non chargée »

Après le mouvement de toutes les fourmis, une **évaporation** globale est appliquée :

$$V_i(s) \leftarrow \beta \cdot V_i(s) \quad \forall s, \forall i$$

où $\beta \in [0,1]$ (proche de 1) est le coefficient d'évaporation.

### Paramètres utilisés dans ce projet

| Paramètre | Symbole | Valeur |
|---|---|---|
| Taille de la grille | $N$ | 513 (ou 1025 selon la config) |
| Nombre de fourmis | $m$ | 5 000 (baseline), 50 000 et 100 000 (tests de scalabilité) |
| Taux d'exploration | $\varepsilon$ | 0.8 |
| Paramètre de bruit | $\alpha$ | 0.7 |
| Coefficient d'évaporation | $\beta$ | 0.999 |
| Position du nid | | $(256, 256)$ |
| Position de la nourriture | | $(500, 500)$ |
| Graine aléatoire | | 2026 |
| Nombre d'itérations | | 500 |

---

## Phase 1 : Audit, Métrologie et Préparation

### Objectifs de la phase

Avant toute optimisation, il est indispensable de :

1. **Valider** que le code fourni compile et s'exécute correctement
2. **Isoler** le calcul pur de l'affichage graphique pour obtenir des mesures de temps fiables
3. **Profiler** le code pour identifier les sections les plus coûteuses (goulots d'étranglement) et orienter la stratégie d'optimisation

Cette phase constitue la **ligne de base (Baseline)** à laquelle toutes les optimisations futures seront comparées.

---

### Étape 1.1 : Compilation et Validation de l'Environnement Initial

#### Objectif

Valider l'environnement de développement, identifier les dépendances logicielles nécessaires et confirmer le bon fonctionnement fonctionnel du code séquentiel fourni.

#### Environnement technique

| Composant | Détail |
|---|---|
| Système d'exploitation | Linux Ubuntu 24.04 (Debian-based) |
| Machine | Dell 16 DC16250, 8 cœurs physiques |
| Compilateur | g++ 13.3.0 (GNU C++ Compiler) |
| Standard C++ | C++17 |
| Bibliothèques externes | SDL2 (Simple DirectMedia Layer) pour la visualisation |
| Outil de build | Makefile fourni avec le projet |

#### Problèmes rencontrés et résolutions

Lors de la première tentative de compilation via l'utilitaire `make`, deux problèmes de configuration ont été identifiés et corrigés :

**Problème 1 — Cible par défaut du Makefile**

La commande `make` sans argument affichait l'aide au lieu de lancer la compilation. Le Makefile ne définissait pas de cible par défaut.

*Correction :* Modification du `Makefile` pour définir la cible `all` comme cible par défaut, de sorte qu'un simple appel à `make` déclenche la compilation complète.

**Problème 2 — Dépendances manquantes**

Une erreur fatale de compilation a révélé l'absence des en-têtes de développement de la bibliothèque SDL2 :

```
fatal error: SDL2/SDL.h: Aucun fichier ou dossier de ce nom
```

*Correction :* Installation du paquet de développement SDL2 :

```bash
sudo apt install libsdl2-dev
```

#### Compilation réussie

Une fois les correctifs appliqués, la compilation s'est déroulée avec succès. La commande finale d'édition de liens :

```bash
g++ -fopenmp -std=c++17 -O3 -march=native -Wall ant.o fractal_land.o renderer.o window.o ant_simu.o -o ant_simu.exe -lpthread -lSDL2
```

Les drapeaux utilisés :
- `-O3 -march=native` : optimisation maximale par le compilateur, ciblant l'architecture matérielle locale
- `-fopenmp` : activation du support OpenMP (préparation pour les phases ultérieures)
- `-Wall` : activation de tous les avertissements pour un code propre
- `-lSDL2 -lpthread` : liaison avec les bibliothèques SDL2 et threads POSIX

> **Note :** À ce stade, le code est structurellement séquentiel. Le drapeau `-fopenmp` n'a pas d'effet visible mais prépare l'environnement pour la Phase 3.

#### Validation fonctionnelle

L'exécutable `ant_simu.exe` a été lancé et l'interface graphique s'est initialisée correctement. La simulation visuelle confirme le comportement attendu :

- Les agents (fourmis) se déplacent sur la grille fractale
- Les phéromones se déposent progressivement (visualisation par des couleurs : rouge pour $V_1$, vert pour $V_2$)
- Les zones claires du terrain ralentissent effectivement les fourmis
- La simulation est fluide et réactive

#### Conclusion de l'étape

Le code source de base est **sain** et l'environnement de développement est **opérationnel**. Nous disposons d'une version de référence fonctionnelle qui servira de point de comparaison pour les mesures de performances futures.

---

### Étape 1.2 : Mise en Place du Protocole de Benchmark

#### Objectif

L'affichage graphique SDL2 introduit un biais significatif dans les mesures de temps : le rafraîchissement de l'écran, la gestion des événements et la synchronisation verticale (VSync) consomment du temps qui n'a rien à voir avec le calcul de la simulation. Pour évaluer objectivement les gains de performance des futures optimisations, il est nécessaire de **dissocier le temps de calcul pur du temps d'affichage**.

L'objectif est de créer un mode d'exécution **non-interactif** et **déterministe** : la simulation s'exécute « à l'aveugle » pendant un nombre fixe d'itérations, sans aucune sortie graphique.

#### Modifications logicielles

Le fichier principal `ant_simu.cpp` a été modifié pour introduire deux modes de fonctionnement contrôlés par une variable booléenne `enable_gui` :

| Aspect | Mode Graphique (`enable_gui = true`) | Mode Calcul (`enable_gui = false`) |
|---|---|---|
| Fenêtre SDL | Initialisée | Non créée |
| Événements clavier/souris | Gérés | Ignorés |
| Rendu à l'écran | `renderer.display()` + `win.blit()` appelés | Aucun appel de rendu |
| Boucle principale | Boucle infinie (stoppée par l'utilisateur) | Boucle `for` bornée à `max_iterations = 500` |
| Sortie | Visuelle | Textuelle (temps mesurés) |

Les modifications concrètes dans le code :

1. Ajout de la variable `bool enable_gui = false;` en début de `main()`
2. Conditionnement de la création de la fenêtre et du renderer par `if (enable_gui)`
3. Conditionnement des appels de rendu (`renderer->display()`, `win->blit()`) par `if (enable_gui)`
4. Remplacement de la boucle `while` infinie par `for (std::size_t it = 0; it < max_iterations; ++it)`
5. Ajout d'un message de fin : `"Simulation terminée après 500 itérations."`

#### Validation

Le programme a été recompilé et exécuté en mode Calcul :

```
$ ./ant_simu.exe
Simulation terminée après 500 itérations.
```

Le programme s'exécute rapidement et se termine proprement sans ouvrir de fenêtre graphique.

#### Conclusion de l'étape

Nous disposons maintenant d'un **candidat stable** pour effectuer des mesures de temps précises. Le mode Benchmark garantit que seul le temps de calcul est mesuré, sans pollution par les opérations d'affichage.

---

### Étape 1.3 : Profiling et Analyse des Performances Initiales

#### Objectif

Identifier les **goulots d'étranglement** du code séquentiel en mesurant le temps passé dans chaque section critique de la simulation. Ces mesures constitueront notre **Baseline** — la référence à laquelle toutes les optimisations seront comparées.

#### Méthodologie

Le code a été instrumenté à l'aide de la bibliothèque `std::chrono` (horloge haute résolution `std::chrono::high_resolution_clock`). Des chronomètres ont été placés autour des trois sections critiques identifiées dans la boucle principale :

1. **Mouvement des fourmis** (`ant::advance`) : pour chaque fourmi, calcul de la direction, déplacement, mise à jour des phéromones, gestion de l'état chargé/non chargé
2. **Évaporation des phéromones** (`phen.do_evaporation()`) : parcours de toute la grille pour multiplier chaque valeur de phéromone par $\beta$
3. **Mise à jour** (`phen.update()`) : échange des buffers (double buffering) et réinitialisation des valeurs fixes (nourriture et nid)

Le temps est mesuré de façon **cumulée** sur les 500 itérations, puis le pourcentage de chaque composant est calculé par rapport au temps total.

#### Résultats

**Conditions de test :** 5 000 fourmis, 500 itérations, grille 513×513, GUI désactivé.

| Composant | Temps cumulé (s) | Pourcentage | Analyse |
|---|---|---|---|
| **Mouvement des fourmis** | 0.285 | **~74.8%** | Goulot d'étranglement principal |
| **Évaporation** | 0.095 | **~24.9%** | Charge significative |
| **Mise à jour (Update)** | 0.001 | ~0.3% | Négligeable |
| **Total** | **0.381** | 100% | |

**Temps moyen par itération :** 0.76 ms

#### Interprétation

La répartition du temps de calcul révèle une situation classique en optimisation :

**1. Le mouvement des fourmis domine largement (75%)**

C'est la section la plus coûteuse car elle implique, pour chaque fourmi et à chaque sous-pas de temps :
- La lecture des phéromones des 4 cellules voisines (accès mémoire potentiellement non contigus)
- Un tirage aléatoire pour la décision exploration vs exploitation
- Un ou plusieurs tirages aléatoires pour la direction (en cas d'exploration)
- La mise à jour de la position et le dépôt de phéromones
- La vérification des conditions de chargement/déchargement

De plus, la structure de données actuelle utilise un `std::vector<ant>` (Array of Structures — AoS). Cela signifie que les données de chaque fourmi (position, état, graine) sont stockées ensemble en mémoire, mais les données du **même attribut** de fourmis différentes sont **dispersées**. Lors du parcours de la boucle `for` sur les fourmis, le processeur charge en cache des données dont seule une partie est utilisée à chaque instruction, ce qui provoque des **défauts de cache** (cache misses).

**2. L'évaporation représente un quart du temps (25%)**

L'évaporation parcourt systématiquement **toute la grille** ($513 \times 513 = 263\,169$ cellules × 2 types de phéromones). C'est une opération très régulière (chaque cellule est traitée indépendamment), ce qui en fait une candidate idéale pour :
- La **vectorisation SIMD** (traitement de plusieurs cellules en parallèle par instruction)
- La **parallélisation OpenMP** (distribution des lignes entre les threads)

**3. La mise à jour est négligeable (0.3%)**

L'opération `update()` se résume à un `swap` de pointeurs entre deux buffers et la réinitialisation de deux cellules (nourriture et nid). Son coût est constant et indépendant de la taille du problème. Il n'est pas pertinent de l'optimiser.

#### Stratégie d'optimisation découlant du profiling

Conformément à la **loi d'Amdahl**, l'accélération maximale théorique est limitée par la fraction séquentielle du code. Avec 75% du temps dans le mouvement des fourmis :

- Si on parallélise **uniquement** le mouvement avec $P$ processeurs : $\text{Speedup}_{\max} = \frac{1}{0.25 + 0.75/P}$
- Avec $P = 8$ : Speedup$_{\max} \approx 2.9\times$
- Avec $P \to \infty$ : Speedup$_{\max} = 4.0\times$

Pour dépasser cette limite, il faudra aussi paralléliser l'évaporation.

**Plan d'action retenu :**

1. **Phase 2 (Vectorisation)** : Restructurer les données des fourmis de AoS vers SoA (Structure of Arrays) pour améliorer la localité de cache et permettre la vectorisation SIMD du compilateur
2. **Phase 3 (OpenMP)** : Paralléliser les boucles de mouvement et d'évaporation avec des directives OpenMP
3. **Phases 4-5 (MPI)** : Distribuer le calcul sur plusieurs processus, soit en dupliquant la carte (Stratégie 1), soit en la découpant (Stratégie 2)

#### Conclusion de l'étape

Le profiling a permis d'identifier clairement le goulot d'étranglement : **le mouvement des fourmis consomme 75% du temps de calcul**. La structure de données actuelle (AoS) est sous-optimale pour le cache. La stratégie d'optimisation est définie et priorisée selon la loi d'Amdahl.

---

### Conclusion de la Phase 1

| Étape | Objectif | Résultat |
|---|---|---|
| 1.1 | Compilation et validation | [OK] Code fonctionnel, environnement opérationnel |
| 1.2 | Mode benchmark sans GUI | [OK] Mesures de temps fiables, 500 itérations déterministes |
| 1.3 | Profiling initial | [OK] Baseline établie : 0.38 s total, 75% mouvement, 25% évaporation |

**Baseline de référence :** 5 000 fourmis, 500 itérations → **0.381 secondes** (0.76 ms/itération)

Nous disposons maintenant de tous les éléments pour passer à l'optimisation du code : une version de référence validée, un protocole de mesure fiable, et une connaissance précise de la répartition du temps de calcul.



## Phase 2 : Vectorisation — Restructuration des Données

### Objectif de la phase

Le profiling de la Phase 1 a révélé que le mouvement des fourmis consomme 75% du temps de calcul. L'analyse de la structure de données a identifié un problème fondamental : le code original utilise un `std::vector<ant>`, c'est-à-dire un tableau d'objets (Array of Structures — AoS). Cette organisation mémoire est sous-optimale pour deux raisons :

1. **Mauvaise localité de cache** : lorsqu'on parcourt les fourmis dans une boucle, le processeur charge en cache l'intégralité de chaque objet `ant` (position, état, graine, padding), alors que chaque instruction n'utilise qu'un seul attribut à la fois. Les lignes de cache sont sous-exploitées.

2. **Obstacle à la vectorisation SIMD** : les instructions SIMD (SSE, AVX) traitent plusieurs données identiques en parallèle (par exemple, 4 ou 8 `double` en une seule instruction). Pour en bénéficier, les données du même type doivent être contiguës en mémoire, ce qui n'est pas le cas avec AoS.

L'objectif de cette phase est de transformer la structure de données en **Structure of Arrays (SoA)** : au lieu d'un tableau de fourmis, on crée des tableaux séparés pour chaque attribut. Cette réorganisation est un **pré-requis technique indispensable** pour la parallélisation efficace (OpenMP et MPI) dans les phases suivantes.

### Contexte technique : AoS vs SoA

#### Organisation AoS (Array of Structures) — Code original

Dans le code initial, chaque fourmi est un objet autonome de la classe `ant` :

```cpp
class ant {
    static double m_eps;    // Coefficient d'exploration (partagé)
    std::size_t m_seed;     // Graine aléatoire propre
    state m_state;          // État : chargée ou non chargée
    position_t m_position;  // Position (x, y)
};

// Utilisation : un vecteur de N objets ant
std::vector<ant> ants(nb_ants);
```

En mémoire, les données sont organisées ainsi (illustration pour 4 fourmis) :

```
Adresse mémoire →
┌──────────┬──────────┬──────────┬──────────┐
│ Fourmi 0 │ Fourmi 1 │ Fourmi 2 │ Fourmi 3 │
│ seed₀    │ seed₁    │ seed₂    │ seed₃    │
│ state₀   │ state₁   │ state₂   │ state₃   │
│ pos_x₀   │ pos_x₁   │ pos_x₂   │ pos_x₃   │
│ pos_y₀   │ pos_y₁   │ pos_y₂   │ pos_y₃   │
│ padding  │ padding  │ padding  │ padding  │
└──────────┴──────────┴──────────┴──────────┘
```

**Problème** : pour traiter toutes les positions X (par exemple, pour le déplacement), le processeur doit charger les données de chaque fourmi entière, en sautant par-dessus les champs `seed`, `state`, `pos_y` et le padding d'alignement. Les accès mémoire sont **non contigus** et les lignes de cache sont gaspillées.

#### Organisation SoA (Structure of Arrays) — Nouvelle version

Avec SoA, on sépare chaque attribut dans son propre tableau :

```cpp
class AntPopulation {
    std::vector<int> m_x;           // Toutes les positions X
    std::vector<int> m_y;           // Toutes les positions Y
    std::vector<int> m_loaded;      // Tous les états (0 ou 1)
    std::vector<std::size_t> m_seeds; // Toutes les graines
    double m_eps;                    // Coefficient d'exploration
};
```

En mémoire :

```
Adresse mémoire →
m_x:      ┌──────┬──────┬──────┬──────┐
           │ x₀   │ x₁   │ x₂   │ x₃   │  ← Contigu !
           └──────┴──────┴──────┴──────┘
m_y:      ┌──────┬──────┬──────┬──────┐
           │ y₀   │ y₁   │ y₂   │ y₃   │  ← Contigu !
           └──────┴──────┴──────┴──────┘
m_loaded: ┌──────┬──────┬──────┬──────┐
           │ l₀   │ l₁   │ l₂   │ l₃   │  ← Contigu !
           └──────┴──────┴──────┴──────┘
m_seeds:  ┌──────┬──────┬──────┬──────┐
           │ s₀   │ s₁   │ s₂   │ s₃   │  ← Contigu !
           └──────┴──────┴──────┴──────┘
```

**Avantage** : lorsqu'une instruction traite toutes les positions X, les données sont **parfaitement contiguës** en mémoire. Le processeur peut charger une ligne de cache entière contenant 16 positions X d'un coup (64 octets / 4 octets par `int`). De plus, le compilateur peut vectoriser la boucle en utilisant des instructions SIMD qui traitent 4 ou 8 positions simultanément.

### Réalisation

#### Fichiers modifiés et créés

| Fichier | Action | Description |
|---|---|---|
| `ant_population.hpp` | **Créé** | Nouvelle classe `AntPopulation` en SoA, contenant les vecteurs `m_x`, `m_y`, `m_loaded`, `m_seeds` et la méthode `advance_all()` |
| `ant_simu.cpp` | **Modifié** | Remplacement du `std::vector<ant>` par un objet `AntPopulation`, adaptation de la boucle principale |
| `renderer.hpp` / `renderer.cpp` | **Modifiés** | Le renderer prend maintenant une référence vers `AntPopulation` au lieu de `std::vector<ant>`. L'affichage boucle sur `m_ref_ants.m_x[i]` et `m_ref_ants.m_y[i]` |
| `ant.cpp` / `ant.hpp` | **Conservés** (non utilisés) | L'ancienne classe `ant` n'est plus instanciée. Les fichiers sont conservés à titre de référence |

#### Détail de la classe `AntPopulation`

La classe `AntPopulation` centralise toute la logique de la simulation des fourmis :

**Constructeur** — Initialisation de la population :
```cpp
AntPopulation(size_t nb_ants, const fractal_land& land, double eps, size_t global_seed)
```
- Alloue les 4 vecteurs à la taille `nb_ants`
- Initialise chaque fourmi avec une position aléatoire sur la grille (via un générateur déterministe dépendant de la graine)
- Attribue une graine unique à chaque fourmi pour garantir des parcours différents
- Toutes les fourmis commencent dans l'état « non chargée » (`m_loaded[i] = 0`)

**Méthode `advance_all()`** — Mouvement de toutes les fourmis :
```cpp
void advance_all(pheronome& phen, const fractal_land& land,
                 const position_t& pos_food, const position_t& pos_nest,
                 std::size_t& cpteur_food)
```
- Boucle sur toutes les fourmis par leur indice `i`
- Pour chaque fourmi, appelle `advance_one_ant(i, ...)` qui reproduit **exactement** la même logique que l'ancienne méthode `ant::advance()` :
  1. Lecture des phéromones des 4 voisins
  2. Décision exploration/exploitation basée sur le tirage aléatoire et le coefficient $\varepsilon$
  3. Déplacement et consommation du budget de mouvement
  4. Dépôt de phéromones sur la nouvelle position
  5. Gestion du chargement (nourriture) et déchargement (nid)

**Différence clé** : au lieu d'accéder aux attributs via `this->m_position`, `this->m_state`, etc., la méthode accède aux tableaux contigus via `m_x[i]`, `m_y[i]`, `m_loaded[i]`, `m_seeds[i]`.

#### Adaptation du programme principal

Dans `ant_simu.cpp`, le changement principal est :

**Avant (AoS) :**
```cpp
std::vector<ant> ants;
for (size_t i = 0; i < nb_ants; ++i)
    ants.emplace_back(/* position, seed */);
// ...
for (auto& a : ants)
    a.advance(phen, land, pos_food, pos_nest, food_quantity);
```

**Après (SoA) :**
```cpp
AntPopulation population(nb_ants, land, eps, seed);
// ...
population.advance_all(phen, land, pos_food, pos_nest, food_quantity);
```

Le code est plus concis et l'objet `AntPopulation` encapsule toute la complexité.

#### Adaptation du renderer

Le renderer a été adapté pour lire les positions directement depuis les vecteurs SoA :

**Avant :**
```cpp
for (auto& a : ants)
    win.pset(a.get_position().x, a.get_position().y);
```

**Après :**
```cpp
for (size_t i = 0; i < m_ref_ants.m_x.size(); ++i)
    win.pset(m_ref_ants.m_x[i], m_ref_ants.m_y[i]);
```

### Résultats

#### Conditions de test

- **Configuration identique à la Phase 1** : 5 000 fourmis, 500 itérations, grille 513×513, GUI désactivé
- **Compilation** : mêmes drapeaux (`-O3 -march=native -fopenmp`)
- **Exécution** : 1 seul thread (séquentiel pur)

#### Tableau comparatif

| Composant | Version Originale (AoS) | Version Vectorisée (SoA) | Variation |
|---|---|---|---|
| **Temps Total** | 0.380 s | 0.382 s | +0.5% (bruit de mesure) |
| **Mouvement des fourmis** | 0.284 s | 0.322 s | **+13.4%** |
| **Évaporation** | 0.094 s | 0.058 s | **−38.3%** |
| **Mise à jour (Update)** | 0.001 s | 0.002 s | Négligeable |

#### Répartition du temps de calcul

```
Version AoS (originale)           Version SoA (vectorisée)
┌──────────────────────────┐     ┌──────────────────────────┐
│▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓ 74.8%│     │▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓84.3%│
│   Mouvement              │     │   Mouvement               │
│                          │     │                           │
│▒▒▒▒▒▒▒ 24.9%            │     │▒▒▒▒ 15.2%                │
│   Évaporation            │     │   Évaporation             │
│                          │     │                           │
│░ 0.3%                    │     │░ 0.5%                     │
│   Update                 │     │   Update                  │
└──────────────────────────┘     └──────────────────────────┘
      Total: 0.380 s                   Total: 0.382 s
```

### Interprétation détaillée

#### 1. Mouvement des fourmis : +13.4% (régression apparente)

La légère augmentation du temps de mouvement en SoA peut paraître surprenante, mais s'explique par plusieurs facteurs :

**Accès à plusieurs tableaux distincts** : en AoS, toutes les données d'une fourmi sont dans le même objet et potentiellement dans la même ligne de cache. En SoA, pour traiter une seule fourmi, il faut accéder à 4 tableaux différents (`m_x[i]`, `m_y[i]`, `m_loaded[i]`, `m_seeds[i]`), ce qui peut générer davantage de défauts de cache **en mode séquentiel**.

**Nature de l'algorithme de mouvement** : la boucle interne de `advance_one_ant()` n'est pas vectorisable par le compilateur car elle contient :
- Des **branchements conditionnels** dépendants des données (choix de direction, collision avec les murs)
- Des **boucles `do-while`** de longueur variable (tant que la fourmi tombe sur un mur)
- Des **accès mémoire indirects** aux phéromones via les coordonnées calculées dynamiquement
- Des **dépendances de données** intra-itération (la position de la fourmi à l'itération $k+1$ dépend de la position à l'itération $k$)

Le compilateur ne peut donc pas exploiter le caractère contigu des tableaux SoA pour cette partie du code. Le bénéfice de SoA se manifestera lors de la **parallélisation** (Phase 3), quand plusieurs fourmis seront traitées simultanément par des threads différents.

**En résumé** : la pénalité de +13% est le prix à payer pour une structure de données qui sera beaucoup plus performante en parallèle.

#### 2. Évaporation : −38.3% (gain spectaculaire)

Ce gain majeur s'explique par la combinaison de deux facteurs :

**Meilleure vectorisation par le compilateur** : la boucle d'évaporation parcourt linéairement un `std::vector<double>` plat :
```cpp
for (std::size_t i = 1; i <= m_dim; ++i)
    for (std::size_t j = 1; j <= m_dim; ++j) {
        m_buffer[idx]     *= m_beta;
        m_buffer[idx + 1] *= m_beta;
    }
```

Ce pattern (parcours séquentiel d'un tableau de `double` avec une multiplication scalaire) est un cas idéal pour la **vectorisation automatique** par le compilateur. Avec le drapeau `-O3 -march=native`, g++ peut utiliser des instructions AVX2 pour traiter 4 `double` par instruction au lieu d'un seul, soit un facteur théorique de 4×.

**Excellente localité de cache** : le stockage plat `[V1₀, V2₀, V1₁, V2₁, ...]` est parfaitement séquentiel en mémoire. Le prefetcher matériel du processeur peut anticiper les accès futurs et pré-charger les données en cache avant qu'elles ne soient nécessaires. Le taux de cache miss est quasi nul.

#### 3. Bilan global : gain nul (+0.5%) mais investissement stratégique

Le temps total reste identique (0.38 s) car le gain sur l'évaporation (−36 ms) est compensé par la pénalité sur le mouvement (+38 ms). Cela confirme que la restructuration SoA n'apporte pas de bénéfice en mode **séquentiel pur** pour un algorithme à branchements complexes.

Cependant, cette transformation est un **investissement stratégique** dont les dividendes se récolteront dans les phases suivantes :

| Bénéfice | Phase concernée |
|---|---|
| Chaque fourmi est identifiée par un simple indice `i` → facile à distribuer entre threads/processus | Phases 3, 4, 5 |
| Les dépendances entre fourmis sont supprimées (chaque fourmi a sa propre graine `m_seeds[i]`) | Phases 3, 4, 5 |
| Les tableaux contigus permettent un partitionnement trivial (fourmis 0 à N/2 → Thread 0, N/2 à N → Thread 1) | Phase 3 |
| La classe `AntPopulation` peut être instanciée localement par chaque processus MPI avec un sous-ensemble de fourmis | Phases 4, 5 |

### Conclusion de la Phase 2

| Critère | Résultat |
|---|---|
| Restructuration AoS → SoA | [OK] Complète |
| Impact sur le temps total (séquentiel) | Neutre (+0.5%) |
| Impact sur le mouvement | +13.4% (attendu en séquentiel) |
| Impact sur l'évaporation | **−38.3%** (vectorisation SIMD automatique) |
| Préparation pour le parallélisme | [OK] Structure prête pour OpenMP et MPI |
| Indépendance des fourmis | [OK] Chaque fourmi a sa propre graine, pas de dépendance inter-fourmis |

**La transformation SoA ne produit pas d'accélération séquentielle notable, mais elle supprime les obstacles structurels à la parallélisation.** C'est une étape fondamentale de l'approche *Data-Oriented Design* : on adapte la structure des données à l'algorithme de traitement (et non l'inverse) pour maximiser la performance sur les architectures modernes.



## Phase 3 : Parallélisme en Mémoire Partagée (OpenMP)

### Objectif de la phase

La Phase 2 a préparé le terrain en restructurant les données au format SoA. L'objectif de cette phase est d'exploiter les **multiples cœurs** de la machine de test pour accélérer les boucles de calcul, en utilisant les directives de compilation **OpenMP** (`#pragma omp`).

OpenMP est un standard de parallélisation en **mémoire partagée** : tous les threads accèdent au même espace d'adressage (même RAM). C'est le modèle le plus simple à mettre en œuvre car il ne nécessite que l'ajout de quelques directives au-dessus des boucles existantes. Cependant, il impose de gérer les **accès concurrents** lorsque plusieurs threads écrivent au même emplacement mémoire.

### Identification des sections parallélisables

D'après le profiling de la Phase 1, trois sections de code sont candidates à la parallélisation :

| Section | % du temps | Parallélisable ? | Difficulté |
|---|---|---|---|
| Mouvement des fourmis (`advance_all`) | 75-84% | [OK] Oui (boucle sur les fourmis) | **Élevée** — écriture concurrente sur la grille de phéromones |
| Évaporation (`do_evaporation`) | 15-25% | [OK] Oui (boucle sur la grille) | **Faible** — chaque cellule est indépendante |
| Mise à jour (`update`) | <1% | [FAIL] Non pertinent | Négligeable |

### Analyse des dépendances de données

Avant de paralléliser, il est crucial d'identifier les **données partagées** et les **conflits d'accès potentiels** :

#### Mouvement des fourmis — Données partagées

```
Thread 0 (fourmis 0..N/4)      Thread 1 (fourmis N/4..N/2)
         │                              │
         ▼                              ▼
    advance_one_ant(i)             advance_one_ant(j)
         │                              │
         ├── Lecture m_x[i], m_y[i]     ├── Lecture m_x[j], m_y[j]
         │   (privé [OK])                  │   (privé [OK])
         │                              │
         ├── Lecture phen(x±1, y±1)     ├── Lecture phen(x±1, y±1)
         │   (partagé, lecture seule [OK]) │   (partagé, lecture seule [OK])
         │                              │
         ├── Écriture phen.mark(pos)    ├── Écriture phen.mark(pos)
         │   ⚠️ CONFLIT si pos == pos    │   ⚠️ CONFLIT si même cellule
         │                              │
         └── Incrémentation cpteur_food └── Incrémentation cpteur_food
             ⚠️ CONFLIT (compteur partagé)   ⚠️ CONFLIT
```

**Deux conflits identifiés :**

1. **Écriture concurrente sur les phéromones** : si deux fourmis gérées par des threads différents se trouvent sur la même cellule (ou des cellules voisines au sein de la même ligne de cache), les écritures dans le buffer de phéromones se chevauchent → **Race condition**

2. **Incrémentation du compteur de nourriture** : variable scalaire incrémentée par n'importe quel thread lorsqu'une fourmi chargée atteint le nid → **Race condition**

#### Évaporation — Aucune dépendance

La boucle d'évaporation est **trivialement parallélisable** : chaque cellule `(i, j)` est traitée indépendamment et l'opération est une simple multiplication `m_buffer[idx] *= m_beta`. Il n'y a aucune dépendance entre les itérations.

### Réalisation

#### 3.1. Protection des données partagées dans `pheronome.hpp`

Pour résoudre le conflit d'écriture sur les phéromones, des directives `#pragma omp atomic write` ont été ajoutées dans la méthode `mark_pheronome()` :

```cpp
void mark_pheronome(const position_t& pos) {
    // ... calcul de val1 et val2 ...

    size_t idx = ((i + 1) * m_stride + (j + 1)) * 2;

    #pragma omp atomic write
    m_buffer[idx] = val1;

    #pragma omp atomic write
    m_buffer[idx + 1] = val2;
}
```

**Explication** : la directive `atomic write` garantit que l'écriture d'un `double` à une adresse donnée est indivisible. Si deux threads tentent d'écrire simultanément à la même adresse, les écritures sont sérialisées : l'une s'exécute avant l'autre, sans corruption de données.

**Choix de `atomic write` vs alternatives :**

| Alternative | Avantage | Inconvénient |
|---|---|---|
| `#pragma omp critical` | Protection globale | **Trop grossier** : sérialise toutes les écritures, même sur des cellules différentes |
| `#pragma omp atomic write` | **Fine granularité** : ne sérialise que les écritures à la même adresse | Coût non nul par écriture |
| Grille privée par thread + réduction | Aucune contention | **Mémoire prohibitive** : $N_{\text{threads}} \times 513^2 \times 2 \times 8$ octets |
| Mutex par cellule | Fine granularité | Surcharge mémoire et d'administration |

Le choix de `atomic write` représente le meilleur compromis pour notre cas d'usage.

#### 3.2. Parallélisation du mouvement des fourmis dans `ant_population.hpp`

La boucle principale de `advance_all()` a été annot��e avec OpenMP :

```cpp
void advance_all(pheronome& phen, const fractal_land& land,
                 const position_t& pos_food, const position_t& pos_nest,
                 std::size_t& cpteur_food)
{
    const size_t nb_ants = m_x.size();

    #pragma omp parallel for schedule(dynamic) reduction(+:cpteur_food)
    for (size_t i = 0; i < nb_ants; ++i) {
        advance_one_ant(i, phen, land, pos_food, pos_nest, cpteur_food);
    }
}
```

**Détail des clauses OpenMP utilisées :**

| Clause | Rôle |
|---|---|
| `parallel for` | Crée une équipe de threads et distribue les itérations de la boucle entre eux |
| `schedule(dynamic)` | Chaque thread prend les itérations une par une (ou par petit bloc) au fur et à mesure qu'il finit les précédentes. Choisi car le temps par fourmi est **variable** (dépend du terrain traversé et du nombre de sous-pas de temps) |
| `reduction(+:cpteur_food)` | Chaque thread possède une copie privée du compteur, initialisée à 0. À la fin de la région parallèle, toutes les copies sont additionnées dans la variable originale. Cela élimine le conflit d'accès au compteur |

#### 3.3. Parallélisation de l'évaporation dans `pheronome.hpp`

```cpp
void do_evaporation() {
    #pragma omp parallel for
    for (std::size_t i = 1; i <= m_dim; ++i) {
        for (std::size_t j = 1; j <= m_dim; ++j) {
            size_t idx = ((i * m_stride) + j) * 2;
            m_buffer[idx]     *= m_beta;
            m_buffer[idx + 1] *= m_beta;
        }
    }
}
```

**Choix de `schedule(static)`** (implicite, par défaut) : la charge de travail est parfaitement uniforme — chaque ligne de la grille demande exactement le même nombre d'opérations. Le scheduler statique divise les lignes en blocs égaux attribués aux threads une seule fois, sans surcharge d'ordonnancement.

#### Récapitulatif des modifications

| Fichier | Modification | Directive OpenMP |
|---|---|---|
| `pheronome.hpp` — `mark_pheronome()` | Protection des écritures concurrentes | `#pragma omp atomic write` |
| `pheronome.hpp` — `do_evaporation()` | Parallélisation du parcours de la grille | `#pragma omp parallel for` |
| `ant_population.hpp` — `advance_all()` | Parallélisation de la boucle sur les fourmis | `#pragma omp parallel for schedule(dynamic) reduction(+:cpteur_food)` |

### Résultats

#### Conditions de test

- **Machine** : Dell 16 DC16250, 8 cœurs physiques (rapportés comme 12 cœurs logiques par le système via Hyper-Threading ou configuration spécifique)
- **Configuration** : 5 000 fourmis, 500 itérations, grille 513×513, GUI désactivé
- **Threads OpenMP** : 12 (valeur par défaut = nombre de cœurs logiques détectés par le système)
- **Compilation** : `g++ -fopenmp -std=c++17 -O3 -march=native`

#### Tableau comparatif (Strong Scaling)

| Composant | Séquentiel (SoA) | OpenMP (12 threads) | Variation |
|---|---|---|---|
| **Temps Total** | 0.383 s | 0.468 s | **+22.2%** (ralentissement) |
| **Mouvement** | 0.323 s | 0.403 s | **+24.8%** |
| **Évaporation** | 0.058 s | 0.062 s | +6.9% |
| **Update** | 0.002 s | 0.003 s | Négligeable |

#### Speedup mesuré

$$\text{Speedup} = \frac{T_{\text{séquentiel}}}{T_{\text{OpenMP}}} = \frac{0.383}{0.468} = 0.82$$

**Le speedup est inférieur à 1** : l'utilisation d'OpenMP a **ralenti** le programme de 22% par rapport à la version séquentielle.

### Analyse critique : pourquoi c'est plus lent ?

Ce résultat, bien que contre-intuitif, est **caractéristique** de certaines situations en parallélisation. Trois phénomènes se conjuguent pour expliquer le ralentissement :

#### 1. Contention atomique et faux partage (*False Sharing*)

Le problème principal réside dans la méthode `mark_pheronome()` appelée à chaque sous-pas de temps de chaque fourmi.

**Contention atomique** : avec un taux d'exploration $\varepsilon = 0.8$, les fourmis se concentrent rapidement autour de deux zones d'intérêt : le **nid** (256, 256) et la **nourriture** (500, 500). Lorsque plusieurs fourmis se trouvent sur des cellules proches, leurs threads respectifs tentent d'exécuter des `atomic write` sur des adresses mémoire voisines. Le processeur doit alors **sérialiser** ces écritures, annulant le bénéfice du parallélisme.

**Faux partage** (*false sharing*) : même si deux fourmis écrivent sur des cellules **différentes**, si ces cellules se trouvent dans la même **ligne de cache** (64 octets, soit 8 `double`), le protocole de cohérence de cache (MESI) invalide la ligne de cache de tous les threads qui y accèdent. Cela provoque un trafic inter-cœurs massif :

```
Thread 0 écrit phen[100]    Thread 1 écrit phen[101]
         │                           │
         ▼                           ▼
   ┌─────────────────────────────────────┐
   │ Ligne de cache (64 octets)          │
   │ phen[96] phen[97] ... phen[103]     │
   └─────────────────────────────────────┘
         │                           │
   INVALIDATION ←──────────────→ INVALIDATION
   (protocole MESI)              (protocole MESI)
```

Chaque invalidation force un rechargement de la ligne de cache depuis un niveau de cache supérieur (L2 ou L3), ce qui coûte typiquement **30 à 100 cycles** par occurrence.

#### 2. Surcharge d'administration des threads (*Fork/Join overhead*)

OpenMP crée et synchronise une équipe de threads **à chaque itération** de la boucle principale (500 fois). Le modèle Fork/Join implique :

1. **Fork** : création (ou réveil) de 12 threads au début de la région parallèle
2. **Distribution** : attribution des itérations aux threads
3. **Calcul** : exécution parallèle
4. **Join** : barrière de synchronisation — tous les threads attendent le plus lent
5. **Réduction** : agrégation du compteur `cpteur_food`

Avec seulement 5 000 fourmis et un temps de mouvement séquentiel de ~0.6 ms par itération, la charge de travail par thread est d'environ **0.6 ms / 12 = 50 µs**. Or, le coût d'un Fork/Join sur OpenMP est typiquement de **10 à 50 µs** selon la plateforme. Le rapport travail/surcharge est donc proche de 1, ce qui signifie que les threads passent autant de temps à se synchroniser qu'à calculer.

#### 3. Déséquilibre de charge (*Load Imbalance*)

Bien que `schedule(dynamic)` ait été utilisé pour gérer la variabilité du temps par fourmi, ce scheduler introduit lui-même une surcharge : chaque thread doit accéder à un compteur atomique partagé pour obtenir sa prochaine itération, ce qui génère de la contention supplémentaire.

De plus, les fourmis situées dans des zones à faible coût de terrain (valeurs proches de 0) effectuent **davantage de sous-pas de temps** que celles dans les zones coûteuses (valeurs proches de 1). Cela crée un déséquilibre naturel que `schedule(dynamic)` atténue partiellement, mais au prix d'une surcharge d'ordonnancement.

#### Résumé des facteurs de ralentissement

| Facteur | Impact estimé | Remédiable ? |
|---|---|---|
| Contention atomique sur `mark_pheronome()` | **Élevé** | Oui, par grilles privées (coût mémoire) |
| Faux partage sur les lignes de cache | **Moyen** | Oui, par padding ou grilles privées |
| Surcharge Fork/Join (500 itérations × 12 threads) | **Moyen** | Partiellement, par augmentation de la charge |
| Surcharge `schedule(dynamic)` | **Faible** | Oui, par `schedule(static)` si charge uniforme |

#### Estimation de la charge minimale pour un speedup positif

Pour que la parallélisation OpenMP soit rentable, il faut que le temps de calcul par thread **dépasse largement** la surcharge d'administration. Avec 12 threads et une surcharge estimée à ~30 µs par Fork/Join :

$$T_{\text{mouvement}} > 12 \times 30\,\mu\text{s} \times N_{\text{itérations}} = 12 \times 30 \times 500 = 180\,\text{ms}$$

Le temps de mouvement séquentiel est de 323 ms, ce qui est juste au-dessus du seuil. Mais les opérations atomiques annulent le gain théorique. Avec **50 000 fourmis** (mouvement séquentiel ~3.2 s), le rapport serait bien plus favorable.

### Conclusion de la Phase 3

| Critère | Résultat |
|---|---|
| Parallélisation OpenMP | [OK] Implémentée (mouvement + évaporation) |
| Protection des données | [OK] `atomic write` + `reduction` |
| Speedup (5 000 fourmis, 12 threads) | **0.82×** [FAIL] (ralentissement de 22%) |
| Cause principale | Contention atomique + surcharge Fork/Join |

**Le parallélisme en mémoire partagée « naïf » atteint ses limites** dans ce contexte. Les accès concurrents en écriture sur la grille de phéromones créent une contention qui annule — et dépasse — le gain de la distribution du calcul. De plus, la charge de travail (5 000 fourmis) est insuffisante pour amortir les coûts d'administration d'OpenMP.

Pour obtenir un speedup réel, deux pistes se présentent :

1. **Augmenter la charge de travail** : avec 50 000+ fourmis, le ratio calcul/overhead deviendrait favorable
2. **Isoler les zones mémoire** : plutôt que de partager la grille de phéromones entre tous les threads (mémoire partagée), on peut donner à chaque processus **sa propre copie** ou **sa propre portion** de la grille. C'est exactement ce que propose le paradigme de **mémoire distribuée (MPI)**, objet des phases suivantes.

Cette expérience illustre un enseignement fondamental en calcul haute performance : **la parallélisation n'est pas toujours synonyme d'accélération**. La performance dépend du rapport entre le travail utile parallélisable et les surcoûts induits (synchronisation, communication, contention mémoire).


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
| Compilateur MPI | `mpic++ --version` | `g++ (Ubuntu 13.3.0-6ubuntu2~24.04) 13.3.0` [OK] |
| OpenMPI (runtime) | `dpkg -l openmpi-bin` | Version 4.1.6-7ubuntu2, déjà installé [OK] |
| MPICH (alternative) | `sudo apt install mpich` | Version 4.2.0-5build3, installé en complément [OK] |

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
| 12 | [FAIL] Échec | — | — | — | — |

**Tableau de speedup :**

| Nb Procs | Temps Total (s) | Speedup ($T_1/T_N$) | Efficacité (Speedup/$N$) |
|:---:|:---:|:---:|:---:|
| 1 | 1.704 | 1.00 | 100% |
| 2 | 1.915 | 0.89 [FAIL] | 44.5% |
| 4 | 6.910 | 0.25 [FAIL] | 6.1% |
| 8 | 8.892 | 0.19 [FAIL] | 2.4% |

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
| 2 | 0.89 [FAIL] | **1.21** [OK] | **1.42** [OK] |
| 4 | 0.25 [FAIL] | 0.70 [FAIL] | 0.91 [FAIL] |
| 8 | 0.19 [FAIL] | 0.64 [FAIL] | 0.91 [FAIL] |
| 12* | — | 0.50 [FAIL] | 0.80 [FAIL] |

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
Speedup: 0.89 [FAIL]            Speedup: 1.21 [OK]           Speedup: 1.42 [OK]
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
| 2 | 0.89 [FAIL] | **1.21** [OK] | **1.42** [OK] |
| 4 | 0.25 [FAIL] | 0.70 [FAIL] | 0.91 [FAIL] |
| 8 | 0.19 [FAIL] | 0.64 [FAIL] | 0.91 [FAIL] |
| 12* | — | 0.50 [FAIL] | 0.80 [FAIL] |

#### Bilan

| Critère | Résultat |
|---|---|
| **Implémentation MPI** | [OK] Fonctionnelle et correcte |
| **Correctesse** | [OK] Résultats cohérents sur toutes les configurations |
| **Meilleur speedup obtenu** | **1.42×** (100k fourmis, 2 processus) |
| **Scalabilité** | [FAIL] Ne passe pas à l'échelle au-delà de 2 processus |
| **Cause principale** | `MPI_Allreduce` sur ~16.8 Mo à chaque itération |
| **Conformité au sujet** | [OK] Les limites annoncées sont confirmées expérimentalement |

#### Leçon retenue

La Stratégie 1 illustre un **anti-pattern classique** en calcul distribué : elle échange **trop de données** par rapport au calcul qu'elle effectue. Le `MPI_Allreduce` sur la grille entière est un coût **incompressible** qui croît avec le nombre de processus, tandis que le gain de calcul (division des fourmis) est **limité** par la taille du problème.

Pour que cette stratégie soit rentable, deux conditions doivent être réunies :
1. Le nombre de fourmis doit être **très grand** (>>100 000) pour que le mouvement domine le temps total
2. La carte doit être **petite** pour que le volume du `MPI_Allreduce` reste gérable

Ces deux conditions sont en tension : une carte petite limite l'intérêt de la simulation, tandis qu'une grande population de fourmis sur une petite carte sature rapidement les phéromones.

La **Stratégie 2 (décomposition de domaine)**, objet de la Phase 5, résout structurellement ce problème en n'échangeant que les **cellules fantômes aux frontières** (~32 Ko au lieu de 16.8 Mo), soit une réduction du volume de communication de **~500×**.


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
| 1 | 100 000 | 100 000 | [OK] Correct |
| 2 | 100 000 | 100 000 | [OK] Correct |
| 4 | 100 000 | 100 000 | [OK] Correct |
| 8 | 100 000 | 100 000 | [OK] Correct |
| 12 | 100 000 | 100 000 | [OK] Correct |

**Aucune fourmi perdue ni dupliquée** lors de la migration. Le protocole de transfert swap-and-pop + `MPI_Sendrecv` est correct et conservatif.

---

### Analyse des résultats

#### Speedup global

| Nb Procs | Temps Total (s) | Speedup Strat. 2 | Speedup Strat. 1 (rappel) |
|:---:|:---:|:---:|:---:|
| 1 | 16.651 | 1.00 | 1.00 |
| 2 | 7.851 | **2.12** [OK] | 1.42 |
| 4 | 8.919 | **1.87** [OK] | 0.91 |
| 8 | 6.820 | **2.44** [OK] | 0.91 |
| 12* | 6.503 | **2.56** [OK] | 0.80 |

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
| **Passe à l'échelle ?** | [FAIL] Non (ralentit dès 4P) | [OK] Oui (accélère jusqu'à 12P) |
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