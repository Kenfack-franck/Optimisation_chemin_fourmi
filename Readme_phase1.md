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
| 1.1 | Compilation et validation | ✅ Code fonctionnel, environnement opérationnel |
| 1.2 | Mode benchmark sans GUI | ✅ Mesures de temps fiables, 500 itérations déterministes |
| 1.3 | Profiling initial | ✅ Baseline établie : 0.38 s total, 75% mouvement, 25% évaporation |

**Baseline de référence :** 5 000 fourmis, 500 itérations → **0.381 secondes** (0.76 ms/itération)

Nous disposons maintenant de tous les éléments pour passer à l'optimisation du code : une version de référence validée, un protocole de mesure fiable, et une connaissance précise de la répartition du temps de calcul.