## Phase 3 : Parallélisme en Mémoire Partagée (OpenMP)

### Objectif de la phase

La Phase 2 a préparé le terrain en restructurant les données au format SoA. L'objectif de cette phase est d'exploiter les **multiples cœurs** de la machine de test pour accélérer les boucles de calcul, en utilisant les directives de compilation **OpenMP** (`#pragma omp`).

OpenMP est un standard de parallélisation en **mémoire partagée** : tous les threads accèdent au même espace d'adressage (même RAM). C'est le modèle le plus simple à mettre en œuvre car il ne nécessite que l'ajout de quelques directives au-dessus des boucles existantes. Cependant, il impose de gérer les **accès concurrents** lorsque plusieurs threads écrivent au même emplacement mémoire.

### Identification des sections parallélisables

D'après le profiling de la Phase 1, trois sections de code sont candidates à la parallélisation :

| Section | % du temps | Parallélisable ? | Difficulté |
|---|---|---|---|
| Mouvement des fourmis (`advance_all`) | 75-84% | ✅ Oui (boucle sur les fourmis) | **Élevée** — écriture concurrente sur la grille de phéromones |
| Évaporation (`do_evaporation`) | 15-25% | ✅ Oui (boucle sur la grille) | **Faible** — chaque cellule est indépendante |
| Mise à jour (`update`) | <1% | ❌ Non pertinent | Négligeable |

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
         │   (privé ✅)                  │   (privé ✅)
         │                              │
         ├── Lecture phen(x±1, y±1)     ├── Lecture phen(x±1, y±1)
         │   (partagé, lecture seule ✅) │   (partagé, lecture seule ✅)
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
| Parallélisation OpenMP | ✅ Implémentée (mouvement + évaporation) |
| Protection des données | ✅ `atomic write` + `reduction` |
| Speedup (5 000 fourmis, 12 threads) | **0.82×** ❌ (ralentissement de 22%) |
| Cause principale | Contention atomique + surcharge Fork/Join |

**Le parallélisme en mémoire partagée « naïf » atteint ses limites** dans ce contexte. Les accès concurrents en écriture sur la grille de phéromones créent une contention qui annule — et dépasse — le gain de la distribution du calcul. De plus, la charge de travail (5 000 fourmis) est insuffisante pour amortir les coûts d'administration d'OpenMP.

Pour obtenir un speedup réel, deux pistes se présentent :

1. **Augmenter la charge de travail** : avec 50 000+ fourmis, le ratio calcul/overhead deviendrait favorable
2. **Isoler les zones mémoire** : plutôt que de partager la grille de phéromones entre tous les threads (mémoire partagée), on peut donner à chaque processus **sa propre copie** ou **sa propre portion** de la grille. C'est exactement ce que propose le paradigme de **mémoire distribuée (MPI)**, objet des phases suivantes.

Cette expérience illustre un enseignement fondamental en calcul haute performance : **la parallélisation n'est pas toujours synonyme d'accélération**. La performance dépend du rapport entre le travail utile parallélisable et les surcoûts induits (synchronisation, communication, contention mémoire).