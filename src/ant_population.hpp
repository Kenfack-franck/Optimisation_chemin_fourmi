#ifndef _ANT_POPULATION_HPP_
#define _ANT_POPULATION_HPP_

#include <vector>
#include <random>
#include <algorithm>
#include "basic_types.hpp"
#include "fractal_land.hpp"
#include "pheronome.hpp"
#include "rand_generator.hpp"

class AntPopulation {
public:
    // Structure of Arrays (SoA) : Les données sont séparées en vecteurs contigus
    std::vector<int> m_x;
    std::vector<int> m_y;
    std::vector<int> m_loaded; // 0 = false, 1 = true (int pour alignement mémoire)
    std::vector<std::size_t> m_seeds;
    
    double m_eps; // Coefficient d'exploration

    // Constructeur
    AntPopulation(size_t nb_ants, const fractal_land& land, double eps, size_t global_seed) 
        : m_eps(eps) 
    {
        m_x.resize(nb_ants);
        m_y.resize(nb_ants);
        m_loaded.resize(nb_ants, 0); // Tout le monde commence non chargé
        m_seeds.resize(nb_ants);

        // Initialisation aléatoire
        RandomGenerator pos_gen(global_seed, 0, land.dimensions()-1);
        // On donne une graine unique à chaque fourmi pour qu'elles aient des parcours différents
        // On utilise un générateur simple pour initialiser les graines
        std::size_t current_seed = global_seed; 
        
        for(size_t i=0; i<nb_ants; ++i) {
            // Position aléatoire initiale
            // Note: on utilise l'index i pour varier la position initiale
            m_x[i] = static_cast<int>(pos_gen(i, 0));
            m_y[i] = static_cast<int>(pos_gen(0, i));
            
            // Génération d'une seed propre à la fourmi
            // Formule simple pour disperser les graines
            current_seed = (1664525 * current_seed + 1013904223) % 0xFFFFFFFF;
            m_seeds[i] = current_seed; 
        }
    }

    // // Méthode vectorisée pour faire avancer TOUTES les fourmis
    // void advance_all(pheronome& phen, const fractal_land& land, 
    //                  const position_t& pos_food, const position_t& pos_nest, 
    //                  std::size_t& cpteur_food) 
    // {
    //     const size_t nb_ants = m_x.size();

    //     // Cette boucle est maintenant candidate idéale pour la vectorisation SIMD
    //     // et plus tard pour OpenMP
    //     for (size_t i = 0; i < nb_ants; ++i) {
    //         advance_one_ant(i, phen, land, pos_food, pos_nest, cpteur_food);
    //     }
    // }

    

     // Méthode vectorisée pour faire avancer TOUTES les fourmis
    void advance_all(pheronome& phen, const fractal_land& land, 
                     const position_t& pos_food, const position_t& pos_nest, 
                     std::size_t& cpteur_food) 
    {
        const size_t nb_ants = m_x.size();

        // Cette boucle est maintenant candidate idéale pour la vectorisation SIMD
        // et plus tard pour OpenMP
        for (size_t i = 0; i < nb_ants; ++i) {
            advance_one_ant(i, phen, land, pos_food, pos_nest, cpteur_food);
        }
    }

private:
    // Logique de mouvement pour UNE fourmi (identifiée par son index i)
    // C'est exactement la même logique que ant.cpp, mais adaptée aux vecteurs
    void advance_one_ant(size_t i, pheronome& phen, const fractal_land& land, 
                         const position_t& pos_food, const position_t& pos_nest, 
                         std::size_t& cpteur_food) 
    {
        // Générateurs locaux utilisant la graine stockée de la fourmi i
        auto ant_choice = [&]() { return rand_double(0., 1., m_seeds[i]); };
        auto dir_choice = [&]() { return rand_int32(1, 4, m_seeds[i]); };

        double consumed_time = 0.;

        while (consumed_time < 1.) {
            int ind_pher = (m_loaded[i] == 1 ? 1 : 0);
            double choix = ant_choice();
            
            int old_x = m_x[i];
            int old_y = m_y[i];
            int new_x = old_x;
            int new_y = old_y;

            // Lecture des phéromones voisins
            // Attention aux accès hors limites potentiels (gérés par pheronome::operator())
            double p_left   = phen(old_x - 1, old_y)[ind_pher];
            double p_right  = phen(old_x + 1, old_y)[ind_pher];
            double p_up     = phen(old_x, old_y - 1)[ind_pher];
            double p_down   = phen(old_x, old_y + 1)[ind_pher];

            double max_phen = std::max({p_left, p_right, p_up, p_down});

            // Décision de mouvement
            if ((choix > m_eps) || (max_phen <= 0.)) {
                // Mouvement aléatoire tant qu'on tombe sur un mur (-1)
                do {
                    new_x = old_x; 
                    new_y = old_y;
                    int d = dir_choice();
                    if (d == 1) new_x -= 1;
                    if (d == 2) new_y -= 1;
                    if (d == 3) new_x += 1;
                    if (d == 4) new_y += 1;
                    // On vérifie le phéromone de la case cible. 
                    // Si -1, c'est un mur ou bordure.
                } while (phen(new_x, new_y)[ind_pher] == -1);
            } else {
                // Mouvement guidé par les phéromones
                if (p_left == max_phen)       new_x -= 1;
                else if (p_right == max_phen) new_x += 1;
                else if (p_up == max_phen)    new_y -= 1;
                else                          new_y += 1;
            }

            // Mise à jour du temps et marquage
            consumed_time += land(new_x, new_y);
            
            // Mise à jour de la position réelle de la fourmi
            m_x[i] = new_x;
            m_y[i] = new_y;

            // Dépôt de phéromones
            phen.mark_pheronome({new_x, new_y});

            // Gestion de la nourriture et du nid
            if (new_x == pos_nest.x && new_y == pos_nest.y) {
                if (m_loaded[i]) {
                    cpteur_food += 1;
                }
                m_loaded[i] = 0; // Déchargée
            }
            if (new_x == pos_food.x && new_y == pos_food.y) {
                m_loaded[i] = 1; // Chargée
            }
        }
    }
};

#endif