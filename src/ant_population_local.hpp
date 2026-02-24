#ifndef _ANT_POPULATION_LOCAL_HPP_
#define _ANT_POPULATION_LOCAL_HPP_

#include <vector>
#include <algorithm>
#include <mpi.h>
#include "basic_types.hpp"
#include "fractal_land.hpp"
#include "pheronome_local.hpp"
#include "domain_decomposition.hpp"
#include "rand_generator.hpp"

/**
 * @brief Population de fourmis pour la décomposition de domaine.
 * 
 * Chaque processus ne gère que les fourmis présentes dans son sous-domaine.
 * Si une fourmi sort du domaine local, elle est envoyée au processus voisin.
 */
class AntPopulationLocal {
public:
    // SoA
    std::vector<int> m_x;
    std::vector<int> m_y;
    std::vector<int> m_loaded;
    std::vector<std::size_t> m_seeds;
    
    double m_eps;
    const DomainDecomposition& m_domain;

    /**
     * @brief Constructeur : initialise des fourmis UNIQUEMENT dans le domaine local
     */
    AntPopulationLocal(size_t total_nb_ants, const fractal_land& land, 
                       double eps, size_t global_seed,
                       const DomainDecomposition& domain)
        : m_eps(eps), m_domain(domain)
    {
        // Chaque processus crée sa part de fourmis dans son domaine
        size_t ants_per_proc = total_nb_ants / domain.nb_procs;
        size_t remainder = total_nb_ants % domain.nb_procs;
        size_t local_ants = ants_per_proc + (static_cast<size_t>(domain.rank) < remainder ? 1 : 0);
        
        m_x.resize(local_ants);
        m_y.resize(local_ants);
        m_loaded.resize(local_ants, 0);
        m_seeds.resize(local_ants);

        // Graine unique par processus
        std::size_t current_seed = global_seed + domain.rank * 100003;
        
        for (size_t i = 0; i < local_ants; ++i) {
            // Position aléatoire DANS le domaine local
            current_seed = (1664525 * current_seed + 1013904223) % 0xFFFFFFFF;
            m_x[i] = static_cast<int>(domain.local_row_start + 
                     (current_seed % domain.local_nrows));
            
            current_seed = (1664525 * current_seed + 1013904223) % 0xFFFFFFFF;
            m_y[i] = static_cast<int>(current_seed % domain.global_dim);
            
            current_seed = (1664525 * current_seed + 1013904223) % 0xFFFFFFFF;
            m_seeds[i] = current_seed;
        }
    }

    size_t size() const { return m_x.size(); }

    /**
     * @brief Fait avancer toutes les fourmis locales.
     * Les fourmis qui sortent du domaine sont marquées pour migration.
     */
    void advance_all(PheronomeLocal& phen, const fractal_land& land,
                     const position_t& pos_food, const position_t& pos_nest,
                     std::size_t& cpteur_food)
    {
        const size_t nb_ants = m_x.size();
        
        for (size_t i = 0; i < nb_ants; ++i) {
            advance_one_ant(i, phen, land, pos_food, pos_nest, cpteur_food);
        }
    }

    /**
     * @brief Migre les fourmis hors domaine vers les processus voisins.
     * 
     * Protocole : 
     *   1. Identifier les fourmis hors domaine (trop haut → rank_top, trop bas → rank_bottom)
     *   2. Envoyer au voisin : [x, y, loaded, seed] par fourmi
     *   3. Recevoir les fourmis entrantes
     *   4. Supprimer les fourmis sorties, ajouter les entrantes
     */
    void migrate_ants() {
        // --- Trier les fourmis à migrer ---
        // Données à envoyer : 4 doubles par fourmi [x, y, loaded, seed]
        std::vector<double> send_top, send_bottom;
        std::vector<size_t> to_remove; // Indices des fourmis à supprimer
        
        for (size_t i = 0; i < m_x.size(); ++i) {
            if (m_x[i] < static_cast<int>(m_domain.local_row_start)) {
                // Fourmi sort par le haut → envoyer au rank_top
                send_top.push_back(static_cast<double>(m_x[i]));
                send_top.push_back(static_cast<double>(m_y[i]));
                send_top.push_back(static_cast<double>(m_loaded[i]));
                send_top.push_back(static_cast<double>(m_seeds[i]));
                to_remove.push_back(i);
            } else if (m_x[i] >= static_cast<int>(m_domain.local_row_end)) {
                // Fourmi sort par le bas → envoyer au rank_bottom
                send_bottom.push_back(static_cast<double>(m_x[i]));
                send_bottom.push_back(static_cast<double>(m_y[i]));
                send_bottom.push_back(static_cast<double>(m_loaded[i]));
                send_bottom.push_back(static_cast<double>(m_seeds[i]));
                to_remove.push_back(i);
            }
        }
        
        // --- Communiquer le nombre de fourmis à envoyer/recevoir ---
        int nb_send_top = static_cast<int>(send_top.size());
        int nb_send_bottom = static_cast<int>(send_bottom.size());
        int nb_recv_top = 0, nb_recv_bottom = 0;
        
        MPI_Status status;
        
        // Échanger les compteurs avec les voisins
        MPI_Sendrecv(&nb_send_top, 1, MPI_INT, m_domain.rank_top, 10,
                     &nb_recv_bottom, 1, MPI_INT, m_domain.rank_bottom, 10,
                     MPI_COMM_WORLD, &status);
        
        MPI_Sendrecv(&nb_send_bottom, 1, MPI_INT, m_domain.rank_bottom, 11,
                     &nb_recv_top, 1, MPI_INT, m_domain.rank_top, 11,
                     MPI_COMM_WORLD, &status);
        
        // --- Échanger les données des fourmis ---
        std::vector<double> recv_top(nb_recv_top);
        std::vector<double> recv_bottom(nb_recv_bottom);
        
        MPI_Sendrecv(send_top.data(), nb_send_top, MPI_DOUBLE, m_domain.rank_top, 20,
                     recv_bottom.data(), nb_recv_bottom, MPI_DOUBLE, m_domain.rank_bottom, 20,
                     MPI_COMM_WORLD, &status);
        
        MPI_Sendrecv(send_bottom.data(), nb_send_bottom, MPI_DOUBLE, m_domain.rank_bottom, 21,
                     recv_top.data(), nb_recv_top, MPI_DOUBLE, m_domain.rank_top, 21,
                     MPI_COMM_WORLD, &status);
        
        // --- Supprimer les fourmis parties (en ordre inverse pour ne pas casser les indices) ---
        std::sort(to_remove.rbegin(), to_remove.rend());
        for (size_t idx : to_remove) {
            // Swap avec le dernier élément et pop_back (O(1))
            size_t last = m_x.size() - 1;
            if (idx != last) {
                std::swap(m_x[idx], m_x[last]);
                std::swap(m_y[idx], m_y[last]);
                std::swap(m_loaded[idx], m_loaded[last]);
                std::swap(m_seeds[idx], m_seeds[last]);
            }
            m_x.pop_back();
            m_y.pop_back();
            m_loaded.pop_back();
            m_seeds.pop_back();
        }
        
        // --- Ajouter les fourmis reçues ---
        auto add_ants = [&](const std::vector<double>& recv) {
            for (size_t k = 0; k + 3 < recv.size(); k += 4) {
                m_x.push_back(static_cast<int>(recv[k]));
                m_y.push_back(static_cast<int>(recv[k + 1]));
                m_loaded.push_back(static_cast<int>(recv[k + 2]));
                m_seeds.push_back(static_cast<std::size_t>(recv[k + 3]));
            }
        };
        add_ants(recv_top);
        add_ants(recv_bottom);
    }

private:
    void advance_one_ant(size_t i, PheronomeLocal& phen, const fractal_land& land,
                         const position_t& pos_food, const position_t& pos_nest,
                         std::size_t& cpteur_food)
    {
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

            // Lecture des phéromones voisins via safe_access (gère les ghost cells)
            const double* p_left  = phen.safe_access(old_x - 1, old_y);
            const double* p_right = phen.safe_access(old_x + 1, old_y);
            const double* p_up    = phen.safe_access(old_x, old_y - 1);
            const double* p_down  = phen.safe_access(old_x, old_y + 1);

            double max_phen = std::max({p_left[ind_pher], p_right[ind_pher],
                                        p_up[ind_pher], p_down[ind_pher]});

            if ((choix > m_eps) || (max_phen <= 0.)) {
                do {
                    new_x = old_x;
                    new_y = old_y;
                    int d = dir_choice();
                    if (d == 1) new_x -= 1;
                    if (d == 2) new_y -= 1;
                    if (d == 3) new_x += 1;
                    if (d == 4) new_y += 1;
                } while (phen.safe_access(new_x, new_y)[ind_pher] == -1);
            } else {
                if (p_left[ind_pher] == max_phen)        new_x -= 1;
                else if (p_right[ind_pher] == max_phen)  new_x += 1;
                else if (p_up[ind_pher] == max_phen)     new_y -= 1;
                else                                      new_y += 1;
            }

            consumed_time += land(new_x, new_y);
            m_x[i] = new_x;
            m_y[i] = new_y;

            // Marquage phéromone seulement si la fourmi est encore dans le domaine local
            if (m_domain.is_local(new_x, new_y)) {
                phen.mark_pheronome({new_x, new_y});
            }

            if (new_x == pos_nest.x && new_y == pos_nest.y) {
                if (m_loaded[i]) {
                    cpteur_food += 1;
                }
                m_loaded[i] = 0;
            }
            if (new_x == pos_food.x && new_y == pos_food.y) {
                m_loaded[i] = 1;
            }
        }
    }
};

#endif