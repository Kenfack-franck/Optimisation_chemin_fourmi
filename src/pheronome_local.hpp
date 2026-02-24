#ifndef _PHERONOME_LOCAL_HPP_
#define _PHERONOME_LOCAL_HPP_

#include <algorithm>
#include <vector>
#include <iostream>
#include <cassert>
#include <omp.h>
#include "basic_types.hpp"
#include "domain_decomposition.hpp"

/**
 * @brief Grille de phéromones locale pour la décomposition de domaine.
 * 
 * Chaque processus ne stocke que ses lignes locales + 1 ghost row en haut/bas.
 * Le stockage est plat : [V1, V2] par cellule, lignes consécutives.
 * 
 * Indexation locale :
 *   local_i = 0                : ghost top
 *   local_i = 1..local_nrows   : données réelles
 *   local_i = local_nrows + 1  : ghost bottom
 *   local_j = 0                : ghost left (bordure)
 *   local_j = 1..global_dim    : données réelles
 *   local_j = global_dim + 1   : ghost right (bordure)
 */
class PheronomeLocal {
public:
    using size_t = unsigned long;

    PheronomeLocal(const DomainDecomposition& domain,
                   const position_t& pos_food, const position_t& pos_nest,
                   double alpha = 0.7, double beta = 0.9999)
        : m_domain(domain),
          m_alpha(alpha), m_beta(beta),
          m_pos_food(pos_food), m_pos_nest(pos_nest),
          m_data(domain.local_phen_size(), 0.0),
          m_buffer(domain.local_phen_size(), 0.0)
    {
        // Initialiser nourriture et nid s'ils sont dans notre domaine
        if (is_food_local()) {
            get_val_local(food_local_i(), food_local_j(), 0) = 1.0;
        }
        if (is_nest_local()) {
            get_val_local(nest_local_i(), nest_local_j(), 1) = 1.0;
        }
        
        cl_update();
        m_buffer = m_data;
    }

    // =================================================================
    // ACCESSEURS avec coordonnées GLOBALES (interface compatible)
    // =================================================================
    
    /**
     * @brief Accès par coordonnées globales. Renvoie pointeur vers [V1, V2].
     * Gère automatiquement la conversion global → local.
     */
    double* operator()(size_t global_i, size_t global_j) {
        size_t local_i = global_i - m_domain.local_row_start + 1;
        size_t local_j = global_j + 1;
        return &m_data[(local_i * m_domain.stride_cols + local_j) * 2];
    }

    const double* operator()(size_t global_i, size_t global_j) const {
        size_t local_i = global_i - m_domain.local_row_start + 1;
        size_t local_j = global_j + 1;
        return &m_data[(local_i * m_domain.stride_cols + local_j) * 2];
    }

    /**
     * @brief Accès pour les cases voisines qui peuvent être dans les ghost cells.
     * Accepte global_i dans [local_row_start-1, local_row_end] (ghost inclus).
     */
    const double* safe_access(int global_i, int global_j) const {
        int local_i = global_i - static_cast<int>(m_domain.local_row_start) + 1;
        int local_j = global_j + 1;
        
        // Vérification des bornes (les ghost cells couvrent local_i = 0 et local_nrows+1)
        if (local_i < 0 || local_i >= static_cast<int>(m_domain.local_stride_rows) ||
            local_j < 0 || local_j >= static_cast<int>(m_domain.stride_cols)) {
            // Hors limites → renvoyer un pointeur vers des valeurs -1 (mur)
            static double wall[2] = {-1.0, -1.0};
            return wall;
        }
        return &m_data[(local_i * m_domain.stride_cols + local_j) * 2];
    }

    // =================================================================
    // ACCESSEURS avec coordonnées LOCALES
    // =================================================================
    
    double& get_val_local(size_t local_i, size_t local_j, int type) {
        return m_data[(local_i * m_domain.stride_cols + local_j) * 2 + type];
    }

    const double& get_val_local(size_t local_i, size_t local_j, int type) const {
        return m_data[(local_i * m_domain.stride_cols + local_j) * 2 + type];
    }

    // =================================================================
    // OPÉRATIONS
    // =================================================================

    void do_evaporation() {
        #pragma omp parallel for
        for (size_t i = 1; i <= m_domain.local_nrows; ++i) {
            for (size_t j = 1; j <= m_domain.global_dim; ++j) {
                size_t idx = (i * m_domain.stride_cols + j) * 2;
                m_buffer[idx]     *= m_beta;
                m_buffer[idx + 1] *= m_beta;
            }
        }
    }

    void mark_pheronome(const position_t& pos) {
        int global_i = pos.x;
        int global_j = pos.y;
        
        // Lecture des voisins (peuvent être dans les ghost cells)
        const double* left  = safe_access(global_i - 1, global_j);
        const double* right = safe_access(global_i + 1, global_j);
        const double* up    = safe_access(global_i, global_j - 1);
        const double* down  = safe_access(global_i, global_j + 1);

        double v1_left  = std::max(left[0],  0.);
        double v1_right = std::max(right[0], 0.);
        double v1_up    = std::max(up[0],    0.);
        double v1_down  = std::max(down[0],  0.);
        
        double v2_left  = std::max(left[1],  0.);
        double v2_right = std::max(right[1], 0.);
        double v2_up    = std::max(up[1],    0.);
        double v2_down  = std::max(down[1],  0.);

        double val1 = m_alpha * std::max({v1_left, v1_right, v1_up, v1_down}) +
                      (1 - m_alpha) * 0.25 * (v1_left + v1_right + v1_up + v1_down);
        
        double val2 = m_alpha * std::max({v2_left, v2_right, v2_up, v2_down}) +
                      (1 - m_alpha) * 0.25 * (v2_left + v2_right + v2_up + v2_down);

        size_t local_i = global_i - m_domain.local_row_start + 1;
        size_t local_j = global_j + 1;
        size_t idx = (local_i * m_domain.stride_cols + local_j) * 2;
        
        #pragma omp atomic write
        m_buffer[idx] = val1;
        
        #pragma omp atomic write
        m_buffer[idx + 1] = val2;
    }

    void update() {
        m_data.swap(m_buffer);
        cl_update();
        
        // Remettre les valeurs fixes
        if (is_food_local()) {
            get_val_local(food_local_i(), food_local_j(), 0) = 1.0;
        }
        if (is_nest_local()) {
            get_val_local(nest_local_i(), nest_local_j(), 1) = 1.0;
        }
    }

    /**
     * @brief Échange les ghost cells avec les processus voisins
     */
    void exchange_ghosts() {
        m_domain.exchange_ghost_cells(m_data);
    }

    // Accesseur pour les données brutes (pour debug/rendu)
    std::vector<double>& get_data() { return m_data; }
    const std::vector<double>& get_data() const { return m_data; }
    size_t get_total_size() const { return m_data.size(); }

private:
    void cl_update() {
        // Bordures gauche/droite à -1
        for (size_t i = 0; i < m_domain.local_stride_rows; ++i) {
            // Colonne gauche (j=0)
            size_t idx_left = (i * m_domain.stride_cols) * 2;
            m_data[idx_left] = -1.; m_data[idx_left + 1] = -1.;
            // Colonne droite (j = global_dim + 1)
            size_t idx_right = (i * m_domain.stride_cols + m_domain.global_dim + 1) * 2;
            m_data[idx_right] = -1.; m_data[idx_right + 1] = -1.;
        }
        
        // Si on est le premier rank, ghost top = -1 (bord supérieur du domaine)
        if (m_domain.rank_top == MPI_PROC_NULL) {
            for (size_t j = 0; j < m_domain.stride_cols; ++j) {
                size_t idx = j * 2;
                m_data[idx] = -1.; m_data[idx + 1] = -1.;
            }
        }
        
        // Si on est le dernier rank, ghost bottom = -1 (bord inférieur)
        if (m_domain.rank_bottom == MPI_PROC_NULL) {
            for (size_t j = 0; j < m_domain.stride_cols; ++j) {
                size_t idx = ((m_domain.local_nrows + 1) * m_domain.stride_cols + j) * 2;
                m_data[idx] = -1.; m_data[idx + 1] = -1.;
            }
        }
    }

    // Helpers pour savoir si nourriture/nid sont dans ce domaine
    bool is_food_local() const {
        return m_domain.is_local(m_pos_food.x, m_pos_food.y);
    }
    bool is_nest_local() const {
        return m_domain.is_local(m_pos_nest.x, m_pos_nest.y);
    }
    size_t food_local_i() const { return m_pos_food.x - m_domain.local_row_start + 1; }
    size_t food_local_j() const { return static_cast<size_t>(m_pos_food.y) + 1; }
    size_t nest_local_i() const { return m_pos_nest.x - m_domain.local_row_start + 1; }
    size_t nest_local_j() const { return static_cast<size_t>(m_pos_nest.y) + 1; }

    const DomainDecomposition& m_domain;
    double m_alpha, m_beta;
    position_t m_pos_food, m_pos_nest;
    std::vector<double> m_data;
    std::vector<double> m_buffer;
};

#endif