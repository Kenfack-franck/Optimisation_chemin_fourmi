#ifndef _PHERONOME_HPP_
#define _PHERONOME_HPP_
#include <algorithm>
#include <vector>
#include <iostream>
#include <cassert>
#include <omp.h>
#include "basic_types.hpp"

class pheronome {
public:
    using size_t = unsigned long;

    pheronome(size_t dim, const position_t& pos_food, const position_t& pos_nest,
              double alpha = 0.7, double beta = 0.9999)
        : m_dim(dim),
          m_stride(dim + 2),
          m_alpha(alpha), m_beta(beta),
          // Taille * 2 car 2 valeurs (V1, V2) par cellule
          m_data(m_stride * m_stride * 2, 0.0), 
          m_buffer(m_stride * m_stride * 2, 0.0),
          m_pos_nest(pos_nest),
          m_pos_food(pos_food) 
    {
        // Init positions spécifiques
        get_val(pos_food.x, pos_food.y, 0) = 1.;
        get_val(pos_nest.x, pos_nest.y, 1) = 1.;
        
        cl_update();
        m_buffer = m_data;
    }

    // Accesseur helper pour simplifier la lecture (i, j, type)
    // type 0 = nourriture (V1), type 1 = nid (V2)
    double& get_val(size_t i, size_t j, int type) {
        return m_data[((i + 1) * m_stride + (j + 1)) * 2 + type];
    }

    const double& get_val(size_t i, size_t j, int type) const {
        return m_data[((i + 1) * m_stride + (j + 1)) * 2 + type];
    }

    // Opérateur () pour compatibilité avec le reste du code
    // Renvoie un pointeur vers la paire de double [v1, v2]
    double* operator()(size_t i, size_t j) {
        return &m_data[((i + 1) * m_stride + (j + 1)) * 2];
    }
    
    const double* operator()(size_t i, size_t j) const {
        return &m_data[((i + 1) * m_stride + (j + 1)) * 2];
    }

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

    void mark_pheronome(const position_t& pos) {
        size_t i = pos.x;
        size_t j = pos.y;
        
        // Lecture voisins (code un peu verbeux pour accès plat)
        double v1_left   = std::max(get_val(i-1, j, 0), 0.);
        double v2_left   = std::max(get_val(i-1, j, 1), 0.);
        double v1_right  = std::max(get_val(i+1, j, 0), 0.);
        double v2_right  = std::max(get_val(i+1, j, 1), 0.);
        double v1_up     = std::max(get_val(i, j-1, 0), 0.);
        double v2_up     = std::max(get_val(i, j-1, 1), 0.);
        double v1_down   = std::max(get_val(i, j+1, 0), 0.);
        double v2_down   = std::max(get_val(i, j+1, 1), 0.);

        double val1 = m_alpha * std::max({v1_left, v1_right, v1_up, v1_down}) +
                      (1 - m_alpha) * 0.25 * (v1_left + v1_right + v1_up + v1_down);
        
        double val2 = m_alpha * std::max({v2_left, v2_right, v2_up, v2_down}) +
                      (1 - m_alpha) * 0.25 * (v2_left + v2_right + v2_up + v2_down);

        // Ecriture dans le buffer
        size_t idx = ((i + 1) * m_stride + (j + 1)) * 2;
        
        #pragma omp atomic write
        m_buffer[idx] = val1;
        
        #pragma omp atomic write
        m_buffer[idx + 1] = val2;
    }

    void update() {
        m_data.swap(m_buffer);
        cl_update();
        get_val(m_pos_food.x, m_pos_food.y, 0) = 1.;
        get_val(m_pos_nest.x, m_pos_nest.y, 1) = 1.;
    }

    // ACCESSEURS POUR MPI
    std::vector<double>& get_data() { return m_data; }
    size_t get_total_size() const { return m_data.size(); }

private:
    void cl_update() {
        // Bords à -1
        for (unsigned long k = 0; k < m_stride; ++k) {
            // Haut et Bas
            size_t idx_top = k * 2;
            size_t idx_bot = (k + m_stride * (m_dim + 1)) * 2;
            m_data[idx_top] = -1.; m_data[idx_top+1] = -1.;
            m_data[idx_bot] = -1.; m_data[idx_bot+1] = -1.;
            
            // Gauche et Droite (k représente ici la ligne)
            size_t idx_left  = (k * m_stride) * 2;
            size_t idx_right = (k * m_stride + m_dim + 1) * 2;
            m_data[idx_left] = -1.; m_data[idx_left+1] = -1.;
            m_data[idx_right] = -1.; m_data[idx_right+1] = -1.;
        }
    }

    size_t index(const position_t& pos) const {
        return ((pos.x + 1) * m_stride + pos.y + 1) * 2;
    }

    unsigned long m_dim, m_stride;
    double m_alpha, m_beta;
    // Vecteur plat : [case0_v1, case0_v2, case1_v1, case1_v2, ...]
    std::vector<double> m_data; 
    std::vector<double> m_buffer;
    position_t m_pos_nest, m_pos_food;
};

#endif