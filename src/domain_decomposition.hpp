#ifndef _DOMAIN_DECOMPOSITION_HPP_
#define _DOMAIN_DECOMPOSITION_HPP_

#include <mpi.h>
#include <vector>
#include <cstddef>
#include <algorithm>
#include <iostream>

/**
 * @brief Gère la décomposition 1D (par lignes) de la grille de phéromones.
 * 
 * Chaque processus possède un sous-ensemble de lignes de la grille,
 * plus une rangée de "ghost cells" en haut et en bas pour les voisins.
 * 
 * Grille globale de dimension `dim` lignes :
 *   Rank 0 : lignes [0, row_end_0)           + ghost bottom
 *   Rank 1 : lignes [row_start_1, row_end_1) + ghost top + ghost bottom
 *   ...
 *   Rank P-1 : lignes [row_start_P-1, dim)   + ghost top
 */
class DomainDecomposition {
public:
    int rank;
    int nb_procs;
    
    // Dimensions globales
    unsigned long global_dim;       // Taille de la grille (ex: 1025)
    
    // Dimensions locales (sans ghost cells)
    unsigned long local_row_start;  // Première ligne locale (dans coords globales)
    unsigned long local_row_end;    // Dernière ligne + 1 (dans coords globales)
    unsigned long local_nrows;      // Nombre de lignes locales
    
    // Voisins MPI
    int rank_top;    // Rang du voisin du haut (-1 si bord)
    int rank_bottom; // Rang du voisin du bas (-1 si bord)
    
    // Stride pour le stockage local (avec ghost cells)
    // local_nrows + 2 ghost rows (top + bottom)
    unsigned long local_stride_rows;
    unsigned long stride_cols;  // = global_dim + 2 (bordures phéromones)

    DomainDecomposition(unsigned long dim, int rank, int nb_procs)
        : rank(rank), nb_procs(nb_procs), global_dim(dim)
    {
        // Découpage équitable des lignes
        unsigned long rows_per_proc = dim / nb_procs;
        unsigned long remainder = dim % nb_procs;
        
        if (static_cast<unsigned long>(rank) < remainder) {
            local_nrows = rows_per_proc + 1;
            local_row_start = rank * local_nrows;
        } else {
            local_nrows = rows_per_proc;
            local_row_start = remainder * (rows_per_proc + 1) + 
                              (rank - remainder) * rows_per_proc;
        }
        local_row_end = local_row_start + local_nrows;
        
        // Voisins
        rank_top    = (rank > 0) ? rank - 1 : MPI_PROC_NULL;
        rank_bottom = (rank < nb_procs - 1) ? rank + 1 : MPI_PROC_NULL;
        
        // Stride avec ghost cells
        local_stride_rows = local_nrows + 2; // +1 ghost top, +1 ghost bottom
        stride_cols = global_dim + 2;        // +1 ghost left, +1 ghost right
        
        if (rank == 0) {
            std::cout << "=== Décomposition de domaine ===" << std::endl;
            std::cout << "Grille globale : " << dim << " x " << dim << std::endl;
        }
        MPI_Barrier(MPI_COMM_WORLD);
        std::cout << "  Rank " << rank << " : lignes [" << local_row_start 
                  << ", " << local_row_end << ") = " << local_nrows 
                  << " lignes, voisins: top=" << rank_top 
                  << " bottom=" << rank_bottom << std::endl;
    }
    
    /**
     * @brief Convertit une coordonnée globale i en indice local (avec offset ghost)
     */
    unsigned long global_to_local_i(unsigned long global_i) const {
        // +1 pour la ghost cell du haut
        return global_i - local_row_start + 1;
    }
    
    /**
     * @brief Vérifie si une position globale (i,j) est dans le domaine local
     */
    bool is_local(int global_i, int global_j) const {
        return global_i >= static_cast<int>(local_row_start) && 
               global_i < static_cast<int>(local_row_end) &&
               global_j >= 0 && 
               global_j < static_cast<int>(global_dim);
    }
    
    /**
     * @brief Taille totale du buffer local de phéromones (avec ghost cells)
     *        2 valeurs par cellule (V1, V2)
     */
    size_t local_phen_size() const {
        return local_stride_rows * stride_cols * 2;
    }

    /**
     * @brief Échange les ghost cells avec les voisins
     * @param local_data Le buffer local de phéromones
     * 
     * Convention de stockage : 
     *   Index dans local_data = ((local_i) * stride_cols + (j+1)) * 2 + type
     *   local_i = 0 : ghost top
     *   local_i = 1..local_nrows : données réelles
     *   local_i = local_nrows+1 : ghost bottom
     */
    void exchange_ghost_cells(std::vector<double>& local_data) const {
        // Nombre de doubles par ligne = stride_cols * 2
        int line_size = static_cast<int>(stride_cols * 2);
        
        MPI_Status status;
        
        // --- Envoi vers le bas, réception du haut ---
        // On envoie notre dernière ligne réelle (local_i = local_nrows) au voisin du bas
        // On reçoit dans notre ghost top (local_i = 0) depuis le voisin du haut
        MPI_Sendrecv(
            &local_data[local_nrows * stride_cols * 2],  // Envoi : dernière ligne réelle
            line_size, MPI_DOUBLE, rank_bottom, 0,
            &local_data[0],                               // Réception : ghost top
            line_size, MPI_DOUBLE, rank_top, 0,
            MPI_COMM_WORLD, &status
        );
        
        // --- Envoi vers le haut, réception du bas ---
        // On envoie notre première ligne réelle (local_i = 1) au voisin du haut
        // On reçoit dans notre ghost bottom (local_i = local_nrows+1) depuis le voisin du bas
        MPI_Sendrecv(
            &local_data[1 * stride_cols * 2],              // Envoi : première ligne réelle
            line_size, MPI_DOUBLE, rank_top, 1,
            &local_data[(local_nrows + 1) * stride_cols * 2], // Réception : ghost bottom
            line_size, MPI_DOUBLE, rank_bottom, 1,
            MPI_COMM_WORLD, &status
        );
    }
};

#endif