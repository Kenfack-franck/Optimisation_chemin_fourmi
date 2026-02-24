#include <vector>
#include <iostream>
#include <chrono>
#include <mpi.h>
#include "fractal_land.hpp"
#include "ant_population_local.hpp"
#include "pheronome_local.hpp"
#include "domain_decomposition.hpp"
#include "rand_generator.hpp"
#include "basic_types.hpp"

int main(int nargs, char* argv[])
{
    MPI_Init(&nargs, &argv);

    int rank, nb_procs;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nb_procs);

    // ============================================================
    // PARAMÈTRES
    // ============================================================
    std::size_t max_iterations = 500;
    std::size_t seed = 2026;
    const int total_nb_ants = 100000; // Modifiable pour les tests
    const double eps = 0.8;
    const double alpha = 0.7;
    const double beta = 0.999;
    position_t pos_nest{256, 256};
    position_t pos_food{500, 500};

    // SDL_Init nécessaire pour SDL_Point
    SDL_Init(SDL_INIT_VIDEO);

    // ============================================================
    // CONSTRUCTION DU TERRAIN (identique sur tous les processus)
    // ============================================================
    fractal_land land(8, 2, 1., 1024);

    double max_val = 0.0, min_val = 0.0;
    for (fractal_land::dim_t i = 0; i < land.dimensions(); ++i)
        for (fractal_land::dim_t j = 0; j < land.dimensions(); ++j) {
            max_val = std::max(max_val, land(i, j));
            min_val = std::min(min_val, land(i, j));
        }
    double delta = max_val - min_val;
    for (fractal_land::dim_t i = 0; i < land.dimensions(); ++i)
        for (fractal_land::dim_t j = 0; j < land.dimensions(); ++j)
            land(i, j) = (land(i, j) - min_val) / delta;

    // ============================================================
    // DÉCOMPOSITION DE DOMAINE
    // ============================================================
    DomainDecomposition domain(land.dimensions(), rank, nb_procs);

    // ============================================================
    // PHÉROMONES LOCALES + POPULATION LOCALE
    // ============================================================
    PheronomeLocal phen(domain, pos_food, pos_nest, alpha, beta);
    AntPopulationLocal population(total_nb_ants, land, eps, seed, domain);

    if (rank == 0) {
        std::cout << "=== Configuration MPI (Stratégie 2 - Domaine) ===" << std::endl;
        std::cout << "Nombre de processus : " << nb_procs << std::endl;
        std::cout << "Fourmis totales     : " << total_nb_ants << std::endl;
        std::cout << "Grille              : " << land.dimensions() << " x " << land.dimensions() << std::endl;
        std::cout << "Mémoire phéromones/proc : " 
                  << (phen.get_total_size() * 8) / 1024 << " Ko (vs " 
                  << (2 * (land.dimensions()+2) * (land.dimensions()+2) * 8) / 1024 / 1024 
                  << " Mo global)" << std::endl;
        std::cout << "=================================================" << std::endl;
    }

    // ============================================================
    // BOUCLE PRINCIPALE
    // ============================================================
    size_t local_food_quantity = 0;
    double total_time_move = 0.0;
    double total_time_comm = 0.0;
    double total_time_evap = 0.0;
    double total_time_update = 0.0;
    double total_time_migrate = 0.0;

    MPI_Barrier(MPI_COMM_WORLD);
    auto start_global = std::chrono::high_resolution_clock::now();

    for (std::size_t it = 0; it < max_iterations; ++it) {

        // 1. ÉCHANGE DES GHOST CELLS (phéromones des voisins)
        auto start_comm = std::chrono::high_resolution_clock::now();
        phen.exchange_ghosts();
        auto end_comm = std::chrono::high_resolution_clock::now();
        total_time_comm += std::chrono::duration<double>(end_comm - start_comm).count();

        // 2. MOUVEMENT DES FOURMIS LOCALES
        auto start_move = std::chrono::high_resolution_clock::now();
        population.advance_all(phen, land, pos_food, pos_nest, local_food_quantity);
        auto end_move = std::chrono::high_resolution_clock::now();
        total_time_move += std::chrono::duration<double>(end_move - start_move).count();

        // 3. MIGRATION DES FOURMIS HORS DOMAINE
        auto start_migrate = std::chrono::high_resolution_clock::now();
        population.migrate_ants();
        auto end_migrate = std::chrono::high_resolution_clock::now();
        total_time_migrate += std::chrono::duration<double>(end_migrate - start_migrate).count();

        // 4. ÉVAPORATION
        auto start_evap = std::chrono::high_resolution_clock::now();
        phen.do_evaporation();
        auto end_evap = std::chrono::high_resolution_clock::now();
        total_time_evap += std::chrono::duration<double>(end_evap - start_evap).count();

        // 5. UPDATE
        auto start_update = std::chrono::high_resolution_clock::now();
        phen.update();
        auto end_update = std::chrono::high_resolution_clock::now();
        total_time_update += std::chrono::duration<double>(end_update - start_update).count();
    }

    auto end_global = std::chrono::high_resolution_clock::now();
    double total_simulation_time = std::chrono::duration<double>(end_global - start_global).count();

    // ============================================================
    // RÉDUCTIONS FINALES
    // ============================================================
    size_t global_food_quantity = 0;
    MPI_Reduce(&local_food_quantity, &global_food_quantity, 1,
               MPI_UNSIGNED_LONG, MPI_SUM, 0, MPI_COMM_WORLD);

    // Nombre total de fourmis restantes (vérification)
    size_t local_ant_count = population.size();
    size_t global_ant_count = 0;
    MPI_Reduce(&local_ant_count, &global_ant_count, 1,
               MPI_UNSIGNED_LONG, MPI_SUM, 0, MPI_COMM_WORLD);

    // ============================================================
    // AFFICHAGE (rank 0)
    // ============================================================
    if (rank == 0) {
        SDL_Quit();

        std::cout << std::endl;
        std::cout << "Simulation terminée après " << max_iterations << " itérations." << std::endl;
        std::cout << "Nourriture totale collectée : " << global_food_quantity << std::endl;
        std::cout << "Fourmis totales restantes   : " << global_ant_count 
                  << " (attendu: " << total_nb_ants << ")" << std::endl;
        std::cout << "------------------------------------------------------------" << std::endl;
        std::cout << "TEMPS TOTAL (Stratégie 2, " << nb_procs << " processus) : " 
                  << total_simulation_time << " s" << std::endl;
        std::cout << "Temps moyen par itération : " 
                  << (total_simulation_time / max_iterations) << " s" << std::endl;
        std::cout << "------------------------------------------------------------" << std::endl;
        std::cout << "DETAIL DES COMPOSANTS (Cumulé sur rank 0) :" << std::endl;
        std::cout << "  - Mouvement fourmis  : " << total_time_move << " s (" 
                  << (total_time_move / total_simulation_time) * 100 << "%)" << std::endl;
        std::cout << "  - Ghost cells (comm) : " << total_time_comm << " s (" 
                  << (total_time_comm / total_simulation_time) * 100 << "%)" << std::endl;
        std::cout << "  - Migration fourmis  : " << total_time_migrate << " s (" 
                  << (total_time_migrate / total_simulation_time) * 100 << "%)" << std::endl;
        std::cout << "  - Evaporation        : " << total_time_evap << " s (" 
                  << (total_time_evap / total_simulation_time) * 100 << "%)" << std::endl;
        std::cout << "  - Update Pheromones  : " << total_time_update << " s (" 
                  << (total_time_update / total_simulation_time) * 100 << "%)" << std::endl;
        std::cout << "------------------------------------------------------------" << std::endl;
    } else {
        SDL_Quit();
    }

    MPI_Finalize();
    return 0;
}