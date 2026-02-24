#include <vector>
#include <iostream>
#include <random>
#include <chrono>
#include <iomanip>
#include <mpi.h>
#include "fractal_land.hpp"
#include "ant_population.hpp"
#include "pheronome.hpp"
#include "renderer.hpp"
#include "window.hpp"
#include "rand_generator.hpp"

int main(int nargs, char* argv[])
{
    // ============================================================
    // 1. INITIALISATION MPI
    // ============================================================
    MPI_Init(&nargs, &argv);

    int rank, nb_procs;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nb_procs);

    // ============================================================
    // 2. PARAMÈTRES DE SIMULATION
    // ============================================================
    bool enable_gui = false;  // GUI désactivé en mode MPI
    std::size_t max_iterations = 500;

    std::size_t seed = 2026;
    
    const int total_nb_ants = 100000;

    const double eps = 0.8;
    const double alpha = 0.7;
    const double beta = 0.999;
    position_t pos_nest{256, 256};
    position_t pos_food{500, 500};

    // ============================================================
    // 3. DIVISION DES FOURMIS ENTRE LES PROCESSUS
    // ============================================================
    // Chaque processus gère un sous-ensemble de fourmis
    int local_nb_ants = total_nb_ants / nb_procs;
    int remainder = total_nb_ants % nb_procs;
    // Les premiers 'remainder' processus prennent une fourmi de plus
    if (rank < remainder) {
        local_nb_ants += 1;
    }

    // Graine unique par processus pour que les fourmis soient différentes
    std::size_t local_seed = seed + rank * 1000;

    if (rank == 0) {
        std::cout << "=== Configuration MPI ===" << std::endl;
        std::cout << "Nombre de processus : " << nb_procs << std::endl;
        std::cout << "Fourmis totales     : " << total_nb_ants << std::endl;
        std::cout << "Fourmis par proc    : ~" << (total_nb_ants / nb_procs) << std::endl;
        std::cout << "=========================" << std::endl;
    }

    // ============================================================
    // 4. CONSTRUCTION DU TERRAIN (identique sur tous les processus)
    // ============================================================
    // SDL_Init nécessaire même sans GUI pour les types SDL_Point
    SDL_Init(SDL_INIT_VIDEO);

    fractal_land land(8, 2, 1., 1024);

    // Normalisation du terrain
    double max_val = 0.0;
    double min_val = 0.0;
    for (fractal_land::dim_t i = 0; i < land.dimensions(); ++i)
        for (fractal_land::dim_t j = 0; j < land.dimensions(); ++j) {
            max_val = std::max(max_val, land(i, j));
            min_val = std::min(min_val, land(i, j));
        }
    double delta = max_val - min_val;
    for (fractal_land::dim_t i = 0; i < land.dimensions(); ++i)
        for (fractal_land::dim_t j = 0; j < land.dimensions(); ++j) {
            land(i, j) = (land(i, j) - min_val) / delta;
        }

    // ============================================================
    // 5. CRÉATION DE LA POPULATION LOCALE ET DES PHÉROMONES
    // ============================================================
    // Chaque processus crée SA population locale avec SA graine
    AntPopulation population(local_nb_ants, land, eps, local_seed);

    // Chaque processus possède une copie COMPLÈTE de la grille de phéromones
    pheronome phen(land.dimensions(), pos_food, pos_nest, alpha, beta);

    // Taille du vecteur de phéromones pour MPI_Allreduce
    size_t phen_total_size = phen.get_total_size();

    // Buffer temporaire pour MPI_Allreduce
    std::vector<double> global_phen_buffer(phen_total_size);

    // GUI uniquement sur le rank 0 (optionnel, désactivé par défaut en MPI)
    Window* win = nullptr;
    Renderer* renderer = nullptr;
    if (enable_gui && rank == 0) {
        win = new Window("Ant Simulation - MPI", 2 * land.dimensions() + 10, land.dimensions() + 266);
        renderer = new Renderer(land, phen, pos_nest, pos_food, population);
    }

    // ============================================================
    // 6. BOUCLE PRINCIPALE DE SIMULATION
    // ============================================================
    size_t local_food_quantity = 0;

    double total_time_move = 0.0;
    double total_time_evap = 0.0;
    double total_time_update = 0.0;
    double total_time_comm = 0.0;

    // Barrière pour synchroniser le démarrage
    MPI_Barrier(MPI_COMM_WORLD);
    auto start_global = std::chrono::high_resolution_clock::now();

    for (std::size_t it = 0; it < max_iterations; ++it) {

        // ---------------------------------------------------------
        // 6.1. MOUVEMENT DES FOURMIS LOCALES
        // ---------------------------------------------------------
        auto start_move = std::chrono::high_resolution_clock::now();

        population.advance_all(phen, land, pos_food, pos_nest, local_food_quantity);

        auto end_move = std::chrono::high_resolution_clock::now();
        total_time_move += std::chrono::duration<double>(end_move - start_move).count();

        // ---------------------------------------------------------
        // 6.2. COMMUNICATION MPI : Fusion des phéromones
        // ---------------------------------------------------------
        // Chaque processus a modifié sa grille locale.
        // On prend le MAX de toutes les grilles (comme dit dans le PDF).
        auto start_comm = std::chrono::high_resolution_clock::now();

        MPI_Allreduce(
            phen.get_data().data(),       // Envoi : données locales
            global_phen_buffer.data(),     // Réception : résultat global
            static_cast<int>(phen_total_size),
            MPI_DOUBLE,
            MPI_MAX,                       // On prend la valeur MAX
            MPI_COMM_WORLD
        );

        // Recopie du résultat global dans la grille locale
        phen.get_data() = global_phen_buffer;

        auto end_comm = std::chrono::high_resolution_clock::now();
        total_time_comm += std::chrono::duration<double>(end_comm - start_comm).count();

        // ---------------------------------------------------------
        // 6.3. ÉVAPORATION (chaque processus fait la même chose)
        // ---------------------------------------------------------
        auto start_evap = std::chrono::high_resolution_clock::now();
        phen.do_evaporation();
        auto end_evap = std::chrono::high_resolution_clock::now();
        total_time_evap += std::chrono::duration<double>(end_evap - start_evap).count();

        // ---------------------------------------------------------
        // 6.4. MISE À JOUR DES PHÉROMONES
        // ---------------------------------------------------------
        auto start_update = std::chrono::high_resolution_clock::now();
        phen.update();
        auto end_update = std::chrono::high_resolution_clock::now();
        total_time_update += std::chrono::duration<double>(end_update - start_update).count();

        // GUI (rank 0 uniquement)
        if (enable_gui && rank == 0) {
            SDL_Event event;
            while (SDL_PollEvent(&event)) {
                if (event.type == SDL_QUIT) max_iterations = 0;
            }
            renderer->display(*win, local_food_quantity);
            win->blit();
        }
    }

    auto end_global = std::chrono::high_resolution_clock::now();
    double total_simulation_time = std::chrono::duration<double>(end_global - start_global).count();

    // ============================================================
    // 7. RÉDUCTION FINALE : Somme de la nourriture collectée
    // ============================================================
    size_t global_food_quantity = 0;
    MPI_Reduce(
        &local_food_quantity,
        &global_food_quantity,
        1,
        MPI_UNSIGNED_LONG,
        MPI_SUM,
        0,                    // Le rank 0 reçoit le résultat
        MPI_COMM_WORLD
    );

    // ============================================================
    // 8. AFFICHAGE DES RÉSULTATS (rank 0 uniquement)
    // ============================================================
    if (rank == 0) {
        if (enable_gui) {
            delete renderer;
            delete win;
        }
        SDL_Quit();

        std::cout << std::endl;
        std::cout << "Simulation terminée après " << max_iterations << " itérations." << std::endl;
        std::cout << "Nourriture totale collectée : " << global_food_quantity << std::endl;
        std::cout << "------------------------------------------------------------" << std::endl;
        std::cout << "TEMPS TOTAL (MPI, " << nb_procs << " processus) : " 
                  << total_simulation_time << " s" << std::endl;
        std::cout << "Temps moyen par itération : " 
                  << (total_simulation_time / max_iterations) << " s" << std::endl;
        std::cout << "------------------------------------------------------------" << std::endl;
        std::cout << "DETAIL DES COMPOSANTS (Cumulé sur rank 0) :" << std::endl;
        std::cout << "  - Mouvement fourmis : " << total_time_move << " s (" 
                  << (total_time_move / total_simulation_time) * 100 << "%)" << std::endl;
        std::cout << "  - Communication MPI : " << total_time_comm << " s (" 
                  << (total_time_comm / total_simulation_time) * 100 << "%)" << std::endl;
        std::cout << "  - Evaporation       : " << total_time_evap << " s (" 
                  << (total_time_evap / total_simulation_time) * 100 << "%)" << std::endl;
        std::cout << "  - Update Pheromones : " << total_time_update << " s (" 
                  << (total_time_update / total_simulation_time) * 100 << "%)" << std::endl;
        std::cout << "------------------------------------------------------------" << std::endl;
    } else {
        SDL_Quit();
    }

    // ============================================================
    // 9. FINALISATION MPI
    // ============================================================
    MPI_Finalize();

    return 0;
}