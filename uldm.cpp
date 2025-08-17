#include <algorithm>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <chrono>

#include "parser.hpp"
#include "solver.hpp"

#include <fftw3.h>
#include <omp.h>

#ifdef FFMPEG
    #include <thread>
    #include "ffmpeg.hpp"
#endif 

int main(int argc, char** argv) {
    if(argc < 2) { 
        std::cerr << "Usage: <config_file:string>";
        return -1;
    }

    std::ifstream file(argv[1]);
    if(!file) { std::cout << "Error opening file!" << std::endl; return -1; }

    parser reader = parser(file);
    settings config = reader.parse();

    std::ofstream en_file("energy.csv");

    double L = config.box_size;
    int N = config.grid_size;
    int NNN = CUBE(N);
    int nthreads = std::min((unsigned int)config.threads, (unsigned int)omp_get_max_threads());

    double vmin =  1000000000.0;
    double vmax = -1000000000.0;

    poisson_schrodinger_state* grid = new poisson_schrodinger_state[NNN];
    double* fourier_grid = new double[NNN];

    fftw_init_threads();
    fftw_plan_with_nthreads(nthreads);

    auto t0 = std::chrono::high_resolution_clock::now();

    en_file << "# t, DE/E_0, E_tot, E_k + E_q, E_self, E_pot\n";
    en_file << 0.0 << ", ";

    generate_fourier_grid(fourier_grid, N, L);
    fftw_import_wisdom_from_filename(config.wisdom_file.c_str());

    auto t_frame_start = std::chrono::high_resolution_clock::now();
    double energy = init_grid(grid, fourier_grid, config, vmin, vmax, en_file);
    auto t_frame_end = std::chrono::high_resolution_clock::now();

    double frame_time = std::chrono::duration<double>(t_frame_end - t_frame_start).count();
    double dx = (double)(L / N);
    double dt = dx * dx / M_PI;
    double max_time = config.max_time;

    std::cout << "[INFO] Simulation settings: [" << N << "x" << N << ", L: " << L 
    << ", max_time: " << max_time << ", threads: " << nthreads << "]\n";

    #ifdef FFMPEG
        std::thread ffmpeg_thread(start_ffmpeg, config);
    #endif

    int nframes = (int)(max_time / dt);

    double mean_frame_time = frame_time;
    int counter = 1;

    for(double t = 0.0; t < max_time; t += dt) {
        t_frame_start = std::chrono::high_resolution_clock::now();
        pseudo_spectral_solver(grid, fourier_grid, config, vmin, vmax, en_file, energy);
        t_frame_end = std::chrono::high_resolution_clock::now();

        frame_time = std::chrono::duration<double>(t_frame_end - t_frame_start).count();
        mean_frame_time *= N++;
        mean_frame_time += frame_time;
        mean_frame_time /= N;

        double eta_s = ((nframes - (int)(t / dt) - 1) * mean_frame_time);
        int eta_m = (int)(eta_s / 60.0);
        eta_s -= eta_m * 60.0;

        std::cout << "\r[INFO] Frames: " 
        << std::setw(4) << std::setfill('0') << (int)(t / dt) + 1
        << "/" << std::setw(4) << std::setfill('0') << nframes 
        << " (" << std::setprecision(4) << frame_time << "s/" << mean_frame_time << "s)" <<
        ". ETA: " << std::setw(2) << std::setfill('0') << eta_m << " min " << std::setw(2) << std::setfill('0') << (int)eta_s << " s   " << std::flush;
        en_file << t << ", ";
    }

    #ifdef FFMPEG
        end_ffmpeg();
    #endif

    auto t1 = std::chrono::high_resolution_clock::now();
    
    double elapsed_s = std::chrono::duration<double>(t1 - t0).count();
    int elapsed_m = (int)(elapsed_s / 60.0);
    elapsed_s -= elapsed_m * 60.0;
    std::cout << "\n[INFO] Execution time: " << std::setw(2) << std::setfill('0') <<
    elapsed_m << " min " << std::setw(2) << std::setfill('0') << (int)elapsed_s << " s\n";

    fftw_export_wisdom_to_filename(config.wisdom_file.c_str());
    fftw_cleanup_threads();
    fftw_cleanup();

    file.close();
    en_file.close();

    delete[] grid;
    delete[] fourier_grid;

    return 0;
}

