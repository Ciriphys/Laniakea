#include "parser.h"

#include <iostream>
#include <fstream>

int main(int argc, char** argv) {
    std::ifstream config_file("mock.cfg");
    parser reader = parser(config_file);

    settings config = reader.parse();

    std::cout << "Retrieved configuration: " << std::endl;
    std::cout << " - grid_size        : " << config.grid_size     << std::endl;
    std::cout << " - box_size         : " << config.box_size      << std::endl;
    std::cout << " - max_time         : " << config.max_time      << std::endl;
    std::cout << " - threads          : " << config.threads       << std::endl;
    std::cout << " - solitons         : " << config.n_sol         << std::endl;
    std::cout << " - ascii            : " << config.ascii         << std::endl;
    std::cout << " - initialize_path  : " << config.soliton_file  << std::endl;
    std::cout << " - wisdom_path      : " << config.wisdom_file   << std::endl;
    std::cout << " - renderer         : " << config.renderer      << std::endl;
    std::cout << " - video_ext        : " << config.video_ext     << std::endl;

    for (soliton sol : config.solitons) {
        std::cout << "\n--- SOLITON ---\n";
        std::cout << " - mass     : " << sol.mass  << std::endl;
        std::cout << " - position : " << sol.pos   << std::endl;
        std::cout << " - velocity : " << sol.vel   << std::endl;
        std::cout << " - phase    : " << sol.phase << std::endl;
        std::cout << "--- END ---\n";
    }

    for (potential pot : config.potentials) {
        std::cout << "\n--- POTENTIAL ---\n";
        std::cout << " - type   : " << pot.type   << std::endl;
        std::cout << " - mass   : " << pot.mass   << std::endl;
        std::cout << " - center : " << pot.center << std::endl;
        std::cout << "--- END ---\n";
    }

    return 0;
}
