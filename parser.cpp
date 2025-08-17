#include "solver.hpp"
#include "parser.hpp"

#include <algorithm>
#include <cctype>
#include <iostream>
#include <sstream>
#include <string>

parser::parser(std::ifstream& handle) : handle(handle) {}
parser::~parser() {
    if(handle) handle.close();
}

settings parser::parse() {
    settings config = {};
    if(!handle) { return config; }

    std::string input;
    soliton current_soliton = {};
    potential current_potential = {};

    bool soliton_flag = false;
    bool potential_flag = false;

    while(std::getline(handle, input)) {
        std::istringstream ss(input);
        std::string key = "";

        if(input.find("---") != std::string::npos) {
            std::transform(input.begin(), input.end(), input.begin(), ::tolower);
            bool soliton_header = input.find("soliton") != std::string::npos;
            bool potential_header = input.find("potential") != std::string::npos;
            bool end_header = input.find("end") != std::string::npos;
            if(soliton_header) {
                if(soliton_flag) { 
                    config.solitons.push_back(current_soliton); 
                    current_soliton = {};
                }
                soliton_flag = true;
            } else if (potential_header) {
                if(soliton_flag) {
                    config.solitons.push_back(current_soliton); 
                    current_soliton = {};
                }
                if (potential_flag) {
                    config.potentials.push_back(current_potential);
                    current_potential = {};
                }
                potential_flag = true;
            } else if (end_header) {
                if(soliton_flag) {
                    config.solitons.push_back(current_soliton); 
                    current_soliton = {};
                }
                if (potential_flag) {
                    config.potentials.push_back(current_potential);
                    current_potential = {};
                }
                soliton_flag = false;
                potential_flag = false;
            } else {
                input = strip(input, "---");
                std::cout << "[Warning] Unknown header section found (" << input << ").\n Parser will ignore and continue.\n";
            }
            } else {
                while(std::getline(ss, input, ':')) {
                input = strip(input);
                if(key == "") { 
                    std::transform(input.begin(), input.end(), input.begin(), ::tolower);
                    key = input;
                } else {
                    if(soliton_flag) {
                        if (key == "mass") {
                            current_soliton.mass = std::stod(input);
                        } else if(key == "pos" || key == "position") {
                            current_soliton.pos = parse_vec3(input);
                        } else if(key == "vel" || key == "velocity") {
                            current_soliton.vel = parse_vec3(input);
                        } else if(key == "phase") {
                            current_soliton.phase = std::stod(input);
                        } else {
                            std::cout << "[Warning] Unknown soliton parameter has been encountered (" << key << ") with value: \"" << input;
                            std::cout << "\". Parser will ignore and continue.\n";
                        }
                    } else if(potential_flag) {
                        if (key == "mass") {
                            current_potential.mass = std::stod(input);
                        } else if (key == "type") {
                            std::string orig = input;
                            input = strip(input, "\"");
                            std::transform(input.begin(), input.end(), input.begin(), ::tolower);
                            current_potential.type = input;
                            if(input != "central" && input != "quadratic") {
                                std::cout << "[WARNING] The specified potential type (" << orig << ") is not supported. Solver will ignore it.\n";
                            }
                        } else if (key == "center" || key == "pos" || key == "position") {
                            current_potential.center = parse_vec3(input);
                        } else {
                            std::cout << "[Warning] Unknown potential parameter has been encountered (" << key << ") with value: \"" << input;
                            std::cout << "\". Parser will ignore and continue.\n";
                        }
                    } else {
                        if(key == "grid_size") {
                            config.grid_size = std::stoi(input);
                        } else if (key == "box_size") {
                            config.box_size = std::stod(input);
                        } else if (key == "max_time") {
                            config.max_time = std::stod(input);
                        } else if (key == "vmin") {
                            config.vmin = std::stod(input);
                        } else if (key == "vmax") {
                            config.vmax = std::stod(input);
                        } else if (key == "threads") {
                            config.threads = std::stoi(input);
                        } else if (key == "framerate") {
                            config.framerate = std::stoi(input);
                        } else if (key == "solitons") {
                            config.n_sol = std::stoi(input);
                        } else if (key == "potentials") {
                            config.n_pot = std::stoi(input);
                        } else if (key == "ascii") {
                            config.ascii = parse_bool(input);
                        } else if (key == "ffmpeg_cmd") {
                            config.ffmpeg_cmd = strip(input, "\"");
                        } else if (key == "initialize_path") {
                            config.soliton_file = strip(input, "\"");
                        } else if (key == "wisdom_path") {
                            config.wisdom_file = strip(input, "\"");
                        } else if (key == "output_filepath" || key == "video_filepath") {
                            config.video_filename = strip(input, "\"");
                        } else if (key == "renderer") {
                            input = strip(input, "\"");
                            std::transform(input.begin(), input.end(), input.begin(), ::tolower);
                            config.renderer = input;
                        } else {
                            std::cout << "[Warning] Unknown key has been encountered (" << key << ") with value: \"" << input;
                            std::cout << "\". Parser will ignore and continue.\n";
                        }
                    }
                }
            }
        }
    }

    if(soliton_flag) config.solitons.push_back(current_soliton); 
    if (potential_flag) config.potentials.push_back(current_potential);

    if(config.n_sol != config.solitons.size()) {
        std::cout << "[Warning] The number of input solitons (" << config.n_sol << ") doesn't match ";
        std::cout << "with the actual solitons provided (" << config.solitons.size() << ").\n";
        std::cout << "If you believe this is a mistake, double check your config file.\n";

        config.n_sol = config.solitons.size();
    }

    if(config.n_pot != config.potentials.size()) {
        std::cout << "[Warning] The number of input potentials (" << config.n_pot << ") doesn't match ";
        std::cout << "with the actual potnentials provided (" << config.potentials.size() << ").\n";
        std::cout << "If you believe this is a mistake, double check your config file.\n";

        config.n_pot = config.potentials.size();
    }

    return config; 
}

vec3 parser::parse_vec3(std::string input) {
    input = strip(input);
    unsigned int back_idx = input.size() - 1;
    unsigned int left = input.find_first_of('(');
    unsigned int right = input.find_last_of(')');

    input = input.substr(left + 1, right - 1);

    std::istringstream vec_input(input);
    vec3 data = {};

    int idx = 0;
    while(std::getline(vec_input, input, ',')) {
        data.data[idx++] = std::stod(input);
    }

    if(left != 0 || right != back_idx) {
        std::cout << "[Warning] Check vec3 input for undesired mispell! Retrieved value: ";
        std::cout << data << std::endl;
    }

    return data;
}

bool parser::parse_bool(std::string input) {
    std::string orig = input;
    input = strip(input, "\"");
    std::transform(input.begin(), input.end(), input.begin(), ::tolower);
    bool result = input == "true" || input == "1";

    if(!result && input != "false" && input != "0") {
        std::cout << "[Warning] Input (" << orig << ") is not true/1, false/0, treated as false.\n";
    }

    return result;
}