#include "parser.h"

#include <algorithm>
#include <cctype>
#include <functional>
#include <iostream>
#include <sstream>
#include <string>

parser::parser(std::ifstream& handle) : handle(handle) {}
parser::~parser() {
    if(handle) handle.close();
}

std::ostream& operator<<(std::ostream& os, const vec3& rhs) {
    os << "(" << rhs.x << ", " << rhs.y << ", " << rhs.z << ")";
    return os;
}

soliton::soliton(double mass, 
                 vec3 pos,
                 vec3 vel, 
                 double phase
                ) : 
                mass(mass),
                pos(pos),
                vel(vel),
                phase(phase)
                {}

vec3::vec3(double x, double y, double z)
    : x(x), y(y), z(z) {}

settings::settings(int grid_size, 
                   double box_size,   
                   double max_time,
                   int threads,
                   int n_sol,
                   int n_pot,
                   bool ascii,
                   std::string soliton_file,
                   std::string wisdom_file,
                   std::string renderer,
                   std::string video_filename,
                   std::vector<soliton> solitons,
                   std::vector<potential> potentials) : 
                   grid_size(grid_size),
                   box_size(box_size),
                   max_time(max_time),
                   threads(threads),
                   n_sol(n_sol),
                   n_pot(n_pot),
                   ascii(ascii),
                   soliton_file(soliton_file),
                   wisdom_file(wisdom_file),
                   video_filename(video_ext),
                   renderer(renderer),
                   solitons(solitons),
                   potentials(potentials) {}

potential::potential(double mass, std::string type, vec3 center)
     : mass(mass), type(type), center(center), cmp_ptr(nullptr) {}

double potential::compute(vec3 pos) {
    return (cmp_ptr ? cmp_ptr(pos) : 0.0);
}

void potential::set_type(std::string type) {
    if(type == "inv_square") {
        cmp_ptr = std::bind(&potential::inv_sqr, this, std::placeholders::_1);
    } else if (type == "quadratic") {
        cmp_ptr = std::bind(&potential::quadratic, this, std::placeholders::_1);
    }

    this->type = type;
}

double potential::inv_sqr(vec3 pos) {
    std::cout << mass << std::endl;
    std::cout << pos << std::endl;
    double radius = DIST(pos.x - center.x, pos.y - center.y, pos.z - center.z);
    return -mass / radius;
}

double potential::quadratic(vec3 pos) {
    double radius = DIST2(pos.x - center.x, pos.y - center.y, pos.z - center.z);
    return 0.5 * radius;
}

std::string strip(const std::string& s, const std::string& chars=" ") {
    size_t begin = 0;
    size_t end = s.size()-1;
    for(; begin < s.size(); begin++)
        if(chars.find_first_of(s[begin]) == std::string::npos)
            break;
    for(; end > begin; end--)
        if(chars.find_first_of(s[end]) == std::string::npos)
            break;
    return s.substr(begin, end-begin+1);
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
                            input = strip(input, "\"");
                            std::transform(input.begin(), input.end(), input.begin(), ::tolower);
                            current_potential.set_type(input);
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
                        } else if (key == "threads") {
                            config.threads = std::stoi(input);
                        } else if (key == "solitons") {
                            config.n_sol = std::stoi(input);
                        } else if (key == "potentials") {
                            config.n_pot = std::stoi(input);
                        } else if (key == "ascii") {
                            std::string orig = input;
                            input = strip(input, "\"");
                            std::transform(input.begin(), input.end(), input.begin(), ::tolower);
                            config.ascii = input == "true" || input == "1";

                            if(!config.ascii && input != "false" && input != "0") {
                                std::cout << "[Warning] Input (" << orig << ") is not true/1, false/0, treated as false.\n";
                            }
                        } else if (key == "initialize_path") {
                            config.soliton_file = strip(input, "\"");
                        } else if (key == "wisdom_path") {
                            config.wisdom_file = strip(input, "\"");
                        } else if (key=="output_filename") {
                            config.video_ext = strip(input, "\"");
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