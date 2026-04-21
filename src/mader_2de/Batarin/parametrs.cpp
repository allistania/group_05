#include "parameters.h"
#include <fstream>
#include <iostream>
#include <string>
#include <map>
#include <vector>
#include <cmath>

// ooieoey aey nicaaiey oanoa
void writeTest() {
    std::string test;
    std::cout << "Iacaaiea oanoa ";
    std::cin >> test;
    std::ofstream file(test);
    std::string rho_L, u_L, p_L, rho_R, u_R, p_R, t_end;
    std::cout << "rho_L= ";
    std::cin >> rho_L;
    std::cout << "u_L= ";
    std::cin >> u_L;
    std::cout << "p_L= ";
    std::cin >> p_L;
    std::cout << "rho_R= ";
    std::cin >> rho_R;
    std::cout << "u_R= ";
    std::cin >> u_R;
    std::cout << "p_R= ";
    std::cin >> p_R;
    std::cout << "t_end= ";
    std::cin >> t_end;

    if (file.is_open()) {
        file << "rho_L= " << rho_L << "\n"; 
        file << "u_L= " << u_L << "\n";
        file << "p_L= " << p_L << "\n";
        file << "rho_R= " << rho_R << "\n"; 
        file << "u_R= " << u_R << "\n";
        file << "p_R= " << p_R << "\n";
        file << "t_end= " << t_end << "\n";
        file << "\n";
        file << "?oeiaianoai e oanoo:\n";
        file << "A oanoao niaa??eony 7 ia?aiao?ia aey oanoia\n";
        file << "1: rho_L e rho_R - yoi ia?aeuiua ieioiinoe aaca neaaa e ni?aaa ia aaneiia?iinoe\n";
        file << "2: p_L e p_R - yoi ia?aeuiua aaaeaiey aaca neaaa e ni?aaa ia aaneiia?iinoe\n";
        file << "3: u_L e u_R - yoi ia?aeuiua nei?inoe aaca neaaa e ni?aaa ia aaneiia?iinoe\n"; 
        file << "4: t_end - a?aiy ?an??oa";

        file.close();
        std::cout << "Aaiiua caienaiu: " << test << std::endl;
    } else {
        std::cerr << "Ioeaea ia oaaeinu caienaou a oaee " << test << std::endl;
    }
}

//eieoeaeecaoey ia?aiao?ia 
void start() {
    std::string start;
    std::cout << "Name of init file: ";
    std::cin >> start;
    std::ofstream file(start);
    
    std::string tau, h, cells, cfl,solver_name, left_bc, right_bc, test_name, count, t_output, key;
    std::vector<std::string> available_solvers = {
        "godunov", "godunov-kolgan", "godunov-kolgan-rodionov", "eno", "weno" , "Fletcher"
    };
    
    std::cout << "Time step: ";
    std::cin >> tau;
    std::cout << "Coordinate step: ";
    std::cin >> h;
    std::cout << "Nubmer of cells: ";
    std::cin >> cells;
    std::cout << "CFK: ";
    std::cin >> cfl;
    std::cout << "Left boundary condition (wall/flow/periodic): ";
    std::cin >> left_bc;
    std::cout << "Right boundary condition (wall/flow/periodic): ";
    std::cin >> right_bc;
    std::cout << "Sod test (sod1/sod2/sod3/sod4/sod5/custom): ";
    std::cin >> test_name;
    std::string t_end_user;
    std::cout << "End time (enter 'default' or specific value): ";
    std::cin >> t_end_user; 
    std::cout << "count: ";
    std::cin >> count;
    std::cout << "Output interval for gif (time between frames): ";
    std::cin >> t_output;
    std::cout << "true - if time steps and talse if count steps: ";
    std::cin >> key;
    std::cout << "\nAvailable solvers: ";
    for (size_t i = 0; i< available_solvers.size(); i++) {
        std::cout << available_solvers[i];
        if (i< available_solvers.size() - 1) std::cout << ", ";
    }
    std::cout << std::endl;

    int num_solvers;
    std::cout << "Numbers of solvers";
    std::cin >> num_solvers;
    
    std::vector<std::string> selected_solvers;
    for (int i = 0; i< num_solvers; i++){
        std::string solver;
        std::cout << "Enter solver number " << i+1<< ": ";
        std::cin >> solver; 
        selected_solvers.push_back(solver);
    }

    std::string solver_names_str; 
    for (size_t i = 0; i<selected_solvers.size(); i++){
        solver_names_str +=selected_solvers[i];
        if (i< selected_solvers.size() -1) solver_names_str+= ",";
    }

    if (file.is_open()) {
        file << "tau = " << tau << "\n"; 
        file << "h = " << h << "\n"; 
        file << "cells = " << cells << "\n";
        file << "cfl = " << cfl << "\n";
        file << "left_bc = " << left_bc << "\n";
        file << "right_bc = " << right_bc << "\n";
        file << "test_name = " << test_name << "\n";
        file << "t_end = " << t_end_user << "\n";
        file << "count = " << count << "\n";
        file << "t_output = " << t_output << "\n";
        file << "key = " << key << "\n";
        file << "solver_names = " << solver_names_str << "\n";
        file << "\n";
        file << "?oeiaianoai e oaeeo eieoeaeecaoee:\n";
        file << "A oaeea niaa??eony 7 ia?aiao?ia aey eieoeaeecaoee\n";
        file << "1: tau e h - yoi ia?aeuiua oaae ii a?aiaie niioaaonoaaiii\n";
        file << "2: cells - yoi eiee?anoai ?an??oiuo y?aae\n";
        file << "3: CFL - ?enei Eo?aioa — O?ea?eona — Eaae\n";
        file << "4: left_dc e right_bc - yoi, niioaaonoaaiii, eaaia e i?aaia a?aie?iia oneiaea. ?aaeeciaaiu o?ao aeaia\n"; 
        file << "wall - io?a?aiea io noaiee \n";
        file << "flow - iioie \n";
        file << "periodic - ia?eiae?aneia AO \n";
        file << "5: test_name - eiy oanoa. Anee io?ai aioiaue oano Niaa, oi ieoai sod{iiia? oanoa}. Iai?eia?, sod1\n";
        file << "6: ?aciinoiay noaia: Aiaoiiaa, Aiaoiiaa-Eieaaia, Aiaoiiaa-Eieaaia-?iaeiiiaa, AII eee AAiI \n";
        file << "7: ia?aiao? count ioaa?aao ca caienu oaeea aey aeoee. Iai?eia?, anee count=100, oi ea?ao? nioo? iia?aoe? aoaao caienu oaeea \n";
        file << "8: a ioee?ea io ia?aiao?a count, t_output ioaa?aao ca caienu oaeea aey aeoee ia ?a?ac count-ue oaa, a ?a?ac t_output-ia a?aiy\n";
        file << "9: True anee ii a?aiaie caienu e False anee ii oaaai:\n";
        file.close();
        std::cout << "Aaiiua caienaiu: " << start << std::endl;
    } else {
        std::cerr << "Ioeaea ia oaaeinu caienaou a oaee " << start << std::endl;
    }
}

bool readTest(const std::string& filename, double& rho_L, double& u_L, double& p_L,
              double& rho_R, double& u_R, double& p_R, double& t_end) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: failed opening of file " << filename << std::endl;
        return false;
    }

    std::string line;
    std::map<std::string, double*> paramMap = {
        {"rho_L", &rho_L},
        {"u_L", &u_L},
        {"p_L", &p_L},
        {"rho_R", &rho_R},
        {"u_R", &u_R},
        {"p_R", &p_R},
        {"t_end", &t_end}
    };
    int paramsFound = 0;

    while (std::getline(file, line)) {
        if (line.empty()) continue;
        size_t equalPos = line.find('=');
        if (equalPos == std::string::npos) continue;
        std::string paramName = line.substr(0, equalPos);
        size_t start = paramName.find_first_not_of(" \t");
        size_t end = paramName.find_last_not_of(" \t");
        if (start != std::string::npos && end != std::string::npos) {
            paramName = paramName.substr(start, end - start + 1);
        }

        std::string valueStr = line.substr(equalPos + 1);
        try {
            double value = std::stod(valueStr);
            if (paramMap.find(paramName) != paramMap.end()) {
                *paramMap[paramName] = value;
                paramsFound++;
                std::cout << "Read parameter: " << paramName << " = " << value << std::endl;
            } else {
                std::cout << "I?aaoi?a?aaiea: iaecaanoiue ia?aiao? '" << paramName << "'" << std::endl;
            }
        } catch (const std::exception& e) {
            std::cerr << "Ioeaea i?aia?aciaaiey cia?aiey aey ia?aiao?a '" << paramName << "'" << std::endl;
        }
    }    
    file.close();
    if (paramsFound == 7) {
        return true;
    } else {
        std::cerr << "I?aaoi?a?aaiea: i?i?eoaii oieuei " << paramsFound << " ec 7 ia?aiao?ia" << std::endl;
        return paramsFound > 0;
    }
}

bool readInit(const std::string& filename, double& tau, double& h, 
              int& cells, double& cfl, std::string& left_bc, 
              std::string& right_bc, std::string& test_name, int& count, std::vector<std::string>& solver_names,
              std::string& t_end_user, double& t_output, std::string& key) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: failed opening of file " << filename << std::endl;
        return false;
    }

    std::string line;
    std::string solver_names_str;
    std::map<std::string, void*> paramMap = {
        {"tau", &tau},
        {"h", &h},
        {"cells", &cells},
        {"cfl", &cfl},
        {"left_bc", &left_bc},
        {"right_bc", &right_bc},
        {"test_name", &test_name},
        {"count", &count},
        {"solver_names", &solver_names_str},
        {"t_end", &t_end_user}, 
        {"t_output", &t_output}, 
        {"key", &key}
    };
    int paramsFound = 0;

    while (std::getline(file, line)) {
        if (line.empty()) continue;
        if (line.empty() || line.find("?oeiaianoai") != std::string::npos) {
        break;  // i?ae?auaai ?oaiea i?e ionoie no?iea eee ia?aea ?oeiaianoaa
        }
        size_t equalPos = line.find('=');
        if (equalPos == std::string::npos) continue;
        std::string paramName = line.substr(0, equalPos);
        size_t start = paramName.find_first_not_of(" \t");
        size_t end = paramName.find_last_not_of(" \t");
        if (start != std::string::npos && end != std::string::npos) {
            paramName = paramName.substr(start, end - start + 1);
        }

        std::string valueStr = line.substr(equalPos + 1);
        // Oae?aai i?iaaeu aie?oa cia?aiey
        start = valueStr.find_first_not_of(" \t");
        end = valueStr.find_last_not_of(" \t");
        if (start != std::string::npos && end != std::string::npos) {
            valueStr = valueStr.substr(start, end - start + 1);
        }

        try {
            if (paramMap.find(paramName) != paramMap.end()) {
                if (paramName == "tau" || paramName == "h" || paramName == "cfl" || paramName == "t_output") {
                    *static_cast<double*>(paramMap[paramName]) = std::stod(valueStr);
                } else if (paramName == "cells" || paramName == "count") {
                    *static_cast<int*>(paramMap[paramName]) = std::stoi(valueStr);
                } else {
                    *static_cast<std::string*>(paramMap[paramName]) = valueStr;
                }
                paramsFound++;
                std::cout << "Read parameter: " << paramName << " = " << valueStr << std::endl;
            } else {
                std::cout << "Warning: unknown parameter '" << paramName << "'" << std::endl;
            }
        } catch (const std::exception& e) {
            std::cerr << "Ioeaea i?aia?aciaaiey cia?aiey aey ia?aiao?a '" << paramName << "'" << std::endl;
        }
    }    
    file.close();


        if (!solver_names_str.empty()) {
        size_t start = 0;
        size_t end = solver_names_str.find(',');
        while (end != std::string::npos) {
            std::string solver = solver_names_str.substr(start, end - start);
            // Oae?aai i?iaaeu aie?oa eiaie noaiu
            size_t solver_start = solver.find_first_not_of(" \t");
            size_t solver_end = solver.find_last_not_of(" \t");
            if (solver_start != std::string::npos && solver_end != std::string::npos) {
                solver = solver.substr(solver_start, solver_end - solver_start + 1);
            }
            solver_names.push_back(solver);
            start = end + 1;
            end = solver_names_str.find(',', start);
        }
        // Aiaaaeyai iineaai?? noaio
        std::string last_solver = solver_names_str.substr(start);
        size_t solver_start = last_solver.find_first_not_of(" \t");
        size_t solver_end = last_solver.find_last_not_of(" \t");
        if (solver_start != std::string::npos && solver_end != std::string::npos) {
            last_solver = last_solver.substr(solver_start, solver_end - solver_start + 1);
        }
        solver_names.push_back(last_solver);
        
        std::cout << "Found " << solver_names.size() << " solvers: ";
        for (size_t i = 0; i < solver_names.size(); i++) {
            std::cout << solver_names[i];
            if (i < solver_names.size() - 1) std::cout << ", ";
        }
        std::cout << std::endl;
    }
    
    if (paramsFound >= 7) {
        return true;
    } else {
        std::cerr << "Warning! Read only " << paramsFound << " 7 params" << std::endl;
        return paramsFound > 0;
    }
}

//aey iieo?aiey oeeo y?aae

int getGhostCellsBySolver (const std::string& solver_name){
    if (solver_name == "godunov-kolgan" || 
        solver_name == "godunov-kolgan-rodionov" || 
        solver_name == "eno" ||
        solver_name == "weno"){return 2;}
    if (solver_name == "godunov"){return 1;}

    return 1;
}

// Ooieoey aey onoaiiaee aioiauo oanoia Niaa
void setupSODTest(const std::string& test_name, 
                  double& rho_L, double& u_L, double& p_L,
                  double& rho_R, double& u_R, double& p_R, 
                  double& t_end,
                  const std::string& t_end_user) {
    
    if (test_name == "sod1") {
        rho_L = 1.0; u_L = 0.0; p_L = 1.0;
        rho_R = 0.125; u_R = 0.0; p_R = 0.1;
        t_end = 0.25;
    }
    else if (test_name == "sod2") {
        rho_L = 1.0; u_L = -2.0; p_L = 0.4;
        rho_R = 1.0; u_R = 2.0; p_R = 0.4;
        t_end = 0.15;
    }
    else if (test_name == "sod3") {
        rho_L = 1.0; u_L = 0.0; p_L = 1000.0;
        rho_R = 1.0; u_R = 0.0; p_R = 0.01;
        t_end = 0.012;
    }
    else if (test_name == "sod4") {
        rho_L = 5.99924; u_L = 19.5975; p_L = 460.894;
        rho_R = 5.99242; u_R = -6.19633; p_R = 46.0950;
        t_end = 0.035;
    }
    else if (test_name == "sod5") {
        rho_L = 1.0; u_L = 0.0; p_L = 0.01;
        rho_R = 1.0; u_R = 0.0; p_R = 100.0;
        t_end = 0.035;
    }
    else {
        // Ii oiie?aie? sod1
        rho_L = 1.0; u_L = 0.0; p_L = 1.0;
        rho_R = 0.125; u_R = 0.0; p_R = 0.1;
        t_end = 0.2;
    }
    
    // I?IAA??AI IIEUCIAAOAEUNEIA A?AI?
    if (t_end_user != "default" && !t_end_user.empty()) {
        try {
            double user_time = std::stod(t_end_user);
            t_end = user_time;  // eniieucoai iieuciaaoaeuneia a?aiy
            std::cout << "Using user-specified time: " << t_end << std::endl;
        } catch (const std::exception& e) {
            std::cout << "Invalid time format, using default: " << t_end << std::endl;
        }
    } else {
        std::cout << "Using default time: " << t_end << std::endl;
    }
    
    std::cout << "Test: " << test_name << std::endl;
    std::cout << "rho_L=" << rho_L << " u_L=" << u_L << " p_L=" << p_L << std::endl;
    std::cout << "rho_R=" << rho_R << " u_R=" << u_R << " p_R=" << p_R << std::endl;
    std::cout << "t_end=" << t_end << std::endl;
}

// Ooieoey aey i?eiaiaiey a?aie?iuo oneiaee
void applyBoundaryConditions(double* rho, double* u, double* p,
                            int cells, int ghost_cells,
                            const std::string& left_bc,
                            const std::string& right_bc) {
    
    // Eiaaenaoey:
    // Aioo?aiiea y?aeee: [ghost_cells, ghost_cells + cells - 1]
    // Eaaua oeeoeaiua:   [0, ghost_cells - 1]
    // I?aaua oeeoeaiua:  [ghost_cells + cells, ghost_cells + cells + ghost_cells - 1]

    // --- Eaaay a?aieoa ---
    if (left_bc == "wall") {
        for (int i = 0; i < ghost_cells; ++i) {
            rho[i] = rho[ghost_cells];
            u[i]   = -u[ghost_cells];
            p[i]   = p[ghost_cells];
        }
    }
    else if (left_bc == "flow") {   // yeno?aiieyoey ioeaaiai ii?yaea
        for (int i = 0; i < ghost_cells; ++i) {
            rho[i] = rho[ghost_cells];
            u[i]   = u[ghost_cells];
            p[i]   = p[ghost_cells];
        }
    }
    else if (left_bc == "periodic") {
        for (int i = 0; i < ghost_cells; ++i) {
            rho[i] = rho[cells + ghost_cells - 1 - i];
            u[i]   = u[cells + ghost_cells - 1 - i];
            p[i]   = p[cells + ghost_cells - 1 - i];
        }
    }
    else if (left_bc == "inlet") {
        // Oeene?iaaiiua ia?aiao?u aoiaa (ii?ii auianoe a a?aoiaiou eee eiinoaiou)
        const double rho_in = 1.0;
        const double u_in   = 0.5;   // iiei?eoaeuiay nei?inou — iioie ai?aai
        const double p_in   = 1.0;
        for (int i = 0; i <= ghost_cells; ++i) {
            rho[i] = rho_in;
            u[i]   = u_in;
            p[i]   = p_in;
        }
    }
    else if (left_bc == "outlet") {
        // Auoia: yeno?aiieyoey ioeaaiai ii?yaea ec ia?aie aioo?aiiae y?aeee
        for (int i = 0; i < ghost_cells; ++i) {
            rho[i] = rho[ghost_cells];
            u[i]   = u[ghost_cells];
            p[i]   = p[ghost_cells];
        }
    }

    // --- I?aaay a?aieoa ---
    if (right_bc == "wall") {
        for (int i = 0; i < ghost_cells; ++i) {
            rho[cells + ghost_cells + i] = rho[cells + ghost_cells - 1];
            u[cells + ghost_cells + i]   = -u[cells + ghost_cells - 1];
            p[cells + ghost_cells + i]   = p[cells + ghost_cells - 1];
        }
    }
    else if (right_bc == "flow") {
        for (int i = 0; i < ghost_cells; ++i) {
            rho[cells + ghost_cells + i] = rho[cells + ghost_cells - 1];
            u[cells + ghost_cells + i]   = u[cells + ghost_cells - 1];
            p[cells + ghost_cells + i]   = p[cells + ghost_cells - 1];
        }
    }
    else if (right_bc == "periodic") {
        for (int i = 0; i < ghost_cells; ++i) {
            rho[cells + ghost_cells + i] = rho[ghost_cells + i];
            u[cells + ghost_cells + i]   = u[ghost_cells + i];
            p[cells + ghost_cells + i]   = p[ghost_cells + i];
        }
    }
    else if (right_bc == "inlet") {
        // Aoia ni?aaa: nei?inou io?eoaoaeuiay (iioie aeaai)
        const double rho_in = 1.0;
        const double u_in   = -0.5;   // io?eoaoaeuiay nei?inou
        const double p_in   = 1.0;
        for (int i = 0; i < ghost_cells; ++i) {
            rho[cells + ghost_cells + i] = rho_in;
            u[cells + ghost_cells + i]   = u_in;
            p[cells + ghost_cells + i]   = p_in;
        }
    }
    else if (right_bc == "outlet") {
        // Auoia: yeno?aiieyoey ioeaaiai ii?yaea ec iineaaiae aioo?aiiae y?aeee
        for (int i = 0; i < ghost_cells; ++i) {
            rho[cells + ghost_cells + i] = rho[cells + ghost_cells - 1];
            u[cells + ghost_cells + i]   = u[cells + ghost_cells - 1];
            p[cells + ghost_cells + i]   = p[cells + ghost_cells - 1];
        }
    }
}