#include "Output.h"

void outputCSV(std::string name, 
               const std::vector<double>& grid_P,
               const std::vector<double>& u,
               const std::vector<double>& P,
               const std::vector<double>& rho,
               const std::vector<double>& I,
               double t) {
    std::ofstream outfile(name);
    std::string tmp;
    tmp = "Time,X,Rho,P,U,I\n";
    outfile << tmp;
    for (int i = 0; i < P.size(); i++) {
        tmp = "";
        tmp += std::to_string(t) + ',';
        tmp += std::to_string(grid_P[i]) + ',';
        tmp += std::to_string(rho[i]) + ',';
        tmp += std::to_string(P[i]) + ',';
        tmp += std::to_string(u[i]) + ',';
        tmp += std::to_string(I[i]) + '\n';
        outfile << tmp;
    }
    outfile.close();
}

// ?aii oeacuaaai aica?auaaiue oei aianoi auto
std::function<bool(double, const std::vector<double>&, const std::vector<double>&,
                   const std::vector<double>&, const std::vector<double>&,
                   const std::vector<double>&)>
createTimeSeriesWriter(double output_interval, std::string base_name) {
    double next_output_time = 0.0;
    
    return [output_interval, base_name, next_output_time](
        double current_time,
        const std::vector<double>& grid_P,
        const std::vector<double>& u,
        const std::vector<double>& P,
        const std::vector<double>& rho,
        const std::vector<double>& I) mutable -> bool {
        
        if (current_time >= next_output_time) {
            std::ostringstream filename;
            filename << base_name << "_" << std::fixed << std::setprecision(1) << current_time << ".csv";
            
            outputCSV(filename.str(), grid_P, u, P, rho, I, current_time);
            
            next_output_time += output_interval;
            return true;
        }
        return false;
    };
}

void clearDirectory(const std::string& directoryPath) {
    try {
        if (std::filesystem::exists(directoryPath)) {
            for (const auto& entry : std::filesystem::directory_iterator(directoryPath)) {
                if (entry.path().extension() == ".csv") {
                    std::filesystem::remove(entry.path());
                }
            }
            std::cout<< "Cleaned directory: " << directoryPath << std::endl;
        } else {
            std::filesystem::create_directories(directoryPath);
            std::cout << "Created directory: " << directoryPath << std::endl;
        }
    } catch (const std::exception& e) {
        std::cerr << "Error cleaning directory " << directoryPath << ": " << e.what() << std::endl;
    }
}


std::tuple<double, double, double> hll(double rhoL, double uL, double pL,
                                       double rhoR, double uR, double pR,
                                       double xi, double gamma )  {
    const double eps = 1e-12;
    
    // 1. Nei?inoe caoea
    double aL = std::sqrt(gamma * pL / rhoL);
    double aR = std::sqrt(gamma * pR / rhoR);
    
    // 2. Ioaiee aieiiauo nei?inoae (Einfeldt, Davis, Toro)
    // Auaa?eoa iaei ec aa?eaioia:
    
    double rho_avg = 0.5 * (rhoL + rhoR);
    double a_avg = 0.5 * (aL + aR);
    
    // Nei?inoe ?io (Roe averages)
    double u_roe = (std::sqrt(rhoL) * uL + std::sqrt(rhoR) * uR) / 
                   (std::sqrt(rhoL) + std::sqrt(rhoR));
    double H_roe = ((std::sqrt(rhoL) * (aL * aL / (gamma - 1) + 0.5 * uL * uL) +
                     std::sqrt(rhoR) * (aR * aR / (gamma - 1) + 0.5 * uR * uR)) /
                    (std::sqrt(rhoL) + std::sqrt(rhoR)));
    double a_roe = std::sqrt((gamma - 1) * (H_roe - 0.5 * u_roe * u_roe));
    
    double SL = u_roe - a_roe;
    double SR = u_roe + a_roe;
    
    // ?anoe?aiea aey neeuiuo ?ac?uaia
    SL = std::min(SL, uL - aL);
    SR = std::max(SR, uR + aR);
    
    // 3. Eiina?aaoeaiua ia?aiaiiua e iioiee
    double EL = pL / (gamma - 1.0) + 0.5 * rhoL * uL * uL;
    double ER = pR / (gamma - 1.0) + 0.5 * rhoR * uR * uR;
    
    std::vector<double> UL = {rhoL, rhoL * uL, EL};
    std::vector<double> UR = {rhoR, rhoR * uR, ER};
    
    std::vector<double> FL = {rhoL * uL, rhoL * uL * uL + pL, uL * (EL + pL)};
    std::vector<double> FR = {rhoR * uR, rhoR * uR * uR + pR, uR * (ER + pR)};
    
    // 4. Au?eneaiea HLL-iioiea aey ? = 0 (ia?ao y?aeeaie)
    std::vector<double> F(3);
    
    if (0.0 <= SL) {
        // Aanu iioie io eaaiai ninoiyiey
        F = FL;
    } 
    else if (0.0 >= SR) {
        // Aanu iioie io i?aaiai ninoiyiey
        F = FR;
    } 
    else {
        // HLL-iioie
        double denom = SR - SL;
        denom = std::max(denom, eps);
        
        for (int i = 0; i < 3; ++i) {
            F[i] = (SR * FL[i] - SL * FR[i] + SL * SR * (UR[i] - UL[i])) / denom;
        }
    }
    
    // 5. Anee aai aaenoaeoaeuii io?ii AINNOAIIAEOU ninoiyiea (ia iioie),
    // eniieucoeoa a?oaie iaoia eee au?eneeoa ec iioiea inoi?i?ii
    
    // Aieiaiea: ainnoaiiaeaiea ninoiyiey ec iioiea HLL iaiaa??ii!
    // Eo?oa eniieuciaaou iaoiau n oi?iui eee i?eaee??iiui ?aoaieai ?eiaia
    
    return std::make_tuple(F[0], F[1], F[2]); // Aica?auaai IIOIEE
}

// Ooieoey HLLC 
std::tuple<double, double, double> hllc(double rhoL, double uL, double pL,
                                        double rhoR, double uR, double pR,
                                        double xi, double gamma) {
    
    // 1. Nei?inoe caoea
    double aL = std::sqrt(gamma * pL / rhoL);
    double aR = std::sqrt(gamma * pR / rhoR);
    
    // 2. Iieiua yia?aee
    double EL = pL / (gamma - 1.0) + 0.5 * rhoL * uL * uL;
    double ER = pR / (gamma - 1.0) + 0.5 * rhoR * uR * uR;
    
    // 3. Ioaiea nei?inoae aiei (eeiaa?eciaaiiay ioaiea ec i?acaioaoee)
    // Au?eneyai n?aaiea cia?aiey ii ?io
    double sqrt_rhoL = std::sqrt(rhoL);
    double sqrt_rhoR = std::sqrt(rhoR);
    double sum_sqrt = sqrt_rhoL + sqrt_rhoR;
    
    double u_tilde = (sqrt_rhoL * uL + sqrt_rhoR * uR) / sum_sqrt;
    double HL = (EL + pL) / rhoL;
    double HR = (ER + pR) / rhoR;
    double H_tilde = (sqrt_rhoL * HL + sqrt_rhoR * HR) / sum_sqrt;
    
    double a_tilde = std::sqrt((gamma - 1.0) * (H_tilde - 0.5 * u_tilde * u_tilde));
    
    // Ioaiee nei?inoae aiei
    double S_L = std::min(uL - aL, u_tilde - a_tilde);
    double S_R = std::max(uR + aR, u_tilde + a_tilde);
    
    // 4. Nei?inou eiioaeoiiai ?ac?uaa (S_*)
    double numerator = pR - pL + rhoL * uL * (S_L - uL) - rhoR * uR * (S_R - uR);
    double denominator = rhoL * (S_L - uL) - rhoR * (S_R - uR);
    
    double S_star;
    if (std::fabs(denominator) < 1e-10) {
        // Anee ciaiaiaoaeu aeecie e ioe?, eniieucoai n?aaiaa
        S_star = 0.5 * (uL + uR);
    } else {
        S_star = numerator / denominator;
    }
    
    // 5. Aaaeaiea a caacaiie iaeanoe (p_*)
    double p_star = pL + rhoL * (uL - S_L) * (uL - S_star);
    // Aeuoa?iaoeaii: p_star = pR + rhoR * (uR - S_R) * (uR - S_star)
    
    // 6. Ieioiinoe a caacaiuo iaeanoyo
    double rho_star_L = rhoL * (S_L - uL) / (S_L - S_star);
    double rho_star_R = rhoR * (S_R - uR) / (S_R - S_star);
    
    // 7. Yia?aee a caacaiuo iaeanoyo
    double E_star_L = ((S_L - uL) * EL + p_star * S_star - pL * uL) / (S_L - S_star);
    double E_star_R = ((S_R - uR) * ER + p_star * S_star - pR * uR) / (S_R - S_star);
    
    // 8. Eiina?aaoeaiua ia?aiaiiua aey ?aciuo iaeanoae
    std::vector<double> UL = {rhoL, rhoL * uL, EL};
    std::vector<double> UR = {rhoR, rhoR * uR, ER};
    std::vector<double> U_star_L = {rho_star_L, rho_star_L * S_star, E_star_L};
    std::vector<double> U_star_R = {rho_star_R, rho_star_R * S_star, E_star_R};
    
    // 9. Auai? ninoiyiey a caaeneiinoe io iiei?aiey xi
    std::vector<double> state(3);
    
    if (xi <= S_L) {
        // Eaaaa eaaie aieiu
        state = UL;
    } 
    else if (xi <= S_star) {
        // Ia?ao eaaie aieiie e eiioaeoiui ?ac?uaii
        state = U_star_L;
    }
    else if (xi < S_R) {
        // Ia?ao eiioaeoiui ?ac?uaii e i?aaie aieiie
        state = U_star_R;
    }
    else {
        // I?aaaa i?aaie aieiu
        state = UR;
    }
    
    // 10. I?aia?aciaaiea eiina?aaoeaiuo ia?aiaiiuo a i?eieoeaiua
    double rho = state[0];
    double u = state[1] / rho;
    double E = state[2];
    double p = (E - 0.5 * rho * u * u) * (gamma - 1.0);
    
    // 11. Ei??aeoey io?eoaoaeuiiai aaaeaiey
    if (p < 0.0) {
        p = 1e-10;
        E = p / (gamma - 1.0) + 0.5 * rho * u * u;
        state[2] = E;
    }
    
    // 12. Ei??aeoey io?eoaoaeuiie ieioiinoe
    if (rho < 0.0) {
        rho = 1e-10;
        // Ia?an?eouaaai nei?inou
        if (rho > 1e-12) {
            u = state[1] / rho;
        }
        p = std::max(p, 1e-10);
    }
    
    return std::make_tuple(rho, u, p);
}
// Aniiiiaaoaeuiua ooieoee aey eiina?aaoeaiuo ia?aiaiiuo e iioieia
std::vector<double> primitive_to_conservative(double rho, double u, double p, double gamma = 1.4) {
    double E = p/(gamma-1.0) + 0.5*rho*u*u;
    return {rho, rho*u, E};
}

std::vector<double> flux_from_primitive(double rho, double u, double p, double gamma = 1.4) {
    double E = p/(gamma-1.0) + 0.5*rho*u*u;
    double F1 = rho*u;
    double F2 = rho*u*u + p;
    double F3 = u*(E + p);
    return {F1, F2, F3};
}

// Iaoia ?onaiiaa
std::tuple<double, double, double> rusanov(double rhoL, double uL, double pL,
                                          double rhoR, double uR, double pR,
                                          double xi) {
    
    double gamma = 1.4;
    
    // Au?eneyai nei?inoe caoea
    double aL = sqrt(gamma * pL / rhoL);
    double aR = sqrt(gamma * pR / rhoR);
    
    // Iaeneiaeuiay nei?inou ?ani?ino?aiaiey aiciouaiee
    double lambda_max = std::max(std::abs(uL) + aL, std::abs(uR) + aR);
    
    // Eiina?aaoeaiua ia?aiaiiua
    std::vector<double> UL(3), UR(3);
    UL[0] = rhoL;  // ianna
    UL[1] = rhoL * uL;  // eiioeun
    UL[2] = pL/(gamma-1.0) + 0.5*rhoL*uL*uL;  // iieiay yia?aey
    
    UR[0] = rhoR;  // ianna
    UR[1] = rhoR * uR;  // eiioeun
    UR[2] = pR/(gamma-1.0) + 0.5*rhoR*uR*uR;  // iieiay yia?aey
    
    // Iioiee ec eaaiai e i?aaiai ninoiyiee
    std::vector<double> FL(3), FR(3);
    
    // Eaaue iioie
    FL[0] = rhoL * uL;  // iioie iannu
    FL[1] = rhoL * uL * uL + pL;  // iioie eiioeuna
    FL[2] = uL * (UL[2] + pL);  // iioie yia?aee
    
    // I?aaue iioie
    FR[0] = rhoR * uR;  // iioie iannu
    FR[1] = rhoR * uR * uR + pR;  // iioie eiioeuna
    FR[2] = uR * (UR[2] + pR);  // iioie yia?aee
    
    // ?eneaiiue iioie ?onaiiaa
    std::vector<double> F(3);
    for (int i = 0; i < 3; ++i) {
        F[i] = 0.5 * (FL[i] + FR[i]) - 0.5 * lambda_max * (UR[i] - UL[i]);
    }
    
    // Aica?auaai o?e eiiiiiaiou iioiea
    return std::make_tuple(F[0], F[1], F[2]);
}
void setDefaultEndTime(const std::string& test_case, double& t_end) {
    // ≈сли врем€ уже задано положительным числом, ничего не мен€ем
    if (t_end != 0.0) return;

    // —тандартные времена дл€ известных тестов —ода
    if (t_end == 0.0) {
        if (test_case == "sod1"){
            t_end = 0.25;
        } else if (test_case == "sod2") {
            t_end = 0.15;
        } else if (test_case == "sod3") {
            t_end = 0.012;
        } else if (test_case == "sod4") {
            t_end = 0.035;
        } else if (test_case == "sod5") {
            t_end = 0.035;
        } else {
        // ≈сли тест не опознан, оставл€ем t_end как есть (может быть 0)
        // или можно установить значение по умолчанию:
             t_end = 0.25;
        }
    }
}
#include <filesystem>  //                               
namespace fs = std::filesystem;

//                                                                
void prepareDirectories() {
    //                            
    if (!fs::exists("plots")) {
        fs::create_directory("plots");
        std::cout << "Created directory: plots" << std::endl;
    }
    //                                GIF
    if (!fs::exists("forGIF/G")) {
        fs::create_directories("forGIF/G");
        std::cout << "Created directory: forGIF/G" << std::endl;
    }
    //               forGIF/G           CSV-      
    std::cout << "Cleaning directory forGIF/G..." << std::endl;
    for (const auto& entry : fs::directory_iterator("forGIF/G")) {
        if (entry.is_regular_file() && entry.path().extension() == ".csv") {
            fs::remove(entry.path());
        }
    }
}