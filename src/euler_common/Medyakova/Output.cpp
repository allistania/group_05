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

// Явно указываем возвращаемый тип вместо auto
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
    
    // 1. Скорости звука
    double aL = std::sqrt(gamma * pL / rhoL);
    double aR = std::sqrt(gamma * pR / rhoR);
    
    // 2. Оценки волновых скоростей (Einfeldt, Davis, Toro)
    // Выберите один из вариантов:
    
    double rho_avg = 0.5 * (rhoL + rhoR);
    double a_avg = 0.5 * (aL + aR);
    
    // Скорости Роу (Roe averages)
    double u_roe = (std::sqrt(rhoL) * uL + std::sqrt(rhoR) * uR) / 
                   (std::sqrt(rhoL) + std::sqrt(rhoR));
    double H_roe = ((std::sqrt(rhoL) * (aL * aL / (gamma - 1) + 0.5 * uL * uL) +
                     std::sqrt(rhoR) * (aR * aR / (gamma - 1) + 0.5 * uR * uR)) /
                    (std::sqrt(rhoL) + std::sqrt(rhoR)));
    double a_roe = std::sqrt((gamma - 1) * (H_roe - 0.5 * u_roe * u_roe));
    
    double SL = u_roe - a_roe;
    double SR = u_roe + a_roe;
    
    // Расширение для сильных разрывов
    SL = std::min(SL, uL - aL);
    SR = std::max(SR, uR + aR);
    
    // 3. Консервативные переменные и потоки
    double EL = pL / (gamma - 1.0) + 0.5 * rhoL * uL * uL;
    double ER = pR / (gamma - 1.0) + 0.5 * rhoR * uR * uR;
    
    std::vector<double> UL = {rhoL, rhoL * uL, EL};
    std::vector<double> UR = {rhoR, rhoR * uR, ER};
    
    std::vector<double> FL = {rhoL * uL, rhoL * uL * uL + pL, uL * (EL + pL)};
    std::vector<double> FR = {rhoR * uR, rhoR * uR * uR + pR, uR * (ER + pR)};
    
    // 4. Вычисление HLL-потока для ? = 0 (между ячейками)
    std::vector<double> F(3);
    
    if (0.0 <= SL) {
        // Весь поток от левого состояния
        F = FL;
    } 
    else if (0.0 >= SR) {
        // Весь поток от правого состояния
        F = FR;
    } 
    else {
        // HLL-поток
        double denom = SR - SL;
        denom = std::max(denom, eps);
        
        for (int i = 0; i < 3; ++i) {
            F[i] = (SR * FL[i] - SL * FR[i] + SL * SR * (UR[i] - UL[i])) / denom;
        }
    }
    
    // 5. Если вам действительно нужно ВОССТАНОВИТЬ состояние (не поток),
    // используйте другой метод или вычислите из потока осторожно
    
    // Внимание: восстановление состояния из потока HLL ненадёжно!
    // Лучше использовать методы с точным или приближённым решением Римана
    
    return std::make_tuple(F[0], F[1], F[2]); // Возвращаем ПОТОКИ
}

// Функция 
// CHECK: HLLC_SOLVER 
std::tuple<double, double, double> hllc(double rhoL, double uL, double pL,
                                        double rhoR, double uR, double pR,
                                        double xi, double gamma) {
    
    // 1. Скорости звука
    double aL = std::sqrt(gamma * pL / rhoL);
    double aR = std::sqrt(gamma * pR / rhoR);
    
    // 2. Полные энергии
    double EL = pL / (gamma - 1.0) + 0.5 * rhoL * uL * uL;
    double ER = pR / (gamma - 1.0) + 0.5 * rhoR * uR * uR;
    
    // 3. Оценка скоростей волн (линеаризованная оценка из презентации)
    // Вычисляем средние значения по Роу
    double sqrt_rhoL = std::sqrt(rhoL);
    double sqrt_rhoR = std::sqrt(rhoR);
    double sum_sqrt = sqrt_rhoL + sqrt_rhoR;
    
    double u_tilde = (sqrt_rhoL * uL + sqrt_rhoR * uR) / sum_sqrt;
    double HL = (EL + pL) / rhoL;
    double HR = (ER + pR) / rhoR;
    double H_tilde = (sqrt_rhoL * HL + sqrt_rhoR * HR) / sum_sqrt;
    
    double a_tilde = std::sqrt((gamma - 1.0) * (H_tilde - 0.5 * u_tilde * u_tilde));
    
    // Оценки скоростей волн
    double S_L = std::min(uL - aL, u_tilde - a_tilde);
    double S_R = std::max(uR + aR, u_tilde + a_tilde);
    
    // 4. Скорость контактного разрыва (S_*)
    double numerator = pR - pL + rhoL * uL * (S_L - uL) - rhoR * uR * (S_R - uR);
    double denominator = rhoL * (S_L - uL) - rhoR * (S_R - uR);
    
    double S_star;
    if (std::fabs(denominator) < 1e-10) {
        // Если знаменатель близок к нулю, используем среднее
        S_star = 0.5 * (uL + uR);
    } else {
        S_star = numerator / denominator;
    }
    
    // 5. Давление в звездной области (p_*)
    double p_star = pL + rhoL * (uL - S_L) * (uL - S_star);
    // Альтернативно: p_star = pR + rhoR * (uR - S_R) * (uR - S_star)
    
    // 6. Плотности в звездных областях
    double rho_star_L = rhoL * (S_L - uL) / (S_L - S_star);
    double rho_star_R = rhoR * (S_R - uR) / (S_R - S_star);
    
    // 7. Энергии в звездных областях
    double E_star_L = ((S_L - uL) * EL + p_star * S_star - pL * uL) / (S_L - S_star);
    double E_star_R = ((S_R - uR) * ER + p_star * S_star - pR * uR) / (S_R - S_star);
    
    // 8. Консервативные переменные для разных областей
    std::vector<double> UL = {rhoL, rhoL * uL, EL};
    std::vector<double> UR = {rhoR, rhoR * uR, ER};
    std::vector<double> U_star_L = {rho_star_L, rho_star_L * S_star, E_star_L};
    std::vector<double> U_star_R = {rho_star_R, rho_star_R * S_star, E_star_R};
    
    // 9. Выбор состояния в зависимости от положения xi
    std::vector<double> state(3);
    
    if (xi <= S_L) {
        // Левее левой волны
        state = UL;
    } 
    else if (xi <= S_star) {
        // Между левой волной и контактным разрывом
        state = U_star_L;
    }
    else if (xi < S_R) {
        // Между контактным разрывом и правой волной
        state = U_star_R;
    }
    else {
        // Правее правой волны
        state = UR;
    }
    
    // 10. Преобразование консервативных переменных в примитивные
    double rho = state[0];
    double u = state[1] / rho;
    double E = state[2];
    double p = (E - 0.5 * rho * u * u) * (gamma - 1.0);
    
    // 11. Коррекция отрицательного давления
    if (p < 0.0) {
        p = 1e-10;
        E = p / (gamma - 1.0) + 0.5 * rho * u * u;
        state[2] = E;
    }
    
    // 12. Коррекция отрицательной плотности
    if (rho < 0.0) {
        rho = 1e-10;
        // Пересчитываем скорость
        if (rho > 1e-12) {
            u = state[1] / rho;
        }
        p = std::max(p, 1e-10);
    }
    
    return std::make_tuple(rho, u, p);
}
// Вспомогательные функции для консервативных переменных и потоков
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

// Метод Русанова
std::tuple<double, double, double> rusanov(double rhoL, double uL, double pL,
                                          double rhoR, double uR, double pR,
                                          double xi) {
    
    double gamma = 1.4;
    
    // Вычисляем скорости звука
    double aL = sqrt(gamma * pL / rhoL);
    double aR = sqrt(gamma * pR / rhoR);
    
    // Максимальная скорость распространения возмущений
    double lambda_max = std::max(std::abs(uL) + aL, std::abs(uR) + aR);
    
    // Консервативные переменные
    std::vector<double> UL(3), UR(3);
    UL[0] = rhoL;  // масса
    UL[1] = rhoL * uL;  // импульс
    UL[2] = pL/(gamma-1.0) + 0.5*rhoL*uL*uL;  // полная энергия
    
    UR[0] = rhoR;  // масса
    UR[1] = rhoR * uR;  // импульс
    UR[2] = pR/(gamma-1.0) + 0.5*rhoR*uR*uR;  // полная энергия
    
    // Потоки из левого и правого состояний
    std::vector<double> FL(3), FR(3);
    
    // Левый поток
    FL[0] = rhoL * uL;  // поток массы
    FL[1] = rhoL * uL * uL + pL;  // поток импульса
    FL[2] = uL * (UL[2] + pL);  // поток энергии
    
    // Правый поток
    FR[0] = rhoR * uR;  // поток массы
    FR[1] = rhoR * uR * uR + pR;  // поток импульса
    FR[2] = uR * (UR[2] + pR);  // поток энергии
    
    // Численный поток Русанова
    std::vector<double> F(3);
    for (int i = 0; i < 3; ++i) {
        F[i] = 0.5 * (FL[i] + FR[i]) - 0.5 * lambda_max * (UR[i] - UL[i]);
    }
    
    // Возвращаем три компоненты потока
    return std::make_tuple(F[0], F[1], F[2]);
}