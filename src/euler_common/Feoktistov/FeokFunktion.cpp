#include <iostream>
#include <cmath>
#include <algorithm>
#include <tuple>
#include <vector>
#include "Batarin/parameters.h"
#include "Medyakova/Output.h" 
using namespace std;
void kabare(int cells, int ghost_cells,
           std::vector<double>& u,
           std::vector<double>& P,
           std::vector<double>& rho,
           std::vector<double>& I,
           double& t_total, double tau, double h,
           const std::string& left_bc,
           const std::string& right_bc) {
               bool use_artificial_viscosity = true;
           double viscosity_coeff = 0.01;
           bool use_anti_diffusion = false;
           double anti_diffusion_coeff = 0.05;
    double gamma = 1.4;
    int total_cells = cells + 2 * ghost_cells;
    
    // Проверка, что массивы достаточно велики
    if (u.size() < total_cells) u.resize(total_cells);
    if (P.size() < total_cells) P.resize(total_cells);
    if (rho.size() < total_cells) rho.resize(total_cells);
    if (I.size() < total_cells) I.resize(total_cells);
    
    // Сохраняем старые значения для временных слоев
    std::vector<double> u_old = u;
    std::vector<double> P_old = P;
    std::vector<double> rho_old = rho;
    
    // Шаг 1: Вычисление консервативных переменных на n+1/2
    std::vector<double> Theta_half(total_cells, 0.0);
    std::vector<double> U_half(total_cells, 0.0);
    std::vector<double> P_half(total_cells, 0.0);
    std::vector<double> C_half(total_cells, 0.0);
    std::vector<double> G_half(total_cells, 0.0);
    std::vector<double> E_half(total_cells, 0.0);
    
    // Текущие консервативные переменные (n)
    std::vector<double> Theta_n(total_cells, 0.0);
    std::vector<double> U_n(total_cells, 0.0);
    std::vector<double> E_n(total_cells, 0.0);
    
    // Вычисление консервативных переменных на n
    for (int i = 0; i < total_cells; ++i) {
        Theta_n[i] = rho_old[i];
        U_n[i] = rho_old[i] * u_old[i];
        E_n[i] = P_old[i] / (gamma - 1.0) + 0.5 * rho_old[i] * u_old[i] * u_old[i];
    }
    
    // Вычисление промежуточных значений (n+1/2) по формуле (1.6)
    for (int i = ghost_cells; i < cells + ghost_cells; ++i) {
        // Индексы для соседних ячеек
        int left = std::max(ghost_cells, i - 1);
        int right = std::min(cells + ghost_cells - 1, i + 1);
        
        // Потоки на левой границе ячейки (i-1/2)
        double rho_left = 0.5 * (rho_old[left] + rho_old[i]);
        double u_left = 0.5 * (u_old[left] + u_old[i]);
        double P_left = 0.5 * (P_old[left] + P_old[i]);
        
        double F1_left = rho_left * u_left;
        double F2_left = rho_left * u_left * u_left + P_left;
        
        double e_left = P_left / ((gamma - 1.0) * rho_left);
        double E_left = e_left + 0.5 * u_left * u_left;
        double F3_left = (E_left + P_left) * u_left;
        
        // Потоки на правой границе ячейки (i+1/2)
        double rho_right = 0.5 * (rho_old[i] + rho_old[right]);
        double u_right = 0.5 * (u_old[i] + u_old[right]);
        double P_right = 0.5 * (P_old[i] + P_old[right]);
        
        double F1_right = rho_right * u_right;
        double F2_right = rho_right * u_right * u_right + P_right;
        
        double e_right = P_right / ((gamma - 1.0) * rho_right);
        double E_right = e_right + 0.5 * u_right * u_right;
        double F3_right = (E_right + P_right) * u_right;
        
        // Обновление на n+1/2 (формула 1.6)
        Theta_half[i] = Theta_n[i] - (tau / (2.0 * h)) * (F1_right - F1_left);
        
        double U_half_cons = U_n[i] - (tau / (2.0 * h)) * (F2_right - F2_left);
        U_half[i] = U_half_cons / Theta_half[i];
        
        double E_half_cons = E_n[i] - (tau / (2.0 * h)) * (F3_right - F3_left);
        E_half[i] = E_half_cons / Theta_half[i];
        
        // Вычисление давления на n+1/2
        double e_int = E_half[i] - 0.5 * U_half[i] * U_half[i];
        P_half[i] = (gamma - 1.0) * Theta_half[i] * e_int;
        
        // Скорость звука и G на n+1/2
        C_half[i] = sqrt(gamma * P_half[i] / std::max(Theta_half[i], 1e-12));
        G_half[i] = 1.0 / (Theta_half[i] * C_half[i]);
    }
    
    // Граничные условия для промежуточных значений
    applyBoundaryConditions(Theta_half.data(), U_half.data(), P_half.data(),
                           cells, ghost_cells, left_bc, right_bc);
    
    // Пересчет C_half и G_half после граничных условий
    for (int i = 0; i < total_cells; ++i) {
        C_half[i] = sqrt(gamma * P_half[i] / std::max(Theta_half[i], 1e-12));
        G_half[i] = 1.0 / (Theta_half[i] * C_half[i]);
    }
    
    // Шаг 2: Вычисление квазиинвариантов на n+1/2
    std::vector<double> R_half(total_cells, 0.0);
    std::vector<double> Q_half(total_cells, 0.0);
    std::vector<double> S_half(total_cells, 0.0);
    
    for (int i = 0; i < total_cells; ++i) {
        R_half[i] = U_half[i] + G_half[i] * P_half[i];
        Q_half[i] = U_half[i] - G_half[i] * P_half[i];
        S_half[i] = P_half[i] - C_half[i] * C_half[i] * Theta_half[i];
    }
    
    // Шаг 3: Экстраполяция квазиинвариантов на n+1
    std::vector<double> r_tilde_left(total_cells, 0.0);
    std::vector<double> q_tilde_left(total_cells, 0.0);
    std::vector<double> s_tilde_left(total_cells, 0.0);
    std::vector<double> r_tilde_right(total_cells, 0.0);
    std::vector<double> q_tilde_right(total_cells, 0.0);
    std::vector<double> s_tilde_right(total_cells, 0.0);
    
    for (int i = ghost_cells; i < cells + ghost_cells; ++i) {
        int left = std::max(ghost_cells, i - 1);
        int right = std::min(cells + ghost_cells - 1, i + 1);
        
        // Экстраполяция из левой ячейки (формула 1.12)
        double r_i_minus_1 = u_old[left] + G_half[left] * P_old[left];
        double q_i_minus_1 = u_old[left] - G_half[left] * P_old[left];
        double s_i_minus_1 = P_old[left] - C_half[left] * C_half[left] * rho_old[left];
        
        r_tilde_left[i] = 2.0 * R_half[left] - r_i_minus_1;
        q_tilde_left[i] = 2.0 * Q_half[left] - q_i_minus_1;
        s_tilde_left[i] = 2.0 * S_half[left] - s_i_minus_1;
        
        // Экстраполяция из правой ячейки
        double r_i_plus_1 = u_old[right] + G_half[right] * P_old[right];
        double q_i_plus_1 = u_old[right] - G_half[right] * P_old[right];
        double s_i_plus_1 = P_old[right] - C_half[right] * C_half[right] * rho_old[right];
        
        r_tilde_right[i] = 2.0 * R_half[i] - r_i_plus_1;
        q_tilde_right[i] = 2.0 * Q_half[i] - q_i_plus_1;
        s_tilde_right[i] = 2.0 * S_half[i] - s_i_plus_1;
    }
    
    // Шаг 4: Коррекция квазиинвариантов (формулы 1.15-1.17)
    for (int i = ghost_cells; i < cells + ghost_cells; ++i) {
        int left = std::max(ghost_cells, i - 1);
        int right = std::min(cells + ghost_cells - 1, i + 1);
        
        // Коррекция для r
        double r_left = u_old[left] + G_half[left] * P_old[left];
        double r_center = R_half[left];
        double r_right = u_old[right] + G_half[right] * P_old[right];
        
        double r_max = std::max({r_left, r_center, r_right});
        double r_min = std::min({r_left, r_center, r_right});
        
        if (r_tilde_left[i] < r_min) r_tilde_left[i] = r_min;
        if (r_tilde_left[i] > r_max) r_tilde_left[i] = r_max;
        if (r_tilde_right[i] < r_min) r_tilde_right[i] = r_min;
        if (r_tilde_right[i] > r_max) r_tilde_right[i] = r_max;
        
        // Коррекция для q
        double q_left = u_old[left] - G_half[left] * P_old[left];
        double q_center = Q_half[left];
        double q_right = u_old[right] - G_half[right] * P_old[right];
        
        double q_max = std::max({q_left, q_center, q_right});
        double q_min = std::min({q_left, q_center, q_right});
        
        if (q_tilde_left[i] < q_min) q_tilde_left[i] = q_min;
        if (q_tilde_left[i] > q_max) q_tilde_left[i] = q_max;
        if (q_tilde_right[i] < q_min) q_tilde_right[i] = q_min;
        if (q_tilde_right[i] > q_max) q_tilde_right[i] = q_max;
        
        // Коррекция для s
        double s_left = P_old[left] - C_half[left] * C_half[left] * rho_old[left];
        double s_center = S_half[left];
        double s_right = P_old[right] - C_half[right] * C_half[right] * rho_old[right];
        
        double s_max = std::max({s_left, s_center, s_right});
        double s_min = std::min({s_left, s_center, s_right});
        
        if (s_tilde_left[i] < s_min) s_tilde_left[i] = s_min;
        if (s_tilde_left[i] > s_max) s_tilde_left[i] = s_max;
        if (s_tilde_right[i] < s_min) s_tilde_right[i] = s_min;
        if (s_tilde_right[i] > s_max) s_tilde_right[i] = s_max;
    }
    
    // Шаг 5: Определение направления характеристик
    std::vector<double> u_tilde(total_cells, 0.0);
    std::vector<double> c_tilde(total_cells, 0.0);
    
    for (int i = ghost_cells; i < cells + ghost_cells; ++i) {
        int left = std::max(ghost_cells, i - 1);
        int right = std::min(cells + ghost_cells - 1, i + 1);
        
        double c_i = sqrt(gamma * P_old[i] / std::max(rho_old[i], 1e-12));
        u_tilde[i] = U_half[left] + U_half[i] - u_old[i];
        c_tilde[i] = C_half[left] + C_half[i] - c_i;
    }
    
    // Шаг 6: Вычисление новых потоковых переменных (u, P, rho) на n+1
    std::vector<double> u_new(total_cells, 0.0);
    std::vector<double> P_new(total_cells, 0.0);
    std::vector<double> rho_new(total_cells, 0.0);
    
    for (int i = ghost_cells; i < cells + ghost_cells; ++i) {
        int left = std::max(ghost_cells, i - 1);
        int right = std::min(cells + ghost_cells - 1, i + 1);
        
        double lambda1 = u_tilde[i] - c_tilde[i];
        double lambda2 = u_tilde[i];
        double lambda3 = u_tilde[i] + c_tilde[i];
        
        // Выбор формул в зависимости от направления характеристик
        if (lambda1 >= 0 && lambda2 >= 0 && lambda3 >= 0) {
            // Сверхзвуковой поток вправо (формула 1.23)
            P_new[i] = (r_tilde_left[i] - q_tilde_left[i]) / (2.0 * G_half[left]);
            u_new[i] = (r_tilde_left[i] + q_tilde_left[i]) / 2.0;
            rho_new[i] = (P_new[i] - s_tilde_left[i]) / (C_half[left] * C_half[left]);
        }
        else if (lambda1 >= 0 && lambda2 >= 0 && lambda3 < 0) {
            // Формула (1.20) первая система
            P_new[i] = (r_tilde_left[i] - q_tilde_right[i]) / (G_half[left] + G_half[i]);
            u_new[i] = (G_half[i] * r_tilde_left[i] + G_half[left] * q_tilde_right[i]) 
                       / (G_half[left] + G_half[i]);
            rho_new[i] = (P_new[i] - s_tilde_left[i]) / (C_half[left] * C_half[left]);
        }
        else if (lambda1 >= 0 && lambda2 < 0 && lambda3 < 0) {
            // Формула (1.22)
            P_new[i] = (r_tilde_left[i] - q_tilde_right[i]) / (G_half[left] + G_half[i]);
            u_new[i] = (G_half[i] * r_tilde_left[i] + G_half[left] * q_tilde_right[i])
                       / (G_half[left] + G_half[i]);
            rho_new[i] = (P_new[i] - s_tilde_right[i]) / (C_half[i] * C_half[i]);
        }
        else if (lambda1 < 0 && lambda2 < 0 && lambda3 < 0) {
            // Сверхзвуковой поток влево (формула после 1.23)
            P_new[i] = (r_tilde_right[i] - q_tilde_right[i]) / (2.0 * G_half[i]);
            u_new[i] = (r_tilde_right[i] + q_tilde_right[i]) / 2.0;
            rho_new[i] = (P_new[i] - s_tilde_right[i]) / (C_half[i] * C_half[i]);
        }
        else {
            // Запасной вариант - сохраняем старые значения
            u_new[i] = u_old[i];
            P_new[i] = P_old[i];
            rho_new[i] = rho_old[i];
        }
        
        // Физические ограничения
        if (rho_new[i] < 1e-8) rho_new[i] = 1e-8;
        if (P_new[i] < 1e-8) P_new[i] = 1e-8;
    }
    
    // === ДОБАВЛЕННАЯ ДИФФУЗИЯ ДЛЯ ПОДАВЛЕНИЯ ОСЦИЛЛЯЦИЙ ===
    // Коэффициенты искусственной диффузии (можно регулировать)
    double diff_coeff_u = 0.01;    // для скорости
    double diff_coeff_p = 0.01;    // для давления
    double diff_coeff_rho = 0.01;  // для плотности
    
    std::vector<double> u_smoothed = u_new;
    std::vector<double> P_smoothed = P_new;
    std::vector<double> rho_smoothed = rho_new;
    
    // Применяем явную схему диффузии (лапласиан)
    for (int i = ghost_cells; i < cells + ghost_cells; ++i) {
        int left = std::max(ghost_cells, i - 1);
        int right = std::min(cells + ghost_cells - 1, i + 1);
        
        // Диффузия для скорости
        double diff_u = diff_coeff_u * (u_new[left] - 2.0 * u_new[i] + u_new[right]);
        u_smoothed[i] = u_new[i] + diff_u;
        
        // Диффузия для давления
        double diff_p = diff_coeff_p * (P_new[left] - 2.0 * P_new[i] + P_new[right]);
        P_smoothed[i] = P_new[i] + diff_p;
        
        // Диффузия для плотности
        double diff_rho = diff_coeff_rho * (rho_new[left] - 2.0 * rho_new[i] + rho_new[right]);
        rho_smoothed[i] = rho_new[i] + diff_rho;
        
        // Ограничения после сглаживания
        if (rho_smoothed[i] < 1e-8) rho_smoothed[i] = 1e-8;
        if (P_smoothed[i] < 1e-8) P_smoothed[i] = 1e-8;
    }
    
    // Используем сглаженные значения
    u_new = u_smoothed;
    P_new = P_smoothed;
    rho_new = rho_smoothed;
    // ======================================================
    
    // Шаг 7: Обновление консервативных переменных на n+1 (формула 1.5)
    std::vector<double> Theta_new(total_cells, 0.0);
    std::vector<double> U_new_cons(total_cells, 0.0);
    std::vector<double> E_new_cons(total_cells, 0.0);
    
    for (int i = ghost_cells; i < cells + ghost_cells; ++i) {
        int left = std::max(ghost_cells, i - 1);
        int right = std::min(cells + ghost_cells - 1, i + 1);
        
        // Потоки на границах ячеек на слое n+1
        double rho_left = 0.5 * (rho_new[left] + rho_new[i]);
        double u_left = 0.5 * (u_new[left] + u_new[i]);
        double P_left = 0.5 * (P_new[left] + P_new[i]);
        
        double F1_left = rho_left * u_left;
        double F2_left = rho_left * u_left * u_left + P_left;
        double e_left = P_left / ((gamma - 1.0) * rho_left);
        double E_left = e_left + 0.5 * u_left * u_left;
        double F3_left = (E_left + P_left) * u_left;
        
        double rho_right = 0.5 * (rho_new[i] + rho_new[right]);
        double u_right = 0.5 * (u_new[i] + u_new[right]);
        double P_right = 0.5 * (P_new[i] + P_new[right]);
        
        double F1_right = rho_right * u_right;
        double F2_right = rho_right * u_right * u_right + P_right;
        double e_right = P_right / ((gamma - 1.0) * rho_right);
        double E_right = e_right + 0.5 * u_right * u_right;
        double F3_right = (E_right + P_right) * u_right;
        
        // Обновление консервативных переменных (формула 1.5)
        Theta_new[i] = Theta_n[i] - (tau / h) * (F1_right - F1_left);
        
        U_new_cons[i] = U_n[i] - (tau / h) * (F2_right - F2_left);
        
        E_new_cons[i] = E_n[i] - (tau / h) * (F3_right - F3_left);
    }
    
    // === ПРИМЕНЕНИЕ ИСКУССТВЕННОЙ ВЯЗКОСТИ ИЗ ФУНКЦИИ ФЕОКТИСТОВА ===
    if (use_artificial_viscosity) {
        // Применяем искусственную вязкость к консервативным переменным
        std::vector<double> Theta_old = Theta_new;
        std::vector<double> U_cons_old = U_new_cons;
        std::vector<double> E_cons_old = E_new_cons;
        
        for (int i = ghost_cells; i < cells + ghost_cells; ++i) {
            if (i > ghost_cells && i < cells + ghost_cells - 1) {
                Theta_new[i] = Theta_old[i] + viscosity_coeff * 
                              (Theta_old[i+1] - 2.0 * Theta_old[i] + Theta_old[i-1]);
                U_new_cons[i] = U_cons_old[i] + viscosity_coeff * 
                               (U_cons_old[i+1] - 2.0 * U_cons_old[i] + U_cons_old[i-1]);
                E_new_cons[i] = E_cons_old[i] + viscosity_coeff * 
                               (E_cons_old[i+1] - 2.0 * E_cons_old[i] + E_cons_old[i-1]);
            }
        }
    }
    
    if (use_anti_diffusion) {
        // Применяем антидиффузию
        std::vector<double> Theta_temp = Theta_new;
        std::vector<double> U_cons_temp = U_new_cons;
        std::vector<double> E_cons_temp = E_new_cons;
        
        for (int i = ghost_cells + 1; i < cells + ghost_cells - 1; ++i) {
            Theta_new[i] = Theta_temp[i] + anti_diffusion_coeff * 
                          (Theta_temp[i+1] - 2.0 * Theta_temp[i] + Theta_temp[i-1]);
            U_new_cons[i] = U_cons_temp[i] + anti_diffusion_coeff * 
                           (U_cons_temp[i+1] - 2.0 * U_cons_temp[i] + U_cons_temp[i-1]);
            E_new_cons[i] = E_cons_temp[i] + anti_diffusion_coeff * 
                           (E_cons_temp[i+1] - 2.0 * E_cons_temp[i] + E_cons_temp[i-1]);
        }
    }
    // ================================================================
    
    // Шаг 8: Преобразование обратно в примитивные переменные и обновление выходных массивов
    for (int i = ghost_cells; i < cells + ghost_cells; ++i) {
        rho[i] = Theta_new[i];
        if (rho[i] < 1e-8) rho[i] = 1e-8;
        
        u[i] = U_new_cons[i] / rho[i];
        
        double e_kin = 0.5 * rho[i] * u[i] * u[i];
        double e_int = E_new_cons[i] - e_kin;
        if (e_int < 1e-8) e_int = 1e-8;
        
        P[i] = (gamma - 1.0) * e_int;
        I[i] = E_new_cons[i];
    }
    
    // Граничные условия
    applyBoundaryConditions(rho.data(), u.data(), P.data(),
                           cells, ghost_cells, left_bc, right_bc);
    
    // Обновление времени
    t_total += tau;
    
    std::cout << "Kabare step: t = " << t_total << " tau = " << tau 
              << " viscosity = " << (use_artificial_viscosity ? "on" : "off")
              << " anti-diffusion = " << (use_anti_diffusion ? "on" : "off") << std::endl;
}
double WaveFunk(double p, double pk, double rho_k){
	       	double gamma = 1.4;

    if (p > pk){
                double Ak = 2.0 / ((gamma + 1) * rho_k);
        	double Bk = (gamma - 1) / (gamma + 1) * pk;
        	return (p - pk) * std::sqrt(Ak / (p + Bk));
    } 
    else{
        return 2.0 * std::sqrt(gamma * pk / rho_k) / (gamma - 1)*(std::pow(p / pk, (gamma - 1) / (2 * gamma)) - 1.0);
    }
}

double ProisvWaveFunk(double p, double pk, double rho_k){
	double gamma = 1.4;

       if (p > pk) {
		double Ak = 2.0 / ((gamma + 1) * rho_k);
		double Bk = (gamma - 1) / (gamma + 1) * pk;
        	return std::sqrt(Ak / (p + Bk)) * (1.0 - (p - pk) / (2.0 * (p + Bk)));
    } else {
	    	return (1.0 / (rho_k * std::sqrt(gamma * pk / rho_k))) *
                std::pow(p / pk, -(gamma + 1) / (2 * gamma));
    }
}

double solve_p_star(double rhoL, double uL, double pL, double rhoR, double uR, double pR){
	double gamma = 1.4, tol = 1e-8;
       	int max_iter = 100;	
	double p0 = 0.5 * (pL + pR) - 0.125 * (uR - uL) * (rhoL + rhoR) *
                (std::sqrt(gamma * pL / rhoL) + std::sqrt(gamma * pR / rhoR));	
	p0 = std::max(tol, p0);
        double fL, fR, f, dfL, dfR, df;
	double p_star = p0;
	
    	for (int i = 0; i < max_iter; i++) {
        	fL = WaveFunk(p_star, pL, rhoL);
        	fR = WaveFunk(p_star, pR, rhoR);
        	f = fL + fR + uR - uL;
        
        	dfL = ProisvWaveFunk(p_star, pL, rhoL);
        	dfR = ProisvWaveFunk(p_star, pR, rhoR);
        	df = dfL + dfR;
        
               if (std::abs(df) < 1e-12) {
            		break;
        }
        
        double dp = -f / df;
        p_star += dp;
        
       
        if (p_star < tol) {
            p_star = tol;
        }
        
        if (std::abs(dp) < tol * p_star) {
            break;
        }
    }
    
    return p_star; 
}



tuple<double, double, double> sample(
    double p_star, double u_star, double rhoL, double uL, double pL,
    double rhoR, double uR, double pR, double x_over_t, double gamma = 1.4
) {
    double cL = sqrt(gamma * pL / rhoL);
    double cR = sqrt(gamma * pR / rhoR);
    
    if (x_over_t < u_star) {
               if (p_star > pL) {
                       double SL = uL - sqrt(((gamma + 1) / (2 * gamma)) * (p_star / pL) + 
                                 (gamma - 1) / (2 * gamma)) * cL;
            if (x_over_t < SL) {
                         return {rhoL, uL, pL};
            } else {
               
                double rho_star = rhoL * ((p_star / pL) + (gamma - 1) / (gamma + 1)) /
                                 (((gamma - 1) / (gamma + 1)) * (p_star / pL) + 1);
                return {rho_star, u_star, p_star};
            }
        } else {
            //                       
            double SHL = uL - cL;  //                        
            double STL = u_star - sqrt(gamma * p_star / 
                                      (rhoL * pow(p_star / pL, 1.0 / gamma)));  //      
            
            if (x_over_t < SHL) {
                //                    
                return {rhoL, uL, pL};
            } else if (x_over_t < STL) {
                //                        
                double u = 2.0 / (gamma + 1) * (cL + (gamma - 1) / 2 * uL + x_over_t);
                double c = 2.0 / (gamma + 1) * (cL + (gamma - 1) / 2 * (uL - x_over_t));
                double rho = rhoL * pow(c / cL, 2.0 / (gamma - 1));
                double p = pL * pow(rho / rhoL, gamma);
                return {rho, u, p};
            } else {
                //                     
                double rho_star = rhoL * pow(p_star / pL, 1.0 / gamma);
                return {rho_star, u_star, p_star};
            }
        }
    } else {
        //                                   
        if (p_star > pR) {
            //                     
            double SR = uR + sqrt(((gamma + 1) / (2 * gamma)) * (p_star / pR) + 
                                 (gamma - 1) / (2 * gamma)) * cR;
            if (x_over_t > SR) {
                //                 
                return {rhoR, uR, pR};
            } else {
                //                  
                double rho_star = rhoR * ((p_star / pR) + (gamma - 1) / (gamma + 1)) /
                                 (((gamma - 1) / (gamma + 1)) * (p_star / pR) + 1);
                return {rho_star, u_star, p_star};
            }
        } else {
            //                        
            double SHR = uR + cR;  //                        
            double STR = u_star + sqrt(gamma * p_star / 
                                      (rhoR * pow(p_star / pR, 1.0 / gamma)));  //      
            
            if (x_over_t > SHR) {
                //                    
                return {rhoR, uR, pR};
            } else if (x_over_t > STR) {
                //                        
                double u = 2.0 / (gamma + 1) * (-cR + (gamma - 1) / 2 * uR + x_over_t);
                double c = 2.0 / (gamma + 1) * (cR - (gamma - 1) / 2 * (uR - x_over_t));
                double rho = rhoR * pow(c / cR, 2.0 / (gamma - 1));
                double p = pR * pow(rho / rhoR, gamma);
                return {rho, u, p};
            } else {
                //                     
                double rho_star = rhoR * pow(p_star / pR, 1.0 / gamma);
                return {rho_star, u_star, p_star};
            }
        }
    }
}
double  u_star(double p_star,double  rhoL,double  uL,double  pL,double  rhoR,double uR,double  pR){
   double gamma=1.4;
    double fL = WaveFunk(p_star, pL, rhoL);
    double fR =WaveFunk (p_star, pR, rhoR);
    return 0.5 * (uL + uR) + 0.5 * (fR - fL);
}

std::tuple<double,double,double> roe(double rhoL,double uL,double pL,double rhoR,double uR,double pR,double xi,double gamma=1.4){
	
	double EL=pL/(gamma-1.0)+0.5*rhoL*uL*uL;
	double ER=pR/(gamma-1.0)+0.5*rhoR*uR*uR;
	
	double sqrt_rhoL=sqrt(rhoL);
	double sqrt_rhoR=sqrt(rhoR);
	double sum_sqrt=sqrt_rhoL+sqrt_rhoR;
	
	double u_tilde=(sqrt_rhoL*uL+sqrt_rhoR*uR)/sum_sqrt;
	double HL=(EL+pL)/rhoL;
	double HR=(ER+pR)/rhoR;
	double H_tilde=(sqrt_rhoL*HL+sqrt_rhoR*HR)/sum_sqrt;
	
	double rho_tilde=sqrt(rhoL*rhoR);
	double a2_tilde=(gamma-1.0)*(H_tilde-0.5*u_tilde*u_tilde);
	double a_tilde=sqrt(a2_tilde);
	
	double delta_p=pR-pL;
	double delta_u=uR-uL;
	double delta_rho=rhoR-rhoL;
	
	double alpha1=(delta_p-rho_tilde*a_tilde*delta_u)/(2.0*a2_tilde);
	double alpha2=delta_rho-delta_p/a2_tilde;
	double alpha3=(delta_p+rho_tilde*a_tilde*delta_u)/(2.0*a2_tilde);
	
	double lambda1=u_tilde-a_tilde;
	double lambda2=u_tilde;
	double lambda3=u_tilde+a_tilde;
	
	double U1=rhoL;
	double U2=rhoL*uL;
	double U3=EL;
	
	if(xi>lambda1){
		U1+=alpha1*1.0;
		U2+=alpha1*(u_tilde-a_tilde);
		U3+=alpha1*(H_tilde-u_tilde*a_tilde);}
	
	if(xi>lambda2){
		U1+=alpha2*1.0;
		U2+=alpha2*u_tilde;
		U3+=alpha2*(0.5*u_tilde*u_tilde);}
	
	if(xi>lambda3){
		U1+=alpha3*1.0;
		U2+=alpha3*(u_tilde+a_tilde);
		U3+=alpha3*(H_tilde+u_tilde*a_tilde);}
	
	double rho=U1;
	double u=U2/rho;
	double e=U3/rho-0.5*u*u;
	
	double p=rho*e*(gamma-1.0);
	if(p<0.0){
		p=1e-10;
		e=p/(rho*(gamma-1.0));
		U3=rho*(e+0.5*u*u);}

	return std::make_tuple(rho,u,p);
}

#include <tuple>
#include <array>
#include <cmath>
#include <algorithm>

// Метод Ошера-Соломона для вычисления численного потока
std::tuple<double, double, double> osher_solomon(
    double rhoL, double uL, double pL,
    double rhoR, double uR, double pR,
    double xi, double gamma) {
    
    const double eps = 1e-12;
    const double gamma_minus_1 = gamma - 1.0;
    
    // Вычисление полной энергии
    double EL = pL / gamma_minus_1 + 0.5 * rhoL * uL * uL;
    double ER = pR / gamma_minus_1 + 0.5 * rhoR * uR * uR;
    
    // Векторы консервативных переменных
    std::array<double, 3> UL = {rhoL, rhoL * uL, EL};
    std::array<double, 3> UR = {rhoR, rhoR * uR, ER};
    
    // Физические потоки
    std::array<double, 3> FL = {
        rhoL * uL,
        rhoL * uL * uL + pL,
        uL * (EL + pL)
    };
    
    std::array<double, 3> FR = {
        rhoR * uR,
        rhoR * uR * uR + pR,
        uR * (ER + pR)
    };
    
    // Квадратура Гаусса с 3 точками
    constexpr int Q = 3;
    const double theta[Q] = {0.11270166537925831148, 0.5, 0.88729833462074168852};
    const double weight[Q] = {5.0/18.0, 4.0/9.0, 5.0/18.0};
    
    // Интеграл от |A| dU/d?
    std::array<double, 3> integral = {0.0, 0.0, 0.0};
    
    // Разности примитивных переменных
    double drho = rhoR - rhoL;
    double du = uR - uL;
    double dp = pR - pL;
    
    for (int k = 0; k < Q; ++k) {
        double psi = theta[k]; // Параметр вдоль пути (0 ? ? ? 1)
        
        // Линейный путь по примитивным переменным
        double rho = rhoL + psi * drho;
        double u = uL + psi * du;
        double p = pL + psi * dp;
        
        // Производные консервативных переменных по ?
        std::array<double, 3> dU_dpsi;
        dU_dpsi[0] = drho;
        dU_dpsi[1] = drho * u + rho * du;
        dU_dpsi[2] = dp / gamma_minus_1 + 0.5 * drho * u * u + rho * u * du;
        
        // Термодинамические величины
        double a = std::sqrt(gamma * p / std::max(rho, eps)); // Скорость звука
        double E = p / gamma_minus_1 + 0.5 * rho * u * u;
        double H = (E + p) / rho; // Энтальпия
        
        // Собственные значения (характеристические скорости)
        double lambda1 = u - a;
        double lambda2 = u;
        double lambda3 = u + a;
        
        // Левая матрица собственных векторов (строки)
        double beta = (gamma - 1.0) / (2.0 * a * a);
        
        // Строка 1 (соответствует ?1)
        double l11 = 0.5 * (beta * u * u + u / a);
        double l12 = -0.5 * (beta * u + 1.0 / a);
        double l13 = 0.5 * beta;
        
        // Строка 2 (соответствует ?2)
        double l21 = 1.0 - beta * u * u;
        double l22 = beta * u;
        double l23 = -beta;
        
        // Строка 3 (соответствует ?3)
        double l31 = 0.5 * (beta * u * u - u / a);
        double l32 = -0.5 * (beta * u - 1.0 / a);
        double l33 = 0.5 * beta;
        
        // Проекции на левые собственные векторы
        double alpha1 = l11 * dU_dpsi[0] + l12 * dU_dpsi[1] + l13 * dU_dpsi[2];
        double alpha2 = l21 * dU_dpsi[0] + l22 * dU_dpsi[1] + l23 * dU_dpsi[2];
        double alpha3 = l31 * dU_dpsi[0] + l32 * dU_dpsi[1] + l33 * dU_dpsi[2];
        
        // Умножение на абсолютные значения собственных чисел
        alpha1 *= std::abs(lambda1);
        alpha2 *= std::abs(lambda2);
        alpha3 *= std::abs(lambda3);
        
        // Правые собственные векторы (столбцы)
        // r1 = [1, u-a, H-u*a]^T
        // r2 = [1, u, 0.5*u^2]^T
        // r3 = [1, u+a, H+u*a]^T
        
        // Вклад в интеграл: R * |?| * L * dU/d?
        std::array<double, 3> term;
        term[0] = alpha1 + alpha2 + alpha3;
        term[1] = (u - a) * alpha1 + u * alpha2 + (u + a) * alpha3;
        term[2] = (H - u * a) * alpha1 + 0.5 * u * u * alpha2 + (H + u * a) * alpha3;
        
        // Добавление с весом квадратуры
        integral[0] += weight[k] * term[0];
        integral[1] += weight[k] * term[1];
        integral[2] += weight[k] * term[2];
    }
    
    // Численный поток по методу Ошера-Соломона
    std::array<double, 3> F_osher;
    for (int i = 0; i < 3; ++i) {
        F_osher[i] = 0.5 * (FL[i] + FR[i]) - 0.5 * integral[i];
    }
    
    return std::make_tuple(F_osher[0], F_osher[1], F_osher[2]);
}std::tuple<double,double,double> superSolve(double rhoL,double uL,double pL,
double rhoR,double uR,double pR,
std::string name){
 double xi = 0.0;
    
    if(name=="tocnoe"){
        auto p_star=solve_p_star(rhoL,uL,pL,rhoR,uR,pR);
        auto u_star_val=u_star(p_star,rhoL,uL,pL,rhoR,uR,pR);
        auto[rho_face,u_face,p_face]=sample(p_star,u_star_val,rhoL,uL,pL,rhoR,uR,pR,xi);
        return{rho_face,u_face,p_face};
    }
    else if(name=="roe"){
        auto[rho_face,u_face,p_face]=roe(rhoL,uL,pL,rhoR,uR,pR,xi);
        return{rho_face,u_face,p_face};
    }
    else if (name == "hll") {
        // Метод HLL
        auto [rho_face, u_face, p_face] = hll(rhoL, uL, pL, rhoR, uR, pR, xi);
        return {rho_face, u_face, p_face};
    }
    else if (name == "hllc") {
        // Метод HLLC
        auto [rho_face, u_face, p_face] = hllc(rhoL, uL, pL, rhoR, uR, pR, xi);
        return {rho_face, u_face, p_face};
    }
    else if (name == "rusanov") {
        auto [rho_face, u_face, p_face] = rusanov(rhoL, uL, pL, rhoR, uR, pR, xi);
        return {rho_face, u_face, p_face};
    }

    else if (name == "osher-solomon") {
        auto [rho_face, u_face, p_face] = osher_solomon(rhoL, uL, pL,rhoR, uR, pR, xi, 1.4);   
        return {rho_face, u_face, p_face};
    }

    
    return{0,0,0};
}
double minmod(double a, double b){
    if (a*a < a*b){return a;}
    else if (b*b < a*b){return b;}
    else if (a*b<0){return 0;}
    return 0;
}
std::tuple<std::vector<double>,std::vector<double>,std::vector<double>>
computeFluxes(int cells,int ghost_cells,
const std::vector<double>& rho,
const std::vector<double>& u,
const std::vector<double>& P,
double h,int Q,double gamma, std::string maneRS){
    
    int total_cells=cells+2*ghost_cells;
    std::vector<double> F1(total_cells,0.0);
    std::vector<double> F2(total_cells,0.0);
    std::vector<double> F3(total_cells,0.0);
    
    double x_over_t=0.0;
    
    for(int i=ghost_cells-1;i<cells+ghost_cells;i++){
        if(Q==0){
            auto[rho_face,u_face,p_face]=superSolve(rho[i],u[i],P[i],
            rho[i+1],u[i+1],P[i+1],maneRS);
         if(maneRS!="rusanov" && maneRS!="hll" &&  maneRS!="osher-solomon"){
            F1[i+1]=rho_face*u_face;
            F2[i+1]=rho_face*u_face*u_face+p_face;
            double E_face=p_face/(gamma-1)+0.5*rho_face*u_face*u_face;
            F3[i+1]=u_face*(E_face+p_face);}
	    else{
	 F1[i+1]=rho_face;
            F2[i+1]=u_face;
            double E_face=p_face/(gamma-1)+0.5*rho_face*u_face*u_face;
            F3[i+1]=p_face;
	 }
        }
        else if(Q==1){
            double drho_i=minmod((rho[i]-rho[i-1]),(rho[i+1]-rho[i]));
            double du_i=minmod((u[i]-u[i-1]),(u[i+1]-u[i]));
            double dP_i=minmod((P[i]-P[i-1]),(P[i+1]-P[i]));
            
            double drho_ip1=minmod((rho[i+1]-rho[i]),(rho[i+2]-rho[i+1]));
            double du_ip1=minmod((u[i+1]-u[i]),(u[i+2]-u[i+1]));
            double dP_ip1=minmod((P[i+1]-P[i]),(P[i+2]-P[i+1]));
            
            double rho_R_i=rho[i]+0.5*drho_i;
            double u_R_i=u[i]+0.5*du_i;
            double P_R_i=P[i]+0.5*dP_i;
            
            double rho_L_ip1=rho[i+1]-0.5*drho_ip1;
            double u_L_ip1=u[i+1]-0.5*du_ip1;
            double P_L_ip1=P[i+1]-0.5*dP_ip1;
            
            auto[rho_face,u_face,p_face]=superSolve(rho_R_i,u_R_i,P_R_i,
            rho_L_ip1,u_L_ip1,P_L_ip1,maneRS);
            
            if(maneRS!="rusanov"){
            F1[i+1]=rho_face*u_face;
            F2[i+1]=rho_face*u_face*u_face+p_face;
            double E_face=p_face/(gamma-1)+0.5*rho_face*u_face*u_face;
            F3[i+1]=u_face*(E_face+p_face);}
	    else{ F1[i+1]=rho_face;
            F2[i+1]=u_face;
            double E_face=p_face/(gamma-1)+0.5*rho_face*u_face*u_face;
            F3[i+1]=p_face;
	 }
        }
    }
    
    return{F1,F2,F3};
}

void GodunovSolve(int cells,int ghost_cells,
std::vector<double>& u,
std::vector<double>& P,
std::vector<double>& rho,
std::vector<double>& I,
double& t_total,double tau,double h,
const std::string& left_bc,
const std::string& right_bc,
int Q, const std::string& time_integrator, std::string nameRS){
        
    double gamma=1.4;
    int total_cells=cells+2*ghost_cells;
    
    std::vector<double> U1(total_cells);
    std::vector<double> U2(total_cells);
    std::vector<double> U3(total_cells);
    
    std::vector<double> U1_new(total_cells);
    std::vector<double> U2_new(total_cells);
    std::vector<double> U3_new(total_cells);
    
    for(int i=0;i<total_cells;i++){
        U1[i]=rho[i];
        U2[i]=rho[i]*u[i];
        U3[i]=P[i]/(gamma-1)+0.5*rho[i]*u[i]*u[i];
    }
    
    std::vector<double> U1_temp,U2_temp,U3_temp;
    std::vector<double> rho_temp,u_temp,P_temp;
    
    if(time_integrator=="rk2"){
        U1_temp=U1;
        U2_temp=U2;
        U3_temp=U3;
        rho_temp=rho;
        u_temp=u;
        P_temp=P;
    }
    
    if(time_integrator=="euler"){
        auto[F1,F2,F3]=computeFluxes(cells,ghost_cells,rho,u,P,h,Q,gamma,nameRS);
        
        for(int i=ghost_cells;i<cells+ghost_cells;i++){
            U1_new[i]=U1[i]+tau/h*(F1[i]-F1[i+1]);
            U2_new[i]=U2[i]+tau/h*(F2[i]-F2[i+1]);
            U3_new[i]=U3[i]+tau/h*(F3[i]-F3[i+1]);
        }
    }
    else{
        auto[F1,F2,F3]=computeFluxes(cells,ghost_cells,rho,u,P,h,Q,gamma,nameRS);
        
        for(int i=ghost_cells;i<cells+ghost_cells;i++){
            U1_new[i]=U1[i]+tau/h*(F1[i]-F1[i+1]);
            U2_new[i]=U2[i]+tau/h*(F2[i]-F2[i+1]);
            U3_new[i]=U3[i]+tau/h*(F3[i]-F3[i+1]);
        }
    }
    
    for(int i=ghost_cells;i<cells+ghost_cells;i++){
        rho[i]=U1_new[i];
        if(rho[i]<1e-5)rho[i]=1e-5;
        
        u[i]=U2_new[i]/rho[i];
        
        double e_in=U3_new[i]-0.5*rho[i]*u[i]*u[i];
        if(e_in<1e-6)e_in=1e-6;
        
        P[i]=e_in*(gamma-1);
        I[i]=U3_new[i];
    }
    
    applyBoundaryConditions(rho.data(),u.data(),P.data(),cells,ghost_cells,left_bc,right_bc);
    
    t_total+=tau;
    
    std::cout<<"t="<<t_total<<" tau="<<tau
    <<" method="<<time_integrator<<std::endl;
}

void newTimeStep(std::vector<double>& u,
               std::vector<double>& P,
               std::vector<double>& rho,
               double& tau, double h, double CFL)
{
    double gamma=1.4;
    std::vector<double> vc;
    for(int i=0; i<u.size(); i++){
    	vc.push_back(std::abs(u[i])+std::sqrt(gamma*P[i]/rho[i]));
    }
    auto maxVC=*std::max_element(vc.begin(), vc.end());
   	
   double tau1=((CFL*h/(maxVC))<tau) ? (CFL*h/(maxVC)) : tau;
   if((CFL*h/(maxVC))<tau){std::cout<<"time was changed"<<tau1<<std::endl;}
   tau=(tau1>1e-10) ? tau1 : 1e-10;
}

std::tuple<std::vector<double>,std::vector<double>,std::vector<double>,
std::vector<double>,std::vector<double>,std::vector<double>>
predictorStep(int cells,int ghost_cells,
const std::vector<double>& rho,
const std::vector<double>& u,
const std::vector<double>& P,
const std::vector<double>& U1,
const std::vector<double>& U2,
const std::vector<double>& U3,
double tau,double h,double gamma){
    
    int total_cells=cells+2*ghost_cells;
    
    std::vector<double> dU1(total_cells,0.0);
    std::vector<double> dU2(total_cells,0.0);
    std::vector<double> dU3(total_cells,0.0);
    
    std::vector<double> rhoHalf(total_cells,0.0);
    std::vector<double> uHalf(total_cells,0.0);
    std::vector<double> pHalf(total_cells,0.0);
    
    for(int i=ghost_cells;i<cells+ghost_cells;i++){
        double drho_i=minmod((rho[i]-rho[i-1]),(rho[i+1]-rho[i]));
        double du_i=minmod((u[i]-u[i-1]),(u[i+1]-u[i]));
        double dP_i=minmod((P[i]-P[i-1]),(P[i+1]-P[i]));
        
        dU1[i]=drho_i/2;
        dU2[i]=du_i/2;
        dU3[i]=dP_i/2;
        
        double rho_R=rho[i]+0.5*drho_i;
        double u_R=u[i]+0.5*du_i;
        double P_R=P[i]+0.5*dP_i;
        
        double rho_L=rho[i]-0.5*drho_i;
        double u_L=u[i]-0.5*du_i;
        double P_L=P[i]-0.5*dP_i;
        
        double rho_star=U1[i]-tau/h*(rho_R*u_R-rho_L*u_L);
        double rhou_star=U2[i]-tau/h*((rho_R*u_R*u_R+P_R)-(rho_L*u_L*u_L+P_L));
        
        double I_L=P_L/(gamma-1)+0.5*rho_L*u_L*u_L;
        double I_R=P_R/(gamma-1)+0.5*rho_R*u_R*u_R;
        
        double F_I_L=rho_L*u_L*I_L;
        double F_I_R=rho_R*u_R*I_R;
        double I_star=U3[i]-tau/h*(F_I_R-F_I_L);
        
        rhoHalf[i]=0.5*(U1[i]+rho_star);
        uHalf[i]=0.5*(U2[i]+rhou_star)/rhoHalf[i];
        double E_half=0.5*(U3[i]+I_star);
        pHalf[i]=(gamma-1)*(E_half-0.5*rhoHalf[i]*uHalf[i]*uHalf[i]);
    }
    
    return{rhoHalf,uHalf,pHalf,dU1,dU2,dU3};
}

std::tuple<std::vector<double>,std::vector<double>,std::vector<double>>
correctorStep(int cells,int ghost_cells,
const std::vector<double>& rhoHalf,
const std::vector<double>& uHalf,
const std::vector<double>& pHalf,
const std::vector<double>& dU1,
const std::vector<double>& dU2,
const std::vector<double>& dU3,
double h,double gamma, std::string nameRS){
    
    int total_cells=cells+2*ghost_cells;
    std::vector<double> F1(total_cells+1,0.0);
    std::vector<double> F2(total_cells+1,0.0);
    std::vector<double> F3(total_cells+1,0.0);
    
    for(int i=ghost_cells-1;i<cells+ghost_cells;i++){
        double rho_R_i=rhoHalf[i]+dU1[i];
        double u_R_i=uHalf[i]+dU2[i];
        double P_R_i=pHalf[i]+dU3[i];
        
        double rho_L_ip1=rhoHalf[i+1]-dU1[i+1];
        double u_L_ip1=uHalf[i+1]-dU2[i+1];
        double P_L_ip1=pHalf[i+1]-dU3[i+1];
        

        auto[rho_face,u_face,p_face]=superSolve(rho_R_i,u_R_i,P_R_i,
        rho_L_ip1,u_L_ip1,P_L_ip1,nameRS);
        
 if(nameRS!="rusanov"){
            F1[i+1]=rho_face*u_face;
            F2[i+1]=rho_face*u_face*u_face+p_face;
            double E_face=p_face/(gamma-1)+0.5*rho_face*u_face*u_face;
            F3[i+1]=u_face*(E_face+p_face);}
	    else{
	 F1[i+1]=rho_face;
            F2[i+1]=u_face;
            double E_face=p_face/(gamma-1)+0.5*rho_face*u_face*u_face;
            F3[i+1]=p_face;
	 }
    }
    
    return{F1,F2,F3};
}

void RodionovSolve(int cells,int ghost_cells,
std::vector<double>& u,
std::vector<double>& P,
std::vector<double>& rho,
std::vector<double>& I,
double& t_total,double tau,double h,
const std::string& left_bc,
const std::string& right_bc, const std::string& time_integrator, std::string nameRS){

    
    double gamma=1.4;
    int total_cells=cells+2*ghost_cells;
    
    std::vector<double> U1(total_cells);
    std::vector<double> U2(total_cells);
    std::vector<double> U3(total_cells);
    
    for(int i=0;i<total_cells;i++){
        U1[i]=rho[i];
        U2[i]=rho[i]*u[i];
        U3[i]=P[i]/(gamma-1)+0.5*rho[i]*u[i]*u[i];
    }
    
    std::vector<double> U1_new(total_cells,0.0);
    std::vector<double> U2_new(total_cells,0.0);
    std::vector<double> U3_new(total_cells,0.0);
    
    if(time_integrator=="euler"){
        std::vector<double> dU1_temp(total_cells,0.0);
        std::vector<double> dU2_temp(total_cells,0.0);
        std::vector<double> dU3_temp(total_cells,0.0);
        
        for(int i=ghost_cells;i<cells+ghost_cells;i++){
            double drho_i=minmod((rho[i]-rho[i-1]),(rho[i+1]-rho[i]));
            double du_i=minmod((u[i]-u[i-1]),(u[i+1]-u[i]));
            double dP_i=minmod((P[i]-P[i-1]),(P[i+1]-P[i]));
            
            dU1_temp[i]=drho_i/2;
            dU2_temp[i]=du_i/2;
            dU3_temp[i]=dP_i/2;
        }
        
        auto[F1,F2,F3]=correctorStep(cells,ghost_cells,rho,u,P,
        dU1_temp,dU2_temp,dU3_temp,h,gamma, nameRS);
        
        for(int i=ghost_cells;i<cells+ghost_cells;i++){
            U1_new[i]=U1[i]+tau/h*(F1[i]-F1[i+1]);
            U2_new[i]=U2[i]+tau/h*(F2[i]-F2[i+1]);
            U3_new[i]=U3[i]+tau/h*(F3[i]-F3[i+1]);
        }
    }
    else if(time_integrator=="rk2"){
        auto[rhoHalf,uHalf,pHalf,dU1,dU2,dU3]=
        predictorStep(cells,ghost_cells,rho,u,P,U1,U2,U3,tau/2,h,gamma);
        
        applyBoundaryConditions(rhoHalf.data(),uHalf.data(),pHalf.data(),
        cells,ghost_cells,left_bc,right_bc);
        
        std::vector<double> U1_half(total_cells);
        std::vector<double> U2_half(total_cells);
        std::vector<double> U3_half(total_cells);
        
        for(int i=0;i<total_cells;i++){
            U1_half[i]=rhoHalf[i];
            U2_half[i]=rhoHalf[i]*uHalf[i];
            U3_half[i]=pHalf[i]/(gamma-1)+0.5*rhoHalf[i]*uHalf[i]*uHalf[i];
        }
        
        auto[F1_half,F2_half,F3_half]=correctorStep(cells,ghost_cells,rhoHalf,uHalf,pHalf,
        dU1,dU2,dU3,h,gamma, nameRS);
        
        for(int i=ghost_cells;i<cells+ghost_cells;i++){
            U1_new[i]=U1[i]+tau/h*(F1_half[i]-F1_half[i+1]);
            U2_new[i]=U2[i]+tau/h*(F2_half[i]-F2_half[i+1]);
            U3_new[i]=U3[i]+tau/h*(F3_half[i]-F3_half[i+1]);
        }
    }
    else{
        auto[rhoHalf,uHalf,pHalf,dU1,dU2,dU3]=
        predictorStep(cells,ghost_cells,rho,u,P,U1,U2,U3,tau,h,gamma);
        
        applyBoundaryConditions(rhoHalf.data(),uHalf.data(),pHalf.data(),
        cells,ghost_cells,left_bc,right_bc);
        
        auto[F1,F2,F3]=correctorStep(cells,ghost_cells,rhoHalf,uHalf,pHalf,
        dU1,dU2,dU3,h,gamma, nameRS);
        
        for(int i=ghost_cells;i<cells+ghost_cells;i++){
            U1_new[i]=U1[i]+tau/h*(F1[i]-F1[i+1]);
            U2_new[i]=U2[i]+tau/h*(F2[i]-F2[i+1]);
            U3_new[i]=U3[i]+tau/h*(F3[i]-F3[i+1]);
        }
    }
    
    for(int i=ghost_cells;i<cells+ghost_cells;i++){
        rho[i]=U1_new[i];
        if(rho[i]<1e-5)rho[i]=1e-5;
        
        u[i]=U2_new[i]/rho[i];
        
        double e_in=U3_new[i]-0.5*rho[i]*u[i]*u[i];
        if(e_in<1e-6)e_in=1e-6;
        
        P[i]=e_in*(gamma-1);
        I[i]=U3_new[i];
    }
    
    applyBoundaryConditions(rho.data(),u.data(),P.data(),
    cells,ghost_cells,left_bc,right_bc);
    
    t_total+=tau;
    
    std::cout<<"t="<<t_total<<" tau="<<tau
    <<" method="<<time_integrator<<std::endl;
}

void makeGIFfile(int cells, int ghost_cells, int g, double current_time, 
                 int save_step_interval, double save_time_interval, 
                 double h, std::vector<double>& u, std::vector<double>& P,
                 std::vector<double>& rho, std::vector<double>& I, 
                 double t_end, const std::string& way, bool use_time_interval = false) {
    
    bool should_save = false;
    
    if (use_time_interval) {
        // Сохранение по временному интервалу
        static double last_saved_time = 0.0;
        if (current_time - last_saved_time >= save_time_interval || g == 0 || current_time - last_saved_time==0) {
            should_save = true;
            last_saved_time = current_time;
        }
    } else {
        // Сохранение по шагам
        if (g % save_step_interval == 0) {
            should_save = true;
        }
    }
    
    if (should_save) {
        std::vector<double> rhoP(cells);
        std::vector<double> uP(cells);
        std::vector<double> pP(cells);
        std::vector<double> IP(cells);
        std::vector<double> grid_P(cells);

        for(int i = 0; i < cells; i++) {
            grid_P[i] = (h * i);
            rhoP[i] = rho[i + ghost_cells];
            pP[i] = P[i + ghost_cells];
            uP[i] = u[i + ghost_cells];
            IP[i] = (I[i + ghost_cells] - 0.5 * rhoP[i] * uP[i] * uP[i]) / (rhoP[i]);
        }
        
        // Добавляем время в имя файла для уникальности
        std::string filename = way + std::to_string(int(current_time * 1000))  + ".csv";
        outputCSV(filename, grid_P, uP, pP, rhoP, IP, current_time);
        
        std::cout << "Saved: " << filename << std::endl;
    }
}

std::tuple<double,double,double> ENOreconstruction_right(int i,
std::vector<double>& u,
std::vector<double>& P,
std::vector<double>& rho){
    
    double Pe,rhoe,ue;

    double Su_L=(u[i]-2*u[i-1]+u[i-2]);
    double Su_M=(u[i+1]-2*u[i]+u[i-1]);
    double Su_R=(u[i+2]-2*u[i+1]+u[i]);

    if(std::abs(Su_L)<std::abs(Su_M)&&std::abs(Su_L)<std::abs(Su_R)){
        ue=(1.0/3.0)*u[i-2]-(7.0/6.0)*u[i-1]+(11.0/6.0)*u[i];
    }
    else if(std::abs(Su_R)<std::abs(Su_L)&&std::abs(Su_R)<std::abs(Su_M)){
        ue=(1.0/3.0)*u[i]+(5.0/6.0)*u[i+1]-(1.0/6.0)*u[i+2];
    }
    else{
        ue=(-1.0/6.0)*u[i-1]+(5.0/6.0)*u[i]+(1.0/3.0)*u[i+1];
    }

    double Sp_L=(P[i]-2*P[i-1]+P[i-2]);
    double Sp_M=(P[i+1]-2*P[i]+P[i-1]);
    double Sp_R=(P[i+2]-2*P[i+1]+P[i]);

    if(std::abs(Sp_L)<std::abs(Sp_M)&&std::abs(Sp_L)<std::abs(Sp_R)){
        Pe=(1.0/3.0)*P[i-2]-(7.0/6.0)*P[i-1]+(11.0/6.0)*P[i];
    }
    else if(std::abs(Sp_R)<std::abs(Sp_L)&&std::abs(Sp_R)<std::abs(Sp_M)){
        Pe=(1.0/3.0)*P[i]+(5.0/6.0)*P[i+1]-(1.0/6.0)*P[i+2];
    }
    else{
        Pe=(-1.0/6.0)*P[i-1]+(5.0/6.0)*P[i]+(1.0/3.0)*P[i+1];
    }
    if(Pe<1e-5)Pe=1e-5;

    double Srho_L=(rho[i]-2*rho[i-1]+rho[i-2]);
    double Srho_M=(rho[i+1]-2*rho[i]+rho[i-1]);
    double Srho_R=(rho[i+2]-2*rho[i+1]+rho[i]);
    
    if(std::abs(Srho_L)<std::abs(Srho_M)&&std::abs(Srho_L)<std::abs(Srho_R)){
        rhoe=(1.0/3.0)*rho[i-2]-(7.0/6.0)*rho[i-1]+(11.0/6.0)*rho[i];
    }
    else if(std::abs(Srho_R)<std::abs(Srho_L)&&std::abs(Srho_R)<std::abs(Srho_M)){
        rhoe=(1.0/3.0)*rho[i]+(5.0/6.0)*rho[i+1]-(1.0/6.0)*rho[i+2];
    }
    else{
        rhoe=(-1.0/6.0)*rho[i-1]+(5.0/6.0)*rho[i]+(1.0/3.0)*rho[i+1];
    }
    if(rhoe<1e-4)rhoe=1e-4;

    return std::make_tuple(ue,Pe,rhoe);
}

std::tuple<double,double,double> ENOreconstruction_left(int i,
std::vector<double>& u,
std::vector<double>& P,
std::vector<double>& rho){
    
    double Pe,rhoe,ue;

    double Su_L=(u[i]-2*u[i-1]+u[i-2]);
    double Su_M=(u[i+1]-2*u[i]+u[i-1]);
    double Su_R=(u[i+2]-2*u[i+1]+u[i]);

    if(std::abs(Su_L)<std::abs(Su_M)&&std::abs(Su_L)<std::abs(Su_R)){
        ue=(-1.0/6.0)*u[i-2]+(5.0/6.0)*u[i-1]+(1.0/3.0)*u[i];
    }
    else if(std::abs(Su_R)<std::abs(Su_L)&&std::abs(Su_R)<std::abs(Su_M)){
        ue=(11.0/6.0)*u[i]-(7.0/6.0)*u[i+1]+(1.0/3.0)*u[i+2];
    }
    else{
        ue=(1.0/3.0)*u[i-1]+(5.0/6.0)*u[i]-(1.0/6.0)*u[i+1];
    }

    double Sp_L=(P[i]-2*P[i-1]+P[i-2]);
    double Sp_M=(P[i+1]-2*P[i]+P[i-1]);
    double Sp_R=(P[i+2]-2*P[i+1]+P[i]);

    if(std::abs(Sp_L)<std::abs(Sp_M)&&std::abs(Sp_L)<std::abs(Sp_R)){
        Pe=(-1.0/6.0)*P[i-2]+(5.0/6.0)*P[i-1]+(1.0/3.0)*P[i];
    }
    else if(std::abs(Sp_R)<std::abs(Sp_L)&&std::abs(Sp_R)<std::abs(Sp_M)){
        Pe=(11.0/6.0)*P[i]-(7.0/6.0)*P[i+1]+(1.0/3.0)*P[i+2];
    }
    else{
        Pe=(1.0/3.0)*P[i-1]+(5.0/6.0)*P[i]-(1.0/6.0)*P[i+1];
    }
    if(Pe<1e-5)Pe=1e-5;

    double Srho_L=(rho[i]-2*rho[i-1]+rho[i-2]);
    double Srho_M=(rho[i+1]-2*rho[i]+rho[i-1]);
    double Srho_R=(rho[i+2]-2*rho[i+1]+rho[i]);
    
    if(std::abs(Srho_L)<std::abs(Srho_M)&&std::abs(Srho_L)<std::abs(Srho_R)){
        rhoe=(-1.0/6.0)*rho[i-2]+(5.0/6.0)*rho[i-1]+(1.0/3.0)*rho[i];
    }
    else if(std::abs(Srho_R)<std::abs(Srho_L)&&std::abs(Srho_R)<std::abs(Srho_M)){
        rhoe=(11.0/6.0)*rho[i]-(7.0/6.0)*rho[i+1]+(1.0/3.0)*rho[i+2];
    }
    else{
        rhoe=(1.0/3.0)*rho[i-1]+(5.0/6.0)*rho[i]-(1.0/6.0)*rho[i+1];
    }
    if(rhoe<1e-4)rhoe=1e-4;

    return std::make_tuple(ue,Pe,rhoe);
}

std::tuple<double,double,double> compute_flux(double rho_face,double u_face,double p_face,double gamma){
    if(rho_face<1e-5)rho_face=1e-5;
    if(p_face<1e-5)p_face=1e-5;
    
    double F1=rho_face*u_face;
    double F2=rho_face*u_face*u_face+p_face;
    double E_face=p_face/(gamma-1)+0.5*rho_face*u_face*u_face;
    double F3=u_face*(E_face+p_face);
    
    return std::make_tuple(F1,F2,F3);
}

std::tuple<std::vector<double>,std::vector<double>,std::vector<double>> 
computeRHS_ENO(int cells,int ghost_cells,
std::vector<double>& u,std::vector<double>& P,std::vector<double>& rho,
double h,double gamma, std::string nameRS){
    
    int total_cells=cells+2*ghost_cells;
    std::vector<double> R1(total_cells,0.0);
    std::vector<double> R2(total_cells,0.0);
    std::vector<double> R3(total_cells,0.0);
    
    double x_over_t=0.0;

    for(int i=ghost_cells;i<cells+ghost_cells;i++){
        auto[ui,pi,rhoi]=ENOreconstruction_right(i,u,P,rho);
        auto[uip,pip,rhoip]=ENOreconstruction_left(i+1,u,P,rho);
        
        auto[rho_face,u_face,p_face]=superSolve(rhoi,ui,pi,rhoip,uip,pip, nameRS);
        
        auto[F1,F2,F3]=compute_flux(rho_face,u_face,p_face,gamma);

        auto[uim,pim,rhoim]=ENOreconstruction_right(i-1,u,P,rho);
        auto[uim1,pim1,rhoim1]=ENOreconstruction_left(i,u,P,rho);
        
        auto[rho_face1,u_face1,p_face1]=superSolve(rhoim,uim,pim,rhoim1,uim1,pim1, nameRS);
        
        auto[F11,F21,F31]=compute_flux(rho_face1,u_face1,p_face1,gamma);

        R1[i]=(F11-F1)/h;
        R2[i]=(F21-F2)/h;
        R3[i]=(F31-F3)/h;
    }

    return std::make_tuple(R1,R2,R3);
}

void updatePrimitiveVariables(int cells,int ghost_cells,
std::vector<double>& u,
std::vector<double>& P,
std::vector<double>& rho,
std::vector<double>& I,
const std::vector<double>& U1,
const std::vector<double>& U2,
const std::vector<double>& U3,
double gamma){
    
    for(int i=ghost_cells;i<cells+ghost_cells;i++){
        rho[i]=U1[i];
        if(rho[i]<1e-5)rho[i]=1e-5;
        
        u[i]=U2[i]/rho[i];
        
        double e_in=U3[i]-0.5*rho[i]*u[i]*u[i];
        if(e_in<1e-6)e_in=1e-6;
        
        P[i]=e_in*(gamma-1);
        I[i]=U3[i];
    }
}

void ENO(int cells,int ghost_cells,
std::vector<double>& u,
std::vector<double>& P,
std::vector<double>& rho,
std::vector<double>& I,
double& t_total,double tau,double h,
const std::string& left_bc,
const std::string& right_bc, const std::string& time_integrator, std::string nameRS){
       
    double gamma=1.4;
    int total_cells=cells+2*ghost_cells;
    
    std::vector<double> U1(total_cells);
    std::vector<double> U2(total_cells);
    std::vector<double> U3(total_cells);
    
    for(int i=0;i<total_cells;i++){
        U1[i]=rho[i];
        U2[i]=rho[i]*u[i];
        U3[i]=P[i]/(gamma-1)+0.5*rho[i]*u[i]*u[i];
    }
    
    std::vector<double> U1_new(total_cells);
    std::vector<double> U2_new(total_cells);
    std::vector<double> U3_new(total_cells);
    
    if(time_integrator=="euler"){
        auto[R1,R2,R3]=computeRHS_ENO(cells,ghost_cells,u,P,rho,h,gamma, nameRS);
        
        for(int i=ghost_cells;i<cells+ghost_cells;i++){
            U1_new[i]=U1[i]+tau*R1[i];
            U2_new[i]=U2[i]+tau*R2[i];
            U3_new[i]=U3[i]+tau*R3[i];
        }
        
        updatePrimitiveVariables(cells,ghost_cells,u,P,rho,I,U1_new,U2_new,U3_new,gamma);
        
        applyBoundaryConditions(rho.data(),u.data(),P.data(),
        cells,ghost_cells,left_bc,right_bc);
    }
    else if(time_integrator=="rk2"){
        auto[R1_1,R2_1,R3_1]=computeRHS_ENO(cells,ghost_cells,u,P,rho,h,gamma, nameRS);
        
        std::vector<double> U1_temp(total_cells);
        std::vector<double> U2_temp(total_cells);
        std::vector<double> U3_temp(total_cells);
        
        for(int i=ghost_cells;i<cells+ghost_cells;i++){
            U1_temp[i]=U1[i]+tau*R1_1[i];
            U2_temp[i]=U2[i]+tau*R2_1[i];
            U3_temp[i]=U3[i]+tau*R3_1[i];
        }
        
        std::vector<double> u_temp=u;
        std::vector<double> P_temp=P;
        std::vector<double> rho_temp=rho;
        
        updatePrimitiveVariables(cells,ghost_cells,u_temp,P_temp,rho_temp,I,
        U1_temp,U2_temp,U3_temp,gamma);
        
        applyBoundaryConditions(rho_temp.data(),u_temp.data(),P_temp.data(),
        cells,ghost_cells,left_bc,right_bc);
        
        auto[R1_2,R2_2,R3_2]=computeRHS_ENO(cells,ghost_cells,u_temp,P_temp,rho_temp,h,gamma, nameRS);
        
        for(int i=ghost_cells;i<cells+ghost_cells;i++){
            U1_new[i]=U1[i]+0.5*tau*(R1_1[i]+R1_2[i]);
            U2_new[i]=U2[i]+0.5*tau*(R2_1[i]+R2_2[i]);
            U3_new[i]=U3[i]+0.5*tau*(R3_1[i]+R3_2[i]);
            
        }        
        updatePrimitiveVariables(cells,ghost_cells,u,P,rho,I,U1_new,U2_new,U3_new,gamma);
        
        applyBoundaryConditions(rho.data(),u.data(),P.data(),
        cells,ghost_cells,left_bc,right_bc);
    }
    else{
        
        auto[R1_1,R2_1,R3_1]=computeRHS_ENO(cells,ghost_cells,u,P,rho,h,gamma, nameRS);
        
        std::vector<double> U1_temp(total_cells);
        std::vector<double> U2_temp(total_cells);
        std::vector<double> U3_temp(total_cells);
        
        for(int i=ghost_cells;i<cells+ghost_cells;i++){
            U1_temp[i]=U1[i]+tau*R1_1[i];
            U2_temp[i]=U2[i]+tau*R2_1[i];
            U3_temp[i]=U3[i]+tau*R3_1[i];
        }
        
        std::vector<double> u_temp=u;
        std::vector<double> P_temp=P;
        std::vector<double> rho_temp=rho;
        
        updatePrimitiveVariables(cells,ghost_cells,u_temp,P_temp,rho_temp,I,
        U1_temp,U2_temp,U3_temp,gamma);
        
        applyBoundaryConditions(rho_temp.data(),u_temp.data(),P_temp.data(),
        cells,ghost_cells,left_bc,right_bc);
        
        auto[R1_2,R2_2,R3_2]=computeRHS_ENO(cells,ghost_cells,u_temp,P_temp,rho_temp,h,gamma, nameRS);
        
        for(int i=ghost_cells;i<cells+ghost_cells;i++){
            U1_new[i]=U1[i]+0.5*tau*(R1_1[i]+R1_2[i]);
            U2_new[i]=U2[i]+0.5*tau*(R2_1[i]+R2_2[i]);
            U3_new[i]=U3[i]+0.5*tau*(R3_1[i]+R3_2[i]);
        }
        
        updatePrimitiveVariables(cells,ghost_cells,u,P,rho,I,U1_new,U2_new,U3_new,gamma);
        applyBoundaryConditions(rho.data(),u.data(),P.data(),
        cells,ghost_cells,left_bc,right_bc);
    }
    
    t_total+=tau;
    
    std::cout<<"t="<<t_total<<" tau="<<tau
    <<" method="<<time_integrator<<std::endl;
}
// CHECK: RECONSTRUCTION
std::tuple<double,double,double> WENO5reconstruction_right(int i,
std::vector<double>& u,std::vector<double>& P,std::vector<double>& rho){
    
    double Pe,rhoe,ue;
    double c0=0.3,c1=0.6,c2=0.1;
    double epsilon=1e-6;

    double u0=(2.0/6.0)*u[i-2]-(7.0/6.0)*u[i-1]+(11.0/6.0)*u[i];
    double u1=(-1.0/6.0)*u[i-1]+(5.0/6.0)*u[i]+(2.0/6.0)*u[i+1];
    double u2=(2.0/6.0)*u[i]+(5.0/6.0)*u[i+1]-(1.0/6.0)*u[i+2];

    double s0_u=(13.0/12.0)*pow(u[i-2]-2*u[i-1]+u[i],2)+
    0.25*pow(u[i-2]-4*u[i-1]+3*u[i],2);
    double s1_u=(13.0/12.0)*pow(u[i-1]-2*u[i]+u[i+1],2)+
    0.25*pow(u[i-1]-u[i+1],2);
    double s2_u=(13.0/12.0)*pow(u[i]-2*u[i+1]+u[i+2],2)+
    0.25*pow(3*u[i]-4*u[i+1]+u[i+2],2);

    double alpha0_u=c0/pow(epsilon+s0_u,2);
    double alpha1_u=c1/pow(epsilon+s1_u,2);
    double alpha2_u=c2/pow(epsilon+s2_u,2);
    double sum_alpha_u=alpha0_u+alpha1_u+alpha2_u;
    double w0_u=alpha0_u/sum_alpha_u;
    double w1_u=alpha1_u/sum_alpha_u;
    double w2_u=alpha2_u/sum_alpha_u;

    ue=w0_u*u0+w1_u*u1+w2_u*u2;

    double P0=(2.0/6.0)*P[i-2]-(7.0/6.0)*P[i-1]+(11.0/6.0)*P[i];
    double P1=(-1.0/6.0)*P[i-1]+(5.0/6.0)*P[i]+(2.0/6.0)*P[i+1];
    double P2=(2.0/6.0)*P[i]+(5.0/6.0)*P[i+1]-(1.0/6.0)*P[i+2];

    double s0_P=(13.0/12.0)*pow(P[i-2]-2*P[i-1]+P[i],2)+
    0.25*pow(P[i-2]-4*P[i-1]+3*P[i],2);
    double s1_P=(13.0/12.0)*pow(P[i-1]-2*P[i]+P[i+1],2)+
    0.25*pow(P[i-1]-P[i+1],2);
    double s2_P=(13.0/12.0)*pow(P[i]-2*P[i+1]+P[i+2],2)+
    0.25*pow(3*P[i]-4*P[i+1]+P[i+2],2);

    double alpha0_P=c0/pow(epsilon+s0_P,2);
    double alpha1_P=c1/pow(epsilon+s1_P,2);
    double alpha2_P=c2/pow(epsilon+s2_P,2);
    double sum_alpha_P=alpha0_P+alpha1_P+alpha2_P;
    double w0_P=alpha0_P/sum_alpha_P;
    double w1_P=alpha1_P/sum_alpha_P;
    double w2_P=alpha2_P/sum_alpha_P;

    Pe=w0_P*P0+w1_P*P1+w2_P*P2;
    if(Pe<1e-5)Pe=1e-5;

    double rho0=(2.0/6.0)*rho[i-2]-(7.0/6.0)*rho[i-1]+(11.0/6.0)*rho[i];
    double rho1=(-1.0/6.0)*rho[i-1]+(5.0/6.0)*rho[i]+(2.0/6.0)*rho[i+1];
    double rho2=(2.0/6.0)*rho[i]+(5.0/6.0)*rho[i+1]-(1.0/6.0)*rho[i+2];

    double s0_rho=(13.0/12.0)*pow(rho[i-2]-2*rho[i-1]+rho[i],2)+
    0.25*pow(rho[i-2]-4*rho[i-1]+3*rho[i],2);
    double s1_rho=(13.0/12.0)*pow(rho[i-1]-2*rho[i]+rho[i+1],2)+
    0.25*pow(rho[i-1]-rho[i+1],2);
    double s2_rho=(13.0/12.0)*pow(rho[i]-2*rho[i+1]+rho[i+2],2)+
    0.25*pow(3*rho[i]-4*rho[i+1]+rho[i+2],2);

    double alpha0_rho=c0/pow(epsilon+s0_rho,2);
    double alpha1_rho=c1/pow(epsilon+s1_rho,2);
    double alpha2_rho=c2/pow(epsilon+s2_rho,2);
    double sum_alpha_rho=alpha0_rho+alpha1_rho+alpha2_rho;
    double w0_rho=alpha0_rho/sum_alpha_rho;
    double w1_rho=alpha1_rho/sum_alpha_rho;
    double w2_rho=alpha2_rho/sum_alpha_rho;

    rhoe=w0_rho*rho0+w1_rho*rho1+w2_rho*rho2;
    if(rhoe<1e-5)rhoe=1e-5;

    return std::make_tuple(ue,Pe,rhoe);
}

std::tuple<double,double,double> WENO5reconstruction_left(int i,
std::vector<double>& u,std::vector<double>& P,std::vector<double>& rho){
    
    double Pe,rhoe,ue;
    double epsilon=1e-6;
    double c0=0.1,c1=0.6,c2=0.3;

    double u0=(-1.0/6.0)*u[i-2]+(5.0/6.0)*u[i-1]+(1.0/3.0)*u[i];
    double u1=(1.0/3.0)*u[i-1]+(5.0/6.0)*u[i]-(1.0/6.0)*u[i+1];
    double u2=(11.0/6.0)*u[i]-(7.0/6.0)*u[i+1]+(1.0/3.0)*u[i+2];

    double s0_u=(13.0/12.0)*pow(u[i-2]-2*u[i-1]+u[i],2)+
    0.25*pow(u[i-2]-4*u[i-1]+3*u[i],2);
    double s1_u=(13.0/12.0)*pow(u[i-1]-2*u[i]+u[i+1],2)+
    0.25*pow(u[i-1]-u[i+1],2);
    double s2_u=(13.0/12.0)*pow(u[i]-2*u[i+1]+u[i+2],2)+
    0.25*pow(3*u[i]-4*u[i+1]+u[i+2],2);

    double alpha0_u=c0/pow(epsilon+s0_u,2);
    double alpha1_u=c1/pow(epsilon+s1_u,2);
    double alpha2_u=c2/pow(epsilon+s2_u,2);
    double sum_alpha_u=alpha0_u+alpha1_u+alpha2_u;
    double w0_u=alpha0_u/sum_alpha_u;
    double w1_u=alpha1_u/sum_alpha_u;
    double w2_u=alpha2_u/sum_alpha_u;

    ue=w0_u*u0+w1_u*u1+w2_u*u2;

    double P0=(-1.0/6.0)*P[i-2]+(5.0/6.0)*P[i-1]+(1.0/3.0)*P[i];
    double P1=(1.0/3.0)*P[i-1]+(5.0/6.0)*P[i]-(1.0/6.0)*P[i+1];
    double P2=(11.0/6.0)*P[i]-(7.0/6.0)*P[i+1]+(1.0/3.0)*P[i+2];

    double s0_P=(13.0/12.0)*pow(P[i-2]-2*P[i-1]+P[i],2)+
    0.25*pow(P[i-2]-4*P[i-1]+3*P[i],2);
    double s1_P=(13.0/12.0)*pow(P[i-1]-2*P[i]+P[i+1],2)+
    0.25*pow(P[i-1]-P[i+1],2);
    double s2_P=(13.0/12.0)*pow(P[i]-2*P[i+1]+P[i+2],2)+
    0.25*pow(3*P[i]-4*P[i+1]+P[i+2],2);

    double alpha0_P=c0/pow(epsilon+s0_P,2);
    double alpha1_P=c1/pow(epsilon+s1_P,2);
    double alpha2_P=c2/pow(epsilon+s2_P,2);
    double sum_alpha_P=alpha0_P+alpha1_P+alpha2_P;
    double w0_P=alpha0_P/sum_alpha_P;
    double w1_P=alpha1_P/sum_alpha_P;
    double w2_P=alpha2_P/sum_alpha_P;

    Pe=w0_P*P0+w1_P*P1+w2_P*P2;
    if(Pe<1e-5)Pe=1e-5;

    double rho0=(-1.0/6.0)*rho[i-2]+(5.0/6.0)*rho[i-1]+(1.0/3.0)*rho[i];
    double rho1=(1.0/3.0)*rho[i-1]+(5.0/6.0)*rho[i]-(1.0/6.0)*rho[i+1];
    double rho2=(11.0/6.0)*rho[i]-(7.0/6.0)*rho[i+1]+(1.0/3.0)*rho[i+2];

    double s0_rho=(13.0/12.0)*pow(rho[i-2]-2*rho[i-1]+rho[i],2)+
    0.25*pow(rho[i-2]-4*rho[i-1]+3*rho[i],2);
    double s1_rho=(13.0/12.0)*pow(rho[i-1]-2*rho[i]+rho[i+1],2)+
    0.25*pow(rho[i-1]-rho[i+1],2);
    double s2_rho=(13.0/12.0)*pow(rho[i]-2*rho[i+1]+rho[i+2],2)+
    0.25*pow(3*rho[i]-4*rho[i+1]+rho[i+2],2);

    double alpha0_rho=c0/pow(epsilon+s0_rho,2);
    double alpha1_rho=c1/pow(epsilon+s1_rho,2);
    double alpha2_rho=c2/pow(epsilon+s2_rho,2);
    double sum_alpha_rho=alpha0_rho+alpha1_rho+alpha2_rho;
    double w0_rho=alpha0_rho/sum_alpha_rho;
    double w1_rho=alpha1_rho/sum_alpha_rho;
    double w2_rho=alpha2_rho/sum_alpha_rho;

    rhoe=w0_rho*rho0+w1_rho*rho1+w2_rho*rho2;
    if(rhoe<1e-5)rhoe=1e-5;

    return std::make_tuple(ue,Pe,rhoe);
}

std::tuple<double,double,double> compute_flux_WENO(double rho_face,double u_face,
double p_face,double gamma){
    if(rho_face<1e-5)rho_face=1e-5;
    if(p_face<1e-5)p_face=1e-5;
    
    double F1=rho_face*u_face;
    double F2=rho_face*u_face*u_face+p_face;
    double E_face=p_face/(gamma-1)+0.5*rho_face*u_face*u_face;
    double F3=u_face*(E_face+p_face);
    
    return std::make_tuple(F1,F2,F3);
}

std::tuple<std::vector<double>,std::vector<double>,std::vector<double>> 
computeRHS_WENO(int cells,int ghost_cells,
std::vector<double>& u,std::vector<double>& P,std::vector<double>& rho,
double h,double gamma, std::string nameRS){
    
    int total_cells=cells+2*ghost_cells;
    std::vector<double> R1(total_cells,0.0);
    std::vector<double> R2(total_cells,0.0);
    std::vector<double> R3(total_cells,0.0);
    
    double x_over_t=0.0;

    for(int i=ghost_cells;i<cells+ghost_cells;i++){
        auto[ui,pi,rhoi]=WENO5reconstruction_right(i,u,P,rho);
        auto[uip,pip,rhoip]=WENO5reconstruction_left(i+1,u,P,rho);
        
        auto[rho_face,u_face,p_face]=superSolve(rhoi,ui,pi,rhoip,uip,pip, nameRS);
        
        auto[F1,F2,F3]=compute_flux_WENO(rho_face,u_face,p_face,gamma);

        auto[uim,pim,rhoim]=WENO5reconstruction_right(i-1,u,P,rho);
        auto[uim1,pim1,rhoim1]=WENO5reconstruction_left(i,u,P,rho);
        
        auto[rho_face1,u_face1,p_face1]=superSolve(rhoim,uim,pim,rhoim1,uim1,pim1,nameRS);
        
        auto[F11,F21,F31]=compute_flux_WENO(rho_face1,u_face1,p_face1,gamma);

        R1[i]=(F11-F1)/h;
        R2[i]=(F21-F2)/h;
        R3[i]=(F31-F3)/h;
    }

    return std::make_tuple(R1,R2,R3);
}

void updatePrimitiveVariables_WENO(int cells,int ghost_cells,
std::vector<double>& u,
std::vector<double>& P,
std::vector<double>& rho,
std::vector<double>& I,
const std::vector<double>& U1,
const std::vector<double>& U2,
const std::vector<double>& U3,
double gamma){
    
    for(int i=ghost_cells;i<cells+ghost_cells;i++){
        rho[i]=U1[i];
        if(rho[i]<1e-5)rho[i]=1e-5;
        
        u[i]=U2[i]/rho[i];
        
        double e_in=U3[i]-0.5*rho[i]*u[i]*u[i];
        if(e_in<1e-6)e_in=1e-6;
        
        P[i]=e_in*(gamma-1);
        I[i]=U3[i];
    }
}

void WENO(int cells,int ghost_cells,
std::vector<double>& u,
std::vector<double>& P,
std::vector<double>& rho,
std::vector<double>& I,
double& t_total,double tau,double h,
const std::string& left_bc,
const std::string& right_bc, const std::string& time_integrator, std::string nameRS){
       
    double gamma=1.4;
    int total_cells=cells+2*ghost_cells;
    
    std::vector<double> U1(total_cells);
    std::vector<double> U2(total_cells);
    std::vector<double> U3(total_cells);
    
    for(int i=0;i<total_cells;i++){
        U1[i]=rho[i];
        U2[i]=rho[i]*u[i];
        U3[i]=P[i]/(gamma-1)+0.5*rho[i]*u[i]*u[i];
    }
    
    std::vector<double> U1_new(total_cells);
    std::vector<double> U2_new(total_cells);
    std::vector<double> U3_new(total_cells);
    
    if(time_integrator=="euler"){
        auto[R1,R2,R3]=computeRHS_WENO(cells,ghost_cells,u,P,rho,h,gamma, nameRS);
        
        for(int i=ghost_cells;i<cells+ghost_cells;i++){
            U1_new[i]=U1[i]+tau*R1[i];
            U2_new[i]=U2[i]+tau*R2[i];
            U3_new[i]=U3[i]+tau*R3[i];
        }
        
        updatePrimitiveVariables_WENO(cells,ghost_cells,u,P,rho,I,U1_new,U2_new,U3_new,gamma);
        
        applyBoundaryConditions(rho.data(),u.data(),P.data(),
        cells,ghost_cells,left_bc,right_bc);
    }
    else if(time_integrator=="rk2"){
        auto[R1_1,R2_1,R3_1]=computeRHS_WENO(cells,ghost_cells,u,P,rho,h,gamma, nameRS);
        
        std::vector<double> U1_temp(total_cells);
        std::vector<double> U2_temp(total_cells);
        std::vector<double> U3_temp(total_cells);
        
        for(int i=ghost_cells;i<cells+ghost_cells;i++){
            U1_temp[i]=U1[i]+tau*R1_1[i];
            U2_temp[i]=U2[i]+tau*R2_1[i];
            U3_temp[i]=U3[i]+tau*R3_1[i];
        }
        
        std::vector<double> u_temp=u;
        std::vector<double> P_temp=P;
        std::vector<double> rho_temp=rho;
        
        updatePrimitiveVariables_WENO(cells,ghost_cells,u_temp,P_temp,rho_temp,I,
        U1_temp,U2_temp,U3_temp,gamma);
        
        applyBoundaryConditions(rho_temp.data(),u_temp.data(),P_temp.data(),
        cells,ghost_cells,left_bc,right_bc);
        
        auto[R1_2,R2_2,R3_2]=computeRHS_WENO(cells,ghost_cells,u_temp,P_temp,rho_temp,h,gamma, nameRS);
        
        for(int i=ghost_cells;i<cells+ghost_cells;i++){
                U1_new[i]=U1[i]+0.5*tau*(R1_1[i]+R1_2[i]);
                U2_new[i]=U2[i]+0.5*tau*(R2_1[i]+R2_2[i]);
                U3_new[i]=U3[i]+0.5*tau*(R3_1[i]+R3_2[i]);
            }
        
        
        updatePrimitiveVariables_WENO(cells,ghost_cells,u,P,rho,I,U1_new,U2_new,U3_new,gamma);
        
        applyBoundaryConditions(rho.data(),u.data(),P.data(),
        cells,ghost_cells,left_bc,right_bc);
    }
    else if(time_integrator=="rk3"){
        auto[R1_1,R2_1,R3_1]=computeRHS_WENO(cells,ghost_cells,u,P,rho,h,gamma, nameRS);
        
        std::vector<double> U1_temp1(total_cells);
        std::vector<double> U2_temp1(total_cells);
        std::vector<double> U3_temp1(total_cells);
        
        for(int i=ghost_cells;i<cells+ghost_cells;i++){
            U1_temp1[i]=U1[i]+tau*R1_1[i];
            U2_temp1[i]=U2[i]+tau*R2_1[i];
            U3_temp1[i]=U3[i]+tau*R3_1[i];
        }
        
        std::vector<double> u_temp1=u;
        std::vector<double> P_temp1=P;
        std::vector<double> rho_temp1=rho;
        
        updatePrimitiveVariables_WENO(cells,ghost_cells,u_temp1,P_temp1,rho_temp1,I,
        U1_temp1,U2_temp1,U3_temp1,gamma);
        
        applyBoundaryConditions(rho_temp1.data(),u_temp1.data(),P_temp1.data(),
        cells,ghost_cells,left_bc,right_bc);
        
        auto[R1_2,R2_2,R3_2]=computeRHS_WENO(cells,ghost_cells,u_temp1,P_temp1,rho_temp1,h,gamma, nameRS);
        
        std::vector<double> U1_temp2(total_cells);
        std::vector<double> U2_temp2(total_cells);
        std::vector<double> U3_temp2(total_cells);
        
        for(int i=ghost_cells;i<cells+ghost_cells;i++){
            U1_temp2[i]=0.75*U1[i]+0.25*(U1_temp1[i]+tau*R1_2[i]);
            U2_temp2[i]=0.75*U2[i]+0.25*(U2_temp1[i]+tau*R2_2[i]);
            U3_temp2[i]=0.75*U3[i]+0.25*(U3_temp1[i]+tau*R3_2[i]);
        }
        
        std::vector<double> u_temp2=u;
        std::vector<double> P_temp2=P;
        std::vector<double> rho_temp2=rho;
        
        updatePrimitiveVariables_WENO(cells,ghost_cells,u_temp2,P_temp2,rho_temp2,I,
        U1_temp2,U2_temp2,U3_temp2,gamma);
        
        applyBoundaryConditions(rho_temp2.data(),u_temp2.data(),P_temp2.data(),
        cells,ghost_cells,left_bc,right_bc);
        
        auto[R1_3,R2_3,R3_3]=computeRHS_WENO(cells,ghost_cells,u_temp2,P_temp2,rho_temp2,h,gamma, nameRS);
        
        for(int i=ghost_cells;i<cells+ghost_cells;i++){
            U1_new[i]=(1.0/3.0)*U1[i]+(2.0/3.0)*(U1_temp2[i]+tau*R1_3[i]);
            U2_new[i]=(1.0/3.0)*U2[i]+(2.0/3.0)*(U2_temp2[i]+tau*R2_3[i]);
            U3_new[i]=(1.0/3.0)*U3[i]+(2.0/3.0)*(U3_temp2[i]+tau*R3_3[i]);
        }
        
        updatePrimitiveVariables_WENO(cells,ghost_cells,u,P,rho,I,U1_new,U2_new,U3_new,gamma);
        
        applyBoundaryConditions(rho.data(),u.data(),P.data(),
        cells,ghost_cells,left_bc,right_bc);
    }
    else{
        std::cout<<"Warning: Unknown time integrator, using RK2 method\n";
        
        auto[R1_1,R2_1,R3_1]=computeRHS_WENO(cells,ghost_cells,u,P,rho,h,gamma, nameRS);
        
        std::vector<double> U1_temp(total_cells);
        std::vector<double> U2_temp(total_cells);
        std::vector<double> U3_temp(total_cells);
        
        for(int i=ghost_cells;i<cells+ghost_cells;i++){
            U1_temp[i]=U1[i]+tau*R1_1[i];
            U2_temp[i]=U2[i]+tau*R2_1[i];
            U3_temp[i]=U3[i]+tau*R3_1[i];
        }
        
        std::vector<double> u_temp=u;
        std::vector<double> P_temp=P;
        std::vector<double> rho_temp=rho;
        
        updatePrimitiveVariables_WENO(cells,ghost_cells,u_temp,P_temp,rho_temp,I,
        U1_temp,U2_temp,U3_temp,gamma);
        
        applyBoundaryConditions(rho_temp.data(),u_temp.data(),P_temp.data(),
        cells,ghost_cells,left_bc,right_bc);
        
        auto[R1_2,R2_2,R3_2]=computeRHS_WENO(cells,ghost_cells,u_temp,P_temp,rho_temp,h,gamma, nameRS);
        
        for(int i=ghost_cells;i<cells+ghost_cells;i++){
            U1_new[i]=U1[i]+0.5*tau*(R1_1[i]+R1_2[i]);
            U2_new[i]=U2[i]+0.5*tau*(R2_1[i]+R2_2[i]);
            U3_new[i]=U3[i]+0.5*tau*(R3_1[i]+R3_2[i]);
        }
        
        updatePrimitiveVariables_WENO(cells,ghost_cells,u,P,rho,I,U1_new,U2_new,U3_new,gamma);
        applyBoundaryConditions(rho.data(),u.data(),P.data(),
        cells,ghost_cells,left_bc,right_bc);
    }
    
    t_total+=tau;
    
    std::cout<<"t="<<t_total<<" tau="<<tau
    <<" method="<<time_integrator<<std::endl;
}
void  Fletcher(int cells,int ghost_cells,
               std::vector<double>& u,
               std::vector<double>& P,
               std::vector<double>& rho,
               std::vector<double>& I,
               double& t_total, double tau, double h, const std::string& left_bc, 
               const std::string& right_bc){

    double gamma=1.4;
	
    int total_cells=cells+2*ghost_cells;   
 
    std::vector<double> U1_newS(total_cells);
    std::vector<double> U2_newS(total_cells);
    std::vector<double> U3_newS(total_cells);
    
    std::vector<double> U1(total_cells);
    std::vector<double> U2(total_cells);
    std::vector<double> U3(total_cells);

    std::vector<double> U1_new(total_cells);
    std::vector<double> U2_new(total_cells);
    std::vector<double> U3_new(total_cells);

    std::vector<double> rhoS(total_cells);
    std::vector<double> US(total_cells);
    std::vector<double> PS(total_cells);
    std::vector<double> IS(total_cells);

    for(int i=0;i<total_cells;i++){
    	U1[i]=rho[i];
        U2[i]=rho[i]*u[i];
    	U3[i]=P[i]/(gamma-1)+0.5*rho[i]*u[i]*u[i];;
    }

    for(int i=ghost_cells;i<cells+ghost_cells;i++){
	
	double F1=rho[i]*u[i];
    	double F2=rho[i]*u[i]*u[i]+P[i];
    	double E=P[i]/(gamma-1)+0.5*rho[i]*u[i]*u[i];
    	double F3=u[i]*(E+P[i]);

    	double F11=rho[i+1]*u[i+1];
    	double F21=rho[i+1]*u[i+1]*u[i+1]+P[i+1];
    	double E1=P[i+1]/(gamma-1)+0.5*rho[i+1]*u[i+1]*u[i+1];
    	double F31=u[i+1]*(E1+P[i+1]);

	
	U1_newS[i]=U1[i]-tau/h*(F11-F1);
    	U2_newS[i]=U2[i]-tau/h*(F21-F2);
    	U3_newS[i]=U3[i]-tau/h*(F31-F3);
	}

    for(int i=ghost_cells;i<cells+ghost_cells;i++){
    	
	rhoS[i]=U1_newS[i];
    	if(rhoS[i]<1e-5)rhoS[i]=1e-5;
    
    	US[i]=U2_newS[i]/rhoS[i];
    
    	double e_in=U3_newS[i]-0.5*rhoS[i]*US[i]*US[i];
    	if(e_in<1e-6)e_in=1e-6;
    
   	PS[i]=e_in*(gamma-1);
    	IS[i]=U3_newS[i];
	}


	applyBoundaryConditions(rhoS.data(), US.data(), PS.data(),cells,ghost_cells,left_bc,right_bc);

	for(int i=ghost_cells;i<cells+ghost_cells;i++){

	double F1=rhoS[i-1]*US[i-1];
    	double F2=rhoS[i-1]*US[i-1]*US[i-1]+PS[i-1];
    	double E=PS[i-1]/(gamma-1)+0.5*rhoS[i-1]*US[i-1]*US[i-1];
    	double F3=US[i-1]*(E+PS[i-1]);

    	double F11=rhoS[i]*US[i];
    	double F21=rhoS[i]*US[i]*US[i]+PS[i];
    	double E1=PS[i]/(gamma-1)+0.5*rhoS[i]*US[i]*US[i];
    	double F31=US[i]*(E1+PS[i]);

	
	U1_new[i]=0.5*(U1_newS[i]+U1[i]-tau/h*(F11-F1));
    	U2_new[i]=0.5*(U2_newS[i]+U2[i]-tau/h*(F21-F2));
    	U3_new[i]=0.5*(U3_newS[i]+U3[i]-tau/h*(F31-F3));
	}
	
	double viscosity_coeff=0.1;
	std::vector<double> U1_viscous(total_cells), U2_viscous(total_cells), U3_viscous(total_cells);

	for(int i=ghost_cells; i<cells+ghost_cells; i++) {
    		if(i > ghost_cells && i<cells+ghost_cells-1) {
        		U1_viscous[i]=U1_new[i]+viscosity_coeff*(U1_new[i+1]-2*U1_new[i]+U1_new[i-1]);
        		U2_viscous[i]=U2_new[i]+viscosity_coeff*(U2_new[i+1]-2*U2_new[i]+U2_new[i-1]);
        		U3_viscous[i]=U3_new[i]+viscosity_coeff*(U3_new[i+1]-2*U3_new[i]+U3_new[i-1]);
    		} 
    		else {
        	U1_viscous[i]=U1_new[i];
        	U2_viscous[i]=U2_new[i];
        	U3_viscous[i]=U3_new[i];
    		}
	}

	double anti_diffusion_coeff=0.1;
	std::vector<double> U1_final(total_cells), U2_final(total_cells), U3_final(total_cells);

	auto applyLimitedAntiDiffusion=[&](const std::vector<double>& U_viscous, 
                                   std::vector<double>& U_final) {
    	int n=total_cells;
    	std::vector<double> phi(n, 0.0); // антидиффузионные потоки
    	std::vector<double> phi_limited(n, 0.0); // ограниченные потоки
    
        for(int i=ghost_cells; i<cells+ghost_cells-1; i++) {
        phi[i]=anti_diffusion_coeff*(U_viscous[i+1]-U_viscous[i]);
    	}
    
    
    	for(int i=ghost_cells+1; i<cells+ghost_cells-2; i++) {
        	double s=(phi[i]>=0) ? 1.0 : -1.0;
        
        	double delta_u_left=U_viscous[i]-U_viscous[i-1];   
        	double delta_u_right=U_viscous[i+2]-U_viscous[i+1]; 
        
        	double abs_phi=std::abs(phi[i]);
        	double term1=s*delta_u_left;
        	double term2=s*delta_u_right;
        
        	double min_val=std::max(0.0, std::min({term1, abs_phi, term2}));
        	phi_limited[i]=s *min_val;
    	}
    
    	for (int i=ghost_cells; i<cells+ghost_cells; i++) {
        	if (i>ghost_cells && i<cells+ghost_cells-1) {
            		U_final[i]=U_viscous[i]-(phi_limited[i]-phi_limited[i-1]);
        }
		else {
            	U_final[i]=U_viscous[i];
        	}
    }
	};


	applyLimitedAntiDiffusion(U1_viscous, U1_final);
	applyLimitedAntiDiffusion(U2_viscous, U2_final);
	applyLimitedAntiDiffusion(U3_viscous, U3_final);
	for(int i = ghost_cells; i < cells + ghost_cells; i++) {
    		rho[i] = U1_final[i];
    	u[i] = U2_final[i] / rho[i];
    	double e_in = U3_final[i] - 0.5 * rho[i] * u[i] * u[i];
    	P[i] = e_in * (gamma-1);
    	if(P[i] < 1e-8) P[i] = 1e-8;
    	I[i] = U3_final[i];
	}
    
    applyBoundaryConditions(rho.data(),u.data(),P.data(),cells,ghost_cells,left_bc,right_bc);
    std::cout<<"t = "<<t_total<<" tau = "<<tau<<std::endl;
    t_total+=tau;
	     
}

void makeNormAnaliz(int count, int N_all, int start) {
double tau, h, cfl, t_output;
    int cells, ghost_cells=3;
    std::string left_bc, right_bc, test_name, t_end_user, key;
    std::vector<std::string> solver_names;
    double rho_L, u_L, p_L, rho_R, u_R, p_R, t_end,  x_over_t ;

    if (readInit("start.txt", tau, h, cells, cfl, left_bc, right_bc, test_name, count, solver_names, t_end_user, t_output, key)) {
        std::cout << "Mesh is read" << std::endl;
        
        if (test_name == "custom") {
            readTest("custom_test.txt", rho_L, u_L, p_L, rho_R, u_R, p_R, t_end);
        } else {
            setupSODTest(test_name, rho_L, u_L, p_L, rho_R, u_R, p_R, t_end, t_end_user);
        }

        double gamma = 1.4;
        std::vector<double> grid_P;
        std::vector<double> errRho, errU, errP, errI;

        for(int N = start; N < N_all; N = N + count) {
            h = 1.0 / N;
            int current_cells = N;
            int total_cells = current_cells + 2 * ghost_cells;
            
            std::cout << "Calculating for N = " << N << ", h = " << h << std::endl;

            grid_P.push_back(h);

            std::vector<double> rho(total_cells);
            std::vector<double> u(total_cells);
            std::vector<double> p(total_cells);
            std::vector<double> I(total_cells);

            int break_point = ghost_cells + current_cells / 2;

            for (int i = 0; i < total_cells; i++) {
                if (i < break_point) {
                    rho[i] = rho_L;
                    u[i] = u_L;
                    p[i] = p_L;
                    I[i] = p[i] / ((gamma - 1) * rho[i]);
                } else {
                    rho[i] = rho_R;
                    u[i] = u_R;
                    p[i] = p_R;
                    I[i] = p[i] / ((gamma - 1) * rho[i]);
                }
            }

            applyBoundaryConditions(rho.data(), u.data(), p.data(), 
                                  current_cells, ghost_cells, left_bc, right_bc);

            double t_total = 0.0;
            while (t_total < t_end) {
                newTimeStep(u, p, rho, tau, h, cfl);
                WENO(current_cells, ghost_cells, u, p, rho, 
                        I, t_total, tau, h, left_bc, right_bc, "rk2", "tocnoe");
                
                if (tau <= 0) {
                    std::cerr << "Error: non-positive tau detected!" << std::endl;
                    break;
                }
                
                t_total += tau;
                
                if (t_total > 10 * t_end) {
                    std::cerr << "Warning: t_total exceeds 10*t_end, breaking loop" << std::endl;
                    break;
                }
            }

            std::vector<double> rhoP(current_cells);
            std::vector<double> uP(current_cells);
            std::vector<double> pP(current_cells);
            std::vector<double> IP(current_cells);
            
            for (int i = 0; i < current_cells; i++) {
                rhoP[i] = rho[i + ghost_cells];
                pP[i] = p[i + ghost_cells];
                uP[i] = u[i + ghost_cells];
                IP[i] = (I[i + ghost_cells] - 0.5 * rhoP[i] * uP[i] * uP[i]) / rhoP[i];
            }

 
            double p_star = solve_p_star(rho_L, u_L, p_L, rho_R, u_R, p_R);
            double u_star1 = u_star(p_star, rho_L, u_L, p_L, rho_R, u_R, p_R);
            
            std::vector<double> rhoAbs, Uabs, Pabs, Iabs;

            for (int i = 0; i < current_cells; i++) {
                x_over_t = ((i + 0.5) * h - 0.5) / t_end;  // Корректный расчет координаты
                auto [rho_analit, u_analit, p_analit] = sample(p_star, u_star1, rho_L, u_L, p_L, 
                                                              rho_R, u_R, p_R, x_over_t);
                double I_analit = p_analit / (rho_analit * (gamma - 1));
                
                rhoAbs.push_back(std::abs(rhoP[i] - rho_analit));
                Uabs.push_back(std::abs(uP[i] - u_analit));
                Pabs.push_back(std::abs(pP[i] - p_analit));
                Iabs.push_back(std::abs(IP[i] - I_analit));
            }
            errRho.push_back(*std::max_element(rhoAbs.begin(), rhoAbs.end()));
            errU.push_back(*std::max_element(Uabs.begin(), Uabs.end()));
            errP.push_back(*std::max_element(Pabs.begin(), Pabs.end()));
            errI.push_back(*std::max_element(Iabs.begin(), Iabs.end()));

            std::cout << "N = " << N << " completed. Errors - rho: " << errRho.back() 
                      << ", u: " << errU.back() << ", p: " << errP.back() << std::endl;
        }

        // Сохранение результатов
        outputCSV("data/errL1.csv", grid_P, errU, errP, errRho, errI, t_end);
        std::cout << "Analysis completed. Results saved to data/errL1.csv" << std::endl;
        
        if (grid_P.size() > 1) {
            std::cout << "\nConvergence rates:" << std::endl;
            for (size_t i = 1; i < grid_P.size(); i++) {
                double rate_rho = log(errRho[i-1] / errRho[i]) / log(grid_P[i-1] / grid_P[i]);
                double rate_u = log(errU[i-1] / errU[i]) / log(grid_P[i-1] / grid_P[i]);
                double rate_p = log(errP[i-1] / errP[i]) / log(grid_P[i-1] / grid_P[i]);
                
                std::cout << "h = " << grid_P[i] << " -> rho: " << rate_rho 
                          << ", u: " << rate_u << ", p: " << rate_p << std::endl;
            }
        }
    } else {
        std::cerr << "Error reading initial parameters!" << std::endl;
    }
}



// Superbee limiter
double superbee(double r){
    if(r<=0.0)return 0.0;
    if(r<=0.5)return 2.0*r;
    if(r<=1.0)return 1.0;
    if(r<=2.0)return r;
    return 2.0;
}

// Van Leer limiter
double vanLeer(double r){
    if(r<=0.0)return 0.0;
    return(r+std::abs(r))/(1.0+std::abs(r));
}

// MC limiter (monotonized central)
double mcLimiter(double r){
    if(r<=0.0)return 0.0;
    double result=std::min(2.0*r,0.5*(1.0+r));
    return std::min(result,2.0);
}

double computeTVDLimiter(double r,const std::string& limiter_type="minmod"){
    if(limiter_type=="minmod"){
        if(r<=0.0)return 0.0;
        return std::min(r,1.0);
    }
    else if(limiter_type=="superbee"){
        return superbee(r);
    }
    else if(limiter_type=="vanleer"){
        return vanLeer(r);
    }
    else if(limiter_type=="mc"){
        return mcLimiter(r);
    }
    else if(limiter_type=="none"){
        return 1.0;  // Без лимитера
    }
    else{
        // По умолчанию minmod
        if(r<=0.0)return 0.0;
        return std::min(r,1.0);
    }
}

std::tuple<double,double,double> computeLimiters(int i,
const std::vector<double>& rho,
const std::vector<double>& u,
const std::vector<double>& P,
const std::string& limiter_type="minmod"){
    
    double epsilon=1e-8;
    
    // For rho
    double drho_up=rho[i+1]-rho[i];
    double drho_down=rho[i]-rho[i-1];
    double r_rho=(std::abs(drho_down)>epsilon)?drho_up/drho_down:0.0;
    
    // For u
    double du_up=u[i+1]-u[i];
    double du_down=u[i]-u[i-1];
    double r_u=(std::abs(du_down)>epsilon)?du_up/du_down:0.0;
    
    // For P
    double dP_up=P[i+1]-P[i];
    double dP_down=P[i]-P[i-1];
    double r_P=(std::abs(dP_down)>epsilon)?dP_up/dP_down:0.0;
    
    double phi_rho=computeTVDLimiter(r_rho,limiter_type);
    double phi_u=computeTVDLimiter(r_u,limiter_type);
    double phi_P=computeTVDLimiter(r_P,limiter_type);
    
    return std::make_tuple(phi_rho,phi_u,phi_P);
}

void applyArtificialViscosity(int cells,int ghost_cells,
                             std::vector<double>& U1,std::vector<double>& U2,std::vector<double>& U3,
                             double viscosity_coeff,double h){
    int total_cells=cells+2*ghost_cells;
    std::vector<double> U1_old=U1;
    std::vector<double> U2_old=U2;
    std::vector<double> U3_old=U3;
    
    for(int i=ghost_cells;i<cells+ghost_cells;i++){
        if(i>ghost_cells&&i<cells+ghost_cells-1){
            U1[i]=U1_old[i]+viscosity_coeff*(U1_old[i+1]-2.0*U1_old[i]+U1_old[i-1]);
            U2[i]=U2_old[i]+viscosity_coeff*(U2_old[i+1]-2.0*U2_old[i]+U2_old[i-1]);
            U3[i]=U3_old[i]+viscosity_coeff*(U3_old[i+1]-2.0*U3_old[i]+U3_old[i-1]);
        }
    }
}

void applyAntiDiffusion(int cells,int ghost_cells,
                       std::vector<double>& U1,std::vector<double>& U2,std::vector<double>& U3,
                       double anti_diffusion_coeff){
    int total_cells=cells+2*ghost_cells;
    
    std::vector<double> U1_temp=U1;
    std::vector<double> U2_temp=U2;
    std::vector<double> U3_temp=U3;
    
       for(int i=ghost_cells+1;i<cells+ghost_cells-1;i++){
        U1[i]=U1_temp[i]+anti_diffusion_coeff*(U1_temp[i+1]-2.0*U1_temp[i]+U1_temp[i-1]);
        U2[i]=U2_temp[i]+anti_diffusion_coeff*(U2_temp[i+1]-2.0*U2_temp[i]+U2_temp[i-1]);
        U3[i]=U3_temp[i]+anti_diffusion_coeff*(U3_temp[i+1]-2.0*U3_temp[i]+U3_temp[i-1]);
    }
}

std::tuple<std::vector<double>,std::vector<double>,std::vector<double>> 
computeHighOrderFluxes(int cells,int ghost_cells,
                      const std::vector<double>& rho,
                      const std::vector<double>& u,
                      const std::vector<double>& P,
                      double h,double gamma,
                      const std::string& high_order_scheme="tvd",
                      const std::string& tvd_scheme="minmod", std::string nameRS="tocnoe"){
    
    int total_cells=cells+2*ghost_cells;
    std::vector<double> F1(total_cells+1,0.0);
    std::vector<double> F2(total_cells+1,0.0);
    std::vector<double> F3(total_cells+1,0.0);
    
    double x_over_t=0.0;
    
    if(high_order_scheme=="tvd"){
        std::vector<double> dU1_temp(total_cells,0.0);
        std::vector<double> dU2_temp(total_cells,0.0);
        std::vector<double> dU3_temp(total_cells,0.0);
        
        for(int i=ghost_cells;i<cells+ghost_cells;i++){
            auto[phi_rho,phi_u,phi_P]=computeLimiters(i,rho,u,P,tvd_scheme);
            
            double drho_L=rho[i]-rho[i-1];
            double drho_R=rho[i+1]-rho[i];
            
            double du_L=u[i]-u[i-1];
            double du_R=u[i+1]-u[i];
            
            double dP_L=P[i]-P[i-1];
            double dP_R=P[i+1]-P[i];
            
            if(tvd_scheme=="minmod"){
                dU1_temp[i]=0.5*minmod(drho_L,drho_R);
                dU2_temp[i]=0.5*minmod(du_L,du_R);
                dU3_temp[i]=0.5*minmod(dP_L,dP_R);
            }else{
                dU1_temp[i]=0.5*phi_rho*drho_L;
                dU2_temp[i]=0.5*phi_u*du_L;
                dU3_temp[i]=0.5*phi_P*dP_L;
            }
        }
        
        auto[F1_high,F2_high,F3_high]=correctorStep(cells,ghost_cells,rho,u,P,
                                                       dU1_temp,dU2_temp,dU3_temp,h,gamma, nameRS);
        return{F1_high,F2_high,F3_high};
    }
    else if(high_order_scheme=="eno"){
        for(int i=ghost_cells-1;i<cells+ghost_cells;i++){
            auto[ui,pi,rhoi]=ENOreconstruction_right(i,const_cast<std::vector<double>&>(u),
                                                         const_cast<std::vector<double>&>(P),
                                                         const_cast<std::vector<double>&>(rho));
            auto[uip,pip,rhoip]=ENOreconstruction_left(i+1,const_cast<std::vector<double>&>(u),
                                                          const_cast<std::vector<double>&>(P),
                                                          const_cast<std::vector<double>&>(rho));
            
            auto[rho_face,u_face,p_face]=superSolve(rhoi,ui,pi,rhoip,uip,pip, nameRS);
            
            F1[i+1]=rho_face*u_face;
            F2[i+1]=rho_face*u_face*u_face+p_face;
            double E_face=p_face/(gamma-1)+0.5*rho_face*u_face*u_face;
            F3[i+1]=u_face*(E_face+p_face);
        }
        return{F1,F2,F3};
    }
    else if(high_order_scheme=="weno"){
        for(int i=ghost_cells-1;i<cells+ghost_cells;i++){
            auto[ui,pi,rhoi]=WENO5reconstruction_right(i,const_cast<std::vector<double>&>(u),
                                                          const_cast<std::vector<double>&>(P),
                                                          const_cast<std::vector<double>&>(rho));
            auto[uip,pip,rhoip]=WENO5reconstruction_left(i+1,const_cast<std::vector<double>&>(u),
                                                           const_cast<std::vector<double>&>(P),
                                                           const_cast<std::vector<double>&>(rho));
            
            auto[rho_face,u_face,p_face]=superSolve(rhoi,ui,pi,rhoip,uip,pip, nameRS);
            
            F1[i+1]=rho_face*u_face;
            F2[i+1]=rho_face*u_face*u_face+p_face;
            double E_face=p_face/(gamma-1)+0.5*rho_face*u_face*u_face;
            F3[i+1]=u_face*(E_face+p_face);
        }
        return{F1,F2,F3};
    }
    else if(high_order_scheme=="rodi"){
        std::vector<double> dU1_temp(total_cells,0.0);
        std::vector<double> dU2_temp(total_cells,0.0);
        std::vector<double> dU3_temp(total_cells,0.0);
        
        for(int i=ghost_cells;i<cells+ghost_cells;i++){
            double drho_i=minmod((rho[i]-rho[i-1]),(rho[i+1]-rho[i]));
            double du_i=minmod((u[i]-u[i-1]),(u[i+1]-u[i]));
            double dP_i=minmod((P[i]-P[i-1]),(P[i+1]-P[i]));
            
            dU1_temp[i]=drho_i/2;
            dU2_temp[i]=du_i/2;
            dU3_temp[i]=dP_i/2;
        }
        
        auto[F1_high,F2_high,F3_high]=correctorStep(cells,ghost_cells,rho,u,P,
                                                       dU1_temp,dU2_temp,dU3_temp,h,gamma, nameRS);
        return{F1_high,F2_high,F3_high};
    }
    else{
        std::vector<double> dU1_temp(total_cells,0.0);
        std::vector<double> dU2_temp(total_cells,0.0);
        std::vector<double> dU3_temp(total_cells,0.0);
        
        for(int i=ghost_cells;i<cells+ghost_cells;i++){
            auto[phi_rho,phi_u,phi_P]=computeLimiters(i,rho,u,P,tvd_scheme);
            
            double drho_L=rho[i]-rho[i-1];
            double drho_R=rho[i+1]-rho[i];
            
            double du_L=u[i]-u[i-1];
            double du_R=u[i+1]-u[i];
            
            double dP_L=P[i]-P[i-1];
            double dP_R=P[i+1]-P[i];
            
            if(tvd_scheme=="minmod"){
                dU1_temp[i]=0.5*minmod(drho_L,drho_R);
                dU2_temp[i]=0.5*minmod(du_L,du_R);
                dU3_temp[i]=0.5*minmod(dP_L,dP_R);
            }else{
                dU1_temp[i]=0.5*phi_rho*drho_L;
                dU2_temp[i]=0.5*phi_u*du_L;
                dU3_temp[i]=0.5*phi_P*dP_L;
            }
        }
        
        auto[F1_high,F2_high,F3_high]=correctorStep(cells,ghost_cells,rho,u,P,
                                                       dU1_temp,dU2_temp,dU3_temp,h,gamma, "tocnoe");
        return{F1_high,F2_high,F3_high};
    }
}

std::tuple<std::vector<double>,std::vector<double>,std::vector<double>>
computeTVDFluxes(int cells,int ghost_cells,
                const std::vector<double>& rho,
                const std::vector<double>& u,
                const std::vector<double>& P,
                double h,double gamma,
                const std::string& high_order_scheme,
                const std::string& reconstruction,
                const std::string& tvd_scheme,
                const std::string& left_bc,
                const std::string& right_bc, std::string nameRS="tocnoe"){
    
    int total_cells=cells+2*ghost_cells;
    

    auto[F1_high,F2_high,F3_high]=computeHighOrderFluxes(cells,ghost_cells,rho,u,P,h,gamma,
                                                            high_order_scheme,tvd_scheme, nameRS);
    
    auto[F1_low,F2_low,F3_low]=computeFluxes(cells,ghost_cells,rho,u,P,h,0,gamma,"tocnoe");
    
    std::vector<double> F1_TVD(total_cells+1,0.0);
    std::vector<double> F2_TVD(total_cells+1,0.0);
    std::vector<double> F3_TVD(total_cells+1,0.0);
    

    for(int i=ghost_cells;i<=cells+ghost_cells;i++){
        double phi=1.0;
        
        if(reconstruction=="high"){

            if(i>ghost_cells&&i<cells+ghost_cells){
                double epsilon=1e-8;
                
                double delta_rho_left=rho[i]-rho[i-1];
                double delta_rho_right=rho[i-1]-rho[i-2];
                double r_rho=(std::abs(delta_rho_right)>epsilon)?
                              delta_rho_left/delta_rho_right:0.0;
                
                double delta_u_left=u[i]-u[i-1];
                double delta_u_right=u[i-1]-u[i-2];
                double r_u=(std::abs(delta_u_right)>epsilon)?
                            delta_u_left/delta_u_right:0.0;
                
                double delta_P_left=P[i]-P[i-1];
                double delta_P_right=P[i-1]-P[i-2];
                double r_P=(std::abs(delta_P_right)>epsilon)?
                            delta_P_left/delta_P_right:0.0;
                
                double phi_rho=computeTVDLimiter(r_rho,tvd_scheme);
                double phi_u=computeTVDLimiter(r_u,tvd_scheme);
                double phi_P=computeTVDLimiter(r_P,tvd_scheme);
                
                phi=std::min(std::min(phi_rho,phi_u),phi_P);
            }
        }else{
            // Чистая схема Годунова
            phi=0.0;
        }
        
        F1_TVD[i]=F1_low[i]+phi*(F1_high[i]-F1_low[i]);
        F2_TVD[i]=F2_low[i]+phi*(F2_high[i]-F2_low[i]);
        F3_TVD[i]=F3_low[i]+phi*(F3_high[i]-F3_low[i]);
    }
    
    return{F1_TVD,F2_TVD,F3_TVD};
}

std::tuple<std::vector<double>,std::vector<double>,std::vector<double>>
eulerStep(int cells,int ghost_cells,
          const std::vector<double>& U1,const std::vector<double>& U2,const std::vector<double>& U3,
          const std::vector<double>& F1_TVD,const std::vector<double>& F2_TVD,const std::vector<double>& F3_TVD,
          double tau,double h){
    
    int total_cells=cells+2*ghost_cells;
    std::vector<double> U1_new=U1;
    std::vector<double> U2_new=U2;
    std::vector<double> U3_new=U3;
    
     for(int i=ghost_cells;i<cells+ghost_cells;i++){
        U1_new[i]=U1[i]+tau/h*(F1_TVD[i]-F1_TVD[i+1]);
        U2_new[i]=U2[i]+tau/h*(F2_TVD[i]-F2_TVD[i+1]);
        U3_new[i]=U3[i]+tau/h*(F3_TVD[i]-F3_TVD[i+1]);
    }
    
    return{U1_new,U2_new,U3_new};
}

std::tuple<std::vector<double>,std::vector<double>,std::vector<double>>
rungeKutta2Step(int cells,int ghost_cells,
               const std::vector<double>& U1,const std::vector<double>& U2,const std::vector<double>& U3,
               const std::vector<double>& rho,const std::vector<double>& u,const std::vector<double>& P,
               double tau,double h,double gamma,
               const std::string& high_order_scheme,const std::string& reconstruction,
               const std::string& tvd_scheme,
               const std::string& left_bc,const std::string& right_bc, std::string nameRS){
    
    int total_cells=cells+2*ghost_cells;
    
    auto[F1_TVD1,F2_TVD1,F3_TVD1]=computeTVDFluxes(cells,ghost_cells,rho,u,P,h,gamma,
                                                      high_order_scheme,reconstruction,tvd_scheme,
                                                      left_bc,right_bc, nameRS);
    
    std::vector<double> k1_U1(total_cells,0.0);
    std::vector<double> k1_U2(total_cells,0.0);
    std::vector<double> k1_U3(total_cells,0.0);
    
    for(int i=ghost_cells;i<cells+ghost_cells;i++){
        k1_U1[i]=(F1_TVD1[i]-F1_TVD1[i+1])/h;
        k1_U2[i]=(F2_TVD1[i]-F2_TVD1[i+1])/h;
        k1_U3[i]=(F3_TVD1[i]-F3_TVD1[i+1])/h;
    }
    
    std::vector<double> U1_temp=U1;
    std::vector<double> U2_temp=U2;
    std::vector<double> U3_temp=U3;
    
    for(int i=ghost_cells;i<cells+ghost_cells;i++){
        U1_temp[i]=U1[i]+tau*k1_U1[i];
        U2_temp[i]=U2[i]+tau*k1_U2[i];
        U3_temp[i]=U3[i]+tau*k1_U3[i];
    }
    
    std::vector<double> rho_temp=rho;
    std::vector<double> u_temp=u;
    std::vector<double> P_temp=P;
    
    for(int i=0;i<total_cells;i++){
        rho_temp[i]=U1_temp[i];
        if(rho_temp[i]<1e-5)rho_temp[i]=1e-5;
        
        u_temp[i]=U2_temp[i]/rho_temp[i];
        
        double e_in=U3_temp[i]-0.5*rho_temp[i]*u_temp[i]*u_temp[i];
        if(e_in<1e-6)e_in=1e-6;
        
        P_temp[i]=e_in*(gamma-1);
    }
    
    applyBoundaryConditions(rho_temp.data(),u_temp.data(),P_temp.data(),
                           cells,ghost_cells,left_bc,right_bc);
    
    auto[F1_TVD2,F2_TVD2,F3_TVD2]=computeTVDFluxes(cells,ghost_cells,rho_temp,u_temp,P_temp,h,gamma,
                                                      high_order_scheme,reconstruction,tvd_scheme,
                                                      left_bc,right_bc, nameRS);
    
    std::vector<double> U1_new=U1;
    std::vector<double> U2_new=U2;
    std::vector<double> U3_new=U3;
    
    for(int i=ghost_cells;i<cells+ghost_cells;i++){
        double k2_U1=(F1_TVD2[i]-F1_TVD2[i+1])/h;
        double k2_U2=(F2_TVD2[i]-F2_TVD2[i+1])/h;
        double k2_U3=(F3_TVD2[i]-F3_TVD2[i+1])/h;
        
        U1_new[i]=U1[i]+tau*0.5*(k1_U1[i]+k2_U1);
        U2_new[i]=U2[i]+tau*0.5*(k1_U2[i]+k2_U2);
        U3_new[i]=U3[i]+tau*0.5*(k1_U3[i]+k2_U3);
    }
    
    return{U1_new,U2_new,U3_new};
}

void feoktistov(int cells,int ghost_cells,
                std::vector<double>& u,std::vector<double>& P,std::vector<double>& rho,
                std::vector<double>& I,double& t_total,double tau,double h,
                const std::string& left_bc,const std::string& right_bc,
                const std::string& high_order_scheme="tvd",
                const std::string& reconstruction="high",
                const std::string& tvd_scheme="minmod",
                bool use_diffusion=false,
                bool use_anti_diffusion=false,
                bool use_artificial_viscosity=false,
                double viscosity_coeff=0.01,
                double anti_diffusion_coeff=0.05,
                const std::string& time_integrator="euler", std::string nameRS="tocnoe")
{      
    double gamma=1.4;
    int total_cells=cells+2*ghost_cells;
    
    std::vector<double> U1(total_cells);
    std::vector<double> U2(total_cells);
    std::vector<double> U3(total_cells);
    
    for(int i=0;i<total_cells;i++){
        U1[i]=rho[i];
        U2[i]=rho[i]*u[i];
        U3[i]=P[i]/(gamma-1)+0.5*rho[i]*u[i]*u[i];
    }
    
    std::vector<double> U1_new(total_cells,0.0);
    std::vector<double> U2_new(total_cells,0.0);
    std::vector<double> U3_new(total_cells,0.0);
    
     if(time_integrator=="euler"){
        // Метод Эйлера (явный)
        auto[F1_TVD,F2_TVD,F3_TVD]=computeTVDFluxes(cells,ghost_cells,rho,u,P,h,gamma,
                                                       high_order_scheme,reconstruction,tvd_scheme,
                                                       left_bc,right_bc, nameRS);
        
        std::tie(U1_new,U2_new,U3_new)=eulerStep(cells,ghost_cells,U1,U2,U3,
                                                    F1_TVD,F2_TVD,F3_TVD,tau,h);
    }
    else if(time_integrator=="rk2"||time_integrator=="runge-kutta"){
        // Метод Рунге-Кутты 2-го порядка
        std::tie(U1_new,U2_new,U3_new)=rungeKutta2Step(cells,ghost_cells,U1,U2,U3,
                                                          rho,u,P,tau,h,gamma,
                                                          high_order_scheme,reconstruction,tvd_scheme,
                                                          left_bc,right_bc, nameRS);
    }
    else{
          std::cout<<"Warning: Unknown time integrator '"<<time_integrator
                  <<"', using Euler method."<<std::endl;
        
        auto[F1_TVD,F2_TVD,F3_TVD]=computeTVDFluxes(cells,ghost_cells,rho,u,P,h,gamma,
                                                       high_order_scheme,reconstruction,tvd_scheme,
                                                       left_bc,right_bc, nameRS);
        
        std::tie(U1_new,U2_new,U3_new)=eulerStep(cells,ghost_cells,U1,U2,U3,
                                                    F1_TVD,F2_TVD,F3_TVD,tau,h);
    }
    
    // Применение искусственной вязкости
    if(use_artificial_viscosity){
        applyArtificialViscosity(cells,ghost_cells,U1_new,U2_new,U3_new,viscosity_coeff,h);
    }
    
    // Применение диффузии
    if(use_diffusion){
  
        for(int i=ghost_cells+1;i<cells+ghost_cells-1;i++){
            double diff_coeff=anti_diffusion_coeff;
            U1_new[i]+=diff_coeff*(U1_new[i+1]-2.0*U1_new[i]+U1_new[i-1]);
            U2_new[i]+=diff_coeff*(U2_new[i+1]-2.0*U2_new[i]+U2_new[i-1]);
            U3_new[i]+=diff_coeff*(U3_new[i+1]-2.0*U3_new[i]+U3_new[i-1]);
        }
    }
    

    if(use_anti_diffusion){
        applyAntiDiffusion(cells,ghost_cells,U1_new,U2_new,U3_new,anti_diffusion_coeff);
    }
    
   
    for(int i=ghost_cells;i<cells+ghost_cells;i++){
        rho[i]=U1_new[i];
        if(rho[i]<1e-5)rho[i]=1e-5;
        
        u[i]=U2_new[i]/rho[i];
        
        double e_in=U3_new[i]-0.5*rho[i]*u[i]*u[i];
        if(e_in<1e-6)e_in=1e-6;
        
        P[i]=e_in*(gamma-1);
        I[i]=U3_new[i];
    }
    
    // Антидиффузионная коррекция для плотности
    for(int i=ghost_cells+1;i<cells+ghost_cells-1;i++){
        double drho_left=rho[i]-rho[i-1];
        double drho_right=rho[i+1]-rho[i];
        
        if(drho_left*drho_right<0.0){
            double rho_avg=0.5*(rho[i-1]+rho[i+1]);
            double limit_factor=0.5;
            
            rho[i]=rho[i]+limit_factor*(rho_avg-rho[i]);
            if(rho[i]<1e-5)rho[i]=1e-5;
            
            U1_new[i]=rho[i];
            U2_new[i]=rho[i]*u[i];
            U3_new[i]=P[i]/(gamma-1)+0.5*rho[i]*u[i]*u[i];
        }
    }
    

    applyBoundaryConditions(rho.data(),u.data(),P.data(),
                           cells,ghost_cells,left_bc,right_bc);
    
    t_total+=tau;
    
    std::cout<<"t="<<t_total<<" tau="<<tau
              <<" time_integrator="<<time_integrator
              <<" high_order="<<high_order_scheme
              <<" TVD="<<tvd_scheme
              <<" reconstruction="<<reconstruction
              <<" viscosity="<<(use_artificial_viscosity?"on":"off")
              <<" anti-diffusion="<<(use_anti_diffusion?"on":"off")<<std::endl;
}

#include <fstream>
#include <string>
#include <sstream>
#include <algorithm>
#include <cctype>
#include <iostream>

void trim(std::string& s) {
    s.erase(s.begin(), std::find_if(s.begin(), s.end(), [](unsigned char ch) {
        return !std::isspace(ch);
    }));
    s.erase(std::find_if(s.rbegin(), s.rend(), [](unsigned char ch) {
        return !std::isspace(ch);
    }).base(), s.end());
}

void readStartFileF(
    std::string& nameRS,
    std::string& time_integrator,
    std::string& high_shem,
    std::string& tipe_reconsrutin,
    std::string& limiter,
    bool& use_diffusion,
    bool& use_anti_diffusion,
    bool& use_artificial_viscosity,
    double& viscosity_coeff,
    double& anti_diffusion_coeff
) {
    // Значения по умолчанию
    nameRS = "tocnoe";
    time_integrator = "rk2";
    high_shem = "weno";
    tipe_reconsrutin = "high";
    limiter = "minmod";
    use_diffusion = false;
    use_anti_diffusion = false;
    use_artificial_viscosity = false;
    viscosity_coeff = 0.001;
    anti_diffusion_coeff = 0.001;
    
    std::ifstream file("start.txt");
    if (!file.is_open()) {
        std::cerr << "Не удалось открыть файл start.txt" << std::endl;
        return;
    }
    
    std::string line;
    while (std::getline(file, line)) {
        // Удаляем комментарии (все что после #)
        size_t comment_pos = line.find('#');
        if (comment_pos != std::string::npos) {
            line = line.substr(0, comment_pos);
        }
        
        trim(line);
        if (line.empty()) continue;
        
        size_t equals_pos = line.find('=');
        if (equals_pos == std::string::npos) continue;
        
        std::string key = line.substr(0, equals_pos);
        std::string value = line.substr(equals_pos + 1);
        
        trim(key);
        trim(value);
        
        // Удаляем кавычки, если они есть
        if (!value.empty() && value.front() == '"' && value.back() == '"') {
            value = value.substr(1, value.size() - 2);
        }
        
        // Преобразуем ключ к нижнему регистру для сравнения
        std::string key_lower = key;
        std::transform(key_lower.begin(), key_lower.end(), key_lower.begin(),
                      [](unsigned char c){ return std::tolower(c); });
        
        if (key_lower == "namers") {
            nameRS = value;
        } else if (key_lower == "time_integrator") {
            time_integrator = value;
        } else if (key_lower == "high_shem") {
            high_shem = value;
        } else if (key_lower == "tipe_reconsrutin") {
            tipe_reconsrutin = value;
        } else if (key_lower == "limiter") {
            limiter = value;
        } else if (key_lower == "use_diffusion") {
            std::string val_lower = value;
            std::transform(val_lower.begin(), val_lower.end(), val_lower.begin(),
                          [](unsigned char c){ return std::tolower(c); });
            use_diffusion = (val_lower == "true" || val_lower == "1" || val_lower == "yes");
        } else if (key_lower == "use_anti_diffusion") {
            std::string val_lower = value;
            std::transform(val_lower.begin(), val_lower.end(), val_lower.begin(),
                          [](unsigned char c){ return std::tolower(c); });
            use_anti_diffusion = (val_lower == "true" || val_lower == "1" || val_lower == "yes");
        } else if (key_lower == "use_artificial_viscosity") {
            std::string val_lower = value;
            std::transform(val_lower.begin(), val_lower.end(), val_lower.begin(),
                          [](unsigned char c){ return std::tolower(c); });
            use_artificial_viscosity = (val_lower == "true" || val_lower == "1" || val_lower == "yes");
        } else if (key_lower == "viscosity_coeff") {
            try {
                viscosity_coeff = std::stod(value);
            } catch (...) {
                std::cerr << "Ошибка чтения viscosity_coeff: " << value << std::endl;
            }
        } else if (key_lower == "anti_diffusion_coeff") {
            try {
                anti_diffusion_coeff = std::stod(value);
            } catch (...) {
                std::cerr << "Ошибка чтения anti_diffusion_coeff: " << value << std::endl;
            }
        }
    }
    
    file.close();
}