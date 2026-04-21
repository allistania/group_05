// ============================================================
// Основной файл: godunov_gmsh.cpp
// Компиляция: g++ -o godunov godunov_gmsh.cpp gmsh_io.cpp -lm
// ============================================================

#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <map>
#include <tuple>
#include <string>
#include <cassert>

// Подключаем библиотеку для чтения MSH 2.2
#include "gmsh_io.hpp"

using namespace std;

// ============================================================
// Константы
// ============================================================
const double GAMMA = 1.4;          // показатель адиабаты
const double CFL   = 0.4;          // число Куранта
const double T_END = 1.0;          // конечное время
const int    SAVE_INTERVAL = 100;   // сохранять каждые 100 итераций

// ============================================================
// Структуры данных
// ============================================================

// Консервативные переменные
struct Conserved {
    double rho;   // плотность
    double rhou;  // x-импульс
    double rhov;  // y-импульс
    double rhoE;  // полная энергия
};

// Примитивные переменные
struct Primitive {
    double rho, u, v, p;
};

// Ячейка (треугольник)
struct Cell {
    int id;                     // индекс (0..ncells-1)
    vector<int> nodes;          // глобальные номера узлов (1-индексация)
    double area;                // площадь
    double xc, yc;              // центр масс
    Conserved U;                // консервативные переменные
    Conserved dU;               // приращение за шаг
    int phys_tag;               // физический тег из MSH (для подобластей)
};

// Грань (ребро между ячейками или граница)
struct Face {
    int id;                     // индекс
    int cellL, cellR;           // индексы ячеек (-1 для границы)
    int type;                   // 0 – внутренняя, 1 – граничная
    int phys_tag;               // физический тег границы (если type==1)
    vector<int> nodes;          // два узла (глобальные номера)
    double nx, ny;              // единичная нормаль, направленная из cellL наружу
    double length;              // длина грани (для 2D)
};

// ============================================================
// Глобальные данные сетки
// ============================================================
vector<double> nodeCoords;      // координаты узлов: 3 * nodeNum
int nodeNum, elemNum;           // общее количество узлов и элементов (всех)
vector<int> elemNodes;          // связность для всех элементов (4*elemNum)
vector<int> elemTags;           // физические теги для всех элементов

vector<Cell> cells;             // только треугольные ячейки (volume)
vector<Face> faces;             // все грани

// ============================================================
// Функции для работы с примитивными/консервативными переменными
// ============================================================
Primitive getPrimitive(const Conserved& U) {
    double rho = U.rho;
    double u = U.rhou / rho;
    double v = U.rhov / rho;
    double e = U.rhoE / rho - 0.5 * (u*u + v*v);
    double p = (GAMMA - 1.0) * rho * e;
    return {rho, u, v, p};
}

Conserved getConserved(const Primitive& prim) {
    double rho = prim.rho;
    double u = prim.u;
    double v = prim.v;
    double p = prim.p;
    double e = p / (GAMMA - 1.0) / rho;
    double E = e + 0.5 * (u*u + v*v);
    return {rho, rho*u, rho*v, rho*E};
}

// ============================================================
// Геометрические вычисления для треугольника
// ============================================================
double triangleArea(const vector<int>& nodes, const vector<double>& coords) {
    double x1 = coords[3*(nodes[0]-1)];
    double y1 = coords[3*(nodes[0]-1)+1];
    double x2 = coords[3*(nodes[1]-1)];
    double y2 = coords[3*(nodes[1]-1)+1];
    double x3 = coords[3*(nodes[2]-1)];
    double y3 = coords[3*(nodes[2]-1)+1];
    return 0.5 * fabs((x2-x1)*(y3-y1) - (x3-x1)*(y2-y1));
}

void triangleCentroid(const vector<int>& nodes, const vector<double>& coords,
                      double& xc, double& yc) {
    double x1 = coords[3*(nodes[0]-1)];
    double y1 = coords[3*(nodes[0]-1)+1];
    double x2 = coords[3*(nodes[1]-1)];
    double y2 = coords[3*(nodes[1]-1)+1];
    double x3 = coords[3*(nodes[2]-1)];
    double y3 = coords[3*(nodes[2]-1)+1];
    xc = (x1 + x2 + x3) / 3.0;
    yc = (y1 + y2 + y3) / 3.0;
}

// ============================================================
// Чтение сетки из MSH 2.2 (через gmsh_io)
// ============================================================
bool readMesh(const string& filename) {
    // Читаем размеры
    gmsh_size_read(filename, nodeNum, elemNum);
    nodeCoords.resize(3 * nodeNum);
    elemNodes.resize(4 * elemNum);   // 4 числа на элемент (3 узла + -1)
    elemTags.resize(elemNum);
    gmsh_data_read(filename, nodeNum, elemNum,
                   nodeCoords.data(), elemNodes.data(), elemTags.data());
    return true;
}

// ============================================================
// Построение ячеек (треугольников)
// ============================================================
void buildCells() {
    cells.clear();
    // Проходим по всем элементам, выбираем треугольники (4-й номер узла == -1)
    for (int i = 0; i < elemNum; ++i) {
        int n3 = elemNodes[4*i+3];
        if (n3 == -1) { // треугольник
            Cell c;
            c.id = cells.size();
            c.nodes = { elemNodes[4*i], elemNodes[4*i+1], elemNodes[4*i+2] };
            c.phys_tag = elemTags[i];
            c.area = triangleArea(c.nodes, nodeCoords);
            triangleCentroid(c.nodes, nodeCoords, c.xc, c.yc);
            // Начальные условия – пока нулевые, заполним позже
            c.U = {0,0,0,0};
            c.dU = {0,0,0,0};
            cells.push_back(c);
        }
    }
    cout << "Построено ячеек: " << cells.size() << endl;
}

// ============================================================
// Построение граней (внутренних и граничных)
// ============================================================
void buildFaces() {
    // Словарь: ребро (упорядоченные номера узлов) -> список индексов ячеек
    map<pair<int,int>, vector<int>> edgeToCells;

    // Заполняем для каждой ячейки её рёбра
    for (size_t cid = 0; cid < cells.size(); ++cid) {
        const auto& cell = cells[cid];
        // Три ребра: (n0,n1), (n1,n2), (n2,n0)
        int n[3] = {cell.nodes[0], cell.nodes[1], cell.nodes[2]};
        for (int i = 0; i < 3; ++i) {
            int a = n[i];
            int b = n[(i+1)%3];
            if (a > b) swap(a,b);
            edgeToCells[{a,b}].push_back(cid);
        }
    }

    // Теперь создаём грани
    faces.clear();
    int faceId = 0;
    for (auto& kv : edgeToCells) {
        const auto& edge = kv.first;
        const auto& cellIndices = kv.second;

        Face f;
        f.id = faceId++;
        f.nodes = {edge.first, edge.second};

        // Координаты концов ребра
        double x1 = nodeCoords[3*(edge.first-1)];
        double y1 = nodeCoords[3*(edge.first-1)+1];
        double x2 = nodeCoords[3*(edge.second-1)];
        double y2 = nodeCoords[3*(edge.second-1)+1];
        double dx = x2 - x1;
        double dy = y2 - y1;
        f.length = hypot(dx, dy);
        // Единичная нормаль (перпендикуляр), позже уточним направление
        f.nx =  dy / f.length;
        f.ny = -dx / f.length;

        if (cellIndices.size() == 2) {
            // Внутренняя грань
            f.type = 0;
            f.cellL = cellIndices[0];
            f.cellR = cellIndices[1];
            f.phys_tag = -1;
            // Уточняем нормаль: она должна быть направлена из cellL в cellR
            // Для этого вычисляем вектор из центра cellL в середину грани
            double cxL = cells[f.cellL].xc, cyL = cells[f.cellL].yc;
            double cxR = cells[f.cellR].xc, cyR = cells[f.cellR].yc;
            double mx = (x1 + x2) / 2.0;
            double my = (y1 + y2) / 2.0;
            double dxL = mx - cxL, dyL = my - cyL;
            double dxR = mx - cxR, dyR = my - cyR;
            double dotL = f.nx * dxL + f.ny * dyL;
            if (dotL < 0) {
                // Нормаль смотрит внутрь cellL, разворачиваем
                f.nx = -f.nx;
                f.ny = -f.ny;
            }
        } else if (cellIndices.size() == 1) {
            // Граничная грань
            f.type = 1;
            f.cellL = cellIndices[0];
            f.cellR = -1;
            // Физический тег пока неизвестен – его нужно получить из граничных элементов MSH
            // Для простоты оставим 0, потом можно будет сопоставить по узлам, но здесь не реализовано.
            // В рабочем коде нужно прочитать линии из файла и сопоставить.
            f.phys_tag = 0;
            // Нормаль направлена наружу от cellL – проверяем и при необходимости разворачиваем
            double cx = cells[f.cellL].xc, cy = cells[f.cellL].yc;
            double mx = (x1 + x2) / 2.0;
            double my = (y1 + y2) / 2.0;
            double dx = mx - cx, dy = my - cy;
            double dot = f.nx * dx + f.ny * dy;
            if (dot < 0) {
                f.nx = -f.nx;
                f.ny = -f.ny;
            }
        } else {
            cerr << "Ошибка: ребро принадлежит более чем двум ячейкам!" << endl;
            exit(1);
        }
        faces.push_back(f);
    }
    cout << "Построено граней: " << faces.size() << endl;
}

// ============================================================
// Функции для решения задачи Римана (из предоставленного кода)
// ============================================================
double WaveFunk(double p, double pk, double rho_k) {
    double gamma = GAMMA;
    if (p > pk) {
        double Ak = 2.0 / ((gamma + 1) * rho_k);
        double Bk = (gamma - 1) / (gamma + 1) * pk;
        return (p - pk) * sqrt(Ak / (p + Bk));
    } else {
        return 2.0 * sqrt(gamma * pk / rho_k) / (gamma - 1) *
               (pow(p / pk, (gamma - 1) / (2 * gamma)) - 1.0);
    }
}

double ProisvWaveFunk(double p, double pk, double rho_k) {
    double gamma = GAMMA;
    if (p > pk) {
        double Ak = 2.0 / ((gamma + 1) * rho_k);
        double Bk = (gamma - 1) / (gamma + 1) * pk;
        return sqrt(Ak / (p + Bk)) * (1.0 - (p - pk) / (2.0 * (p + Bk)));
    } else {
        return (1.0 / (rho_k * sqrt(gamma * pk / rho_k))) *
               pow(p / pk, -(gamma + 1) / (2 * gamma));
    }
}

double solve_p_star(double rhoL, double uL, double pL, double rhoR, double uR, double pR) {
    double gamma = GAMMA, tol = 1e-8;
    int max_iter = 100;
    double p0 = 0.5 * (pL + pR) - 0.125 * (uR - uL) * (rhoL + rhoR) *
                (sqrt(gamma * pL / rhoL) + sqrt(gamma * pR / rhoR));
    p0 = max(tol, p0);
    double p_star = p0;
    for (int i = 0; i < max_iter; ++i) {
        double fL = WaveFunk(p_star, pL, rhoL);
        double fR = WaveFunk(p_star, pR, rhoR);
        double f = fL + fR + uR - uL;
        double dfL = ProisvWaveFunk(p_star, pL, rhoL);
        double dfR = ProisvWaveFunk(p_star, pR, rhoR);
        double df = dfL + dfR;
        if (fabs(df) < 1e-12) break;
        double dp = -f / df;
        p_star += dp;
        if (p_star < tol) p_star = tol;
        if (fabs(dp) < tol * p_star) break;
    }
    return p_star;
}

double u_star(double p_star, double rhoL, double uL, double pL, double rhoR, double uR, double pR) {
    double fL = WaveFunk(p_star, pL, rhoL);
    double fR = WaveFunk(p_star, pR, rhoR);
    return 0.5 * (uL + uR) + 0.5 * (fR - fL);
}

tuple<double, double, double> sample(
    double p_star, double u_star, double rhoL, double uL, double pL,
    double rhoR, double uR, double pR, double x_over_t) {
    double gamma = GAMMA;
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
            double SHL = uL - cL;
            double STL = u_star - sqrt(gamma * p_star / (rhoL * pow(p_star / pL, 1.0 / gamma)));
            if (x_over_t < SHL) {
                return {rhoL, uL, pL};
            } else if (x_over_t < STL) {
                double u = 2.0 / (gamma + 1) * (cL + (gamma - 1) / 2 * uL + x_over_t);
                double c = 2.0 / (gamma + 1) * (cL + (gamma - 1) / 2 * (uL - x_over_t));
                double rho = rhoL * pow(c / cL, 2.0 / (gamma - 1));
                double p = pL * pow(rho / rhoL, gamma);
                return {rho, u, p};
            } else {
                double rho_star = rhoL * pow(p_star / pL, 1.0 / gamma);
                return {rho_star, u_star, p_star};
            }
        }
    } else {
        if (p_star > pR) {
            double SR = uR + sqrt(((gamma + 1) / (2 * gamma)) * (p_star / pR) +
                                 (gamma - 1) / (2 * gamma)) * cR;
            if (x_over_t > SR) {
                return {rhoR, uR, pR};
            } else {
                double rho_star = rhoR * ((p_star / pR) + (gamma - 1) / (gamma + 1)) /
                                 (((gamma - 1) / (gamma + 1)) * (p_star / pR) + 1);
                return {rho_star, u_star, p_star};
            }
        } else {
            double SHR = uR + cR;
            double STR = u_star + sqrt(gamma * p_star / (rhoR * pow(p_star / pR, 1.0 / gamma)));
            if (x_over_t > SHR) {
                return {rhoR, uR, pR};
            } else if (x_over_t > STR) {
                double u = 2.0 / (gamma + 1) * (-cR + (gamma - 1) / 2 * uR + x_over_t);
                double c = 2.0 / (gamma + 1) * (cR - (gamma - 1) / 2 * (uR - x_over_t));
                double rho = rhoR * pow(c / cR, 2.0 / (gamma - 1));
                double p = pR * pow(rho / rhoR, gamma);
                return {rho, u, p};
            } else {
                double rho_star = rhoR * pow(p_star / pR, 1.0 / gamma);
                return {rho_star, u_star, p_star};
            }
        }
    }
}

// Точное решение Римана
tuple<double, double, double> exactRiemann(double rhoL, double uL, double pL,
                                           double rhoR, double uR, double pR, double xi) {
    double p_star = solve_p_star(rhoL, uL, pL, rhoR, uR, pR);
    double u_star_val = u_star(p_star, rhoL, uL, pL, rhoR, uR, pR);
    return sample(p_star, u_star_val, rhoL, uL, pL, rhoR, uR, pR, xi);
}

// Решатель Роу
tuple<double, double, double> roeRiemann(double rhoL, double uL, double pL,
                                         double rhoR, double uR, double pR, double xi) {
    double gamma = GAMMA;
    double EL = pL/(gamma-1.0) + 0.5*rhoL*uL*uL;
    double ER = pR/(gamma-1.0) + 0.5*rhoR*uR*uR;
    double sqrt_rhoL = sqrt(rhoL);
    double sqrt_rhoR = sqrt(rhoR);
    double sum_sqrt = sqrt_rhoL + sqrt_rhoR;
    double u_tilde = (sqrt_rhoL*uL + sqrt_rhoR*uR) / sum_sqrt;
    double HL = (EL + pL) / rhoL;
    double HR = (ER + pR) / rhoR;
    double H_tilde = (sqrt_rhoL*HL + sqrt_rhoR*HR) / sum_sqrt;
    double rho_tilde = sqrt(rhoL*rhoR);
    double a2_tilde = (gamma-1.0)*(H_tilde - 0.5*u_tilde*u_tilde);
    double a_tilde = sqrt(a2_tilde);
    double delta_p = pR - pL;
    double delta_u = uR - uL;
    double delta_rho = rhoR - rhoL;
    double alpha1 = (delta_p - rho_tilde*a_tilde*delta_u) / (2.0*a2_tilde);
    double alpha2 = delta_rho - delta_p / a2_tilde;
    double alpha3 = (delta_p + rho_tilde*a_tilde*delta_u) / (2.0*a2_tilde);
    double lambda1 = u_tilde - a_tilde;
    double lambda2 = u_tilde;
    double lambda3 = u_tilde + a_tilde;
    double U1 = rhoL;
    double U2 = rhoL*uL;
    double U3 = EL;
    if (xi > lambda1) {
        U1 += alpha1 * 1.0;
        U2 += alpha1 * (u_tilde - a_tilde);
        U3 += alpha1 * (H_tilde - u_tilde*a_tilde);
    }
    if (xi > lambda2) {
        U1 += alpha2 * 1.0;
        U2 += alpha2 * u_tilde;
        U3 += alpha2 * (0.5*u_tilde*u_tilde);
    }
    if (xi > lambda3) {
        U1 += alpha3 * 1.0;
        U2 += alpha3 * (u_tilde + a_tilde);
        U3 += alpha3 * (H_tilde + u_tilde*a_tilde);
    }
    double rho = U1;
    double u = U2 / rho;
    double e = U3 / rho - 0.5*u*u;
    double p = rho*e*(gamma-1.0);
    if (p < 0.0) p = 1e-10;
    return {rho, u, p};
}

// Решатель Осер-Соломон (из вашего кода)
tuple<double, double, double> osherSolomonRiemann(double rhoL, double uL, double pL,
                                                  double rhoR, double uR, double pR,
                                                  double xi, double gamma) {
    const double eps = 1e-12;
    const double gamma_minus_1 = gamma - 1.0;
    double EL = pL/gamma_minus_1 + 0.5*rhoL*uL*uL;
    double ER = pR/gamma_minus_1 + 0.5*rhoR*uR*uR;
    array<double,3> UL = {rhoL, rhoL*uL, EL};
    array<double,3> UR = {rhoR, rhoR*uR, ER};
    array<double,3> FL = {rhoL*uL, rhoL*uL*uL + pL, uL*(EL + pL)};
    array<double,3> FR = {rhoR*uR, rhoR*uR*uR + pR, uR*(ER + pR)};
    constexpr int Q = 3;
    const double theta[Q] = {0.11270166537925831148, 0.5, 0.88729833462074168852};
    const double weight[Q] = {5.0/18.0, 4.0/9.0, 5.0/18.0};
    array<double,3> integral = {0.0,0.0,0.0};
    double drho = rhoR - rhoL;
    double du = uR - uL;
    double dp = pR - pL;
    for (int k=0; k<Q; ++k) {
        double psi = theta[k];
        double rho = rhoL + psi*drho;
        double u = uL + psi*du;
        double p = pL + psi*dp;
        array<double,3> dU_dpsi;
        dU_dpsi[0] = drho;
        dU_dpsi[1] = drho*u + rho*du;
        dU_dpsi[2] = dp/gamma_minus_1 + 0.5*drho*u*u + rho*u*du;
        double a = sqrt(gamma*p / max(rho, eps));
        double E = p/gamma_minus_1 + 0.5*rho*u*u;
        double H = (E + p) / rho;
        double lambda1 = u - a;
        double lambda2 = u;
        double lambda3 = u + a;
        double beta = (gamma-1.0) / (2.0*a*a);
        double l11 = 0.5*(beta*u*u + u/a);
        double l12 = -0.5*(beta*u + 1.0/a);
        double l13 = 0.5*beta;
        double l21 = 1.0 - beta*u*u;
        double l22 = beta*u;
        double l23 = -beta;
        double l31 = 0.5*(beta*u*u - u/a);
        double l32 = -0.5*(beta*u - 1.0/a);
        double l33 = 0.5*beta;
        double alpha1 = l11*dU_dpsi[0] + l12*dU_dpsi[1] + l13*dU_dpsi[2];
        double alpha2 = l21*dU_dpsi[0] + l22*dU_dpsi[1] + l23*dU_dpsi[2];
        double alpha3 = l31*dU_dpsi[0] + l32*dU_dpsi[1] + l33*dU_dpsi[2];
        alpha1 *= fabs(lambda1);
        alpha2 *= fabs(lambda2);
        alpha3 *= fabs(lambda3);
        array<double,3> term;
        term[0] = alpha1 + alpha2 + alpha3;
        term[1] = (u-a)*alpha1 + u*alpha2 + (u+a)*alpha3;
        term[2] = (H - u*a)*alpha1 + 0.5*u*u*alpha2 + (H + u*a)*alpha3;
        integral[0] += weight[k] * term[0];
        integral[1] += weight[k] * term[1];
        integral[2] += weight[k] * term[2];
    }
    array<double,3> F_osher;
    for (int i=0; i<3; ++i) {
        F_osher[i] = 0.5*(FL[i] + FR[i]) - 0.5*integral[i];
    }
    // Для xi=0 возвращаем поток F_osher
    return {F_osher[0], F_osher[1], F_osher[2]};
}

// Решатель HLL
tuple<double, double, double> hllRiemann(double rhoL, double uL, double pL,
                                         double rhoR, double uR, double pR, double xi) {
    double gamma = GAMMA;
    double cL = sqrt(gamma * pL / rhoL);
    double cR = sqrt(gamma * pR / rhoR);
    double SL = uL - cL;
    double SR = uR + cR;
    double fluxL_m = rhoL * uL;
    double fluxL_mom = rhoL * uL * uL + pL;
    double fluxL_E = uL * (pL/(gamma-1) + 0.5*rhoL*uL*uL + pL);
    double fluxR_m = rhoR * uR;
    double fluxR_mom = rhoR * uR * uR + pR;
    double fluxR_E = uR * (pR/(gamma-1) + 0.5*rhoR*uR*uR + pR);
    double f_m, f_mom, f_E;
    if (xi <= SL) {
        f_m = fluxL_m; f_mom = fluxL_mom; f_E = fluxL_E;
    } else if (xi >= SR) {
        f_m = fluxR_m; f_mom = fluxR_mom; f_E = fluxR_E;
    } else {
        double UL_m = rhoL;
        double UL_mom = rhoL * uL;
        double UL_E = pL/(gamma-1) + 0.5*rhoL*uL*uL;
        double UR_m = rhoR;
        double UR_mom = rhoR * uR;
        double UR_E = pR/(gamma-1) + 0.5*rhoR*uR*uR;
        f_m = (SR*fluxL_m - SL*fluxR_m + SL*SR*(UR_m - UL_m)) / (SR - SL);
        f_mom = (SR*fluxL_mom - SL*fluxR_mom + SL*SR*(UR_mom - UL_mom)) / (SR - SL);
        f_E = (SR*fluxL_E - SL*fluxR_E + SL*SR*(UR_E - UL_E)) / (SR - SL);
    }
    return {f_m, f_mom, f_E};
}

// Решатель HLLC
tuple<double, double, double> hllcRiemann(double rhoL, double uL, double pL,
                                          double rhoR, double uR, double pR, double xi) {
    double gamma = GAMMA;
    double cL = sqrt(gamma * pL / rhoL);
    double cR = sqrt(gamma * pR / rhoR);
    double SL = uL - cL;
    double SR = uR + cR;
    // Оценка для S* (по Batten et al.)
    double rho_avg = 0.5*(rhoL+rhoR);
    double c_avg = 0.5*(cL+cR);
    double Sstar = (rhoR*uR*(SR-uR) + rhoL*uL*(uL-SL) + pL - pR) /
                   (rhoR*(SR-uR) + rhoL*(uL-SL));
    double fluxL_m = rhoL * uL;
    double fluxL_mom = rhoL * uL*uL + pL;
    double fluxL_E = uL * (pL/(gamma-1) + 0.5*rhoL*uL*uL + pL);
    double fluxR_m = rhoR * uR;
    double fluxR_mom = rhoR * uR*uR + pR;
    double fluxR_E = uR * (pR/(gamma-1) + 0.5*rhoR*uR*uR + pR);
    double f_m, f_mom, f_E;
    if (xi <= SL) {
        f_m = fluxL_m; f_mom = fluxL_mom; f_E = fluxL_E;
    } else if (xi >= SR) {
        f_m = fluxR_m; f_mom = fluxR_mom; f_E = fluxR_E;
    } else if (xi <= Sstar) {
        // левый промежуток
        double Ustar_m = rhoL * (SL - uL) / (SL - Sstar);
        double Ustar_mom = Ustar_m * Sstar;
        double Ustar_E = Ustar_m * ( (pL/(gamma-1) + 0.5*rhoL*uL*uL)/rhoL +
                                    (Sstar - uL)*(Sstar + pL/(rhoL*(SL-uL))) );
        f_m = fluxL_m + SL*(Ustar_m - rhoL);
        f_mom = fluxL_mom + SL*(Ustar_mom - rhoL*uL);
        f_E = fluxL_E + SL*(Ustar_E - (pL/(gamma-1) + 0.5*rhoL*uL*uL));
    } else {
        // правый промежуток
        double Ustar_m = rhoR * (SR - uR) / (SR - Sstar);
        double Ustar_mom = Ustar_m * Sstar;
        double Ustar_E = Ustar_m * ( (pR/(gamma-1) + 0.5*rhoR*uR*uR)/rhoR +
                                    (Sstar - uR)*(Sstar + pR/(rhoR*(SR-uR))) );
        f_m = fluxR_m + SR*(Ustar_m - rhoR);
        f_mom = fluxR_mom + SR*(Ustar_mom - rhoR*uR);
        f_E = fluxR_E + SR*(Ustar_E - (pR/(gamma-1) + 0.5*rhoR*uR*uR));
    }
    return {f_m, f_mom, f_E};
}

// Решатель Русанова
tuple<double, double, double> rusanovRiemann(double rhoL, double uL, double pL,
                                             double rhoR, double uR, double pR, double xi) {
    double gamma = GAMMA;
    double cL = sqrt(gamma * pL / rhoL);
    double cR = sqrt(gamma * pR / rhoR);
    double Smax = max(fabs(uL)+cL, fabs(uR)+cR);
    double fluxL_m = rhoL * uL;
    double fluxL_mom = rhoL * uL*uL + pL;
    double fluxL_E = uL * (pL/(gamma-1) + 0.5*rhoL*uL*uL + pL);
    double fluxR_m = rhoR * uR;
    double fluxR_mom = rhoR * uR*uR + pR;
    double fluxR_E = uR * (pR/(gamma-1) + 0.5*rhoR*uR*uR + pR);
    double f_m = 0.5*(fluxL_m + fluxR_m) - 0.5*Smax*(rhoR - rhoL);
    double f_mom = 0.5*(fluxL_mom + fluxR_mom) - 0.5*Smax*(rhoR*uR - rhoL*uL);
    double f_E = 0.5*(fluxL_E + fluxR_E) - 0.5*Smax*((pR/(gamma-1)+0.5*rhoR*uR*uR) -
                                                    (pL/(gamma-1)+0.5*rhoL*uL*uL));
    return {f_m, f_mom, f_E};
}

// Универсальная функция: возвращает поток (масса, импульс, энергия) через грань
// с нормалью n (единичная), для xi=0 (т.е. поток в направлении нормали)
tuple<double, double, double> riemannSolver(const Primitive& left, const Primitive& right,
                                            double nx, double ny, const string& method) {
    // Приводим к одномерной задаче: нормальная скорость
    double uL_n = left.u * nx + left.v * ny;
    double uR_n = right.u * nx + right.v * ny;
    double rhoL = left.rho, pL = left.p;
    double rhoR = right.rho, pR = right.p;
    double xi = 0.0; // для потока в нормальном направлении

    tuple<double, double, double> res;
    if (method == "exact") {
        res = exactRiemann(rhoL, uL_n, pL, rhoR, uR_n, pR, xi);
    } else if (method == "roe") {
        res = roeRiemann(rhoL, uL_n, pL, rhoR, uR_n, pR, xi);
    } else if (method == "osher") {
        res = osherSolomonRiemann(rhoL, uL_n, pL, rhoR, uR_n, pR, xi, GAMMA);
    } else if (method == "hll") {
        res = hllRiemann(rhoL, uL_n, pL, rhoR, uR_n, pR, xi);
    } else if (method == "hllc") {
        res = hllcRiemann(rhoL, uL_n, pL, rhoR, uR_n, pR, xi);
    } else if (method == "rusanov") {
        res = rusanovRiemann(rhoL, uL_n, pL, rhoR, uR_n, pR, xi);
    } else {
        cerr << "Unknown Riemann solver: " << method << endl;
        exit(1);
    }
    // res: (поток массы, поток нормального импульса, поток энергии)
    // поток импульса возвращается как скаляр в нормальном направлении,
    // его нужно разложить на компоненты x,y, умножив на nx, ny
    double flux_mass = get<0>(res);
    double flux_mom_n = get<1>(res);
    double flux_energy = get<2>(res);
    double flux_mom_x = flux_mom_n * nx;
    double flux_mom_y = flux_mom_n * ny;
    return {flux_mass, flux_mom_x, flux_mom_y, flux_energy};
}

// ============================================================
// Граничные условия
// ============================================================
Primitive getBoundaryState(const Primitive& interior, int phys_tag, double nx, double ny) {
    // Для простоты реализуем несколько типов:
    // phys_tag = 10 (inflow) – зададим фиксированные значения (например, сверхзвуковой поток)
    // phys_tag = 11 (outflow) – экстраполяция нулевой производной
    // phys_tag = 12 (wall) – зеркальное отражение по нормальной скорости
    Primitive bcState = interior;
    if (phys_tag == 10) { // inflow: задаём значения (например, набегающий поток)
        // Здесь можно задать параметры свободного потока
        bcState.rho = 1.0;
        bcState.u = 1.0;  // горизонтальная скорость
        bcState.v = 0.0;
        bcState.p = 1.0 / GAMMA; // чтобы звуковая скорость была 1
    } else if (phys_tag == 11) { // outflow: экстраполяция
        // ничего не меняем, оставляем interior
    } else if (phys_tag == 12) { // wall: непроницаемая стенка
        // Нормальная скорость меняет знак
        double u_n = interior.u*nx + interior.v*ny;
        double u_tx = interior.u - u_n*nx;
        double u_ty = interior.v - u_n*ny;
        bcState.u = u_tx - u_n*nx;  // отражаем нормальную компоненту
        bcState.v = u_ty - u_n*ny;
        // плотность и давление не меняются (адиабатическая стенка)
    }
    return bcState;
}

// ============================================================
// Обновление решения (один шаг по времени)
// ============================================================
void updateSolution(double dt, const string& method) {
    // Обнуляем приращения
    for (auto& cell : cells) {
        cell.dU = {0,0,0,0};
    }

    // Проходим по всем граням
    for (const auto& face : faces) {
        Conserved flux; // поток массы, x-имп, y-имп, энергии

        // Левые и правые примитивные переменные
        Primitive left = getPrimitive(cells[face.cellL].U);
        Primitive right;

        if (face.type == 0) { // внутренняя грань
            right = getPrimitive(cells[face.cellR].U);
        } else { // граничная грань
            right = getBoundaryState(left, face.phys_tag, face.nx, face.ny);
        }

        auto [f_m, f_mx, f_my, f_E] = riemannSolver(left, right, face.nx, face.ny, method);
        flux = {f_m, f_mx, f_my, f_E};

        // Добавляем в левую ячейку: -flux * length / area
        double invVolL = 1.0 / cells[face.cellL].area;
        cells[face.cellL].dU.rho   -= flux.rho   * face.length * invVolL;
        cells[face.cellL].dU.rhou  -= flux.rhou  * face.length * invVolL;
        cells[face.cellL].dU.rhov  -= flux.rhov  * face.length * invVolL;
        cells[face.cellL].dU.rhoE  -= flux.rhoE  * face.length * invVolL;

        if (face.type == 0) { // внутренняя: добавить с противоположным знаком в правую ячейку
            double invVolR = 1.0 / cells[face.cellR].area;
            cells[face.cellR].dU.rho   += flux.rho   * face.length * invVolR;
            cells[face.cellR].dU.rhou  += flux.rhou  * face.length * invVolR;
            cells[face.cellR].dU.rhov  += flux.rhov  * face.length * invVolR;
            cells[face.cellR].dU.rhoE  += flux.rhoE  * face.length * invVolR;
        }
    }

    // Обновляем консервативные переменные
    for (auto& cell : cells) {
        cell.U.rho   += dt * cell.dU.rho;
        cell.U.rhou  += dt * cell.dU.rhou;
        cell.U.rhov  += dt * cell.dU.rhov;
        cell.U.rhoE  += dt * cell.dU.rhoE;

        // Проверка на положительность
        if (cell.U.rho <= 0.0) {
            cerr << "Negative density in cell " << cell.id << endl;
            cell.U.rho = 1e-10;
        }
        Primitive prim = getPrimitive(cell.U);
        if (prim.p <= 0.0) {
            // Восстанавливаем давление
            double e = cell.U.rhoE / cell.U.rho - 0.5*(prim.u*prim.u + prim.v*prim.v);
            if (e <= 0.0) e = 1e-10;
            double p = (GAMMA-1.0) * cell.U.rho * e;
            cell.U.rhoE = cell.U.rho * (e + 0.5*(prim.u*prim.u + prim.v*prim.v));
        }
    }
}

// ============================================================
// Вычисление максимального шага по времени (CFL)
// ============================================================
double computeTimeStep(double cfl) {
    double max_lambda = 0.0;
    for (const auto& cell : cells) {
        Primitive prim = getPrimitive(cell.U);
        double a = sqrt(GAMMA * prim.p / prim.rho);
        double u_abs = sqrt(prim.u*prim.u + prim.v*prim.v);
        double lambda = u_abs + a;
        // Характерный размер ячейки: sqrt(площади)
        double h = sqrt(cell.area);
        double local = lambda / h;
        if (local > max_lambda) max_lambda = local;
    }
    double dt = cfl / max_lambda;
    return dt;
}

// ============================================================
// Запись результатов в VTK
// ============================================================
void writeVTK(int iter, double time) {
    string filename = "solution_" + to_string(iter) + ".vtk";
    ofstream vtk(filename);
    vtk << "# vtk DataFile Version 3.0\n";
    vtk << "Godunov solution, time = " << time << "\n";
    vtk << "ASCII\n";
    vtk << "DATASET UNSTRUCTURED_GRID\n";

    // Точки (все узлы)
    vtk << "POINTS " << nodeNum << " float\n";
    for (int i = 0; i < nodeNum; ++i) {
        vtk << nodeCoords[3*i] << " " << nodeCoords[3*i+1] << " " << nodeCoords[3*i+2] << "\n";
    }

    // Ячейки (треугольники)
    vtk << "CELLS " << cells.size() << " " << 4 * cells.size() << "\n";
    for (const auto& cell : cells) {
        vtk << "3 " << cell.nodes[0]-1 << " " << cell.nodes[1]-1 << " " << cell.nodes[2]-1 << "\n";
    }
    vtk << "CELL_TYPES " << cells.size() << "\n";
    for (size_t i = 0; i < cells.size(); ++i) vtk << "5\n";

    // Данные на ячейках
    vtk << "CELL_DATA " << cells.size() << "\n";
    vtk << "SCALARS density float 1\n";
    vtk << "LOOKUP_TABLE default\n";
    for (const auto& cell : cells) vtk << cell.U.rho << "\n";

    vtk << "VECTORS velocity float\n";
    for (const auto& cell : cells) {
        double u = cell.U.rhou / cell.U.rho;
        double v = cell.U.rhov / cell.U.rho;
        vtk << u << " " << v << " 0\n";
    }

    vtk << "SCALARS pressure float 1\n";
    vtk << "LOOKUP_TABLE default\n";
    for (const auto& cell : cells) {
        Primitive prim = getPrimitive(cell.U);
        vtk << prim.p << "\n";
    }

    vtk.close();
    cout << "Saved " << filename << endl;
}

// ============================================================
// Начальные условия
// ============================================================
void setInitialConditions() {
    // Пример: однородный поток с параметрами, соответствующими inflow
    // (можно задать разные условия для подобластей, если есть)
    for (auto& cell : cells) {
        // Задаём примитивные переменные
        double rho0 = 1.0;
        double u0 = 1.0;
        double v0 = 0.0;
        double p0 = 1.0 / GAMMA; // чтобы c = 1
        Primitive prim0 = {rho0, u0, v0, p0};
        cell.U = getConserved(prim0);
    }
}

// ============================================================
// Основная функция
// ============================================================
int main(int argc, char** argv) {
    if (argc < 2) {
        cerr << "Usage: " << argv[0] << " mesh.msh [method]" << endl;
        cerr << "  method: exact, roe, osher, hll, hllc, rusanov (default roe)" << endl;
        return 1;
    }

    string meshFile = argv[1];
    string method = "roe";
    if (argc >= 3) method = argv[2];

    // 1. Чтение сетки
    if (!readMesh(meshFile)) {
        cerr << "Failed to read mesh file: " << meshFile << endl;
        return 1;
    }

    // 2. Построение ячеек и граней
    buildCells();
    buildFaces();

    // 3. Начальные условия
    setInitialConditions();

    // 4. Временной цикл
    double t = 0.0;
    int iter = 0;
    while (t < T_END) {
        double dt = computeTimeStep(CFL);
        if (t + dt > T_END) dt = T_END - t;

        updateSolution(dt, method);
        t += dt;
        iter++;

        if (iter % SAVE_INTERVAL == 0 || t >= T_END) {
            writeVTK(iter, t);
        }

        cout << "Iter " << iter << ", t = " << t << ", dt = " << dt << endl;
    }

    cout << "Simulation finished." << endl;
    return 0;
}