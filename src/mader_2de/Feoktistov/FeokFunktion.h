#ifndef FEOKFUNCTIONS_H 
#define FEOKFUNCTIONS_H
#include <iostream>
#include <cmath>
#include <algorithm>
#include <tuple>
#include <vector>
using namespace std;
double WaveFunk(double p, double pk, double rho_k);
double ProisvWaveFunk(double p, double pk, double rho_k);
double solve_p_star(double rhoL, double uL, double pL, double rhoR, double uR, double pR);
tuple<double, double, double> sample(
    double p_star, double u_star, double rhoL, double uL, double pL,
    double rhoR, double uR, double pR, double x_over_t, double gamma = 1.4
);
double  u_star(double p_star,double  rhoL,double  uL,double  pL,double  rhoR,double uR,double  pR);

void GodunovSolve2D(
    int cells_x, int cells_y, int ghost_cells,
    std::vector<std::vector<double>>& u,
    std::vector<std::vector<double>>& v,
    std::vector<std::vector<double>>& P,
    std::vector<std::vector<double>>& rho,
    std::vector<std::vector<double>>& I,
    double& t_total, double tau, double hx, double hy,
    const std::string& left_bc,   const std::string& right_bc,
    const std::string& bottom_bc, const std::string& top_bc,
    int Q, const std::string& time_integrator, const std::string& nameRS);

std::tuple<double, double, double> roe(
    double rhoL, double uL, double pL,
    double rhoR, double uR, double pR,
    double xi, double gamma = 1.4);

std::tuple<double,double,double> superSolve(double rhoL,double uL,double pL,
double rhoR,double uR,double pR,
std::string name);

void RodionovSolve2D(
    int cells_x, int cells_y, int ghost_cells,
    std::vector<std::vector<double>>& u,
    std::vector<std::vector<double>>& v,
    std::vector<std::vector<double>>& P,
    std::vector<std::vector<double>>& rho,
    std::vector<std::vector<double>>& I,
    double& t_total, double tau, double hx, double hy,
    const std::string& left_bc, const std::string& right_bc,
    const std::string& bottom_bc, const std::string& top_bc,
    const std::string& time_integrator, const std::string& nameRS);

void ENO2DSolve(
    int cells_x, int cells_y, int ghost_cells,
    std::vector<std::vector<double>>& u,
    std::vector<std::vector<double>>& v,
    std::vector<std::vector<double>>& P,
    std::vector<std::vector<double>>& rho,
    std::vector<std::vector<double>>& I,
    double& t_total, double tau, double hx, double hy,
    const std::string& left_bc, const std::string& right_bc,
    const std::string& bottom_bc, const std::string& top_bc,
    const std::string& time_integrator,
    const std::string& nameRS);

void WENO2DSolve(
    int cells_x, int cells_y, int ghost_cells,
    std::vector<std::vector<double>>& u,
    std::vector<std::vector<double>>& v,
    std::vector<std::vector<double>>& P,
    std::vector<std::vector<double>>& rho,
    std::vector<std::vector<double>>& I,
    double& t_total, double tau, double hx, double hy,
    const std::string& left_bc, const std::string& right_bc,
    const std::string& bottom_bc, const std::string& top_bc,
    const std::string& time_integrator,
    const std::string& nameRS);

void Fletcher2DSolve(
    int cells_x, int cells_y, int ghost_cells,
    std::vector<std::vector<double>>& u,
    std::vector<std::vector<double>>& v,
    std::vector<std::vector<double>>& P,
    std::vector<std::vector<double>>& rho,
    std::vector<std::vector<double>>& I,
    double& t_total, double tau, double hx, double hy,
    const std::string& left_bc, const std::string& right_bc,
    const std::string& bottom_bc, const std::string& top_bc,
    double viscosity_coeff, double anti_diffusion_coeff);
void saveToCSV(
    const std::string& filename,
    const std::vector<std::vector<double>>& rho,
    const std::vector<std::vector<double>>& u,
    const std::vector<std::vector<double>>& v,
    const std::vector<std::vector<double>>& P,
    const std::vector<std::vector<double>>& I,
    double hx, double hy,
    int cells_x, int cells_y, int ghost_cells,
    double t);
void newTimeStep2D(
    const std::vector<std::vector<double>>& u,
    const std::vector<std::vector<double>>& v,
    const std::vector<std::vector<double>>& P,
    const std::vector<std::vector<double>>& rho,
    double& tau,
    double hx,
    double hy,
    double CFL,
    int ghost_cells);
void setInitialConditions(
    std::vector<std::vector<double>>& rho,
    std::vector<std::vector<double>>& u,
    std::vector<std::vector<double>>& v,
    std::vector<std::vector<double>>& P,
    std::vector<std::vector<double>>& I,
    double Lx, double Ly,
    double hx, double hy,
    int ghost_cells,
    double gamma,
    const std::string& test_case);
bool readConfigFromFile(const std::string& filename,
                        double& Lx, double& Ly,
                        int& cells_x, int& cells_y,
                        int& ghost_cells,
                        double& CFL,double& tau, double& t_end,
                        std::string& time_integrator,
                        std::string& nameRS,
                        std::string& solver,
                        double& viscosity_coeff,
                        double& anti_diffusion_coeff,
                        std::string& left_bc,
                        std::string& right_bc,
                        std::string& bottom_bc,
                        std::string& top_bc,
                        int& save_every,
                        std::string& test_case);
void FLIC(
    int cells_x, int cells_y, int ghost_cells,
    std::vector<std::vector<double>>& u,
    std::vector<std::vector<double>>& v,
    std::vector<std::vector<double>>& P,
    std::vector<std::vector<double>>& rho,
    std::vector<std::vector<double>>& I,
    double& t_total, double tau, double hx, double hy,
    const std::string& left_bc, const std::string& right_bc,
    const std::string& bottom_bc, const std::string& top_bc,
    const std::string& time_integrator,
    const std::string& nameRS,
    bool order_xy);

void FLICCF(
    int cells_x, int cells_y, int ghost_cells,
    std::vector<std::vector<double>>& u,
    std::vector<std::vector<double>>& v,
    std::vector<std::vector<double>>& P,
    std::vector<std::vector<double>>& rho,
    std::vector<std::vector<double>>& I,
    double& t_total, double tau, double hx, double hy,
    const std::string& left_bc, const std::string& right_bc,
    const std::string& bottom_bc, const std::string& top_bc,
    const std::string& time_integrator,
    const std::string& nameRS);

void Mader(
    int cells_x, int cells_y, int ghost_cells,
    std::vector<std::vector<double>>& u,
    std::vector<std::vector<double>>& v,
    std::vector<std::vector<double>>& P,
    std::vector<std::vector<double>>& rho,
    std::vector<std::vector<double>>& I,
    std::vector<std::vector<double>>& W,          
    double& t_total, double tau, double hx, double hy,
    const std::string& left_bc, const std::string& right_bc,
    const std::string& bottom_bc, const std::string& top_bc,
    double viscosity_coeff,
    double MINWT, double GASW, int VCNT,
    double Z_freq, double E_act_over_R, double R_gas);
void SIMPLEPISO2DSolve(
    const std::string& solver_type,
    int cells_x, int cells_y, int ghost_cells,
    std::vector<std::vector<double>>& u,
    std::vector<std::vector<double>>& v,
    std::vector<std::vector<double>>& P,
    std::vector<std::vector<double>>& rho,
    std::vector<std::vector<double>>& I,
    double& t_total, double tau, double hx, double hy,
    const std::string& left_bc, const std::string& right_bc,
    const std::string& bottom_bc, const std::string& top_bc,
    double alpha_u, double alpha_p, int nCorr,
    double U_lid);
#endif