#ifndef FEOKFUNCTIONS_H 
#define FEOKFUNCTIONS_H
#include <iostream>
#include <cmath>
#include <algorithm>
#include <tuple>
#include <vector>
using namespace std;
void kabare(int cells, int ghost_cells,
std::vector<double>& u,
std::vector<double>& P,
std::vector<double>& rho,
std::vector<double>& I,
double& t_total, double tau, double h,
const std::string& left_bc,
const std::string& right_bc);
double WaveFunk(double p, double pk, double rho_k);
double ProisvWaveFunk(double p, double pk, double rho_k);
double solve_p_star(double rhoL, double uL, double pL, double rhoR, double uR, double pR);
tuple<double, double, double> sample(
    double p_star, double u_star, double rhoL, double uL, double pL,
    double rhoR, double uR, double pR, double x_over_t, double gamma = 1.4
);
double  u_star(double p_star,double  rhoL,double  uL,double  pL,double  rhoR,double uR,double  pR);
double  GodunovSolve(int cells,int ghost_cells,
               std::vector<double>& u,
               std::vector<double>& P,
               std::vector<double>& rho,
               std::vector<double>& I,
               double& t, double tau, double h, const std::string& left_bc, 
               const std::string& right_bc, int Q, const std::string& time_integrator, std::string nameRS);
void newTimeStep(std::vector<double>& u,
               std::vector<double>& P,
               std::vector<double>& rho,
               double& tau, double h, double CFL);

void  RodionovSolve(int cells,int ghost_cells,
               std::vector<double>& u,
               std::vector<double>& P,
               std::vector<double>& rho,
               std::vector<double>& I,
               double& t, double tau, double h, const std::string& left_bc, 
               const std::string& right_bc, const std::string& time_integrator, std::string nameRS);

void makeGIFfile(int cells, int ghost_cells, int g, double current_time, 
                 int save_step_interval, double save_time_interval, 
                 double h, std::vector<double>& u, std::vector<double>& P,
                 std::vector<double>& rho, std::vector<double>& I, 
                 double t_end, const std::string& way, bool use_time_interval = false);
double minmod(double a, double b);

tuple<double, double, double> ENOreconstruction(int i, std::vector<double>& u,
               std::vector<double>& P,
               std::vector<double>& rho, const std::string& time_integrator);

void  ENO(int cells,int ghost_cells,
               std::vector<double>& u,
               std::vector<double>& P,
               std::vector<double>& rho,
               std::vector<double>& I,
               double& t, double tau, double h, const std::string& left_bc, 
               const std::string& right_bc, const std::string& time_integrator, std::string nameRS);

void  WENO(int cells,int ghost_cells,
               std::vector<double>& u,
               std::vector<double>& P,
               std::vector<double>& rho,
               std::vector<double>& I,
               double& t, double tau, double h, const std::string& left_bc, 
               const std::string& right_bc, const std::string& time_integrator, std::string nameRS);

void  Fletcher(int cells,int ghost_cells,
               std::vector<double>& u,
               std::vector<double>& P,
               std::vector<double>& rho,
               std::vector<double>& I,
               double& t, double tau, double h, const std::string& left_bc, 
               const std::string& right_bc);

void feoktistov(int cells, int ghost_cells, 
                std::vector<double>& u, std::vector<double>& P, std::vector<double>& rho,
                std::vector<double>& I, double& t_total, double tau, double h,
                const std::string& left_bc, const std::string& right_bc,
                const std::string& high_order_scheme = "tvd",
                const std::string& reconstruction = "high",
                const std::string& tvd_scheme = "minmod",
                bool use_diffusion = false,
                bool use_anti_diffusion = false,
                bool use_artificial_viscosity = false,
                double viscosity_coeff = 0.01,
                double anti_diffusion_coeff = 0.05,
                const std::string& time_integrator = "euler", std::string nameRS="tocnoe");

void makeNormAnaliz(int count, int N_all, int start);

std::tuple<double, double, double> roe(
    double rhoL, double uL, double pL,
    double rhoR, double uR, double pR,
    double xi, double gamma = 1.4);

std::tuple<double,double,double> superSolve(double rhoL,double uL,double pL,
double rhoR,double uR,double pR,
std::string name);

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
    double& anti_diffusion_coeff); 

#endif