#include "Feoktistov/FeokFunktion.h"
#include "Batarin/parameters.h"
#include "Medyakova/Output.h"

int main()
{
    std::string config_file = "start.txt";

    double Lx, Ly;
    int cells_x, cells_y;
    int ghost_cells, save_every;
    double CFL, t_end, tau;
    std::string time_integrator, nameRS, solver;
    double viscosity_coeff, anti_diffusion_coeff;
    std::string left_bc, right_bc, bottom_bc, top_bc;
    const double GAMMA = 1.4; 
    std::string test_case; 
if (!readConfigFromFile(config_file,
                        Lx, Ly,
                        cells_x, cells_y,
                        ghost_cells,
                        CFL, tau,t_end,
                        time_integrator,
                        nameRS,
                        solver,
                        viscosity_coeff,
                        anti_diffusion_coeff,
                        left_bc,
                        right_bc,
                        bottom_bc,
                        top_bc,
                        save_every,
                        test_case)) {
    return 1;
} 

    prepareDirectories();
    setDefaultEndTime(test_case, t_end);
    int Nx = cells_x + 2 * ghost_cells;
    int Ny = cells_y + 2 * ghost_cells;
    double hx = Lx / cells_x;
    double hy = Ly / cells_y;

    std::vector<std::vector<double>> rho(Ny, std::vector<double>(Nx));
    std::vector<std::vector<double>> u(Ny, std::vector<double>(Nx));
    std::vector<std::vector<double>> v(Ny, std::vector<double>(Nx));
    std::vector<std::vector<double>> P(Ny, std::vector<double>(Nx));
    std::vector<std::vector<double>> I(Ny, std::vector<double>(Nx));

setInitialConditions(rho, u, v, P, I, Lx, Ly, hx, hy, ghost_cells, GAMMA, test_case);
    for (int j = 0; j < Ny; ++j) {
        applyBoundaryConditions(rho[j].data(), u[j].data(), P[j].data(),
                                cells_x, ghost_cells, left_bc, right_bc);
    }

    for (int i = 0; i < Nx; ++i) {
        std::vector<double> rho_col(Ny), u_col(Ny), v_col(Ny), P_col(Ny);
        for (int j = 0; j < Ny; ++j) {
            rho_col[j] = rho[j][i];
            u_col[j]   = u[j][i];
            v_col[j]   = v[j][i];
            P_col[j]   = P[j][i];
        }

        applyBoundaryConditions(rho_col.data(), u_col.data(), P_col.data(),
                                cells_y, ghost_cells, bottom_bc, top_bc);


        for (int j = 0; j < Ny; ++j) {
            rho[j][i] = rho_col[j];
            u[j][i]   = u_col[j];
            v[j][i]   = v_col[j];
            P[j][i]   = P_col[j];
        }
    }


    double t_total = 0.0;
    int step = 0;
    while (t_total < t_end - 1e-12) {
        double tau_step = std::min(tau, t_end - t_total);

        if (solver == "godunov") {
	    int Q=0;
             GodunovSolve2D(
                cells_x, cells_y, ghost_cells,
                u, v, P, rho, I,
                t_total, tau_step, hx, hy,
                left_bc, right_bc, bottom_bc, top_bc,Q,
                time_integrator, nameRS);
        }
        else if (solver == "Kolgan") {
	    int Q=1;
             GodunovSolve2D(
                cells_x, cells_y, ghost_cells,
                u, v, P, rho, I,
                t_total, tau_step, hx, hy,
                left_bc, right_bc, bottom_bc, top_bc,Q,
                time_integrator, nameRS);
        }
        else if (solver == "Rodionov") {
            RodionovSolve2D(
                cells_x, cells_y, ghost_cells,
                u, v, P, rho, I,
                t_total, tau_step, hx, hy,
                left_bc, right_bc, bottom_bc, top_bc,
                time_integrator, nameRS);
        }

        else if (solver == "eno") {
            ENO2DSolve(
                cells_x, cells_y, ghost_cells,
                u, v, P, rho, I,
                t_total, tau_step, hx, hy,
                left_bc, right_bc, bottom_bc, top_bc,
                time_integrator, nameRS);
        }
        else if (solver == "weno") {
            WENO2DSolve(
                cells_x, cells_y, ghost_cells,
                u, v, P, rho, I,
                t_total, tau_step, hx, hy,
                left_bc, right_bc, bottom_bc, top_bc,
                time_integrator, nameRS);
        }
        else if (solver == "fletcher") {
            Fletcher2DSolve(
                cells_x, cells_y, ghost_cells,
                u, v, P, rho, I,
                t_total, tau_step, hx, hy,
                left_bc, right_bc, bottom_bc, top_bc,
                viscosity_coeff, anti_diffusion_coeff);
        }
	else if (solver == "FLIC") {
            FLIC(
                cells_x, cells_y, ghost_cells,
                u, v, P, rho, I,
                t_total, tau_step, hx, hy,
                left_bc, right_bc, bottom_bc, top_bc,
                time_integrator, nameRS);
        }
        else {
            return 1;
        }
	if(step%save_every==0){
	std::string filename = "forGIF/G/frame_" + std::to_string(step) + ".csv";
        saveToCSV(filename, rho, u, v, P, I, hx, hy, cells_x, cells_y, ghost_cells, t_total);
	}
// CHECK: FLIC_CFL
	newTimeStep2D(u, v, P, rho, tau, hx, hy, CFL, ghost_cells);

        step++;
    }


    std::string filename = "output"".csv";
    saveToCSV(filename, rho, u, v, P, I, hx, hy, cells_x, cells_y, ghost_cells, t_total);

    return 0;
}