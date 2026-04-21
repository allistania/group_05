#ifndef PARAMETERS_H
#define PARAMETERS_H

#include <string>
#include <vector>

// Функции для работы с тестовыми параметрами
void writeTest();
bool readTest(const std::string& filename, 
              double& rho_L, double& u_L, double& p_L,
              double& rho_R, double& u_R, double& p_R, 
              double& t_end);

// Функции для работы с параметрами инициализации
void start();

bool readInit(const std::string& filename, double& tau, double& h, 
              int& cells, double& cfl, std::string& left_bc, 
              std::string& right_bc, std::string& test_name, int& count, 
              std::vector<std::string>& solver_names, std::string& t_end_user, 
              double& t_output, std::string& key);
// Функции для готовых тестов Сода
void setupSODTest(const std::string& test_name, 
                  double& rho_L, double& u_L, double& p_L,
                  double& rho_R, double& u_R, double& p_R, 
                  double& t_end, const std::string& t_end_user);

// Функции для граничных условий
void applyBoundaryConditions(double* rho, double* u, double* p, 
                            int cells, int ghost_cells,
                            const std::string& left_bc, 
                            const std::string& right_bc);

int getGhostCellsBySolver(const std::string& solver_name);
#endif // PARAMETERS_H