#include "Feoktistov/FeokFunktion.h"
#include "Medyakova/Output.h"
#include "Batarin/parameters.h"
int main() {
     
    std::string numerical_dir = "forGIF/G/";
    std::string exact_dir = "forGIF/tocnoe/";
    std::string data_dir = "data/";

    clearDirectory(numerical_dir);
    clearDirectory(exact_dir);   
    clearDirectory(data_dir);  

    double rho_L, u_L, p_L, rho_R, u_R, p_R, t_output;
    double x_over_t, t_end;
    double t, tau, CFL;
    double h;

    std::string nameRS, time_integrator, high_shem, tipe_reconsrutin, limiter;
    bool use_diffusion, use_anti_diffusion, use_artificial_viscosity;
    double viscosity_coeff, anti_diffusion_coeff;
    
    readStartFileF(
        nameRS,
        time_integrator,
        high_shem,
        tipe_reconsrutin,
        limiter,
        use_diffusion,
        use_anti_diffusion,
        use_artificial_viscosity,
        viscosity_coeff,
        anti_diffusion_coeff
    );

    double gamma=1.4;
    std::string filename = "init_params.txt", filename1 = "test_config.txt",
    left_bc, right_bc, test_name, solver_name;
    int cells, ghost_cells, count;
    std::string  t_end_user, key;
    std::vector<std::string> solver_names;

    if (readInit("start.txt", tau, h, cells, CFL, left_bc, right_bc, test_name, count, solver_names, t_end_user, t_output, key)) {
        std::cout << "Mesh is read" << std::endl;
        
        //                   
        if (test_name == "custom") {
            //                               
            readTest("custom_test.txt", rho_L, u_L, p_L, rho_R, u_R, p_R, t_end);
        } else {
            //                                
            setupSODTest(test_name, rho_L, u_L, p_L, rho_R, u_R, p_R, t_end, t_end_user);
        }
        
        int ghost_cells = getGhostCellsBySolver(solver_name);
        std::cout << "Solver: " << solver_name << " -> fict cells: " << ghost_cells << std::endl;
        //                                                
  ghost_cells=5;

        int total_cells = cells + 2 * ghost_cells;
        std::vector<double> rho(total_cells);
        std::vector<double> u(total_cells);
        std::vector<double> p(total_cells);
        std::vector<double> I(total_cells);
  int cells1 = cells+1;
      std::vector<double> rho1(cells), grid_P(cells);
      std::vector<double> u1(cells);
      std::vector<double> p1(cells);
      std::vector<double> I1(cells);
  //                      (               )
        int break_point = ghost_cells + cells / 2;
        for(auto m : solver_names){
  ghost_cells=5;
  int break_point = ghost_cells + cells / 2;

  double t_total=0;
        for (int i = 0; i < total_cells; i++) {
            if (i < break_point) {
    rho[i] = rho_L;
                u[i] = u_L;
                p[i] = p_L;
    I[i] = p[i]/(gamma-1)/rho[i];
    }
      else {
                rho[i] = rho_R;
                u[i] = u_R;
                p[i] = p_R;
    I[i] = p[i]/(gamma-1)/rho[i];

            }
           }
   
        //                            
       applyBoundaryConditions(rho.data(), u.data(), p.data(), 
                              cells, ghost_cells, left_bc, right_bc);
          
    int g=0;
    
    while(t_total<t_end){
    newTimeStep(u,p,rho,tau,h, CFL);
    if(m=="godunov"){ 
  GodunovSolve(cells, ghost_cells,u,p,rho, I,t_total, tau, h,left_bc, right_bc, 0,time_integrator, nameRS);
  }
    
    if(m=="godunov-kolgan"){ 
  GodunovSolve(cells, ghost_cells,u,p,rho, I,t_total, tau, h,left_bc, right_bc, 1,time_integrator, nameRS);
     }
    if(m=="godunov-kolgan-rodionov"){
  RodionovSolve(cells, ghost_cells,u,p,rho, I,t_total, tau, h,left_bc, right_bc,time_integrator, nameRS);
  }
    if(m=="eno"){
    ENO(cells, ghost_cells,u,p,rho, I,t_total, tau, h,left_bc, right_bc,time_integrator, nameRS);
  }
    if(m=="weno"){ 

  std::string A =  "forGIF/G/hislennoereshG", B = "forGIF/tocnoe/tochnoe" ;
      WENO(cells, ghost_cells,u,p,rho, I,t_total, tau, h,left_bc, right_bc,time_integrator, nameRS);

  makeGIFfile(cells,ghost_cells,g,t_total, count, t_output, h, u,p, rho, I, t_total, A, true);

      double p_star1=solve_p_star(rho_L, u_L, p_L, rho_R, u_R, p_R);
      double u_star1=u_star(p_star1,  rho_L, u_L,  p_L,  rho_R, u_R, p_R);

      for(int i=0; i<cells; i++){

        grid_P[i]=(h*i);
      x_over_t=(i*h-0.5)/t_total;
      auto [rho, u, p] = sample(p_star1, u_star1, rho_L, u_L, p_L, rho_R, u_R, p_R, x_over_t);
      rho1[i]=rho;
      p1[i]=p;
      u1[i]=u;
      I1[i]=p/(gamma-1)+0.5*rho*u*u;

    }
     makeGIFfile(cells,0,g,t_total, count, t_output, h, u1,p1, rho1, I1, t_total, B, true);     g++;

   }
if(m=="Fletcher"){
kabare(cells, ghost_cells,u,p,rho, I,t_total, tau, h,left_bc, right_bc);
  }
if(m=="customShem"){
feoktistov(cells, ghost_cells, u, p, rho, I, t_total, tau, h,
           left_bc, right_bc,
           high_shem, tipe_reconsrutin, limiter,
           use_diffusion, use_anti_diffusion, use_artificial_viscosity, viscosity_coeff, anti_diffusion_coeff, time_integrator, nameRS);
 }

    }
    std::vector<double> rhoP(cells);
    std::vector<double> uP(cells);
    std::vector<double> pP(cells);
    std::vector<double> IP(cells);
    for(int i=0; i<cells; i++){
  grid_P[i]=(h*i);
  rhoP[i]=rho[i+ghost_cells];
  pP[i]=p[i+ghost_cells];
  uP[i]=u[i+ghost_cells];
  IP[i]=(I[i+ghost_cells]-0.5*rhoP[i]*uP[i]*uP[i])/(rhoP[i]);

    }
    if(m=="godunov") outputCSV("data/hislennoereshG.csv", grid_P, uP, pP, rhoP, IP, t_end);
    if(m=="godunov-kolgan") outputCSV("data/hislennoereshGK.csv", grid_P, uP, pP, rhoP, IP, t_end);
    if(m=="godunov-kolgan-rodionov") outputCSV("data/hislennoereshGKR.csv", grid_P, uP, pP, rhoP, IP, t_end);
    if(m=="eno") outputCSV("data/hislennoereshENO.csv", grid_P, uP, pP, rhoP, IP, t_end);
    if(m=="weno") outputCSV("data/hislennoereshWENO.csv", grid_P, uP, pP, rhoP, IP, t_end);
    if(m=="Fletcher") outputCSV("data/hislennoereshKabare.csv", grid_P, uP, pP, rhoP, IP, t_end);
    if(m=="customShem") outputCSV("data/hislennoereshCustom.csv", grid_P, uP, pP, rhoP, IP, t_end);

    }
    
    double p_star=solve_p_star(rho_L, u_L, p_L, rho_R, u_R, p_R);
    double u_star1=u_star(p_star,  rho_L, u_L,  p_L,  rho_R, u_R, p_R);

    for(int i=0; i<cells; i++){

      grid_P[i]=(h*i);
      x_over_t=(i*h-0.5)/t_end;
      auto [rho, u, p] = sample(p_star, u_star1, rho_L, u_L, p_L, rho_R, u_R, p_R, x_over_t);
      rho1[i]=rho;
      p1[i]=p;
      u1[i]=u;
      I1[i]=(p/(rho*(gamma-1)));

    }
    outputCSV("data/tochnoeresh.csv", grid_P, u1,p1,rho1,I1, t_end);

    }
    else {
        std::cerr << "Error of reading mesh file" << std::endl;
        return 1;
    }


    return 0;
}