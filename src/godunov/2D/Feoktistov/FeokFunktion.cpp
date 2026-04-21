#include <iostream>
#include <cmath>
#include <algorithm>
#include <tuple>
#include <vector>
#include "Batarin/parameters.h"
#include "Medyakova/Output.h" 
using namespace std;

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

std::tuple<double,double,double> osher_solomon(
    double rhoL,double uL,double pL,
    double rhoR,double uR,double pR,
    double xi,double gamma){
    
    const double eps=1e-12;
    const double gamma_minus_1=gamma-1.0;
    
    double EL=pL/gamma_minus_1+0.5*rhoL*uL*uL;
    double ER=pR/gamma_minus_1+0.5*rhoR*uR*uR;
    
    std::array<double,3> UL={rhoL,rhoL*uL,EL};
    std::array<double,3> UR={rhoR,rhoR*uR,ER};
    
    std::array<double,3> FL={
        rhoL*uL,
        rhoL*uL*uL+pL,
        uL*(EL+pL)
    };
    
    std::array<double,3> FR={
        rhoR*uR,
        rhoR*uR*uR+pR,
        uR*(ER+pR)
    };
    
    constexpr int Q=3;
    const double theta[Q]={0.11270166537925831148,0.5,0.88729833462074168852};
    const double weight[Q]={5.0/18.0,4.0/9.0,5.0/18.0};
    
    //            |A| dU/d?
    std::array<double,3> integral={0.0,0.0,0.0};
    
    double drho=rhoR-rhoL;
    double du=uR-uL;
    double dp=pR-pL;
    
    for(int k=0;k<Q;++k){
        double psi=theta[k];//                   (0???1)
        
        double rho=rhoL+psi*drho;
        double u=uL+psi*du;
        double p=pL+psi*dp;
        
        std::array<double,3> dU_dpsi;
        dU_dpsi[0]=drho;
        dU_dpsi[1]=drho*u+rho*du;
        dU_dpsi[2]=dp/gamma_minus_1+0.5*drho*u*u+rho*u*du;
        
        
        double a=std::sqrt(gamma*p/std::max(rho,eps));//              
        double E=p/gamma_minus_1+0.5*rho*u*u;
        double H=(E+p)/rho;        
       
	double lambda1=u-a;
        double lambda2=u;
        double lambda3=u+a;
        
        double beta=(gamma-1.0)/(2.0*a*a);
        

        double l11=0.5*(beta*u*u+u/a);
        double l12=-0.5*(beta*u+1.0/a);
        double l13=0.5*beta;
        
        double l21=1.0-beta*u*u;
        double l22=beta*u;
        double l23=-beta;
        
        double l31=0.5*(beta*u*u-u/a);
        double l32=-0.5*(beta*u-1.0/a);
        double l33=0.5*beta;
        
        double alpha1=l11*dU_dpsi[0]+l12*dU_dpsi[1]+l13*dU_dpsi[2];
        double alpha2=l21*dU_dpsi[0]+l22*dU_dpsi[1]+l23*dU_dpsi[2];
        double alpha3=l31*dU_dpsi[0]+l32*dU_dpsi[1]+l33*dU_dpsi[2];
        
        alpha1*=std::abs(lambda1);
        alpha2*=std::abs(lambda2);
        alpha3*=std::abs(lambda3);
        
        //r1=[1,u-a,H-u*a]^T
        //r2=[1,u,0.5*u^2]^T
        //r3=[1,u+a,H+u*a]^T
        
     
        std::array<double,3> term;
        term[0]=alpha1+alpha2+alpha3;
        term[1]=(u-a)*alpha1+u*alpha2+(u+a)*alpha3;
        term[2]=(H-u*a)*alpha1+0.5*u*u*alpha2+(H+u*a)*alpha3;
        
        //                             
        integral[0]+=weight[k]*term[0];
        integral[1]+=weight[k]*term[1];
        integral[2]+=weight[k]*term[2];
    }
    
    //                               -        
    std::array<double,3> F_osher;
    for(int i=0;i<3;++i){
        F_osher[i]=0.5*(FL[i]+FR[i])-0.5*integral[i];
    }
    
    return std::make_tuple(F_osher[0],F_osher[1],F_osher[2]);
}

std::tuple<double,double,double> superSolve(double rhoL,double uL,double pL,
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
        //       HLL
        auto [rho_face, u_face, p_face] = hll(rhoL, uL, pL, rhoR, uR, pR, xi);
        return {rho_face, u_face, p_face};
    }
    else if (name == "hllc") {
        //       HLLC
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


std::tuple<
    std::vector<std::vector<double>>, 
    std::vector<std::vector<double>>, 
    std::vector<std::vector<double>>, 
    std::vector<std::vector<double>>,
    std::vector<std::vector<double>>,
    std::vector<std::vector<double>>,
    std::vector<std::vector<double>>,
    std::vector<std::vector<double>>
>
computeFluxes2D(
    int cells_x, int cells_y, int ghost_cells,
    const std::vector<std::vector<double>>& rho,
    const std::vector<std::vector<double>>& u,
    const std::vector<std::vector<double>>& v,
    const std::vector<std::vector<double>>& P,
    double hx, double hy,
    int Q,
    double gamma,
    const std::string& nameRS)
{
    int Nx = cells_x + 2 * ghost_cells;
    int Ny = cells_y + 2 * ghost_cells;

    std::vector<std::vector<double>> F1(Ny, std::vector<double>(Nx - 1, 0.0));
    std::vector<std::vector<double>> F2(Ny, std::vector<double>(Nx - 1, 0.0));
    std::vector<std::vector<double>> F3(Ny, std::vector<double>(Nx - 1, 0.0));
    std::vector<std::vector<double>> F4(Ny, std::vector<double>(Nx - 1, 0.0));

    std::vector<std::vector<double>> G1(Ny - 1, std::vector<double>(Nx, 0.0));
    std::vector<std::vector<double>> G2(Ny - 1, std::vector<double>(Nx, 0.0));
    std::vector<std::vector<double>> G3(Ny - 1, std::vector<double>(Nx, 0.0));
    std::vector<std::vector<double>> G4(Ny - 1, std::vector<double>(Nx, 0.0));

    if(Q==0){
    for (int iy = 0; iy < Ny; ++iy) {
        for (int ix = 0; ix < Nx - 1; ++ix) {
            double rho_L = rho[iy][ix];
            double u_L   = u[iy][ix];
            double v_L   = v[iy][ix];
            double p_L   = P[iy][ix];

            double rho_R = rho[iy][ix + 1];
            double u_R   = u[iy][ix + 1];
            double v_R   = v[iy][ix + 1];
            double p_R   = P[iy][ix + 1];

            auto [rho_face, u_face, p_face] = superSolve(rho_L, u_L, p_L,
                                                         rho_R, u_R, p_R,
                                                         nameRS);

            double v_star = (u_face >= 0.0) ? v_L : v_R;

            double E_face = p_face / (gamma - 1.0)
                          + 0.5 * rho_face * (u_face * u_face + v_star * v_star);
	    if(nameRS!="rusanov"){
            F1[iy][ix] = rho_face * u_face;
            F2[iy][ix] = rho_face * u_face * u_face + p_face;
            F3[iy][ix] = rho_face * u_face * v_star;
            F4[iy][ix] = u_face * (E_face + p_face);}
	    
	    else{
	    F1[iy][ix] = rho_face;
            F2[iy][ix] = u_face;
            F3[iy][ix] = rho_face * v_star;
            F4[iy][ix] = p_face;
	}
        }
    }

    for (int iy = 0; iy < Ny - 1; ++iy) {
        for (int ix = 0; ix < Nx; ++ix) {
            double rho_B = rho[iy][ix];
            double u_B   = u[iy][ix];
            double v_B   = v[iy][ix];
            double p_B   = P[iy][ix];

            double rho_T = rho[iy + 1][ix];
            double u_T   = u[iy + 1][ix];
            double v_T   = v[iy + 1][ix];
            double p_T   = P[iy + 1][ix];

            auto [rho_face, v_face, p_face] = superSolve(rho_B, v_B, p_B,
                                                         rho_T, v_T, p_T,
                                                         nameRS);

            double u_star = (v_face >= 0.0) ? u_B : u_T;

            double E_face = p_face / (gamma - 1.0)
                          + 0.5 * rho_face * (u_star * u_star + v_face * v_face);

            if(nameRS!="rusanov"){
	    G1[iy][ix] = rho_face * v_face;
            G2[iy][ix] = rho_face * v_face * u_star;
            G3[iy][ix] = rho_face * v_face * v_face + p_face;
            G4[iy][ix] = v_face * (E_face + p_face);
	   }
	else{
	    G1[iy][ix] = rho_face;
            G2[iy][ix] = rho_face * u_star;
            G3[iy][ix] =  v_face;
            G4[iy][ix] =  p_face;

	}

        }
    }
    }
   if (Q == 1) {
    for (int iy = 0; iy < Ny; ++iy) {
        for (int ix = 1; ix <= Nx-3; ++ix) {
            double drho_i  = minmod(rho[iy][ix] - rho[iy][ix-1], rho[iy][ix+1] - rho[iy][ix]);
            double du_i    = minmod(u[iy][ix]   - u[iy][ix-1],   u[iy][ix+1] - u[iy][ix]);
            double dP_i    = minmod(P[iy][ix]   - P[iy][ix-1],   P[iy][ix+1] - P[iy][ix]);

            double drho_ip1 = minmod(rho[iy][ix+1] - rho[iy][ix], rho[iy][ix+2] - rho[iy][ix+1]);
            double du_ip1   = minmod(u[iy][ix+1]   - u[iy][ix],   u[iy][ix+2] - u[iy][ix+1]);
            double dP_ip1   = minmod(P[iy][ix+1]   - P[iy][ix],   P[iy][ix+2] - P[iy][ix+1]);

            double rho_L = rho[iy][ix]   + 0.5 * drho_i;
            double u_L   = u[iy][ix]     + 0.5 * du_i;
            double p_L   = P[iy][ix]     + 0.5 * dP_i;
            double rho_R = rho[iy][ix+1] - 0.5 * drho_ip1;
            double u_R   = u[iy][ix+1]   - 0.5 * du_ip1;
            double p_R   = P[iy][ix+1]   - 0.5 * dP_ip1;

            rho_L = std::max(rho_L, 1e-10);
            p_L   = std::max(p_L,   1e-10);
            rho_R = std::max(rho_R, 1e-10);
            p_R   = std::max(p_R,   1e-10);

            double v_L = v[iy][ix];
            double v_R = v[iy][ix+1];

            auto [rho_face, u_face, p_face] = superSolve(rho_L, u_L, p_L,
                                                          rho_R, u_R, p_R,
                                                          nameRS);

            double v_star = (u_face >= 0.0) ? v_L : v_R;
            double E_face = p_face / (gamma - 1.0)
                          + 0.5 * rho_face * (u_face * u_face + v_star * v_star);

            F1[iy][ix] = rho_face * u_face;
            F2[iy][ix] = rho_face * u_face * u_face + p_face;
            F3[iy][ix] = rho_face * u_face * v_star;
            F4[iy][ix] = u_face * (E_face + p_face);
        }
        for (int ix : {0, Nx-2}) {
            if (ix < 0 || ix >= Nx-1) continue;
            double rho_L = rho[iy][ix];
            double u_L   = u[iy][ix];
            double v_L   = v[iy][ix];
            double p_L   = P[iy][ix];
            double rho_R = rho[iy][ix+1];
            double u_R   = u[iy][ix+1];
            double v_R   = v[iy][ix+1];
            double p_R   = P[iy][ix+1];

            auto [rho_face, u_face, p_face] = superSolve(rho_L, u_L, p_L,
                                                          rho_R, u_R, p_R,
                                                          nameRS);
            double v_star = (u_face >= 0.0) ? v_L : v_R;
            double E_face = p_face / (gamma - 1.0)
                          + 0.5 * rho_face * (u_face * u_face + v_star * v_star);
            F1[iy][ix] = rho_face * u_face;
            F2[iy][ix] = rho_face * u_face * u_face + p_face;
            F3[iy][ix] = rho_face * u_face * v_star;
            F4[iy][ix] = u_face * (E_face + p_face);
        }
    }

    for (int iy = 1; iy <= Ny-3; ++iy) {
        for (int ix = 0; ix < Nx; ++ix) {
            double drho_i  = minmod(rho[iy][ix] - rho[iy-1][ix], rho[iy+1][ix] - rho[iy][ix]);
            double dv_i    = minmod(v[iy][ix]   - v[iy-1][ix],   v[iy+1][ix] - v[iy][ix]);
            double dP_i    = minmod(P[iy][ix]   - P[iy-1][ix],   P[iy+1][ix] - P[iy][ix]);

            double drho_ip1 = minmod(rho[iy+1][ix] - rho[iy][ix], rho[iy+2][ix] - rho[iy+1][ix]);
            double dv_ip1   = minmod(v[iy+1][ix]   - v[iy][ix],   v[iy+2][ix] - v[iy+1][ix]);
            double dP_ip1   = minmod(P[iy+1][ix]   - P[iy][ix],   P[iy+2][ix] - P[iy+1][ix]);

            double rho_B = rho[iy][ix]   + 0.5 * drho_i;
            double v_B   = v[iy][ix]     + 0.5 * dv_i;
            double p_B   = P[iy][ix]     + 0.5 * dP_i;
            double rho_T = rho[iy+1][ix] - 0.5 * drho_ip1;
            double v_T   = v[iy+1][ix]   - 0.5 * dv_ip1;
            double p_T   = P[iy+1][ix]   - 0.5 * dP_ip1;

            rho_B = std::max(rho_B, 1e-10);
            p_B   = std::max(p_B,   1e-10);
            rho_T = std::max(rho_T, 1e-10);
            p_T   = std::max(p_T,   1e-10);

            double u_B = u[iy][ix];
            double u_T = u[iy+1][ix];

            auto [rho_face, v_face, p_face] = superSolve(rho_B, v_B, p_B,
                                                          rho_T, v_T, p_T,
                                                          nameRS);

            double u_star = (v_face >= 0.0) ? u_B : u_T;
            double E_face = p_face / (gamma - 1.0)
                          + 0.5 * rho_face * (u_star * u_star + v_face * v_face);

            G1[iy][ix] = rho_face * v_face;
            G2[iy][ix] = rho_face * v_face * u_star;
            G3[iy][ix] = rho_face * v_face * v_face + p_face;
            G4[iy][ix] = v_face * (E_face + p_face);
        }
    }
    for (int iy : {0, Ny-2}) {
        if (iy < 0 || iy >= Ny-1) continue;
        for (int ix = 0; ix < Nx; ++ix) {
            double rho_B = rho[iy][ix];
            double v_B   = v[iy][ix];
            double p_B   = P[iy][ix];
            double rho_T = rho[iy+1][ix];
            double v_T   = v[iy+1][ix];
            double p_T   = P[iy+1][ix];
            double u_B   = u[iy][ix];
            double u_T   = u[iy+1][ix];

            auto [rho_face, v_face, p_face] = superSolve(rho_B, v_B, p_B,
                                                          rho_T, v_T, p_T,
                                                          nameRS);
            double u_star = (v_face >= 0.0) ? u_B : u_T;
            double E_face = p_face / (gamma - 1.0)
                          + 0.5 * rho_face * (u_star * u_star + v_face * v_face);
            G1[iy][ix] = rho_face * v_face;
            G2[iy][ix] = rho_face * v_face * u_star;
            G3[iy][ix] = rho_face * v_face * v_face + p_face;
            G4[iy][ix] = v_face * (E_face + p_face);
        }
    }
}

    return {F1, F2, F3, F4, G1, G2, G3, G4};
}
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
    int Q, const std::string& time_integrator, const std::string& nameRS)
{
    const double gamma = 1.4;
    const double eps_rho = 1e-5;
    const double eps_e   = 1e-6;

    int Nx = cells_x + 2 * ghost_cells;
    int Ny = cells_y + 2 * ghost_cells;

    std::vector<std::vector<double>> U1(Ny, std::vector<double>(Nx));
    std::vector<std::vector<double>> U2(Ny, std::vector<double>(Nx));
    std::vector<std::vector<double>> U3(Ny, std::vector<double>(Nx));
    std::vector<std::vector<double>> U4(Ny, std::vector<double>(Nx));

    std::vector<std::vector<double>> U1_new(Ny, std::vector<double>(Nx));
    std::vector<std::vector<double>> U2_new(Ny, std::vector<double>(Nx));
    std::vector<std::vector<double>> U3_new(Ny, std::vector<double>(Nx));
    std::vector<std::vector<double>> U4_new(Ny, std::vector<double>(Nx));

    for (int j = 0; j < Ny; ++j) {
        for (int i = 0; i < Nx; ++i) {
            U1[j][i] = rho[j][i];
            U2[j][i] = rho[j][i] * u[j][i];
            U3[j][i] = rho[j][i] * v[j][i];
            U4[j][i] = P[j][i] / (gamma - 1.0) + 0.5 * rho[j][i] * (u[j][i]*u[j][i] + v[j][i]*v[j][i]);
        }
    }

    for (int j = 0; j < Ny; ++j) {
        applyBoundaryConditions(rho[j].data(), u[j].data(), P[j].data(), cells_x, ghost_cells, left_bc, right_bc);
            }

    for (int i = 0; i < Nx; ++i) {
        std::vector<double> rho_col(Ny), u_col(Ny), v_col(Ny), P_col(Ny);
        for (int j = 0; j < Ny; ++j) {
            rho_col[j] = rho[j][i];
            u_col[j]   = u[j][i];
            v_col[j]   = v[j][i];
            P_col[j]   = P[j][i];
        }

        applyBoundaryConditions(rho_col.data(), u_col.data(), P_col.data(), cells_y, ghost_cells, bottom_bc, top_bc);

        
        for (int j = 0; j < Ny; ++j) {
            rho[j][i] = rho_col[j];
            u[j][i]   = u_col[j];
            v[j][i]   = v_col[j];
            P[j][i]   = P_col[j];
        }
    }

    auto [F1, F2, F3, F4, G1, G2, G3, G4] = computeFluxes2D(
        cells_x, cells_y, ghost_cells,
        rho, u, v, P, hx, hy, Q, gamma, nameRS);

    if (time_integrator == "euler") {
        for (int j = ghost_cells; j < Ny - ghost_cells; ++j) {
            for (int i = ghost_cells; i < Nx - ghost_cells; ++i) {
                double dFx1 = F1[j][i] - F1[j][i-1];
                double dFx2 = F2[j][i] - F2[j][i-1];
                double dFx3 = F3[j][i] - F3[j][i-1];
                double dFx4 = F4[j][i] - F4[j][i-1];

                double dGy1 = G1[j][i] - G1[j-1][i];
                double dGy2 = G2[j][i] - G2[j-1][i];
                double dGy3 = G3[j][i] - G3[j-1][i];
                double dGy4 = G4[j][i] - G4[j-1][i];

                U1_new[j][i] = U1[j][i] - tau/hx * dFx1 - tau/hy * dGy1;
                U2_new[j][i] = U2[j][i] - tau/hx * dFx2 - tau/hy * dGy2;
                U3_new[j][i] = U3[j][i] - tau/hx * dFx3 - tau/hy * dGy3;
                U4_new[j][i] = U4[j][i] - tau/hx * dFx4 - tau/hy * dGy4;

                U1_new[j][i] = std::max(U1_new[j][i], eps_rho);
            }
        }

        for (int j = ghost_cells; j < Ny - ghost_cells; ++j) {
            for (int i = ghost_cells; i < Nx - ghost_cells; ++i) {
                rho[j][i] = U1_new[j][i];
                u[j][i] = U2_new[j][i] / rho[j][i];
                v[j][i] = U3_new[j][i] / rho[j][i];
                double Ek = 0.5 * rho[j][i] * (u[j][i]*u[j][i] + v[j][i]*v[j][i]);
                double e_in = (U4_new[j][i] - Ek) / rho[j][i];
                if (e_in < eps_e) e_in = eps_e;
                P[j][i] = e_in * rho[j][i] * (gamma - 1.0);
                I[j][i] = e_in;
            }
        }

        U1.swap(U1_new);
        U2.swap(U2_new);
        U3.swap(U3_new);
        U4.swap(U4_new);
    }
    else {
        std::cerr << "Time integrator " << time_integrator << " not implemented.\n";
    }

    t_total += tau;
    std::cout << "t = " << t_total << ", tau = " << tau
              << ", method = " << time_integrator
              << ", RS = " << nameRS << std::endl;
}

std::tuple<
    std::vector<std::vector<double>>, // drho_x
    std::vector<std::vector<double>>, // du_x
    std::vector<std::vector<double>>, // dv_x
    std::vector<std::vector<double>>, // dP_x
    std::vector<std::vector<double>>, // drho_y
    std::vector<std::vector<double>>, // du_y
    std::vector<std::vector<double>>, // dv_y
    std::vector<std::vector<double>>  // dP_y
>
computeSlopes2D(
    int cells_x, int cells_y, int ghost_cells,
    const std::vector<std::vector<double>>& rho,
    const std::vector<std::vector<double>>& u,
    const std::vector<std::vector<double>>& v,
    const std::vector<std::vector<double>>& P)
{
    int Nx = cells_x + 2 * ghost_cells;
    int Ny = cells_y + 2 * ghost_cells;

    std::vector<std::vector<double>> drho_x(Ny, std::vector<double>(Nx, 0.0));
    std::vector<std::vector<double>> du_x(Ny, std::vector<double>(Nx, 0.0));
    std::vector<std::vector<double>> dv_x(Ny, std::vector<double>(Nx, 0.0));
    std::vector<std::vector<double>> dP_x(Ny, std::vector<double>(Nx, 0.0));
    std::vector<std::vector<double>> drho_y(Ny, std::vector<double>(Nx, 0.0));
    std::vector<std::vector<double>> du_y(Ny, std::vector<double>(Nx, 0.0));
    std::vector<std::vector<double>> dv_y(Ny, std::vector<double>(Nx, 0.0));
    std::vector<std::vector<double>> dP_y(Ny, std::vector<double>(Nx, 0.0));

    for (int j = 0; j < Ny; ++j) {
        for (int i = ghost_cells; i < Nx - ghost_cells; ++i) {
            drho_x[j][i] = minmod(rho[j][i] - rho[j][i-1], rho[j][i+1] - rho[j][i]);
            du_x[j][i]   = minmod(u[j][i] - u[j][i-1],     u[j][i+1] - u[j][i]);
            dv_x[j][i]   = minmod(v[j][i] - v[j][i-1],     v[j][i+1] - v[j][i]);
            dP_x[j][i]   = minmod(P[j][i] - P[j][i-1],     P[j][i+1] - P[j][i]);
        }
    }

    for (int j = ghost_cells; j < Ny - ghost_cells; ++j) {
        for (int i = 0; i < Nx; ++i) {
            drho_y[j][i] = minmod(rho[j][i] - rho[j-1][i], rho[j+1][i] - rho[j][i]);
            du_y[j][i]   = minmod(u[j][i] - u[j-1][i],     u[j+1][i] - u[j][i]);
            dv_y[j][i]   = minmod(v[j][i] - v[j-1][i],     v[j+1][i] - v[j][i]);
            dP_y[j][i]   = minmod(P[j][i] - P[j-1][i],     P[j+1][i] - P[j][i]);
        }
    }

    return {drho_x, du_x, dv_x, dP_x, drho_y, du_y, dv_y, dP_y};
}

std::tuple<
    std::vector<std::vector<double>>, // F1 (ernnr, x-dinie)
    std::vector<std::vector<double>>, // F2 (x-cedoeun, x-dinie)
    std::vector<std::vector<double>>, // F3 (y-cedoeun, x-dinie)
    std::vector<std::vector<double>>, // F4 (yildac?, x-dinie)
    std::vector<std::vector<double>>, // G1 (ernnr, y-dinie)
    std::vector<std::vector<double>>, // G2 (x-cedoeun, y-dinie)
    std::vector<std::vector<double>>, // G3 (y-cedoeun, y-dinie)
    std::vector<std::vector<double>>  // G4 (yildac?, y-dinie)
>
computeFluxesWithGivenSlopes(
    int cells_x, int cells_y, int ghost_cells,
    const std::vector<std::vector<double>>& rho,
    const std::vector<std::vector<double>>& u,
    const std::vector<std::vector<double>>& v,
    const std::vector<std::vector<double>>& P,
    const std::vector<std::vector<double>>& drho_x,
    const std::vector<std::vector<double>>& du_x,
    const std::vector<std::vector<double>>& dv_x,
    const std::vector<std::vector<double>>& dP_x,
    const std::vector<std::vector<double>>& drho_y,
    const std::vector<std::vector<double>>& du_y,
    const std::vector<std::vector<double>>& dv_y,
    const std::vector<std::vector<double>>& dP_y,
    double gamma,
    const std::string& nameRS)
{
    int Nx = cells_x + 2 * ghost_cells;
    int Ny = cells_y + 2 * ghost_cells;

    std::vector<std::vector<double>> F1(Ny, std::vector<double>(Nx - 1, 0.0));
    std::vector<std::vector<double>> F2(Ny, std::vector<double>(Nx - 1, 0.0));
    std::vector<std::vector<double>> F3(Ny, std::vector<double>(Nx - 1, 0.0));
    std::vector<std::vector<double>> F4(Ny, std::vector<double>(Nx - 1, 0.0));

    std::vector<std::vector<double>> G1(Ny - 1, std::vector<double>(Nx, 0.0));
    std::vector<std::vector<double>> G2(Ny - 1, std::vector<double>(Nx, 0.0));
    std::vector<std::vector<double>> G3(Ny - 1, std::vector<double>(Nx, 0.0));
    std::vector<std::vector<double>> G4(Ny - 1, std::vector<double>(Nx, 0.0));

    for (int j = 0; j < Ny; ++j) {
        for (int i = 0; i < Nx - 1; ++i) {
            double rho_L = rho[j][i] + 0.5 * drho_x[j][i];
            double u_L   = u[j][i]   + 0.5 * du_x[j][i];
            double v_L   = v[j][i]   + 0.5 * dv_x[j][i];
            double p_L   = P[j][i]   + 0.5 * dP_x[j][i];

            double rho_R = rho[j][i+1] - 0.5 * drho_x[j][i+1];
            double u_R   = u[j][i+1]   - 0.5 * du_x[j][i+1];
            double v_R   = v[j][i+1]   - 0.5 * dv_x[j][i+1];
            double p_R   = P[j][i+1]   - 0.5 * dP_x[j][i+1];

            auto [rho_face, u_face, p_face] = superSolve(rho_L, u_L, p_L,
                                                         rho_R, u_R, p_R,
                                                         nameRS);

            double v_star = (u_face >= 0.0) ? v_L : v_R;
            double E_face = p_face / (gamma - 1.0)
                          + 0.5 * rho_face * (u_face * u_face + v_star * v_star);

            F1[j][i] = rho_face * u_face;
            F2[j][i] = rho_face * u_face * u_face + p_face;
            F3[j][i] = rho_face * u_face * v_star;
            F4[j][i] = u_face * (E_face + p_face);
        }
    }

    for (int j = 0; j < Ny - 1; ++j) {
        for (int i = 0; i < Nx; ++i) {
            double rho_B = rho[j][i] + 0.5 * drho_y[j][i];
            double u_B   = u[j][i]   + 0.5 * du_y[j][i];
            double v_B   = v[j][i]   + 0.5 * dv_y[j][i];
            double p_B   = P[j][i]   + 0.5 * dP_y[j][i];

            double rho_T = rho[j+1][i] - 0.5 * drho_y[j+1][i];
            double u_T   = u[j+1][i]   - 0.5 * du_y[j+1][i];
            double v_T   = v[j+1][i]   - 0.5 * dv_y[j+1][i];
            double p_T   = P[j+1][i]   - 0.5 * dP_y[j+1][i];

            auto [rho_face, v_face, p_face] = superSolve(rho_B, v_B, p_B,
                                                         rho_T, v_T, p_T,
                                                         nameRS);

            double u_star = (v_face >= 0.0) ? u_B : u_T;
            double E_face = p_face / (gamma - 1.0)
                          + 0.5 * rho_face * (u_star * u_star + v_face * v_face);

            G1[j][i] = rho_face * v_face;
            G2[j][i] = rho_face * v_face * u_star;
            G3[j][i] = rho_face * v_face * v_face + p_face;
            G4[j][i] = v_face * (E_face + p_face);
        }
    }

    return {F1, F2, F3, F4, G1, G2, G3, G4};
}

std::tuple<
    std::vector<std::vector<double>>, // rhoHalf
    std::vector<std::vector<double>>, // uHalf
    std::vector<std::vector<double>>, // vHalf
    std::vector<std::vector<double>>, // pHalf
    std::vector<std::vector<double>>, // drho_x
    std::vector<std::vector<double>>, // du_x
    std::vector<std::vector<double>>, // dv_x
    std::vector<std::vector<double>>, // dP_x
    std::vector<std::vector<double>>, // drho_y
    std::vector<std::vector<double>>, // du_y
    std::vector<std::vector<double>>, // dv_y
    std::vector<std::vector<double>>  // dP_y
>
predictorStep2D(
    int cells_x, int cells_y, int ghost_cells,
    const std::vector<std::vector<double>>& rho,
    const std::vector<std::vector<double>>& u,
    const std::vector<std::vector<double>>& v,
    const std::vector<std::vector<double>>& P,
    const std::vector<std::vector<double>>& U1,
    const std::vector<std::vector<double>>& U2,
    const std::vector<std::vector<double>>& U3,
    const std::vector<std::vector<double>>& U4,
    double tau, double hx, double hy,
    double gamma, const std::string& nameRS)
{
    int Nx = cells_x + 2 * ghost_cells;
    int Ny = cells_y + 2 * ghost_cells;

    auto [drho_x, du_x, dv_x, dP_x, drho_y, du_y, dv_y, dP_y] =
        computeSlopes2D(cells_x, cells_y, ghost_cells, rho, u, v, P);

    auto [F1, F2, F3, F4, G1, G2, G3, G4] =
        computeFluxesWithGivenSlopes(cells_x, cells_y, ghost_cells,
                                     rho, u, v, P,
                                     drho_x, du_x, dv_x, dP_x,
                                     drho_y, du_y, dv_y, dP_y,
                                     gamma, nameRS);

    std::vector<std::vector<double>> U1_temp(Ny, std::vector<double>(Nx, 0.0));
    std::vector<std::vector<double>> U2_temp(Ny, std::vector<double>(Nx, 0.0));
    std::vector<std::vector<double>> U3_temp(Ny, std::vector<double>(Nx, 0.0));
    std::vector<std::vector<double>> U4_temp(Ny, std::vector<double>(Nx, 0.0));

    for (int j = ghost_cells; j < Ny - ghost_cells; ++j) {
        for (int i = ghost_cells; i < Nx - ghost_cells; ++i) {
            double dFx1 = F1[j][i] - F1[j][i-1];
            double dFx2 = F2[j][i] - F2[j][i-1];
            double dFx3 = F3[j][i] - F3[j][i-1];
            double dFx4 = F4[j][i] - F4[j][i-1];

            double dGy1 = G1[j][i] - G1[j-1][i];
            double dGy2 = G2[j][i] - G2[j-1][i];
            double dGy3 = G3[j][i] - G3[j-1][i];
            double dGy4 = G4[j][i] - G4[j-1][i];

            U1_temp[j][i] = U1[j][i] - 0.5 * tau / hx * dFx1 - 0.5 * tau / hy * dGy1;
            U2_temp[j][i] = U2[j][i] - 0.5 * tau / hx * dFx2 - 0.5 * tau / hy * dGy2;
            U3_temp[j][i] = U3[j][i] - 0.5 * tau / hx * dFx3 - 0.5 * tau / hy * dGy3;
            U4_temp[j][i] = U4[j][i] - 0.5 * tau / hx * dFx4 - 0.5 * tau / hy * dGy4;
        }
    }

    const double eps_rho = 1e-5;
    const double eps_e   = 1e-6;

    std::vector<std::vector<double>> rhoHalf(Ny, std::vector<double>(Nx, 0.0));
    std::vector<std::vector<double>> uHalf(Ny, std::vector<double>(Nx, 0.0));
    std::vector<std::vector<double>> vHalf(Ny, std::vector<double>(Nx, 0.0));
    std::vector<std::vector<double>> pHalf(Ny, std::vector<double>(Nx, 0.0));

    for (int j = ghost_cells; j < Ny - ghost_cells; ++j) {
        for (int i = ghost_cells; i < Nx - ghost_cells; ++i) {
            rhoHalf[j][i] = std::max(U1_temp[j][i], eps_rho);
            uHalf[j][i]   = U2_temp[j][i] / rhoHalf[j][i];
            vHalf[j][i]   = U3_temp[j][i] / rhoHalf[j][i];
            double Ek = 0.5 * rhoHalf[j][i] * (uHalf[j][i]*uHalf[j][i] + vHalf[j][i]*vHalf[j][i]);
            double e_in = (U4_temp[j][i] - Ek) / rhoHalf[j][i];
            if (e_in < eps_e) e_in = eps_e;
            pHalf[j][i] = e_in * rhoHalf[j][i] * (gamma - 1.0);
        }
    }

    return {rhoHalf, uHalf, vHalf, pHalf,
            drho_x, du_x, dv_x, dP_x,
            drho_y, du_y, dv_y, dP_y};
}

std::tuple<
    std::vector<std::vector<double>>,
    std::vector<std::vector<double>>,
    std::vector<std::vector<double>>,
    std::vector<std::vector<double>>,
    std::vector<std::vector<double>>,
    std::vector<std::vector<double>>,
    std::vector<std::vector<double>>,
    std::vector<std::vector<double>>>
correctorStep2D(
    int cells_x, int cells_y, int ghost_cells,
    const std::vector<std::vector<double>>& rhoHalf,
    const std::vector<std::vector<double>>& uHalf,
    const std::vector<std::vector<double>>& vHalf,
    const std::vector<std::vector<double>>& pHalf,
    const std::vector<std::vector<double>>& drho_x,
    const std::vector<std::vector<double>>& du_x,
    const std::vector<std::vector<double>>& dv_x,
    const std::vector<std::vector<double>>& dP_x,
    const std::vector<std::vector<double>>& drho_y,
    const std::vector<std::vector<double>>& du_y,
    const std::vector<std::vector<double>>& dv_y,
    const std::vector<std::vector<double>>& dP_y,
    double gamma, const std::string& nameRS)
{
    return computeFluxesWithGivenSlopes(cells_x, cells_y, ghost_cells,
                                        rhoHalf, uHalf, vHalf, pHalf,
                                        drho_x, du_x, dv_x, dP_x,
                                        drho_y, du_y, dv_y, dP_y,
                                        gamma, nameRS);
}
void applyBCs2D(
    std::vector<std::vector<double>>& rho,
    std::vector<std::vector<double>>& u,
    std::vector<std::vector<double>>& v,
    std::vector<std::vector<double>>& P,
    int Nx, int Ny,
    int cells_x, int cells_y, int ghost_cells,
    const std::string& left_bc, const std::string& right_bc,
    const std::string& bottom_bc, const std::string& top_bc)
{
    // Ddceli?le oneiac? di x (nndiec)
    for (int j = 0; j < Ny; ++j) {
        applyBoundaryConditions(rho[j].data(), u[j].data(), P[j].data(),
                                cells_x, ghost_cells, left_bc, right_bc);
    }
    // Ddceli?le oneiac? di y (nnieaou)
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
}
void RodionovSolve2D(
    int cells_x, int cells_y, int ghost_cells,
    std::vector<std::vector<double>>& u,
    std::vector<std::vector<double>>& v,
    std::vector<std::vector<double>>& P,
    std::vector<std::vector<double>>& rho,
    std::vector<std::vector<double>>& I,
    double& t_total, double tau, double hx, double hy,
    const std::string& left_bc,   const std::string& right_bc,
    const std::string& bottom_bc, const std::string& top_bc,
    const std::string& time_integrator, const std::string& nameRS)
{
    const double gamma = 1.4;
    const double eps_rho = 1e-5;
    const double eps_e   = 1e-6;

    int Nx = cells_x + 2 * ghost_cells;
    int Ny = cells_y + 2 * ghost_cells;

    std::vector<std::vector<double>> U1(Ny, std::vector<double>(Nx));
    std::vector<std::vector<double>> U2(Ny, std::vector<double>(Nx));
    std::vector<std::vector<double>> U3(Ny, std::vector<double>(Nx));
    std::vector<std::vector<double>> U4(Ny, std::vector<double>(Nx));

    for (int j = 0; j < Ny; ++j) {
        for (int i = 0; i < Nx; ++i) {
            U1[j][i] = rho[j][i];
            U2[j][i] = rho[j][i] * u[j][i];
            U3[j][i] = rho[j][i] * v[j][i];
            U4[j][i] = P[j][i] / (gamma - 1.0)
                     + 0.5 * rho[j][i] * (u[j][i]*u[j][i] + v[j][i]*v[j][i]);
        }
    }

    applyBCs2D(rho, u, v, P, Nx, Ny, cells_x, cells_y, ghost_cells,
               left_bc, right_bc, bottom_bc, top_bc);

     if (time_integrator == "euler") {
        auto [rhoHalf, uHalf, vHalf, pHalf,
              drho_x, du_x, dv_x, dP_x,
              drho_y, du_y, dv_y, dP_y] =
            predictorStep2D(cells_x, cells_y, ghost_cells,
                            rho, u, v, P, U1, U2, U3, U4,
                            tau, hx, hy, gamma, nameRS); 

        applyBCs2D(rhoHalf, uHalf, vHalf, pHalf, Nx, Ny, cells_x, cells_y, ghost_cells,
                   left_bc, right_bc, bottom_bc, top_bc);

        auto [F1_half, F2_half, F3_half, F4_half,
              G1_half, G2_half, G3_half, G4_half] =
            correctorStep2D(cells_x, cells_y, ghost_cells,
                            rhoHalf, uHalf, vHalf, pHalf,
                            drho_x, du_x, dv_x, dP_x,
                            drho_y, du_y, dv_y, dP_y,
                            gamma, nameRS);

        std::vector<std::vector<double>> U1_new(Ny, std::vector<double>(Nx, 0.0));
        std::vector<std::vector<double>> U2_new(Ny, std::vector<double>(Nx, 0.0));
        std::vector<std::vector<double>> U3_new(Ny, std::vector<double>(Nx, 0.0));
        std::vector<std::vector<double>> U4_new(Ny, std::vector<double>(Nx, 0.0));

        for (int j = ghost_cells; j < Ny - ghost_cells; ++j) {
            for (int i = ghost_cells; i < Nx - ghost_cells; ++i) {
                double dFx1 = F1_half[j][i] - F1_half[j][i-1];
                double dFx2 = F2_half[j][i] - F2_half[j][i-1];
                double dFx3 = F3_half[j][i] - F3_half[j][i-1];
                double dFx4 = F4_half[j][i] - F4_half[j][i-1];

                double dGy1 = G1_half[j][i] - G1_half[j-1][i];
                double dGy2 = G2_half[j][i] - G2_half[j-1][i];
                double dGy3 = G3_half[j][i] - G3_half[j-1][i];
                double dGy4 = G4_half[j][i] - G4_half[j-1][i];

                U1_new[j][i] = U1[j][i] - tau/hx * dFx1 - tau/hy * dGy1;
                U2_new[j][i] = U2[j][i] - tau/hx * dFx2 - tau/hy * dGy2;
                U3_new[j][i] = U3[j][i] - tau/hx * dFx3 - tau/hy * dGy3;
                U4_new[j][i] = U4[j][i] - tau/hx * dFx4 - tau/hy * dGy4;

                U1_new[j][i] = std::max(U1_new[j][i], eps_rho);
            }
        }

        for (int j = ghost_cells; j < Ny - ghost_cells; ++j) {
            for (int i = ghost_cells; i < Nx - ghost_cells; ++i) {
                rho[j][i] = U1_new[j][i];
                u[j][i]   = U2_new[j][i] / rho[j][i];
                v[j][i]   = U3_new[j][i] / rho[j][i];
                double Ek = 0.5 * rho[j][i] * (u[j][i]*u[j][i] + v[j][i]*v[j][i]);
                double e_in = (U4_new[j][i] - Ek) / rho[j][i];
                if (e_in < eps_e) e_in = eps_e;
                P[j][i] = e_in * rho[j][i] * (gamma - 1.0);
                I[j][i] = U4_new[j][i];
            }
        }
        applyBCs2D(rho, u, v, P, Nx, Ny, cells_x, cells_y, ghost_cells,
                   left_bc, right_bc, bottom_bc, top_bc);
    }
    else {
        std::cerr << "Time integrator " << time_integrator << " not implemented in RodionovSolve2D.\n";
        return;
    }

    t_total += tau;
    std::cout << "t = " << t_total << ", tau = " << tau
              << ", method = " << time_integrator
              << ", RS = " << nameRS << std::endl;
}

std::tuple<double, double, double, double> ENOreconstruction_x_right(
    int i, int j,
    const std::vector<std::vector<double>>& u,
    const std::vector<std::vector<double>>& v,
    const std::vector<std::vector<double>>& P,
    const std::vector<std::vector<double>>& rho)
{
    double Su_L = u[j][i] - 2*u[j][i-1] + u[j][i-2];
    double Su_M = u[j][i+1] - 2*u[j][i] + u[j][i-1];
    double Su_R = u[j][i+2] - 2*u[j][i+1] + u[j][i];
    double ue;
    if (std::abs(Su_L) < std::abs(Su_M) && std::abs(Su_L) < std::abs(Su_R))
        ue = (1.0/3.0)*u[j][i-2] - (7.0/6.0)*u[j][i-1] + (11.0/6.0)*u[j][i];
    else if (std::abs(Su_R) < std::abs(Su_L) && std::abs(Su_R) < std::abs(Su_M))
        ue = (1.0/3.0)*u[j][i] + (5.0/6.0)*u[j][i+1] - (1.0/6.0)*u[j][i+2];
    else
        ue = (-1.0/6.0)*u[j][i-1] + (5.0/6.0)*u[j][i] + (1.0/3.0)*u[j][i+1];

    double Sv_L = v[j][i] - 2*v[j][i-1] + v[j][i-2];
    double Sv_M = v[j][i+1] - 2*v[j][i] + v[j][i-1];
    double Sv_R = v[j][i+2] - 2*v[j][i+1] + v[j][i];
    double ve;
    if (std::abs(Sv_L) < std::abs(Sv_M) && std::abs(Sv_L) < std::abs(Sv_R))
        ve = (1.0/3.0)*v[j][i-2] - (7.0/6.0)*v[j][i-1] + (11.0/6.0)*v[j][i];
    else if (std::abs(Sv_R) < std::abs(Sv_L) && std::abs(Sv_R) < std::abs(Sv_M))
        ve = (1.0/3.0)*v[j][i] + (5.0/6.0)*v[j][i+1] - (1.0/6.0)*v[j][i+2];
    else
        ve = (-1.0/6.0)*v[j][i-1] + (5.0/6.0)*v[j][i] + (1.0/3.0)*v[j][i+1];

    double Sp_L = P[j][i] - 2*P[j][i-1] + P[j][i-2];
    double Sp_M = P[j][i+1] - 2*P[j][i] + P[j][i-1];
    double Sp_R = P[j][i+2] - 2*P[j][i+1] + P[j][i];
    double Pe;
    if (std::abs(Sp_L) < std::abs(Sp_M) && std::abs(Sp_L) < std::abs(Sp_R))
        Pe = (1.0/3.0)*P[j][i-2] - (7.0/6.0)*P[j][i-1] + (11.0/6.0)*P[j][i];
    else if (std::abs(Sp_R) < std::abs(Sp_L) && std::abs(Sp_R) < std::abs(Sp_M))
        Pe = (1.0/3.0)*P[j][i] + (5.0/6.0)*P[j][i+1] - (1.0/6.0)*P[j][i+2];
    else
        Pe = (-1.0/6.0)*P[j][i-1] + (5.0/6.0)*P[j][i] + (1.0/3.0)*P[j][i+1];
    if (Pe < 1e-5) Pe = 1e-5;

    double Srho_L = rho[j][i] - 2*rho[j][i-1] + rho[j][i-2];
    double Srho_M = rho[j][i+1] - 2*rho[j][i] + rho[j][i-1];
    double Srho_R = rho[j][i+2] - 2*rho[j][i+1] + rho[j][i];
    double rhoe;
    if (std::abs(Srho_L) < std::abs(Srho_M) && std::abs(Srho_L) < std::abs(Srho_R))
        rhoe = (1.0/3.0)*rho[j][i-2] - (7.0/6.0)*rho[j][i-1] + (11.0/6.0)*rho[j][i];
    else if (std::abs(Srho_R) < std::abs(Srho_L) && std::abs(Srho_R) < std::abs(Srho_M))
        rhoe = (1.0/3.0)*rho[j][i] + (5.0/6.0)*rho[j][i+1] - (1.0/6.0)*rho[j][i+2];
    else
        rhoe = (-1.0/6.0)*rho[j][i-1] + (5.0/6.0)*rho[j][i] + (1.0/3.0)*rho[j][i+1];
    if (rhoe < 1e-4) rhoe = 1e-4;

    return std::make_tuple(ue, ve, Pe, rhoe);
}
std::tuple<double, double, double, double> ENOreconstruction_x_left(
    int i, int j,
    const std::vector<std::vector<double>>& u,
    const std::vector<std::vector<double>>& v,
    const std::vector<std::vector<double>>& P,
    const std::vector<std::vector<double>>& rho)
{
    // u
    double Su_L = u[j][i] - 2*u[j][i-1] + u[j][i-2];
    double Su_M = u[j][i+1] - 2*u[j][i] + u[j][i-1];
    double Su_R = u[j][i+2] - 2*u[j][i+1] + u[j][i];
    double ue;
    if (std::abs(Su_L) < std::abs(Su_M) && std::abs(Su_L) < std::abs(Su_R))
        ue = (-1.0/8.0)*u[j][i-2] + (3.0/4.0)*u[j][i-1] + (3.0/8.0)*u[j][i];
    else if (std::abs(Su_R) < std::abs(Su_L) && std::abs(Su_R) < std::abs(Su_M))
        ue = (15.0/8.0)*u[j][i] - (5.0/4.0)*u[j][i+1] + (3.0/8.0)*u[j][i+2];
    else
        ue = (3.0/8.0)*u[j][i-1] + (3.0/4.0)*u[j][i] - (1.0/8.0)*u[j][i+1];

    // v
    double Sv_L = v[j][i] - 2*v[j][i-1] + v[j][i-2];
    double Sv_M = v[j][i+1] - 2*v[j][i] + v[j][i-1];
    double Sv_R = v[j][i+2] - 2*v[j][i+1] + v[j][i];
    double ve;
    if (std::abs(Sv_L) < std::abs(Sv_M) && std::abs(Sv_L) < std::abs(Sv_R))
        ve = (-1.0/8.0)*v[j][i-2] + (3.0/4.0)*v[j][i-1] + (3.0/8.0)*v[j][i];
    else if (std::abs(Sv_R) < std::abs(Sv_L) && std::abs(Sv_R) < std::abs(Sv_M))
        ve = (15.0/8.0)*v[j][i] - (5.0/4.0)*v[j][i+1] + (3.0/8.0)*v[j][i+2];
    else
        ve = (3.0/8.0)*v[j][i-1] + (3.0/4.0)*v[j][i] - (1.0/8.0)*v[j][i+1];

    // P
    double Sp_L = P[j][i] - 2*P[j][i-1] + P[j][i-2];
    double Sp_M = P[j][i+1] - 2*P[j][i] + P[j][i-1];
    double Sp_R = P[j][i+2] - 2*P[j][i+1] + P[j][i];
    double Pe;
    if (std::abs(Sp_L) < std::abs(Sp_M) && std::abs(Sp_L) < std::abs(Sp_R))
        Pe = (-1.0/8.0)*P[j][i-2] + (3.0/4.0)*P[j][i-1] + (3.0/8.0)*P[j][i];
    else if (std::abs(Sp_R) < std::abs(Sp_L) && std::abs(Sp_R) < std::abs(Sp_M))
        Pe = (15.0/8.0)*P[j][i] - (5.0/4.0)*P[j][i+1] + (3.0/8.0)*P[j][i+2];
    else
        Pe = (3.0/8.0)*P[j][i-1] + (3.0/4.0)*P[j][i] - (1.0/8.0)*P[j][i+1];
    if (Pe < 1e-5) Pe = 1e-5;

    // rho
    double Srho_L = rho[j][i] - 2*rho[j][i-1] + rho[j][i-2];
    double Srho_M = rho[j][i+1] - 2*rho[j][i] + rho[j][i-1];
    double Srho_R = rho[j][i+2] - 2*rho[j][i+1] + rho[j][i];
    double rhoe;
    if (std::abs(Srho_L) < std::abs(Srho_M) && std::abs(Srho_L) < std::abs(Srho_R))
        rhoe = (-1.0/8.0)*rho[j][i-2] + (3.0/4.0)*rho[j][i-1] + (3.0/8.0)*rho[j][i];
    else if (std::abs(Srho_R) < std::abs(Srho_L) && std::abs(Srho_R) < std::abs(Srho_M))
        rhoe = (15.0/8.0)*rho[j][i] - (5.0/4.0)*rho[j][i+1] + (3.0/8.0)*rho[j][i+2];
    else
        rhoe = (3.0/8.0)*rho[j][i-1] + (3.0/4.0)*rho[j][i] - (1.0/8.0)*rho[j][i+1];
    if (rhoe < 1e-4) rhoe = 1e-4;

    return std::make_tuple(ue, ve, Pe, rhoe);
}
std::tuple<double, double, double, double> ENOreconstruction_y_top(
    int i, int j,
    const std::vector<std::vector<double>>& u,
    const std::vector<std::vector<double>>& v,
    const std::vector<std::vector<double>>& P,
    const std::vector<std::vector<double>>& rho)
{
    double Su_L = u[j][i] - 2*u[j-1][i] + u[j-2][i];
    double Su_M = u[j+1][i] - 2*u[j][i] + u[j-1][i];
    double Su_R = u[j+2][i] - 2*u[j+1][i] + u[j][i];
    double ue;
    if (std::abs(Su_L) < std::abs(Su_M) && std::abs(Su_L) < std::abs(Su_R))
        ue = (1.0/3.0)*u[j-2][i] - (7.0/6.0)*u[j-1][i] + (11.0/6.0)*u[j][i];
    else if (std::abs(Su_R) < std::abs(Su_L) && std::abs(Su_R) < std::abs(Su_M))
        ue = (1.0/3.0)*u[j][i] + (5.0/6.0)*u[j+1][i] - (1.0/6.0)*u[j+2][i];
    else
        ue = (-1.0/6.0)*u[j-1][i] + (5.0/6.0)*u[j][i] + (1.0/3.0)*u[j+1][i];

    double Sv_L = v[j][i] - 2*v[j-1][i] + v[j-2][i];
    double Sv_M = v[j+1][i] - 2*v[j][i] + v[j-1][i];
    double Sv_R = v[j+2][i] - 2*v[j+1][i] + v[j][i];
    double ve;
    if (std::abs(Sv_L) < std::abs(Sv_M) && std::abs(Sv_L) < std::abs(Sv_R))
        ve = (1.0/3.0)*v[j-2][i] - (7.0/6.0)*v[j-1][i] + (11.0/6.0)*v[j][i];
    else if (std::abs(Sv_R) < std::abs(Sv_L) && std::abs(Sv_R) < std::abs(Sv_M))
        ve = (1.0/3.0)*v[j][i] + (5.0/6.0)*v[j+1][i] - (1.0/6.0)*v[j+2][i];
    else
        ve = (-1.0/6.0)*v[j-1][i] + (5.0/6.0)*v[j][i] + (1.0/3.0)*v[j+1][i];

    // P
    double Sp_L = P[j][i] - 2*P[j-1][i] + P[j-2][i];
    double Sp_M = P[j+1][i] - 2*P[j][i] + P[j-1][i];
    double Sp_R = P[j+2][i] - 2*P[j+1][i] + P[j][i];
    double Pe;
    if (std::abs(Sp_L) < std::abs(Sp_M) && std::abs(Sp_L) < std::abs(Sp_R))
        Pe = (1.0/3.0)*P[j-2][i] - (7.0/6.0)*P[j-1][i] + (11.0/6.0)*P[j][i];
    else if (std::abs(Sp_R) < std::abs(Sp_L) && std::abs(Sp_R) < std::abs(Sp_M))
        Pe = (1.0/3.0)*P[j][i] + (5.0/6.0)*P[j+1][i] - (1.0/6.0)*P[j+2][i];
    else
        Pe = (-1.0/6.0)*P[j-1][i] + (5.0/6.0)*P[j][i] + (1.0/3.0)*P[j+1][i];
    if (Pe < 1e-5) Pe = 1e-5;

    // rho
    double Srho_L = rho[j][i] - 2*rho[j-1][i] + rho[j-2][i];
    double Srho_M = rho[j+1][i] - 2*rho[j][i] + rho[j-1][i];
    double Srho_R = rho[j+2][i] - 2*rho[j+1][i] + rho[j][i];
    double rhoe;
    if (std::abs(Srho_L) < std::abs(Srho_M) && std::abs(Srho_L) < std::abs(Srho_R))
        rhoe = (1.0/3.0)*rho[j-2][i] - (7.0/6.0)*rho[j-1][i] + (11.0/6.0)*rho[j][i];
    else if (std::abs(Srho_R) < std::abs(Srho_L) && std::abs(Srho_R) < std::abs(Srho_M))
        rhoe = (1.0/3.0)*rho[j][i] + (5.0/6.0)*rho[j+1][i] - (1.0/6.0)*rho[j+2][i];
    else
        rhoe = (-1.0/6.0)*rho[j-1][i] + (5.0/6.0)*rho[j][i] + (1.0/3.0)*rho[j+1][i];
    if (rhoe < 1e-4) rhoe = 1e-4;

    return std::make_tuple(ue, ve, Pe, rhoe);
}

std::tuple<double, double, double, double> ENOreconstruction_y_bottom(
    int i, int j,
    const std::vector<std::vector<double>>& u,
    const std::vector<std::vector<double>>& v,
    const std::vector<std::vector<double>>& P,
    const std::vector<std::vector<double>>& rho)
{
    // u
    double Su_L = u[j][i] - 2*u[j-1][i] + u[j-2][i];
    double Su_M = u[j+1][i] - 2*u[j][i] + u[j-1][i];
    double Su_R = u[j+2][i] - 2*u[j+1][i] + u[j][i];
    double ue;
    if (std::abs(Su_L) < std::abs(Su_M) && std::abs(Su_L) < std::abs(Su_R))
        ue = (-1.0/8.0)*u[j-2][i] + (3.0/4.0)*u[j-1][i] + (3.0/8.0)*u[j][i];
    else if (std::abs(Su_R) < std::abs(Su_L) && std::abs(Su_R) < std::abs(Su_M))
        ue = (15.0/8.0)*u[j][i] - (5.0/4.0)*u[j+1][i] + (3.0/8.0)*u[j+2][i];
    else
        ue = (3.0/8.0)*u[j-1][i] + (3.0/4.0)*u[j][i] - (1.0/8.0)*u[j+1][i];

    // v
    double Sv_L = v[j][i] - 2*v[j-1][i] + v[j-2][i];
    double Sv_M = v[j+1][i] - 2*v[j][i] + v[j-1][i];
    double Sv_R = v[j+2][i] - 2*v[j+1][i] + v[j][i];
    double ve;
    if (std::abs(Sv_L) < std::abs(Sv_M) && std::abs(Sv_L) < std::abs(Sv_R))
        ve = (-1.0/8.0)*v[j-2][i] + (3.0/4.0)*v[j-1][i] + (3.0/8.0)*v[j][i];
    else if (std::abs(Sv_R) < std::abs(Sv_L) && std::abs(Sv_R) < std::abs(Sv_M))
        ve = (15.0/8.0)*v[j][i] - (5.0/4.0)*v[j+1][i] + (3.0/8.0)*v[j+2][i];
    else
        ve = (3.0/8.0)*v[j-1][i] + (3.0/4.0)*v[j][i] - (1.0/8.0)*v[j+1][i];

    // P
    double Sp_L = P[j][i] - 2*P[j-1][i] + P[j-2][i];
    double Sp_M = P[j+1][i] - 2*P[j][i] + P[j-1][i];
    double Sp_R = P[j+2][i] - 2*P[j+1][i] + P[j][i];
    double Pe;
    if (std::abs(Sp_L) < std::abs(Sp_M) && std::abs(Sp_L) < std::abs(Sp_R))
        Pe = (-1.0/8.0)*P[j-2][i] + (3.0/4.0)*P[j-1][i] + (3.0/8.0)*P[j][i];
    else if (std::abs(Sp_R) < std::abs(Sp_L) && std::abs(Sp_R) < std::abs(Sp_M))
        Pe = (15.0/8.0)*P[j][i] - (5.0/4.0)*P[j+1][i] + (3.0/8.0)*P[j+2][i];
    else
        Pe = (3.0/8.0)*P[j-1][i] + (3.0/4.0)*P[j][i] - (1.0/8.0)*P[j+1][i];
    if (Pe < 1e-5) Pe = 1e-5;

    // rho
    double Srho_L = rho[j][i] - 2*rho[j-1][i] + rho[j-2][i];
    double Srho_M = rho[j+1][i] - 2*rho[j][i] + rho[j-1][i];
    double Srho_R = rho[j+2][i] - 2*rho[j+1][i] + rho[j][i];
    double rhoe;
    if (std::abs(Srho_L) < std::abs(Srho_M) && std::abs(Srho_L) < std::abs(Srho_R))
        rhoe = (-1.0/8.0)*rho[j-2][i] + (3.0/4.0)*rho[j-1][i] + (3.0/8.0)*rho[j][i];
    else if (std::abs(Srho_R) < std::abs(Srho_L) && std::abs(Srho_R) < std::abs(Srho_M))
        rhoe = (15.0/8.0)*rho[j][i] - (5.0/4.0)*rho[j+1][i] + (3.0/8.0)*rho[j+2][i];
    else
        rhoe = (3.0/8.0)*rho[j-1][i] + (3.0/4.0)*rho[j][i] - (1.0/8.0)*rho[j+1][i];
    if (rhoe < 1e-4) rhoe = 1e-4;

    return std::make_tuple(ue, ve, Pe, rhoe);
}
// CHECK: STATE_VECTOR_2D
std::tuple<
    std::vector<std::vector<double>>, // F1
    std::vector<std::vector<double>>, // F2
    std::vector<std::vector<double>>, // F3
    std::vector<std::vector<double>>, // F4
    std::vector<std::vector<double>>, // G1
    std::vector<std::vector<double>>, // G2
    std::vector<std::vector<double>>, // G3
    std::vector<std::vector<double>>  // G4
>computeFluxes2D_ENO(
    int cells_x, int cells_y, int ghost_cells,
    const std::vector<std::vector<double>>& rho,
    const std::vector<std::vector<double>>& u,
    const std::vector<std::vector<double>>& v,
    const std::vector<std::vector<double>>& P,
    double gamma,
    const std::string& nameRS)
{
    int Nx = cells_x + 2 * ghost_cells;
    int Ny = cells_y + 2 * ghost_cells;

    std::vector<std::vector<double>> F1(Ny, std::vector<double>(Nx - 1, 0.0));
    std::vector<std::vector<double>> F2(Ny, std::vector<double>(Nx - 1, 0.0));
    std::vector<std::vector<double>> F3(Ny, std::vector<double>(Nx - 1, 0.0));
    std::vector<std::vector<double>> F4(Ny, std::vector<double>(Nx - 1, 0.0));

    std::vector<std::vector<double>> G1(Ny - 1, std::vector<double>(Nx, 0.0));
    std::vector<std::vector<double>> G2(Ny - 1, std::vector<double>(Nx, 0.0));
    std::vector<std::vector<double>> G3(Ny - 1, std::vector<double>(Nx, 0.0));
    std::vector<std::vector<double>> G4(Ny - 1, std::vector<double>(Nx, 0.0));

    for (int j = 0; j < Ny; ++j) {
        for (int i = ghost_cells-1; i < Nx - ghost_cells; ++i) {
            auto [u_L, v_L, p_L, rho_L] = ENOreconstruction_x_right(i, j, u, v, P, rho);
            auto [u_R, v_R, p_R, rho_R] = ENOreconstruction_x_left(i+1, j, u, v, P, rho);

            auto [rho_face, u_face, p_face] = superSolve(rho_L, u_L, p_L,
                                                         rho_R, u_R, p_R,
                                                         nameRS);
            double v_star = (u_face >= 0.0) ? v_L : v_R;
            double E_face = p_face / (gamma - 1.0) + 0.5 * rho_face * (u_face*u_face + v_star*v_star);

            F1[j][i] = rho_face * u_face;
            F2[j][i] = rho_face * u_face * u_face + p_face;
            F3[j][i] = rho_face * u_face * v_star;
            F4[j][i] = u_face * (E_face + p_face);
        }
    }

    for (int j = ghost_cells-1; j < Ny - ghost_cells; ++j) {
        for (int i = 0; i < Nx; ++i) {
            auto [u_B, v_B, p_B, rho_B] = ENOreconstruction_y_top(i, j, u, v, P, rho);
            auto [u_T, v_T, p_T, rho_T] = ENOreconstruction_y_bottom(i, j+1, u, v, P, rho);

            auto [rho_face, v_face, p_face] = superSolve(rho_B, v_B, p_B,
                                                         rho_T, v_T, p_T,
                                                         nameRS);
            double u_star = (v_face >= 0.0) ? u_B : u_T;
            double E_face = p_face / (gamma - 1.0) + 0.5 * rho_face * (u_star*u_star + v_face*v_face);

            G1[j][i] = rho_face * v_face;
            G2[j][i] = rho_face * v_face * u_star;
            G3[j][i] = rho_face * v_face * v_face + p_face;
            G4[j][i] = v_face * (E_face + p_face);
        }
    }

    return {F1, F2, F3, F4, G1, G2, G3, G4};
}
std::tuple<
    std::vector<std::vector<double>>, // R1 (ernnr)
    std::vector<std::vector<double>>, // R2 (x-cedoeun)
    std::vector<std::vector<double>>, // R3 (y-cedoeun)
    std::vector<std::vector<double>>  // R4 (yildac?)
>
computeRHS_ENO_2D(
    int cells_x, int cells_y, int ghost_cells,
    const std::vector<std::vector<double>>& rho,
    const std::vector<std::vector<double>>& u,
    const std::vector<std::vector<double>>& v,
    const std::vector<std::vector<double>>& P,
    double hx, double hy,
    double gamma,
    const std::string& nameRS)
{
    int Nx = cells_x + 2 * ghost_cells;
    int Ny = cells_y + 2 * ghost_cells;

    auto [F1, F2, F3, F4, G1, G2, G3, G4] = computeFluxes2D_ENO(
        cells_x, cells_y, ghost_cells,
        rho, u, v, P, gamma, nameRS);

    std::vector<std::vector<double>> R1(Ny, std::vector<double>(Nx, 0.0));
    std::vector<std::vector<double>> R2(Ny, std::vector<double>(Nx, 0.0));
    std::vector<std::vector<double>> R3(Ny, std::vector<double>(Nx, 0.0));
    std::vector<std::vector<double>> R4(Ny, std::vector<double>(Nx, 0.0));

    for (int j = ghost_cells; j < Ny - ghost_cells; ++j) {
        for (int i = ghost_cells; i < Nx - ghost_cells; ++i) {
            double dFx1 = F1[j][i] - F1[j][i-1];
            double dFx2 = F2[j][i] - F2[j][i-1];
            double dFx3 = F3[j][i] - F3[j][i-1];
            double dFx4 = F4[j][i] - F4[j][i-1];

            double dGy1 = G1[j][i] - G1[j-1][i];
            double dGy2 = G2[j][i] - G2[j-1][i];
            double dGy3 = G3[j][i] - G3[j-1][i];
            double dGy4 = G4[j][i] - G4[j-1][i];

            R1[j][i] = -(dFx1 / hx + dGy1 / hy);
            R2[j][i] = -(dFx2 / hx + dGy2 / hy);
            R3[j][i] = -(dFx3 / hx + dGy3 / hy);
            R4[j][i] = -(dFx4 / hx + dGy4 / hy);
        }
    }

    return {R1, R2, R3, R4};
}

void updatePrimitiveVariables2D(
    int cells_x, int cells_y, int ghost_cells,
    std::vector<std::vector<double>>& u,
    std::vector<std::vector<double>>& v,
    std::vector<std::vector<double>>& P,
    std::vector<std::vector<double>>& rho,
    std::vector<std::vector<double>>& I,
    const std::vector<std::vector<double>>& U1,
    const std::vector<std::vector<double>>& U2,
    const std::vector<std::vector<double>>& U3,
    const std::vector<std::vector<double>>& U4,
    double gamma)
{
    const double eps_rho = 1e-5;
    const double eps_e   = 1e-6;

    int Nx = cells_x + 2 * ghost_cells;
    int Ny = cells_y + 2 * ghost_cells;

    for (int j = ghost_cells; j < Ny - ghost_cells; ++j) {
        for (int i = ghost_cells; i < Nx - ghost_cells; ++i) {
            rho[j][i] = std::max(U1[j][i], eps_rho);
            u[j][i]   = U2[j][i] / rho[j][i];
            v[j][i]   = U3[j][i] / rho[j][i];
            double Ek = 0.5 * rho[j][i] * (u[j][i]*u[j][i] + v[j][i]*v[j][i]);
            double e_in = (U4[j][i] - Ek) / rho[j][i];
            if (e_in < eps_e) e_in = eps_e;
            P[j][i] = e_in * rho[j][i] * (gamma - 1.0);
            I[j][i] = U4[j][i];
        }
    }
}
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
    const std::string& nameRS)
{
    const double gamma = 1.4;
    int Nx = cells_x + 2 * ghost_cells;
    int Ny = cells_y + 2 * ghost_cells;

    std::vector<std::vector<double>> U1(Ny, std::vector<double>(Nx));
    std::vector<std::vector<double>> U2(Ny, std::vector<double>(Nx));
    std::vector<std::vector<double>> U3(Ny, std::vector<double>(Nx));
    std::vector<std::vector<double>> U4(Ny, std::vector<double>(Nx));

    for (int j = 0; j < Ny; ++j) {
        for (int i = 0; i < Nx; ++i) {
            U1[j][i] = rho[j][i];
            U2[j][i] = rho[j][i] * u[j][i];
            U3[j][i] = rho[j][i] * v[j][i];
            U4[j][i] = P[j][i] / (gamma - 1.0) + 0.5 * rho[j][i] * (u[j][i]*u[j][i] + v[j][i]*v[j][i]);
        }
    }

    applyBCs2D(rho, u, v, P, Nx, Ny, cells_x, cells_y, ghost_cells,
               left_bc, right_bc, bottom_bc, top_bc);

    std::vector<std::vector<double>> U1_new(Ny, std::vector<double>(Nx, 0.0));
    std::vector<std::vector<double>> U2_new(Ny, std::vector<double>(Nx, 0.0));
    std::vector<std::vector<double>> U3_new(Ny, std::vector<double>(Nx, 0.0));
    std::vector<std::vector<double>> U4_new(Ny, std::vector<double>(Nx, 0.0));

    if (time_integrator == "euler") {
        auto [R1, R2, R3, R4] = computeRHS_ENO_2D(
            cells_x, cells_y, ghost_cells,
            rho, u, v, P, hx, hy, gamma, nameRS);

        for (int j = ghost_cells; j < Ny - ghost_cells; ++j) {
            for (int i = ghost_cells; i < Nx - ghost_cells; ++i) {
                U1_new[j][i] = U1[j][i] + tau * R1[j][i];
                U2_new[j][i] = U2[j][i] + tau * R2[j][i];
                U3_new[j][i] = U3[j][i] + tau * R3[j][i];
                U4_new[j][i] = U4[j][i] + tau * R4[j][i];
            }
        }

        updatePrimitiveVariables2D(cells_x, cells_y, ghost_cells,
                                   u, v, P, rho, I,
                                   U1_new, U2_new, U3_new, U4_new, gamma);
    }
    else if (time_integrator == "rk2") {
        auto [R1_1, R2_1, R3_1, R4_1] = computeRHS_ENO_2D(
            cells_x, cells_y, ghost_cells,
            rho, u, v, P, hx, hy, gamma, nameRS);

        std::vector<std::vector<double>> U1_temp(Ny, std::vector<double>(Nx, 0.0));
        std::vector<std::vector<double>> U2_temp(Ny, std::vector<double>(Nx, 0.0));
        std::vector<std::vector<double>> U3_temp(Ny, std::vector<double>(Nx, 0.0));
        std::vector<std::vector<double>> U4_temp(Ny, std::vector<double>(Nx, 0.0));

        for (int j = ghost_cells; j < Ny - ghost_cells; ++j) {
            for (int i = ghost_cells; i < Nx - ghost_cells; ++i) {
                U1_temp[j][i] = U1[j][i] + tau * R1_1[j][i];
                U2_temp[j][i] = U2[j][i] + tau * R2_1[j][i];
                U3_temp[j][i] = U3[j][i] + tau * R3_1[j][i];
                U4_temp[j][i] = U4[j][i] + tau * R4_1[j][i];
            }
        }

        std::vector<std::vector<double>> u_temp = u;
        std::vector<std::vector<double>> v_temp = v;
        std::vector<std::vector<double>> P_temp = P;
        std::vector<std::vector<double>> rho_temp = rho;
        std::vector<std::vector<double>> I_temp = I; // il cndieucolnn?, ii ae? niaelnnceinnc

        updatePrimitiveVariables2D(cells_x, cells_y, ghost_cells,
                                   u_temp, v_temp, P_temp, rho_temp, I_temp,
                                   U1_temp, U2_temp, U3_temp, U4_temp, gamma);

        applyBCs2D(rho_temp, u_temp, v_temp, P_temp, Nx, Ny,
                   cells_x, cells_y, ghost_cells,
                   left_bc, right_bc, bottom_bc, top_bc);

        auto [R1_2, R2_2, R3_2, R4_2] = computeRHS_ENO_2D(
            cells_x, cells_y, ghost_cells,
            rho_temp, u_temp, v_temp, P_temp, hx, hy, gamma, nameRS);

        for (int j = ghost_cells; j < Ny - ghost_cells; ++j) {
            for (int i = ghost_cells; i < Nx - ghost_cells; ++i) {
                U1_new[j][i] = U1[j][i] + 0.5 * tau * (R1_1[j][i] + R1_2[j][i]);
                U2_new[j][i] = U2[j][i] + 0.5 * tau * (R2_1[j][i] + R2_2[j][i]);
                U3_new[j][i] = U3[j][i] + 0.5 * tau * (R3_1[j][i] + R3_2[j][i]);
                U4_new[j][i] = U4[j][i] + 0.5 * tau * (R4_1[j][i] + R4_2[j][i]);
            }
        }

        updatePrimitiveVariables2D(cells_x, cells_y, ghost_cells,
                                   u, v, P, rho, I,
                                   U1_new, U2_new, U3_new, U4_new, gamma);
    }
    else {
        std::cerr << "Time integrator " << time_integrator << " not implemented in ENO2DSolve.\n";
        return;
    }

    applyBCs2D(rho, u, v, P, Nx, Ny, cells_x, cells_y, ghost_cells,
               left_bc, right_bc, bottom_bc, top_bc);

    t_total += tau;
    std::cout << "t = " << t_total << ", tau = " << tau
              << ", method = " << time_integrator
              << ", RS = " << nameRS << std::endl;
}


std::tuple<double, double, double, double> WENO5reconstruction_x_right(
    int i, int j,
    const std::vector<std::vector<double>>& u,
    const std::vector<std::vector<double>>& v,
    const std::vector<std::vector<double>>& P,
    const std::vector<std::vector<double>>& rho)
{
    double epsilon = 1e-6;
    double c0 = 0.3, c1 = 0.6, c2 = 0.1;

    double u0 = (2.0/6.0)*u[j][i-2] - (7.0/6.0)*u[j][i-1] + (11.0/6.0)*u[j][i];
    double u1 = (-1.0/6.0)*u[j][i-1] + (5.0/6.0)*u[j][i] + (2.0/6.0)*u[j][i+1];
    double u2 = (2.0/6.0)*u[j][i] + (5.0/6.0)*u[j][i+1] - (1.0/6.0)*u[j][i+2];

    double s0_u = (13.0/12.0)*pow(u[j][i-2] - 2*u[j][i-1] + u[j][i], 2) +
                  0.25*pow(u[j][i-2] - 4*u[j][i-1] + 3*u[j][i], 2);
    double s1_u = (13.0/12.0)*pow(u[j][i-1] - 2*u[j][i] + u[j][i+1], 2) +
                  0.25*pow(u[j][i-1] - u[j][i+1], 2);
    double s2_u = (13.0/12.0)*pow(u[j][i] - 2*u[j][i+1] + u[j][i+2], 2) +
                  0.25*pow(3*u[j][i] - 4*u[j][i+1] + u[j][i+2], 2);

    double alpha0_u = c0 / pow(epsilon + s0_u, 2);
    double alpha1_u = c1 / pow(epsilon + s1_u, 2);
    double alpha2_u = c2 / pow(epsilon + s2_u, 2);
    double sum_alpha_u = alpha0_u + alpha1_u + alpha2_u;
    double w0_u = alpha0_u / sum_alpha_u;
    double w1_u = alpha1_u / sum_alpha_u;
    double w2_u = alpha2_u / sum_alpha_u;

    double ue = w0_u * u0 + w1_u * u1 + w2_u * u2;

    double v0 = (2.0/6.0)*v[j][i-2] - (7.0/6.0)*v[j][i-1] + (11.0/6.0)*v[j][i];
    double v1 = (-1.0/6.0)*v[j][i-1] + (5.0/6.0)*v[j][i] + (2.0/6.0)*v[j][i+1];
    double v2 = (2.0/6.0)*v[j][i] + (5.0/6.0)*v[j][i+1] - (1.0/6.0)*v[j][i+2];

    double s0_v = (13.0/12.0)*pow(v[j][i-2] - 2*v[j][i-1] + v[j][i], 2) +
                  0.25*pow(v[j][i-2] - 4*v[j][i-1] + 3*v[j][i], 2);
    double s1_v = (13.0/12.0)*pow(v[j][i-1] - 2*v[j][i] + v[j][i+1], 2) +
                  0.25*pow(v[j][i-1] - v[j][i+1], 2);
    double s2_v = (13.0/12.0)*pow(v[j][i] - 2*v[j][i+1] + v[j][i+2], 2) +
                  0.25*pow(3*v[j][i] - 4*v[j][i+1] + v[j][i+2], 2);

    double alpha0_v = c0 / pow(epsilon + s0_v, 2);
    double alpha1_v = c1 / pow(epsilon + s1_v, 2);
    double alpha2_v = c2 / pow(epsilon + s2_v, 2);
    double sum_alpha_v = alpha0_v + alpha1_v + alpha2_v;
    double w0_v = alpha0_v / sum_alpha_v;
    double w1_v = alpha1_v / sum_alpha_v;
    double w2_v = alpha2_v / sum_alpha_v;

    double ve = w0_v * v0 + w1_v * v1 + w2_v * v2;

    double P0 = (2.0/6.0)*P[j][i-2] - (7.0/6.0)*P[j][i-1] + (11.0/6.0)*P[j][i];
    double P1 = (-1.0/6.0)*P[j][i-1] + (5.0/6.0)*P[j][i] + (2.0/6.0)*P[j][i+1];
    double P2 = (2.0/6.0)*P[j][i] + (5.0/6.0)*P[j][i+1] - (1.0/6.0)*P[j][i+2];

    double s0_P = (13.0/12.0)*pow(P[j][i-2] - 2*P[j][i-1] + P[j][i], 2) +
                  0.25*pow(P[j][i-2] - 4*P[j][i-1] + 3*P[j][i], 2);
    double s1_P = (13.0/12.0)*pow(P[j][i-1] - 2*P[j][i] + P[j][i+1], 2) +
                  0.25*pow(P[j][i-1] - P[j][i+1], 2);
    double s2_P = (13.0/12.0)*pow(P[j][i] - 2*P[j][i+1] + P[j][i+2], 2) +
                  0.25*pow(3*P[j][i] - 4*P[j][i+1] + P[j][i+2], 2);

    double alpha0_P = c0 / pow(epsilon + s0_P, 2);
    double alpha1_P = c1 / pow(epsilon + s1_P, 2);
    double alpha2_P = c2 / pow(epsilon + s2_P, 2);
    double sum_alpha_P = alpha0_P + alpha1_P + alpha2_P;
    double w0_P = alpha0_P / sum_alpha_P;
    double w1_P = alpha1_P / sum_alpha_P;
    double w2_P = alpha2_P / sum_alpha_P;

    double Pe = w0_P * P0 + w1_P * P1 + w2_P * P2;
    if (Pe < 1e-5) Pe = 1e-5;

    double rho0 = (2.0/6.0)*rho[j][i-2] - (7.0/6.0)*rho[j][i-1] + (11.0/6.0)*rho[j][i];
    double rho1 = (-1.0/6.0)*rho[j][i-1] + (5.0/6.0)*rho[j][i] + (2.0/6.0)*rho[j][i+1];
    double rho2 = (2.0/6.0)*rho[j][i] + (5.0/6.0)*rho[j][i+1] - (1.0/6.0)*rho[j][i+2];

    double s0_rho = (13.0/12.0)*pow(rho[j][i-2] - 2*rho[j][i-1] + rho[j][i], 2) +
                    0.25*pow(rho[j][i-2] - 4*rho[j][i-1] + 3*rho[j][i], 2);
    double s1_rho = (13.0/12.0)*pow(rho[j][i-1] - 2*rho[j][i] + rho[j][i+1], 2) +
                    0.25*pow(rho[j][i-1] - rho[j][i+1], 2);
    double s2_rho = (13.0/12.0)*pow(rho[j][i] - 2*rho[j][i+1] + rho[j][i+2], 2) +
                    0.25*pow(3*rho[j][i] - 4*rho[j][i+1] + rho[j][i+2], 2);

    double alpha0_rho = c0 / pow(epsilon + s0_rho, 2);
    double alpha1_rho = c1 / pow(epsilon + s1_rho, 2);
    double alpha2_rho = c2 / pow(epsilon + s2_rho, 2);
    double sum_alpha_rho = alpha0_rho + alpha1_rho + alpha2_rho;
    double w0_rho = alpha0_rho / sum_alpha_rho;
    double w1_rho = alpha1_rho / sum_alpha_rho;
    double w2_rho = alpha2_rho / sum_alpha_rho;

    double rhoe = w0_rho * rho0 + w1_rho * rho1 + w2_rho * rho2;
    if (rhoe < 1e-5) rhoe = 1e-5;

    return std::make_tuple(ue, ve, Pe, rhoe);
}

std::tuple<double, double, double, double> WENO5reconstruction_x_left(
    int i, int j,
    const std::vector<std::vector<double>>& u,
    const std::vector<std::vector<double>>& v,
    const std::vector<std::vector<double>>& P,
    const std::vector<std::vector<double>>& rho)
{
    double epsilon = 1e-6;
    double c0 = 0.1, c1 = 0.6, c2 = 0.3;

    double u0 = (-1.0/6.0)*u[j][i-2] + (5.0/6.0)*u[j][i-1] + (1.0/3.0)*u[j][i];
    double u1 = (1.0/3.0)*u[j][i-1] + (5.0/6.0)*u[j][i] - (1.0/6.0)*u[j][i+1];
    double u2 = (11.0/6.0)*u[j][i] - (7.0/6.0)*u[j][i+1] + (1.0/3.0)*u[j][i+2];

    double s0_u = (13.0/12.0)*pow(u[j][i-2] - 2*u[j][i-1] + u[j][i], 2) +
                  0.25*pow(u[j][i-2] - 4*u[j][i-1] + 3*u[j][i], 2);
    double s1_u = (13.0/12.0)*pow(u[j][i-1] - 2*u[j][i] + u[j][i+1], 2) +
                  0.25*pow(u[j][i-1] - u[j][i+1], 2);
    double s2_u = (13.0/12.0)*pow(u[j][i] - 2*u[j][i+1] + u[j][i+2], 2) +
                  0.25*pow(3*u[j][i] - 4*u[j][i+1] + u[j][i+2], 2);

    double alpha0_u = c0 / pow(epsilon + s0_u, 2);
    double alpha1_u = c1 / pow(epsilon + s1_u, 2);
    double alpha2_u = c2 / pow(epsilon + s2_u, 2);
    double sum_alpha_u = alpha0_u + alpha1_u + alpha2_u;
    double w0_u = alpha0_u / sum_alpha_u;
    double w1_u = alpha1_u / sum_alpha_u;
    double w2_u = alpha2_u / sum_alpha_u;

    double ue = w0_u * u0 + w1_u * u1 + w2_u * u2;

    double v0 = (-1.0/6.0)*v[j][i-2] + (5.0/6.0)*v[j][i-1] + (1.0/3.0)*v[j][i];
    double v1 = (1.0/3.0)*v[j][i-1] + (5.0/6.0)*v[j][i] - (1.0/6.0)*v[j][i+1];
    double v2 = (11.0/6.0)*v[j][i] - (7.0/6.0)*v[j][i+1] + (1.0/3.0)*v[j][i+2];

    double s0_v = (13.0/12.0)*pow(v[j][i-2] - 2*v[j][i-1] + v[j][i], 2) +
                  0.25*pow(v[j][i-2] - 4*v[j][i-1] + 3*v[j][i], 2);
    double s1_v = (13.0/12.0)*pow(v[j][i-1] - 2*v[j][i] + v[j][i+1], 2) +
                  0.25*pow(v[j][i-1] - v[j][i+1], 2);
    double s2_v = (13.0/12.0)*pow(v[j][i] - 2*v[j][i+1] + v[j][i+2], 2) +
                  0.25*pow(3*v[j][i] - 4*v[j][i+1] + v[j][i+2], 2);

    double alpha0_v = c0 / pow(epsilon + s0_v, 2);
    double alpha1_v = c1 / pow(epsilon + s1_v, 2);
    double alpha2_v = c2 / pow(epsilon + s2_v, 2);
    double sum_alpha_v = alpha0_v + alpha1_v + alpha2_v;
    double w0_v = alpha0_v / sum_alpha_v;
    double w1_v = alpha1_v / sum_alpha_v;
    double w2_v = alpha2_v / sum_alpha_v;

    double ve = w0_v * v0 + w1_v * v1 + w2_v * v2;

    double P0 = (-1.0/6.0)*P[j][i-2] + (5.0/6.0)*P[j][i-1] + (1.0/3.0)*P[j][i];
    double P1 = (1.0/3.0)*P[j][i-1] + (5.0/6.0)*P[j][i] - (1.0/6.0)*P[j][i+1];
    double P2 = (11.0/6.0)*P[j][i] - (7.0/6.0)*P[j][i+1] + (1.0/3.0)*P[j][i+2];

    double s0_P = (13.0/12.0)*pow(P[j][i-2] - 2*P[j][i-1] + P[j][i], 2) +
                  0.25*pow(P[j][i-2] - 4*P[j][i-1] + 3*P[j][i], 2);
    double s1_P = (13.0/12.0)*pow(P[j][i-1] - 2*P[j][i] + P[j][i+1], 2) +
                  0.25*pow(P[j][i-1] - P[j][i+1], 2);
    double s2_P = (13.0/12.0)*pow(P[j][i] - 2*P[j][i+1] + P[j][i+2], 2) +
                  0.25*pow(3*P[j][i] - 4*P[j][i+1] + P[j][i+2], 2);

    double alpha0_P = c0 / pow(epsilon + s0_P, 2);
    double alpha1_P = c1 / pow(epsilon + s1_P, 2);
    double alpha2_P = c2 / pow(epsilon + s2_P, 2);
    double sum_alpha_P = alpha0_P + alpha1_P + alpha2_P;
    double w0_P = alpha0_P / sum_alpha_P;
    double w1_P = alpha1_P / sum_alpha_P;
    double w2_P = alpha2_P / sum_alpha_P;

    double Pe = w0_P * P0 + w1_P * P1 + w2_P * P2;
    if (Pe < 1e-5) Pe = 1e-5;

    double rho0 = (-1.0/6.0)*rho[j][i-2] + (5.0/6.0)*rho[j][i-1] + (1.0/3.0)*rho[j][i];
    double rho1 = (1.0/3.0)*rho[j][i-1] + (5.0/6.0)*rho[j][i] - (1.0/6.0)*rho[j][i+1];
    double rho2 = (11.0/6.0)*rho[j][i] - (7.0/6.0)*rho[j][i+1] + (1.0/3.0)*rho[j][i+2];

    double s0_rho = (13.0/12.0)*pow(rho[j][i-2] - 2*rho[j][i-1] + rho[j][i], 2) +
                    0.25*pow(rho[j][i-2] - 4*rho[j][i-1] + 3*rho[j][i], 2);
    double s1_rho = (13.0/12.0)*pow(rho[j][i-1] - 2*rho[j][i] + rho[j][i+1], 2) +
                    0.25*pow(rho[j][i-1] - rho[j][i+1], 2);
    double s2_rho = (13.0/12.0)*pow(rho[j][i] - 2*rho[j][i+1] + rho[j][i+2], 2) +
                    0.25*pow(3*rho[j][i] - 4*rho[j][i+1] + rho[j][i+2], 2);

    double alpha0_rho = c0 / pow(epsilon + s0_rho, 2);
    double alpha1_rho = c1 / pow(epsilon + s1_rho, 2);
    double alpha2_rho = c2 / pow(epsilon + s2_rho, 2);
    double sum_alpha_rho = alpha0_rho + alpha1_rho + alpha2_rho;
    double w0_rho = alpha0_rho / sum_alpha_rho;
    double w1_rho = alpha1_rho / sum_alpha_rho;
    double w2_rho = alpha2_rho / sum_alpha_rho;

    double rhoe = w0_rho * rho0 + w1_rho * rho1 + w2_rho * rho2;
    if (rhoe < 1e-5) rhoe = 1e-5;

    return std::make_tuple(ue, ve, Pe, rhoe);
}

std::tuple<double, double, double, double> WENO5reconstruction_y_top(
    int i, int j,
    const std::vector<std::vector<double>>& u,
    const std::vector<std::vector<double>>& v,
    const std::vector<std::vector<double>>& P,
    const std::vector<std::vector<double>>& rho)
{
    double epsilon = 1e-6;
    double c0 = 0.3, c1 = 0.6, c2 = 0.1;  // alnr ae? aldoile adric (ere ae? ddraie)

    double u0 = (2.0/6.0)*u[j-2][i] - (7.0/6.0)*u[j-1][i] + (11.0/6.0)*u[j][i];
    double u1 = (-1.0/6.0)*u[j-1][i] + (5.0/6.0)*u[j][i] + (2.0/6.0)*u[j+1][i];
    double u2 = (2.0/6.0)*u[j][i] + (5.0/6.0)*u[j+1][i] - (1.0/6.0)*u[j+2][i];

    double s0_u = (13.0/12.0)*pow(u[j-2][i] - 2*u[j-1][i] + u[j][i], 2) +
                  0.25*pow(u[j-2][i] - 4*u[j-1][i] + 3*u[j][i], 2);
    double s1_u = (13.0/12.0)*pow(u[j-1][i] - 2*u[j][i] + u[j+1][i], 2) +
                  0.25*pow(u[j-1][i] - u[j+1][i], 2);
    double s2_u = (13.0/12.0)*pow(u[j][i] - 2*u[j+1][i] + u[j+2][i], 2) +
                  0.25*pow(3*u[j][i] - 4*u[j+1][i] + u[j+2][i], 2);

    double alpha0_u = c0 / pow(epsilon + s0_u, 2);
    double alpha1_u = c1 / pow(epsilon + s1_u, 2);
    double alpha2_u = c2 / pow(epsilon + s2_u, 2);
    double sum_alpha_u = alpha0_u + alpha1_u + alpha2_u;
    double w0_u = alpha0_u / sum_alpha_u;
    double w1_u = alpha1_u / sum_alpha_u;
    double w2_u = alpha2_u / sum_alpha_u;

    double ue = w0_u * u0 + w1_u * u1 + w2_u * u2;

    double v0 = (2.0/6.0)*v[j-2][i] - (7.0/6.0)*v[j-1][i] + (11.0/6.0)*v[j][i];
    double v1 = (-1.0/6.0)*v[j-1][i] + (5.0/6.0)*v[j][i] + (2.0/6.0)*v[j+1][i];
    double v2 = (2.0/6.0)*v[j][i] + (5.0/6.0)*v[j+1][i] - (1.0/6.0)*v[j+2][i];

    double s0_v = (13.0/12.0)*pow(v[j-2][i] - 2*v[j-1][i] + v[j][i], 2) +
                  0.25*pow(v[j-2][i] - 4*v[j-1][i] + 3*v[j][i], 2);
    double s1_v = (13.0/12.0)*pow(v[j-1][i] - 2*v[j][i] + v[j+1][i], 2) +
                  0.25*pow(v[j-1][i] - v[j+1][i], 2);
    double s2_v = (13.0/12.0)*pow(v[j][i] - 2*v[j+1][i] + v[j+2][i], 2) +
                  0.25*pow(3*v[j][i] - 4*v[j+1][i] + v[j+2][i], 2);

    double alpha0_v = c0 / pow(epsilon + s0_v, 2);
    double alpha1_v = c1 / pow(epsilon + s1_v, 2);
    double alpha2_v = c2 / pow(epsilon + s2_v, 2);
    double sum_alpha_v = alpha0_v + alpha1_v + alpha2_v;
    double w0_v = alpha0_v / sum_alpha_v;
    double w1_v = alpha1_v / sum_alpha_v;
    double w2_v = alpha2_v / sum_alpha_v;

    double ve = w0_v * v0 + w1_v * v1 + w2_v * v2;

    double P0 = (2.0/6.0)*P[j-2][i] - (7.0/6.0)*P[j-1][i] + (11.0/6.0)*P[j][i];
    double P1 = (-1.0/6.0)*P[j-1][i] + (5.0/6.0)*P[j][i] + (2.0/6.0)*P[j+1][i];
    double P2 = (2.0/6.0)*P[j][i] + (5.0/6.0)*P[j+1][i] - (1.0/6.0)*P[j+2][i];

    double s0_P = (13.0/12.0)*pow(P[j-2][i] - 2*P[j-1][i] + P[j][i], 2) +
                  0.25*pow(P[j-2][i] - 4*P[j-1][i] + 3*P[j][i], 2);
    double s1_P = (13.0/12.0)*pow(P[j-1][i] - 2*P[j][i] + P[j+1][i], 2) +
                  0.25*pow(P[j-1][i] - P[j+1][i], 2);
    double s2_P = (13.0/12.0)*pow(P[j][i] - 2*P[j+1][i] + P[j+2][i], 2) +
                  0.25*pow(3*P[j][i] - 4*P[j+1][i] + P[j+2][i], 2);

    double alpha0_P = c0 / pow(epsilon + s0_P, 2);
    double alpha1_P = c1 / pow(epsilon + s1_P, 2);
    double alpha2_P = c2 / pow(epsilon + s2_P, 2);
    double sum_alpha_P = alpha0_P + alpha1_P + alpha2_P;
    double w0_P = alpha0_P / sum_alpha_P;
    double w1_P = alpha1_P / sum_alpha_P;
    double w2_P = alpha2_P / sum_alpha_P;

    double Pe = w0_P * P0 + w1_P * P1 + w2_P * P2;
    if (Pe < 1e-5) Pe = 1e-5;

    double rho0 = (2.0/6.0)*rho[j-2][i] - (7.0/6.0)*rho[j-1][i] + (11.0/6.0)*rho[j][i];
    double rho1 = (-1.0/6.0)*rho[j-1][i] + (5.0/6.0)*rho[j][i] + (2.0/6.0)*rho[j+1][i];
    double rho2 = (2.0/6.0)*rho[j][i] + (5.0/6.0)*rho[j+1][i] - (1.0/6.0)*rho[j+2][i];

    double s0_rho = (13.0/12.0)*pow(rho[j-2][i] - 2*rho[j-1][i] + rho[j][i], 2) +
                    0.25*pow(rho[j-2][i] - 4*rho[j-1][i] + 3*rho[j][i], 2);
    double s1_rho = (13.0/12.0)*pow(rho[j-1][i] - 2*rho[j][i] + rho[j+1][i], 2) +
                    0.25*pow(rho[j-1][i] - rho[j+1][i], 2);
    double s2_rho = (13.0/12.0)*pow(rho[j][i] - 2*rho[j+1][i] + rho[j+2][i], 2) +
                    0.25*pow(3*rho[j][i] - 4*rho[j+1][i] + rho[j+2][i], 2);

    double alpha0_rho = c0 / pow(epsilon + s0_rho, 2);
    double alpha1_rho = c1 / pow(epsilon + s1_rho, 2);
    double alpha2_rho = c2 / pow(epsilon + s2_rho, 2);
    double sum_alpha_rho = alpha0_rho + alpha1_rho + alpha2_rho;
    double w0_rho = alpha0_rho / sum_alpha_rho;
    double w1_rho = alpha1_rho / sum_alpha_rho;
    double w2_rho = alpha2_rho / sum_alpha_rho;

    double rhoe = w0_rho * rho0 + w1_rho * rho1 + w2_rho * rho2;
    if (rhoe < 1e-5) rhoe = 1e-5;

    return std::make_tuple(ue, ve, Pe, rhoe);
}

std::tuple<double, double, double, double> WENO5reconstruction_y_bottom(
    int i, int j,
    const std::vector<std::vector<double>>& u,
    const std::vector<std::vector<double>>& v,
    const std::vector<std::vector<double>>& P,
    const std::vector<std::vector<double>>& rho)
{
    double epsilon = 1e-6;
    double c0 = 0.1, c1 = 0.6, c2 = 0.3;  // alnr ae? iccile adric (ere ae? elaie)

    double u0 = (-1.0/6.0)*u[j-2][i] + (5.0/6.0)*u[j-1][i] + (1.0/3.0)*u[j][i];
    double u1 = (1.0/3.0)*u[j-1][i] + (5.0/6.0)*u[j][i] - (1.0/6.0)*u[j+1][i];
    double u2 = (11.0/6.0)*u[j][i] - (7.0/6.0)*u[j+1][i] + (1.0/3.0)*u[j+2][i];

    double s0_u = (13.0/12.0)*pow(u[j-2][i] - 2*u[j-1][i] + u[j][i], 2) +
                  0.25*pow(u[j-2][i] - 4*u[j-1][i] + 3*u[j][i], 2);
    double s1_u = (13.0/12.0)*pow(u[j-1][i] - 2*u[j][i] + u[j+1][i], 2) +
                  0.25*pow(u[j-1][i] - u[j+1][i], 2);
    double s2_u = (13.0/12.0)*pow(u[j][i] - 2*u[j+1][i] + u[j+2][i], 2) +
                  0.25*pow(3*u[j][i] - 4*u[j+1][i] + u[j+2][i], 2);

    double alpha0_u = c0 / pow(epsilon + s0_u, 2);
    double alpha1_u = c1 / pow(epsilon + s1_u, 2);
    double alpha2_u = c2 / pow(epsilon + s2_u, 2);
    double sum_alpha_u = alpha0_u + alpha1_u + alpha2_u;
    double w0_u = alpha0_u / sum_alpha_u;
    double w1_u = alpha1_u / sum_alpha_u;
    double w2_u = alpha2_u / sum_alpha_u;

    double ue = w0_u * u0 + w1_u * u1 + w2_u * u2;

    double v0 = (-1.0/6.0)*v[j-2][i] + (5.0/6.0)*v[j-1][i] + (1.0/3.0)*v[j][i];
    double v1 = (1.0/3.0)*v[j-1][i] + (5.0/6.0)*v[j][i] - (1.0/6.0)*v[j+1][i];
    double v2 = (11.0/6.0)*v[j][i] - (7.0/6.0)*v[j+1][i] + (1.0/3.0)*v[j+2][i];

    double s0_v = (13.0/12.0)*pow(v[j-2][i] - 2*v[j-1][i] + v[j][i], 2) +
                  0.25*pow(v[j-2][i] - 4*v[j-1][i] + 3*v[j][i], 2);
    double s1_v = (13.0/12.0)*pow(v[j-1][i] - 2*v[j][i] + v[j+1][i], 2) +
                  0.25*pow(v[j-1][i] - v[j+1][i], 2);
    double s2_v = (13.0/12.0)*pow(v[j][i] - 2*v[j+1][i] + v[j+2][i], 2) +
                  0.25*pow(3*v[j][i] - 4*v[j+1][i] + v[j+2][i], 2);

    double alpha0_v = c0 / pow(epsilon + s0_v, 2);
    double alpha1_v = c1 / pow(epsilon + s1_v, 2);
    double alpha2_v = c2 / pow(epsilon + s2_v, 2);
    double sum_alpha_v = alpha0_v + alpha1_v + alpha2_v;
    double w0_v = alpha0_v / sum_alpha_v;
    double w1_v = alpha1_v / sum_alpha_v;
    double w2_v = alpha2_v / sum_alpha_v;

    double ve = w0_v * v0 + w1_v * v1 + w2_v * v2;

    double P0 = (-1.0/6.0)*P[j-2][i] + (5.0/6.0)*P[j-1][i] + (1.0/3.0)*P[j][i];
    double P1 = (1.0/3.0)*P[j-1][i] + (5.0/6.0)*P[j][i] - (1.0/6.0)*P[j+1][i];
    double P2 = (11.0/6.0)*P[j][i] - (7.0/6.0)*P[j+1][i] + (1.0/3.0)*P[j+2][i];

    double s0_P = (13.0/12.0)*pow(P[j-2][i] - 2*P[j-1][i] + P[j][i], 2) +
                  0.25*pow(P[j-2][i] - 4*P[j-1][i] + 3*P[j][i], 2);
    double s1_P = (13.0/12.0)*pow(P[j-1][i] - 2*P[j][i] + P[j+1][i], 2) +
                  0.25*pow(P[j-1][i] - P[j+1][i], 2);
    double s2_P = (13.0/12.0)*pow(P[j][i] - 2*P[j+1][i] + P[j+2][i], 2) +
                  0.25*pow(3*P[j][i] - 4*P[j+1][i] + P[j+2][i], 2);

    double alpha0_P = c0 / pow(epsilon + s0_P, 2);
    double alpha1_P = c1 / pow(epsilon + s1_P, 2);
    double alpha2_P = c2 / pow(epsilon + s2_P, 2);
    double sum_alpha_P = alpha0_P + alpha1_P + alpha2_P;
    double w0_P = alpha0_P / sum_alpha_P;
    double w1_P = alpha1_P / sum_alpha_P;
    double w2_P = alpha2_P / sum_alpha_P;

    double Pe = w0_P * P0 + w1_P * P1 + w2_P * P2;
    if (Pe < 1e-5) Pe = 1e-5;

    double rho0 = (-1.0/6.0)*rho[j-2][i] + (5.0/6.0)*rho[j-1][i] + (1.0/3.0)*rho[j][i];
    double rho1 = (1.0/3.0)*rho[j-1][i] + (5.0/6.0)*rho[j][i] - (1.0/6.0)*rho[j+1][i];
    double rho2 = (11.0/6.0)*rho[j][i] - (7.0/6.0)*rho[j+1][i] + (1.0/3.0)*rho[j+2][i];

    double s0_rho = (13.0/12.0)*pow(rho[j-2][i] - 2*rho[j-1][i] + rho[j][i], 2) +
                    0.25*pow(rho[j-2][i] - 4*rho[j-1][i] + 3*rho[j][i], 2);
    double s1_rho = (13.0/12.0)*pow(rho[j-1][i] - 2*rho[j][i] + rho[j+1][i], 2) +
                    0.25*pow(rho[j-1][i] - rho[j+1][i], 2);
    double s2_rho = (13.0/12.0)*pow(rho[j][i] - 2*rho[j+1][i] + rho[j+2][i], 2) +
                    0.25*pow(3*rho[j][i] - 4*rho[j+1][i] + rho[j+2][i], 2);

    double alpha0_rho = c0 / pow(epsilon + s0_rho, 2);
    double alpha1_rho = c1 / pow(epsilon + s1_rho, 2);
    double alpha2_rho = c2 / pow(epsilon + s2_rho, 2);
    double sum_alpha_rho = alpha0_rho + alpha1_rho + alpha2_rho;
    double w0_rho = alpha0_rho / sum_alpha_rho;
    double w1_rho = alpha1_rho / sum_alpha_rho;
    double w2_rho = alpha2_rho / sum_alpha_rho;

    double rhoe = w0_rho * rho0 + w1_rho * rho1 + w2_rho * rho2;
    if (rhoe < 1e-5) rhoe = 1e-5;

    return std::make_tuple(ue, ve, Pe, rhoe);
}

std::tuple<
    std::vector<std::vector<double>>, // F1
    std::vector<std::vector<double>>, // F2
    std::vector<std::vector<double>>, // F3
    std::vector<std::vector<double>>, // F4
    std::vector<std::vector<double>>, // G1
    std::vector<std::vector<double>>, // G2
    std::vector<std::vector<double>>, // G3
    std::vector<std::vector<double>>  // G4
>
computeFluxes2D_WENO(
    int cells_x, int cells_y, int ghost_cells,
    const std::vector<std::vector<double>>& rho,
    const std::vector<std::vector<double>>& u,
    const std::vector<std::vector<double>>& v,
    const std::vector<std::vector<double>>& P,
    double gamma,
    const std::string& nameRS)
{
    int Nx = cells_x + 2 * ghost_cells;
    int Ny = cells_y + 2 * ghost_cells;

    std::vector<std::vector<double>> F1(Ny, std::vector<double>(Nx - 1, 0.0));
    std::vector<std::vector<double>> F2(Ny, std::vector<double>(Nx - 1, 0.0));
    std::vector<std::vector<double>> F3(Ny, std::vector<double>(Nx - 1, 0.0));
    std::vector<std::vector<double>> F4(Ny, std::vector<double>(Nx - 1, 0.0));

    std::vector<std::vector<double>> G1(Ny - 1, std::vector<double>(Nx, 0.0));
    std::vector<std::vector<double>> G2(Ny - 1, std::vector<double>(Nx, 0.0));
    std::vector<std::vector<double>> G3(Ny - 1, std::vector<double>(Nx, 0.0));
    std::vector<std::vector<double>> G4(Ny - 1, std::vector<double>(Nx, 0.0));

    for (int j = 0; j < Ny; ++j) {
        for (int i = ghost_cells - 1; i < Nx - ghost_cells; ++i) {
            // Elail ninni?icl (cc ??leec i) ? cir?licl ir ddraie nnidiil ??leec i
            auto [u_L, v_L, p_L, rho_L] = WENO5reconstruction_x_right(i, j, u, v, P, rho);
            // Ddrail ninni?icl (cc ??leec i+1) ? cir?licl ir elaie nnidiil ??leec i+1
            auto [u_R, v_R, p_R, rho_R] = WENO5reconstruction_x_left(i+1, j, u, v, P, rho);

            auto [rho_face, u_face, p_face] = superSolve(rho_L, u_L, p_L,
                                                         rho_R, u_R, p_R,
                                                         nameRS);
            double v_star = (u_face >= 0.0) ? v_L : v_R;
            double E_face = p_face / (gamma - 1.0) + 0.5 * rho_face * (u_face*u_face + v_star*v_star);

            F1[j][i] = rho_face * u_face;
            F2[j][i] = rho_face * u_face * u_face + p_face;
            F3[j][i] = rho_face * u_face * v_star;
            F4[j][i] = u_face * (E_face + p_face);
        }
    }

    for (int j = ghost_cells - 1; j < Ny - ghost_cells; ++j) {
        for (int i = 0; i < Nx; ++i) {
            // Iccill ninni?icl (cc ??leec j) ? cir?licl ir aldoile nnidiil ??leec j
            auto [u_B, v_B, p_B, rho_B] = WENO5reconstruction_y_top(i, j, u, v, P, rho);
            // Aldoill ninni?icl (cc ??leec j+1) ? cir?licl ir iccile nnidiil ??leec j+1
            auto [u_T, v_T, p_T, rho_T] = WENO5reconstruction_y_bottom(i, j+1, u, v, P, rho);

            auto [rho_face, v_face, p_face] = superSolve(rho_B, v_B, p_B,
                                                         rho_T, v_T, p_T,
                                                         nameRS);
            double u_star = (v_face >= 0.0) ? u_B : u_T;
            double E_face = p_face / (gamma - 1.0) + 0.5 * rho_face * (u_star*u_star + v_face*v_face);

            G1[j][i] = rho_face * v_face;
            G2[j][i] = rho_face * v_face * u_star;
            G3[j][i] = rho_face * v_face * v_face + p_face;
            G4[j][i] = v_face * (E_face + p_face);
        }
    }

    return {F1, F2, F3, F4, G1, G2, G3, G4};
}

std::tuple<
    std::vector<std::vector<double>>, // R1
    std::vector<std::vector<double>>, // R2
    std::vector<std::vector<double>>, // R3
    std::vector<std::vector<double>>  // R4
>
computeRHS_WENO_2D(
    int cells_x, int cells_y, int ghost_cells,
    const std::vector<std::vector<double>>& rho,
    const std::vector<std::vector<double>>& u,
    const std::vector<std::vector<double>>& v,
    const std::vector<std::vector<double>>& P,
    double hx, double hy,
    double gamma,
    const std::string& nameRS)
{
    int Nx = cells_x + 2 * ghost_cells;
    int Ny = cells_y + 2 * ghost_cells;

    auto [F1, F2, F3, F4, G1, G2, G3, G4] = computeFluxes2D_WENO(
        cells_x, cells_y, ghost_cells,
        rho, u, v, P, gamma, nameRS);

    std::vector<std::vector<double>> R1(Ny, std::vector<double>(Nx, 0.0));
    std::vector<std::vector<double>> R2(Ny, std::vector<double>(Nx, 0.0));
    std::vector<std::vector<double>> R3(Ny, std::vector<double>(Nx, 0.0));
    std::vector<std::vector<double>> R4(Ny, std::vector<double>(Nx, 0.0));

    for (int j = ghost_cells; j < Ny - ghost_cells; ++j) {
        for (int i = ghost_cells; i < Nx - ghost_cells; ++i) {
            double dFx1 = F1[j][i] - F1[j][i-1];
            double dFx2 = F2[j][i] - F2[j][i-1];
            double dFx3 = F3[j][i] - F3[j][i-1];
            double dFx4 = F4[j][i] - F4[j][i-1];

            double dGy1 = G1[j][i] - G1[j-1][i];
            double dGy2 = G2[j][i] - G2[j-1][i];
            double dGy3 = G3[j][i] - G3[j-1][i];
            double dGy4 = G4[j][i] - G4[j-1][i];

            R1[j][i] = -(dFx1 / hx + dGy1 / hy);
            R2[j][i] = -(dFx2 / hx + dGy2 / hy);
            R3[j][i] = -(dFx3 / hx + dGy3 / hy);
            R4[j][i] = -(dFx4 / hx + dGy4 / hy);
        }
    }

    return {R1, R2, R3, R4};
}

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
    const std::string& nameRS)
{
    const double gamma = 1.4;
    int Nx = cells_x + 2 * ghost_cells;
    int Ny = cells_y + 2 * ghost_cells;

    std::vector<std::vector<double>> U1(Ny, std::vector<double>(Nx));
    std::vector<std::vector<double>> U2(Ny, std::vector<double>(Nx));
    std::vector<std::vector<double>> U3(Ny, std::vector<double>(Nx));
    std::vector<std::vector<double>> U4(Ny, std::vector<double>(Nx));

    for (int j = 0; j < Ny; ++j) {
        for (int i = 0; i < Nx; ++i) {
            U1[j][i] = rho[j][i];
            U2[j][i] = rho[j][i] * u[j][i];
            U3[j][i] = rho[j][i] * v[j][i];
            U4[j][i] = P[j][i] / (gamma - 1.0) + 0.5 * rho[j][i] * (u[j][i]*u[j][i] + v[j][i]*v[j][i]);
        }
    }

    applyBCs2D(rho, u, v, P, Nx, Ny, cells_x, cells_y, ghost_cells,
               left_bc, right_bc, bottom_bc, top_bc);

    std::vector<std::vector<double>> U1_new(Ny, std::vector<double>(Nx, 0.0));
    std::vector<std::vector<double>> U2_new(Ny, std::vector<double>(Nx, 0.0));
    std::vector<std::vector<double>> U3_new(Ny, std::vector<double>(Nx, 0.0));
    std::vector<std::vector<double>> U4_new(Ny, std::vector<double>(Nx, 0.0));

    if (time_integrator == "euler") {
        auto [R1, R2, R3, R4] = computeRHS_WENO_2D(
            cells_x, cells_y, ghost_cells,
            rho, u, v, P, hx, hy, gamma, nameRS);

        for (int j = ghost_cells; j < Ny - ghost_cells; ++j) {
            for (int i = ghost_cells; i < Nx - ghost_cells; ++i) {
                U1_new[j][i] = U1[j][i] + tau * R1[j][i];
                U2_new[j][i] = U2[j][i] + tau * R2[j][i];
                U3_new[j][i] = U3[j][i] + tau * R3[j][i];
                U4_new[j][i] = U4[j][i] + tau * R4[j][i];
            }
        }

        updatePrimitiveVariables2D(cells_x, cells_y, ghost_cells,
                                   u, v, P, rho, I,
                                   U1_new, U2_new, U3_new, U4_new, gamma);
    }
    else if (time_integrator == "rk2") {
        auto [R1_1, R2_1, R3_1, R4_1] = computeRHS_WENO_2D(
            cells_x, cells_y, ghost_cells,
            rho, u, v, P, hx, hy, gamma, nameRS);

        std::vector<std::vector<double>> U1_temp(Ny, std::vector<double>(Nx, 0.0));
        std::vector<std::vector<double>> U2_temp(Ny, std::vector<double>(Nx, 0.0));
        std::vector<std::vector<double>> U3_temp(Ny, std::vector<double>(Nx, 0.0));
        std::vector<std::vector<double>> U4_temp(Ny, std::vector<double>(Nx, 0.0));

        for (int j = ghost_cells; j < Ny - ghost_cells; ++j) {
            for (int i = ghost_cells; i < Nx - ghost_cells; ++i) {
                U1_temp[j][i] = U1[j][i] + tau * R1_1[j][i];
                U2_temp[j][i] = U2[j][i] + tau * R2_1[j][i];
                U3_temp[j][i] = U3[j][i] + tau * R3_1[j][i];
                U4_temp[j][i] = U4[j][i] + tau * R4_1[j][i];
            }
        }

        std::vector<std::vector<double>> u_temp = u;
        std::vector<std::vector<double>> v_temp = v;
        std::vector<std::vector<double>> P_temp = P;
        std::vector<std::vector<double>> rho_temp = rho;
        std::vector<std::vector<double>> I_temp = I;

        updatePrimitiveVariables2D(cells_x, cells_y, ghost_cells,
                                   u_temp, v_temp, P_temp, rho_temp, I_temp,
                                   U1_temp, U2_temp, U3_temp, U4_temp, gamma);

        applyBCs2D(rho_temp, u_temp, v_temp, P_temp, Nx, Ny,
                   cells_x, cells_y, ghost_cells,
                   left_bc, right_bc, bottom_bc, top_bc);

        auto [R1_2, R2_2, R3_2, R4_2] = computeRHS_WENO_2D(
            cells_x, cells_y, ghost_cells,
            rho_temp, u_temp, v_temp, P_temp, hx, hy, gamma, nameRS);

        for (int j = ghost_cells; j < Ny - ghost_cells; ++j) {
            for (int i = ghost_cells; i < Nx - ghost_cells; ++i) {
                U1_new[j][i] = U1[j][i] + 0.5 * tau * (R1_1[j][i] + R1_2[j][i]);
                U2_new[j][i] = U2[j][i] + 0.5 * tau * (R2_1[j][i] + R2_2[j][i]);
                U3_new[j][i] = U3[j][i] + 0.5 * tau * (R3_1[j][i] + R3_2[j][i]);
                U4_new[j][i] = U4[j][i] + 0.5 * tau * (R4_1[j][i] + R4_2[j][i]);
            }
        }

        updatePrimitiveVariables2D(cells_x, cells_y, ghost_cells,
                                   u, v, P, rho, I,
                                   U1_new, U2_new, U3_new, U4_new, gamma);
    }
    else {
        std::cerr << "Time integrator " << time_integrator << " not implemented in WENO2DSolve.\n";
        return;
    }

    applyBCs2D(rho, u, v, P, Nx, Ny, cells_x, cells_y, ghost_cells,
               left_bc, right_bc, bottom_bc, top_bc);

    t_total += tau;
    std::cout << "t = " << t_total << ", tau = " << tau
              << ", method = " << time_integrator
              << ", RS = " << nameRS << std::endl;
}


std::tuple<
    std::vector<std::vector<double>>, // R1 (ernnr)
    std::vector<std::vector<double>>, // R2 (x-cedoeun)
    std::vector<std::vector<double>>, // R3 (y-cedoeun)
    std::vector<std::vector<double>>  // R4 (yildac?)
>
computeRHS_Fletcher_2D(
    int cells_x, int cells_y, int ghost_cells,
    const std::vector<std::vector<double>>& rho,
    const std::vector<std::vector<double>>& u,
    const std::vector<std::vector<double>>& v,
    const std::vector<std::vector<double>>& P,
    double hx, double hy, double gamma)
{
    int Nx = cells_x + 2 * ghost_cells;
    int Ny = cells_y + 2 * ghost_cells;

    std::vector<std::vector<double>> R1(Ny, std::vector<double>(Nx, 0.0));
    std::vector<std::vector<double>> R2(Ny, std::vector<double>(Nx, 0.0));
    std::vector<std::vector<double>> R3(Ny, std::vector<double>(Nx, 0.0));
    std::vector<std::vector<double>> R4(Ny, std::vector<double>(Nx, 0.0));


    for (int j = ghost_cells; j < Ny - ghost_cells; ++j) {
        for (int i = ghost_cells; i < Nx - ghost_cells; ++i) {
            double E_ij = P[j][i]/(gamma-1.0) + 0.5*rho[j][i]*(u[j][i]*u[j][i] + v[j][i]*v[j][i]);
            double F1_ij = rho[j][i] * u[j][i];
            double F2_ij = rho[j][i] * u[j][i] * u[j][i] + P[j][i];
            double F3_ij = rho[j][i] * u[j][i] * v[j][i];
            double F4_ij = u[j][i] * (E_ij + P[j][i]);

            double G1_ij = rho[j][i] * v[j][i];
            double G2_ij = rho[j][i] * v[j][i] * u[j][i];
            double G3_ij = rho[j][i] * v[j][i] * v[j][i] + P[j][i];
            double G4_ij = v[j][i] * (E_ij + P[j][i]);

            double E_ip1 = P[j][i+1]/(gamma-1.0) + 0.5*rho[j][i+1]*(u[j][i+1]*u[j][i+1] + v[j][i+1]*v[j][i+1]);
            double F1_ip1 = rho[j][i+1] * u[j][i+1];
            double F2_ip1 = rho[j][i+1] * u[j][i+1] * u[j][i+1] + P[j][i+1];
            double F3_ip1 = rho[j][i+1] * u[j][i+1] * v[j][i+1];
            double F4_ip1 = u[j][i+1] * (E_ip1 + P[j][i+1]);

            double E_jp1 = P[j+1][i]/(gamma-1.0) + 0.5*rho[j+1][i]*(u[j+1][i]*u[j+1][i] + v[j+1][i]*v[j+1][i]);
            double G1_jp1 = rho[j+1][i] * v[j+1][i];
            double G2_jp1 = rho[j+1][i] * v[j+1][i] * u[j+1][i];
            double G3_jp1 = rho[j+1][i] * v[j+1][i] * v[j+1][i] + P[j+1][i];
            double G4_jp1 = v[j+1][i] * (E_jp1 + P[j+1][i]);

            double dFx1 = F1_ip1 - F1_ij;
            double dFx2 = F2_ip1 - F2_ij;
            double dFx3 = F3_ip1 - F3_ij;
            double dFx4 = F4_ip1 - F4_ij;

            double dGy1 = G1_jp1 - G1_ij;
            double dGy2 = G2_jp1 - G2_ij;
            double dGy3 = G3_jp1 - G3_ij;
            double dGy4 = G4_jp1 - G4_ij;

            R1[j][i] = -(dFx1 / hx + dGy1 / hy);
            R2[j][i] = -(dFx2 / hx + dGy2 / hy);
            R3[j][i] = -(dFx3 / hx + dGy3 / hy);
            R4[j][i] = -(dFx4 / hx + dGy4 / hy);
        }
    }

    return {R1, R2, R3, R4};
}

void applyFCT2D(
    std::vector<std::vector<double>>& U1,
    std::vector<std::vector<double>>& U2,
    std::vector<std::vector<double>>& U3,
    std::vector<std::vector<double>>& U4,
    int ghost_cells,
    double viscosity_coeff,
    double anti_diffusion_coeff)
{
    int Ny = U1.size();
    int Nx = U1[0].size();

    std::vector<std::vector<double>> U1_viscous(Ny, std::vector<double>(Nx, 0.0));
    std::vector<std::vector<double>> U2_viscous(Ny, std::vector<double>(Nx, 0.0));
    std::vector<std::vector<double>> U3_viscous(Ny, std::vector<double>(Nx, 0.0));
    std::vector<std::vector<double>> U4_viscous(Ny, std::vector<double>(Nx, 0.0));

    for (int j = ghost_cells; j < Ny - ghost_cells; ++j) {
        for (int i = ghost_cells; i < Nx - ghost_cells; ++i) {
            double lap1 = (U1[j][i+1] - 2.0*U1[j][i] + U1[j][i-1])
                        + (U1[j+1][i] - 2.0*U1[j][i] + U1[j-1][i]);
            double lap2 = (U2[j][i+1] - 2.0*U2[j][i] + U2[j][i-1])
                        + (U2[j+1][i] - 2.0*U2[j][i] + U2[j-1][i]);
            double lap3 = (U3[j][i+1] - 2.0*U3[j][i] + U3[j][i-1])
                        + (U3[j+1][i] - 2.0*U3[j][i] + U3[j-1][i]);
            double lap4 = (U4[j][i+1] - 2.0*U4[j][i] + U4[j][i-1])
                        + (U4[j+1][i] - 2.0*U4[j][i] + U4[j-1][i]);

            U1_viscous[j][i] = U1[j][i] + viscosity_coeff * lap1;
            U2_viscous[j][i] = U2[j][i] + viscosity_coeff * lap2;
            U3_viscous[j][i] = U3[j][i] + viscosity_coeff * lap3;
            U4_viscous[j][i] = U4[j][i] + viscosity_coeff * lap4;
        }
    }

    std::vector<std::vector<double>> U1_final = U1_viscous;
    std::vector<std::vector<double>> U2_final = U2_viscous;
    std::vector<std::vector<double>> U3_final = U3_viscous;
    std::vector<std::vector<double>> U4_final = U4_viscous;

    auto limit_flux = [&](int dir, 
                          const std::vector<std::vector<double>>& U,
                          std::vector<std::vector<double>>& U_out)
    {
        int nx = Nx, ny = Ny;
        std::vector<std::vector<double>> phi(ny, std::vector<double>(nx, 0.0));
        std::vector<std::vector<double>> phi_limited(ny, std::vector<double>(nx, 0.0));

        if (dir == 0) { // di x
            for (int j = ghost_cells; j < ny - ghost_cells; ++j) {
                for (int i = ghost_cells; i < nx - ghost_cells - 1; ++i) {
                    phi[j][i] = anti_diffusion_coeff * (U[j][i+1] - U[j][i]);
                }
            }
            for (int j = ghost_cells; j < ny - ghost_cells; ++j) {
                for (int i = ghost_cells + 1; i < nx - ghost_cells - 2; ++i) {
                    double s = (phi[j][i] >= 0) ? 1.0 : -1.0;
                    double delta_left  = U[j][i] - U[j][i-1];
                    double delta_right = U[j][i+2] - U[j][i+1];
                    double abs_phi = std::abs(phi[j][i]);
                    double min_val = std::max(0.0, std::min({s * delta_left, abs_phi, s * delta_right}));
                    phi_limited[j][i] = s * min_val;
                }
            }
            for (int j = ghost_cells; j < ny - ghost_cells; ++j) {
                for (int i = ghost_cells + 1; i < nx - ghost_cells - 1; ++i) {
                    U_out[j][i] -= (phi_limited[j][i] - phi_limited[j][i-1]);
                }
            }
        }
        else { 
            for (int j = ghost_cells; j < ny - ghost_cells - 1; ++j) {
                for (int i = ghost_cells; i < nx - ghost_cells; ++i) {
                    phi[j][i] = anti_diffusion_coeff * (U[j+1][i] - U[j][i]);
                }
            }
            for (int j = ghost_cells + 1; j < ny - ghost_cells - 2; ++j) {
                for (int i = ghost_cells; i < nx - ghost_cells; ++i) {
                    double s = (phi[j][i] >= 0) ? 1.0 : -1.0;
                    double delta_bottom = U[j][i] - U[j-1][i];
                    double delta_top    = U[j+2][i] - U[j+1][i];
                    double abs_phi = std::abs(phi[j][i]);
                    double min_val = std::max(0.0, std::min({s * delta_bottom, abs_phi, s * delta_top}));
                    phi_limited[j][i] = s * min_val;
                }
            }
            for (int j = ghost_cells + 1; j < ny - ghost_cells - 1; ++j) {
                for (int i = ghost_cells; i < nx - ghost_cells; ++i) {
                    U_out[j][i] -= (phi_limited[j][i] - phi_limited[j-1][i]);
                }
            }
        }
    };

    limit_flux(0, U1_viscous, U1_final);
    limit_flux(1, U1_viscous, U1_final);
    limit_flux(0, U2_viscous, U2_final);
    limit_flux(1, U2_viscous, U2_final);
    limit_flux(0, U3_viscous, U3_final);
    limit_flux(1, U3_viscous, U3_final);
    limit_flux(0, U4_viscous, U4_final);
    limit_flux(1, U4_viscous, U4_final);

    U1 = U1_final;
    U2 = U2_final;
    U3 = U3_final;
    U4 = U4_final;
}

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
    double viscosity_coeff, double anti_diffusion_coeff)
{
    const double gamma = 1.4;
    int Nx = cells_x + 2 * ghost_cells;
    int Ny = cells_y + 2 * ghost_cells;

    std::vector<std::vector<double>> U1(Ny, std::vector<double>(Nx));
    std::vector<std::vector<double>> U2(Ny, std::vector<double>(Nx));
    std::vector<std::vector<double>> U3(Ny, std::vector<double>(Nx));
    std::vector<std::vector<double>> U4(Ny, std::vector<double>(Nx));

    for (int j = 0; j < Ny; ++j) {
        for (int i = 0; i < Nx; ++i) {
            U1[j][i] = rho[j][i];
            U2[j][i] = rho[j][i] * u[j][i];
            U3[j][i] = rho[j][i] * v[j][i];
            U4[j][i] = P[j][i] / (gamma - 1.0) + 0.5 * rho[j][i] * (u[j][i]*u[j][i] + v[j][i]*v[j][i]);
        }
    }

    applyBCs2D(rho, u, v, P, Nx, Ny, cells_x, cells_y, ghost_cells,
               left_bc, right_bc, bottom_bc, top_bc);

    auto [R1_1, R2_1, R3_1, R4_1] = computeRHS_Fletcher_2D(
        cells_x, cells_y, ghost_cells, rho, u, v, P, hx, hy, gamma);

    std::vector<std::vector<double>> U1_pred(Ny, std::vector<double>(Nx, 0.0));
    std::vector<std::vector<double>> U2_pred(Ny, std::vector<double>(Nx, 0.0));
    std::vector<std::vector<double>> U3_pred(Ny, std::vector<double>(Nx, 0.0));
    std::vector<std::vector<double>> U4_pred(Ny, std::vector<double>(Nx, 0.0));

    for (int j = ghost_cells; j < Ny - ghost_cells; ++j) {
        for (int i = ghost_cells; i < Nx - ghost_cells; ++i) {
            U1_pred[j][i] = U1[j][i] + tau * R1_1[j][i];
            U2_pred[j][i] = U2[j][i] + tau * R2_1[j][i];
            U3_pred[j][i] = U3[j][i] + tau * R3_1[j][i];
            U4_pred[j][i] = U4[j][i] + tau * R4_1[j][i];
        }
    }

    std::vector<std::vector<double>> rho_pred = rho;
    std::vector<std::vector<double>> u_pred = u;
    std::vector<std::vector<double>> v_pred = v;
    std::vector<std::vector<double>> P_pred = P;
    std::vector<std::vector<double>> I_pred = I;

    updatePrimitiveVariables2D(cells_x, cells_y, ghost_cells,
                               u_pred, v_pred, P_pred, rho_pred, I_pred,
                               U1_pred, U2_pred, U3_pred, U4_pred, gamma);

    applyBCs2D(rho_pred, u_pred, v_pred, P_pred, Nx, Ny,
               cells_x, cells_y, ghost_cells,
               left_bc, right_bc, bottom_bc, top_bc);

    auto [R1_2, R2_2, R3_2, R4_2] = computeRHS_Fletcher_2D(
        cells_x, cells_y, ghost_cells, rho_pred, u_pred, v_pred, P_pred, hx, hy, gamma);

    std::vector<std::vector<double>> U1_new(Ny, std::vector<double>(Nx, 0.0));
    std::vector<std::vector<double>> U2_new(Ny, std::vector<double>(Nx, 0.0));
    std::vector<std::vector<double>> U3_new(Ny, std::vector<double>(Nx, 0.0));
    std::vector<std::vector<double>> U4_new(Ny, std::vector<double>(Nx, 0.0));

    for (int j = ghost_cells; j < Ny - ghost_cells; ++j) {
        for (int i = ghost_cells; i < Nx - ghost_cells; ++i) {
            U1_new[j][i] = 0.5 * (U1_pred[j][i] + U1[j][i] + tau * R1_2[j][i]);
            U2_new[j][i] = 0.5 * (U2_pred[j][i] + U2[j][i] + tau * R2_2[j][i]);
            U3_new[j][i] = 0.5 * (U3_pred[j][i] + U3[j][i] + tau * R3_2[j][i]);
            U4_new[j][i] = 0.5 * (U4_pred[j][i] + U4[j][i] + tau * R4_2[j][i]);
        }
    }

    applyFCT2D(U1_new, U2_new, U3_new, U4_new, ghost_cells,
               viscosity_coeff, anti_diffusion_coeff);

    updatePrimitiveVariables2D(cells_x, cells_y, ghost_cells,
                               u, v, P, rho, I,
                               U1_new, U2_new, U3_new, U4_new, gamma);

    applyBCs2D(rho, u, v, P, Nx, Ny, cells_x, cells_y, ghost_cells,
               left_bc, right_bc, bottom_bc, top_bc);

    t_total += tau;
    std::cout << "t = " << t_total << ", tau = " << tau
              << ", method = Fletcher2D" << std::endl;
}

void saveToCSV(
    const std::string& filename,
    const std::vector<std::vector<double>>& rho,
    const std::vector<std::vector<double>>& u,
    const std::vector<std::vector<double>>& v,
    const std::vector<std::vector<double>>& P,
    const std::vector<std::vector<double>>& I,
    double hx, double hy,
    int cells_x, int cells_y, int ghost_cells,
    double t)
{

    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Ircaer inedunc? oreer: " << filename << std::endl;
        return;
    }

    file << "# \n";
    file << "x,y,rho,u,v,P,e_internal\n";

    for (int j = ghost_cells; j < cells_y + ghost_cells; ++j) {
        for (int i = ghost_cells; i < cells_x + ghost_cells; ++i) {
            double x = (i - ghost_cells + 0.5) * hx;
            double y = (j - ghost_cells + 0.5) * hy;
            file << x << "," << y << ","
                 << rho[j][i] << "," << u[j][i] << "," << v[j][i] << ","
                 << P[j][i] << "," << I[j][i] << "\n";
        }
    }

    file.close();
}
void newTimeStep2D(
    const std::vector<std::vector<double>>& u,
    const std::vector<std::vector<double>>& v,
    const std::vector<std::vector<double>>& P,
    const std::vector<std::vector<double>>& rho,
    double& tau,
    double hx,
    double hy,
    double CFL,
    int ghost_cells)
{
    const double gamma = 1.4;
    double min_dt = 1e20;

    int ny = static_cast<int>(u.size());
    int nx = static_cast<int>(u[0].size());
// CHECK: CFL_2D
    for (int j = ghost_cells; j < ny - ghost_cells; ++j) {
        for (int i = ghost_cells; i < nx - ghost_cells; ++i) {
            double c = std::sqrt(gamma * P[j][i] / rho[j][i]);
            double dt_x = hx / (std::abs(u[j][i]) + c);
            double dt_y = hy / (std::abs(v[j][i]) + c);
            double dt_local = std::min(dt_x, dt_y);
            if (dt_local < min_dt) {
                min_dt = dt_local;
            }
        }
    }

    double tau_new = CFL * min_dt;
    if (tau_new < 1e-10) tau_new = 1e-10;

    // Lnec tau lu? il crari (dldaue rra) cec iiail cir?licl eliurl nleoulai,
    // ni iaiiae?le tau. A ddincaiie neo?rl innrae?le alc ccelilice.
    if (tau == 0.0 || tau_new < tau) {
        if (tau != 0.0) {
            std::cout << "Time step reduced from " << tau << " to " << tau_new << std::endl;
        }
        tau = tau_new;
    }
}
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
    const std::string& test_case)
{
    int Ny = rho.size();
    int Nx = rho[0].size();

    if (test_case == "sod1" || test_case.empty()) {
        for (int j = 0; j < Ny; ++j) {
            for (int i = 0; i < Nx; ++i) {
                double x = (i - ghost_cells + 0.5) * hx;
                if (x < 0.5 * Lx) {
                    rho[j][i] = 1.0;
                    u[j][i]   = 0.0;
                    v[j][i]   = 0.0;
                    P[j][i]   = 1.0;
                } else {
                    rho[j][i] = 0.125;
                    u[j][i]   = 0.0;
                    v[j][i]   = 0.0;
                    P[j][i]   = 0.1;
                }
                I[j][i] = P[j][i] / (rho[j][i] * (gamma - 1.0));
            }
        }
    }
    else     if (test_case == "sod2" || test_case.empty()) {
        for (int j = 0; j < Ny; ++j) {
            for (int i = 0; i < Nx; ++i) {
                double x = (i - ghost_cells + 0.5) * hx;
                if (x < 0.5 * Lx) {
                    rho[j][i] = 1.0;
                    u[j][i]   = -2.0;
                    v[j][i]   = 0.0;
                    P[j][i]   = 0.4;
                } else {
                    rho[j][i] = 1.0;
                    u[j][i]   = 2.0;
                    v[j][i]   = 0.0;
                    P[j][i]   = 0.4;
                }
                I[j][i] = P[j][i] / (rho[j][i] * (gamma - 1.0));
            }
        }
    }
    else     if (test_case == "sod3" || test_case.empty()) {
        for (int j = 0; j < Ny; ++j) {
            for (int i = 0; i < Nx; ++i) {
                double x = (i - ghost_cells + 0.5) * hx;
                if (x < 0.5 * Lx) {
                    rho[j][i] = 1.0;
                    u[j][i]   = 0.0;
                    v[j][i]   = 0.0;
                    P[j][i]   = 1000.0;
                } else {
                    rho[j][i] = 1.0;
                    u[j][i]   = 0.0;
                    v[j][i]   = 0.0;
                    P[j][i]   = 0.01;
                }
                I[j][i] = P[j][i] / (rho[j][i] * (gamma - 1.0));
            }
        }
    }
    else     if (test_case == "sod4" || test_case.empty()) {
        for (int j = 0; j < Ny; ++j) {
            for (int i = 0; i < Nx; ++i) {
                double x = (i - ghost_cells + 0.5) * hx;
                if (x < 0.5 * Lx) {
                    rho[j][i] = 1.0;
                    u[j][i]   = 0.0;
                    v[j][i]   = 0.0;
                    P[j][i]   = 0.01;
                } else {
                    rho[j][i] = 1.0;
                    u[j][i]   = 0.0;
                    v[j][i]   = 0.0;
                    P[j][i]   = 100.0;
                }
                I[j][i] = P[j][i] / (rho[j][i] * (gamma - 1.0));
            }
        }
    }
    else   if (test_case == "sod5" || test_case.empty()) {
        for (int j = 0; j < Ny; ++j) {
            for (int i = 0; i < Nx; ++i) {
                double x = (i - ghost_cells + 0.5) * hx;
                if (x < 0.5 * Lx) {
                    rho[j][i] = 5.99924;
                    u[j][i]   = 19.5975;
                    v[j][i]   = 0.0;
                    P[j][i]   = 460.894;
                } else {
                    rho[j][i] = 5.99242;
                    u[j][i]   = -6.19633;
                    v[j][i]   = 0.0;
                    P[j][i]   = 46.0950;
                }
                I[j][i] = P[j][i] / (rho[j][i] * (gamma - 1.0));
            }
        }
    }
     
    else   if (test_case == "soll" || test_case.empty()) {
        for (int j = 0; j < Ny; ++j) {
            for (int i = 0; i < Nx; ++i) {
                if ((i*hx-0.5)*(i*hx-0.5)+(j*hy-0.5)*(j*hy-0.5) < 0.125 * Lx) {
                    rho[j][i] = 2;
                    u[j][i]   = 0.1;
                    v[j][i]   = 0.1;
                    P[j][i]   = 2;
                } else {
                    rho[j][i] = 1;
                    u[j][i]   = 0;
                    v[j][i]   = 0.0;
                    P[j][i]   = 1;
                }
                I[j][i] = P[j][i] / (rho[j][i] * (gamma - 1.0));
            }
        }
    }
    else   if (test_case == "forflow" || test_case.empty()) {
        for (int j = 0; j < Ny; ++j) {
            for (int i = 0; i < Nx; ++i) {
                    rho[j][i] = 1;
                    u[j][i]   = 0.0;
                    v[j][i]   = 0.0;
                    P[j][i]   = 1;
                    I[j][i] = P[j][i] / (rho[j][i] * (gamma - 1.0));
            }
        }
    }


    else {
        setInitialConditions(rho, u, v, P, I, Lx, Ly, hx, hy, ghost_cells, gamma, "sod");
    }
}
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
                        std::string& test_case)   // iiaue drdrelnd
{
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Ircaer: il oareinu inedunu oree " << filename << std::endl;
        return false;
    }

    std::string line;
    while (std::getline(file, line)) {
        if (line.empty() || line[0] == '#') continue;

        std::istringstream iss(line);
        std::string key;
        iss >> key;

        if (key == "Lx") iss >> Lx;
        else if (key == "Ly") iss >> Ly;
        else if (key == "cells_x") iss >> cells_x;
        else if (key == "cells_y") iss >> cells_y;
        else if (key == "ghost_cells") iss >> ghost_cells;
        else if (key == "CFL") iss >> CFL;
	else if (key == "tau") iss >> tau;
        else if (key == "t_end") iss >> t_end;
        else if (key == "time_integrator") iss >> time_integrator;
        else if (key == "nameRS") iss >> nameRS;
        else if (key == "solver") iss >> solver;
        else if (key == "viscosity_coeff") iss >> viscosity_coeff;
        else if (key == "anti_diffusion_coeff") iss >> anti_diffusion_coeff;
        else if (key == "left_bc") iss >> left_bc;
        else if (key == "right_bc") iss >> right_bc;
        else if (key == "bottom_bc") iss >> bottom_bc;
        else if (key == "top_bc") iss >> top_bc;
        else if (key == "save_every") iss >> save_every;
        else if (key == "test_case") iss >> test_case;   // iiaue eet?
        else {
        }
    }

    file.close();
    return true;
}


// Вспомогательная функция для одного направления (X или Y) с заданным шагом dt
void flicDirection(
    int cells_x, int cells_y, int ghost_cells,
    std::vector<std::vector<double>>& rho,
    std::vector<std::vector<double>>& u,
    std::vector<std::vector<double>>& v,
    std::vector<std::vector<double>>& P,
    std::vector<std::vector<double>>& I,   // полная энергия на объём
    double dt, double h,                    // шаг по времени и пространству для данного направления
    const std::string& dir,                  // "x" или "y"
    const std::string& left_bc, const std::string& right_bc,
    const std::string& bottom_bc, const std::string& top_bc,
    double gamma)
{
    int Nx = cells_x + 2 * ghost_cells;
    int Ny = cells_y + 2 * ghost_cells;

    // ---- Сохраняем консервативные переменные до начала этапа ----
    std::vector<std::vector<double>> U1(Ny, std::vector<double>(Nx));
    std::vector<std::vector<double>> U2(Ny, std::vector<double>(Nx));
    std::vector<std::vector<double>> U3(Ny, std::vector<double>(Nx));
    std::vector<std::vector<double>> U4(Ny, std::vector<double>(Nx));
    for (int j = 0; j < Ny; ++j) {
        for (int i = 0; i < Nx; ++i) {
            U1[j][i] = rho[j][i];
            U2[j][i] = rho[j][i] * u[j][i];
            U3[j][i] = rho[j][i] * v[j][i];
            U4[j][i] = I[j][i];
        }
    }

    // ---- ЛАГРАНЖЕВ ЭТАП ----
    std::vector<std::vector<double>> U1_tilde = U1;   // плотность не меняется
    std::vector<std::vector<double>> U2_tilde = U2;
    std::vector<std::vector<double>> U3_tilde = U3;
    std::vector<std::vector<double>> U4_tilde = U4;

    if (dir == "x") {
        // Обновление u и I под действием давления по x
        for (int j = ghost_cells; j < Ny - ghost_cells; ++j) {
            for (int i = ghost_cells; i < Nx - ghost_cells; ++i) {
                double rho_cur = U1[j][i];

                // Давление на левой и правой гранях ячейки (i)
                double p_left  = 0.5 * (P[j][i-1] + P[j][i]);
                double p_right = 0.5 * (P[j][i] + P[j][i+1]);

                // Скорость u на гранях (средняя для потока работы давления)
                double u_left  = 0.5 * (u[j][i-1] + u[j][i]);
                double u_right = 0.5 * (u[j][i] + u[j][i+1]);

                // Обновление скорости u (лагранжев шаг)
                double u_new = u[j][i] - dt / rho_cur * (p_right - p_left) / h;
                U2_tilde[j][i] = rho_cur * u_new;   // новый импульс

                // Обновление полной энергии (работа давления)
                double pu_left  = p_left * u_left;
                double pu_right = p_right * u_right;
                double I_new = I[j][i] - dt / rho_cur * (pu_right - pu_left) / h;
                U4_tilde[j][i] = I_new;

                // v не меняется
                U3_tilde[j][i] = U3[j][i];
            }
        }
    } else { // dir == "y"
        // Обновление v и I под действием давления по y
        for (int j = ghost_cells; j < Ny - ghost_cells; ++j) {
            for (int i = ghost_cells; i < Nx - ghost_cells; ++i) {
                double rho_cur = U1[j][i];

                // Давление на нижней и верхней гранях ячейки (j)
                double p_bottom = 0.5 * (P[j-1][i] + P[j][i]);
                double p_top    = 0.5 * (P[j][i] + P[j+1][i]);

                // Скорость v на гранях
                double v_bottom = 0.5 * (v[j-1][i] + v[j][i]);
                double v_top    = 0.5 * (v[j][i] + v[j+1][i]);

                // Обновление скорости v
                double v_new = v[j][i] - dt / rho_cur * (p_top - p_bottom) / h;
                U3_tilde[j][i] = rho_cur * v_new;

                // Обновление полной энергии
                double pv_bottom = p_bottom * v_bottom;
                double pv_top    = p_top * v_top;
                double I_new = I[j][i] - dt / rho_cur * (pv_top - pv_bottom) / h;
                U4_tilde[j][i] = I_new;

                // u не меняется
                U2_tilde[j][i] = U2[j][i];
            }
        }
    }

    // ---- ЭЙЛЕРОВ ЭТАП (перенос через границы) ----
    std::vector<std::vector<double>> U1_new = U1_tilde;
    std::vector<std::vector<double>> U2_new = U2_tilde;
    std::vector<std::vector<double>> U3_new = U3_tilde;
    std::vector<std::vector<double>> U4_new = U4_tilde;

    if (dir == "x") {
        // Вычисление потоков через вертикальные грани (i+1/2)
        std::vector<std::vector<double>> flux_mass(Ny, std::vector<double>(Nx-1, 0.0));
        std::vector<std::vector<double>> flux_imp_x(Ny, std::vector<double>(Nx-1, 0.0));
        std::vector<std::vector<double>> flux_imp_y(Ny, std::vector<double>(Nx-1, 0.0));
        std::vector<std::vector<double>> flux_energy(Ny, std::vector<double>(Nx-1, 0.0));

        for (int j = ghost_cells; j < Ny - ghost_cells; ++j) {
            for (int i = ghost_cells; i < Nx - ghost_cells - 1; ++i) {
                // Скорость на грани (средняя после лагранжева этапа)
                double u_face = 0.5 * (U2_tilde[j][i] / U1_tilde[j][i] + U2_tilde[j][i+1] / U1_tilde[j][i+1]);

                // Определяем донорскую ячейку
                int donor = (u_face > 0.0) ? i : i+1;

                // Значения из донора
                double rho_d  = U1_tilde[j][donor];
                double u_d    = U2_tilde[j][donor] / rho_d;
                double v_d    = U3_tilde[j][donor] / rho_d;
                double I_d    = U4_tilde[j][donor];

                // Потоки
                double factor = u_face * dt;
                flux_mass[j][i]   = rho_d * factor;
                flux_imp_x[j][i]  = (rho_d * u_d) * factor;
                flux_imp_y[j][i]  = (rho_d * v_d) * factor;
                flux_energy[j][i] = I_d * factor;
            }
        }

        // Обновление консервативных переменных в ячейках
        for (int j = ghost_cells; j < Ny - ghost_cells; ++j) {
            for (int i = ghost_cells; i < Nx - ghost_cells; ++i) {
                U1_new[j][i] = U1_tilde[j][i] - (flux_mass[j][i] - flux_mass[j][i-1]) / h;
                U2_new[j][i] = U2_tilde[j][i] - (flux_imp_x[j][i] - flux_imp_x[j][i-1]) / h;
                U3_new[j][i] = U3_tilde[j][i] - (flux_imp_y[j][i] - flux_imp_y[j][i-1]) / h;
                U4_new[j][i] = U4_tilde[j][i] - (flux_energy[j][i] - flux_energy[j][i-1]) / h;
            }
        }
    } else { // dir == "y"
        // Вычисление потоков через горизонтальные грани (j+1/2)
        std::vector<std::vector<double>> flux_mass(Ny-1, std::vector<double>(Nx, 0.0));
        std::vector<std::vector<double>> flux_imp_x(Ny-1, std::vector<double>(Nx, 0.0));
        std::vector<std::vector<double>> flux_imp_y(Ny-1, std::vector<double>(Nx, 0.0));
        std::vector<std::vector<double>> flux_energy(Ny-1, std::vector<double>(Nx, 0.0));

        for (int j = ghost_cells; j < Ny - ghost_cells - 1; ++j) {
            for (int i = ghost_cells; i < Nx - ghost_cells; ++i) {
                // Скорость на грани (средняя после лагранжева этапа)
                double v_face = 0.5 * (U3_tilde[j][i] / U1_tilde[j][i] + U3_tilde[j+1][i] / U1_tilde[j+1][i]);

                int donor = (v_face > 0.0) ? j : j+1;

                double rho_d  = U1_tilde[donor][i];
                double u_d    = U2_tilde[donor][i] / rho_d;
                double v_d    = U3_tilde[donor][i] / rho_d;
                double I_d    = U4_tilde[donor][i];

                double factor = v_face * dt;
                flux_mass[j][i]   = rho_d * factor;
                flux_imp_x[j][i]  = (rho_d * u_d) * factor;
                flux_imp_y[j][i]  = (rho_d * v_d) * factor;
                flux_energy[j][i] = I_d * factor;
            }
        }

        // Обновление консервативных переменных в ячейках
        for (int j = ghost_cells; j < Ny - ghost_cells; ++j) {
            for (int i = ghost_cells; i < Nx - ghost_cells; ++i) {
                U1_new[j][i] = U1_tilde[j][i] - (flux_mass[j][i] - flux_mass[j-1][i]) / h;
                U2_new[j][i] = U2_tilde[j][i] - (flux_imp_x[j][i] - flux_imp_x[j-1][i]) / h;
                U3_new[j][i] = U3_tilde[j][i] - (flux_imp_y[j][i] - flux_imp_y[j-1][i]) / h;
                U4_new[j][i] = U4_tilde[j][i] - (flux_energy[j][i] - flux_energy[j-1][i]) / h;
            }
        }
    }

    // ---- Восстановление примитивных переменных ----
    const double eps_rho = 1e-5, eps_e = 1e-6;
    for (int j = ghost_cells; j < Ny - ghost_cells; ++j) {
        for (int i = ghost_cells; i < Nx - ghost_cells; ++i) {
            rho[j][i] = std::max(U1_new[j][i], eps_rho);
            u[j][i]   = U2_new[j][i] / rho[j][i];
            v[j][i]   = U3_new[j][i] / rho[j][i];
            double Ek = 0.5 * rho[j][i] * (u[j][i]*u[j][i] + v[j][i]*v[j][i]);
            double e_int = (U4_new[j][i] - Ek) / rho[j][i];
            if (e_int < eps_e) e_int = eps_e;
            P[j][i] = (gamma - 1.0) * rho[j][i] * e_int;
            I[j][i] = U4_new[j][i];
        }
    }

    // Применение граничных условий
    applyBCs2D(rho, u, v, P, Nx, Ny, cells_x, cells_y, ghost_cells,
               left_bc, right_bc, bottom_bc, top_bc);
}

// Основная функция двумерного метода FLIC с расщеплением Strang
void FLIC2DSolve(
    int cells_x, int cells_y, int ghost_cells,
    std::vector<std::vector<double>>& u,
    std::vector<std::vector<double>>& v,
    std::vector<std::vector<double>>& P,
    std::vector<std::vector<double>>& rho,
    std::vector<std::vector<double>>& I,
    double& t_total, double tau, double hx, double hy,
    const std::string& left_bc, const std::string& right_bc,
    const std::string& bottom_bc, const std::string& top_bc,
    double gamma,
    int step)  // номер текущего шага
{
    double tau_half = tau * 0.5;

    // Чередуем порядок направлений в зависимости от чётности шага
    if (step % 2 == 0) {
        // Чётный шаг: X - Y 
        flicDirection(cells_x, cells_y, ghost_cells,
                      rho, u, v, P, I,
                      tau_half, hx, "x",
                      left_bc, right_bc, bottom_bc, top_bc, gamma);

        flicDirection(cells_x, cells_y, ghost_cells,
                      rho, u, v, P, I,
                      tau, hy, "y",
                      left_bc, right_bc, bottom_bc, top_bc, gamma);

        flicDirection(cells_x, cells_y, ghost_cells,
                      rho, u, v, P, I,
                      tau_half, hx, "x",
                      left_bc, right_bc, bottom_bc, top_bc, gamma);

    } else {
        // Нечётный шаг: Y - X 
        flicDirection(cells_x, cells_y, ghost_cells,
                      rho, u, v, P, I,
                      tau_half, hy, "y",
                      left_bc, right_bc, bottom_bc, top_bc, gamma);

        flicDirection(cells_x, cells_y, ghost_cells,
                      rho, u, v, P, I,
                      tau, hx, "x",
                      left_bc, right_bc, bottom_bc, top_bc, gamma);
        flicDirection(cells_x, cells_y, ghost_cells,
                      rho, u, v, P, I,
                      tau_half, hx, "x",
                      left_bc, right_bc, bottom_bc, top_bc, gamma);

    }

    t_total += tau;
    std::cout << "t = " << t_total << ", tau = " << tau
              << ", method = FLIC2D" << std::endl;
}