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
    std::vector<std::vector<double>>, // F1 (масса, x-поток)
    std::vector<std::vector<double>>, // F2 (x-импульс, x-поток)
    std::vector<std::vector<double>>, // F3 (y-импульс, x-поток)
    std::vector<std::vector<double>>, // F4 (энергия, x-поток)
    std::vector<std::vector<double>>, // G1 (масса, y-поток)
    std::vector<std::vector<double>>, // G2 (x-импульс, y-поток)
    std::vector<std::vector<double>>, // G3 (y-импульс, y-поток)
    std::vector<std::vector<double>>  // G4 (энергия, y-поток)
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
    // Применяем условия по x (строки)
    for (int j = 0; j < Ny; ++j) {
        applyBoundaryConditions(rho[j].data(), u[j].data(), P[j].data(),
                                cells_x, ghost_cells, left_bc, right_bc);
    }
    // Применяем условия по y (столбцы)
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
    std::vector<std::vector<double>>, // R1 (масса)
    std::vector<std::vector<double>>, // R2 (x-импульс)
    std::vector<std::vector<double>>, // R3 (y-импульс)
    std::vector<std::vector<double>>  // R4 (энергия)
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
        std::vector<std::vector<double>> I_temp = I; // не используется, но для совместимости

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
    double c0 = 0.3, c1 = 0.6, c2 = 0.1;  // веса для верхней грани (как для правой)

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
    double c0 = 0.1, c1 = 0.6, c2 = 0.3;  // веса для нижней грани (как для левой)

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
            // Левое состояние (из ячейки i) – значение на правой стороне ячейки i
            auto [u_L, v_L, p_L, rho_L] = WENO5reconstruction_x_right(i, j, u, v, P, rho);
            // Правое состояние (из ячейки i+1) – значение на левой стороне ячейки i+1
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
            // Нижнее состояние (из ячейки j) – значение на верхней стороне ячейки j
            auto [u_B, v_B, p_B, rho_B] = WENO5reconstruction_y_top(i, j, u, v, P, rho);
            // Верхнее состояние (из ячейки j+1) – значение на нижней стороне ячейки j+1
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
    std::vector<std::vector<double>>, // R1 (масса)
    std::vector<std::vector<double>>, // R2 (x-импульс)
    std::vector<std::vector<double>>, // R3 (y-импульс)
    std::vector<std::vector<double>>  // R4 (энергия)
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

        if (dir == 0) { // по x
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
        std::cerr << "Ошибка открытия файла: " << filename << std::endl;
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

    // Если tau ещё не задан (первый шаг) или новое значение меньше текущего,
    // то обновляем tau. В противном случае оставляем без изменений.
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
    else   if (test_case == "Teilor" || test_case.empty()) {
        for (int j = 0; j < Ny; ++j) {
            for (int i = 0; i < Nx; ++i) {
                    rho[j][i] = 1;
                    u[j][i]   = -std::cos(3.1415*i*hx)*std::sin(3.1415*i*hy);
                    v[j][i]   = std::sin(3.1415*i*hx)*std::cos(3.1415*i*hy);
                    P[j][i]   = -1/4*(std::cos(2*3.1415*i*hx)+std::cos(2*3.1415*i*hy));
                    I[j][i] = P[j][i] / (rho[j][i] * (gamma - 1.0));
            }
        }
    }

else if (test_case == "detonation") {
    // Параметры идеализированного ВВ (? = 3, данные из теста)
    const double gamma = 3.0;                // показатель адиабаты продуктов
    const double rho0 = 1.84;                 // начальная плотность [г/см?]
    const double D = 0.88;                     // скорость детонации Чепмена–Жуге [см/мкс]
    const double Q = D*D / (2.0 * (gamma*gamma - 1.0)); // энергия взрыва [Мбар·см?/г]
    const double P_CJ = rho0 * D*D / (gamma + 1.0);      // давление на плоскости CJ
    const double rho_CJ = rho0 * (gamma + 1.0) / gamma;  // плотность за фронтом
    const double u_CJ = D / (gamma + 1.0);                // скорость вещества за фронтом
    const double I_thermal_CJ = P_CJ / ((gamma - 1.0) * rho_CJ); // тепловая часть энергии
    const double I_CJ = I_thermal_CJ + Q;                 // полная внутренняя энергия

    // Ширина детонатора (например, 5 ячеек)
    int det_cells = 1;
    double det_length = det_cells * hx;

    for (int j = 0; j < Ny; ++j) {
        for (int i = 0; i < Nx; ++i) {
            double x = (i - ghost_cells + 0.5) * hx; // координата центра ячейки
            if (x < det_length) {
                // Детонатор (продукты)
                rho[j][i] = rho_CJ;
                u[j][i]   = u_CJ;
                v[j][i]   = 0.0;
                P[j][i]   = P_CJ;
                I[j][i]   = I_CJ;
            } else {
                // Невзорвавшееся ВВ
                rho[j][i] = rho0;
                u[j][i]   = 0.0;
                v[j][i]   = 0.0;
                P[j][i]   = 1e-6; // очень маленькое давление, чтобы избежать деления на ноль
                I[j][i]   = 0.0;
            }
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
                        std::string& test_case)   // новый параметр
{
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Ошибка: не удалось открыть файл " << filename << std::endl;
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
        else if (key == "test_case") iss >> test_case;   // новый ключ
        else {
        }
    }

    file.close();
    return true;
}

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
    bool order_xy){ 

    const double gamma = 1.4;
    int Nx = cells_x + 2 * ghost_cells;
    int Ny = cells_y + 2 * ghost_cells;

    std::vector<std::vector<double>> E(Ny, std::vector<double>(Nx));
    for (int j = 0; j < Ny; ++j){
        for (int i = 0; i < Nx; ++i){
            E[j][i] = I[j][i] + 0.5 * (u[j][i]*u[j][i] + v[j][i]*v[j][i]);}}

    for (int j = 0; j < Ny; ++j) {
        applyBoundaryConditions(rho[j].data(), u[j].data(), P[j].data(),
                                cells_x, ghost_cells, left_bc, right_bc);
    }
    for (int i = 0; i < Nx; ++i) {
        std::vector<double> rho_col(Ny), v_col(Ny), P_col(Ny);
        for (int j = 0; j < Ny; ++j) {
            rho_col[j] = rho[j][i];
            v_col[j]   = v[j][i];
            P_col[j]   = P[j][i];
        }
        applyBoundaryConditions(rho_col.data(), v_col.data(), P_col.data(),
                                cells_y, ghost_cells, bottom_bc, top_bc);
        for (int j = 0; j < Ny; ++j) {
            rho[j][i] = rho_col[j];
            v[j][i]   = v_col[j];
            P[j][i]   = P_col[j];
        }
    }

    if (order_xy) {
        std::vector<std::vector<double>> u_tilde(Ny, std::vector<double>(Nx));
        std::vector<std::vector<double>> v_tilde(Ny, std::vector<double>(Nx));
        std::vector<std::vector<double>> E_tilde(Ny, std::vector<double>(Nx));

        for (int j = ghost_cells; j < Ny - ghost_cells; ++j) {
            for (int i = ghost_cells; i < Nx - ghost_cells; ++i) {
                double pL = 0.5 * (P[j][i] + P[j][i-1]);
                double pR = 0.5 * (P[j][i] + P[j][i+1]);
                double uL = 0.5 * (u[j][i] + u[j][i-1]);
                double uR = 0.5 * (u[j][i] + u[j][i+1]);

                u_tilde[j][i] = u[j][i] - (tau / rho[j][i]) * (pR - pL) / hx;
                v_tilde[j][i] = v[j][i];
                E_tilde[j][i] = E[j][i] - (tau / rho[j][i]) * (pR * uR - pL * uL) / hx;
            }
        }

        for (int j = 0; j < Ny; ++j) {
            applyBoundaryConditions(rho[j].data(), u_tilde[j].data(), P[j].data(),
                                    cells_x, ghost_cells, left_bc, right_bc);
        }
        for (int i = 0; i < Nx; ++i) {
            std::vector<double> rho_col(Ny), u_col(Ny), v_col(Ny), E_col(Ny), P_col(Ny);
            for (int j = 0; j < Ny; ++j) {
                rho_col[j] = rho[j][i];
                u_col[j]   = u_tilde[j][i];
                v_col[j]   = v_tilde[j][i];
                E_col[j]   = E_tilde[j][i];
                P_col[j]   = P[j][i];
            }
            applyBoundaryConditions(rho_col.data(), u_col.data(), P_col.data(),
                                    cells_y, ghost_cells, bottom_bc, top_bc);
            applyBoundaryConditions(rho_col.data(), v_col.data(), P_col.data(),
                                    cells_y, ghost_cells, bottom_bc, top_bc);
            applyBoundaryConditions(rho_col.data(), E_col.data(), P_col.data(),
                                    cells_y, ghost_cells, bottom_bc, top_bc);
            for (int j = 0; j < Ny; ++j) {
                u_tilde[j][i] = u_col[j];
                v_tilde[j][i] = v_col[j];
                E_tilde[j][i] = E_col[j];
            }
        }

        std::vector<std::vector<double>> rho_u_tilde(Ny, std::vector<double>(Nx));
        std::vector<std::vector<double>> rho_v_tilde(Ny, std::vector<double>(Nx));
        std::vector<std::vector<double>> rho_E_tilde(Ny, std::vector<double>(Nx));
        for (int j = 0; j < Ny; ++j)
            for (int i = 0; i < Nx; ++i) {
                rho_u_tilde[j][i] = rho[j][i] * u_tilde[j][i];
                rho_v_tilde[j][i] = rho[j][i] * v_tilde[j][i];
                rho_E_tilde[j][i] = rho[j][i] * E_tilde[j][i];
            }

        std::vector<std::vector<double>> rho_new(Ny, std::vector<double>(Nx));
        std::vector<std::vector<double>> rho_u_new(Ny, std::vector<double>(Nx));
        std::vector<std::vector<double>> rho_v_new(Ny, std::vector<double>(Nx));
        std::vector<std::vector<double>> rho_E_new(Ny, std::vector<double>(Nx));
        for (int j = 0; j < Ny; ++j)
            for (int i = 0; i < Nx; ++i) {
                rho_new[j][i]    = rho[j][i];
                rho_u_new[j][i]  = rho_u_tilde[j][i];
                rho_v_new[j][i]  = rho_v_tilde[j][i];
                rho_E_new[j][i]  = rho_E_tilde[j][i];
            }

        for (int j = ghost_cells; j < Ny - ghost_cells; ++j) {
            for (int i = ghost_cells; i < Nx - ghost_cells + 1; ++i) {
                double u_face = 0.5 * (u_tilde[j][i-1] + u_tilde[j][i]);
                int donor_i = (u_face >= 0) ? i-1 : i;

                double flux_mass   = rho[j][donor_i] * u_face;
                double flux_rho_u  = rho_u_tilde[j][donor_i] * u_face;
                double flux_rho_v  = rho_v_tilde[j][donor_i] * u_face;
                double flux_rho_E  = rho_E_tilde[j][donor_i] * u_face;

                if (i >= ghost_cells && i < Nx - ghost_cells) {
                    rho_new[j][i]    += (tau / hx) * flux_mass;
                    rho_u_new[j][i]  += (tau / hx) * flux_rho_u;
                    rho_v_new[j][i]  += (tau / hx) * flux_rho_v;
                    rho_E_new[j][i]  += (tau / hx) * flux_rho_E;
                }
                if (i-1 >= ghost_cells && i-1 < Nx - ghost_cells) {
                    rho_new[j][i-1]    -= (tau / hx) * flux_mass;
                    rho_u_new[j][i-1]  -= (tau / hx) * flux_rho_u;
                    rho_v_new[j][i-1]  -= (tau / hx) * flux_rho_v;
                    rho_E_new[j][i-1]  -= (tau / hx) * flux_rho_E;
                }
            }
        }

        for (int j = ghost_cells; j < Ny - ghost_cells; ++j)
            for (int i = ghost_cells; i < Nx - ghost_cells; ++i) {
                double rho_inv = 1.0 / rho_new[j][i];
                u[j][i] = rho_u_new[j][i] * rho_inv;
                v[j][i] = rho_v_new[j][i] * rho_inv;
                E[j][i] = rho_E_new[j][i] * rho_inv;
                I[j][i] = E[j][i] - 0.5*(u[j][i]*u[j][i] + v[j][i]*v[j][i]);
                P[j][i] = (gamma - 1.0) * rho_new[j][i] * I[j][i];
                rho[j][i] = rho_new[j][i];
            }

        for (int j = 0; j < Ny; ++j)
            applyBoundaryConditions(rho[j].data(), u[j].data(), P[j].data(),
                                    cells_x, ghost_cells, left_bc, right_bc);
        for (int i = 0; i < Nx; ++i) {
            std::vector<double> rho_col(Ny), v_col(Ny), P_col(Ny);
            for (int j = 0; j < Ny; ++j) {
                rho_col[j] = rho[j][i];
                v_col[j]   = v[j][i];
                P_col[j]   = P[j][i];
            }
            applyBoundaryConditions(rho_col.data(), v_col.data(), P_col.data(),
                                    cells_y, ghost_cells, bottom_bc, top_bc);
            for (int j = 0; j < Ny; ++j) {
                rho[j][i] = rho_col[j];
                v[j][i]   = v_col[j];
                P[j][i]   = P_col[j];
            }
        }
    for (int j = 0; j < Ny; ++j){
        for (int i = 0; i < Nx; ++i){
            E[j][i] = I[j][i] + 0.5 * (u[j][i]*u[j][i] + v[j][i]*v[j][i]);}}

        std::vector<std::vector<double>> u_tilde2(Ny, std::vector<double>(Nx));
        std::vector<std::vector<double>> v_tilde2(Ny, std::vector<double>(Nx));
        std::vector<std::vector<double>> E_tilde2(Ny, std::vector<double>(Nx));

        for (int j = ghost_cells; j < Ny - ghost_cells; ++j)
            for (int i = ghost_cells; i < Nx - ghost_cells; ++i) {
                double pB = 0.5 * (P[j][i] + P[j-1][i]);
                double pT = 0.5 * (P[j][i] + P[j+1][i]);
                double vB = 0.5 * (v[j][i] + v[j-1][i]);
                double vT = 0.5 * (v[j][i] + v[j+1][i]);

                u_tilde2[j][i] = u[j][i];
                v_tilde2[j][i] = v[j][i] - (tau / rho[j][i]) * (pT - pB) / hy;
                E_tilde2[j][i] = E[j][i] - (tau / rho[j][i]) * (pT * vT - pB * vB) / hy;
            }

        for (int j = 0; j < Ny; ++j) {
            applyBoundaryConditions(rho[j].data(), u_tilde2[j].data(), P[j].data(),
                                    cells_x, ghost_cells, left_bc, right_bc);
        }
        for (int i = 0; i < Nx; ++i) {
            std::vector<double> rho_col(Ny), u_col(Ny), v_col(Ny), E_col(Ny), P_col(Ny);
            for (int j = 0; j < Ny; ++j) {
                rho_col[j] = rho[j][i];
                u_col[j]   = u_tilde2[j][i];
                v_col[j]   = v_tilde2[j][i];
                E_col[j]   = E_tilde2[j][i];
                P_col[j]   = P[j][i];
            }
            applyBoundaryConditions(rho_col.data(), u_col.data(), P_col.data(),
                                    cells_y, ghost_cells, bottom_bc, top_bc);
            applyBoundaryConditions(rho_col.data(), v_col.data(), P_col.data(),
                                    cells_y, ghost_cells, bottom_bc, top_bc);
            applyBoundaryConditions(rho_col.data(), E_col.data(), P_col.data(),
                                    cells_y, ghost_cells, bottom_bc, top_bc);
            for (int j = 0; j < Ny; ++j) {
                u_tilde2[j][i] = u_col[j];
                v_tilde2[j][i] = v_col[j];
                E_tilde2[j][i] = E_col[j];
            }
        }

        for (int j = 0; j < Ny; ++j)
            for (int i = 0; i < Nx; ++i) {
                rho_u_tilde[j][i] = rho[j][i] * u_tilde2[j][i];
                rho_v_tilde[j][i] = rho[j][i] * v_tilde2[j][i];
                rho_E_tilde[j][i] = rho[j][i] * E_tilde2[j][i];
            }

        for (int j = 0; j < Ny; ++j)
            for (int i = 0; i < Nx; ++i) {
                rho_new[j][i]    = rho[j][i];
                rho_u_new[j][i]  = rho_u_tilde[j][i];
                rho_v_new[j][i]  = rho_v_tilde[j][i];
                rho_E_new[j][i]  = rho_E_tilde[j][i];
            }

        for (int j = ghost_cells; j < Ny - ghost_cells + 1; ++j)
            for (int i = ghost_cells; i < Nx - ghost_cells; ++i) {
                double v_face = 0.5 * (v_tilde2[j-1][i] + v_tilde2[j][i]);
                int donor_j = (v_face >= 0) ? j-1 : j;

                double flux_mass   = rho[donor_j][i] * v_face;
                double flux_rho_u  = rho_u_tilde[donor_j][i] * v_face;
                double flux_rho_v  = rho_v_tilde[donor_j][i] * v_face;
                double flux_rho_E  = rho_E_tilde[donor_j][i] * v_face;

                if (j >= ghost_cells && j < Ny - ghost_cells) {
                    rho_new[j][i]    += (tau / hy) * flux_mass;
                    rho_u_new[j][i]  += (tau / hy) * flux_rho_u;
                    rho_v_new[j][i]  += (tau / hy) * flux_rho_v;
                    rho_E_new[j][i]  += (tau / hy) * flux_rho_E;
                }
                if (j-1 >= ghost_cells && j-1 < Ny - ghost_cells) {
                    rho_new[j-1][i]    -= (tau / hy) * flux_mass;
                    rho_u_new[j-1][i]  -= (tau / hy) * flux_rho_u;
                    rho_v_new[j-1][i]  -= (tau / hy) * flux_rho_v;
                    rho_E_new[j-1][i]  -= (tau / hy) * flux_rho_E;
                }
            }

        for (int j = ghost_cells; j < Ny - ghost_cells; ++j)
            for (int i = ghost_cells; i < Nx - ghost_cells; ++i) {
                double rho_inv = 1.0 / rho_new[j][i];
                u[j][i] = rho_u_new[j][i] * rho_inv;
                v[j][i] = rho_v_new[j][i] * rho_inv;
                E[j][i] = rho_E_new[j][i] * rho_inv;
                I[j][i] = E[j][i] - 0.5*(u[j][i]*u[j][i] + v[j][i]*v[j][i]);
                P[j][i] = (gamma - 1.0) * rho_new[j][i] * I[j][i];
                rho[j][i] = rho_new[j][i];
            }
    } else {
        std::vector<std::vector<double>> u_tildeY(Ny, std::vector<double>(Nx));
        std::vector<std::vector<double>> v_tildeY(Ny, std::vector<double>(Nx));
        std::vector<std::vector<double>> E_tildeY(Ny, std::vector<double>(Nx));
        std::vector<std::vector<double>> rho_u_tilde(Ny, std::vector<double>(Nx));
        std::vector<std::vector<double>> rho_v_tilde(Ny, std::vector<double>(Nx));
        std::vector<std::vector<double>> rho_E_tilde(Ny, std::vector<double>(Nx));
        std::vector<std::vector<double>> rho_new(Ny, std::vector<double>(Nx));
        std::vector<std::vector<double>> rho_u_new(Ny, std::vector<double>(Nx));
        std::vector<std::vector<double>> rho_v_new(Ny, std::vector<double>(Nx));
        std::vector<std::vector<double>> rho_E_new(Ny, std::vector<double>(Nx));

        for (int j = ghost_cells; j < Ny - ghost_cells; ++j)
            for (int i = ghost_cells; i < Nx - ghost_cells; ++i) {
                double pB = 0.5 * (P[j][i] + P[j-1][i]);
                double pT = 0.5 * (P[j][i] + P[j+1][i]);
                double vB = 0.5 * (v[j][i] + v[j-1][i]);
                double vT = 0.5 * (v[j][i] + v[j+1][i]);

                u_tildeY[j][i] = u[j][i];
                v_tildeY[j][i] = v[j][i] - (tau / rho[j][i]) * (pT - pB) / hy;
                E_tildeY[j][i] = E[j][i] - (tau / rho[j][i]) * (pT * vT - pB * vB) / hy;
            }

        for (int j = 0; j < Ny; ++j) {
            applyBoundaryConditions(rho[j].data(), u_tildeY[j].data(), P[j].data(),
                                    cells_x, ghost_cells, left_bc, right_bc);
         }
        for (int i = 0; i < Nx; ++i) {
            std::vector<double> rho_col(Ny), u_col(Ny), v_col(Ny), E_col(Ny), P_col(Ny);
            for (int j = 0; j < Ny; ++j) {
                rho_col[j] = rho[j][i];
                u_col[j]   = u_tildeY[j][i];
                v_col[j]   = v_tildeY[j][i];
                E_col[j]   = E_tildeY[j][i];
                P_col[j]   = P[j][i];
            }
            applyBoundaryConditions(rho_col.data(), u_col.data(), P_col.data(),
                                    cells_y, ghost_cells, bottom_bc, top_bc);
            applyBoundaryConditions(rho_col.data(), v_col.data(), P_col.data(),
                                    cells_y, ghost_cells, bottom_bc, top_bc);
            applyBoundaryConditions(rho_col.data(), E_col.data(), P_col.data(),
                                    cells_y, ghost_cells, bottom_bc, top_bc);
            for (int j = 0; j < Ny; ++j) {
                u_tildeY[j][i] = u_col[j];
                v_tildeY[j][i] = v_col[j];
                E_tildeY[j][i] = E_col[j];
            }
        }

        for (int j = 0; j < Ny; ++j)
            for (int i = 0; i < Nx; ++i) {
                rho_u_tilde[j][i] = rho[j][i] * u_tildeY[j][i];
                rho_v_tilde[j][i] = rho[j][i] * v_tildeY[j][i];
                rho_E_tilde[j][i] = rho[j][i] * E_tildeY[j][i];
            }

        for (int j = 0; j < Ny; ++j)
            for (int i = 0; i < Nx; ++i) {
                rho_new[j][i]    = rho[j][i];
                rho_u_new[j][i]  = rho_u_tilde[j][i];
                rho_v_new[j][i]  = rho_v_tilde[j][i];
                rho_E_new[j][i]  = rho_E_tilde[j][i];
            }

        for (int j = ghost_cells; j < Ny - ghost_cells + 1; ++j)
            for (int i = ghost_cells; i < Nx - ghost_cells; ++i) {
                double v_face = 0.5 * (v_tildeY[j-1][i] + v_tildeY[j][i]);
                int donor_j = (v_face >= 0) ? j-1 : j;

                double flux_mass   = rho[donor_j][i] * v_face;
                double flux_rho_u  = rho_u_tilde[donor_j][i] * v_face;
                double flux_rho_v  = rho_v_tilde[donor_j][i] * v_face;
                double flux_rho_E  = rho_E_tilde[donor_j][i] * v_face;

                if (j >= ghost_cells && j < Ny - ghost_cells) {
                    rho_new[j][i]    += (tau / hy) * flux_mass;
                    rho_u_new[j][i]  += (tau / hy) * flux_rho_u;
                    rho_v_new[j][i]  += (tau / hy) * flux_rho_v;
                    rho_E_new[j][i]  += (tau / hy) * flux_rho_E;
                }
                if (j-1 >= ghost_cells && j-1 < Ny - ghost_cells) {
                    rho_new[j-1][i]    -= (tau / hy) * flux_mass;
                    rho_u_new[j-1][i]  -= (tau / hy) * flux_rho_u;
                    rho_v_new[j-1][i]  -= (tau / hy) * flux_rho_v;
                    rho_E_new[j-1][i]  -= (tau / hy) * flux_rho_E;
                }
            }

        for (int j = ghost_cells; j < Ny - ghost_cells; ++j)
            for (int i = ghost_cells; i < Nx - ghost_cells; ++i) {
                double rho_inv = 1.0 / rho_new[j][i];
                u[j][i] = rho_u_new[j][i] * rho_inv;
                v[j][i] = rho_v_new[j][i] * rho_inv;
                E[j][i] = rho_E_new[j][i] * rho_inv;
                I[j][i] = E[j][i] - 0.5*(u[j][i]*u[j][i] + v[j][i]*v[j][i]);
                P[j][i] = (gamma - 1.0) * rho_new[j][i] * I[j][i];
                rho[j][i] = rho_new[j][i];
            }

        for (int j = 0; j < Ny; ++j)
            applyBoundaryConditions(rho[j].data(), u[j].data(), P[j].data(),
                                    cells_x, ghost_cells, left_bc, right_bc);
        for (int i = 0; i < Nx; ++i) {
            std::vector<double> rho_col(Ny), v_col(Ny), P_col(Ny);
            for (int j = 0; j < Ny; ++j) {
                rho_col[j] = rho[j][i];
                v_col[j]   = v[j][i];
                P_col[j]   = P[j][i];
            }
            applyBoundaryConditions(rho_col.data(), v_col.data(), P_col.data(),
                                    cells_y, ghost_cells, bottom_bc, top_bc);
            for (int j = 0; j < Ny; ++j) {
                rho[j][i] = rho_col[j];
                v[j][i]   = v_col[j];
                P[j][i]   = P_col[j];
            }
        }
        for (int j = 0; j < Ny; ++j){
        for (int i = 0; i < Nx; ++i){
            E[j][i] = I[j][i] + 0.5 * (u[j][i]*u[j][i] + v[j][i]*v[j][i]);}}

        std::vector<std::vector<double>> u_tildeX(Ny, std::vector<double>(Nx));
        std::vector<std::vector<double>> v_tildeX(Ny, std::vector<double>(Nx));
        std::vector<std::vector<double>> E_tildeX(Ny, std::vector<double>(Nx));

        for (int j = ghost_cells; j < Ny - ghost_cells; ++j)
            for (int i = ghost_cells; i < Nx - ghost_cells; ++i) {
                double pL = 0.5 * (P[j][i] + P[j][i-1]);
                double pR = 0.5 * (P[j][i] + P[j][i+1]);
                double uL = 0.5 * (u[j][i] + u[j][i-1]);
                double uR = 0.5 * (u[j][i] + u[j][i+1]);

                u_tildeX[j][i] = u[j][i] - (tau / rho[j][i]) * (pR - pL) / hx;
                v_tildeX[j][i] = v[j][i];
                E_tildeX[j][i] = E[j][i] - (tau / rho[j][i]) * (pR * uR - pL * uL) / hx;
            }

        for (int j = 0; j < Ny; ++j) {
            applyBoundaryConditions(rho[j].data(), u_tildeX[j].data(), P[j].data(),
                                    cells_x, ghost_cells, left_bc, right_bc);
        }
        for (int i = 0; i < Nx; ++i) {
            std::vector<double> rho_col(Ny), u_col(Ny), v_col(Ny), E_col(Ny), P_col(Ny);
            for (int j = 0; j < Ny; ++j) {
                rho_col[j] = rho[j][i];
                u_col[j]   = u_tildeX[j][i];
                v_col[j]   = v_tildeX[j][i];
                E_col[j]   = E_tildeX[j][i];
                P_col[j]   = P[j][i];
            }
            applyBoundaryConditions(rho_col.data(), u_col.data(), P_col.data(),
                                    cells_y, ghost_cells, bottom_bc, top_bc);
            applyBoundaryConditions(rho_col.data(), v_col.data(), P_col.data(),
                                    cells_y, ghost_cells, bottom_bc, top_bc);
            applyBoundaryConditions(rho_col.data(), E_col.data(), P_col.data(),
                                    cells_y, ghost_cells, bottom_bc, top_bc);
            for (int j = 0; j < Ny; ++j) {
                u_tildeX[j][i] = u_col[j];
                v_tildeX[j][i] = v_col[j];
                E_tildeX[j][i] = E_col[j];
            }
        }

        for (int j = 0; j < Ny; ++j)
            for (int i = 0; i < Nx; ++i) {
                rho_u_tilde[j][i] = rho[j][i] * u_tildeX[j][i];
                rho_v_tilde[j][i] = rho[j][i] * v_tildeX[j][i];
                rho_E_tilde[j][i] = rho[j][i] * E_tildeX[j][i];
            }

        for (int j = 0; j < Ny; ++j)
            for (int i = 0; i < Nx; ++i) {
                rho_new[j][i]    = rho[j][i];
                rho_u_new[j][i]  = rho_u_tilde[j][i];
                rho_v_new[j][i]  = rho_v_tilde[j][i];
                rho_E_new[j][i]  = rho_E_tilde[j][i];
            }

        for (int j = ghost_cells; j < Ny - ghost_cells; ++j)
            for (int i = ghost_cells; i < Nx - ghost_cells + 1; ++i) {
                double u_face = 0.5 * (u_tildeX[j][i-1] + u_tildeX[j][i]);
                int donor_i = (u_face >= 0) ? i-1 : i;

                double flux_mass   = rho[j][donor_i] * u_face;
                double flux_rho_u  = rho_u_tilde[j][donor_i] * u_face;
                double flux_rho_v  = rho_v_tilde[j][donor_i] * u_face;
                double flux_rho_E  = rho_E_tilde[j][donor_i] * u_face;

                if (i >= ghost_cells && i < Nx - ghost_cells) {
                    rho_new[j][i]    += (tau / hx) * flux_mass;
                    rho_u_new[j][i]  += (tau / hx) * flux_rho_u;
                    rho_v_new[j][i]  += (tau / hx) * flux_rho_v;
                    rho_E_new[j][i]  += (tau / hx) * flux_rho_E;
                }
                if (i-1 >= ghost_cells && i-1 < Nx - ghost_cells) {
                    rho_new[j][i-1]    -= (tau / hx) * flux_mass;
                    rho_u_new[j][i-1]  -= (tau / hx) * flux_rho_u;
                    rho_v_new[j][i-1]  -= (tau / hx) * flux_rho_v;
                    rho_E_new[j][i-1]  -= (tau / hx) * flux_rho_E;
                }
            }

        for (int j = ghost_cells; j < Ny - ghost_cells; ++j)
            for (int i = ghost_cells; i < Nx - ghost_cells; ++i) {
                double rho_inv = 1.0 / rho_new[j][i];
                u[j][i] = rho_u_new[j][i] * rho_inv;
                v[j][i] = rho_v_new[j][i] * rho_inv;
                E[j][i] = rho_E_new[j][i] * rho_inv;
                I[j][i] = E[j][i] - 0.5*(u[j][i]*u[j][i] + v[j][i]*v[j][i]);
                P[j][i] = (gamma - 1.0) * rho_new[j][i] * I[j][i];
                rho[j][i] = rho_new[j][i];
            }
    }

    for (int j = 0; j < Ny; ++j) {
        applyBoundaryConditions(rho[j].data(), u[j].data(), P[j].data(),
                                cells_x, ghost_cells, left_bc, right_bc);
    }
    for (int i = 0; i < Nx; ++i) {
        std::vector<double> rho_col(Ny), v_col(Ny), P_col(Ny);
        for (int j = 0; j < Ny; ++j) {
            rho_col[j] = rho[j][i];
            v_col[j]   = v[j][i];
            P_col[j]   = P[j][i];
        }
        applyBoundaryConditions(rho_col.data(), v_col.data(), P_col.data(),
                                cells_y, ghost_cells, bottom_bc, top_bc);
        for (int j = 0; j < Ny; ++j) {
            rho[j][i] = rho_col[j];
            v[j][i]   = v_col[j];
            P[j][i]   = P_col[j];
        }
    }

    t_total += tau;
}

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
    const std::string& nameRS)
{
    bool order_xy = true; 

    const double gamma = 1.4;
    int Nx = cells_x + 2 * ghost_cells;
    int Ny = cells_y + 2 * ghost_cells;

    std::vector<std::vector<double>> E(Ny, std::vector<double>(Nx));
    for (int j = 0; j < Ny; ++j){
        for (int i = 0; i < Nx; ++i){
            E[j][i] = I[j][i] + 0.5 * (u[j][i]*u[j][i] + v[j][i]*v[j][i]);}}

    for (int j = 0; j < Ny; ++j) {
        applyBoundaryConditions(rho[j].data(), u[j].data(), P[j].data(),
                                cells_x, ghost_cells, left_bc, right_bc);
    }
    for (int i = 0; i < Nx; ++i) {
        std::vector<double> rho_col(Ny), v_col(Ny), P_col(Ny);
        for (int j = 0; j < Ny; ++j) {
            rho_col[j] = rho[j][i];
            v_col[j]   = v[j][i];
            P_col[j]   = P[j][i];
        }
        applyBoundaryConditions(rho_col.data(), v_col.data(), P_col.data(),
                                cells_y, ghost_cells, bottom_bc, top_bc);
        for (int j = 0; j < Ny; ++j) {
            rho[j][i] = rho_col[j];
            v[j][i]   = v_col[j];
            P[j][i]   = P_col[j];
        }
    }

    if (order_xy) {
        std::vector<std::vector<double>> u_tilde(Ny, std::vector<double>(Nx));
        std::vector<std::vector<double>> v_tilde(Ny, std::vector<double>(Nx));
        std::vector<std::vector<double>> E_tilde(Ny, std::vector<double>(Nx));

        for (int j = ghost_cells; j < Ny - ghost_cells; ++j) {
            for (int i = ghost_cells; i < Nx - ghost_cells; ++i) {
                double pL = 0.5 * (P[j][i] + P[j][i-1]);
                double pR = 0.5 * (P[j][i] + P[j][i+1]);
                double uL = 0.5 * (u[j][i] + u[j][i-1]);
                double uR = 0.5 * (u[j][i] + u[j][i+1]);

                u_tilde[j][i] = u[j][i] - (tau / rho[j][i]) * (pR - pL) / hx;
                v_tilde[j][i] = v[j][i];
                E_tilde[j][i] = E[j][i] - (tau / rho[j][i]) * (pR * uR - pL * uL) / hx;
            }
        }

        for (int j = 0; j < Ny; ++j) {
            applyBoundaryConditions(rho[j].data(), u_tilde[j].data(), P[j].data(),
                                    cells_x, ghost_cells, left_bc, right_bc);
        }
        for (int i = 0; i < Nx; ++i) {
            std::vector<double> rho_col(Ny), u_col(Ny), v_col(Ny), E_col(Ny), P_col(Ny);
            for (int j = 0; j < Ny; ++j) {
                rho_col[j] = rho[j][i];
                u_col[j]   = u_tilde[j][i];
                v_col[j]   = v_tilde[j][i];
                E_col[j]   = E_tilde[j][i];
                P_col[j]   = P[j][i];
            }
            applyBoundaryConditions(rho_col.data(), u_col.data(), P_col.data(),
                                    cells_y, ghost_cells, bottom_bc, top_bc);
            applyBoundaryConditions(rho_col.data(), v_col.data(), P_col.data(),
                                    cells_y, ghost_cells, bottom_bc, top_bc);
            applyBoundaryConditions(rho_col.data(), E_col.data(), P_col.data(),
                                    cells_y, ghost_cells, bottom_bc, top_bc);
            for (int j = 0; j < Ny; ++j) {
                u_tilde[j][i] = u_col[j];
                v_tilde[j][i] = v_col[j];
                E_tilde[j][i] = E_col[j];
            }
        }

        std::vector<std::vector<double>> rho_u_tilde(Ny, std::vector<double>(Nx));
        std::vector<std::vector<double>> rho_v_tilde(Ny, std::vector<double>(Nx));
        std::vector<std::vector<double>> rho_E_tilde(Ny, std::vector<double>(Nx));
        for (int j = 0; j < Ny; ++j)
            for (int i = 0; i < Nx; ++i) {
                rho_u_tilde[j][i] = rho[j][i] * u_tilde[j][i];
                rho_v_tilde[j][i] = rho[j][i] * v_tilde[j][i];
                rho_E_tilde[j][i] = rho[j][i] * E_tilde[j][i];
            }

        std::vector<std::vector<double>> rho_new(Ny, std::vector<double>(Nx));
        std::vector<std::vector<double>> rho_u_new(Ny, std::vector<double>(Nx));
        std::vector<std::vector<double>> rho_v_new(Ny, std::vector<double>(Nx));
        std::vector<std::vector<double>> rho_E_new(Ny, std::vector<double>(Nx));
        for (int j = 0; j < Ny; ++j)
            for (int i = 0; i < Nx; ++i) {
                rho_new[j][i]    = rho[j][i];
                rho_u_new[j][i]  = rho_u_tilde[j][i];
                rho_v_new[j][i]  = rho_v_tilde[j][i];
                rho_E_new[j][i]  = rho_E_tilde[j][i];
            }

        for (int j = ghost_cells; j < Ny - ghost_cells; ++j) {
            for (int i = ghost_cells; i < Nx - ghost_cells + 1; ++i) {
                double u_face = 0.5 * (u_tilde[j][i-1] + u_tilde[j][i]);
                int donor_i = (u_face >= 0) ? i-1 : i;

                double flux_mass   = rho[j][donor_i] * u_face;
                double flux_rho_u  = rho_u_tilde[j][donor_i] * u_face;
                double flux_rho_v  = rho_v_tilde[j][donor_i] * u_face;
                double flux_rho_E  = rho_E_tilde[j][donor_i] * u_face;

                if (i >= ghost_cells && i < Nx - ghost_cells) {
                    rho_new[j][i]    += (tau / hx) * flux_mass;
                    rho_u_new[j][i]  += (tau / hx) * flux_rho_u;
                    rho_v_new[j][i]  += (tau / hx) * flux_rho_v;
                    rho_E_new[j][i]  += (tau / hx) * flux_rho_E;
                }
                if (i-1 >= ghost_cells && i-1 < Nx - ghost_cells) {
                    rho_new[j][i-1]    -= (tau / hx) * flux_mass;
                    rho_u_new[j][i-1]  -= (tau / hx) * flux_rho_u;
                    rho_v_new[j][i-1]  -= (tau / hx) * flux_rho_v;
                    rho_E_new[j][i-1]  -= (tau / hx) * flux_rho_E;
                }
            }
        }

        for (int j = ghost_cells; j < Ny - ghost_cells; ++j)
            for (int i = ghost_cells; i < Nx - ghost_cells; ++i) {
                double rho_inv = 1.0 / rho_new[j][i];
                u[j][i] = rho_u_new[j][i] * rho_inv;
                v[j][i] = rho_v_new[j][i] * rho_inv;
                E[j][i] = rho_E_new[j][i] * rho_inv;
                I[j][i] = E[j][i] - 0.5*(u[j][i]*u[j][i] + v[j][i]*v[j][i]);
                P[j][i] = (gamma - 1.0) * rho_new[j][i] * I[j][i];
                rho[j][i] = rho_new[j][i];
            }

        for (int j = 0; j < Ny; ++j)
            applyBoundaryConditions(rho[j].data(), u[j].data(), P[j].data(),
                                    cells_x, ghost_cells, left_bc, right_bc);
        for (int i = 0; i < Nx; ++i) {
            std::vector<double> rho_col(Ny), v_col(Ny), P_col(Ny);
            for (int j = 0; j < Ny; ++j) {
                rho_col[j] = rho[j][i];
                v_col[j]   = v[j][i];
                P_col[j]   = P[j][i];
            }
            applyBoundaryConditions(rho_col.data(), v_col.data(), P_col.data(),
                                    cells_y, ghost_cells, bottom_bc, top_bc);
            for (int j = 0; j < Ny; ++j) {
                rho[j][i] = rho_col[j];
                v[j][i]   = v_col[j];
                P[j][i]   = P_col[j];
            }
        }
    for (int j = 0; j < Ny; ++j){
        for (int i = 0; i < Nx; ++i){
            E[j][i] = I[j][i] + 0.5 * (u[j][i]*u[j][i] + v[j][i]*v[j][i]);}}

        std::vector<std::vector<double>> u_tilde2(Ny, std::vector<double>(Nx));
        std::vector<std::vector<double>> v_tilde2(Ny, std::vector<double>(Nx));
        std::vector<std::vector<double>> E_tilde2(Ny, std::vector<double>(Nx));

        for (int j = ghost_cells; j < Ny - ghost_cells; ++j)
            for (int i = ghost_cells; i < Nx - ghost_cells; ++i) {
                double pB = 0.5 * (P[j][i] + P[j-1][i]);
                double pT = 0.5 * (P[j][i] + P[j+1][i]);
                double vB = 0.5 * (v[j][i] + v[j-1][i]);
                double vT = 0.5 * (v[j][i] + v[j+1][i]);

                u_tilde2[j][i] = u[j][i];
                v_tilde2[j][i] = v[j][i] - (tau / rho[j][i]) * (pT - pB) / hy;
                E_tilde2[j][i] = E[j][i] - (tau / rho[j][i]) * (pT * vT - pB * vB) / hy;
            }

        for (int j = 0; j < Ny; ++j) {
            applyBoundaryConditions(rho[j].data(), u_tilde2[j].data(), P[j].data(),
                                    cells_x, ghost_cells, left_bc, right_bc);
        }
        for (int i = 0; i < Nx; ++i) {
            std::vector<double> rho_col(Ny), u_col(Ny), v_col(Ny), E_col(Ny), P_col(Ny);
            for (int j = 0; j < Ny; ++j) {
                rho_col[j] = rho[j][i];
                u_col[j]   = u_tilde2[j][i];
                v_col[j]   = v_tilde2[j][i];
                E_col[j]   = E_tilde2[j][i];
                P_col[j]   = P[j][i];
            }
            applyBoundaryConditions(rho_col.data(), u_col.data(), P_col.data(),
                                    cells_y, ghost_cells, bottom_bc, top_bc);
            applyBoundaryConditions(rho_col.data(), v_col.data(), P_col.data(),
                                    cells_y, ghost_cells, bottom_bc, top_bc);
            applyBoundaryConditions(rho_col.data(), E_col.data(), P_col.data(),
                                    cells_y, ghost_cells, bottom_bc, top_bc);
            for (int j = 0; j < Ny; ++j) {
                u_tilde2[j][i] = u_col[j];
                v_tilde2[j][i] = v_col[j];
                E_tilde2[j][i] = E_col[j];
            }
        }

        for (int j = 0; j < Ny; ++j)
            for (int i = 0; i < Nx; ++i) {
                rho_u_tilde[j][i] = rho[j][i] * u_tilde2[j][i];
                rho_v_tilde[j][i] = rho[j][i] * v_tilde2[j][i];
                rho_E_tilde[j][i] = rho[j][i] * E_tilde2[j][i];
            }

        for (int j = 0; j < Ny; ++j)
            for (int i = 0; i < Nx; ++i) {
                rho_new[j][i]    = rho[j][i];
                rho_u_new[j][i]  = rho_u_tilde[j][i];
                rho_v_new[j][i]  = rho_v_tilde[j][i];
                rho_E_new[j][i]  = rho_E_tilde[j][i];
            }

        for (int j = ghost_cells; j < Ny - ghost_cells + 1; ++j)
            for (int i = ghost_cells; i < Nx - ghost_cells; ++i) {
                double v_face = 0.5 * (v_tilde2[j-1][i] + v_tilde2[j][i]);
                int donor_j = (v_face >= 0) ? j-1 : j;

                double flux_mass   = rho[donor_j][i] * v_face;
                double flux_rho_u  = rho_u_tilde[donor_j][i] * v_face;
                double flux_rho_v  = rho_v_tilde[donor_j][i] * v_face;
                double flux_rho_E  = rho_E_tilde[donor_j][i] * v_face;

                if (j >= ghost_cells && j < Ny - ghost_cells) {
                    rho_new[j][i]    += (tau / hy) * flux_mass;
                    rho_u_new[j][i]  += (tau / hy) * flux_rho_u;
                    rho_v_new[j][i]  += (tau / hy) * flux_rho_v;
                    rho_E_new[j][i]  += (tau / hy) * flux_rho_E;
                }
                if (j-1 >= ghost_cells && j-1 < Ny - ghost_cells) {
                    rho_new[j-1][i]    -= (tau / hy) * flux_mass;
                    rho_u_new[j-1][i]  -= (tau / hy) * flux_rho_u;
                    rho_v_new[j-1][i]  -= (tau / hy) * flux_rho_v;
                    rho_E_new[j-1][i]  -= (tau / hy) * flux_rho_E;
                }
            }

        for (int j = ghost_cells; j < Ny - ghost_cells; ++j)
            for (int i = ghost_cells; i < Nx - ghost_cells; ++i) {
                double rho_inv = 1.0 / rho_new[j][i];
                u[j][i] = rho_u_new[j][i] * rho_inv;
                v[j][i] = rho_v_new[j][i] * rho_inv;
                E[j][i] = rho_E_new[j][i] * rho_inv;
                I[j][i] = E[j][i] - 0.5*(u[j][i]*u[j][i] + v[j][i]*v[j][i]);
                P[j][i] = (gamma - 1.0) * rho_new[j][i] * I[j][i];
                rho[j][i] = rho_new[j][i];
            }
    } else {
        std::vector<std::vector<double>> u_tildeY(Ny, std::vector<double>(Nx));
        std::vector<std::vector<double>> v_tildeY(Ny, std::vector<double>(Nx));
        std::vector<std::vector<double>> E_tildeY(Ny, std::vector<double>(Nx));
        std::vector<std::vector<double>> rho_u_tilde(Ny, std::vector<double>(Nx));
        std::vector<std::vector<double>> rho_v_tilde(Ny, std::vector<double>(Nx));
        std::vector<std::vector<double>> rho_E_tilde(Ny, std::vector<double>(Nx));
        std::vector<std::vector<double>> rho_new(Ny, std::vector<double>(Nx));
        std::vector<std::vector<double>> rho_u_new(Ny, std::vector<double>(Nx));
        std::vector<std::vector<double>> rho_v_new(Ny, std::vector<double>(Nx));
        std::vector<std::vector<double>> rho_E_new(Ny, std::vector<double>(Nx));

        for (int j = ghost_cells; j < Ny - ghost_cells; ++j)
            for (int i = ghost_cells; i < Nx - ghost_cells; ++i) {
                double pB = 0.5 * (P[j][i] + P[j-1][i]);
                double pT = 0.5 * (P[j][i] + P[j+1][i]);
                double vB = 0.5 * (v[j][i] + v[j-1][i]);
                double vT = 0.5 * (v[j][i] + v[j+1][i]);

                u_tildeY[j][i] = u[j][i];
                v_tildeY[j][i] = v[j][i] - (tau / rho[j][i]) * (pT - pB) / hy;
                E_tildeY[j][i] = E[j][i] - (tau / rho[j][i]) * (pT * vT - pB * vB) / hy;
            }

        for (int j = 0; j < Ny; ++j) {
            applyBoundaryConditions(rho[j].data(), u_tildeY[j].data(), P[j].data(),
                                    cells_x, ghost_cells, left_bc, right_bc);
         }
        for (int i = 0; i < Nx; ++i) {
            std::vector<double> rho_col(Ny), u_col(Ny), v_col(Ny), E_col(Ny), P_col(Ny);
            for (int j = 0; j < Ny; ++j) {
                rho_col[j] = rho[j][i];
                u_col[j]   = u_tildeY[j][i];
                v_col[j]   = v_tildeY[j][i];
                E_col[j]   = E_tildeY[j][i];
                P_col[j]   = P[j][i];
            }
            applyBoundaryConditions(rho_col.data(), u_col.data(), P_col.data(),
                                    cells_y, ghost_cells, bottom_bc, top_bc);
            applyBoundaryConditions(rho_col.data(), v_col.data(), P_col.data(),
                                    cells_y, ghost_cells, bottom_bc, top_bc);
            applyBoundaryConditions(rho_col.data(), E_col.data(), P_col.data(),
                                    cells_y, ghost_cells, bottom_bc, top_bc);
            for (int j = 0; j < Ny; ++j) {
                u_tildeY[j][i] = u_col[j];
                v_tildeY[j][i] = v_col[j];
                E_tildeY[j][i] = E_col[j];
            }
        }

        for (int j = 0; j < Ny; ++j)
            for (int i = 0; i < Nx; ++i) {
                rho_u_tilde[j][i] = rho[j][i] * u_tildeY[j][i];
                rho_v_tilde[j][i] = rho[j][i] * v_tildeY[j][i];
                rho_E_tilde[j][i] = rho[j][i] * E_tildeY[j][i];
            }

        for (int j = 0; j < Ny; ++j)
            for (int i = 0; i < Nx; ++i) {
                rho_new[j][i]    = rho[j][i];
                rho_u_new[j][i]  = rho_u_tilde[j][i];
                rho_v_new[j][i]  = rho_v_tilde[j][i];
                rho_E_new[j][i]  = rho_E_tilde[j][i];
            }

        for (int j = ghost_cells; j < Ny - ghost_cells + 1; ++j)
            for (int i = ghost_cells; i < Nx - ghost_cells; ++i) {
                double v_face = 0.5 * (v_tildeY[j-1][i] + v_tildeY[j][i]);
                int donor_j = (v_face >= 0) ? j-1 : j;

                double flux_mass   = rho[donor_j][i] * v_face;
                double flux_rho_u  = rho_u_tilde[donor_j][i] * v_face;
                double flux_rho_v  = rho_v_tilde[donor_j][i] * v_face;
                double flux_rho_E  = rho_E_tilde[donor_j][i] * v_face;

                if (j >= ghost_cells && j < Ny - ghost_cells) {
                    rho_new[j][i]    += (tau / hy) * flux_mass;
                    rho_u_new[j][i]  += (tau / hy) * flux_rho_u;
                    rho_v_new[j][i]  += (tau / hy) * flux_rho_v;
                    rho_E_new[j][i]  += (tau / hy) * flux_rho_E;
                }
                if (j-1 >= ghost_cells && j-1 < Ny - ghost_cells) {
                    rho_new[j-1][i]    -= (tau / hy) * flux_mass;
                    rho_u_new[j-1][i]  -= (tau / hy) * flux_rho_u;
                    rho_v_new[j-1][i]  -= (tau / hy) * flux_rho_v;
                    rho_E_new[j-1][i]  -= (tau / hy) * flux_rho_E;
                }
            }

        for (int j = ghost_cells; j < Ny - ghost_cells; ++j)
            for (int i = ghost_cells; i < Nx - ghost_cells; ++i) {
                double rho_inv = 1.0 / rho_new[j][i];
                u[j][i] = rho_u_new[j][i] * rho_inv;
                v[j][i] = rho_v_new[j][i] * rho_inv;
                E[j][i] = rho_E_new[j][i] * rho_inv;
                I[j][i] = E[j][i] - 0.5*(u[j][i]*u[j][i] + v[j][i]*v[j][i]);
                P[j][i] = (gamma - 1.0) * rho_new[j][i] * I[j][i];
                rho[j][i] = rho_new[j][i];
            }

        for (int j = 0; j < Ny; ++j)
            applyBoundaryConditions(rho[j].data(), u[j].data(), P[j].data(),
                                    cells_x, ghost_cells, left_bc, right_bc);
        for (int i = 0; i < Nx; ++i) {
            std::vector<double> rho_col(Ny), v_col(Ny), P_col(Ny);
            for (int j = 0; j < Ny; ++j) {
                rho_col[j] = rho[j][i];
                v_col[j]   = v[j][i];
                P_col[j]   = P[j][i];
            }
            applyBoundaryConditions(rho_col.data(), v_col.data(), P_col.data(),
                                    cells_y, ghost_cells, bottom_bc, top_bc);
            for (int j = 0; j < Ny; ++j) {
                rho[j][i] = rho_col[j];
                v[j][i]   = v_col[j];
                P[j][i]   = P_col[j];
            }
        }
    for (int j = 0; j < Ny; ++j){
        for (int i = 0; i < Nx; ++i){
            E[j][i] = I[j][i] + 0.5 * (u[j][i]*u[j][i] + v[j][i]*v[j][i]);}}

        std::vector<std::vector<double>> u_tildeX(Ny, std::vector<double>(Nx));
        std::vector<std::vector<double>> v_tildeX(Ny, std::vector<double>(Nx));
        std::vector<std::vector<double>> E_tildeX(Ny, std::vector<double>(Nx));

        for (int j = ghost_cells; j < Ny - ghost_cells; ++j)
            for (int i = ghost_cells; i < Nx - ghost_cells; ++i) {
                double pL = 0.5 * (P[j][i] + P[j][i-1]);
                double pR = 0.5 * (P[j][i] + P[j][i+1]);
                double uL = 0.5 * (u[j][i] + u[j][i-1]);
                double uR = 0.5 * (u[j][i] + u[j][i+1]);

                u_tildeX[j][i] = u[j][i] - (tau / rho[j][i]) * (pR - pL) / hx;
                v_tildeX[j][i] = v[j][i];
                E_tildeX[j][i] = E[j][i] - (tau / rho[j][i]) * (pR * uR - pL * uL) / hx;
            }

        for (int j = 0; j < Ny; ++j) {
            applyBoundaryConditions(rho[j].data(), u_tildeX[j].data(), P[j].data(),
                                    cells_x, ghost_cells, left_bc, right_bc);
        }
        for (int i = 0; i < Nx; ++i) {
            std::vector<double> rho_col(Ny), u_col(Ny), v_col(Ny), E_col(Ny), P_col(Ny);
            for (int j = 0; j < Ny; ++j) {
                rho_col[j] = rho[j][i];
                u_col[j]   = u_tildeX[j][i];
                v_col[j]   = v_tildeX[j][i];
                E_col[j]   = E_tildeX[j][i];
                P_col[j]   = P[j][i];
            }
            applyBoundaryConditions(rho_col.data(), u_col.data(), P_col.data(),
                                    cells_y, ghost_cells, bottom_bc, top_bc);
            applyBoundaryConditions(rho_col.data(), v_col.data(), P_col.data(),
                                    cells_y, ghost_cells, bottom_bc, top_bc);
            applyBoundaryConditions(rho_col.data(), E_col.data(), P_col.data(),
                                    cells_y, ghost_cells, bottom_bc, top_bc);
            for (int j = 0; j < Ny; ++j) {
                u_tildeX[j][i] = u_col[j];
                v_tildeX[j][i] = v_col[j];
                E_tildeX[j][i] = E_col[j];
            }
        }

        for (int j = 0; j < Ny; ++j)
            for (int i = 0; i < Nx; ++i) {
                rho_u_tilde[j][i] = rho[j][i] * u_tildeX[j][i];
                rho_v_tilde[j][i] = rho[j][i] * v_tildeX[j][i];
                rho_E_tilde[j][i] = rho[j][i] * E_tildeX[j][i];
            }

        for (int j = 0; j < Ny; ++j)
            for (int i = 0; i < Nx; ++i) {
                rho_new[j][i]    = rho[j][i];
                rho_u_new[j][i]  = rho_u_tilde[j][i];
                rho_v_new[j][i]  = rho_v_tilde[j][i];
                rho_E_new[j][i]  = rho_E_tilde[j][i];
            }

        for (int j = ghost_cells; j < Ny - ghost_cells; ++j)
            for (int i = ghost_cells; i < Nx - ghost_cells + 1; ++i) {
                double u_face = 0.5 * (u_tildeX[j][i-1] + u_tildeX[j][i]);
                int donor_i = (u_face >= 0) ? i-1 : i;

                double flux_mass   = rho[j][donor_i] * u_face;
                double flux_rho_u  = rho_u_tilde[j][donor_i] * u_face;
                double flux_rho_v  = rho_v_tilde[j][donor_i] * u_face;
                double flux_rho_E  = rho_E_tilde[j][donor_i] * u_face;

                if (i >= ghost_cells && i < Nx - ghost_cells) {
                    rho_new[j][i]    += (tau / hx) * flux_mass;
                    rho_u_new[j][i]  += (tau / hx) * flux_rho_u;
                    rho_v_new[j][i]  += (tau / hx) * flux_rho_v;
                    rho_E_new[j][i]  += (tau / hx) * flux_rho_E;
                }
                if (i-1 >= ghost_cells && i-1 < Nx - ghost_cells) {
                    rho_new[j][i-1]    -= (tau / hx) * flux_mass;
                    rho_u_new[j][i-1]  -= (tau / hx) * flux_rho_u;
                    rho_v_new[j][i-1]  -= (tau / hx) * flux_rho_v;
                    rho_E_new[j][i-1]  -= (tau / hx) * flux_rho_E;
                }
            }

        for (int j = ghost_cells; j < Ny - ghost_cells; ++j)
            for (int i = ghost_cells; i < Nx - ghost_cells; ++i) {
                double rho_inv = 1.0 / rho_new[j][i];
                u[j][i] = rho_u_new[j][i] * rho_inv;
                v[j][i] = rho_v_new[j][i] * rho_inv;
                E[j][i] = rho_E_new[j][i] * rho_inv;
                I[j][i] = E[j][i] - 0.5*(u[j][i]*u[j][i] + v[j][i]*v[j][i]);
                P[j][i] = (gamma - 1.0) * rho_new[j][i] * I[j][i];
                rho[j][i] = rho_new[j][i];
            }
    }

    for (int j = 0; j < Ny; ++j) {
        applyBoundaryConditions(rho[j].data(), u[j].data(), P[j].data(),
                                cells_x, ghost_cells, left_bc, right_bc);
    }
    for (int i = 0; i < Nx; ++i) {
        std::vector<double> rho_col(Ny), v_col(Ny), P_col(Ny);
        for (int j = 0; j < Ny; ++j) {
            rho_col[j] = rho[j][i];
            v_col[j]   = v[j][i];
            P_col[j]   = P[j][i];
        }
        applyBoundaryConditions(rho_col.data(), v_col.data(), P_col.data(),
                                cells_y, ghost_cells, bottom_bc, top_bc);
        for (int j = 0; j < Ny; ++j) {
            rho[j][i] = rho_col[j];
            v[j][i]   = v_col[j];
            P[j][i]   = P_col[j];
        }
    }
   
std::vector<std::vector<bool>> isSolid(Ny, std::vector<bool>(Nx, false));
for (int j = ghost_cells; j < Ny-ghost_cells; ++j) {
    for (int i = ghost_cells; i < Nx-ghost_cells; ++i) {
	
	if((i*hx - 0.4)*(i*hx - 0.4)+(j*hy - 0.5)*(j*hy - 0.5)<0.01 ){
	   /* rho[j][i] = 10.0;
            v[j][i]   = 0;
	    u[j][i]   = 0;
            P[j][i]   = 0.1;*/
	    isSolid[j][i] = true;
    }    
    }
}
for (int j = ghost_cells; j < Ny - ghost_cells; ++j) {
    for (int i = ghost_cells; i < Nx - ghost_cells; ++i) {
        if (!isSolid[j][i]) { 
            if (isSolid[j][i+1]) {
                rho[j][i+1] = rho[j][i];
                P[j][i+1]   = P[j][i];
                u[j][i+1]   = -u[j][i];
                v[j][i+1]   = v[j][i];
            }
            if (isSolid[j][i-1]) {
                rho[j][i-1] = rho[j][i];
                P[j][i-1]   = P[j][i];
                u[j][i-1]   = -u[j][i];
                v[j][i-1]   = v[j][i];
            }
            // Сверху
            if (isSolid[j+1][i]) {
                rho[j+1][i] = rho[j][i];
                P[j+1][i]   = P[j][i];
                u[j+1][i]   = u[j][i];
                v[j+1][i]   = -v[j][i];
            }
            // Снизу
            if (isSolid[j-1][i]) {
                rho[j-1][i] = rho[j][i];
                P[j-1][i]   = P[j][i];
                u[j-1][i]   = u[j][i];
                v[j-1][i]   = -v[j][i];
            }
        }
    }
}

    t_total += tau;   
    std::cout << "t = " << t_total << ", tau = " << tau
              << ", method = FLICCF" << std::endl;

}


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
    double Z_freq, double E_act_over_R, double R_gas)
{
    const double gamma = 1.4;
    const double MINGRHO = 0.0;

    static int step_counter = 0;
    step_counter++;

    int Nx = cells_x + 2 * ghost_cells;
    int Ny = cells_y + 2 * ghost_cells;

    int i_start = ghost_cells;
    int i_end = ghost_cells + cells_x - 1;
    int j_start = ghost_cells;
    int j_end = ghost_cells + cells_y - 1;

    std::vector<std::vector<double>> u_tilde(Ny, std::vector<double>(Nx, 0.0));
    std::vector<std::vector<double>> v_tilde(Ny, std::vector<double>(Nx, 0.0));
    std::vector<std::vector<double>> rho_tilde(Ny, std::vector<double>(Nx, 0.0));
    std::vector<std::vector<double>> I_tilde(Ny, std::vector<double>(Nx, 0.0));
    std::vector<std::vector<double>> E_total(Ny, std::vector<double>(Nx, 0.0));

    std::vector<std::vector<double>> q1(Ny, std::vector<double>(Nx, 0.0));
    std::vector<std::vector<double>> q2(Ny, std::vector<double>(Nx, 0.0));
    std::vector<std::vector<double>> q3(Ny, std::vector<double>(Nx, 0.0));
    std::vector<std::vector<double>> q4(Ny, std::vector<double>(Nx, 0.0));

    std::vector<std::vector<double>> DM(Ny, std::vector<double>(Nx, 0.0));
    std::vector<std::vector<double>> DE(Ny, std::vector<double>(Nx, 0.0));
    std::vector<std::vector<double>> DPU(Ny, std::vector<double>(Nx, 0.0));
    std::vector<std::vector<double>> DPV(Ny, std::vector<double>(Nx, 0.0));
    std::vector<std::vector<double>> DW(Ny, std::vector<double>(Nx, 0.0)); 

double dt=tau;
const double C = 0.25;
const double S = 1.5;
const double gamma_s = 1.2;
const double rho0 = 1.84;
const double V0 = 1.0 / rho0;
const double Q = 0.08;
 R_gas = 1e-5; 

for (int j = j_start; j <= j_end; ++j) {
    for (int i = i_start; i <= i_end; ++i) {
        double V = 1.0 / rho[j][i];
        double dV = V0 - V;
        // Защита от dV < 0 
        if (dV < 0.0) {
            double I_thermal = I[j][i] - (1.0 - W[j][i]) * Q;
            if (I_thermal < 0.0) I_thermal = 0.0;
// CHECK: MADER_EOS
            P[j][i] = (gamma_s / V) * I_thermal; // так как I_H=0
            if (P[j][i] < 0.0) P[j][i] = 0.0;
        } else {
            double denom = V0 - S * dV;
            if (fabs(denom) < 1e-12) denom = 1e-12;
            double P_H = (C * C * dV) / (denom * denom);
            double I_H = 0.5 * P_H * dV;
            double I_thermal = I[j][i] - (1.0 - W[j][i]) * Q;
            if (I_thermal < 0.0) I_thermal = 0.0;
            P[j][i] = P_H + (gamma_s / V) * (I_thermal - I_H);
            if (P[j][i] < 0.0) P[j][i] = 0.0;
        }
        double T = P[j][i] / (rho[j][i] * R_gas);
    }
}  
// CHECK: MADER_ARRHENIUS  
    for (int j = j_start; j <= j_end; ++j) {
        for (int i = i_start; i <= i_end; ++i) {
            double T = P[j][i] / (rho[j][i] * R_gas); 
            if (T > MINWT && W[j][i] > GASW && step_counter > VCNT) {
                double rate = Z_freq * exp(-E_act_over_R / T); // exp(-E*/RT)
                double dW = dt * rate * W[j][i];
                W[j][i] -= dW;
                if (W[j][i] < 0.0) W[j][i] = 0.0;
                if (W[j][i] < GASW) W[j][i] = 0.0; 
            }
        }
    }

    for (int j = 0; j < Ny; ++j) {
        applyBoundaryConditions(rho[j].data(), u[j].data(), P[j].data(),
                                cells_x, ghost_cells, left_bc, right_bc);
    }
// CHECK: MADER_VELOCITY
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


/*for (int j = j_start; j <= j_end; ++j) {
    for (int i = i_start; i <= i_end; ++i) {
        if (v[j][i] >= v[j+1][i]) {
            double dv = v[j][i] - v[j+1][i];
            q1[j][i] = viscosity_coeff * rho[j][i] * dv * dv;  
        } else {
            q1[j][i] = 0.0;
        }

        if (u[j][i] >= u[j][i+1]) {
            double du = u[j][i] - u[j][i+1];
            q2[j][i] = viscosity_coeff * rho[j][i] * du * du;  
        } else {
            q2[j][i] = 0.0;
        }
    }
} */
// CHECK: MADER_VISC
for (int j = j_start; j <= j_end; ++j) {
    for (int i = i_start; i <= i_end; ++i) {
        if (j > j_start && v[j-1][i] > v[j][i]) { 
            double dv = v[j-1][i] - v[j][i];
            q1[j][i] = viscosity_coeff * 0.5 * (rho[j-1][i] + rho[j][i]) * dv * dv;
        } else {
            q1[j][i] = 0.0;
        }

        if (i > i_start && u[j][i-1] > u[j][i]) { 
            double du = u[j][i-1] - u[j][i];
            q2[j][i] = viscosity_coeff * 0.5 * (rho[j][i-1] + rho[j][i]) * du * du;
        } else {
            q2[j][i] = 0.0;
        }
    }
}

for (int j = j_start; j <= j_end; ++j) {
    for (int i = i_start; i <= i_end; ++i) {
        q3[j][i] = q1[j+1][i]; 
        q4[j][i] = q2[j][i+1]; 
    }
}

    for (int j = j_start; j <= j_end; ++j) {
        for (int i = i_start; i <= i_end; ++i) {
            double P1 = P[j-1][i];
            double P2 = P[j][i-1];
            double P3 = P[j+1][i];
            double P4 = P[j][i+1];

            u_tilde[j][i] = u[j][i] - (dt / (rho[j][i] * hx)) * ((P4 + q4[j][i]) - (P2 + q2[j][i]));
            v_tilde[j][i] = v[j][i] - (dt / (rho[j][i] * hy)) * ((P3 + q3[j][i]) - (P1 + q1[j][i]));
        }
    }
// CHECK: MADER_ZIP
    for (int j = j_start; j <= j_end; ++j) {
        for (int i = i_start; i <= i_end; ++i) {
            double u_left  = 0.5 * (u_tilde[j][i-1] + u_tilde[j][i]);
            double u_right = 0.5 * (u_tilde[j][i] + u_tilde[j][i+1]);
            double v_bottom = 0.5 * (v_tilde[j-1][i] + v_tilde[j][i]);
            double v_top    = 0.5 * (v_tilde[j][i] + v_tilde[j+1][i]);

            double divU = (u_right - u_left) / hx + (v_top - v_bottom) / hy;
            rho_tilde[j][i] = rho[j][i] - rho[j][i] * dt * divU;
            if (rho_tilde[j][i] < MINGRHO) rho_tilde[j][i] = MINGRHO;

            double work = 0.0;
            work += P[j][i] * ((u_right - u_left)/hx + (v_top - v_bottom)/hy) * dt;
	    work += q4[j][i] * (u_tilde[j][i+1] - u_tilde[j][i]) / hx * dt;
	    work += q2[j][i] * (u_tilde[j][i] - u_tilde[j][i-1]) / hx * dt;
	    work += q3[j][i] * (v_tilde[j+1][i] - v_tilde[j][i]) / hy * dt;
	    work += q1[j][i] * (v_tilde[j][i] - v_tilde[j-1][i]) / hy * dt;
            I_tilde[j][i] = I[j][i] - work / rho[j][i];
            if (I_tilde[j][i] < 0.0) I_tilde[j][i] = 0.0;

            E_total[j][i] = I_tilde[j][i] + 0.5 * (u_tilde[j][i]*u_tilde[j][i] + v_tilde[j][i]*v_tilde[j][i]);
        }
    }

    for (int j = 0; j < Ny; ++j) {
        for (int i = 0; i < Nx; ++i) {
            DM[j][i] = 0.0; DE[j][i] = 0.0; DPU[j][i] = 0.0; DPV[j][i] = 0.0; DW[j][i] = 0.0;
        }
    }
// CHECK: MADER_DONOR
    for (int j = j_start; j <= j_end-1; ++j) {
        for (int i = i_start; i <= i_end; ++i) {
            double V_avg = 0.5 * (v_tilde[j][i] + v_tilde[j+1][i]);
            double beta = (V_avg * dt / hy) / (1.0 + (v_tilde[j][i] - v_tilde[j+1][i]) * dt / hy + 1e-12);
            beta  = std::max(-1.0, std::min(1.0, beta));
	    int donor_j, donor_i, acc_j, acc_i;
            double abs_beta;
            if (beta >= 0) {
                donor_j = j; donor_i = i; acc_j = j+1; acc_i = i; abs_beta = beta;
            } else {
                donor_j = j+1; donor_i = i; acc_j = j; acc_i = i; abs_beta = -beta;
            }
            double DMASS = rho_tilde[donor_j][donor_i] * abs_beta;

            DM[donor_j][donor_i] -= DMASS;
            DM[acc_j][acc_i] += DMASS;
            DE[donor_j][donor_i] -= DMASS * E_total[donor_j][donor_i];
            DE[acc_j][acc_i] += DMASS * E_total[donor_j][donor_i];
            DPU[donor_j][donor_i] -= DMASS * u_tilde[donor_j][donor_i];
            DPU[acc_j][acc_i] += DMASS * u_tilde[donor_j][donor_i];
            DPV[donor_j][donor_i] -= DMASS * v_tilde[donor_j][donor_i];
            DPV[acc_j][acc_i] += DMASS * v_tilde[donor_j][donor_i];
            // Перенос массовой доли W (используем текущее значение W из основного массива)
            double W_donor = W[donor_j][donor_i];
            DW[donor_j][donor_i] -= DMASS * W_donor;
            DW[acc_j][acc_i] += DMASS * W_donor;
        }
    }

    for (int j = j_start; j <= j_end; ++j) {
        for (int i = i_start; i <= i_end-1; ++i) {
            double U_avg = 0.5 * (u_tilde[j][i] + u_tilde[j][i+1]);
            double alpha = (U_avg * dt / hx) / (1.0 + (u_tilde[j][i] - u_tilde[j][i+1]) * dt / hx + 1e-12);
            alpha = std::max(-1.0, std::min(1.0, alpha));
	    int donor_j, donor_i, acc_j, acc_i;
            double abs_alpha;
            if (alpha >= 0) {
                donor_j = j; donor_i = i; acc_j = j; acc_i = i+1; abs_alpha = alpha;
            } else {
                donor_j = j; donor_i = i+1; acc_j = j; acc_i = i; abs_alpha = -alpha;
            }
            double DMASS = rho_tilde[donor_j][donor_i] * abs_alpha;

            DM[donor_j][donor_i] -= DMASS;
            DM[acc_j][acc_i] += DMASS;
            DE[donor_j][donor_i] -= DMASS * E_total[donor_j][donor_i];
            DE[acc_j][acc_i] += DMASS * E_total[donor_j][donor_i];
            DPU[donor_j][donor_i] -= DMASS * u_tilde[donor_j][donor_i];
            DPU[acc_j][acc_i] += DMASS * u_tilde[donor_j][donor_i];
            DPV[donor_j][donor_i] -= DMASS * v_tilde[donor_j][donor_i];
            DPV[acc_j][acc_i] += DMASS * v_tilde[donor_j][donor_i];
            double W_donor = W[donor_j][donor_i];
            DW[donor_j][donor_i] -= DMASS * W_donor;
            DW[acc_j][acc_i] += DMASS * W_donor;
        }
    }
// CHECK: MADER_REPARTITION
    for (int j = j_start; j <= j_end; ++j) {
        for (int i = i_start; i <= i_end; ++i) {
            double rho_new = rho_tilde[j][i] + DM[j][i];
            if (rho_new < MINGRHO) rho_new = MINGRHO;

            double u_new = (rho_tilde[j][i] * u_tilde[j][i] + DPU[j][i]) / rho_new;
            double v_new = (rho_tilde[j][i] * v_tilde[j][i] + DPV[j][i]) / rho_new;
            double E_new = (rho_tilde[j][i] * E_total[j][i] + DE[j][i]) / rho_new;
            double I_new = E_new - 0.5 * (u_new * u_new + v_new * v_new);
            if (I_new < 0.0) I_new = 0.0;

            // Обновление W
            double W_new = (rho_tilde[j][i] * W[j][i] + DW[j][i]) / rho_new;
            if (W_new < 0.0) W_new = 0.0;
            if (W_new > 1.0) W_new = 1.0;

            rho[j][i] = rho_new;
            u[j][i] = u_new;
            v[j][i] = v_new;
            I[j][i] = I_new;
            W[j][i] = W_new;
        }
    }

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
/*
std::vector<std::vector<bool>> isSolid(Ny, std::vector<bool>(Nx, false));
for (int j = ghost_cells; j < Ny-ghost_cells; ++j) {
    for (int i = ghost_cells; i < Nx-ghost_cells; ++i) {
	
	if(i*hx  < 0.01 && j*hy < 0.02 ){
	    rho[j][i] = 4.0;
            v[j][i]   = 0;
	    u[j][i]   = 0;
            P[j][i]   = 0.1;
	    isSolid[j][i] = true;
    }    
    }
}
for (int j = ghost_cells; j < Ny - ghost_cells; ++j) {
    for (int i = ghost_cells; i < Nx - ghost_cells; ++i) {
        if (!isSolid[j][i]) { 
            if (isSolid[j][i+1]) {
                rho[j][i+1] = rho[j][i];
                P[j][i+1]   = P[j][i];
                u[j][i+1]   = -u[j][i];
                v[j][i+1]   = v[j][i];
            }
            if (isSolid[j][i-1]) {
                rho[j][i-1] = rho[j][i];
                P[j][i-1]   = P[j][i];
                u[j][i-1]   = -u[j][i];
                v[j][i-1]   = v[j][i];
            }
            // Сверху
            if (isSolid[j+1][i]) {
                rho[j+1][i] = rho[j][i];
                P[j+1][i]   = P[j][i];
                u[j+1][i]   = u[j][i];
                v[j+1][i]   = -v[j][i];
            }
            // Снизу
            if (isSolid[j-1][i]) {
                rho[j-1][i] = rho[j][i];
                P[j-1][i]   = P[j][i];
                u[j-1][i]   = u[j][i];
                v[j-1][i]   = -v[j][i];
            }
        }
    }
}
*/

    std::cout << t_total << " " << dt << std::endl;
    t_total += 2*dt;
}



const double RHO = 1.0;      
const double NU = 0.01;       

inline double interp_x(const std::vector<std::vector<double>>& f, int i, int j) {
    return 0.5 * (f[j][i] + f[j][i+1]);
}
inline double interp_y(const std::vector<std::vector<double>>& f, int i, int j) {
    return 0.5 * (f[j][i] + f[j+1][i]);
}


void applyBoundaryConditions(
    std::vector<std::vector<double>>& u,   
    std::vector<std::vector<double>>& v,   
    std::vector<std::vector<double>>& P,   
    std::vector<std::vector<double>>& rho,
    int cells_x, int cells_y, int ghost_cells,
    const std::string& left_bc, const std::string& right_bc,
    const std::string& bottom_bc, const std::string& top_bc,
    double U_lid,
    double t_current){    
    double nu=NU;        

    int Nx = cells_x + 2 * ghost_cells;
    int Ny = cells_y + 2 * ghost_cells;
    int Nu = Nx + 1;
    int Nv = Ny + 1;
    double dx = 1.0 / cells_x;   // область [0,1]x[0,1]
    double dy = 1.0 / cells_y;
    (void)rho;
    double decay = exp(-2.0 * M_PI * M_PI * nu * t_current);

    if (left_bc == "outflow") {
        for (int j = 0; j < Ny; ++j) {
            for (int g = 0; g < ghost_cells; ++g) {
                u[j][g] = u[j][ghost_cells];
                v[j][g] = v[j][ghost_cells];
                P[j][g] = P[j][ghost_cells];
            }
        }
    } else if (left_bc == "wall") {
        for (int j = 0; j < Ny; ++j) {
            for (int g = 0; g < ghost_cells; ++g) {
                u[j][g] = -u[j][ghost_cells];
                v[j][g] = -v[j][ghost_cells];
                P[j][g] = P[j][ghost_cells];
            }
        }
    } else if (left_bc == "inlet") {
        for (int j = 0; j < Ny; ++j) {
            for (int g = 0; g < ghost_cells; ++g) {
                u[j][g] = 0.0;
                v[j][g] = 1.0;
                P[j][g] = P[j][ghost_cells];
            }
        }
    } else if (left_bc == "taylorGreen") {
        for (int j = 0; j < Ny; ++j) {
            for (int g = 0; g < ghost_cells; ++g) {
                double y = (j - ghost_cells + 0.5) * dy;
                u[j][g] = -cos(0.0) * sin(M_PI * y) * decay;
                y = (j - ghost_cells) * dy;
                v[j][g] = sin(0.0) * cos(M_PI * y) * decay; 
                P[j][g] = P[j][ghost_cells];  
            }
        }
    } else if (left_bc == "inlet_poiseuille") {
        for (int j = 0; j < Ny; ++j) {
	    double h = 0.5, H = 1, U_avg = 5.5;
            double y = (j - ghost_cells) * dy;
            double u_in = 0.0;
            if (y > h && y < H) {
                u_in = 6.0 * U_avg * (y - h) * (H - y) / ((H - h) * (H - h));
            }
            for (int g = 0; g < ghost_cells; ++g) {
                u[j][g] = u_in;
                v[j][g] = 0.0;
                P[j][g] = P[j][ghost_cells];
            }
        }}
    if (right_bc == "outflow") {
        for (int j = 0; j < Ny; ++j) {
            for (int g = 0; g < ghost_cells; ++g) {
                int idx = Nx - ghost_cells + g;
                int src = Nx - ghost_cells - 1;
                u[j][idx] = u[j][src];
                v[j][idx] = v[j][src];
                P[j][idx] = P[j][src];
            }
        }
    } else if (right_bc == "wall") {
        for (int j = 0; j < Ny; ++j) {
            for (int g = 0; g < ghost_cells; ++g) {
                int idx = Nx - ghost_cells + g;
                int src = Nx - ghost_cells - 1;
                u[j][idx] = -u[j][src];
                v[j][idx] = -v[j][src];
                P[j][idx] = P[j][src];
            }
        }
    } else if (right_bc == "taylorGreen") {
        for (int j = 0; j < Ny; ++j) {
            for (int g = 0; g < ghost_cells; ++g) {
                int idx = Nx - ghost_cells + g;
                double y = (j - ghost_cells + 0.5) * dy;
                u[j][idx] = -cos(M_PI) * sin(M_PI * y) * decay; 
                y = (j - ghost_cells) * dy;
                v[j][idx] = sin(M_PI) * cos(M_PI * y) * decay; 
                P[j][idx] = P[j][Nx - ghost_cells - 1];
            }
        }
    }

    if (bottom_bc == "outflow") {
        for (int i = 0; i < Nx; ++i) {
            for (int g = 0; g < ghost_cells; ++g) {
                u[g][i] = u[ghost_cells][i];
                v[g][i] = v[ghost_cells][i];
                P[g][i] = P[ghost_cells][i];
            }
        }
    } else if (bottom_bc == "wall") {
        for (int i = 0; i < Nx; ++i) {
            for (int g = 0; g < ghost_cells; ++g) {
                u[g][i] = -u[ghost_cells][i];
                v[g][i] = -v[ghost_cells][i];
                P[g][i] = P[ghost_cells][i];
            }
        }
    } else if (bottom_bc == "taylorGreen") {
        // Нижняя граница: y = 0
        for (int i = 0; i < Nx; ++i) {
            for (int g = 0; g < ghost_cells; ++g) {
                double x = (i - ghost_cells + 0.5) * dx;
                u[g][i] = -cos(M_PI * x) * sin(0.0) * decay; 
                x = (i - ghost_cells) * dx;
                v[g][i] = sin(M_PI * x) * cos(0.0) * decay;
                P[g][i] = P[ghost_cells][i];
            }
        }
    }

    if (top_bc == "outflow") {
        for (int i = 0; i < Nx; ++i) {
            for (int g = 0; g < ghost_cells; ++g) {
                int idx = Ny - ghost_cells + g;
                int src = Ny - ghost_cells - 1;
                u[idx][i] = u[src][i];
                v[idx][i] = v[src][i];
                P[idx][i] = P[src][i];
            }
        }
    } else if (top_bc == "wall") {
        for (int i = 0; i < Nx; ++i) {
            for (int g = 0; g < ghost_cells; ++g) {
                int idx = Ny - ghost_cells + g;
                int src = Ny - ghost_cells - 1;
                u[idx][i] = -u[src][i];
                v[idx][i] = -v[src][i];
                P[idx][i] = P[src][i];
            }
        }
    } else if (top_bc == "movingWall") {
        for (int i = 0; i < Nx; ++i) {
            for (int g = 0; g < ghost_cells; ++g) {
                int idx = Ny - ghost_cells + g;
                int src = Ny - ghost_cells - 1;
                u[idx][i] = 2.0 * U_lid - u[src][i];
                v[idx][i] = -v[src][i];
                P[idx][i] = P[src][i];
            }
        }
    } else if (top_bc == "taylorGreen") {
        for (int i = 0; i < Nx; ++i) {
            for (int g = 0; g < ghost_cells; ++g) {
                int idx = Ny - ghost_cells + g;
                double x = (i - ghost_cells + 0.5) * dx;
                u[idx][i] = -cos(M_PI * x) * sin(M_PI) * decay; 
                x = (i - ghost_cells) * dx;
                v[idx][i] = sin(M_PI * x) * cos(M_PI) * decay;
                P[idx][i] = P[Ny - ghost_cells - 1][i];
            }
        }
    }
}

void momentumCoeffsU(
    const std::vector<std::vector<double>>& u,
    const std::vector<std::vector<double>>& v,
    const std::vector<std::vector<double>>& rho,
    int i, int j,
    double hx, double hy, double dt,
    double& ap, double& ae, double& aw, double& an, double& as, double& su) {

    double u_e = 0.5 * (u[j][i] + u[j][i+1]);
    double u_w = 0.5 * (u[j][i] + u[j][i-1]);
    double v_n = 0.5 * (v[j+1][i-1] + v[j+1][i]);
    double v_s = 0.5 * (v[j][i-1]   + v[j][i]);

    double rho_e = interp_x(rho, i, j);
    double rho_w = interp_x(rho, i-1, j);
    double rho_n = interp_y(rho, i, j);
    double rho_s = interp_y(rho, i, j-1);

    double Fe = rho_e * u_e * hy;
    double Fw = rho_w * u_w * hy;
    double Fn = rho_n * v_n * hx;
    double Fs = rho_s * v_s * hx;

    double De = NU * hy / hx;
    double Dw = NU * hy / hx;
    double Dn = NU * hx / hy;
    double Ds = NU * hx / hy;

    ae = std::max(-Fe, 0.0) + De;
    aw = std::max( Fw, 0.0) + Dw;
    an = std::max(-Fn, 0.0) + Dn;
    as = std::max( Fs, 0.0) + Ds;

    ap = ae + aw + an + as;
    if (dt > 0) ap += rho[j][i] * hx * hy / dt;   

    double H = ae * u[j][i+1] + aw * u[j][i-1] + an * u[j+1][i] + as * u[j-1][i];
    if (dt > 0) H += (rho[j][i] * hx * hy / dt) * u[j][i];

    su = H;
}

void momentumCoeffsV(
    const std::vector<std::vector<double>>& u,
    const std::vector<std::vector<double>>& v,
    const std::vector<std::vector<double>>& rho,
    int i, int j,
    double hx, double hy, double dt,
    double& ap, double& ae, double& aw, double& an, double& as, double& su) {

    double u_e = 0.5 * (u[j-1][i+1] + u[j][i+1]);
    double u_w = 0.5 * (u[j-1][i]   + u[j][i]);
    double v_n = 0.5 * (v[j][i] + v[j+1][i]);
    double v_s = 0.5 * (v[j][i] + v[j-1][i]);

    double rho_e = interp_x(rho, i, j);
    double rho_w = interp_x(rho, i-1, j);
    double rho_n = interp_y(rho, i, j);
    double rho_s = interp_y(rho, i, j-1);

    double Fe = rho_e * u_e * hy;
    double Fw = rho_w * u_w * hy;
    double Fn = rho_n * v_n * hx;
    double Fs = rho_s * v_s * hx;

    double De = NU * hy / hx;
    double Dw = NU * hy / hx;
    double Dn = NU * hx / hy;
    double Ds = NU * hx / hy;

    ae = std::max(-Fe, 0.0) + De;
    aw = std::max( Fw, 0.0) + Dw;
    an = std::max(-Fn, 0.0) + Dn;
    as = std::max( Fs, 0.0) + Ds;

    ap = ae + aw + an + as;
    if (dt > 0) ap += rho[j][i] * hx * hy / dt;

    double H = ae * v[j][i+1] + aw * v[j][i-1] + an * v[j+1][i] + as * v[j-1][i];
    if (dt > 0) H += (rho[j][i] * hx * hy / dt) * v[j][i];

    su = H;
}

void updateDensityEnergy(
    std::vector<std::vector<double>>& rho,
    std::vector<std::vector<double>>& I,
    const std::vector<std::vector<double>>& P,
    int cells_x, int cells_y, int ghost_cells) {
    (void)rho; (void)I; (void)P; (void)cells_x; (void)cells_y; (void)ghost_cells;
}

void solvePressureCorrection(
    const std::vector<std::vector<double>>& u_star,
    const std::vector<std::vector<double>>& v_star,
    const std::vector<std::vector<double>>& aP_u_eff,   
    const std::vector<std::vector<double>>& aP_v_eff,  
    std::vector<std::vector<double>>& p_corr,
    int cells_x, int cells_y, int ghost_cells,
    double hx, double hy, double omega) {

    int Nx = cells_x + 2 * ghost_cells;
    int Ny = cells_y + 2 * ghost_cells;

    std::vector<std::vector<double>> aE(Ny, std::vector<double>(Nx, 0.0));
    std::vector<std::vector<double>> aW(Ny, std::vector<double>(Nx, 0.0));
    std::vector<std::vector<double>> aN(Ny, std::vector<double>(Nx, 0.0));
    std::vector<std::vector<double>> aS(Ny, std::vector<double>(Nx, 0.0));
    std::vector<std::vector<double>> aP(Ny, std::vector<double>(Nx, 0.0));
    std::vector<std::vector<double>> b(Ny, std::vector<double>(Nx, 0.0));

    for (int j = ghost_cells; j < cells_y + ghost_cells; ++j) {
        for (int i = ghost_cells; i < cells_x + ghost_cells; ++i) {
            if (i < cells_x + ghost_cells - 1)
                aE[j][i] = hy * hy / aP_u_eff[j][i+1];
            if (i > ghost_cells)
                aW[j][i] = hy * hy / aP_u_eff[j][i];
            if (j < cells_y + ghost_cells - 1)
                aN[j][i] = hx * hx / aP_v_eff[j+1][i];
            if (j > ghost_cells)
                aS[j][i] = hx * hx / aP_v_eff[j][i];

            aP[j][i] = aE[j][i] + aW[j][i] + aN[j][i] + aS[j][i];

            double div = (u_star[j][i+1] - u_star[j][i]) * hy
                       + (v_star[j+1][i] - v_star[j][i]) * hx;
            b[j][i] = -div;
        }
    }

    int fix_j = ghost_cells;
    int fix_i = ghost_cells;
    aP[fix_j][fix_i] = 1e30;
    b[fix_j][fix_i] = 0.0;

    const int max_iter = 500;
    const double tol = 1e-10;
    double residual = 0.0;

    for (int iter = 0; iter < max_iter; ++iter) {
        residual = 0.0;
        for (int j = ghost_cells; j < cells_y + ghost_cells; ++j) {
            for (int i = ghost_cells; i < cells_x + ghost_cells; ++i) {
                double sum = aE[j][i] * p_corr[j][i+1]
                           + aW[j][i] * p_corr[j][i-1]
                           + aN[j][i] * p_corr[j+1][i]
                           + aS[j][i] * p_corr[j-1][i];
                double new_p = (1.0 - omega) * p_corr[j][i] + omega * (sum + b[j][i]) / aP[j][i];
                double diff = new_p - p_corr[j][i];
                residual += diff * diff;
                p_corr[j][i] = new_p;
            }
        }
        if (std::sqrt(residual) < tol) break;
    }
}

void correctVelocityPressure(
    std::vector<std::vector<double>>& u,
    std::vector<std::vector<double>>& v,
    std::vector<std::vector<double>>& P,
    const std::vector<std::vector<double>>& p_corr,
    const std::vector<std::vector<double>>& aP_u_eff,
    const std::vector<std::vector<double>>& aP_v_eff,
    int cells_x, int cells_y, int ghost_cells,
    double hx, double hy, double alpha_p) {

    int Nx = cells_x + 2 * ghost_cells;
    int Ny = cells_y + 2 * ghost_cells;

    for (int j = ghost_cells; j < cells_y + ghost_cells; ++j) {
        for (int i = ghost_cells; i < cells_x + ghost_cells; ++i) {
            double dpdx = (p_corr[j][i] - p_corr[j][i-1]) * hy;
            u[j][i] -= (hy / aP_u_eff[j][i]) * dpdx;
        }
    }

    for (int j = ghost_cells; j < cells_y + ghost_cells; ++j) {
        for (int i = ghost_cells; i < cells_x + ghost_cells; ++i) {
            double dpdy = (p_corr[j][i] - p_corr[j-1][i]) * hx;
            v[j][i] -= (hx / aP_v_eff[j][i]) * dpdy;
        }
    }

    for (int j = ghost_cells; j < cells_y + ghost_cells; ++j) {
        for (int i = ghost_cells; i < cells_x + ghost_cells; ++i) {
            P[j][i] += alpha_p * p_corr[j][i];
        }
    }
}

void simpleStep(
    std::vector<std::vector<double>>& u,
    std::vector<std::vector<double>>& v,
    std::vector<std::vector<double>>& P,
    std::vector<std::vector<double>>& rho,
    int cells_x, int cells_y, int ghost_cells,
    double hx, double hy,
    double alpha_u, double alpha_p,
    const std::string& left_bc, const std::string& right_bc,
    const std::string& bottom_bc, const std::string& top_bc,
    double U_lid, double t_total) {

    int Nx = cells_x + 2 * ghost_cells;
    int Ny = cells_y + 2 * ghost_cells;
    int Nu = Nx + 1;
    int Nv = Ny + 1;

    auto u_old = u;
    auto v_old = v;
    auto P_old = P;

    std::vector<std::vector<double>> aP_u(Ny, std::vector<double>(Nu, 1.0));
    std::vector<std::vector<double>> aP_v(Nv, std::vector<double>(Nx, 1.0));
    std::vector<std::vector<double>> H_u(Ny, std::vector<double>(Nu, 0.0));
    std::vector<std::vector<double>> H_v(Nv, std::vector<double>(Nx, 0.0));

    for (int j = ghost_cells; j < cells_y + ghost_cells; ++j) {
        for (int i = ghost_cells; i < cells_x + ghost_cells; ++i) {
            double ap, ae, aw, an, as, su;
            momentumCoeffsU(u_old, v_old, rho, i, j, hx, hy, 0.0,
                            ap, ae, aw, an, as, su);
            aP_u[j][i] = ap;
            H_u[j][i] = su;

            momentumCoeffsV(u_old, v_old, rho, i, j, hx, hy, 0.0,
                            ap, ae, aw, an, as, su);
            aP_v[j][i] = ap;
            H_v[j][i] = su;
        }
    }

    std::vector<std::vector<double>> u_star(Ny, std::vector<double>(Nu, 0.0));
    std::vector<std::vector<double>> v_star(Nv, std::vector<double>(Nx, 0.0));

    for (int j = ghost_cells; j < cells_y + ghost_cells; ++j) {
        for (int i = ghost_cells; i < cells_x + ghost_cells; ++i) {
            double dpdx = (P_old[j][i] - P_old[j][i-1]) * hy;
            u_star[j][i] = (1 - alpha_u) * u_old[j][i]
                         + (alpha_u / aP_u[j][i]) * (H_u[j][i] - dpdx);
        }
    }

    for (int j = ghost_cells; j < cells_y + ghost_cells; ++j) {
        for (int i = ghost_cells; i < cells_x + ghost_cells; ++i) {
            double dpdy = (P_old[j][i] - P_old[j-1][i]) * hx;
            v_star[j][i] = (1 - alpha_u) * v_old[j][i]
                         + (alpha_u / aP_v[j][i]) * (H_v[j][i] - dpdy);
        }
    }

    applyBoundaryConditions(u_star, v_star, P_old, rho, cells_x, cells_y, ghost_cells,
                            left_bc, right_bc, bottom_bc, top_bc, U_lid, t_total);

    std::vector<std::vector<double>> aP_u_eff = aP_u;
    std::vector<std::vector<double>> aP_v_eff = aP_v;
    for (int j = 0; j < Ny; ++j)
        for (int i = 0; i < Nu; ++i)
            aP_u_eff[j][i] /= alpha_u;
    for (int j = 0; j < Nv; ++j)
        for (int i = 0; i < Nx; ++i)
            aP_v_eff[j][i] /= alpha_u;

    std::vector<std::vector<double>> p_corr(Ny, std::vector<double>(Nx, 0.0));
    solvePressureCorrection(u_star, v_star, aP_u_eff, aP_v_eff, p_corr,
                            cells_x, cells_y, ghost_cells, hx, hy, 1.5);

    u = u_star;
    v = v_star;
    correctVelocityPressure(u, v, P, p_corr, aP_u_eff, aP_v_eff,
                            cells_x, cells_y, ghost_cells, hx, hy, alpha_p);

    applyBoundaryConditions(u, v, P, rho, cells_x, cells_y, ghost_cells,
                            left_bc, right_bc, bottom_bc, top_bc, U_lid, t_total);
}

void pisoStep(
    std::vector<std::vector<double>>& u,
    std::vector<std::vector<double>>& v,
    std::vector<std::vector<double>>& P,
    std::vector<std::vector<double>>& rho,
    int cells_x, int cells_y, int ghost_cells,
    double hx, double hy, double dt, int nCorr,
    const std::string& left_bc, const std::string& right_bc,
    const std::string& bottom_bc, const std::string& top_bc,
    double U_lid, double t_total) {

    int Nx = cells_x + 2 * ghost_cells;
    int Ny = cells_y + 2 * ghost_cells;
    int Nu = Nx + 1;
    int Nv = Ny + 1;

    auto u_n = u;
    auto v_n = v;
    auto P_n = P;

    std::vector<std::vector<double>> aP_u(Ny, std::vector<double>(Nu, 1.0));
    std::vector<std::vector<double>> aP_v(Nv, std::vector<double>(Nx, 1.0));
    std::vector<std::vector<double>> H_u(Ny, std::vector<double>(Nu, 0.0));
    std::vector<std::vector<double>> H_v(Nv, std::vector<double>(Nx, 0.0));

    for (int j = ghost_cells; j < cells_y + ghost_cells; ++j) {
        for (int i = ghost_cells; i < cells_x + ghost_cells; ++i) {
            double ap, ae, aw, an, as, su;
            momentumCoeffsU(u_n, v_n, rho, i, j, hx, hy, dt,
                            ap, ae, aw, an, as, su);
            aP_u[j][i] = ap;
            H_u[j][i] = su;

            momentumCoeffsV(u_n, v_n, rho, i, j, hx, hy, dt,
                            ap, ae, aw, an, as, su);
            aP_v[j][i] = ap;
            H_v[j][i] = su;
        }
    }

    std::vector<std::vector<double>> u_star(Ny, std::vector<double>(Nu, 0.0));
    std::vector<std::vector<double>> v_star(Nv, std::vector<double>(Nx, 0.0));

    for (int j = ghost_cells; j < cells_y + ghost_cells; ++j) {
        for (int i = ghost_cells; i < cells_x + ghost_cells; ++i) {
            double dpdx = (P_n[j][i] - P_n[j][i-1]) * hy;
            u_star[j][i] = (H_u[j][i] - dpdx) / aP_u[j][i];
        }
    }

    for (int j = ghost_cells; j < cells_y + ghost_cells; ++j) {
        for (int i = ghost_cells; i < cells_x + ghost_cells; ++i) {
            double dpdy = (P_n[j][i] - P_n[j-1][i]) * hx;
            v_star[j][i] = (H_v[j][i] - dpdy) / aP_v[j][i];
        }
    }

    applyBoundaryConditions(u_star, v_star, P_n, rho, cells_x, cells_y, ghost_cells,
                            left_bc, right_bc, bottom_bc, top_bc, U_lid, t_total);

    auto aP_u_eff = aP_u;
    auto aP_v_eff = aP_v;

    auto u_m = u_star;
    auto v_m = v_star;
    auto P_m = P_n;

    for (int corr = 0; corr < nCorr; ++corr) {
        std::vector<std::vector<double>> p_corr(Ny, std::vector<double>(Nx, 0.0));
        solvePressureCorrection(u_m, v_m, aP_u_eff, aP_v_eff, p_corr,
                                cells_x, cells_y, ghost_cells, hx, hy, 1.5);

        for (int j = ghost_cells; j < cells_y + ghost_cells; ++j) {
            for (int i = ghost_cells; i < cells_x + ghost_cells; ++i) {
                double dpdx = (p_corr[j][i] - p_corr[j][i-1]) * hy;
                double dpdy = (p_corr[j][i] - p_corr[j-1][i]) * hx;
                u_m[j][i] -= (hy / aP_u_eff[j][i]) * dpdx;
                v_m[j][i] -= (hx / aP_v_eff[j][i]) * dpdy;
                P_m[j][i] += p_corr[j][i];
            }
        }
        applyBoundaryConditions(u_m, v_m, P_m, rho, cells_x, cells_y, ghost_cells,
                                left_bc, right_bc, bottom_bc, top_bc, U_lid, t_total);
    }

    u = u_m;
    v = v_m;
    P = P_m;
}

void pimpleStep(
    std::vector<std::vector<double>>& u,
    std::vector<std::vector<double>>& v,
    std::vector<std::vector<double>>& P,
    std::vector<std::vector<double>>& rho,
    int cells_x, int cells_y, int ghost_cells,
    double hx, double hy, double dt,
    int nOuter, int nCorr,
    double alpha_u, double alpha_p,
    const std::string& left_bc, const std::string& right_bc,
    const std::string& bottom_bc, const std::string& top_bc,
    double U_lid, double t_total) {

    int Nx = cells_x + 2 * ghost_cells;
    int Ny = cells_y + 2 * ghost_cells;
    int Nu = Nx + 1;
    int Nv = Ny + 1;

    auto u_n = u;
    auto v_n = v;
    auto P_n = P;

    std::vector<std::vector<double>> aP_u(Ny, std::vector<double>(Nu, 1.0));
    std::vector<std::vector<double>> aP_v(Nv, std::vector<double>(Nx, 1.0));
    std::vector<std::vector<double>> H_u(Ny, std::vector<double>(Nu, 0.0));
    std::vector<std::vector<double>> H_v(Nv, std::vector<double>(Nx, 0.0));

    for (int j = ghost_cells; j < cells_y + ghost_cells; ++j) {
        for (int i = ghost_cells; i < cells_x + ghost_cells; ++i) {
            double ap, ae, aw, an, as, su;
            momentumCoeffsU(u_n, v_n, rho, i, j, hx, hy, dt,
                            ap, ae, aw, an, as, su);
            aP_u[j][i] = ap;
            H_u[j][i] = su;

            momentumCoeffsV(u_n, v_n, rho, i, j, hx, hy, dt,
                            ap, ae, aw, an, as, su);
            aP_v[j][i] = ap;
            H_v[j][i] = su;
        }
    }

    for (int outer = 0; outer < nOuter; ++outer) {
        double au = (outer < nOuter - 1) ? alpha_u : 1.0;
        double ap_rel = (outer < nOuter - 1) ? alpha_p : 1.0;

        std::vector<std::vector<double>> u_star(Ny, std::vector<double>(Nu, 0.0));
        std::vector<std::vector<double>> v_star(Nv, std::vector<double>(Nx, 0.0));

        for (int j = ghost_cells; j < cells_y + ghost_cells; ++j) {
            for (int i = ghost_cells; i < cells_x + ghost_cells; ++i) {
                double dpdx = (P_n[j][i] - P_n[j][i-1]) * hy;
                u_star[j][i] = (1 - au) * u_n[j][i]
                             + (au / aP_u[j][i]) * (H_u[j][i] - dpdx);
            }
        }

        for (int j = ghost_cells; j < cells_y + ghost_cells; ++j) {
            for (int i = ghost_cells; i < cells_x + ghost_cells; ++i) {
                double dpdy = (P_n[j][i] - P_n[j-1][i]) * hx;
                v_star[j][i] = (1 - au) * v_n[j][i]
                             + (au / aP_v[j][i]) * (H_v[j][i] - dpdy);
            }
        }

        applyBoundaryConditions(u_star, v_star, P_n, rho, cells_x, cells_y, ghost_cells,
                                left_bc, right_bc, bottom_bc, top_bc, U_lid, t_total);

        std::vector<std::vector<double>> aP_u_eff = aP_u;
        std::vector<std::vector<double>> aP_v_eff = aP_v;
        for (int j = 0; j < Ny; ++j)
            for (int i = 0; i < Nu; ++i)
                aP_u_eff[j][i] /= au;
        for (int j = 0; j < Nv; ++j)
            for (int i = 0; i < Nx; ++i)
                aP_v_eff[j][i] /= au;

        auto u_m = u_star;
        auto v_m = v_star;
        auto P_m = P_n;

        for (int corr = 0; corr < nCorr; ++corr) {
            std::vector<std::vector<double>> p_corr(Ny, std::vector<double>(Nx, 0.0));
            solvePressureCorrection(u_m, v_m, aP_u_eff, aP_v_eff, p_corr,
                                    cells_x, cells_y, ghost_cells, hx, hy, 1.5);

            for (int j = ghost_cells; j < cells_y + ghost_cells; ++j) {
                for (int i = ghost_cells; i < cells_x + ghost_cells; ++i) {
                    double dpdx = (p_corr[j][i] - p_corr[j][i-1]) * hy;
                    double dpdy = (p_corr[j][i] - p_corr[j-1][i]) * hx;
                    u_m[j][i] -= (hy / aP_u_eff[j][i]) * dpdx;
                    v_m[j][i] -= (hx / aP_v_eff[j][i]) * dpdy;
                    P_m[j][i] += ap_rel * p_corr[j][i];
                }
            }
            applyBoundaryConditions(u_m, v_m, P_m, rho, cells_x, cells_y, ghost_cells,
                                    left_bc, right_bc, bottom_bc, top_bc, U_lid, t_total);
        }

        u_n = u_m;
        v_n = v_m;
        P_n = P_m;
    }

    u = u_n;
    v = v_n;
    P = P_n;
}

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
    double U_lid) {

    bool transient = (tau > 0);
    double dt = transient ? tau : 1.0;

    applyBoundaryConditions(u, v, P, rho, cells_x, cells_y, ghost_cells,
                            left_bc, right_bc, bottom_bc, top_bc, U_lid, t_total);

    if (solver_type == "SIMPLE") {
        simpleStep(u, v, P, rho, cells_x, cells_y, ghost_cells,
                   hx, hy, alpha_u, alpha_p,
                   left_bc, right_bc, bottom_bc, top_bc, U_lid, t_total);
        t_total += dt;
	std::cout<<t_total<<" "<<dt<<std::endl;
    }
    else if (solver_type == "PISO") {
        pisoStep(u, v, P, rho, cells_x, cells_y, ghost_cells,
                 hx, hy, dt, nCorr,
                 left_bc, right_bc, bottom_bc, top_bc, U_lid, t_total);
        t_total += dt;
	std::cout<<t_total<<" "<<dt<<std::endl;

    }
    else if (solver_type == "PIMPLE") {
        int nOuter = 3;
        int nCorrInner = 1;
        pimpleStep(u, v, P, rho, cells_x, cells_y, ghost_cells,
                   hx, hy, dt, nOuter, nCorrInner,
                   alpha_u, alpha_p,
                   left_bc, right_bc, bottom_bc, top_bc, U_lid, t_total);
        t_total += dt;
	std::cout<<t_total<<" "<<dt<<std::endl;

    }

    updateDensityEnergy(rho, I, P, cells_x, cells_y, ghost_cells);

    applyBoundaryConditions(u, v, P, rho, cells_x, cells_y, ghost_cells,
                            left_bc, right_bc, bottom_bc, top_bc, U_lid, t_total);
    int Nx = cells_x + 2 * ghost_cells;
    int Ny = cells_y + 2 * ghost_cells;
    double dx = 15.0 / cells_x;
    double dy = 1.0 / cells_y;
std::vector<std::vector<bool>> isSolid(Ny, std::vector<bool>(Nx, false));
for (int j = ghost_cells; j < Ny-ghost_cells; ++j) {
    for (int i = ghost_cells; i < Nx-ghost_cells; ++i) {
	double y = j  * dy;
	double x = i  * dx;

	if(y < 0.5 && x<2.5){
	    rho[j][i] = 2.0;
            v[j][i]   = 0;
	    u[j][i]   = 0;
            P[j][i]   = 0.1;
	    isSolid[j][i] = true;
    }    
    }
}
for (int j = ghost_cells; j < Ny - ghost_cells; ++j) {
    for (int i = ghost_cells; i < Nx - ghost_cells; ++i) {

        if (!isSolid[j][i]) { 
            if (isSolid[j][i+1]) {
                rho[j][i+1] = rho[j][i];
                P[j][i+1]   = P[j][i];
                u[j][i+1]   = -u[j][i];
                v[j][i+1]   = -v[j][i];
            }
            if (isSolid[j][i-1]) {
                rho[j][i-1] = rho[j][i];
                P[j][i-1]   = P[j][i];
                u[j][i-1]   = -u[j][i];
                v[j][i-1]   = -v[j][i];
            }
            // Сверху
            if (isSolid[j+1][i]) {
                rho[j+1][i] = rho[j][i];
                P[j+1][i]   = P[j][i];
                u[j+1][i]   = -u[j][i];
                v[j+1][i]   = -v[j][i];
            }
            // Снизу
            if (isSolid[j-1][i]) {
                rho[j-1][i] = rho[j][i];
                P[j-1][i]   = P[j][i];
                u[j-1][i]   = -u[j][i];
                v[j-1][i]   = -v[j][i];
            }
        }
    }
}
}
