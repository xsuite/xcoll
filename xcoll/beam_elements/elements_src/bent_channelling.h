// copyright ############################### #
// This file is part of the Xcoll package.   #
// Copyright (c) CERN, 2025.                 #
// ######################################### #

#ifndef XCOLL_BENT_CHANNELLINGDEV_H
#define XCOLL_BENT_CHANNELLINGDEV_H


#include <math.h>
#include <stdio.h>



// -------- Constants (SI) --------
double aTF = 1.94e-11; //m
double uT = 7.5e-12; //m
double dp = 1.92e-10; //m
double pc = 4.0e11; //eV
double U0mol = 23.9037; //eV
double U0 = 21.34; //eV
double alpha_i = 0.722452;
double beta_i = 0.573481;
double R = 10.0; // m

// -------- More Constants --------
/*gpufun*/
double xc(double dp, double uT){  
    return dp*0.5 - 0.6565*uT; //m
    }

/*gpufun*/
double eta(double dp, double uT){
    double x_c = xc(dp, uT);
    double t =  x_c/(dp*0.5);
    return t*t;
    } 


// we need the root as well
//double sqrt_eta = sqrt(eta(dp, uT));

/*gpufun*/
double d(double beta_i, double aTF, double dp){
    return (beta_i/aTF) * (dp*0.5);
    }

/*gpufun*/
double UN(double U0mol, double alpha_i, double beta_i, double aTF, double dp){
    return 2.0*U0mol*alpha_i/beta_i * exp(-(beta_i/aTF)*(dp*0.5)); //eV
    }

/*gpufun*/ 
double U(double dp, double uT, double U0mol, double alpha_i, double beta_i,
                                     double aTF, double pc){
    double x_c = xc(dp, uT);
    double eta_ = eta(dp, uT);
    double d_ = d(beta_i, aTF, dp);
    double UN_ = UN(U0mol, alpha_i, beta_i, aTF, dp);
    // 1/m^2
    return (UN_/pc) * eta_ * (d_/x_c) * (d_/x_c);
}



// we need the root as well
//double sqrt_U = sqrt(U(dp, uT, U0mol, alpha_i, beta_i, aTF, dp));


/*gpufun*/
double Rc(double dp, double uT,double U0, double pc){
    double eta_ = eta(dp, uT);
    double x_c = xc(dp, uT);
    return (pc/(2.0*eta_*U0)) * x_c; //m
    }


/*gpufun*/
double lambda(double dp, double uT,double U0, double pc){
    double x_c = xc(dp, uT);
    double R_c = Rc(dp, uT, U0, pc);
    return 2.0*PI*sqrt(x_c*R_c); //m
    }
/*gpufun*/
double L(double dp, double uT, double U0, double pc){
    return 6.0*lambda(dp, uT, U0, pc); //m
    }

/*gpufun*/    
double theta_c(double dp, double uT, double U0, double pc){
    return sqrt((2.0*eta(dp, uT)*U0)/pc) * (1.0 - Rc(dp, uT, U0, pc)/R); //rad
    }


// -------- Suzuki – Yoshida  --------
// Order-4 (n=3)
double chi1_4  =  1.351207191959658;
double chi0_4  = -1.702414383919315;

// Order-6 (n=5)
double chi1_6  =  1.174671758089364;
double chi0_6  = -1.349343516178727;

// Order-8 (n=7)
double chi1_8  =  1.116182939325386;
double chi0_8  = -1.232365878650771;

// Order-10 (n=9)
double chi1_10 =  1.087027106299171;
double chi0_10 = -1.174054212598341;

// Order-12 (n=11)
double chi1_12 =  1.069565719632538;
double chi0_12 = -1.139131439265076;

//================================================
//-------Harmonic + Nonlinear and Bend------------
//================================================

/*gpufun*/
void fM2H0(double x0,double theta0, double s, double dp, double uT, double U0mol, double alpha_i, double beta_i,
                                     double aTF, double pc, double* x, double* theta) {
    double sqrt_U = sqrt(U(dp, uT, U0mol, alpha_i, beta_i, aTF, pc));
    double u=sqrt_U*s;
    *x = x0*cos(u)+theta0/sqrt_U*sin(u);
    *theta = -x0*sqrt_U*sin(u)+theta0*cos(u);
    }

/*gpufun*/
void fM2H1(double x0,double theta0, double s, double dp, double uT, double beta_i, double aTF, double R, double* x, double* theta){
    double d_=d(beta_i, aTF, dp);
    double sqrt_eta = sqrt(eta(dp, uT));
    double D = sqrt_eta*d_;
    double U_ = U(dp, uT, U0mol, alpha_i, beta_i, aTF, pc);
    double x_c = xc(dp, uT);
    *x = x0;
    *theta = theta0 + U_*x0*s-U_*x_c/D*sinh(D*x0/x_c)*s-s/R;
    }

/*gpufun*/
void fM2O2v1(double x0, double theta0, double s, double dp, double uT, double U0mol, double alpha_i, double beta_i, double aTF, double pc, double R, double* x, double* theta){
    
    double half = 0.5 * s;

    // 1) drift: (x,θ) <- H0(s/2)(x0,θ0)
    double x1, th1;
    fM2H0(x0, theta0, half, dp, uT, U0mol, alpha_i, beta_i, aTF, pc, &x1, &th1);

    // 2) kick: (x,θ) <- H1(s)(x1,θ1)
    double x2, th2;
    fM2H1(x1, th1, s, dp, uT, beta_i, aTF, R, &x2, &th2);

    // 3) drift: (x,θ) <- H0(s/2)(x2,θ2)
    fM2H0(x2, th2, half, dp, uT, U0mol, alpha_i, beta_i, aTF, pc, x, theta);
    }


/*gpufun*/ 
void fM2O2v2(double x0, double theta0, double s, double dp, double uT, double U0mol, double alpha_i, double beta_i, double aTF, double pc, double R, double* x, double* theta){
    double half = 0.5 * s;

    // 1) kick
    double x1, th1;
    fM2H1(x0, theta0, half, dp, uT, beta_i, aTF, R, &x1, &th1);

    // 2) drift
    double x2, th2;
    fM2H0(x1, th1, s,dp, uT, U0mol, alpha_i, beta_i, aTF, pc, &x2, &th2);

    // 3) kick
    fM2H1(x2, th2, half, dp, uT, beta_i, aTF, R, x, theta); 
    }


// ===================== Version 1 - Higher Orders ====================== 

/*gpufun*/  
void fM2O4v1(double x0, double theta0, double s, double dp, double uT, double U0mol, double alpha_i, double beta_i, double aTF, double pc, double R, double* x, double* theta){
    
    // 1) drift-kick-drift
    double x1, th1;
    fM2O2v1(x0, theta0, chi1_4*s, dp,  uT, U0mol, alpha_i, beta_i, aTF, pc, R, &x1, &th1); 
    // 2) drift-kick-drift
    double x2, th2;
    fM2O2v1(x1, th1, chi0_4*s, dp,  uT, U0mol, alpha_i, beta_i, aTF, pc, R, &x2, &th2); 
    // 3) drift-kick-drift
    fM2O2v1(x2, th2, chi1_4*s, dp,  uT, U0mol, alpha_i, beta_i, aTF, pc, R, x, theta); 
    }


/*gpufun*/ 
void fM2O6v1(double x0, double theta0, double s, double dp, double uT, double U0mol, double alpha_i, double beta_i, double aTF, double pc, double R, double* x, double* theta){
    
    // 1)
    double x1, th1;
    fM2O4v1(x0, theta0, chi1_6 * s, dp, uT, U0mol, alpha_i, beta_i, aTF, pc, R, &x1, &th1);  

    // 2)
    double x2, th2;
    fM2O4v1(x1, th1, chi0_6 * s, dp, uT, U0mol, alpha_i, beta_i, aTF, pc, R, &x2, &th2);

    // 3)
    fM2O4v1(x2, th2, chi1_6 * s, dp, uT, U0mol, alpha_i, beta_i, aTF, pc, R, x, theta);
    }


/*gpufun*/  
void fM2O8v1(double x0, double theta0, double s, double dp, double uT, double U0mol, double alpha_i, double beta_i, double aTF, double pc, double R, double* x, double* theta){
    
    double x1, th1;
    fM2O6v1(x0, theta0, chi1_8 * s, dp, uT, U0mol, alpha_i, beta_i, aTF, pc, R, &x1, &th1);

    double x2, th2;
    fM2O6v1(x1, th1, chi0_8 * s, dp, uT, U0mol, alpha_i, beta_i, aTF, pc, R, &x2, &th2);

    fM2O6v1(x2, th2, chi1_8 * s, dp, uT, U0mol, alpha_i, beta_i, aTF, pc, R, x, theta);
    }


/*gpufun*/  
void fM2O10v1(double x0, double theta0, double s, double dp, double uT, double U0mol, double alpha_i, double beta_i, double aTF, double pc, double R, double* x, double* theta){
    
    double x1, th1;
    fM2O8v1(x0, theta0, chi1_10 * s, dp, uT, U0mol, alpha_i, beta_i, aTF, pc, R, &x1, &th1);

    double x2, th2;
    fM2O8v1(x1, th1, chi0_10 * s, dp, uT, U0mol, alpha_i, beta_i, aTF, pc, R, &x2, &th2);

    fM2O8v1(x2, th2, chi1_10 * s, dp, uT, U0mol, alpha_i, beta_i, aTF, pc, R, x, theta);
    }


/*gpufun*/  
void fM2O12v1(double x0, double theta0, double s, double dp, double uT, double U0mol, double alpha_i, double beta_i, double aTF, double pc, double R, double* x, double* theta){

    double x1, th1;
    fM2O10v1(x0, theta0, chi1_12 * s, dp, uT, U0mol, alpha_i, beta_i, aTF, pc, R, &x1, &th1);

    double x2, th2;
    fM2O10v1(x1, th1, chi0_12 * s, dp, uT, U0mol, alpha_i, beta_i, aTF, pc, R, &x2, &th2);

    fM2O10v1(x2, th2, chi1_12 * s, dp, uT, U0mol, alpha_i, beta_i, aTF, pc, R, x, theta);
    }


// ====================== Version 2 - Higher Orders ======================

/*gpufun*/  
void fM2O4v2(double x0, double theta0, double s, double dp, double uT, double U0mol, double alpha_i, double beta_i, double aTF, double pc, double R, double* x, double* theta){

    double x1, th1;
    fM2O2v2(x0, theta0, chi1_4 * s, dp, uT, U0mol, alpha_i, beta_i, aTF, pc, R, &x1, &th1);  // kick-drift-kick

    double x2, th2;
    fM2O2v2(x1, th1, chi0_4 * s, dp, uT, U0mol, alpha_i, beta_i, aTF, pc, R, &x2, &th2);

    fM2O2v2(x2, th2, chi1_4 * s, dp, uT, U0mol, alpha_i, beta_i, aTF, pc, R, x, theta);
    }


/*gpufun*/  
void fM2O6v2(double x0, double theta0, double s, double dp, double uT, double U0mol, double alpha_i, double beta_i, double aTF, double pc, double R, double* x, double* theta){
    
    double x1, th1;
    fM2O4v2(x0, theta0, chi1_6 * s, dp, uT, U0mol, alpha_i, beta_i, aTF, pc, R, &x1, &th1);

    double x2, th2;
    fM2O4v2(x1, th1, chi0_6 * s, dp, uT, U0mol, alpha_i, beta_i, aTF, pc, R, &x2, &th2);

    fM2O4v2(x2, th2, chi1_6 * s, dp, uT, U0mol, alpha_i, beta_i, aTF, pc, R, x, theta);
    }

/*gpufun*/ 
void fM2O8v2(double x0, double theta0, double s, double dp, double uT, double U0mol, double alpha_i, double beta_i, double aTF, double pc, double R, double* x, double* theta){

    double x1, th1;
    fM2O6v2(x0, theta0, chi1_8 * s, dp, uT, U0mol, alpha_i, beta_i, aTF, pc, R, &x1, &th1);

    double x2, th2;
    fM2O6v2(x1, th1, chi0_8 * s, dp, uT, U0mol, alpha_i, beta_i, aTF, pc, R, &x2, &th2);

    fM2O6v2(x2, th2, chi1_8 * s, dp, uT, U0mol, alpha_i, beta_i, aTF, pc, R, x, theta);
    }


/*gpufun*/  
void fM2O10v2(double x0, double theta0, double s, double dp, double uT, double U0mol, double alpha_i, double beta_i, double aTF, double pc, double R, double* x, double* theta){

    double x1, th1;
    fM2O8v2(x0, theta0, chi1_10 * s, dp, uT, U0mol, alpha_i, beta_i, aTF, pc, R, &x1, &th1);

    double x2, th2;
    fM2O8v2(x1, th1, chi0_10 * s, dp, uT, U0mol, alpha_i, beta_i, aTF, pc, R, &x2, &th2);

    fM2O8v2(x2, th2, chi1_10 * s, dp, uT, U0mol, alpha_i, beta_i, aTF, pc, R, x, theta);
    }

/*gpufun*/
void fM2O12v2(double x0, double theta0, double s, double dp, double uT, double U0mol, double alpha_i, double beta_i, double aTF, double pc, double R, double* x, double* theta){

    double x1, th1;
    fM2O10v2(x0, theta0, chi1_12 * s, dp, uT, U0mol, alpha_i, beta_i, aTF, pc, R, &x1, &th1);

    double x2, th2;
    fM2O10v2(x1, th1, chi0_12 * s, dp, uT, U0mol, alpha_i, beta_i, aTF, pc, R, &x2, &th2);

    fM2O10v2(x2, th2, chi1_12 * s, dp, uT, U0mol, alpha_i, beta_i, aTF, pc, R, x, theta);
    }




//================================================
//--------Harmonic and Bend + Nonlinear-----------
//================================================

/*gpufun*/
void fM3H0(double x0, double theta0, double s,
           double dp, double uT, double U0mol, double alpha_i, double beta_i,
           double aTF, double pc, double R,
           double* x, double* theta)
{
    double sqrt_U = sqrt(U(dp, uT, U0mol, alpha_i, beta_i, aTF, pc));
    double u = sqrt_U * s;

    double inv_UR  = 1.0 / ( (sqrt_U*sqrt_U) * R );
    double inv_sUR = 1.0 / ( sqrt_U * R );

    *x     = (x0 + inv_UR) * cos(u) + (theta0 / sqrt_U) * sin(u) - inv_UR;
    *theta = (-x0 * sqrt_U - inv_sUR) * sin(u) + theta0 * cos(u);
}

/*gpufun*/
void fM3H1(double x0, double theta0, double s,
           double dp, double uT, double U0mol, double alpha_i, double beta_i,
           double aTF, double pc,
           double* x, double* theta)
{
    double d_        = d(beta_i, aTF, dp);
    double sqrt_eta_ = sqrt(eta(dp, uT));
    double D         = sqrt_eta_ * d_;
    double U_        = U(dp, uT, U0mol, alpha_i, beta_i, aTF, pc);
    double x_c       = xc(dp, uT);

    *x     = x0;
    *theta = theta0 + ( U_ * x0 - U_ * x_c / D * sinh(D * x0 / x_c) ) * s;
}

/*gpufun*/
void fM3O2v1(double x0, double theta0, double s,
             double dp, double uT, double U0mol, double alpha_i, double beta_i,
             double aTF, double pc, double R,
             double* x, double* theta)
{
    double half = 0.5 * s;

    double x1, th1;
    fM3H0(x0, theta0, half, dp, uT, U0mol, alpha_i, beta_i, aTF, pc, R, &x1, &th1);

    double x2, th2;
    fM3H1(x1, th1, s,    dp, uT, U0mol, alpha_i, beta_i, aTF, pc,      &x2, &th2);

    fM3H0(x2, th2, half, dp, uT, U0mol, alpha_i, beta_i, aTF, pc, R,     x,  theta);
}

/*gpufun*/
void fM3O2v2(double x0, double theta0, double s,
             double dp, double uT, double U0mol, double alpha_i, double beta_i,
             double aTF, double pc, double R,
             double* x, double* theta)
{
    double half = 0.5 * s;

    double x1, th1;
    fM3H1(x0, theta0, half, dp, uT, U0mol, alpha_i, beta_i, aTF, pc,      &x1, &th1);

    double x2, th2;
    fM3H0(x1, th1, s,    dp, uT, U0mol, alpha_i, beta_i, aTF, pc, R,      &x2, &th2);

    fM3H1(x2, th2, half, dp, uT, U0mol, alpha_i, beta_i, aTF, pc,           x,  theta);
}

/*gpufun*/
void fM3O4v1(double x0, double theta0, double s,
             double dp, double uT, double U0mol, double alpha_i, double beta_i,
             double aTF, double pc, double R,
             double* x, double* theta)
{
    double x1, th1;
    fM3O2v1(x0, theta0, chi1_4 * s, dp, uT, U0mol, alpha_i, beta_i, aTF, pc, R, &x1, &th1);

    double x2, th2;
    fM3O2v1(x1, th1, chi0_4 * s, dp, uT, U0mol, alpha_i, beta_i, aTF, pc, R, &x2, &th2);

    fM3O2v1(x2, th2, chi1_4 * s, dp, uT, U0mol, alpha_i, beta_i, aTF, pc, R,  x,  theta);
}

/*gpufun*/
void fM3O6v1(double x0, double theta0, double s,
             double dp, double uT, double U0mol, double alpha_i, double beta_i,
             double aTF, double pc, double R,
             double* x, double* theta)
{
    double x1, th1;
    fM3O4v1(x0, theta0, chi1_6 * s, dp, uT, U0mol, alpha_i, beta_i, aTF, pc, R, &x1, &th1);

    double x2, th2;
    fM3O4v1(x1, th1, chi0_6 * s, dp, uT, U0mol, alpha_i, beta_i, aTF, pc, R, &x2, &th2);

    fM3O4v1(x2, th2, chi1_6 * s, dp, uT, U0mol, alpha_i, beta_i, aTF, pc, R,  x,  theta);
}

/*gpufun*/
void fM3O8v1(double x0, double theta0, double s,
             double dp, double uT, double U0mol, double alpha_i, double beta_i,
             double aTF, double pc, double R,
             double* x, double* theta)
{
    double x1, th1;
    fM3O6v1(x0, theta0, chi1_8 * s, dp, uT, U0mol, alpha_i, beta_i, aTF, pc, R, &x1, &th1);

    double x2, th2;
    fM3O6v1(x1, th1, chi0_8 * s, dp, uT, U0mol, alpha_i, beta_i, aTF, pc, R, &x2, &th2);

    fM3O6v1(x2, th2, chi1_8 * s, dp, uT, U0mol, alpha_i, beta_i, aTF, pc, R,  x,  theta);
}

/*gpufun*/
void fM3O10v1(double x0, double theta0, double s,
              double dp, double uT, double U0mol, double alpha_i, double beta_i,
              double aTF, double pc, double R,
              double* x, double* theta)
{
    double x1, th1;
    fM3O8v1(x0, theta0, chi1_10 * s, dp, uT, U0mol, alpha_i, beta_i, aTF, pc, R, &x1, &th1);

    double x2, th2;
    fM3O8v1(x1, th1, chi0_10 * s, dp, uT, U0mol, alpha_i, beta_i, aTF, pc, R, &x2, &th2);

    fM3O8v1(x2, th2, chi1_10 * s, dp, uT, U0mol, alpha_i, beta_i, aTF, pc, R,  x,  theta);
}

/*gpufun*/
void fM3O12v1(double x0, double theta0, double s,
              double dp, double uT, double U0mol, double alpha_i, double beta_i,
              double aTF, double pc, double R,
              double* x, double* theta)
{
    double x1, th1;
    fM3O10v1(x0, theta0, chi1_12 * s, dp, uT, U0mol, alpha_i, beta_i, aTF, pc, R, &x1, &th1);

    double x2, th2;
    fM3O10v1(x1, th1, chi0_12 * s, dp, uT, U0mol, alpha_i, beta_i, aTF, pc, R, &x2, &th2);

    fM3O10v1(x2, th2, chi1_12 * s, dp, uT, U0mol, alpha_i, beta_i, aTF, pc, R,  x,  theta);
}

/*gpufun*/
void fM3O4v2(double x0, double theta0, double s,
             double dp, double uT, double U0mol, double alpha_i, double beta_i,
             double aTF, double pc, double R,
             double* x, double* theta)
{
    double x1, th1;
    fM3O2v2(x0, theta0, chi1_4 * s, dp, uT, U0mol, alpha_i, beta_i, aTF, pc, R, &x1, &th1);

    double x2, th2;
    fM3O2v2(x1, th1, chi0_4 * s, dp, uT, U0mol, alpha_i, beta_i, aTF, pc, R, &x2, &th2);

    fM3O2v2(x2, th2, chi1_4 * s, dp, uT, U0mol, alpha_i, beta_i, aTF, pc, R,  x,  theta);
}

/*gpufun*/
void fM3O6v2(double x0, double theta0, double s,
             double dp, double uT, double U0mol, double alpha_i, double beta_i,
             double aTF, double pc, double R,
             double* x, double* theta)
{
    double x1, th1;
    fM3O4v2(x0, theta0, chi1_6 * s, dp, uT, U0mol, alpha_i, beta_i, aTF, pc, R, &x1, &th1);

    double x2, th2;
    fM3O4v2(x1, th1, chi0_6 * s, dp, uT, U0mol, alpha_i, beta_i, aTF, pc, R, &x2, &th2);

    fM3O4v2(x2, th2, chi1_6 * s, dp, uT, U0mol, alpha_i, beta_i, aTF, pc, R,  x,  theta);
}

/*gpufun*/
void fM3O8v2(double x0, double theta0, double s,
             double dp, double uT, double U0mol, double alpha_i, double beta_i,
             double aTF, double pc, double R,
             double* x, double* theta)
{
    double x1, th1;
    fM3O6v2(x0, theta0, chi1_8 * s, dp, uT, U0mol, alpha_i, beta_i, aTF, pc, R, &x1, &th1);

    double x2, th2;
    fM3O6v2(x1, th1, chi0_8 * s, dp, uT, U0mol, alpha_i, beta_i, aTF, pc, R, &x2, &th2);

    fM3O6v2(x2, th2, chi1_8 * s, dp, uT, U0mol, alpha_i, beta_i, aTF, pc, R,  x,  theta);
}

/*gpufun*/
void fM3O10v2(double x0, double theta0, double s,
              double dp, double uT, double U0mol, double alpha_i, double beta_i,
              double aTF, double pc, double R,
              double* x, double* theta)
{
    double x1, th1;
    fM3O8v2(x0, theta0, chi1_10 * s, dp, uT, U0mol, alpha_i, beta_i, aTF, pc, R, &x1, &th1);

    double x2, th2;
    fM3O8v2(x1, th1, chi0_10 * s, dp, uT, U0mol, alpha_i, beta_i, aTF, pc, R, &x2, &th2);

    fM3O8v2(x2, th2, chi1_10 * s, dp, uT, U0mol, alpha_i, beta_i, aTF, pc, R,  x,  theta);
}

/*gpufun*/
void fM3O12v2(double x0, double theta0, double s,
              double dp, double uT, double U0mol, double alpha_i, double beta_i,
              double aTF, double pc, double R,
              double* x, double* theta)
{
    double x1, th1;
    fM3O10v2(x0, theta0, chi1_12 * s, dp, uT, U0mol, alpha_i, beta_i, aTF, pc, R, &x1, &th1);

    double x2, th2;
    fM3O10v2(x1, th1, chi0_12 * s, dp, uT, U0mol, alpha_i, beta_i, aTF, pc, R, &x2, &th2);

    fM3O10v2(x2, th2, chi1_12 * s, dp, uT, U0mol, alpha_i, beta_i, aTF, pc, R,  x,  theta);
}

//================================================
//----------Simplified Moliere + Bend-------------
//================================================

//------Functions from channelling.h--------------

/*gpufun*/
double U_simplemoliere(double x, double U_N, double beta_i_over_a_TF) {
    //return U_N*(cosh(beta_i_over_a_TF*x) - 1.0);
    // In order to be stable near zero, I use the other form
    double sinh_;
    sinh_=sinh(beta_i_over_a_TF*x/2);
    return 2*U_N*sinh_*sinh_;
    }

/*gpufun*/
double E_T_simplemoliere(double x, double theta, double bpc, double U_N, double beta_i_over_a_TF) {
    return bpc*theta*theta/2.0 + U_simplemoliere(x, U_N, beta_i_over_a_TF);
    }

/*gpufun*/
double m_simplemoliere(double U_N, double E_T ) {
    return 1.0 + 2.0*U_N /E_T;
    }

/*gpufun*/
double mp_simplemoliere(double U_N, double E_T ) {
    return E_T/(E_T + 2.0*U_N);
    }

/*gpufun*/
double nu_simplemoliere(double theta, double bpc, double beta_i_over_a_TF, double E_T) {
    // in unit 1.E4
    double sign = -1;
    if (theta < 0) {
        sign = 1;
        }
    return sign*beta_i_over_a_TF*sqrt(E_T/(2.0*bpc));
}

/*gpufun*/
double phi_simplemoliere(double x, double theta, double nu, double bpc, double U_N,
                         double beta_i_over_a_TF, double m, double mp, double sqrt_mp) {
    double U = U_simplemoliere(x, U_N, beta_i_over_a_TF);
    double sign = -1;
    if (x < 0) {
        sign = 1;
    }
    //if (fabs(theta)< 1e-14) {
        // If theta is very close to 0, we can avoid numerical issues with asin(1)
        //return sign*sqrt_mp*ellpk(mp);
    //}
    double phi_amplitude = asin(sqrt(m*U/(U+2*U_N)));
    // double alpha = 1/sqrt(m*pow(sinh(x*beta_i / (2*a_TF)), 2.) + pow(theta*beta_i /(2*a_TF*nu), 2.));
    // double phi_amplitude = asin(sqrt((1 - alpha*alpha)/mp));
    return sign*sqrt_mp*ellik(phi_amplitude, mp);
}

/*gpufun*/
void motion_parameters(double x0, double theta0, double z, double nu, double E_T, double U_N, double twotimes_a_TF_over_beta_i, double phi, double m, double mp, double sqrt_mp,
     double* x, double* theta) {
    if (E_T == 0) {
        *x = 0;
        *theta = 0;
    }
    else {
        //double u = sqrt(m) * (nu * 1.e4 * z + phi);
         double u = sqrt(m)*(nu*z + phi);
        double sn, cn, dn, am;
        ellpj(u, mp, &sn, &cn, &dn, &am);

        *x = -twotimes_a_TF_over_beta_i * asinh(sqrt_mp * sn / dn);
        *theta = -twotimes_a_TF_over_beta_i * nu * cn / dn;
    }
}
//-----------------------------------------------------------


/*gpufun*/
void fM4H0(double x0, double theta0, double s,
           double dp, double uT, double U0mol, double alpha_i, double beta_i,
           double aTF, double pc, double R, // R and pc unused
           double* x, double* theta)
{
    // Precompute constants in the exact same style as your working code
    const double bpc = 150.0e9;                     // 150 GeV in eV
    const double beta_i_over_a_TF = beta_i / aTF;
    const double U_N = UN(U0mol, alpha_i, beta_i, aTF, dp);

    const double E_T = E_T_simplemoliere(x0, theta0, bpc, U_N, beta_i_over_a_TF);
    
    // m, mp, sqrt_mp as in your helpers
    const double m       = m_simplemoliere(U_N, E_T);       // = 1 + 2U_N/E_T
    const double mp      = mp_simplemoliere(U_N, E_T);      // = E_T/(E_T+2U_N) = 1/m
    const double sqrt_mp = sqrt(mp);

    // ν(θ) with your sign convention and units
    const double nu = nu_simplemoliere(theta0, bpc, beta_i_over_a_TF, E_T);

    // ϕ(x,θ) via your ellik-based helper
    const double phi = phi_simplemoliere(x0, theta0, nu, bpc, U_N,
                                         beta_i_over_a_TF, m, mp, sqrt_mp);

    // 2*aTF/beta_i (you use this in motion_parameters)
    const double twotimes_a_TF_over_beta_i = 2.0 * aTF / beta_i;

    // Exact nonlinear evolution (sn/cn/dn via ellpj inside)
    motion_parameters(x0, theta0, s,               // z = s (meters)
                      nu, E_T, U_N,                // energy & frequency
                      twotimes_a_TF_over_beta_i,
                      phi, m, mp, sqrt_mp,
                       x,  theta);
}

/*gpufun*/
void fM4H1(double x0, double theta0, double s,
           double dp, double uT, double U0mol, double alpha_i, double beta_i,
           double aTF, double pc, double R,
           double* x, double* theta)
{
    // Pure bending kick for Method-4
    *x     = x0;
    *theta = theta0 - s / R;
}


/*gpufun*/
void fM4O2v2(double x0, double theta0, double s,
             double dp, double uT, double U0mol, double alpha_i, double beta_i,
             double aTF, double pc, double R,
             double* x, double* theta)
{
    double half = 0.5 * s;

    // 1) kick
    double x1, th1;
    fM4H1(x0, theta0, half, dp, uT, U0mol, alpha_i, beta_i, aTF, pc, R, &x1, &th1);

    // 2) drift
    double x2, th2;
    fM4H0(x1, th1, s,    dp, uT, U0mol, alpha_i, beta_i, aTF, pc, R, &x2, &th2);

    // 3) kick
    fM4H1(x2, th2, half, dp, uT, U0mol, alpha_i, beta_i, aTF, pc, R,  x,   theta);
}

/* ===================== Higher Orders v1 (drift–kick–drift inside) ===================== */

/*gpufun*/
void fM4O2v1(double x0, double theta0, double s,
             double dp, double uT, double U0mol, double alpha_i, double beta_i,
             double aTF, double pc, double R,
             double* x, double* theta)
{
    double half = 0.5 * s;

    double x1, th1;
    fM4H0(x0, theta0, half, dp, uT, U0mol, alpha_i, beta_i, aTF, pc, R, &x1, &th1);

    double x2, th2;
    fM4H1(x1, th1, s,   dp, uT, U0mol, alpha_i, beta_i, aTF, pc, R, &x2, &th2);

    fM4H0(x2, th2, half, dp, uT, U0mol, alpha_i, beta_i, aTF, pc, R,  x,   theta);
}

/*gpufun*/
void fM4O4v1(double x0, double theta0, double s,
             double dp, double uT, double U0mol, double alpha_i, double beta_i,
             double aTF, double pc, double R,
             double* x, double* theta)
{
    double x1, th1; fM4O2v1(x0, theta0, chi1_4 * s, dp, uT, U0mol, alpha_i, beta_i, aTF, pc, R, &x1, &th1);
    double x2, th2; fM4O2v1(x1, th1, chi0_4 * s, dp, uT, U0mol, alpha_i, beta_i, aTF, pc, R, &x2, &th2);
                    fM4O2v1(x2, th2, chi1_4 * s, dp, uT, U0mol, alpha_i, beta_i, aTF, pc, R,  x,   theta);
}

/*gpufun*/
void fM4O6v1(double x0, double theta0, double s,
             double dp, double uT, double U0mol, double alpha_i, double beta_i,
             double aTF, double pc, double R,
             double* x, double* theta)
{
    double x1, th1; fM4O4v1(x0, theta0, chi1_6 * s, dp, uT, U0mol, alpha_i, beta_i, aTF, pc, R, &x1, &th1);
    double x2, th2; fM4O4v1(x1, th1, chi0_6 * s, dp, uT, U0mol, alpha_i, beta_i, aTF, pc, R, &x2, &th2);
                    fM4O4v1(x2, th2, chi1_6 * s, dp, uT, U0mol, alpha_i, beta_i, aTF, pc, R,  x,   theta);
}

/*gpufun*/
void fM4O8v1(double x0, double theta0, double s,
             double dp, double uT, double U0mol, double alpha_i, double beta_i,
             double aTF, double pc, double R,
             double* x, double* theta)
{
    double x1, th1; fM4O6v1(x0, theta0, chi1_8 * s, dp, uT, U0mol, alpha_i, beta_i, aTF, pc, R, &x1, &th1);
    double x2, th2; fM4O6v1(x1, th1, chi0_8 * s, dp, uT, U0mol, alpha_i, beta_i, aTF, pc, R, &x2, &th2);
                    fM4O6v1(x2, th2, chi1_8 * s, dp, uT, U0mol, alpha_i, beta_i, aTF, pc, R,  x,   theta);
}

/*gpufun*/
void fM4O10v1(double x0, double theta0, double s,
              double dp, double uT, double U0mol, double alpha_i, double beta_i,
              double aTF, double pc, double R,
              double* x, double* theta)
{
    double x1, th1; fM4O8v1 (x0, theta0, chi1_10 * s, dp, uT, U0mol, alpha_i, beta_i, aTF, pc, R, &x1, &th1);
    double x2, th2; fM4O8v1 (x1, th1, chi0_10 * s, dp, uT, U0mol, alpha_i, beta_i, aTF, pc, R, &x2, &th2);
                    fM4O8v1 (x2, th2, chi1_10 * s, dp, uT, U0mol, alpha_i, beta_i, aTF, pc, R,  x,   theta);
}

/*gpufun*/
void fM4O12v1(double x0, double theta0, double s,
              double dp, double uT, double U0mol, double alpha_i, double beta_i,
              double aTF, double pc, double R,
              double* x, double* theta)
{
    double x1, th1; fM4O10v1(x0, theta0, chi1_12 * s, dp, uT, U0mol, alpha_i, beta_i, aTF, pc, R, &x1, &th1);
    double x2, th2; fM4O10v1(x1, th1, chi0_12 * s, dp, uT, U0mol, alpha_i, beta_i, aTF, pc, R, &x2, &th2);
                    fM4O10v1(x2, th2, chi1_12 * s, dp, uT, U0mol, alpha_i, beta_i, aTF, pc, R,  x,   theta);
}

/* ====================== Higher Orders v2 (kick–drift–kick inside) ====================== */

/*gpufun*/
void fM4O4v2(double x0, double theta0, double s,
             double dp, double uT, double U0mol, double alpha_i, double beta_i,
             double aTF, double pc, double R,
             double* x, double* theta)
{
    double x1, th1; fM4O2v2(x0, theta0, chi1_4 * s, dp, uT, U0mol, alpha_i, beta_i, aTF, pc, R, &x1, &th1);
    double x2, th2; fM4O2v2(x1, th1, chi0_4 * s, dp, uT, U0mol, alpha_i, beta_i, aTF, pc, R, &x2, &th2);
                    fM4O2v2(x2, th2, chi1_4 * s, dp, uT, U0mol, alpha_i, beta_i, aTF, pc, R,  x,   theta);
}

/*gpufun*/
void fM4O6v2(double x0, double theta0, double s,
             double dp, double uT, double U0mol, double alpha_i, double beta_i,
             double aTF, double pc, double R,
             double* x, double* theta)
{
    double x1, th1; fM4O4v2(x0, theta0, chi1_6 * s, dp, uT, U0mol, alpha_i, beta_i, aTF, pc, R, &x1, &th1);
    double x2, th2; fM4O4v2(x1, th1, chi0_6 * s, dp, uT, U0mol, alpha_i, beta_i, aTF, pc, R, &x2, &th2);
                    fM4O4v2(x2, th2, chi1_6 * s, dp, uT, U0mol, alpha_i, beta_i, aTF, pc, R,  x,   theta);
}

/*gpufun*/
void fM4O8v2(double x0, double theta0, double s,
             double dp, double uT, double U0mol, double alpha_i, double beta_i,
             double aTF, double pc, double R,
             double* x, double* theta)
{
    double x1, th1; fM4O6v2(x0, theta0, chi1_8 * s, dp, uT, U0mol, alpha_i, beta_i, aTF, pc, R, &x1, &th1);
    double x2, th2; fM4O6v2(x1, th1, chi0_8 * s, dp, uT, U0mol, alpha_i, beta_i, aTF, pc, R, &x2, &th2);
                    fM4O6v2(x2, th2, chi1_8 * s, dp, uT, U0mol, alpha_i, beta_i, aTF, pc, R,  x,   theta);
}

/*gpufun*/
void fM4O10v2(double x0, double theta0, double s,
              double dp, double uT, double U0mol, double alpha_i, double beta_i,
              double aTF, double pc, double R,
              double* x, double* theta)
{
    double x1, th1; fM4O8v2 (x0, theta0, chi1_10 * s, dp, uT, U0mol, alpha_i, beta_i, aTF, pc, R, &x1, &th1);
    double x2, th2; fM4O8v2 (x1, th1, chi0_10 * s, dp, uT, U0mol, alpha_i, beta_i, aTF, pc, R, &x2, &th2);
                    fM4O8v2 (x2, th2, chi1_10 * s, dp, uT, U0mol, alpha_i, beta_i, aTF, pc, R,  x,   theta);
}

/*gpufun*/
void fM4O12v2(double x0, double theta0, double s,
              double dp, double uT, double U0mol, double alpha_i, double beta_i,
              double aTF, double pc, double R,
              double* x, double* theta)
{
    double x1, th1; fM4O10v2(x0, theta0, chi1_12 * s, dp, uT, U0mol, alpha_i, beta_i, aTF, pc, R, &x1, &th1);
    double x2, th2; fM4O10v2(x1, th1, chi0_12 * s, dp, uT, U0mol, alpha_i, beta_i, aTF, pc, R, &x2, &th2);
                    fM4O10v2(x2, th2, chi1_12 * s, dp, uT, U0mol, alpha_i, beta_i, aTF, pc, R,  x,   theta);
}

/*gpufun*/
void choose_method(
    double x0, double theta0, double s,
    double dp, double uT, double U0mol, double alpha_i, double beta_i,
    double aTF, double pc, double R,
    int meth_sel, int var_sel, int ord_sel,
    double* x, double* theta)
{
    #ifdef method
    #undef method
    #endif
    #ifdef variant
    #undef variant
    #endif
    #ifdef order
    #undef order
    #endif

    if (meth_sel == 2) {

        if (var_sel == 1) {
            switch (ord_sel) {
                case 12: fM2O12v1(x0,theta0,s,dp,uT,U0mol,alpha_i,beta_i,aTF,pc,R,x,theta); break;
                case 10: fM2O10v1(x0,theta0,s,dp,uT,U0mol,alpha_i,beta_i,aTF,pc,R,x,theta); break;
                case  8: fM2O8v1 (x0,theta0,s,dp,uT,U0mol,alpha_i,beta_i,aTF,pc,R,x,theta); break;
                case  6: fM2O6v1 (x0,theta0,s,dp,uT,U0mol,alpha_i,beta_i,aTF,pc,R,x,theta); break;
                case  4: fM2O4v1 (x0,theta0,s,dp,uT,U0mol,alpha_i,beta_i,aTF,pc,R,x,theta); break;
                default: fM2O2v1 (x0,theta0,s,dp,uT,U0mol,alpha_i,beta_i,aTF,pc,R,x,theta); break;
            }
        } else {
            switch (ord_sel) {
                case 12: fM2O12v2(x0,theta0,s,dp,uT,U0mol,alpha_i,beta_i,aTF,pc,R,x,theta); break;
                case 10: fM2O10v2(x0,theta0,s,dp,uT,U0mol,alpha_i,beta_i,aTF,pc,R,x,theta); break;
                case  8: fM2O8v2 (x0,theta0,s,dp,uT,U0mol,alpha_i,beta_i,aTF,pc,R,x,theta); break;
                case  6: fM2O6v2 (x0,theta0,s,dp,uT,U0mol,alpha_i,beta_i,aTF,pc,R,x,theta); break;
                case  4: fM2O4v2 (x0,theta0,s,dp,uT,U0mol,alpha_i,beta_i,aTF,pc,R,x,theta); break;
                default: fM2O2v2 (x0,theta0,s,dp,uT,U0mol,alpha_i,beta_i,aTF,pc,R,x,theta); break;
            }
        }

    } else if (meth_sel == 3) {

        if (var_sel == 1) {
            switch (ord_sel) {
                case 12: fM3O12v1(x0,theta0,s,dp,uT,U0mol,alpha_i,beta_i,aTF,pc,R,x,theta); break;
                case 10: fM3O10v1(x0,theta0,s,dp,uT,U0mol,alpha_i,beta_i,aTF,pc,R,x,theta); break;
                case  8: fM3O8v1 (x0,theta0,s,dp,uT,U0mol,alpha_i,beta_i,aTF,pc,R,x,theta); break;
                case  6: fM3O6v1 (x0,theta0,s,dp,uT,U0mol,alpha_i,beta_i,aTF,pc,R,x,theta); break;
                case  4: fM3O4v1 (x0,theta0,s,dp,uT,U0mol,alpha_i,beta_i,aTF,pc,R,x,theta); break;
                default: fM3O2v1 (x0,theta0,s,dp,uT,U0mol,alpha_i,beta_i,aTF,pc,R,x,theta); break;
            }
        } else {
            switch (ord_sel) {
                case 12: fM3O12v2(x0,theta0,s,dp,uT,U0mol,alpha_i,beta_i,aTF,pc,R,x,theta); break;
                case 10: fM3O10v2(x0,theta0,s,dp,uT,U0mol,alpha_i,beta_i,aTF,pc,R,x,theta); break;
                case  8: fM3O8v2 (x0,theta0,s,dp,uT,U0mol,alpha_i,beta_i,aTF,pc,R,x,theta); break;
                case  6: fM3O6v2 (x0,theta0,s,dp,uT,U0mol,alpha_i,beta_i,aTF,pc,R,x,theta); break;
                case  4: fM3O4v2 (x0,theta0,s,dp,uT,U0mol,alpha_i,beta_i,aTF,pc,R,x,theta); break;
                default: fM3O2v2 (x0,theta0,s,dp,uT,U0mol,alpha_i,beta_i,aTF,pc,R,x,theta); break;
            }
        }

    } else { // meth_sel == 4

        if (var_sel == 1) {
            switch (ord_sel) {
                case 12: fM4O12v1(x0,theta0,s,dp,uT,U0mol,alpha_i,beta_i,aTF,pc,R,x,theta); break;
                case 10: fM4O10v1(x0,theta0,s,dp,uT,U0mol,alpha_i,beta_i,aTF,pc,R,x,theta); break;
                case  8: fM4O8v1 (x0,theta0,s,dp,uT,U0mol,alpha_i,beta_i,aTF,pc,R,x,theta); break;
                case  6: fM4O6v1 (x0,theta0,s,dp,uT,U0mol,alpha_i,beta_i,aTF,pc,R,x,theta); break;
                case  4: fM4O4v1 (x0,theta0,s,dp,uT,U0mol,alpha_i,beta_i,aTF,pc,R,x,theta); break;
                default: fM4O2v1 (x0,theta0,s,dp,uT,U0mol,alpha_i,beta_i,aTF,pc,R,x,theta); break;
            }
        } else {
            switch (ord_sel) {
                case 12: fM4O12v2(x0,theta0,s,dp,uT,U0mol,alpha_i,beta_i,aTF,pc,R,x,theta); break;
                case 10: fM4O10v2(x0,theta0,s,dp,uT,U0mol,alpha_i,beta_i,aTF,pc,R,x,theta); break;
                case  8: fM4O8v2 (x0,theta0,s,dp,uT,U0mol,alpha_i,beta_i,aTF,pc,R,x,theta); break;
                case  6: fM4O6v2 (x0,theta0,s,dp,uT,U0mol,alpha_i,beta_i,aTF,pc,R,x,theta); break;
                case  4: fM4O4v2 (x0,theta0,s,dp,uT,U0mol,alpha_i,beta_i,aTF,pc,R,x,theta); break;
                default: fM4O2v2 (x0,theta0,s,dp,uT,U0mol,alpha_i,beta_i,aTF,pc,R,x,theta); break;
            }
        }
    }
}



// -------- Initial conditions (example) --------
//double x0 = 0.97*xc; //m
//double theta0 = 0.0; //rad

/*gpufun*/
void BentChannellingDev_track_local_particle(BentChannellingDevData el, LocalParticle* part0) {
    
    double length = BentChannellingDevData_get_length(el);
    int method = BentChannellingDevData_get_method(el);  
    int order = BentChannellingDevData_get_order(el);
    int variant = BentChannellingDevData_get_variant(el);
    //start_per_particle_block (part0->part)
        double x0 = LocalParticle_get_x(part);     // m
        double theta0 = LocalParticle_get_xp(part); // rad
        double x;
        double theta;
        double s = length;

        //fM2O12v1(x0, theta0, s, dp, uT, U0mol, alpha_i, beta_i, aTF, pc, R, &x, &theta);
        //fM4O12v1(x0, theta0, s, dp, uT, U0mol, alpha_i, beta_i, aTF, pc, R, &x, &theta);
        choose_method(x0, theta0, s, dp, uT, U0mol, alpha_i, beta_i, aTF, pc, R, method, variant, order, &x, &theta);

        LocalParticle_set_x(part, x);
        LocalParticle_set_xp(part, theta);
    
    //end_per_particle_block
    }











#endif /* XCOLL_BENT_CHANNELLINGDEV_H */
