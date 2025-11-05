// copyright ############################### #
// This file is part of the Xcoll package.   #
// Copyright (c) CERN, 2025.                 #
// ######################################### #

#ifndef XCOLL_BENT_CHANNELLINGDEV_M3V1O12_H
#define XCOLL_BENT_CHANNELLINGDEV_M3V1O12_H


#include <math.h>
#include <stdio.h>



// -------- Constants (SI) --------
double aTF = 1.94e-11; //m
double uT = 7.5e-12; //m
double dp = 1.92e-10; //m
double pc = 4.0e11; //eV
double bpc = 150.0e9; // 150 GeV 
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
                                     double aTF, double bpc){
    double x_c = xc(dp, uT);
    double eta_ = eta(dp, uT);
    double d_ = d(beta_i, aTF, dp);
    double UN_ = UN(U0mol, alpha_i, beta_i, aTF, dp);
    // 1/m^2
    return (UN_/bpc) * eta_ * (d_/x_c) * (d_/x_c);
}



// we need the root as well
//double sqrt_U = sqrt(U(dp, uT, U0mol, alpha_i, beta_i, aTF, dp));


/*gpufun*/
double Rc(double dp, double uT,double U0, double bpc){
    double eta_ = eta(dp, uT);
    double x_c = xc(dp, uT);
    return (bpc/(2.0*eta_*U0)) * x_c; //m
    }


/*gpufun*/
double lambda(double dp, double uT,double U0, double bpc){
    double x_c = xc(dp, uT);
    double R_c = Rc(dp, uT, U0, bpc);
    return 2.0*PI*sqrt(x_c*R_c); //m
    }
/*gpufun*/
double L(double dp, double uT, double U0, double bpc){
    return 6.0*lambda(dp, uT, U0, bpc); //m
    }

/*gpufun*/    
double theta_c(double dp, double uT, double U0, double bpc){
    return sqrt((2.0*eta(dp, uT)*U0)/bpc) * (1.0 - Rc(dp, uT, U0, bpc)/R); //rad
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
//--------Harmonic and Bend + Nonlinear-----------
//================================================

/*gpufun*/
void fM3H0(double x0, double theta0, double s,
           double dp, double uT, double U0mol, double alpha_i, double beta_i,
           double aTF, double bpc, double R,
           double* x, double* theta)
{
    double sqrt_U = sqrt(U(dp, uT, U0mol, alpha_i, beta_i, aTF, bpc));
    double u = sqrt_U * s;

    double inv_UR  = 1.0 / ( (sqrt_U*sqrt_U) * R );
    double inv_sUR = 1.0 / ( sqrt_U * R );

    *x     = (x0 + inv_UR)*cos(u) + (theta0/sqrt_U)*sin(u) - inv_UR;
    *theta = (-x0*sqrt_U - inv_sUR)*sin(u) + theta0*cos(u);
}

/*gpufun*/
void fM3H1(double x0, double theta0, double s,
           double dp, double uT, double U0mol, double alpha_i, double beta_i,
           double aTF, double bpc,
           double* x, double* theta)
{
    double d_        = d(beta_i, aTF, dp);
    double sqrt_eta_ = sqrt(eta(dp, uT));
    double D         = sqrt_eta_ * d_;
    double U_        = U(dp, uT, U0mol, alpha_i, beta_i, aTF, bpc);
    double x_c       = xc(dp, uT);

    *x     = x0;
    *theta = theta0 + (U_*x0 - U_*x_c/D*sinh(D*x0/x_c))*s;
}





//----------Integrators for method 3--------------
/*gpufun*/
void fM3O2v1(double x0, double theta0, double s,
             double dp, double uT, double U0mol, double alpha_i, double beta_i,
             double aTF, double bpc, double R,
             double* x, double* theta)
{
    double half = 0.5 * s;

    double x1, th1;
    fM3H0(x0, theta0, half, dp, uT, U0mol, alpha_i, beta_i, aTF, bpc, R, &x1, &th1);

    double x2, th2;
    fM3H1(x1, th1, s,    dp, uT, U0mol, alpha_i, beta_i, aTF, bpc, &x2, &th2);

    fM3H0(x2, th2, half, dp, uT, U0mol, alpha_i, beta_i, aTF, bpc, R, x, theta);
}


/*gpufun*/
void fM3O4v1(double x0, double theta0, double s,
             double dp, double uT, double U0mol, double alpha_i, double beta_i,
             double aTF, double bpc, double R,
             double* x, double* theta)
{
    double x1, th1;
    fM3O2v1(x0, theta0, chi1_4 * s, dp, uT, U0mol, alpha_i, beta_i, aTF, bpc, R, &x1, &th1);

    double x2, th2;
    fM3O2v1(x1, th1, chi0_4 * s, dp, uT, U0mol, alpha_i, beta_i, aTF, bpc, R, &x2, &th2);

    fM3O2v1(x2, th2, chi1_4 * s, dp, uT, U0mol, alpha_i, beta_i, aTF, bpc, R,  x,  theta);
}

/*gpufun*/
void fM3O6v1(double x0, double theta0, double s,
             double dp, double uT, double U0mol, double alpha_i, double beta_i,
             double aTF, double bpc, double R,
             double* x, double* theta)
{
    double x1, th1;
    fM3O4v1(x0, theta0, chi1_6 * s, dp, uT, U0mol, alpha_i, beta_i, aTF, bpc, R, &x1, &th1);

    double x2, th2;
    fM3O4v1(x1, th1, chi0_6 * s, dp, uT, U0mol, alpha_i, beta_i, aTF, bpc, R, &x2, &th2);

    fM3O4v1(x2, th2, chi1_6 * s, dp, uT, U0mol, alpha_i, beta_i, aTF, bpc, R,  x,  theta);
}

/*gpufun*/
void fM3O8v1(double x0, double theta0, double s,
             double dp, double uT, double U0mol, double alpha_i, double beta_i,
             double aTF, double bpc, double R,
             double* x, double* theta)
{
    double x1, th1;
    fM3O6v1(x0, theta0, chi1_8 * s, dp, uT, U0mol, alpha_i, beta_i, aTF, bpc, R, &x1, &th1);

    double x2, th2;
    fM3O6v1(x1, th1, chi0_8 * s, dp, uT, U0mol, alpha_i, beta_i, aTF, bpc, R, &x2, &th2);

    fM3O6v1(x2, th2, chi1_8 * s, dp, uT, U0mol, alpha_i, beta_i, aTF, bpc, R,  x,  theta);
}

/*gpufun*/
void fM3O10v1(double x0, double theta0, double s,
              double dp, double uT, double U0mol, double alpha_i, double beta_i,
              double aTF, double bpc, double R,
              double* x, double* theta)
{
    double x1, th1;
    fM3O8v1(x0, theta0, chi1_10 * s, dp, uT, U0mol, alpha_i, beta_i, aTF, bpc, R, &x1, &th1);

    double x2, th2;
    fM3O8v1(x1, th1, chi0_10 * s, dp, uT, U0mol, alpha_i, beta_i, aTF, bpc, R, &x2, &th2);

    fM3O8v1(x2, th2, chi1_10 * s, dp, uT, U0mol, alpha_i, beta_i, aTF, bpc, R,  x,  theta);
}

/*gpufun*/
void fM3O12v1(double x0, double theta0, double s,
              double dp, double uT, double U0mol, double alpha_i, double beta_i,
              double aTF, double bpc, double R,
              double* x, double* theta)
{
    double x1, th1;
    fM3O10v1(x0, theta0, chi1_12 * s, dp, uT, U0mol, alpha_i, beta_i, aTF, bpc, R, &x1, &th1);

    double x2, th2;
    fM3O10v1(x1, th1, chi0_12 * s, dp, uT, U0mol, alpha_i, beta_i, aTF, bpc, R, &x2, &th2);

    fM3O10v1(x2, th2, chi1_12 * s, dp, uT, U0mol, alpha_i, beta_i, aTF, bpc, R,  x,  theta);
}





//--- Method 3 ---//
//==== Version 1 ===//



/*gpufun*/
void BentChannellingDevM3V1o12_track_local_particle(BentChannellingDevM3V1o12Data el, LocalParticle* part0) {

    double length = BentChannellingDevM3V1o12Data_get_length(el);

    //start_per_particle_block (part0->part)
        double x0     = LocalParticle_get_x(part);
        double theta0 = LocalParticle_get_xp(part);
        double x;
        double theta;
        double s = length;

        fM3O12v1(x0, theta0, s, dp, uT, U0mol, alpha_i, beta_i, aTF, bpc, R, &x, &theta);

        LocalParticle_set_x(part, x);
        LocalParticle_set_xp(part, theta);
    //end_per_particle_block
}





#endif /* XCOLL_BENT_CHANNELLINGDEV_M3V1O12_H */


