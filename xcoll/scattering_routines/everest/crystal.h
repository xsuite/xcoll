#include <math.h>
#include <stdio.h>
#include <stdlib.h>

double tlcut_cry = 0.0009982;

double aTF = 0.194e-10; // Screening function [m]
double dP  = 1.920e-10; // Distance between planes (110) [m]
double u1  = 0.075e-10; // Thermal vibrations amplitude

// pp cross-sections and parameters for energy dependence
double pptref_cry = 0.040;
double freeco_cry = 1.618;

// Processes
int proc_out         =  -1;     // Crystal not hit
int proc_AM          =   1;     // Amorphous
int proc_VR          =   2;     // Volume reflection
int proc_CH          =   3;     // Channeling
int proc_VC          =   4;     // Volume capture
int proc_absorbed    =   5;     // Absorption
int proc_DC          =   6;     // Dechanneling
int proc_pne         =   7;     // Proton-neutron elastic interaction
int proc_ppe         =   8;     // Proton-proton elastic interaction
int proc_diff        =   9;     // Single diffractive
int proc_ruth        =  10;     // Rutherford scattering
int proc_ch_absorbed =  15;     // Channeling followed by absorption
int proc_ch_pne      =  17;     // Channeling followed by proton-neutron elastic interaction
int proc_ch_ppe      =  18;     // Channeling followed by proton-proton elastic interaction
int proc_ch_diff     =  19;     // Channeling followed by single diffractive
int proc_ch_ruth     =  20;     // Channeling followed by Rutherford scattering
int proc_TRVR        = 100;     // Volume reflection in VR-AM transition region
int proc_TRAM        = 101;     // Amorphous in VR-AM transition region

int temp = 0;

double* movech(double nam, double dz, double x, double xp, double yp, double pc, double r, double rc, double rho, double anuc, double zatom, double emr, double hcut, double bnref, double csect, double csref0, double csref1, double csref5, double eUm, double collnt, double iProc) {

    //from .random import get_random, set_rutherford_parameters, get_random_ruth, get_random_gauss
    static double result[5];
    
    double pmae = 0.51099890;
    double pmap = 938.271998;

    double xp_in = xp;
    double yp_in = yp;
    double pc_in = pc;

    double cs[6] = {0.0,0.0,0.0,0.0,0.0,0.0};
    double cprob[6] = {0.0,0.0,0.0,0.0,0.0,0.0};
    double xran_cry[1] = {0.0};
    set_rutherford_parameters(zatom, emr, hcut);

    //New treatment of scattering routine based on standard sixtrack routine

    //Useful calculations for cross-section and event topology calculation
    double ecmsq  = ((2.*pmap)*1.0e-3)*pc;
    double xln15s = log(0.15*ecmsq);

    //New models, see Claudia's thesis
    double pptot = (0.041084 - 0.0023302*log(ecmsq)) + 0.00031514 * pow(log(ecmsq),2.);
    double ppel  = (11.7 - 1.59*log(ecmsq) + 0.134 * pow(log(ecmsq),2.))/1.0e3;
    double ppsd  = (4.3 + 0.3*log(ecmsq))/1.0e3;
    double bpp   = 7.156 + 1.439*log(sqrt(ecmsq));

    xran_cry[0] = get_random_ruth();

    //Rescale the total and inelastic cross-section accordigly to the average density seen
    double x_i = x;
    int np  = x_i/dP;    //Calculate in which crystalline plane the particle enters
    x_i = x_i - np*dP;    //Rescale the incoming x at the left crystalline plane
    x_i = x_i - (dP/2.); //Rescale the incoming x in the middle of crystalline planes

    double pv   = pow(pc,2.)/sqrt(pow(pc,2.) + pow((pmap*1.0e-3),2.))*1.0e9;          //Calculate pv=P/E
    double Ueff = eUm*((2.*x_i)/dP)*((2.*x_i)/dP) + pv*x_i/r; //Calculate effective potential
    double Et   = (pv*pow(xp,2.))/2. + Ueff;                            //Calculate transverse energy
    double Ec   = (eUm*(1.-rc/r))*(1.-rc/r);                  //Calculate critical energy in bent crystals

    //To avoid negative Et
    double xminU = ((-pow(dP,2.)*pc)*1.0e9)/(8.*eUm*r);
    double Umin  = fabs((eUm*((2.*xminU)/dP))*((2.*xminU)/dP) + pv*xminU/r);
    Et    = Et + Umin;
    Ec    = Ec + Umin;

    //Calculate min e max of the trajectory between crystalline planes
    double x_min = (-(dP/2.)*rc)/r - (dP/2.)*sqrt(Et/Ec);
    double x_max = (-(dP/2.)*rc)/r + (dP/2.)*sqrt(Et/Ec);

    //Change ref. frame and go back with 0 on the crystalline plane on the left
    x_min = x_min - dP/2.;
    x_max = x_max - dP/2.;

    //Calculate the "normal density" in m^-3
    double N_am  = ((rho*6.022e23)*1.0e6)/anuc;

    //Calculate atomic density at min and max of the trajectory oscillation
    // erf returns the error function of complex argument
    double rho_max = ((N_am*dP)/2.)*(erf(x_max/sqrt(2*pow(u1,2.))) - erf((dP-x_max)/sqrt(2*pow(u1,2.))));
    double rho_min = ((N_am*dP)/2.)*(erf(x_min/sqrt(2*pow(u1,2.))) - erf((dP-x_min)/sqrt(2*pow(u1,2.))));

    //"zero-approximation" of average nuclear density seen along the trajectory
    double avrrho  = (rho_max - rho_min)/(x_max - x_min);
    avrrho  = (2.*avrrho)/N_am;

    double csref_tot_rsc  = csref0*avrrho; //Rescaled total ref cs
    double csref_inel_rsc = csref1*avrrho; //Rescaled inelastic ref cs

    //Cross-section calculation
    double freep = freeco_cry * pow(anuc,1./3.);

    //compute pp and pn el+single diff contributions to cross-section (both added : quasi-elastic or qel later)
    cs[3] = freep*ppel;
    cs[4] = freep*ppsd;

    //correct TOT-CSec for energy dependence of qel
    //TOT CS is here without a Coulomb contribution
    cs[0] = csref_tot_rsc + freep*(pptot - pptref_cry);

    //Also correct inel-CS
    if(csref_tot_rsc == 0) {
        cs[1] = 0;
    }
        
    else {
        cs[1] = (csref_inel_rsc*cs[0])/csref_tot_rsc;
    }
        
    //Nuclear Elastic is TOT-inel-qel ( see definition in RPP)
    cs[2] = ((cs[0] - cs[1]) - cs[3]) - cs[4];
    cs[5] = csref5;

    //Now add Coulomb
    cs[0] = cs[0] + cs[5];

    //Calculate cumulative probability
    //cprob[6] = {0.0,0.0,0.0,0.0,0.0,0.0};
    cprob[5] = 1;
    
    if (cs[0] == 0) {
        int i;
        for (i = 1; i < 5; ++i) {
            cprob[i] = cprob[i-1];
        }
    }
        
    else {
        int i;
        for (i = 1; i < 5; ++i) {
            cprob[i] = cprob[i-1] + cs[i]/cs[0];
        }
    }

    //Multiple Coulomb Scattering
    xp = xp*1.0e3;
    yp = yp*1.0e3;

    //Turn on/off nuclear interactions
    if (nam == 0) {
        result[0] = x;
        result[1] = xp;
        result[2] = yp;
        result[3] = pc;
        result[4] = iProc;

        return result;
    }
    
    double nuc_cl_l;
    //Can nuclear interaction happen?
    //Rescaled nuclear collision length
    if (avrrho == 0) {
        nuc_cl_l = 1.0e6;
    }  
    else {
        nuc_cl_l = collnt/avrrho;
    }

    double zlm = -nuc_cl_l*log(get_random());

    //write(889,*) x_i,pv,Ueff,Et,Ec,N_am,avrrho,csref_tot_rsc,csref_inel_rsc,nuc_cl_l

    if (zlm < dz) {

        //Choose nuclear interaction
        double aran = get_random();
        int i = 1;
        double bn;
        double xm2;
        double bsd = 0.0;
        double teta;
        double tx;
        double tz;

        while (aran > cprob[i]) {
            i=i+1;
        }
        
        int ichoix = i;

        //Do the interaction
        double t = 0 ; //default value to cover ichoix=1
        
        if (ichoix==1) {
            iProc = proc_ch_absorbed; //deep inelastic, impinging p disappeared
        } 
            
        else if (ichoix==2) { //p-n elastic
            iProc = proc_ch_pne;
            bn    = (bnref*cs[0])/csref_tot_rsc;
            t     = -log(get_random())/bn;
        }

        else if (ichoix==3) { //p-p elastic
            iProc = proc_ch_ppe;
            t     = -log(get_random())/bpp;
        }

        else if (ichoix==4) { //Single diffractive
            iProc = proc_ch_diff;
            xm2   = exp(get_random()*xln15s);
            pc    = pc*(1 - xm2/ecmsq);

            if (xm2 < 2.) {
                bsd = 2*bpp;
            }
            else if (xm2 >= 2. && xm2 <= 5.) {
                bsd = ((106.0 - 17.0*xm2)*bpp)/36.0;
            }
            else if (xm2 > 5.) {
                bsd = (7*bpp)/12.0;
            }
            //end if
            t = -log(get_random())/bsd;
        }

        else { //(ichoix==5)
            iProc      = proc_ch_ruth;
            xran_cry[0] = get_random_ruth();
            t = xran_cry[0];
        }


        //Calculate the related kick -----------
        if (ichoix == 4) {
            teta = sqrt(t)/pc_in; //DIFF has changed PC!!!
        }
        else {
            teta = sqrt(t)/pc;
        }

        tx = teta*get_random_gauss()*1.0e3;
        tz = teta*get_random_gauss()*1.0e3;

        //Change p angle
        xp = xp + tx;
        yp = yp + tz;
    }

    xp = xp/1.0e3;
    yp = yp/1.0e3;

    result[0] = x;
    result[1] = xp;
    result[2] = yp;
    result[3] = pc;
    result[4] = iProc;

    return result;

}


double* moveam(double nam, double dz, double dei, double dlr, double xp, double yp, double pc, double anuc, double zatom, double emr, double hcut, double bnref, double csref0, double csref1, double csref5, double collnt, double iProc) {

    // from .random import get_random, set_rutherford_parameters, get_random_ruth, get_random_gauss
    static double result[4];

    double xp_in = xp;
    double yp_in = yp;
    double pc_in = pc;

    double pmae = 0.51099890;
    double pmap = 938.271998;

    double cs[6] = {0.0,0.0,0.0,0.0,0.0,0.0};
    double cprob[6] = {0.0,0.0,0.0,0.0,0.0,0.0};

    // New treatment of scattering routine based on standard sixtrack routine
    // useful calculations for cross-section and event topology calculation
    double ecmsq  = 2.*pmap*1.0e-3*pc;
    double xln15s = log(0.15*ecmsq);

    // New models, see Claudia's thesis
    double pptot = 0.041084 - 0.0023302*log(ecmsq) + 0.00031514 * pow(log(ecmsq),2.);
    double ppel  = (11.7 - 1.59*log(ecmsq) + 0.134 * pow(log(ecmsq),2.))/1.0e3;
    double ppsd  = (4.3 + 0.3*log(ecmsq))/1.0e3;
    double bpp   = 7.156 + 1.439*log(sqrt(ecmsq));

    // // Distribution for Ruth. scatt.
    set_rutherford_parameters(zatom, emr, hcut);

    // Cross-section calculation
    // freep: number of nucleons involved in single scattering
    double freep = freeco_cry * pow(anuc,1./3.);

    // Compute pp and pn el+single diff contributions to cross-section (both added : quasi-elastic or qel later)
    cs[3] = freep*ppel;
    cs[4] = freep*ppsd;

    // Correct TOT-CSec for energy dependence of qel
    // TOT CS is here without a Coulomb contribution
    cs[0] = csref0 + freep*(pptot - pptref_cry);
    double bn = (bnref*cs[0])/csref0;

    // Also correct inel-CS
    cs[1] = (csref1*cs[0])/csref0;

    // Nuclear Elastic is TOT-inel-qel ( see definition in RPP)
    cs[2] = ((cs[0] - cs[1]) - cs[3]) - cs[4];
    cs[5] = csref5;

    // Now add Coulomb
    cs[0] = cs[0] + cs[5];

    // Calculate cumulative probability
    cprob[5] = 1;

    int i;
    for (i = 1; i < 5; ++i) {
        cprob[i] = cprob[i-1] + cs[i]/cs[0];
    }

    // Multiple Coulomb Scattering
    xp  = xp*1.0e3;
    yp  = yp*1.0e3;
    pc  = pc - dei*dz; // Energy lost because of ionization process[GeV]

    double dya   = (13.6/pc)*sqrt(dz/dlr); // RMS of coloumb scattering MCS (mrad)
    double kxmcs = dya*get_random_gauss();
    double kymcs = dya*get_random_gauss();

    xp = xp+kxmcs;
    yp = yp+kymcs;

    if (nam == 0) {

        result[0] = xp;
        result[1] = yp;
        result[2] = pc;
        result[3] = iProc;
        return result; // Turn on/off nuclear interactions
    }

    // Can nuclear interaction happen?
    double zlm = -collnt*log(get_random());

    if (zlm < dz) {
        // Choose nuclear interaction
        double aran = get_random();
        int i=1;

        while (aran > cprob[i]) {
            i = i+1;
            //goto 10
        }

        int ichoix = i;

        // Do the interaction
        double t = 0 ;// default value to cover ichoix=1
        if (ichoix==1) {
        //case(1) // Deep inelastic, impinging p disappeared
            iProc = proc_absorbed;
        }

        else if (ichoix==2) { // p-n elastic
            iProc = proc_pne;
            t     = -log(get_random())/bn;
        }

        else if (ichoix==3) { // p-p elastic
            iProc = proc_ppe;
            t     = -log(get_random())/bpp;
        }

        else if (ichoix==4) { // Single diffractive
            iProc = proc_diff;
            double xm2 = exp(get_random()*xln15s);
            double bsd = 0.0;
            pc    = pc*(1 - xm2/ecmsq);

            if(xm2 < 2) {
                bsd = 2*bpp;
            }
            else if(xm2 >= 2 && xm2 <= 5) {
                bsd = ((106.0 - 17.0*xm2)*bpp)/36.0;
            }
            else if(xm2 > 5) {
                bsd = 7.0*bpp/12.0;
            }
        
            t = -log(get_random())/bsd;
        }

        else { //(ichoix==5)
            iProc = proc_ruth;
        // in python: t = get_random_ruth()
            //length_cry = 1
            t = get_random_ruth();
        }

        // end select

        double teta;

        // Calculate the related kick
        if(ichoix == 4) {
            teta = sqrt(t)/pc_in ;// DIFF has changed PC
        }

        else {
            teta = sqrt(t)/pc;
        }

        double tx = teta*get_random_gauss()*1.0e3;
        double tz = teta*get_random_gauss()*1.0e3;

        // Change p angle
        xp = xp + tx;
        yp = yp + tz;
        
    }

    xp = xp/1.0e3;
    yp = yp/1.0e3;

    result[0] = xp;
    result[1] = yp;
    result[2] = pc;
    result[3] = iProc;

    return result; // Turn on/off nuclear interactions

}