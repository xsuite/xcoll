#include <math.h>

double iterat(double a, double b, double dh, double s) {

    double ds = s;

    while (1) {
        
        ds = ds*0.5;

        if (pow(s,3) < pow((a+b*s),2)) {
            s = s+ds;
        } else {
            s = s-ds; 
        }

        if (ds < dh) {
            break;
        } else { 
            continue;
        }
    }

    return s;

}
  

double soln3(double a, double b, double dh, double smax, double s) {


    if (b == 0) {
        s = pow(a,0.6666666666666667);

        if (s > smax) {
            s = smax;
        }
        return s;
    }

    if (a == 0) {    
        if (b > 0) {
            s = pow(b,2);
        } else {
            s = 0;
        }
        if (s > smax) {
            s = smax;
        }
        return s;
    }
        
    if (b > 0) {

        if (pow(smax,3) <= pow((a + b*smax),2)) {
            s = smax;
            return s;
        } else {
            s = smax*0.5;
            s = iterat(a,b,dh,s);
        }
    } else {
        double c = (-1*a)/b;
        if (smax < c) {
            if ((pow(smax,3)) <= pow((a + b*smax),2)) {
                s = smax;
                return s;
            } else {
                s = smax*0.5;
                s = iterat(a,b,dh,s);
            }
        } else {
            s = c*0.5;
            s = iterat(a,b,dh,s);
        }
    }
   
    return s;

}


// double mcs(double s, double mc_radl, double mc_zlm1, double mc_p0, double mc_x, double mc_xp, double mc_z, double mc_zp, double mc_dpop) {

//     double theta    = 13.6e-3/(mc_p0*(1+mc_dpop)); // dpop   = (p - p0)/p0;
//     double h   = 0.001;
//     double dh  = 0.0001;
//     double bn0 = 0.4330127019;
//     double rlen0 = mc_zlm1/mc_radl;
//     double rlen  = rlen0;

//     mc_x     = (mc_x/theta)/mc_radl;
//     mc_xp    = mc_xp/theta;
//     mc_z     = (mc_z/theta)/mc_radl;
//     mc_zp    = mc_zp/theta;


//     while (1) {
        
//         double ae = bn0*mc_x;
//         double be = bn0*mc_xp;

//         // #######################################
//         // ae = np.array(ae, dtype=np.double64)
//         // be = np.array(be, dtype=np.double64)
//         // dh = np.array(dh, dtype=np.double64)
//         // rlen = np.array(rlen, dtype=np.double64)
//         // s = np.array(s, dtype=np.double64)
//         // #######################################
//         double s = soln3(ae,be,dh,rlen,s);

//         if (s < h) {
//             s = h;
//         }

//         mc_x, mc_xp = scamcs(mc_x,mc_xp,s);

//         if (mc_x <= 0) {
//             s = (rlen0-rlen)+s;
//             break; // go to 20
//         }

//         if ((s+dh) >= rlen) {
//             s = rlen0;
//             break; // go to 20
//         }
//         // go to 10
//         rlen = rlen-s;
//     }

//     mc_z, mc_zp = scamcs(mc_z,mc_zp,s);

//     s  = s*mc_radl;
//     mc_x  = (mc_x*theta)*mc_radl;
//     mc_xp = mc_xp*theta;
//     mc_z  = (mc_z*theta)*mc_radl;
//     mc_zp = mc_zp*theta;

//     return s, mc_x, mc_xp, mc_z, mc_zp, mc_dpop ;


// }


// int my_square(int i) {

//     return i*i;
// }