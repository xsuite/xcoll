float mcs(float s, float mc_radl, float mc_zlm1, float mc_p0, float mc_x, float mc_xp, float mc_z, float mc_zp, float mc_dpop) {

    float theta    = 13.6e-3/(mc_p0*(1+mc_dpop)); // dpop   = (p - p0)/p0;
    float h   = 0.001;
    float dh  = 0.0001;
    float bn0 = 0.4330127019;
    float rlen0 = mc_zlm1/mc_radl;
    float rlen  = rlen0;

    mc_x     = (mc_x/theta)/mc_radl;
    mc_xp    = mc_xp/theta;
    mc_z     = (mc_z/theta)/mc_radl;
    mc_zp    = mc_zp/theta;


    while (1) {
        
        float ae = bn0*mc_x;
        float be = bn0*mc_xp;
        

        // #######################################
        // ae = np.array(ae, dtype=np.float64)
        // be = np.array(be, dtype=np.float64)
        // dh = np.array(dh, dtype=np.float64)
        // rlen = np.array(rlen, dtype=np.float64)
        // s = np.array(s, dtype=np.float64)
        // #######################################
        float s = soln3(ae,be,dh,rlen,s);

        if (s < h) {
            s = h;
        }

        mc_x, mc_xp = scamcs(mc_x,mc_xp,s);

        if (mc_x <= 0) {
            s = (rlen0-rlen)+s;
            break; // go to 20
        }

        if ((s+dh) >= rlen) {
            s = rlen0;
            break; // go to 20
        }
        // go to 10
        rlen = rlen-s;
    }

    mc_z, mc_zp = scamcs(mc_z,mc_zp,s);

    s  = s*mc_radl;
    mc_x  = (mc_x*theta)*mc_radl;
    mc_xp = mc_xp*theta;
    mc_z  = (mc_z*theta)*mc_radl;
    mc_zp = mc_zp*theta;

    return s, mc_x, mc_xp, mc_z, mc_zp, mc_dpop;


}

int square(int i) {
	return i * i;
}
