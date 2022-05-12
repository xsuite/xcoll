import numpy as np

def crystal(*, val_part_hit, val_part_abs, val_part_impact, val_part_indiv, run_exenergy, run_anuc, run_zatom, run_emr, run_rho, run_hcut, run_bnref, run_csref0, run_csref1, run_csref4, run_csref5, run_dlri, run_dlyi, run_eum, run_ai, run_collnt, run_bn, c_length, nhit, nabs, isimp, s, zlm, x, xp, xp_in0, z, zp, p, x_in0):

    try:
        import xcoll.beam_elements.pyk2 as pyk2
    except ImportError:
        raise Exception("Error: Failed importing pyK2 (did you compile?). Cannot track.")
    
    val_part_hit = np.array(val_part_hit, dtype=np.int64)
    val_part_abs = np.array(val_part_abs, dtype=np.int64)
    val_part_impact = np.array(val_part_impact, dtype=np.float64)
    val_part_indiv = np.array(val_part_indiv, dtype=np.float64)
    run_exenergy = np.array(run_exenergy, dtype=np.float64)
    run_bn = np.array(run_bn, dtype=np.float64)
    nhit = np.array(nhit, dtype=np.int64)
    nabs = np.array(nabs, dtype=np.int64)
    isimp = np.array(isimp)
    s = np.array(s, dtype=np.float64)
    zlm = np.array(zlm, dtype=np.float64)
    x = np.array(x, dtype=np.float64)
    xp = np.array(xp, dtype=np.float64)
    xp_in0 = np.array(xp_in0, dtype=np.float64)
    z = np.array(z, dtype=np.float64)
    zp = np.array(zp, dtype=np.float64)
    p = np.array(p, dtype=np.float64)
    x_in = np.array(x_in, dtype=np.float64)

            ##########################################

    pyk2.pyk2_crystal(
                    val_part_hit=val_part_hit,
                    val_part_abs=val_part_abs,
                    val_part_impact=val_part_impact,
                    val_part_indiv=val_part_indiv,
                    run_exenergy=run_exenergy,
                    run_anuc=run_anuc,
                    run_zatom=run_zatom,
                    run_emr=run_emr,
                    run_rho=run_rho,
                    run_hcut=run_hcut,
                    run_bnref=run_bnref,
                    run_csref0=run_csref0,
                    run_csref1=run_csref1,
                    run_csref4=run_csref4,
                    run_csref5=run_csref5,
                    run_dlri=run_dlri,
                    run_dlyi=run_dlyi,
                    run_eum=run_eum,
                    run_ai=run_ai,
                    run_collnt=run_collnt,
                    run_bn=run_bn,
                    c_length=c_length,
                    nhit=nhit,
                    nabs=nabs,
                    isimp=isimp,
                    s=s,
                    zlm=zlm,
                    x=x,
                    xp=xp,
                    xp_in0=xp_in0,
                    z=z,
                    zp=zp,
                    p=p,
                    x_in0=x_in0
                    )

    return val_part_hit, val_part_abs, val_part_impact, val_part_indiv, run_exenergy, run_bn, nhit, nabs, isimp, s, zlm, x, xp, xp_in0, z, zp, p, x_in