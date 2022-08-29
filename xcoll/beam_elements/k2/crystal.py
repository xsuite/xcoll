import numpy as np

def crystal(*,x,xp,z,zp,s,p,x0,xp0,zlm,simp,isimp,val_part_hit,val_part_abs,val_part_impact,val_part_indiv,c_length,exenergy,rho,anuc,zatom,emr,
            dlri,dlyi,ai,eum,collnt,hcut,bnref,csref0,csref1,csref4,csref5,csect,nhit,nabs,nmat):


    from .k2_random import get_random
    
    try:
        import xcoll.beam_elements.pyk2 as pyk2
    except ImportError:
        raise Exception("Error: Failed importing pyK2 (did you compile?). Cannot track.") 

    val_part_hit = np.array(val_part_hit)
    val_part_abs = np.array(val_part_abs)
    val_part_impact = np.array(val_part_impact)
    val_part_indiv = np.array(val_part_indiv)

    exenergy = np.array(exenergy)
    anuc = np.array(anuc)
    zatom = np.array(zatom)
    emr = np.array(emr)
    rho = np.array(rho)
    hcut = np.array(hcut)
    bnref = np.array(bnref)
    simp = np.array(simp)

    csref0 = np.array(csref0)
    csref1 = np.array(csref1)
    csref4 = np.array(csref4)
    csref5 = np.array(csref5)

    csect = np.array(csect)
    dlri = np.array(dlri)
    dlyi = np.array(dlyi)
    eum = np.array(eum)
    ai = np.array(ai)
    collnt = np.array(collnt)
    c_length = np.array(c_length)

    nhit = np.array(nhit)
    nabs = np.array(nabs)
    isimp = np.array(isimp)
    s = np.array(s)
    zlm = np.array(zlm)
    x0 = np.array(x0)
    xp0 = np.array(xp0)

    x = np.array(x)
    xp = np.array(xp)
    z = np.array(z)
    zp = np.array(zp)
    p = np.array(p)
    nmat = np.array(nmat)
    csect = np.array(csect)

    ########################################################

    pyk2.pyk2_docrystal(x=x,
                        xp=xp,
                        z=z,
                        zp=zp,
                        s=s,
                        p=p,
                        x0=x0,
                        xp0=xp0, 
                        zlm=zlm,
                        s_imp=simp,
                        isimp=isimp,
                        nhit=nhit,
                        nabs=nabs,
                        lhit=val_part_hit,
                        part_abs=val_part_abs,
                        impact=val_part_impact,
                        indiv=val_part_indiv,
                        c_length=c_length,
                        exenergy=exenergy,
                        rho=rho,
                        anuc=anuc,
                        zatom=zatom,
                        emr=emr,
                        dlri=dlri,
                        dlyi=dlyi,
                        ai=ai,
                        eum=eum,
                        collnt=collnt,
                        hcut=hcut,
                        csref0=csref0,
                        csref1=csref1,
                        csref4=csref4,
                        csref5=csref5,
                        nmat=nmat,
                        bnref=bnref,
                        csect=csect
                        )

    ########################################################
    
    
    return val_part_hit, val_part_abs, val_part_impact, val_part_indiv, nhit, nabs, simp, isimp, s, zlm, x, xp, x0, xp0, z, zp, p
