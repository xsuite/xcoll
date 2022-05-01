import numpy as np
import pandas as pd
from scipy.constants import c as clight
import xpart as xp


# Note: SixTrack initial.dat is with respect to the closed orbit when using TRAC,
#       but in the lab frame when using SIMU
def particles_to_sixtrack_initial(part, filename):
    with open(filename, 'w') as pf:
        part_xp = np.array(part.px*part.rpp)
        part_yp = np.array(part.py*part.rpp)
        charge = [ round(q) for q in part.charge_ratio*part.q0 ]
        mass_ratio = [ round(m) for m in part.charge_ratio/part.chi ]
        mass = np.array(part.mass0*part.charge_ratio/part.chi)
        part_p = np.array((1+part.delta)*part.p0c)
        # sigma = - beta0*c*dt
        part_dt = - np.array(part.zeta/part.rvv/part.beta0) / clight
        data = pd.DataFrame({
            'particle ID':   list(range(1, len(part.state)+1)),
            'parent ID':     list(range(1, len(part.state)+1)),
            'weight':        1,      # unused
            'x [m]':         np.array(part.x),
            'y [m]':         np.array(part.y),
            'z':             1,      # unused
            'xp [rad]':      part_xp,
            'yp [rad]':      part_yp,
            'zp':            1,      # unused
            'mass number':   mass_ratio,        # This is not correct! If the parent particle
                                                # is an ion, mass_ratio will not approximate
                                                # the mass number. Works for protons though
            'atomic number': charge,
            'mass [GeV/c2]': mass*1e-9,
            'p [GeV/c2]':    part_p*1e-9,
            'delta t [ms]':  part_dt*1e3
        })
        data.to_csv(pf, sep=' ', header=False, index=False)
        
        
        
def sixtrack_initial_to_particles(sixtrack_input_file, p0c, *, mass0=xp.PROTON_MASS_EV, q0=1):
    data = pd.read_csv(sixtrack_input_file, delim_whitespace=True, header=None,
                   names=['pid','parent','weight','x','y','z','xp','yp','zp', 'A', 'Z', 'm', 'p', 'dt'])

    charge_ratio = np.array( data.Z/q0 )
    chi = np.array( charge_ratio*mass0/data.m/1e9 )
    x = np.array(data.x)
    y = np.array(data.y)
    px = np.array(data.xp*data.p*1e9/p0c)
    py = np.array(data.yp*data.p*1e9/p0c)
    beta = data.p / np.sqrt( data.p**2 + data.m**2 )
    zeta = np.array(-data.dt*1e-3*clight * beta)
    delta = np.array( (data.p*1e9 - p0c) / p0c )

    return xp.Particles(
        mass0=mass0, q0=q0, p0c=p0c, charge_ratio=charge_ratio, chi=chi,
        x=x, y=y, px=px, py=py, zeta=zeta, delta=delta
    )


def sixtrack_dump2_to_particles(sixtrack_dump_file, p0c, *, mass0=xp.PROTON_MASS_EV):
    header = pd.read_csv(sixtrack_dump_file, delim_whitespace=True, header=None, nrows=2).loc[1][1:11].values
    data = pd.read_csv(sixtrack_dump_file, delim_whitespace=True, header=None, skiprows=2, names=header)
    data['state'] = 1
    mask = (data['x[mm]']==0).values & (data['xp[mrad]']==0).values & (data['y[mm]']==0).values & \
            (data['yp[mrad]']==0).values & (data['sigma[mm]']==0).values & (data['(E-E0)/E0[1]']==0).values
    data.loc[mask,'state'] = 0

    # delta+1 = pc/p0c = sqrt(e**2 + m**2)/p0c = sqrt[(e/e0)**2 * (p0c**2 + m0**2) + m**2] / p0c
    #         = sqrt[(e/e0)**2 * (1 + (m0/p0c)**2) + (m/p0c)**2]
    e_over_e0 = (data['(E-E0)/E0[1]'] + 1)
    m0_over_p0 = mass0 / p0c
    delta_plus_1 = np.sqrt(e_over_e0**2 * (1+m0_over_p0**2) - m0_over_p0**2)

    x = np.array(data['x[mm]'])*1e-3
    y = np.array(data['y[mm]'])*1e-3
    px = np.array(data['xp[mrad]']*1e-3*delta_plus_1)
    py = np.array(data['yp[mrad]']*1e-3*delta_plus_1)
    # zeta = beta / beta0 * sigma = (pc/e) / (p0c/e0) * sigma = (pc/p0c) / (e/e0) * sigma
    sigma_to_zeta = np.sqrt(1 + m0_over_p0**2 - m0_over_p0**2 / e_over_e0**2)
    zeta = np.array(sigma_to_zeta * data['sigma[mm]']*1e-3)
    delta = np.array(delta_plus_1 - 1)

    return xp.Particles(
        mass0=mass0, p0c=p0c, x=x, y=y, px=px, py=py, zeta=zeta, delta=delta, state=data['state'].values
    )
