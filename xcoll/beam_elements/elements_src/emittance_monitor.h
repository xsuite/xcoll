// copyright ############################### #
// This file is part of the Xcoll package.   #
// Copyright (c) CERN, 2024.                 #
// ######################################### #

#ifndef XCOLL_EMITTANCE_MONITOR_H
#define XCOLL_EMITTANCE_MONITOR_H

#ifdef XO_CONTEXT_CPU
#include <math.h>    // for round
#include <stdint.h>  // for int64_t etc
#endif  // XO_CONTEXT_CPU

#ifndef C_LIGHT
#define C_LIGHT 299792458.0
#endif


/*gpufun*/
void EmittanceMonitor_track_local_particle(EmittanceMonitorData el, LocalParticle* part0){
    int64_t const start_at_turn = EmittanceMonitorData_get_start_at_turn(el);
    int64_t const part_id_start = EmittanceMonitorData_get_part_id_start(el);
    int64_t const part_id_end   = EmittanceMonitorData_get_part_id_end(el);
    double const frev = EmittanceMonitorData_get_frev(el);
    double const sampling_frequency = EmittanceMonitorData_get_sampling_frequency(el);

    EmittanceMonitorRecord record = EmittanceMonitorData_getp_data(el);

    int64_t max_slot = EmittanceMonitorRecord_len_count(record);

    int const plane_selector = EmittanceMonitorData_get__plane_selector(el);
    int const monitor_horizontal   =  plane_selector % 2;
    int const monitor_vertical     = (plane_selector >> 1) % 2;
    int const monitor_longitudinal = (plane_selector >> 2) % 2;

    //start_per_particle_block(part0->part)
        int64_t particle_id = LocalParticle_get_particle_id(part);
        if (part_id_end < 0 || (part_id_start <= particle_id && particle_id < part_id_end)){

            // zeta is the absolute path length deviation from the reference particle: zeta = (s - beta0*c*t)
            // but without limits, i.e. it can exceed the circumference (for coasting beams)
            // as the particle falls behind or overtakes the reference particle
            double const zeta = LocalParticle_get_zeta(part);
            double const beta0 = LocalParticle_get_beta0(part);
            double const at_turn = LocalParticle_get_at_turn(part);

            double x =0;
            double px = 0;
            double y = 0;
            double py = 0;
            double pzeta = 0;

            // compute sample index
            int64_t slot = round(sampling_frequency * ( (at_turn-start_at_turn)/frev - zeta/beta0/C_LIGHT ));

            if (slot >= 0 && slot < max_slot){
                /*gpuglmem*/ double *count = EmittanceMonitorRecord_getp1_count(record, slot); atomicAdd(count, 1);
                EmittanceMonitorRecord_set_cached(record, slot, 0);
                EmittanceMonitorRecord_set_cached_modes(record, slot, 0);

                if (monitor_horizontal){
                    x  = LocalParticle_get_x(part);
                    px = LocalParticle_get_px(part);
                    /*gpuglmem*/ double *x_sum1     = EmittanceMonitorRecord_getp1_x_sum1(record, slot);     atomicAdd(x_sum1, x);
                    /*gpuglmem*/ double *px_sum1    = EmittanceMonitorRecord_getp1_px_sum1(record, slot);    atomicAdd(px_sum1, px);
                    /*gpuglmem*/ double *x_x_sum2   = EmittanceMonitorRecord_getp1_x_x_sum2(record, slot);   atomicAdd(x_x_sum2, x*x);
                    /*gpuglmem*/ double *x_px_sum2  = EmittanceMonitorRecord_getp1_x_px_sum2(record, slot);  atomicAdd(x_px_sum2, x*px);
                    /*gpuglmem*/ double *px_px_sum2 = EmittanceMonitorRecord_getp1_px_px_sum2(record, slot); atomicAdd(px_px_sum2, px*px);
                }

                if (monitor_vertical){
                    y  = LocalParticle_get_y(part);
                    py = LocalParticle_get_py(part);
                    /*gpuglmem*/ double *y_sum1     = EmittanceMonitorRecord_getp1_y_sum1(record, slot);     atomicAdd(y_sum1, y);
                    /*gpuglmem*/ double *py_sum1    = EmittanceMonitorRecord_getp1_py_sum1(record, slot);    atomicAdd(py_sum1, py);
                    /*gpuglmem*/ double *y_y_sum2   = EmittanceMonitorRecord_getp1_y_y_sum2(record, slot);   atomicAdd(y_y_sum2, y*y);
                    /*gpuglmem*/ double *y_py_sum2  = EmittanceMonitorRecord_getp1_y_py_sum2(record, slot);  atomicAdd(y_py_sum2, y*py);
                    /*gpuglmem*/ double *py_py_sum2 = EmittanceMonitorRecord_getp1_py_py_sum2(record, slot); atomicAdd(py_py_sum2, py*py);
                }

                if (monitor_longitudinal){
                    double const ptau = LocalParticle_get_ptau(part);
                    pzeta = ptau/beta0;
                    /*gpuglmem*/ double *zeta_sum1        = EmittanceMonitorRecord_getp1_zeta_sum1(record, slot);        atomicAdd(zeta_sum1, zeta);
                    /*gpuglmem*/ double *pzeta_sum1       = EmittanceMonitorRecord_getp1_pzeta_sum1(record, slot);       atomicAdd(pzeta_sum1, pzeta);
                    /*gpuglmem*/ double *zeta_zeta_sum2   = EmittanceMonitorRecord_getp1_zeta_zeta_sum2(record, slot);   atomicAdd(zeta_zeta_sum2, zeta*zeta);
                    /*gpuglmem*/ double *zeta_pzeta_sum2  = EmittanceMonitorRecord_getp1_zeta_pzeta_sum2(record, slot);  atomicAdd(zeta_pzeta_sum2, zeta*pzeta);
                    /*gpuglmem*/ double *pzeta_pzeta_sum2 = EmittanceMonitorRecord_getp1_pzeta_pzeta_sum2(record, slot); atomicAdd(pzeta_pzeta_sum2, pzeta*pzeta);
                }

                if (monitor_horizontal && monitor_vertical){
                    /*gpuglmem*/ double *x_y_sum2   = EmittanceMonitorRecord_getp1_x_y_sum2(record, slot);   atomicAdd(x_y_sum2, x*y);
                    /*gpuglmem*/ double *x_py_sum2  = EmittanceMonitorRecord_getp1_x_py_sum2(record, slot);  atomicAdd(x_py_sum2, x*py);
                    /*gpuglmem*/ double *px_y_sum2  = EmittanceMonitorRecord_getp1_px_y_sum2(record, slot);  atomicAdd(px_y_sum2, px*y);
                    /*gpuglmem*/ double *px_py_sum2 = EmittanceMonitorRecord_getp1_px_py_sum2(record, slot); atomicAdd(px_py_sum2, px*py);
                }

                if (monitor_horizontal && monitor_longitudinal){
                    /*gpuglmem*/ double *x_zeta_sum2   = EmittanceMonitorRecord_getp1_x_zeta_sum2(record, slot);   atomicAdd(x_zeta_sum2, x*zeta);
                    /*gpuglmem*/ double *x_pzeta_sum2  = EmittanceMonitorRecord_getp1_x_pzeta_sum2(record, slot);  atomicAdd(x_pzeta_sum2, x*pzeta);
                    /*gpuglmem*/ double *px_zeta_sum2  = EmittanceMonitorRecord_getp1_px_zeta_sum2(record, slot);  atomicAdd(px_zeta_sum2, px*zeta);
                    /*gpuglmem*/ double *px_pzeta_sum2 = EmittanceMonitorRecord_getp1_px_pzeta_sum2(record, slot); atomicAdd(px_pzeta_sum2, px*pzeta);
                }

                if (monitor_vertical && monitor_longitudinal){
                    /*gpuglmem*/ double *y_zeta_sum2   = EmittanceMonitorRecord_getp1_y_zeta_sum2(record, slot);   atomicAdd(y_zeta_sum2, y*zeta);
                    /*gpuglmem*/ double *y_pzeta_sum2  = EmittanceMonitorRecord_getp1_y_pzeta_sum2(record, slot);  atomicAdd(y_pzeta_sum2, y*pzeta);
                    /*gpuglmem*/ double *py_zeta_sum2  = EmittanceMonitorRecord_getp1_py_zeta_sum2(record, slot);  atomicAdd(py_zeta_sum2, py*zeta);
                    /*gpuglmem*/ double *py_pzeta_sum2 = EmittanceMonitorRecord_getp1_py_pzeta_sum2(record, slot); atomicAdd(py_pzeta_sum2, py*pzeta);
                }
            }
        }
	//end_per_particle_block
}

#endif /* XCOLL_EMITTANCE_MONITOR_H */